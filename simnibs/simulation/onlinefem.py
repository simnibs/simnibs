import os
import time
import copy
import logging
import platform
import warnings

import numpy as np
from scipy import sparse

from simnibs import __version__
from simnibs.simulation.tms_coil.tms_coil import TmsCoil
from simnibs.mesh_tools import mesh_io
from simnibs.mesh_tools.mesh_io import Msh, read_msh
from simnibs.utils import simnibs_logger
from simnibs.utils.simnibs_logger import logger
from simnibs.utils.file_finder import Templates, SubjectFiles

from .fem import get_dirichlet_node_index_cog, TDCSFEMNeumann, TMSFEM
from .sim_struct import SimuList

# TODO: import here as numba interacts badly with pyqt (GUI). remove comment
# here once this is resolved.
# from simnibs.simulation.numba_fem_utils import sumf, sumf2, node2elmf, sumf3, postp, postp_mag, spmatmul

class OnlineFEM:
    """
    OnlineFEM class for fast FEM calculations using the Pardiso solver of the MKL library

    !!! NOTE !!!
    TMS:
        The mesh used in the OnlineFem and in the FemTargetPointCloud have to have the same node coordinates and the same node ordering.
        They also have to have the same elements and the same element ordering.
    TES:
        The mesh used in the OnlineFem and in the FemTargetPointCloud have to have the same node coordinates and the same node ordering.

    Parameters
    ----------
    mesh : Msh object or str
        Head mesh or path to it.
    method : str
        Specify simulation type ('TMS' or 'TES')
    roi : list of FemTargetPointCloud object [n_roi]
        Target coordinates for the E-field calculation.
    anisotropy_type : str
        Type of anisotropy for simulation ('scalar', 'vn', 'mc')
    solver_options : str
        Options for the FEM solver (Node or "pardiso))
    fn_logger : str
        Filename of log file (optional)
    useElements : bool
        True: interpolate the dadt field using positions (coordinates) of elements centroids
        False: interpolate the dadt field using positions (coordinates) of nodes
    fn_coil : str
        Path to coil file (.tcd or .ccd or .nii)
    dataType : list of int [n_roi], optional, default: 0
        Calc. magn. of e-field for dataType=0 otherwise return Ex, Ey, Ez. Defined for each ROI.
    coil : Coil instance
        TMS coil
    electrode : ElectrodeArrayPair or CircularArray instance
        TES electrode
    dirichlet_node : int
        Index of dirichlet node (indexing starting with 1)
    """
    def __init__(self, mesh, method, roi, anisotropy_type="scalar", solver_options="pardiso", fn_logger=None,
                 useElements=True, fn_coil=None, dataType=0, coil=None, electrode=None, dirichlet_node=None,
                 solver_loglevel=logging.DEBUG, cpus=None, cond=None):
        """
        Constructor of the OnlineFEM class
        """

        self.method = method                        # 'TES' or 'TMS'
        self.dataType = dataType                    # calc. magn. of e-field for dataType=0 otherwise return Ex, Ey, Ez
        self.cond = cond
        self.anisotropy_type = anisotropy_type      # 'scalar', 'dir', 'vn', 'mc'
        self.A = None                               # stiffness matrix
        self.b = None                               # rhs
        self.v = None                               # electric potential
        self.solver = None                          # solver
        self.fn_coil = fn_coil                      # filename of TMS coil (.nii)
        self.roi = roi                              # list of ROI instances
        self.fn_logger = fn_logger                  # either name of log file or True (using standard simnibs logger in later case)
        self.ff_templates = Templates()             # initialize file finder templates
        self.force_integrals = None                 # precomputed force integrals for rhs
        self.coil = coil                            # TMS coil
        self.electrode = electrode                  # TES electrode
        self.n_iter_dirichlet_correction = 0        # number of iterations required for dirichlet correction
        self.aniso_maxratio = 10                    # not used (included to use methods of the SimuList class)
        self.name = None                            # not used (included to use methods of the SimuList class)
        self.aniso_maxcond = 2                      # not used (included to use methods of the SimuList class)

        if type(self.roi) is not list:
            self.roi = [self.roi]

        self.n_roi = len(self.roi)
        self.e = [0 for _ in range(self.n_roi)]
        self.dadt = None

        # True: interpolate the dadt field using positions (coordinates) of elements centroids
        # False: interpolate the dadt field using positions (coordinates) of nodes
        self.useElements = useElements

        if platform.system() == "Darwin" and solver_options == "pardiso":
            warnings.warn("solver `pardiso` was specified but MKL PARDISO is not available for Darwin systems; falling back to MUMPS.")
            solver_options = "mumps"

        self.solver_options = solver_options

        # creating logger
        if fn_logger is not None:
            self._logging = True
            self._set_logger(fn_logger)
        else:
            self._logging = False

        self.solver_loglevel = solver_loglevel

        if not cpus is None:
            logger.info(f"Attempting to limit number of threads to {cpus}")
            if solver_options == "pardiso":
                os.environ["MKL_DOMAIN_NUM_THREADS"] = (
                    f"MKL_DOMAIN_PARDISO={str(int(cpus))}"
                )
                # Note: Using MKL_DOMAIN_PARDISO means that only PARDISO is affected, VML, BLAS, LAPACK is potentially still more threads otherwise change to MKL_NUM_THREADS OR MKL_DOMAIN_ALL
            else: #mainly for solver_options == "mumps" but PETsc can also use MPI
                os.environ["MKL_DOMAIN_NUM_THREADS"] = f"OMP_NUM_THREADS={str(int(cpus))}"
            from numba import set_num_threads

            #below is to limit numba threads
            set_num_threads(int(cpus))
            from numba import get_num_threads
            logger.info(f"Numba reports {get_num_threads()} threads available")

        # read mesh or store in self
        if type(mesh) is str:
            self.mesh = mesh_io.read_msh(mesh)
        elif type(mesh) == Msh:
            self.mesh = mesh
        else:
            raise TypeError("'mesh' parameter has to be either path to .msh file or SimNIBS mesh object.")

        # get subject specific filenames
        if self.mesh.fn:
            ff_subject = SubjectFiles(fnamehead=self.mesh.fn)
            self.fn_tensor_nifti = ff_subject.tensor_file
        else:
            self.fn_tensor_nifti = None

        # get Dirichlet node index where V=0 (close to center of gravity of mesh but not on a surface)
        if dirichlet_node is None:
            self.dirichlet_node = get_dirichlet_node_index_cog(mesh=self.mesh)  # self.mesh.nodes.node_number[self.mesh.nodes.node_coord[:, 2].argmin()]
        else:
            self.dirichlet_node = dirichlet_node

        # For TMS we use only tetrahedra (elm_type=4). The triangles (elm_type=3) are removed.
        if self.method == "TMS":
            node_nr_before = self.mesh.nodes.nr
            nodes_before = np.copy(self.mesh.nodes.node_coord)
            self.mesh = self.mesh.crop_mesh(elm_type=4)
            # ensure that the nodes did not change
            assert self.mesh.nodes.nr == node_nr_before, 'The mesh nodes changed when removing triangles! Please make sure the mesh is valid.'
            assert (self.mesh.nodes.node_coord == nodes_before).all(), 'The mesh nodes changed when removing triangles! Please make sure the mesh is valid.'

        # prepare solver and set self.solver
        self._set_matrices_and_prepare_solver()

        # prepare coil for TMS
        ################################################################################################################
        if method == "TMS":
            if (self.coil is not None and self.fn_coil is not None) or (self.coil is None and self.fn_coil is None):
                raise AssertionError("Provide either filename of TMS coil or Coil object.")

            if coil is None:
                self.coil = TmsCoil.from_file(self.fn_coil)

                if self._logging:
                    logger.info(f'Loaded coil from file: {self.fn_coil}')

        # prepare electrode for TES
        ################################################################################################################
        if method == "TES" and electrode is None:
            raise AssertionError("Please provide TES electrode object for TES simulations.")

    def _set_logger(self,log_fn):
        """
        Set-up loggger to write to a file

        Parameters
        ----------
        log_fn: str, or bool
            if True, the standard simnibs logger is used
            if a filename is given, the imnibs logger is set up to save to this file
        """
        if type(log_fn) == bool:
            self._logging = log_fn
            return
        if type(log_fn) == str:
            fh = logging.FileHandler(log_fn, mode='w')
            formatter = logging.Formatter(
                f'[ %(name)s {__version__} - %(asctime)s - %(process)d ]%(levelname)s: %(message)s')
            fh.setFormatter(formatter)
            fh.setLevel(logging.DEBUG)
            logger.addHandler(fh)
            simnibs_logger.register_excepthook(logger)
            return
        raise ValueError(f"log_fn {log_fn}: has to be string or bool")

    def update_field(self, electrode=None, matsimnibs=None, didt=1e6, fn_electrode_txt=None,
                     dirichlet_correction=False):
        """
        Calculating and updating electric field for given coil position (matsimnibs) for TMS or electrode position (TES)

        Parameters
        ----------
        electrode : list of list of np.ndarray [n_sims][n_electrodes]
            List of list of the surface tags or nodes where the currents to be applied for multiple simulations.
        currents : list of list of np.ndarray [n_sims][n_electrodes]
            List of list of the currents in each surface for multiple simulations.
            If Dirichlet correction has to be performed, the input currents are overwritten by an initial guess
            obtained from previous runs to accelerate convergence.
        matsimnibs : np.array of float [4 x 4 x n_sim]
            Tensor containing the coil positions and orientations in SimNIBS space for multiple simulations.
        didt : float
            Rate of change of coil current (A/s) (e.g. 1 A/us = 1e6 A/s)
        fn_electrode_txt : str, optional, default: None
            Filename of .txt file containing the electrode node coords and the optimized currents
            Note: Will be only used if dirichlet_correction == True
        dirichlet_correction : bool, optional, default: False
            Apply iterative Dirichlet correction to ensure same voltage over the whole electrode when solving.

        Returns
        -------
        e : list of lists [n_sim][n_roi] containing np.arrays of float
            Electric field for queried simulations in ROIs.
        """

        if type(electrode) is not list:
            electrode = [electrode]

        if self.method == "TMS":
            if matsimnibs.ndim < 3:
                matsimnibs = matsimnibs[:, :, np.newaxis]
            n_sim = matsimnibs.shape[2]
        else:
            n_sim = len(electrode)

        self.e = [[0 for _ in range(self.n_roi)] for _ in range(n_sim)]

        # loop over simulation conditions (multiple coil positions or separate electrode configurations)
        for i_sim in range(n_sim):
            start = time.time()

            # determine RHS
            ############################################################################################################
            if self.method == "TMS":
                self.b = self.set_rhs(matsimnibs=matsimnibs[:, :, i_sim])

            elif self.method == "TES":
                self.b = self.set_rhs(electrode=electrode[i_sim])

            # solve for potential
            ############################################################################################################
            if self.method == "TES" and dirichlet_correction:
                self.v = self.solve_dirichlet_correction(electrode=electrode[i_sim],
                                                         fn_electrode_txt=fn_electrode_txt)
            else:
                self.v = self.solve(b=self.b)

            # calculate e-field
            ############################################################################################################
            for i_roi, r in enumerate(self.roi):
                if self.v is None:
                    self.e[i_sim][i_roi] = None
                    if self._logging:
                        logger.warning("Warning! Simulation failed! Returning e-field: None!")
                else:
                    self.e[i_sim][i_roi] = r.calc_fields(v=self.v, dadt=self.dadt, dataType=self.dataType[i_roi])

                if self.method == "TMS":
                    self.e[i_sim][i_roi] *= didt

            stop = time.time()

            if self._logging:
                logger.info(f"Finished simulation #{i_sim + 1}/{n_sim} (time: {(stop-start):.3f}s).")

        return self.e

    def set_rhs(self, electrode=None, matsimnibs=None):
        """
        Set up right hand side (force vector) of equation system.

        Parameters
        ----------
        electrode : ElectrodeArrayPair or CircularArray object instance
            Electrode
        matsimnibs : np.array of float [4 x 4]
            Matrix containing the coil position and orientation in SimNIBS space.

        Returns
        -------
        b : np.array of float [n_nodes - 1]
            Right hand side of equation system (without Dirichlet node)
        """
        # TODO: import here as numba interacts badly with pyqt (GUI). remove
        # here once this is resolved.
        from simnibs.simulation.numba_fem_utils import sumf2, node2elmf, sumf3

        if self.method == "TES":
            # gather electrode currents and associated node indices
            electrodes = []
            currents = []
            for _electrode_array in electrode._electrode_arrays:
                for _ele in _electrode_array.electrodes:
                    if _ele.node_current is not None:
                        currents.append(_ele.node_current)
                        self.fem.weigh_by_area = False
                    else:
                        currents.append(_ele.ele_current)
                        self.fem.weigh_by_area = True
                    electrodes.append(_ele.node_idx + 1)

            b = self.fem.assemble_rhs(electrodes=electrodes,       # list of node indices of electrodes
                                      currents=currents)           # list of electrode currents

        elif self.method == "TMS":
            # determine magnetic vector potential
            if self.useElements:
                self.dadt = self.coil.get_da_dt_at_coordinates(self.coordinates.T, matsimnibs)
            else:
                da_dt = self.coil.get_da_dt_at_coordinates(self.coordinates.T, matsimnibs).T
                #self.dadt = NodeData(da_dt, mesh=self.mesh).node_data2elm_data()[:]

                self.dadt = node2elmf(da_dt, self.reshaped_node_numbersT)

            #reshaped_node_numbers=self.reshaped_node_numbersT,
            #useElements=self.useElements

            # isotropic
            if self.cond.ndim == 1:
                b = sumf2(x=self.force_integrals, y=self.dadt, w=self.reshaped_node_numbers)

            # anisotropic
            elif self.cond.ndim == 3:
                b = sumf3(v=self.volume, dadt=self.dadt, g=self.gradient, nn=(self.node_numbers-1).reshape(-1), c=self.cond)

        else:
            raise NotImplementedError("Simulation method not implemented yet. Method is either 'TMS' or 'TES'.")

        return b

    def normalize_solution(self, v_nodes, channel_id, ele_id, node_area, currents, electrode):
        """
        Normalize electric potential and input currents for later correction step.

        Parameters
        ----------
        v_nodes : np.ndarray of float [n_nodes]
            Electric potential in electrode nodes for all channels where current is applied
        channel_id : np.ndarray of int [n_nodes]
            Channel ID of nodes where current is applied
        ele_id : np.ndarray of int [n_nodes]
            Ele ID of nodes where current is applied
        node_area : np.ndarray of int [n_nodes]
            Area of nodes where current is applied
        currents : list of np.ndarray of float [n_channel][n_ele or n_nodes_channel]
            Applied currents to nodes
        electrode : CircularArray or ElectrodeArrayPair instance
            Electrode

        Returns
        -------
        v_ele_norm : list of np.array of float [n_channel][n_ele or n_nodes_channel]
            Normalized difference of electrode voltages for each channel
        currents_ele_norm : np.array of float [n_ele or n_nodes_channel]
            Normalized difference of electrode currents for each channel
        """
        channel_id_unique = np.unique(channel_id)

        # average and weight potential with node area over electrodes and whole channel (separately)
        v_mean_ele = []         # [n_channel][n_ele]
        v_mean_channel = []     # [n_channel]

        for i, _channel_id in enumerate(channel_id_unique):
            channel_mask = _channel_id == channel_id

            # average area weighted potential over whole channel
            v_mean_channel.append(np.sum(node_area[channel_mask] * v_nodes[channel_mask]) / \
                                  np.sum(node_area[channel_mask]))

            ele_id_unique = np.unique(ele_id[channel_mask])

            if electrode.dirichlet_correction_detailed:
                # do not average voltages over electrode for detailed Dirichlet correction for this channel
                v_mean_ele.append(v_nodes[channel_mask])
            else:
                # determine average potentials over electrodes for this channel
                v_mean_ele.append([])
                for _ele_id in ele_id_unique:
                    ele_mask = _ele_id == ele_id
                    mask = ele_mask * channel_mask
                    v_mean_ele[i].append(np.sum(node_area[mask] * v_nodes[mask]) / np.sum(node_area[mask]))

                v_mean_ele[i] = np.hstack(v_mean_ele[i])

        # print(f"v_mean_ele: {v_mean_ele}")

        v_ele_norm = [0 for _ in range(electrode._n_channel)]
        currents_ele_norm = [0 for _ in range(electrode._n_channel)]

        # calculate rel. difference of electrode potentials to average electrode potential
        for i_channel in range(electrode._n_channel):
            std_scale = np.std(v_mean_channel)

            # calculate differences of electrode potentials to channel potential
            v_ele_norm[i_channel] = (v_mean_ele[i_channel] - np.mean(v_mean_ele[i_channel])) / std_scale
            currents_ele_norm[i_channel] = (currents[i_channel] - np.mean(currents[i_channel])) / np.mean(currents[i_channel])

            if (np.isnan(v_ele_norm[i_channel])).any() or (np.isnan(currents_ele_norm[i_channel])).any():
                te = 1

        return v_ele_norm, currents_ele_norm

    def solve_dirichlet_correction(self, electrode, fn_electrode_txt=None):
        """
        Solve system of equations Ax=b and corrects input currents such that the electrodes with the same channel ID
        have the same potential. Finally, add Dirichlet node (V=0) to solution.

        Parameters
        ----------
        electrode : CircularArray or ElectrodeArrayPair instance
            Electrode
        fn_electrode_txt : str, optional, default: None
            Filename of .txt file containing the electrode node coords and the optimized currents

        Returns
        -------
        v : np.array of float [n_nodes]
            Corrected solution (including the Dirichlet node at the right position)
        """

        electrode.compile_node_arrays()

        n_nodes_total = electrode._node_current.shape[0]

        # update rhs
        b = self.set_rhs(electrode=electrode)

        th_maxrelerr = 0.02
        th_wrong_current_sign = 0.02

        maxiter = 1000

        # create nodal arrays
        # v_nodes | v_node_idx | channel_id_nodes | ele_id_nodes | node_area
        # ------------------------------------------------------------------
        #   ...         4               0                 0           ...
        #   ...         1               0                 0           ...
        #   ...         15              0                 1           ...
        #   ...         3               0                 1           ...
        #   ...         9               1                 0           ...
        #   ...         8               1                 0           ...
        #   ...         22              1                 0           ...
        #   ...         0               1                 1           ...
        #   ...         10              1                 1           ...
        #   ...         11              1                 1           ...

        n_channel = electrode._n_channel

        I_mean = np.zeros(len(electrode._channel_id_unique))

        # read currents from electrode
        I = [[] for _ in range(n_channel)]
        I_sign = [[] for _ in range(n_channel)]
        I_total = [[] for _ in range(n_channel)]

        for i_channel, _channel_id in enumerate(electrode._channel_id_unique):
            # total channel current
            I_total[i_channel] = electrode._current_channel[i_channel]

            if electrode.dirichlet_correction_detailed:
                # take nodal current
                I[i_channel] = electrode._node_current[electrode._node_channel_id == _channel_id]
                I_sign[i_channel] = electrode._node_current_sign[electrode._node_channel_id == _channel_id]

                # mean current over electrodes [n_channel]
                I_mean[i_channel] = electrode._current_channel[i_channel] / np.sum(electrode._node_channel_id == _channel_id)

            else:
                # take total electrode current (sum nodal current up)
                for _ele_id in np.unique(electrode._node_ele_id[electrode._node_channel_id == _channel_id]):
                    mask = (electrode._node_channel_id == _channel_id) * (electrode._node_ele_id == _ele_id)
                    I[i_channel].append(np.sum(electrode._node_current[mask]))
                    I_sign[i_channel].append(electrode._node_current_sign[mask][0])
                I[i_channel] = np.hstack(I[i_channel])
                I_sign[i_channel] = np.hstack(I_sign[i_channel])

                # mean current over electrodes [n_channel]
                I_mean[i_channel] = electrode._current_mean[i_channel]

        # Solve iteratively until maxrelerr is reached
        # Iteration: 0

        # solve
        v = self.solve(b)
        v_nodes = v[electrode._node_idx]

        # v_norm: [n_channel][n_ele], I_norm: [n_channel][n_ele]
        v_norm, I_norm = self.normalize_solution(v_nodes=v_nodes,
                                                 channel_id=electrode._node_channel_id,
                                                 ele_id=electrode._node_ele_id,
                                                 node_area=electrode._node_area,
                                                 currents=I,
                                                 electrode=electrode)

        for i_channel in range(n_channel):
            v_norm[i_channel] = v_norm[i_channel][np.newaxis, :]
            I_norm[i_channel] = I_norm[i_channel][np.newaxis, :]

        j = 0
        maxrelerr = np.max([np.max(np.abs(v_norm[k][-1])) for k in range(n_channel)])
        denom_factor = 100
        n_wrong_current_signs = [0]

        while maxrelerr > th_maxrelerr or n_wrong_current_signs[-1] / n_nodes_total > th_wrong_current_sign:
            if self._logging:
                if j == 0:
                    logger.debug(
                        f"iter: {j:03}, err: {maxrelerr:.3f}, df={denom_factor:.3f}")
                else:
                    logger.debug(f"iter: {j:03}, err: {maxrelerr:.3f}, wcs: {n_wrong_current_signs[-1] / n_nodes_total:.3f}, df={denom_factor:.3f}")

            j += 1
            if j == maxiter:
                break

            if j > maxiter:
                if self._logging:
                    logger.warning('warning: did not converge after ' + str(maxiter) + ' iterations')

                # reset to original currents
                for _electrode_array in electrode._electrode_arrays:
                    for _ele in _electrode_array.electrodes:
                        _ele.ele_current = _ele.ele_current_init

                electrode.compile_node_arrays()

                # return None

            if j == 1:
                # It 1: make small step hopefully in the right direction
                I = []
                all_pos = False
                step_size_factor = -0.1

                while not all_pos:
                    for i_channel in range(n_channel):
                        if len(I_norm[i_channel][-1, :]) == 1:
                            I.append(np.array([I_mean[i_channel]]))
                            all_pos = True
                        else:
                            I.append(step_size_factor * np.sign(I_mean[i_channel]) * v_norm[i_channel][0, :] / np.max(np.abs(v_norm[i_channel][0, :])))

                            if ((I[i_channel] - np.mean(I[i_channel]) + 1) < 0).any():
                                step_size_factor /= 10
                                if self._logging:
                                    logger.debug(f"Decreasing step size factor to {step_size_factor}")
                                all_pos = False
                                break

                            all_pos = True

            else:
                # It 2 ... use gradient descent
                n_wrong_current_signs.append(0)
                for i_channel in range(n_channel):
                    if len(I_norm[i_channel][-1, :]) == 1:
                        I[i_channel] = I[i_channel]
                    else:
                        denom = (v_norm[i_channel][-1, :] - v_norm[i_channel][-2, :])
                        denom += np.sign(denom) * maxrelerr / denom_factor  # regularize: avoid too small denominators to gain stability
                        denom[denom == 0] = 1e-12
                        I[i_channel] = I_norm[i_channel][-1] - (I_norm[i_channel][-1, :] - I_norm[i_channel][-2, :]) / denom * v_norm[i_channel][-1, :]
                        n_wrong_current_signs[-1] += np.sum(((I[i_channel] - np.mean(I[i_channel]) + 1) < 0))
                        # print(f"Channel: {i_channel}: {n_wrong_current_signs[i_channel]}/{len(I[i_channel])} wrong current signs.")

                if (np.array(n_wrong_current_signs[-3:]) == n_wrong_current_signs[-1]).all():
                    denom_factor = np.random.rand(1)[0] * 20
                elif n_wrong_current_signs[-1] == 0:
                    denom_factor = 15
                else:
                    denom_factor = 100 * n_wrong_current_signs[-1] / n_nodes_total

            # convert back from I_norm to I
            I = [I_mean[i_channel] * (I[i_channel] - np.mean(I[i_channel]) + 1) for i_channel in range(n_channel)]

            # print(f"sum current: sum(I[0]) = {np.sum(I[0])}   sum(I[1]) = {np.sum(I[1])}")

            # write currents in electrodes
            for i_channel, _channel_id in enumerate(electrode._channel_id_unique):
                if electrode.dirichlet_correction_detailed:
                    # write nodal current
                    electrode._node_current[electrode._node_channel_id == _channel_id] = I[i_channel]
                else:
                    # calculate nodal current from total electrode current
                    for i_ele, _ele_id in enumerate(np.unique(electrode._node_ele_id[electrode._node_channel_id == _channel_id])):
                        mask = (electrode._node_channel_id == _channel_id) * (electrode._node_ele_id == _ele_id)
                        electrode._node_current[mask] = I[i_channel][i_ele] * electrode._node_area[mask] / np.sum(electrode._node_area[mask])

            # update electrodes from node arrays
            electrode.update_electrode_from_node_arrays()

            # update rhs
            b = self.set_rhs(electrode=electrode)

            # solve
            v = self.solve(b)
            v_nodes = v[electrode._node_idx]

            # normalize
            _v_norm, _I_norm = self.normalize_solution(v_nodes=v_nodes,
                                                       channel_id=electrode._node_channel_id,
                                                       ele_id=electrode._node_ele_id,
                                                       node_area=electrode._node_area,
                                                       currents=I,
                                                       electrode=electrode)

            # append
            for i_channel in range(n_channel):
                v_norm[i_channel] = np.vstack((v_norm[i_channel], _v_norm[i_channel][np.newaxis, :]))
                I_norm[i_channel] = np.vstack((I_norm[i_channel], _I_norm[i_channel][np.newaxis, :]))

            # update error
            maxrelerr = np.max([np.max(np.abs(v_norm[i_channel][-1])) for i_channel in range(n_channel)])

        # final number of iterations to determine optimal currents
        self.n_iter_dirichlet_correction = j

        masks_sign = [(np.sign(I[i_channel]) == I_sign[i_channel]) for i_channel in range(n_channel)]

        # test if signs are correct return no solution
        if not np.array([(masks_sign[i_channel]).all() for i_channel in range(n_channel)]).all():
            masks_sign_not = copy.copy(masks_sign)
            masks_sign_not[0] = np.logical_not(masks_sign[0])
            masks_sign_not[1] = np.logical_not(masks_sign[1])

            # ensure correct sign
            I = [np.abs(I[i_channel]) * I_sign[i_channel] for i_channel in range(n_channel)]

            # normalize currents again to match maximal channel current (we flipped currents before)
            I = [I[i_channel] / np.sum(np.abs(I[i_channel])) * np.abs(I_total[i_channel]) for i_channel in range(n_channel)]

            # write final currents in electrode
            for i_channel, _channel_id in enumerate(electrode._channel_id_unique):
                if electrode.dirichlet_correction_detailed:
                    # write nodal current
                    electrode._node_current[electrode._node_channel_id == _channel_id] = I[i_channel]
                else:
                    # calculate nodal current from total electrode current
                    for i_ele, _ele_id in enumerate(
                            np.unique(electrode._node_ele_id[electrode._node_channel_id == _channel_id])):
                        mask = (electrode._node_channel_id == _channel_id) * (electrode._node_ele_id == _ele_id)
                        electrode._node_current[mask] = I[i_channel][i_ele] * electrode._node_area[mask] / np.sum(
                            electrode._node_area[mask])

            # update electrodes from node arrays
            electrode.update_electrode_from_node_arrays()

            # apply current outlier correction
            if electrode.current_outlier_correction:
                electrode.apply_current_outlier_correction()

            # update rhs
            b = self.set_rhs(electrode=electrode)

            # solve
            v = self.solve(b)
            v_nodes = v[electrode._node_idx]

            # normalize
            _v_norm, _ = self.normalize_solution(v_nodes=v_nodes,
                                                 channel_id=electrode._node_channel_id,
                                                 ele_id=electrode._node_ele_id,
                                                 node_area=electrode._node_area,
                                                 currents=I,
                                                 electrode=electrode)

            # final error
            maxrelerr = np.max([np.max(np.abs(_v_norm[i_channel])) for i_channel in range(n_channel)])

            if self._logging:
                logger.debug(f"Correcting current signs (pos: {np.sum(masks_sign_not[0])}/{len(masks_sign_not[0])}, "
                      f"neg: {np.sum(masks_sign_not[1])}/{len(masks_sign_not[1])}), final maxrelerr: {maxrelerr:3f}")

        # add optimal currents to CurrentEstimator (training data) to improve estimation in further iterations
        # this is done electrode wise and not node wise because the number of nodes changes depending on the
        # electrode position
        if electrode._current_estimator is not None:
            electrode_pos = [_electrode_array.electrode_pos for _electrode_array in electrode._electrode_arrays]

            I_ele = []
            for i_channel, _channel_id in enumerate(electrode._channel_id_unique):
                I_ele.append([])
                for i_ele, _ele_id in enumerate(np.unique(electrode._node_ele_id[electrode._node_channel_id == _channel_id])):
                    mask = (electrode._node_channel_id == _channel_id) * (electrode._node_ele_id == _ele_id)
                    I_ele.append(np.sum(electrode._node_current[mask] * electrode._node_area[mask] / np.sum(electrode._node_area[mask])))

            electrode._current_estimator.add_training_data(electrode_pos=np.hstack(electrode_pos),
                                                           current=np.hstack(I_ele))

        if fn_electrode_txt is not None:
            np.savetxt(fn_electrode_txt, np.hstack((electrode._node_coords, electrode._node_current[:, np.newaxis])))

        return np.squeeze(v)

    def solve(self, b):
        """
        Solve system of equations Ax=b and add Dirichlet node (V=0) to solution.

        Parameters
        ----------
        b : np.array of float [n_nodes]
            Right hand side of equation system (with Dirichlet node)

        Returns
        -------
        v : np.array of float [n_nodes]
            Solution (including the Dirichlet node at the right position)
        """

        # remove dirichlet node
        b_reduced = np.delete(b, self.dirichlet_node-1)

        # solve equation system
        x = self.solver.solve(b_reduced)

        # add Dirichlet node to solution
        v = np.insert(x, self.dirichlet_node-1, 0)

        return np.squeeze(v)

    def _set_matrices_and_prepare_solver(self):
        """
        Set matrices and initialize the pardiso solver for update_position in self.solver.
        """

        if self.cond is None:
            # prepare FEM
            self.simulist = SimuList(mesh=self.mesh)

            # prepare conductivity (scalar or anisotropic)
            self.simulist.anisotropy_type = self.anisotropy_type
            self.simulist.fn_tensor_nifti = self.fn_tensor_nifti
            self.cond = self.simulist.cond2elmdata()

        if self.method == "TMS":
            self.cond = self.cond.value.squeeze()

            if self.cond.ndim == 2:
                self.cond = self.cond.reshape(-1, 3, 3)

            self.node_numbers = self.mesh.elm.node_number_list  # self.node_numbers
            useElements = self.useElements
            node_coordinates = self.mesh.nodes.node_coord.T

            # get the coordinates of nodal points/element centers to prepare for the calculation of dadt.
            # Use self.coordinates to interpolate the field in calculate_dadt().
            self.coordinates = get_coordinates(node_coordinates, self.node_numbers, useElements)

            # set the reshaped_node_numbers (there are two different ways to reshape it and sometimes
            # it is more efficient to use one or the other
            self.reshaped_node_numbers = (self.node_numbers - 1).T.ravel()
            self.reshaped_node_numbersT = (self.node_numbers - 1).ravel()

            # get the distances between local node [0] and nodes [1], [2] and [3] in each element
            local_dist, _ = _get_local_distances(self.node_numbers, node_coordinates)

            # calculate the volume of a tetrahedron given the coordinates of its four nodal points
            self.volume = np.abs(np.linalg.det(local_dist)) / 6.

            # get gradient operator
            self.gradient = _get_gradient(local_dist)

            # set the force integrals. We use the force integrals to assemble the right hand side force vector
            if self.cond.ndim == 1:
                self.force_integrals = get_force_integrals(self.volume, self.gradient, self.cond)

            self.fem = TMSFEM(mesh=self.mesh,
                              cond=self.cond,
                              solver_options=self.solver_options,
                              solver_loglevel=self.solver_loglevel,
                              dirichlet_node=self.dirichlet_node
                              )

            self.fem.prepare_solver()
            self.solver = self.fem._solver

        elif self.method == "TES":
            self.fem = TDCSFEMNeumann(mesh=self.mesh,
                                      cond=self.cond,
                                      ground_electrode=self.dirichlet_node,
                                      input_type="nodes", #this option is actually not supported fem.py claims only 'tag' is supported
                                      solver_options=self.solver_options,
                                      solver_loglevel=self.solver_loglevel
                                      )

            self.fem.prepare_solver()
            self.solver = self.fem._solver


def assemble_force_vector(force_integrals, reshaped_node_numbers, dadt):
    """
    Assembly of the force vector in a system of linear equations stiffmat * x = forcevec. for TMS.

    Parameters
    ----------
    force_integrals : np.array of float [4, number_of_elements, 3]
        Force_integrals for rhs calculation derived by volume * conductivity * gradient
    reshaped_node_numbers : np.array of int [4 * n_elements + 1]
        Flattened node number list (connectivity matrix)
        (node_numbers - 1).T.reshape(-1)
    dadt: NodeData or ElementData
        dA/dt field at each node or element

    Returns
    -------
    forcevec: np.array [n_nodes - 1]
        Right-hand side (without the Dirichlet node)
    """
    # TODO: import here as numba interacts badly with pyqt (GUI). remove once
    # this is resolved.
    from simnibs.simulation.numba_fem_utils import sumf

    # integrate in each node of each element, the value for repeated nodes will be summed
    # together later
    node_integrals = np.zeros(force_integrals.shape[:2], order='C')

    # force_integrals: (4, number_of_elements, 3); dadt: (number_of_elements, 3); node_integrals: (4, number_of_elements)
    sumf(force_integrals, dadt, node_integrals)

    # Assembles the right hand side for TMS.
    # forcevec: (number_of_nodes,), reshaped_node_numbers: (number_of_elements*4,), node_integrals: (4, number_of_elements)

    # keep np.bincount to make the testing more stable. May change it back to bincount_nb() for performance reason.
    # forcevec = bincount_nb(reshaped_node_numbers, node_integrals.reshape(-1)).reshape(-1, 1)
    forcevec = np.bincount(reshaped_node_numbers, node_integrals.reshape(-1)).reshape(-1, 1)

    return forcevec

def calculate_element_centers(node_coordinates, node_numbers):
    """
    Calculate the center of each element. The center is the average of the vertices:
    Centroid = (a + b + c + d) / 4.

    Parameters
    ----------
    node_coordinates : np.array of size [3 x n_nodes]
        Coordinates of the nodes (x, y, z)
    node_numbers : np.array of size [n_elements x 4]
        Node number list (connectivity list)

    Returns
    -------
    element_centers : np.array of float [3 x n_elements]
        Coordinates (x, y, z) of the element centers
    """
    # node_coordinates: (number_of_nodes, 3)
    # node_numbers: (number_of_elements, 4)
    index = node_numbers - 1

    # calculate a/4, b/4, c/4 and d/4 first
    coordinates_normalized = node_coordinates * (1./index.shape[1])

    # locate the coordinates of nodes in each element
    # (3, number_of_elements, 4)
    pos = coordinates_normalized[:, index]

    # then calculate the sum of a/4, b/4, c/4 and d/4
    # (3, number_of_elements)
    element_centers = np.einsum('ijk->ij', pos)

    return element_centers


def get_force_integrals(volume, gradient, conductivity):
    """
    Assembly of the integration part in the force vector without dadt for TMS and conductivity.ndim == 1
    in physical space (unit: cubic metre)

    Parameters
    ----------
    volume : np.array of float [n_elements]
        Volume of the tetrahedra in (mm³)
    gradient : np.array of size [n_elements, 4, 3]
        Gradient in each tetrahedra
    conductivity : np.array of float [n_elements]
        Electrical conductivity value (in S/m) assigned to each tetrahedra

    Returns
    -------
    force_integrals : np.array of float [4, number_of_elements, 3]
        Force_integrals for rhs calculation derived by volume * conductivity * gradient
    """

    if conductivity.ndim != 1:
        raise ValueError('Invalid conductivity array')

    # volume and conductivity are 1D array shape=(M,)
    vol_cond = volume * conductivity * 1e-6  # 1e-6 is to convert 'mm' to 'm'

    # calculate the force_integrals = volume * conductivity * gradient, and rearrange the dimension to 4xMx3
    force_integrals = np.swapaxes(-vol_cond[:, None, None] * gradient, 0, 1)  # (4, number_of_elements, 3)

    return force_integrals


def get_coordinates(node_coordinates, node_numbers, useElements):
    """
    Calculate the centers of the elements.
    If useElements is true, get the coordinates of element centers; otherwise get the coordinates of nodes

    Parameters
    ----------
    node_coordinates : np.array of size [3 x n_nodes]
        Coordinates of the nodes (x, y, z)
    node_numbers : np.array of size [n_elements x 4]
        Node number list (connectivity list)
    useElements : bool
        Get coordinates in the element centers (True) or of the nodes (False)

    Returns
    -------
    coordinates : np.array of float of size [n_element] or [n_nodes]
        Centers of the elements or node coordinates.
    """

    # useElements = True: get the coordinates of the element centers
    if useElements:
        # calculate the coordinates of the element centers
        # Notice coordinates: input and output are the same
        coordinates = calculate_element_centers(node_coordinates, node_numbers)

    else:
        # get the coordinates of the nodes
        coordinates = node_coordinates

    return coordinates


def delete_row_csr(mat, i):
    """
    Delete rows from a sparse matrix in Compressed Sparse Row (CSR) format.

    Parameters
    ----------
    mat : scipy.sparse matrix in csr format
        Sparse matrix (CSR)
    i : int or list of int
        Indices of rows to delete from mat.

    Returns
    -------
    mat : scipy.sparse matrix in csr format
        Sparse matrix with deleted rows.
    """
    if not isinstance(mat, sparse.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    n = mat.indptr[i+1] - mat.indptr[i]
    if n > 0:
        mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
        mat.data = mat.data[:-n]
        mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
        mat.indices = mat.indices[:-n]
    mat.indptr[i:-1] = mat.indptr[i+1:]
    mat.indptr[i:] -= n
    mat.indptr = mat.indptr[:-1]
    mat._shape = (mat._shape[0]-1, mat._shape[1])

    return mat


def delete_col_csr(mat, i):
    """
    Delete columns from a sparse matrix in Compressed Sparse Row (CSR) format.

    Parameters
    ----------
    mat : scipy.sparse matrix in csr format
        Sparse matrix (CSR)
    i : int or list of int
        Indices of columns to delete from mat.

    Returns
    -------
    mat : scipy.sparse matrix in csr format
        Sparse matrix with deleted columns.
    """
    if not isinstance(mat, sparse.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")

    if type(i) is int:
        i = [i]

    mask = np.ones(mat.shape[1], dtype=bool)
    mask[i] = False

    return mat[:, mask]


class FemTargetPointCloud:
    """
    Class containing methods to compute the electric field at point cloud points from the electric potential at nodes (fast).

    !!! NOTE !!!
    TMS:
        The mesh used in the OnlineFem and in the FemTargetPointCloud have to have the same node coordinates and the same node ordering.
        They also have to have the same elements and the same element ordering.
    TES:
        The mesh used in the OnlineFem and in the FemTargetPointCloud have to have the same node coordinates and the same node ordering.

    Parameters
    ----------
    mesh : Msh object instance
        Head mesh
    center : np.ndarray of float [n_roi_center x 3]
        The point coordinates where the e-field is calculated, e.g. element center of tetrahedra.
    gradient : np.array of float [n_tet_mesh_required x 4 x 3]
        Pre-calculated gradient operator of the tetrahedral edges.
    nearest_neighbor : bool
        Whether to use SPR interpolation or nearest neighbor
    fill_nearest : boolean
        Whether to fill values outside the mesh with the nearest value or assign zeros to these points.
        If set to True then use nearest neighbor assigns the nearest value; otherwise assign to 0

    Attributes
    ----------
    center : np.ndarray of float [n_roi_center x 3]
        The point coordinates where the e-field is calculated, e.g. element center of triangles or tetrahedra.
    n_center : int
        Number of center points in the ROI.
    gradient : np.array of float [n_tet_mesh_required x 4 x 3]
        Gradient operator of the tetrahedral edges.
    node_index_list :  [n_tet_mesh_required x 4]
        Connectivity list of the head model (only using the required tetrahedra) (0 indexed)
    sF : sparse matrix of float [n_points_ROI x n_tet_mesh_required ]
        Sparse matrix for SPR interpolation.
    inside : np.array of bool [n_points]
        Indicator if points are lying inside the model.
    idx : np.array of int [n_tet_mesh_required]
        Indices of tetrahedra, in whcih the E-field is calculated
    n_tet_mesh : int
        Number of tetrahedra in the whole head mesh
    """
    def __init__(self, mesh, center=None, gradient=None, nearest_neighbor=False, fill_nearest=False):

        self.center = center
        self.sF = None
        self.n_center = None
        self.gradient = gradient
        self.fill_nearest = fill_nearest

        if type(self.center) is list:
            self.center = np.array(self.center)

        # crop mesh that only tetrahedra are included
        mesh_cropped: Msh = mesh.crop_mesh(elm_type=4)
        # ensure that the nodes did not change
        assert mesh_cropped.nodes.nr == mesh.nodes.nr, 'The mesh nodes changed when removing triangles! Please make sure the mesh is valid.'
        assert (mesh_cropped.nodes.node_coord == mesh.nodes.node_coord).all(), 'The mesh nodes changed when removing triangles! Please make sure the mesh is valid.'

        self.n_tet_mesh = mesh_cropped.elm.node_number_list.shape[0]

        # compute gradient
        if not gradient:
            # get the lengths of the tetrahedral edges _get_local_distances(node_numbers, node_coordinates)
            local_dist, _ = _get_local_distances(node_numbers=mesh_cropped.elm.node_number_list,
                                                 node_coordinates=mesh_cropped.nodes.node_coord.T)
            # get gradient operator
            gradient = _get_gradient(local_dist)

        if nearest_neighbor:
            th_with_points, bar = mesh_cropped.find_tetrahedron_with_points(center, compute_baricentric=True)
            self.inside = th_with_points != -1
            self.any_outside = np.any(~self.inside)
            self.idx = th_with_points - 1
            if self.any_outside:
                if fill_nearest == 1:  # fill == 'nearest'
                    _, nearest = mesh_cropped.find_closest_element(center[~self.inside], return_index=True)
                    self.idx[~self.inside] = nearest - 1

            self.gradient = gradient[self.idx]
            self.node_index_list = mesh_cropped.elm.node_number_list[self.idx] - 1
            self.n_center = self.center.shape[0]
        else:
            # determine sF matrix for fast interpolation
            # compute sF matrix
            self._get_sF_matrix(mesh_cropped, self.center, fill_nearest)
            self.gradient = gradient[self.idx]
            self.node_index_list = mesh_cropped.elm.node_number_list[self.idx] - 1
            self.n_center = self.center.shape[0]


    def calc_fields(self, v, dadt=None, dataType=0):
        """
        Calculate electric field on target points from v (and A)

        !!! NOTE !!!
        TMS:
            OnlineFem - Tetrahedra mesh
            TargetPointCloud - Tetrahedra mesh
            E calculated from:
                1. v on the mesh nodes
                2. da_dt on the mesh tetrahedra
        TES:
            OnlineFem - Tetrahedra and Triangle (surface) mesh
            TargetPointCloud - Tetrahedra mesh
            E calculated from:
                1. v on the mesh nodes

        Parameters
        ----------
        v : np.ndarray of float [n_nodes_total]
            Electric potential in each node in the whole head model
        dadt : np.ndarray of float, optional, default: None
            Magnetic vector potential in each element in the whole head model (for TMS)
        dataType : int, optional, default: 0
            Return magnitude of electric field (dataType = 0) otherwise return x, y, z components

        Returns
        -------
        e : np.ndarray of float [center]
            Electric field in the target positions
        """
        # TODO: import here as numba interacts badly with pyqt (GUI). remove
        # here once this is resolved.
        from simnibs.simulation.numba_fem_utils import postp, postp_mag, spmatmul

        # get the E field in all tetrahedra (v: (number_of_nodes, 1))
        ################################################################################################################
        if dadt is None:
            # TES
            ############################################################################################################
            if dataType == 0:
                fields = postp_mag(self.gradient, v, np.zeros((self.n_tet_mesh, 3)), self.node_index_list, self.idx)
            else:
                fields = postp(self.gradient, v, np.zeros((self.n_tet_mesh, 3)), self.node_index_list, self.idx)

            # fields = np.einsum('ijk,ij->ik', self.gradient, - (v * 1e3)[self.node_index_list])
            # if dataType == 0:
            #     fields = np.linalg.norm(fields, axis=1)

        else:
            # TMS (dadt should be in elements)
            ############################################################################################################
            if dataType == 0:
                fields = postp_mag(self.gradient, v, dadt, self.node_index_list, self.idx)
            else:
                fields = postp(self.gradient, v, dadt, self.node_index_list, self.idx)

            # fields = np.einsum('ijk,ij->ik', self.gradient, - (v * 1e3)[self.node_index_list]) - dadt[self.idx]

        # Calculate field in ROI
        ############################################################################################################
        if self.sF is not None:
            # interpolate to ROI using sF matrix
            # e = self.sF @ fields
            e = np.zeros((self.n_center, fields.shape[1]))
            spmatmul(self.sF.data, self.sF.indptr, self.sF.indices, fields, e)
        else:
            e = fields
            if not self.fill_nearest and self.any_outside:
                e[~self.inside] = 0

        return e

    def _get_sF_matrix(self, msh, center, fill_nearest, tags=None):
        """
        Create a sparse matrix for SPR interpolation from element data to arbitrary positions (here: the surface nodes)

        Sets the following object variables:
           sF       sparse.csr_matrix: (number_of_center_points, number_of_kept_tetrahedra)
           idx      index of the tetrahedra included in sF as columns
           inside   indices of the positions inside the mesh

        Parameters
        ----------
        msh : Msh object
            Loaded mesh.
        center : np.array of float [n_center_ROI x 3]
            The coordinates of the points that we want to interpolate.
        fill_nearest : boolean
            Whether to fill values outside the mesh with the nearest value or assign zeros to these points.
            If set to True then use nearest neighbor assigns the nearest value; otherwise assign to 0
        tags : list or None, optional, default: None
            The tissue type defines the selected volume in the loaded mesh. Defaults to None.
        """

        assert fill_nearest in [False, True]

        # Set the volume to be GM. The interpolation will use only the tetrahedra in the volume.
        if tags is None:
            th_indices = msh.elm.elm_number
        else:
            th_indices = msh.elm.elm_number[np.isin(msh.elm.tag1, tags)]

        th_with_points, bar = msh.find_tetrahedron_with_points(center, compute_baricentric=True)
        inside = np.isin(th_with_points, th_indices)
        self.inside = inside

        # if any points are inside
        if np.any(inside):
            # get the 'elm_number' of the tetrahedra in 'msh' which contain 'points' in 'points' order
            # assert where_inside.shape == th.shape
            th = th_with_points[inside]

            # interpolate the E field to the points using SPR
            sF = self._get_sF_inside_tissues(msh, th, np.where(inside)[0], bar, center.shape[0])

        # Finally, fill in the unassigned values
        if np.any(~inside):
            if fill_nearest:  # fill == 'nearest'
                if tags is not None:
                    is_in = np.isin(msh.elm.elm_number, th_indices)
                    elm_in_volume = msh.elm.elm_number[is_in]
                    m_in_volume = msh.crop_mesh(elements=elm_in_volume)

                    _, nearest = m_in_volume.find_closest_element(center[~inside], return_index=True)

                    sF[np.where(~inside)[0], elm_in_volume[nearest - 1] - 1] = 1
                else:
                    _, nearest = msh.find_closest_element(center[~inside], return_index=True)

                    sF[np.where(~inside)[0], nearest - 1] = 1

        # convert to csc matrix for fast column indexing
        sF = sF.tocsc()
        self.idx = np.nonzero(sF.indptr[:-1] != sF.indptr[1:])[0]

        # delete the all-zero columns in sF
        sF = sF[:, self.idx]

        # convert to csr matrix for fast row indexing in the matrix multiplication
        self.sF = sparse.csr_matrix(sF)

    def _get_sF_inside_tissues(self, msh, th, w, bar, n_center):
        """
        Create a sparse matrix to interpolate from element data to arbitrary positions using the
        superconvergent patch recovery (SPR) approach.

        Parameters
        ----------
        msh : Msh object
            Loaded mesh
        th : np.array of int [n_center_ROI]
            Indices of the elements in the global mesh (start from 0) that contains the ROI center points.
        w : np.array of int [n_center_ROI_in]
            Indices of the center points that are inside the mesh
        bar : np.array of float [n_points_ROI x 4]
            Barycentric coordinates of the ROI center points of the tetrahedra nodes.
        n_center : int
            Number of ROI center points

        Returns
        -------
        sF : sparse matrix of float [n_center_ROI x n_tet_mesh]
            Sparse matrix for interpolation
        """

        # initialize the sparse matrix
        sF = sparse.dok_matrix((n_center, msh.elm.nr))

        # get the 'tag1' from 'msh' for every element in 'th' in 'center' order
        tag1_inside = msh.elm.tag1[th - 1]

        for t in np.unique(tag1_inside):
            # find the elements in 'tag1_inside' which equals to 't'
            is_t = tag1_inside == t

            if np.any(is_t):
                logger.debug('points inside volume with tag1 = {}: {}'.format(t, np.sum(is_t)))

                # 'm_tag' contains only the tetrahedra with 'elm_number == th_with_t'
                m_tag = msh.crop_mesh(tags=t)

                # 'm_with_t' is sorted because 'elm_number' is always sorted
                m_with_t = msh.elm.elm_number[msh.elm.tag1 == t]

                # the 'elm_number' of elements in 'msh'. These elements contain points and 'tag1 == t'
                th_with_t = th[is_t]

                # get the indices of elements in 'm_tag'. The 'elm_number' of the same elements are 'th_with_t' in 'msh'. 'idx' starts from 0, not 1.
                idx = np.searchsorted(m_with_t, th_with_t)

                # convert the 'elmdata' to 'nodedata'
                sM = self._elm2point_SPR(m_tag, bar[w[is_t]], idx)

                i, j, v = sparse.find(sM)
                sF[(w[is_t])[i], m_with_t[j] - 1] = v

        return sF

    def _elm2point_SPR(self, msh, bar, idx):
        """
        Create the sparse matrix to interpolate from element data to arbitrary positions
        using superconvergent patch recovery

        Zienkiewicz, Olgierd Cecil, and Jian Zhong Zhu. "The superconvergent patch recovery and a posteriori error
        estimates. Part 1: The recovery technique." International Journal for Numerical Methods in Engineering
        33.7 (1992): 1331-1364.

        Parameters
        ----------
        msh : Msh object
            Loaded mesh
        bar : np.array of float [n_center_ROI x 4]
            Barycentric coordinates of the ROI center points of the tetrahedra nodes.
        idx : np.array of int [n_center_ROI]
            Indices of the elements in mesh (start from 0) of ROI center points that need to calculate SPR.

        Returns
        -------
        sF : sparse matrix (CSR) of float [bar.shape[0], n_tetrahedra_mesh]
            Sparse matrix for interpolation
        """

        if len(msh.elm.tetrahedra) == 0:
            raise ValueError("Can only transform volume data")

        # Get the point in the outside surface
        points_outside = np.unique(msh.elm.get_outside_faces())
        outside_points_mask = np.isin(msh.elm[msh.elm.tetrahedra],
                                      points_outside).reshape(-1, 4)

        th_indices = msh.elm.tetrahedra
        th_nodes = msh.elm[th_indices] - 1
        masked_th_nodes = np.copy(th_nodes)
        masked_th_nodes[outside_points_mask] = -1
        masked_th_nodes += 1

        # Calculates the quantities needed for the superconvergent patch recovery
        uq_in, th_nodes = np.unique(masked_th_nodes, return_inverse=True)
        baricenters = msh.elements_baricenters()[th_indices]

        volumes = msh.elements_volumes_and_areas()[th_indices]
        baricenters = np.hstack([np.ones((baricenters.shape[0], 1)), baricenters])

        A = np.empty((msh.nodes.nr + 1, 4, 4))
        for i in range(4):
            for j in range(i, 4):
                A[:, i, j] = np.bincount(masked_th_nodes.reshape(-1),
                                         np.repeat(baricenters[:, i], 4) *
                                         np.repeat(baricenters[:, j], 4),
                                         minlength=msh.nodes.nr + 1)

        # This here only ensures we can invert
        outside = np.isclose(A[:, 0, 0], 0)
        for i in range(4):
            A[outside, i, i] = 1

        A[:, 1, 0] = A[:, 0, 1]
        A[:, 2, 0] = A[:, 0, 2]
        A[:, 3, 0] = A[:, 0, 3]
        A[:, 2, 1] = A[:, 1, 2]
        A[:, 3, 1] = A[:, 1, 3]
        A[:, 3, 2] = A[:, 2, 3]

        Ainv = np.linalg.inv(A)

        node_pos = np.hstack(
            [np.ones((msh.nodes.nr, 1)), msh.nodes.node_coord])
        # Added a dummy to the first position
        node_pos = np.vstack([np.ones((1, 4)), node_pos])

        M = sparse.csr_matrix((msh.nodes.nr + 1, msh.elm.nr))
        for i in range(4):
            M += sparse.csr_matrix(
                (np.einsum(
                    'bi, bij, bj -> b',
                    node_pos[masked_th_nodes[:, i]],
                    Ainv[masked_th_nodes[:, i]],
                    baricenters),
                 (masked_th_nodes[:, i], th_indices - 1)),
                shape=M.shape)

        # Assigns the average value to the points in the outside surface
        th_nodes = msh.elm[th_indices] - 1
        masked_th_nodes = np.copy(th_nodes)
        masked_th_nodes[~outside_points_mask] = -1
        masked_th_nodes += 1

        for i in range(4):
            M += sparse.csr_matrix(
                (volumes, (masked_th_nodes[:, i], th_indices - 1)),
                shape=M.shape)

        M = M[1:]

        node_vols = np.bincount(
            th_nodes.reshape(-1),
            np.repeat(volumes, 4),
            minlength=msh.nodes.nr + 1)

        normalization = np.ones(msh.nodes.nr)
        normalization[points_outside - 1] = 1 / (node_vols[points_outside - 1] + np.finfo(float).eps)

        D = sparse.dia_matrix(
            (normalization, 0), shape=(msh.nodes.nr, msh.nodes.nr))
        M = D.dot(M)

        # interpolate from the nodal values to the point values using linear interpolation
        P = bar.shape[0]
        sS = sparse.dok_matrix((P, msh.nodes.nr))
        node_idx = msh.elm.node_number_list[idx] - 1
        sS[np.arange(P)[:, None], node_idx] = bar
        sS = sparse.csr_matrix(sS)
        res = sS @ M

        return res

def _get_local_distances(node_numbers, node_coordinates):
    """
    Get the distances between local node [0] and nodes [1], [2] and [3] in each element

    Parameters
    ----------
    node_numbers : np.array of size [number_of_elements, 4]
        Node number list (connectivity list of tetrahedra)
    node_coordinates : np.array of size [3 x n_nodes]
        Coordinates of the nodes (x, y, z)

    Returns
    -------
    local_dist : np.array of size [number_of_elements, 3, 3]
        Distances between local node [0] and nodes [1], [2] and [3] in each element
    local_node_coords : np.array of size [number_of_elements, 3, 3]
        Node coordinates
    """

    # get the local node coordinates in each element
    # node_coordinates: (3, number_of_nodes)
    # node_numbers: (number_of_elements, 4)
    # local_node_coords: (3, number_of_elements, 4)
    local_node_coords = node_coordinates[:, node_numbers - 1]

    # get the distances between local node [0] and nodes [1], [2] and [3] in each element
    # local_dist: (3, number_of_elements, 3)
    local_dist = local_node_coords[:, :, 1:] - local_node_coords[:, :, [0]]

    # out: (number_of_elements, 3, 3)
    return np.moveaxis(local_dist, 0, -1), (local_node_coords[:, :, 0].T)[:, :, None]


def _get_gradient(local_dist):
    """
    Calculate the gradient of a function in each tetrahedra.

    local_dist * gradient = project_matrix

    project_matrix is a projection matrix
    project_matrix = [-1, 1, 0, 0]
                     [-1, 0, 1, 0]
                     [-1, 0, 0, 1]
    and local_dist is the distances between local node [0] and nodes [1], [2] and [3] in each element.

    Parameters
    ----------
    local_dist : np.array of size [n_elements, 3, 3]
        Distances between local node [0] and nodes [1], [2] and [3] in each element (derived from get_local_distances)

    Returns
    -------
    gradient : np.array of size [n_elements, 4, 3]
        Gradient in each tetrahedra
    """

    # define projection matrix
    project_matrix = np.hstack([-np.ones((3, 1)), np.eye(3)])

    # solve local_dist * gradient = project_matrix.
    gradient = np.linalg.solve(local_dist, project_matrix[None, :, :])

    # swapaxes
    return np.swapaxes(gradient, 1, 2)
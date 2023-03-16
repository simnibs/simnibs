import os
import time
import h5py
import copy
import logging
import numpy as np
import nibabel as nib

from numpy.linalg import pinv
from scipy import sparse
from scipy.ndimage import zoom
from scipy.optimize import minimize

from .fem import FEMSystem, get_dirichlet_node_index_cog, DirichletBC, dofMap
from .sim_struct import SimuList
from ..mesh_tools import Msh, mesh_io, read_msh
from ..utils.simnibs_logger import logger
from ..utils.file_finder import Templates, SubjectFiles
from ..utils.utils_numba import sumf, sumf2, map_coord_lin_trans, node2elmf, sumf3
from .region_of_interest import _get_gradient, _get_local_distances
from simnibs.simulation import pardiso


class OnlineFEM:
    """
    OnlineFEM class for fast FEM calculations using the Pardiso solver of the MKL library

    Parameters
    ----------
    mesh : Msh object or str
        Head mesh or path to it.
    method : str
        Specify simulation type ('TMS' or 'TES')
    roi : RegionOfInterest object
        Region of interest.
    anisotropy_type : str
        Type of anisotropy for simulation ('scalar', 'vn', 'mc')
    solver_options : str
        Options for theF EM solver (Node or "pardiso))
    fn_results : str
        Filename of results file the data is saved in
    useElements : bool
        True: interpolate the dadt field using positions (coordinates) of elements centroids
        False: interpolate the dadt field using positions (coordinates) of nodes
    grid_spacing : float, optional, default: None
        Grid spacing/interval for A field interpolation
    order : int, optional, default: 1
        Interpolation order of magnetic vector potential of coil (0: nearest neighbor, 1: linear)
    fn_coil : str
        Path to coil file (.ccd or .nii)
    didtmax : float, optional, default: 1e6
        Scaling factor of magnetic vector potential (standard: 1 A/us -> 1e6)
    dataType : int, optional, default: 0
        Calc. magn. of e-field for dataType=0 otherwise return Ex, Ey, Ez
    coil :

    electrode :

    """
    def __init__(self, mesh, method, roi, anisotropy_type="scalar", solver_options=None, fn_results=None,
                 useElements=True, fn_coil=None, dataType=0, coil=None, electrode=None):
        """
        Constructor of the OnlineFEM class
        """
        self.method = method                        # 'TES' or 'TMS'
        self.dataType = dataType                    # calc. magn. of e-field for dataType=0 otherwise return Ex, Ey, Ez
        self.anisotropy_type = anisotropy_type      # 'scalar', 'dir', 'vn', 'mc'
        self.A = None                               # stiffness matrix
        self.b = None                               # rhs
        self.v = None                               # electric potential
        self.solver = None                          # solver
        self.fn_coil = fn_coil                      # filename of TMS coil (.nii)
        self.roi = roi                              # list of ROI instances
        self.fn_results = fn_results                # name of output results file (.hdf5)
        self.ff_templates = Templates()             # initialize file finder templates
        self.force_integrals = None                 # precomputed force integrals for rhs
        self.currents_guess = None                  # guess of electrode currents to accelerate further iterations
        self.coil = coil                            # TMS coil
        self.electrode = electrode                  # TES electrode
        self.amplitude = 1
        self.scale = 1

        if self.electrode is not None:
            if type(electrode) is list:
                self.dirichlet_correction = []
                for _electrode in electrode:
                    self.dirichlet_correction.append(_electrode.dirichlet_correction)
            else:
                self.dirichlet_correction = [electrode.dirichlet_correction]
        else:
            self.dirichlet_correction = False

        if type(self.roi) is not list:
            self.roi = [self.roi]

        self.n_roi = len(self.roi)
        self.e = [0 for _ in range(self.n_roi)]
        self.dadt = None

        # True: interpolate the dadt field using positions (coordinates) of elements centroids
        # False: interpolate the dadt field using positions (coordinates) of nodes
        self.useElements = useElements

        if solver_options is None:
            self.solver_options = "pardiso"
        else:
            self.solver_options = solver_options

        # creating logger
        self.logger = setup_logger(os.path.join(os.path.split(fn_results)[0], "simnibs_simulation_" + time.strftime("%Y%m%d-%H%M%S")))

        # read mesh or store in self
        if type(mesh) is str:
            self.mesh = mesh_io.read_msh(mesh)
        elif type(mesh) == Msh:
            self.mesh = mesh
        else:
            raise TypeError("'mesh' parameter has to be either path to .msh file or SimNIBS mesh object.")

        # get subject specific filenames
        self.ff_subject = SubjectFiles(fnamehead=self.mesh.fn)
        self.fn_tensor_nifti = self.ff_subject.tensor_file

        # get Dirichlet node index where V=0 (close to center of gravity of mesh but not on a surface)
        self.dirichlet_node = get_dirichlet_node_index_cog(mesh=self.mesh)  # self.mesh.nodes.node_number[self.mesh.nodes.node_coord[:, 2].argmin()] #   # mesh.nodes.nr

        # For TMS we use only tetrahedra (elm_type=4). The triangles (elm_type=3) are removed.
        if self.method == "TMS":
            self.mesh = remove_triangles_from_mesh(self.mesh)

        # prepare solver and set self.solver
        self._set_matrices_and_prepare_solver()

        # prepare coil for TMS
        ################################################################################################################
        if method == "TMS":
            if (self.coil is not None and self.fn_coil is not None) or (self.coil is None and self.fn_coil is None):
                raise AssertionError("Provide either filename of TMS coil or Coil object.")

            if coil is None:
                self.coil = Coil(filename=self.fn_coil,
                                 logger=self.logger)
                # scale field to make it more appropriate for matrix solve
                self.coil.a_field *= 1e9
                self.logger.info(f'Loaded coil from file: {self.fn_coil}')

        # prepare electrode for TES
        ################################################################################################################
        if method == "TES" and electrode is None:
            raise AssertionError("Please provide TES electrode object for TES simulations.")

    def update_field(self, electrodes=None, matsimnibs=None, sim_idx_hdf5=0, didt=1e6):
        """
        Calculating and updating electric field for given coil position (matsimnibs) for TMS or electrode position (TES)

        Parameters
        ----------
        electrodes : list of list of np.ndarray [n_sims][n_electrodes]
            List of list of the surface tags or nodes where the currents to be applied for multiple simulations.
        currents : list of list of np.ndarray [n_sims][n_electrodes]
            List of list of the currents in each surface for multiple simulations.
            If Dirichlet correction has to be performed, the input currents are overwritten by an initial guess
            obtained from previous runs to accelerate convergence.
        matsimnibs : np.array of float [4 x 4 x n_sim]
            Tensor containing the coil positions and orientations in SimNIBS space for multiple simulations.
        didt : float
            Rate of change of coil current (A/s) (e.g. 1 A/us = 1e6 A/s)
        sim_idx_hdf5 : int
            Simulation index, to continue to write in .hdf5 file

        Returns
        -------
        e : list of lists [n_sim][n_roi] containing np.arrays of float
            Electric field for queried simulations in ROIs.
        """

        if self.method == "TMS":
            if matsimnibs.ndim < 3:
                matsimnibs = matsimnibs[:, :, np.newaxis]
            n_sim = matsimnibs.shape[2]
        else:
            n_sim = len(electrodes)

        self.e = [[0 for _ in range(self.n_roi)] for _ in range(n_sim)]

        # loop over simulation conditions (multiple coil positions or separate electrode configurations)
        for i_sim in range(n_sim):
            start = time.time()

            # determine RHS
            ############################################################################################################
            if self.method == "TMS":
                self.b = self.set_rhs(matsimnibs=matsimnibs[:, :, i_sim])

            elif self.method == "TES":
                self.b = self.set_rhs(electrode=electrodes[i_sim])

            # solve for potential
            ############################################################################################################
            if self.method == "TES" and self.dirichlet_correction[i_sim]:
                self.v = self.solve_dirichlet_correction(b=self.b,
                                                         electrode=self.electrode[i_sim])
            else:
                self.v = self.solve(b=self.b)

            # calculate e-field
            ############################################################################################################
            for i_roi, r in enumerate(self.roi):
                if self.v is None:
                    self.e[i_sim][i_roi] = None
                    self.logger.log(20, "Warning! Simulation failed! Returning e-field: None!")
                else:
                    self.e[i_sim][i_roi] = r.calc_fields(v=self.v, dadt=self.dadt, dataType=self.dataType)

                if self.method == "TMS":
                    self.e[i_sim][i_roi] *= didt * 1e-9

                # TODO save e-fields outside this function?
                # store results (overwrite if existing)
                if self.fn_results:
                    with h5py.File(self.fn_results, "a") as f:
                        try:
                            if f"{sim_idx_hdf5:04}" in f[f"e/roi_{i_roi}"].keys():
                                del f[f"e/roi_{i_roi}/{sim_idx_hdf5:04}"]
                        except KeyError:
                            pass
                        f.create_dataset(name=f"e/roi_{i_roi}/{sim_idx_hdf5:04}", data=self.e[i_sim][i_roi])

            sim_idx_hdf5 += 1

            stop = time.time()

            self.logger.info(f"Finished simulation #{i_sim + 1}/{n_sim} (time: {(stop-start):.3f}s).")

        return self.e

    def set_rhs(self, electrode=None, matsimnibs=None):
        """
        Set up right hand side (force vector) of equation system.

        Parameters
        ----------
        electrode :

        matsimnibs : np.array of float [4 x 4]
            Matrix containing the coil position and orientation in SimNIBS space.

        Returns
        -------
        b : np.array of float [n_nodes - 1]
            Right hand side of equation system (without Dirichlet node)
        """
        if self.method == "TES":
            # gather electrode currents and associated node indices
            electrodes = []
            currents = []
            for _electrode_array in electrode.electrode_arrays:
                for _ele in _electrode_array.electrodes:
                    electrodes.append(_ele.node_idx + 1)
                    currents.append(_ele.current)

            b = self.fem.assemble_tdcs_neumann_rhs(electrodes=electrodes,       # list of node indices of electrodes
                                                   currents=currents,           # list of electrode currents
                                                   input_type='node',
                                                   areas=self.msh_nodes_areas)

        elif self.method == "TMS":
            # determine magnetic vector potential
            self.dadt = calculate_dadt(a_affine=self.coil.affine,
                                       matsimnibs=matsimnibs,
                                       coordinates=self.coordinates,
                                       a_field=self.coil.a_field,
                                       reshaped_node_numbers=self.reshaped_node_numbersT,
                                       useElements=self.useElements)
            # isotropic
            if self.cond.ndim == 1:
                # b = assemble_force_vector(force_integrals=self.force_integrals,
                #                           reshaped_node_numbers=self.reshaped_node_numbers,
                #                           dadt=self.dadt)

                b = sumf2(x=self.force_integrals, y=self.dadt, w=self.reshaped_node_numbers)

            # anisotropic
            elif self.cond.ndim == 3:
                # # integrate in each node of each element, the value for repeated nodes will be summed together later
                # elm_node_integral = np.zeros((len(self.node_numbers), 4), dtype=np.float64)
                # sigma_dadt = np.einsum('aij, aj -> ai', self.cond, self.dadt)
                #
                # for i in range(4):
                #     elm_node_integral[:, i] = -self.volume * (sigma_dadt * self.gradient[:, i, :]).sum(axis=1)
                #
                # elm_node_integral *= 1e-6  # from m to mm
                #
                # b = np.bincount((self.node_numbers-1).reshape(-1), elm_node_integral.reshape(-1))

                b = sumf3(v=self.volume, dadt=self.dadt, g=self.gradient, nn=(self.node_numbers-1).reshape(-1), c=self.cond)

        else:
            raise NotImplementedError("Simulation method not implemented yet. Method is either 'TMS' or 'TES'.")

        # dof_map_copy = copy.deepcopy(self.dof_map)
        # b, _ = self.bc.apply_to_rhs(self.solver._A, b, dof_map_copy)

        return b

    def normalize_solution(self, v_nodes, channel_id, ele_id, node_area, currents, electrode):
        """
        Normalize electric potential and input currents for later correction step.

        Parameters
        ----------
        v_nodes : list of list of np.arrays of float [n_channel][n_ele][n_nodes]
            Electric potential in electrode nodes for all channels
        channel_id :

        ele_id :

        node_area :

        currents :

        electrode :


        Returns
        -------
        v_ele_norm : list of np.array of float [n_channel][n_ele]
            Normalized difference of electrode voltages for each channel
        currents_ele_norm : np.array of float [n_ele]
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

            # for this channel determine average potentials over electrodes
            v_mean_ele.append([])
            for _ele_id in ele_id_unique:
                ele_mask = _ele_id == ele_id
                mask = ele_mask * channel_mask
                v_mean_ele[i].append(np.sum(node_area[mask] * v_nodes[mask]) / \
                                     np.sum(node_area[mask]))

            v_mean_ele[i] = np.array(v_mean_ele[i])

        print(f"v_mean_ele: {v_mean_ele}")

        v_ele_norm = [np.zeros(electrode.n_ele_per_channel[i]) for i in range(electrode.n_channel)]
        currents_ele_norm = [np.zeros(electrode.n_ele_per_channel[i]) for i in range(electrode.n_channel)]

        # calculate rel. difference of electrode potentials to average electrode potential
        for i_channel in range(electrode.n_channel):
            if len(v_mean_ele[i_channel]) == 1:
                std_scale = np.std(v_mean_channel)
            else:
                std_scale = np.std(v_mean_channel)
            # calculate differences of electrode potentials to channel potential
            v_ele_norm[i_channel] = (v_mean_ele[i_channel] - np.mean(v_mean_ele[i_channel])) / std_scale
            currents_ele_norm[i_channel] = (currents[i_channel] - np.mean(currents[i_channel])) / np.mean(currents[i_channel])

        return v_ele_norm, currents_ele_norm

    def optimize_currents(self, x, electrode, v_node_idx, channel_id_nodes, ele_id_nodes, node_area, signs):
        """

        :return:
        """
        channel_id_unique = np.unique(channel_id_nodes)
        n_channel = len(channel_id_unique)

        # write currents in electrode
        ele_counter = 0
        for _electrode_array in electrode.electrode_arrays:
            for _ele in _electrode_array.electrodes:
                _ele.current = x[ele_counter]
                ele_counter += 1

        # set rhs
        b = self.set_rhs(electrode=electrode)

        # solve
        self.v = self.solve(b)
        v_nodes = self.v[v_node_idx]

        # extract and average voltages of channels and electrodes
        # average and weight potential with node area over electrodes and whole channel (separately)
        v_mean_ele = []  # [n_channel][n_ele]
        v_mean_channel = []  # [n_channel]
        v_var_channel = np.zeros(len(channel_id_unique))   # [n_channel]

        for i, _channel_id in enumerate(channel_id_unique):
            channel_mask = _channel_id == channel_id_nodes

            # average area weighted potential over whole channel
            v_mean_channel.append(np.sum(node_area[channel_mask] * v_nodes[channel_mask]) / \
                                  np.sum(node_area[channel_mask]))

            ele_id_unique = np.unique(ele_id_nodes[channel_mask])

            # for this channel determine average potentials over electrodes
            v_mean_ele.append([])
            for _ele_id in ele_id_unique:
                ele_mask = _ele_id == ele_id_nodes
                mask = ele_mask * channel_mask
                v_mean_ele[i].append(np.sum(node_area[mask] * v_nodes[mask]) / \
                                     np.sum(node_area[mask]))

            v_mean_ele[i] = np.array(v_mean_ele[i])
            v_var_channel[i] = np.var(v_mean_ele[i])

        # evaluate goal function (variance of channel voltages over electrodes)
        var_v = np.linalg.norm(v_var_channel)

        print(x,var_v)

        return var_v

    def current_constraint(self, x, electrode, v_node_idx, channel_id_nodes, ele_id_nodes, node_area, signs):
        current_pos_mask = np.where(signs > 0)[0]
        current_neg_mask = np.where(signs < 0)[0]

        c = np.abs(np.sum(x[current_pos_mask]) - electrode.current_total) + \
            np.abs(np.sum(np.abs(x[current_neg_mask])) - electrode.current_total) + \
            np.abs(np.sum(x))

        # c = np.sum(np.abs(x)) - 2*electrode.current_total

        return c

    # def solve_dirichlet_correction(self, b, electrode):
    #     """
    #
    #     :param electrode:
    #     :return:
    #     """
    #     # create arrays
    #     # v_nodes | v_node_idx | channel_id_nodes | ele_id_nodes | node_area
    #     # ------------------------------------------------------------------
    #     #   ...         4               0                 0           ...
    #     #   ...         1               0                 0           ...
    #     #   ...         15              0                 1           ...
    #     #   ...         3               0                 1           ...
    #     #   ...         9               1                 0           ...
    #     #   ...         8               1                 0           ...
    #     #   ...         22              1                 0           ...
    #     #   ...         0               1                 1           ...
    #     #   ...         10              1                 1           ...
    #     #   ...         11              1                 1           ...
    #
    #     n_channel = electrode.n_channel
    #     channel_id_nodes = []
    #     ele_id_nodes = []
    #     v_node_idx = []
    #     node_area = []
    #
    #     for _electrode_array in electrode.electrode_arrays:
    #         for i_ele, _ele in enumerate(_electrode_array.electrodes):
    #             channel_id_nodes.append(np.array([_ele.channel_id] * _ele.n_nodes))
    #             ele_id_nodes.append(np.array([i_ele] * _ele.n_nodes))
    #             node_area.append(_ele.node_area)
    #             v_node_idx.append(_ele.node_idx)
    #
    #     channel_id_nodes = np.hstack(channel_id_nodes)
    #     ele_id_nodes = np.hstack(ele_id_nodes)
    #     node_area = np.hstack(node_area)
    #     v_node_idx = np.hstack(v_node_idx)
    #
    #     # set initial values and bounds
    #     x0 = [np.zeros(electrode.n_ele_per_channel[i_channel]) for i_channel in range(n_channel)]
    #     bounds = [[] for i_channel in range(n_channel)]
    #
    #     ele_counter_channel = [0] * n_channel
    #     for _electrode_array in electrode.electrode_arrays:
    #         for _ele in _electrode_array.electrodes:
    #             x0[_ele.channel_id][ele_counter_channel[_ele.channel_id]] = _ele.current_init
    #             bounds[_ele.channel_id].append(np.sort(
    #                 [_ele.current_init * 0.1,
    #                  np.sign(electrode.current_mean[_ele.channel_id]) * electrode.current_total]))
    #             ele_counter_channel[_ele.channel_id] += 1
    #     bounds = sum(bounds, [])
    #     x0 = np.hstack(x0).flatten()
    #     signs = np.sign(x0)
    #     self.logger.log(20, f"Starting current optimization with x0: {x0}")
    #
    #     # setup optimization algorithm
    #     res = minimize(self.optimize_currents,
    #                    x0=x0,
    #                    bounds=bounds,
    #                    args=(electrode, v_node_idx, channel_id_nodes, ele_id_nodes, node_area, signs),
    #                    method="trust-constr",  #"SLSQP",
    #                    constraints={"fun": self.current_constraint, "type": "eq", "args": (electrode, v_node_idx, channel_id_nodes, ele_id_nodes, node_area, signs)},
    #                    jac=None, hess=None, hessp=None, tol=None, callback=None, options={"gtol": 1e-5, "xtol":1e-5})
    #     if res.success:
    #         self.logger.log(20, f"Optimized currents: {res.x}")
    #         return self.v
    #     else:
    #         return None

    def solve_dirichlet_correction(self, b, electrode):
        """
        Solve system of equations Ax=b and corrects input currents such that the electrodes with the same channel ID
        have the same potential. Finally, add Dirichlet node (V=0) to solution.

        Parameters
        ----------
        b : np.array of float [n_nodes]
            Right hand side of equation system (with Dirichlet node)
        electrode :


        Returns
        -------
        v : np.array of float [n_nodes]
            Corrected solution (including the Dirichlet node at the right position)
        """
        th_maxrelerr = 0.01
        maxiter = 40

        # create arrays
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

        n_channel = electrode.n_channel
        channel_id_nodes = []
        ele_id_nodes = []
        v_node_idx = []
        node_area = []

        for _electrode_array in electrode.electrode_arrays:
            for i_ele, _ele in enumerate(_electrode_array.electrodes):
                channel_id_nodes.append(np.array([_ele.channel_id] * _ele.n_nodes))
                ele_id_nodes.append(np.array([i_ele] * _ele.n_nodes))
                node_area.append(_ele.node_area)
                v_node_idx.append(_ele.node_idx)

        channel_id_nodes = np.hstack(channel_id_nodes)
        ele_id_nodes = np.hstack(ele_id_nodes)
        node_area = np.hstack(node_area)
        v_node_idx = np.hstack(v_node_idx)

        # mean current over electrodes [n_channel]
        I_mean = electrode.current_mean

        # read currents from electrode
        I = [np.zeros(electrode.n_ele_per_channel[i_channel]) for i_channel in range(n_channel)]
        I_sign = [np.zeros(electrode.n_ele_per_channel[i_channel]) for i_channel in range(n_channel)]
        bounds = [[] for i_channel in range(n_channel)]
        ele_counter_channel = [0] * n_channel
        for _electrode_array in electrode.electrode_arrays:
            for _ele in _electrode_array.electrodes:
                I[_ele.channel_id][ele_counter_channel[_ele.channel_id]] = _ele.current
                I_sign[_ele.channel_id][ele_counter_channel[_ele.channel_id]] = _ele.current_sign
                bounds[_ele.channel_id].append(np.sort(
                                    [_ele.current * 0.1,
                                     np.sign(electrode.current_mean[_ele.channel_id]) * electrode.current_total]))
                ele_counter_channel[_ele.channel_id] += 1

        # Solve iteratively until maxrelerr is reached
        # Iteration: 0

        # solve
        v = self.solve(b)
        v_nodes = v[v_node_idx]

        # v_norm: [n_channel][n_ele], I_norm: [n_channel][n_ele]
        v_norm, I_norm = self.normalize_solution(v_nodes=v_nodes,
                                                 channel_id=channel_id_nodes,
                                                 ele_id=ele_id_nodes,
                                                 node_area=node_area,
                                                 currents=I,
                                                 electrode=electrode)

        for i_channel in range(n_channel):
            v_norm[i_channel] = v_norm[i_channel][np.newaxis, :]
            I_norm[i_channel] = I_norm[i_channel][np.newaxis, :]

        j = 0
        maxrelerr = np.max([np.max(np.abs(v_norm[k][-1])) for k in range(n_channel)])
        while maxrelerr > th_maxrelerr:
            print(f"iter: {j} current: {I}, error: {maxrelerr}")
            j += 1

            if j > maxiter:
                print('warning: did not converge after ' + str(maxiter) + ' iterations')
                return None

            if j == 1:
                # It 1: make small step hopefully in the right direction
                I = []
                for i_channel in range(n_channel):
                    if electrode.n_ele_per_channel[i_channel] == 1:
                        I.append(np.array([I_mean[i_channel]]))
                    else:
                        I.append(-0.1 * np.sign(I_mean[i_channel]) * v_norm[i_channel][0, :] / np.max(np.abs(v_norm[i_channel][0, :])))

            else:
                # It 2 ... use gradient descent
                for i_channel in range(n_channel):
                    if electrode.n_ele_per_channel[i_channel] == 1:
                        I[i_channel] = I[i_channel]
                    else:
                        denom = (v_norm[i_channel][-1, :] - v_norm[i_channel][-2, :])
                        denom += np.sign(denom) * maxrelerr / 15  # regularize: avoid too small denominators to gain stability
                        I[i_channel] = I_norm[i_channel][-1] - (I_norm[i_channel][-1, :] - I_norm[i_channel][-2, :]) / denom * v_norm[i_channel][-1, :]

            # convert back from I_norm to I
            I = [I_mean[i_channel] * (I[i_channel] - np.mean(I[i_channel]) + 1) for i_channel in range(n_channel)]

            # write currents in electrode
            ele_counter_channel = [0] * n_channel
            for _electrode_array in electrode.electrode_arrays:
                for _ele in _electrode_array.electrodes:
                    _ele.current = I[_ele.channel_id][ele_counter_channel[_ele.channel_id]]
                    ele_counter_channel[_ele.channel_id] += 1

            # update rhs
            b = self.set_rhs(electrode=electrode)

            # solve
            v = self.solve(b)
            v_nodes = v[v_node_idx]

            # normalize
            _v_norm, _I_norm = self.normalize_solution(v_nodes=v_nodes,
                                                       channel_id=channel_id_nodes,
                                                       ele_id=ele_id_nodes,
                                                       node_area=node_area,
                                                       currents=I,
                                                       electrode=electrode)

            # append
            for i_channel in range(n_channel):
                v_norm[i_channel] = np.vstack((v_norm[i_channel], _v_norm[i_channel][np.newaxis, :]))
                I_norm[i_channel] = np.vstack((I_norm[i_channel], _I_norm[i_channel][np.newaxis, :]))

            # update error
            maxrelerr = np.max([np.max(np.abs(v_norm[i_channel][-1])) for i_channel in range(n_channel)])

        # test if signs are correct otherwise restart
        if not np.array([(np.sign(I[i_channel]) == I_sign[i_channel]).all() for i_channel in range(n_channel)]).all():
            print(f"Final currents have wrong signs!")
            print(f"I: { *I, }")
            return None

        # add optimal currents to CurrentEstimator to improve estimation in further iterations
        electrode_pos = [_electrode_array.electrode_pos for _electrode_array in electrode.electrode_arrays]
        ele_counter_channel = [0] * n_channel

        for _electrode_array in electrode.electrode_arrays:
            for _ele in _electrode_array.electrodes:
                _ele.current_estimator.add_training_data(electrode_pos=np.hstack(electrode_pos),
                                                         current=I[_ele.channel_id][ele_counter_channel[_ele.channel_id]])
                ele_counter_channel[_ele.channel_id] += 1

        self.logger.log(20, f"Optimal currents: { *I, }")

        # reset to original currents
        for _electrode_array in electrode.electrode_arrays:
            for _ele in _electrode_array.electrodes:
                _ele.current = _ele.current_init

        # import matplotlib.pyplot as plt
        # import matplotlib
        # matplotlib.use("Qt5Agg")
        #
        # t = np.arange(j + 1)
        # for i_channel in range(n_channel):
        #     v_normP = v_norm[i_channel]
        #     I_normP = I_norm[i_channel]
        #     I_normP = I_norm[i_channel]
        #
        #     fig, axs = plt.subplots(5, 1)
        #     for i in range(electrode.n_ele_per_channel[i_channel]):
        #         axs[0].plot(t, v_normP[:, i])
        #         axs[1].plot(t, I_normP[:, i])
        #     axs[2].plot(t, np.log10(np.mean(np.abs(v_normP), axis=1)))
        #     axs[2].plot(t, np.log10(np.max(np.abs(v_normP), axis=1)))
        #
        #     axs[3].plot(t, np.log10(np.max(np.abs(v_normP), axis=1)))
        #
        #     axs[0].set_title('v_norm')
        #     axs[1].set_title('I_norm')
        #     axs[2].set_title('max and mean relative error of v_norm (log10)')
        #
        #     fig.tight_layout()
        #     plt.show()

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

        logger.disabled = True

        # remove dirichlet node
        b_reduced = np.delete(b, self.dirichlet_node-1)

        # solve equation system
        x = self.solver.solve(b_reduced)

        # add Dirichlet node to solution
        v = np.insert(x, self.dirichlet_node-1, 0)

        logger.disabled = False

        return np.squeeze(v)

    def _set_matrices_and_prepare_solver(self):
        """
        Set matrices and initialize the pardiso solver for update_position in self.solver.
        """

        # TODO: TES works same as TMS (use assemble stiffnes matrix function here as well and only rhs precomputation (force integral) is different)
        if self.method == "TMS":
            # prepare conductivity (scalar or anisotropic)
            self.simulist = SimuList(mesh=self.mesh)
            self.simulist.anisotropy_type = self.anisotropy_type
            self.simulist.fn_tensor_nifti = self.fn_tensor_nifti
            self.cond = self.simulist.cond2elmdata()
            self.cond = self.cond.value.squeeze()

            if self.cond.ndim == 2:
                self.cond = self.cond.reshape(-1, 3, 3)

            # self.fem = FEMSystem.tms(mesh=self.mesh, cond=cond, solver_options=self.solver_options)
            # self.fem.prepare_solver()
            # self.solver = self.fem._solver

            self.node_numbers = self.mesh.elm.node_number_list  # self.node_numbers
            useElements = self.useElements
            node_coordinates = self.mesh.nodes.node_coord.T
            tag1 = self.mesh.elm.tag1

            # get the number of nodes
            number_of_nodes = node_coordinates.shape[1]

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

            # get the conductivity values of each tetrahedron
            # tag1 indexes starting from 1, the array cond_table indexes starting from 0
            # Use isotropic conductivities
            # conductivity = np.array([float(c.value) for c in standard_cond()])[tag1 - 1]

            # set the force integrals. We use the force integrals to assemble the right hand side force vector
            if self.cond.ndim == 1:
                self.force_integrals = get_force_integrals(self.volume, self.gradient, self.cond)

            #  assemble the left hand side stiffness matrix
            # (volume, gradient, conductivity, node_numbers, number_of_nodes, dirichlet_node
            self.A = assemble_stiffness_matrix(volume=self.volume,
                                               gradient=self.gradient,
                                               conductivity=self.cond,
                                               node_numbers=self.node_numbers,
                                               number_of_nodes=number_of_nodes)

            # self.A_reduced = copy.deepcopy(self.A)
            # dof_map_copy = copy.deepcopy(self.dof_map)
            # self.A_reduced, _ = self.bc.apply_to_matrix(self.A_reduced, dof_map_copy)

            self.A = delete_row_csr(self.A, self.dirichlet_node-1)
            self.A = delete_col_csr(self.A, self.dirichlet_node-1)
            self.solver = pardiso.Solver(self.A)

        elif self.method == "TES":
            # Calculate node areas for whole mesh
            self.msh_nodes_areas = self.mesh.nodes_areas()

            # prepare FEM
            self.simulist = SimuList(mesh=self.mesh)
            self.simulist.anisotropy_type = self.anisotropy_type
            self.simulist.fn_tensor_nifti = self.fn_tensor_nifti
            cond = self.simulist.cond2elmdata()

            self.fem = FEMSystem.tdcs_neumann(mesh=self.mesh,
                                              cond=cond,
                                              ground_electrode=self.dirichlet_node,
                                              solver_options=self.solver_options,
                                              input_type='node')
            self.fem.prepare_solver()
            self.solver = self.fem._solver


class Coil:
    """
    Coil class for TMS calculations. Contains precalculated magnetic vector potential and resolution information.

    Parameters
    ----------
    filename : str
        Path to the coil file (.ccd or .nii format)
    didtmax : float
        Maximum rate of change of coil current (A/us)
    logger : logger object
        Logger

    Attributes
    ----------
    affine : np.array of float [4 x 4]
        Affine matrix describing resolution, location and orientation of untransformed magnetic vector potential
    a_field : np.array of float [3 x N_x x N_y x N_z]
        Magnetic vector potential (A_x, A_y, A_z)
    """
    def __init__(self, filename, logger=None):
        """
        Constructor of Coil class
        """
        self.filename = filename

        if not isinstance(filename, str):
            raise NameError(
                'Failed to parse input volume (not string or nibabel nifti1 volume)')

        # load the A field and affine matrix
        if filename.endswith('.nii.gz') or filename.endswith('.nii'):
            # update the A field and affine matrix
            coil = nib.load(filename)
            A = np.moveaxis(np.asanyarray(coil.dataobj), -1, 0)
            affine = coil.affine

            if logger is not None:
                logger.debug('Load A from coil file: ' + filename)
        else:
            raise ValueError('Coil file must be a nifti file')

        self.a_field = np.asfortranarray(A)
        self.affine = affine


def calculate_dadt(a_affine, matsimnibs, coordinates, a_field, reshaped_node_numbers, useElements):
    """
    Interpolates the dadt field from a coil file.

    Parameters
    ----------
    a_affine : np.array of float [4 x 4]
        Affine matrix describing resolution, location and orientation of untransformed magnetic vector potential
    matsimnibs : np.array of float [4 x 4]
        Matrix containing the coil positions and orientations in SimNIBS space.
    coordinates :

    a_field : np.array of float [3 x N_x x N_y x N_z]
        Magnetic vector potential (A_x, A_y, A_z)
    reshaped_node_numbers : np.array of int [4 * n_elements + 1]
        Flattened node number list (connectivity matrix)
        node_numbers.reshape(-1)
    useElements : bool
        Get coordinates in the element centers (True) or of the nodes (False)

    Returns
    -------
    dadt : np.array of float [3 x N_x x N_y x N_z]
        Magnetic vector potential (A_x, A_y, A_z)
    """

    # Get the affine transformation from the "coordinates" to the coil space (defined in coil file)
    trans = np.dot(pinv(a_affine), pinv(matsimnibs))
    M1 = np.ascontiguousarray(trans[:3, :3])
    t = np.ascontiguousarray(trans[:3, 3])
    M2 = np.ascontiguousarray(matsimnibs[:3, :3])

    if not useElements:
        M2 *= 0.25
        dadt = map_coord_lin_trans(a_field, coordinates, M1, t, M2)
        dadt = node2elmf(dadt, reshaped_node_numbers)
    else:
        dadt = map_coord_lin_trans(a_field, coordinates, M1, t, M2).T

    return dadt


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

    # integrate in each node of each element, the value for repeated nodes will be summed
    # together later
    node_integrals = np.zeros(force_integrals.shape[:2], order='C')

    # force_integrals: (4, number_of_elements, 3); dadt: (number_of_elements, 3); node_integrals: (4, number_of_elements)
    # node_integrals = np.einsum('ijk,jk->ij', force_integrals, dadt)
    sumf(force_integrals, dadt, node_integrals)

    # Assembles the right hand side for TMS.
    # forcevec: (number_of_nodes,), reshaped_node_numbers: (number_of_elements*4,), node_integrals: (4, number_of_elements)
    # forcevec = bincount_nb(reshaped_node_numbers, node_integrals.reshape(-1)).reshape(-1, 1)

    # keep np.bincount to make the testing more stable. May change it back to bincount_nb() for performance reason.
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


def assemble_stiffness_matrix(volume, gradient, conductivity, node_numbers, number_of_nodes):
    """
    Assembly of the l.h.s stiffness matrix. Based in the OptVS algorithm in Cuvelier et. al. 2016.

    Cuvelier, F., Japhet, C., & Scarella, G. (2016). An efficient way to assemble finite element matrices in
    vector languages. BIT Numerical Mathematics, 56(3), 833-864.

    Parameters
    ----------
    volume : np.array of float [n_elements]
        Volume of the tetrahedra
    gradient : np.array of size [n_elements, 4, 3]
        Gradient in each tetrahedra
    conductivity : np.array of float [n_elements]
        Electrical conductivity value (in S/m) assigned to each tetrahedra
    node_numbers : np.array of size [n_elements x 4]
        Node number list (connectivity list)
    number_of_nodes : int
        Number of nodes
    bc : DirichletBC object
        Dirichlet boundary condition object

    Returns
    -------
    stiffmat : scipy.sparse [(n_nodes-1) x (n_nodes-1)]
        Stiffness matrix in sparse (CSR) format (without Dirichlet node)
    """

    # gradient: (number_of_elements, 4, 3),
    # volume: (number_of_elements,),
    # conductivity: (number_of_elements,),
    # node_numbers: (number_of_elements, 4),
    # number_of_nodes: number_of_nodes
    # (number_of_elements, 4)
    index = node_numbers - 1

    dim = np.arange(index.shape[1])
    idx = np.repeat(dim, len(dim))
    idy = np.tile(dim, len(dim))

    # Simplify the integration using commutative law of dot product (elementary-wise product in this case
    # units == 'mm': * 1e6 from the gradient operator, 1e-9 from the volume
    if conductivity.ndim == 1:
        factor = (volume * conductivity * 1e-3)[:, None, None] * gradient
    elif conductivity.ndim == 3:
        factor = volume[:, None, None] * np.einsum('aij, ajk -> aik', gradient, conductivity) * 1e-3

    # if cond.ndim == 1:
    #     vGc = vols[:, None, None]*G*cond[:, None, None]
    # elif cond.ndim == 3:
    #     vGc = vols[:, None, None]*np.einsum('aij, ajk -> aik', G, cond)

    stiffmat = sparse.coo_matrix((number_of_nodes, number_of_nodes), dtype='float64')

    stiffmat.data = (factor @ np.swapaxes(gradient, 1, 2)).flatten()
    stiffmat.row = (index[:, idx]).flatten()
    stiffmat.col = (index[:, idy]).flatten()

    stiffmat = stiffmat.tocsr()
    stiffmat.eliminate_zeros()

    # Make stiffmat symmetric if necessary
    # stiffmat = (stiffmat + stiffmat.T) * 0.5
    stiffmat = stiffmat.sorted_indices()

    # remove row and column for the dirichlet boundary condition
    # stiffmat = delete_row_csr(stiffmat, dirichlet_node-1)
    # stiffmat = delete_cols_csr(stiffmat, dirichlet_node-1)

    return stiffmat


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


def remove_triangles_from_mesh(mesh):
    """
    Load the mesh file. For the E field calculation, we use only tetrahedra (elm_type=4).
    The triangles (elm_type=3) are removed.

    Parameters
    ----------
    mesh : Msh object
        Mesh object

    Returns
    -------
    mesh : Msh object
        Mesh object without triangles
    """

    # keep the tetrahedra in the mesh file (triangles are removed)
    mesh = mesh.crop_mesh(elm_type=4)

    assert mesh.elm.node_number_list.max() == mesh.nodes.nr

    # get the nodal index and value with dirichlet boundary condition
    index_naught = mesh.nodes.node_coord[:, 2].argmin()

    # switch the node with dirichlet bc and the last node
    mesh.nodes.node_coord[[index_naught, -1]] = mesh.nodes.node_coord[[-1, index_naught]]

    # update the node_number_list
    mask1 = np.nonzero(mesh.elm.node_number_list == index_naught + 1)
    mask2 = np.nonzero(mesh.elm.node_number_list == mesh.nodes.nr)
    mesh.elm.node_number_list[mask1] = mesh.nodes.nr
    mesh.elm.node_number_list[mask2] = index_naught + 1

    return mesh


def get_mesh_with_tetrahedra(mesh_file, logger=None):
    """
    Load the mesh file. For the E field calculation, we use only tetrahedra (elm_type=4).
    The triangles (elm_type=3) are removed.

    Parameters
    ----------
    mesh_file : str
        Filename (incl. path) to .msh file.
    logger : logger object
        Logger

    Returns
    -------
    mesh : Msh object
        Mesh object
    """
    # read_mesh: load the mesh file; crop_mesh: keep the tetrahedra in the mesh file (triangles are removed)
    mesh = read_msh(mesh_file).crop_mesh(elm_type=4)

    assert mesh.elm.node_number_list.max() == mesh.nodes.nr

    # get the nodal index and value with dirichlet boundary condition
    index_naught = mesh.nodes.node_coord[:, 2].argmin()

    # switch the node with dirichlet bc and the last node
    mesh.nodes.node_coord[[index_naught, -1]] = mesh.nodes.node_coord[[-1, index_naught]]

    # update the node_number_list
    mask1 = np.nonzero(mesh.elm.node_number_list == index_naught + 1)
    mask2 = np.nonzero(mesh.elm.node_number_list == mesh.nodes.nr)
    mesh.elm.node_number_list[mask1] = mesh.nodes.nr
    mesh.elm.node_number_list[mask2] = index_naught + 1

    if logger is not None:
        logger.info('Loaded mesh file: ' + mesh_file)

    return mesh


def setup_logger(logname, filemode='w', format='[ %(name)s ] %(levelname)s: %(message)s', datefmt='%H:%M:%S'):
    """
    Setup logger.

    Parameters
    ----------
    logname : str
        Filename of logfile
    filemode : str, optional, default: 'w'
        'a' append or 'w' overwrite existing logfile.
    format : str, optional, default: '[ %(name)s ] %(levelname)s: %(message)s'
        format of the output message.
    datefmt : str, optional, default: '%H:%M:%S'
        Format of the output time.

    Returns
    -------
    logger : logger instance
        Logger using loging module
    """

    logging.basicConfig(filename=logname,
                        filemode=filemode,
                        format=format,
                        datefmt=datefmt,
                        level=logging.INFO)

    logger = logging.getLogger("simnibs")

    return logger


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

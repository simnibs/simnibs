import os
import copy
import time
import copy
import h5py
import logging
import datetime
import numpy as np
import nibabel as nib
import scipy.ndimage.morphology as mrph

from numpy.linalg import eig
from scipy.optimize import direct, Bounds

from ..mesh_tools import Msh
from ..mesh_tools import mesh_io
from ..mesh_tools import surface
from ..simulation.sim_struct import ELECTRODE, SimuList
from ..simulation.fem import FEMSystem, get_dirichlet_node_index_cog
from ..simulation.onlinefem import OnlineFEM
from ..utils.file_finder import Templates, SubjectFiles
from ..utils.transformations import subject2mni_coords
from ..utils.ellipsoid import Ellipsoid, subject2ellipsoid, ellipsoid2subject

SIMNIBSDIR = os.path.split(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))[0]

# TODO: adapt matlab_tools/opt_struct.m and include new class
# 1:1 Mapping der Elektroden -> von Neumann Kanal 1 -> Sink 1, Kanal 2 -> Sink 2
# 1 Kanal -> mehrere sink elektroden -> fake Dirichlet -> Kanal 1 -> Sink 1,2,3 (zusammengeschaltet) gleiche Spannung


class TESoptimize():
    """
    Defines a TES optimization problem using a direct approach

    Parameters
    --------------
    electrode : Electrode Object
        Electrode object containing ElectrodeArray instances
        (see /simulation/array_layout.py for pre-implemented examples)
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float (optional)
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int (optional)
        Maximum number of active electrodes. Default: no maximum
    init_pos : str or list of str or list of str and np.ndarray of float [3]
        Initial positions of movable Electrode arrays (for each movable array)
    fn_eeg_cap : str, optional, default: 'EEG10-10_UI_Jurak_2007.csv'
        Filename of EEG cap to use for initial position (without path)
        - 'EEG10-10_UI_Jurak_2007.csv'
        - 'easycap_BC_TMS64_X21.csv'
        - 'EEG10-10_Cutini_2011.csv'
        - 'EEG10-10_Neuroelectrics.csv'
        - 'EEG10-20_extended_SPM12.csv'
        - 'EEG10-20_Okamoto_2004.csv'
    roi : list of RegionOfInterest class instances
        Region of interest(s) the field is evaluated in.
    anisotropy_type : str
        Specify type of anisotropy for simulation ('scalar', 'vn' or 'mc')
    weights : np.array of float [n_roi]
        Weights for optimizer for ROI specific goal function weighting
    min_electrode_distance : float, optional, default: None
        Minimum electrode distance to ensure during optimization (in mm).


    plot : bool, optional, default: False
        Plot configurations in output folder for visualization and control
    name: str (optional)
        Name of optimization problem. Default: optimization
    target: list of TDCStarget objects (optional)
        Targets for the optimization. Default: no target
    avoid: list of TDCSavoid objects
        list of TDCSavoid objects defining regions to avoid

    Attributes
    --------------
    electrode : Electrode Object
        Electrode object containing ElectrodeArray instances
        (see /simulation/array_layout.py for pre-implemented examples)
    ellipsoid : Ellipsoid Object
        Best fitting ellipsoid to valid skin reagion (used for coordinate system definition)
    msh_nodes_areas : np.ndarray of float [n_nodes_msh]
        Areas of nodes
    node_idx_msh : np.ndarray of int [n_nodes_skin]
        Indices of skin surface nodes in global msh
    """

    def __init__(self,
                 mesh,
                 electrode,
                 roi,
                 skin_mask=None,
                 init_pos=None,
                 fn_eeg_cap=None,
                 solver_options=None,
                 optimizer_options=None,
                 weights=None,
                 anisotropy_type=None,
                 min_electrode_distance=None,
                 plot=False,
                 goal="mean",
                 optimizer="direct",
                 output_folder=None):
        """
        Constructor of TESoptimize class instance
        """
        # setup output folders, logging and IO
        ################################################################################################################
        assert type(output_folder) is str, "Please provide an output folder to save optimization results in."
        self.output_folder = output_folder
        self.plot_folder = os.path.join(self.output_folder, "plots")
        self.plot = plot
        self.fn_results_hdf5 = os.path.join(self.output_folder, "opt.hdf5")

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        if not os.path.exists(self.plot_folder):
            os.makedirs(self.plot_folder)

        # setup logger
        self.logger = setup_logger(os.path.join(output_folder, "simnibs_simulation_" + time.strftime("%Y%m%d-%H%M%S")))
        self.logger.log(20, "Setting up output folders, logging and IO ...")

        # setup headmodel
        ################################################################################################################
        self.logger.log(20, "Setting up headmodel ...")

        if fn_eeg_cap is None:
            fn_eeg_cap = 'EEG10-20_Okamoto_2004.csv'  # 'EEG10-10_UI_Jurak_2007.csv'
        else:
            fn_eeg_cap = os.path.split(fn_eeg_cap)[1]

        # read mesh or store in self
        if type(mesh) is str:
            self.mesh = mesh_io.read_msh(mesh)
        elif type(mesh) == Msh:
            self.mesh = mesh
        else:
            raise TypeError("msh has to be either path to .msh file or SimNIBS mesh object.")

        # get subject specific filenames
        self.ff_templates = Templates()
        self.ff_subject = SubjectFiles(fnamehead=self.mesh.fn)

        # set dirichlet node to closest node of center of gravity of head model (indexing starting with 1)
        self.dirichlet_node = get_dirichlet_node_index_cog(mesh=self.mesh)

        # Calculate node areas for whole mesh
        self.mesh_nodes_areas = self.mesh.nodes_areas()

        self.ellipsoid = Ellipsoid()
        self.fn_eeg_cap = fn_eeg_cap

        # set skin mask of upper head
        if skin_mask is None:
            self.fn_electrode_mask = self.ff_templates.mni_volume_upper_head_mask
        else:
            self.fn_electrode_mask = skin_mask

        # relabel internal air
        self.mesh_relabel = relabel_internal_air(m=self.mesh,
                                                 subpath=os.path.split(self.mesh.fn)[0],
                                                 label_skin=1005,
                                                 label_new=1099,
                                                 label_internal_air=501)

        # create skin surface
        self.skin_surface = surface.Surface(mesh=self.mesh_relabel, labels=1005)

        # determine point indices where the electrodes may be applied during optimization
        self.skin_surface = self.valid_skin_region(skin_surface=self.skin_surface, mesh=self.mesh_relabel)

        # get mapping between skin_surface node indices and global mesh nodes
        self.node_idx_msh = np.where(np.isin(self.mesh.nodes.node_coord, self.skin_surface.nodes).all(axis=1))[0]

        # fit optimal ellipsoid to valid skin points
        self.ellipsoid.fit(points=self.skin_surface.nodes)

        # setup ROI
        ################################################################################################################
        self.logger.log(20, "Setting up ROI ...")
        if type(roi) is not list:
            roi = [roi]

        self.roi = roi
        self.n_roi = len(roi)

        # setup electrode
        ################################################################################################################
        self.logger.log(20, "Setting up electrodes ...")
        if type(electrode) is not list:
            electrode = [electrode]

        self.electrode = electrode
        self.electrode_pos_opt = None
        self.min_electrode_distance = min_electrode_distance
        self.dirichlet_correction = False  # will be checked later if required

        # number of independent stimulation channels (on after another)
        self.n_channel_stim = len(electrode)

        # list containing the number of freely movable arrays for each channel [i_channel_stim]
        self.n_ele_free = [len(ele.electrode_arrays) for ele in self.electrode]

        # list containing beta, lambda, alpha for each freely movable array and for each stimulation channel
        self.electrode_pos = [[np.zeros(3) for _ in range(n_ele_free)] for n_ele_free in self.n_ele_free]

        # set initial positions
        if self.n_channel_stim > 1 and type(init_pos) is str:
            raise AssertionError("Please provide a list of initial positions for each stimulation channel"
                                 "containing a list of initial positions for each freely movable array.")

        if self.n_channel_stim == 1 and type(init_pos) is str:
            init_pos = [[init_pos]]

        if self.n_channel_stim == 1 and type(init_pos) is np.ndarray:
            init_pos = [[init_pos]]

        self.init_pos = init_pos
        self.init_pos_subject_coords = [[] for _ in range(self.n_channel_stim)]

        init_pos_list = [["C3", "C4", "F3", "P4"],          # init defaults for first i_channel_stim
                         ["Fz", "Pz", "P3", "F4"]]          # init defaults for second i_channel_stim

        # collect all electrode currents in list of np.array [i_channel_stim][i_channel_ele]
        self.current = [np.zeros(len(np.unique(self.electrode[i_channel_stim].channel_id)))
                        for i_channel_stim in range(self.n_channel_stim)]

        self.dirichlet_correction = [False] * self.n_channel_stim
        for i_channel_stim in range(self.n_channel_stim):
            for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
                if _electrode_array.dirichlet_correction:
                    self.dirichlet_correction[i_channel_stim] = True

                for _electrode in _electrode_array.electrodes:
                    self.current[i_channel_stim][_electrode.channel_id] += _electrode.current

        # set initial positions of electrodes if nothing is provided
        assert self.n_channel_stim <= len(init_pos_list), "Please provide initial electrode positions."

        if init_pos is None:
            self.init_pos = [0 for _ in range(self.n_channel_stim)]
            for i_channel_stim in range(self.n_channel_stim):
                if self.n_ele_free[i_channel_stim] > len(init_pos_list[i_channel_stim]):
                    raise NotImplementedError("Please specify initial coordinates or EEG electrode positions for each"
                                              "freely movable electrode array (init_pos)!")
                self.init_pos[i_channel_stim] = [init_pos_list[i_channel_stim][i]
                                                 for i in range(self.n_ele_free[i_channel_stim])]

        # get subject coordinates of initial positions
        for i_channel_stim in range(self.n_channel_stim):
            # user provided EEG electrode position as str (e.g. "C3", ...)
            if type(self.init_pos[i_channel_stim][0]) is str:
                for eeg_pos in self.init_pos[i_channel_stim]:
                    tmp = ELECTRODE()
                    tmp.centre = eeg_pos
                    tmp.substitute_positions_from_cap(cap=self.ff_subject.get_eeg_cap(cap_name=self.fn_eeg_cap))
                    self.init_pos_subject_coords[i_channel_stim].append(tmp.centre)
            # user provided coordinates in subject space as np.array
            else:
                self.init_pos_subject_coords[i_channel_stim] = self.init_pos[i_channel_stim]

            # transform initial positions from subject to ellipsoid space
            for i_ele_free, coords in enumerate(self.init_pos_subject_coords[i_channel_stim]):
                # get closest point idx on subject surface
                point_idx = np.argmin(np.linalg.norm(coords-self.skin_surface.nodes, axis=1))

                # electrode positon in ellipsoid space
                self.electrode_pos[i_channel_stim][i_ele_free][:2] = self.ellipsoid.cartesian2jacobi(
                    coords=self.ellipsoid.ellipsoid2cartesian(
                        coords=subject2ellipsoid(
                            coords=self.skin_surface.nodes[point_idx, :],
                            normals=self.skin_surface.nodes_normals[point_idx, :],
                            ellipsoid=self.ellipsoid)))

                # set initial orientation alpha to zero
                self.electrode_pos[i_channel_stim][i_ele_free][2] = 0.

        # setup optimization
        ################################################################################################################
        self.logger.log(20, "Setting up optimization algorithm ...")

        # equal ROI weighting if None is provided
        if weights is None:
            weights = np.ones(len(self.roi)) / len(self.roi)

        assert len(weights) == len(self.roi), "Number of weights has to match the number ROIs"

        self.goal = goal
        self.optimizer = optimizer
        self.weights = weights

        # parameter bounds for optimizer
        bounds = Bounds(lb=[-np.pi / 2, -np.pi, -np.pi] * np.sum(self.n_ele_free),
                        ub=[np.pi / 2, np.pi, np.pi] * np.sum(self.n_ele_free))

        # set standard options for optimizer
        if self.optimizer == "direct":
            self.optimizer_options = {"bounds": bounds,
                                      "vol_tol": 1. / 3600000000. * 3 * np.sum(self.n_ele_free),
                                      "len_tol": 1. / 3600000000.,
                                      "f_min_rtol": 1e-12,
                                      "maxiter": 1000}

        # insert user specific options
        if optimizer_options is not None:
            for key in optimizer_options:
                self.optimizer_options[key] = optimizer_options[key]

        # self.max_total_current = max_total_current
        # self.max_individual_current = max_individual_current
        # self.max_active_electrodes = max_active_electrodes

        # setup FEM
        ################################################################################################################
        if anisotropy_type is None:
            anisotropy_type = "scalar"

        if solver_options is None:
            solver_options = "pardiso"

        # prepare FEM
        self.ofem = OnlineFEM(mesh=self.mesh,
                              electrode=electrode,
                              method="TES",
                              roi=self.roi,
                              anisotropy_type=anisotropy_type,
                              solver_options=solver_options,
                              fn_results=self.fn_results_hdf5,
                              useElements=True,
                              dataType=0)

        # self.logger.log(20, 'Preparing FEM')
        # self.simulist = SimuList(mesh=self.mesh)
        # self.simulist.anisotropy_type = self.anisotropy_type
        # cond = self.simulist.cond2elmdata()
        # self.fem = FEMSystem.tdcs_neumann(mesh=self.mesh,
        #                                   cond=cond,
        #                                   ground_electrode=self.dirichlet_node,
        #                                   solver_options=self.solver_options,
        #                                   input_type='node')
        # self.fem.prepare_solver()

        # log summary
        ################################################################################################################
        self.logger.log(25, f"="*100)
        self.logger.log(25, f"headmodel:             {self.mesh.fn}")
        self.logger.log(25, f"n_roi:                 {self.n_roi}")
        self.logger.log(25, f"anisotropy type:       {self.ofem.anisotropy_type}")
        self.logger.log(25, f"n_channel_stim:        {self.n_channel_stim}")
        self.logger.log(25, f"fn_eeg_cap:            {self.fn_eeg_cap}")
        self.logger.log(25, f"fn_electrode_mask:     {self.fn_electrode_mask}")
        self.logger.log(25, f"FEM solver options:    {self.ofem.solver_options}")
        self.logger.log(25, f"dirichlet_correction:  {self.dirichlet_correction}")
        self.logger.log(25, f"optimizer:             {self.optimizer}")
        self.logger.log(25, f"goal:                  {self.goal}")
        self.logger.log(25, f"weights:               {self.weights}")
        self.logger.log(25, f"output_folder:         {self.output_folder}")
        self.logger.log(25, f"fn_results_hdf5:       {self.fn_results_hdf5}")

        if self.optimizer_options is not None:
            for key in self.optimizer_options:
                if key != "bounds":
                    self.logger.log(25, f"{key}:              {self.optimizer_options[key]}")

        for i_channel_stim in range(self.n_channel_stim):
            self.logger.log(25, f"Stimulation: {i_channel_stim} (n_ele_free: {self.n_ele_free[i_channel_stim]})")
            self.logger.log(25, f"---------------------------------------------")

            for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
                self.logger.log(25, f"Electrode array [i_channel_stim][i_array]: [{i_channel_stim}][{i_array}]")
                self.logger.log(25, f"\tn_ele: {_electrode_array.n_ele}")
                self.logger.log(25, f"\tinit_pos: {self.init_pos[i_channel_stim][i_array]}")
                # self.logger.log(25, f"\tcenter: {_electrode_array.center}")
                # self.logger.log(25, f"\tradius: {_electrode_array.radius}")
                # self.logger.log(25, f"\tlength_x: {_electrode_array.length_x}")
                # self.logger.log(25, f"\tlength_y: {_electrode_array.length_y}")

        self.logger.log(25, f"=" * 100)

        # perform optimization
        ################################################################################################################
        self.optimize()

        # temporary plot functions
        ################################################################################################################
        if self.plot:
            import matplotlib
            matplotlib.use('Qt5Agg')
            import matplotlib.pyplot as plt

            theta = np.linspace(0, np.pi, 180)
            phi = np.linspace(0, 2 * np.pi, 360)

            beta = np.linspace(-np.pi/2, np.pi/2, 180)
            lam = np.linspace(0, 2 * np.pi, 360)

            coords_sphere = np.array(np.meshgrid(theta, phi)).T.reshape(-1, 2)
            coords_sphere_jac = np.array(np.meshgrid(beta, lam)).T.reshape(-1, 2)
            eli_coords = self.ellipsoid.ellipsoid2cartesian(coords=coords_sphere, return_normal=False)
            eli_coords_jac = self.ellipsoid.jacobi2cartesian(coords=coords_sphere_jac, return_normal=False)
            # eli_coords_rot = (self.ellipsoid.rotmat.T @ (eli_coords - self.ellipsoid.center).T).T
            # coords_sphere_test = self.ellipsoid.cartesian2ellipsoid(coords=eli_coords)
            # np.isclose(coords_sphere, coords_sphere_test)
            # fig = plt.figure()
            # ax = fig.add_subplot(projection='3d')
            # ax.scatter(eli_coords[:, 0], eli_coords[:, 1], eli_coords[:, 2])
            # ax.scatter(eli_coords_rot[:, 0], eli_coords_rot[:, 1], eli_coords_rot[:, 2])

            np.savetxt(os.path.join(self.output_folder, "plots", "fitted_ellipsoid.txt"), eli_coords)
            np.savetxt(os.path.join(self.output_folder, "plots", "fitted_ellipsoid_jacobian.txt"), eli_coords_jac)
            # np.savetxt(os.path.join("/data/pt_01756/studies/ttf/fitted_ellipsoid.txt"), eli_coords)
            # np.savetxt(os.path.join(self.output_folder, "plots", "fitted_ellipsoid_rot.txt"), eli_coords_rot)

            import pynibs
            # write hdf5 _geo files for visualization in paraview
            pynibs.write_geo_hdf5_surf(out_fn=os.path.join(self.output_folder, "plots", "upper_head_region_geo.hdf5"),
                                       points=self.skin_surface.nodes,
                                       con=self.skin_surface.tr_nodes,
                                       replace=True,
                                       hdf5_path='/mesh')

            pynibs.write_data_hdf5_surf(data=[np.zeros(self.skin_surface.tr_nodes.shape[0])],
                                        data_names=["domain"],
                                        data_hdf_fn_out=os.path.join(self.output_folder, "plots", "upper_head_region_data.hdf5"),
                                        geo_hdf_fn=os.path.join(self.output_folder, "plots", "upper_head_region_geo.hdf5"),
                                        replace=True)

    def valid_skin_region(self, skin_surface, mesh):
        """
        Determine the nodes of the scalp surface where the electrode can be applied (not ears and face etc.)

        Parameters
        ----------
        skin_surface : Surface object
            Surface of the mesh (mesh_tools/surface.py)
        mesh : Msh object
            Mesh object created by SimNIBS (mesh_tools/mesh_io.py)
        """
        # load mask of valid electrode positions (in MNI space)
        mask_img = nib.load(self.fn_electrode_mask)
        mask_img_data = mask_img.get_fdata()

        # transform skin surface points to MNI space
        skin_nodes_mni_ras = subject2mni_coords(coordinates=skin_surface.nodes,
                                                m2m_folder=os.path.split(mesh.fn)[0],
                                                transformation_type='nonl')

        # transform coordinates to voxel space
        skin_nodes_mni_voxel = np.floor(np.linalg.inv(mask_img.affine) @
                                        np.hstack((skin_nodes_mni_ras,
                                                   np.ones(skin_nodes_mni_ras.shape[0])[:, np.newaxis]))
                                        .transpose())[:3, :].transpose().astype(int)
        skin_nodes_mni_voxel[skin_nodes_mni_voxel[:, 0] >= mask_img.shape[0], 0] = mask_img.shape[0] - 1
        skin_nodes_mni_voxel[skin_nodes_mni_voxel[:, 1] >= mask_img.shape[1], 1] = mask_img.shape[1] - 1
        skin_nodes_mni_voxel[skin_nodes_mni_voxel[:, 2] >= mask_img.shape[2], 2] = mask_img.shape[2] - 1
        # skin_nodes_mni_voxel[skin_nodes_mni_voxel < 0] = 0

        # get boolean mask of valid skin points
        skin_surface.mask_valid_nodes = mask_img_data[skin_nodes_mni_voxel[:, 0],
                                                      skin_nodes_mni_voxel[:, 1],
                                                      skin_nodes_mni_voxel[:, 2]].astype(bool)

        # remove points outside of MNI space (lower neck)
        skin_surface.mask_valid_nodes[(skin_nodes_mni_voxel < 0).any(axis=1)] = False

        skin_surface.mask_valid_tr = np.zeros(skin_surface.tr_centers.shape).astype(bool)

        unique_points = np.unique(skin_surface.tr_nodes[skin_surface.mask_valid_nodes[skin_surface.tr_nodes].all(axis=1), :])
        for point in unique_points:
            idx_where = np.where(skin_surface.tr_nodes == point)
            skin_surface.mask_valid_tr[idx_where[0], idx_where[1]] = True
        skin_surface.mask_valid_tr = skin_surface.mask_valid_tr.all(axis=1)

        # determine connectivity list of valid skin region (creates new node and connectivity list)
        skin_surface.nodes, skin_surface.tr_nodes = create_new_connectivity_list_point_mask(
            points=skin_surface.nodes,
            con=skin_surface.tr_nodes,
            point_mask=skin_surface.mask_valid_nodes)

        # identify spurious skin patches inside head and remove them
        tri_domain = np.ones(skin_surface.tr_nodes.shape[0]).astype(int) * -1
        point_domain = np.ones(skin_surface.nodes.shape[0]).astype(int) * -1

        domain = 0
        while (tri_domain == -1).any():
            nodes_idx_of_domain = np.array([])
            tri_idx_of_domain = np.where(tri_domain == -1)[0][0]

            n_current = -1
            n_last = 0
            while n_last != n_current:
                n_last = copy.deepcopy(n_current)
                nodes_idx_of_domain = np.unique(np.append(nodes_idx_of_domain, skin_surface.tr_nodes[tri_idx_of_domain, :])).astype(int)
                tri_idx_of_domain = np.isin(skin_surface.tr_nodes, nodes_idx_of_domain).any(axis=1)
                n_current = np.sum(tri_idx_of_domain)
                # print(f"domain: {domain}, n_current: {n_current}")

            tri_domain[tri_idx_of_domain] = domain
            point_domain[nodes_idx_of_domain] = domain
            domain += 1

        domain_idx_main = np.argmax([np.sum(point_domain == d) for d in range(domain)])

        skin_surface.nodes, skin_surface.tr_nodes = create_new_connectivity_list_point_mask(
            points=skin_surface.nodes,
            con=skin_surface.tr_nodes,
            point_mask=point_domain == domain_idx_main)

        skin_surface.nodes_areas = skin_surface.nodes_areas[skin_surface.mask_valid_nodes]
        skin_surface.nodes_normals = skin_surface.nodes_normals[skin_surface.mask_valid_nodes, :]
        skin_surface.surf2msh_nodes = skin_surface.surf2msh_nodes[skin_surface.mask_valid_nodes]
        skin_surface.surf2msh_triangles = skin_surface.surf2msh_triangles[skin_surface.mask_valid_tr]
        skin_surface.tr_areas = skin_surface.tr_areas[skin_surface.mask_valid_tr]
        skin_surface.tr_centers = skin_surface.tr_centers[skin_surface.mask_valid_tr, :]
        skin_surface.tr_normals = skin_surface.tr_normals[skin_surface.mask_valid_tr, :]

        return skin_surface

    def get_nodes_electrode(self, electrode_pos, plot=False):
        """
        Assigns the skin points of the electrodes in electrode array and writes the points in
        electrode_array.electrodes[i].nodes and electrode_array.electrodes[i].node_area.
        Estimate optimal electrode currents based on previous simulations.

        Parameters
        ----------
        electrode_pos : list of list of np.ndarray of float [3] of length [n_channel_stim][n_ele_free]
            Spherical coordinates (beta, lambda) and orientation angle (alpha) for each electrode array.
                      electrode array 1                        electrode array 2
            [ np.array([beta_1, lambda_1, alpha_1]),   np.array([beta_2, lambda_2, alpha_2]) ]
        plot : bool, optional, default: False
            Generate output files for plotting (debugging)

        Returns
        -------
        node_idx_dict : list of dict
            List [n_channel_stim] containing dicts with electrode channel IDs as keys and node indices.
        """

        electrode_coords_subject = [0 for _ in range(self.n_channel_stim)]

        for i_channel_stim in range(self.n_channel_stim):
            # collect all parameters
            start = np.zeros((self.n_ele_free[i_channel_stim], 3))
            a = np.zeros((self.n_ele_free[i_channel_stim], 3))
            b = np.zeros((self.n_ele_free[i_channel_stim], 3))
            cx = np.zeros((self.n_ele_free[i_channel_stim], 3))
            cy = np.zeros((self.n_ele_free[i_channel_stim], 3))
            n = np.zeros((self.n_ele_free[i_channel_stim], 3))
            start_shifted_ = np.zeros((len(electrode_pos[i_channel_stim]), 3))
            distance = []
            alpha = []

            for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
                start[i_array, :], n[i_array, :] = self.ellipsoid.jacobi2cartesian(
                    coords=electrode_pos[i_channel_stim][i_array][:2],
                    return_normal=True)

                c0, n[i_array, :] = self.ellipsoid.jacobi2cartesian(coords=electrode_pos[i_channel_stim][i_array][:2], return_normal=True)
                a[i_array, :] = self.ellipsoid.jacobi2cartesian(coords=np.array([electrode_pos[i_channel_stim][i_array][0] - 1e-2, electrode_pos[i_channel_stim][i_array][1]])) - c0
                b[i_array, :] = self.ellipsoid.jacobi2cartesian(coords=np.array([electrode_pos[i_channel_stim][i_array][0], electrode_pos[i_channel_stim][i_array][1] - 1e-2])) - c0
                a[i_array, :] /= np.linalg.norm(a[i_array, :])
                b[i_array, :] /= np.linalg.norm(b[i_array, :])

                start_shifted_[i_array, :] = c0 + (1e-3 * ((a[i_array, :]) * np.cos(electrode_pos[i_channel_stim][i_array][2]) +
                                                           (b[i_array, :]) * np.sin(electrode_pos[i_channel_stim][i_array][2])))

                cy[i_array, :] = start_shifted_[i_array, :] - start[i_array, :]
                cy[i_array, :] /= np.linalg.norm(cy[i_array, :])
                cx[i_array, :] = np.cross(cy[i_array, :], -n[i_array, :])
                cx[i_array, :] /= np.linalg.norm(cx[i_array, :])

                distance.append(_electrode_array.distance)
                alpha.append(electrode_pos[i_channel_stim][i_array][2] + _electrode_array.angle)
                # alpha.append(get_array_direction(electrode_pos=electrode_pos[i_array], ellipsoid=self.ellipsoid) +
                #              _electrode_array.angle)

            distance = np.array(distance).flatten()
            alpha = np.array(alpha).flatten()

            for i_a, _alpha in enumerate(alpha):
                if _alpha > np.pi:
                    alpha[i_a] = _alpha - 2 * np.pi
                elif _alpha < -np.pi:
                    alpha[i_a] = _alpha + 2 * np.pi

            start = np.vstack([np.tile(start[i_array, :], (_electrode_array.n_ele, 1))
                               for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays)])

            # determine electrode center on ellipsoid
            electrode_coords_eli_cart = self.ellipsoid.get_geodesic_destination(start=start,
                                                                                distance=distance,
                                                                                alpha=alpha,
                                                                                n_steps=400)

            n = self.ellipsoid.get_normal(coords=electrode_coords_eli_cart)

            # transform to ellipsoidal coordinates
            electrode_coords_eli_eli = self.ellipsoid.cartesian2ellipsoid(coords=electrode_coords_eli_cart)

            # project coordinates to subject
            ele_idx, electrode_coords_subject[i_channel_stim] = ellipsoid2subject(coords=electrode_coords_eli_eli,
                                                                                  ellipsoid=self.ellipsoid,
                                                                                  surface=self.skin_surface)

            if len(ele_idx) != len(alpha):
                # return "Electrode position: invalid (not all electrodes in valid skin region)"
                print("Electrode position: invalid (not all electrodes in valid skin region)")

        # i_ele = 0
        # ele_idx_rect = []
        # start_shifted = []
        #
        # for i_array, _electrode_array in enumerate(self.electrode.electrode_arrays):
        #     for _electrode in _electrode_array.electrodes:
        #         if _electrode.type == "rectangular":
        #             ele_idx_rect.append(i_ele)
        #             start_shifted.append(start_shifted_[i_array])
        #         i_ele += 1

        # loop over electrodes and determine node indices
        node_idx_dict = [dict() for _ in range(self.n_channel_stim)]
        node_coords_list = [[] for _ in range(np.sum(self.n_ele_free))]
        i_array_global = 0

        for i_channel_stim in range(self.n_channel_stim):
            i_ele = 0
            for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
                for _electrode in _electrode_array.electrodes:
                    if _electrode.type == "spherical":
                        # mask with a sphere
                        mask = np.linalg.norm(self.skin_surface.nodes - electrode_coords_subject[i_channel_stim][i_ele, :], axis=1) < _electrode.radius

                        # save position of electrode in subject space to posmat field
                        _electrode.posmat[:3, 3] = electrode_coords_subject[i_channel_stim][i_ele, :]

                    elif _electrode.type == "rectangular":
                        cx_local = np.cross(n[i_ele, :], cy[i_array, :])

                        # rotate skin nodes to normalized electrode space
                        rotmat = np.array([[cx_local[0], cy[i_array, 0], n[i_ele, 0]],
                                           [cx_local[1], cy[i_array, 1], n[i_ele, 1]],
                                           [cx_local[2], cy[i_array, 2], n[i_ele, 2]]])
                        center = np.array([electrode_coords_subject[i_channel_stim][i_ele, 0],
                                           electrode_coords_subject[i_channel_stim][i_ele, 1],
                                           electrode_coords_subject[i_channel_stim][i_ele, 2]])

                        # save position of electrode in subject space to posmat field
                        _electrode.posmat = np.vstack((np.hstack((rotmat, center[:, np.newaxis])), np.array([0, 0, 0, 1])))

                        skin_nodes_rotated = (self.skin_surface.nodes - center) @ rotmat

                        # mask with a box
                        mask_x = np.logical_and(skin_nodes_rotated[:, 0] > -_electrode.length_x / 2,
                                                skin_nodes_rotated[:, 0] < +_electrode.length_x / 2)
                        mask_y = np.logical_and(skin_nodes_rotated[:, 1] > -_electrode.length_y / 2,
                                                skin_nodes_rotated[:, 1] < +_electrode.length_y / 2)
                        mask_z = np.logical_and(skin_nodes_rotated[:, 2] > -20,
                                                skin_nodes_rotated[:, 2] < +20)
                        mask = np.logical_and(np.logical_and(mask_x, mask_y), mask_z)
                    else:
                        raise AssertionError("Electrodes have to be either 'spherical' or 'rectangular'")

                    # node areas
                    _electrode.node_area = self.skin_surface.nodes_areas[mask]

                    # total effective area of all nodes
                    _electrode.area_skin = np.sum(_electrode.node_area)

                    # electrode position is invalid if it overlaps with invalid skin region and area is not "complete"
                    if _electrode.area_skin < 0.8 * _electrode.area:
                        return "Electrode position: invalid (partly overlaps with invalid skin region)"

                    # save node indices (refering to global mesh)
                    _electrode.node_idx = self.node_idx_msh[mask]

                    # save node coords (refering to global mesh)
                    _electrode.node_coords = self.skin_surface.nodes[mask]

                    # save number of nodes assigned to this electrode
                    _electrode.n_nodes = len(_electrode.node_idx)

                    node_coords_list[i_array_global].append(_electrode.node_coords)

                    # group node indices of same channel IDs
                    if _electrode.channel_id in node_idx_dict[i_channel_stim].keys():
                        node_idx_dict[i_channel_stim][_electrode.channel_id] = np.append(node_idx_dict[i_channel_stim][_electrode.channel_id], _electrode.node_idx)
                    else:
                        node_idx_dict[i_channel_stim][_electrode.channel_id] = _electrode.node_idx

                    i_ele += 1

                # gather all electrode node coords of freely movable arrays
                node_coords_list[i_array_global] = np.vstack(node_coords_list[i_array_global])
                i_array_global += 1

        # check if electrode distance is sufficient
        if self.min_electrode_distance is not None and self.min_electrode_distance > 0:
            i_array_test_start = 1
            # start with first array and test if all node coords are too close to other arrays
            for i_array_global in range(np.sum(self.n_ele_free)):
                for node_coord in node_coords_list[i_array_global]:
                    for i_array_test in range(i_array_test_start, np.sum(self.n_ele_free)):
                        # calculate euclidic distance between node coords
                        min_dist = np.min(np.linalg.norm(node_coords_list[i_array_test] - node_coord, axis=1))
                        # stop testing if an electrode is too close
                        if min_dist < self.min_electrode_distance:
                            return "Electrode position: invalid (minimal distance between electrodes too small)"

                i_array_test_start += 1

        if plot:
            # np.savetxt(os.path.join(self.plot_folder, "electrode_coords_center_ellipsoid.txt"), electrode_coords_eli_cart)
            # np.savetxt(os.path.join(self.plot_folder, "electrode_coords_center_subject.txt"), electrode_coords_subject[0])
            # np.savetxt(os.path.join(self.plot_folder, "electrode_coords_nodes_subject.txt"), node_coords_list[0])
            # np.savetxt(os.path.join(self.plot_folder, "electrode_coords_nodes_subject.txt"), self.mesh.nodes.node_coord[self.node_idx_msh])
            # np.savetxt(os.path.join(self.plot_folder, "electrode_coords_nodes_subject.txt"), self.mesh.nodes.node_coord[self.electrode[0].electrode_arrays[0].electrodes[0].node_idx])
            for i_channel_stim in range(self.n_channel_stim):
                node_idx_all = np.hstack([np.hstack(node_idx_dict[i_channel_stim][id])
                                          for id in node_idx_dict[i_channel_stim].keys()])
                points_nodes = self.mesh.nodes.node_coord[node_idx_all, :]
                np.savetxt(os.path.join(self.plot_folder, f"electrode_coords_nodes_subject_{i_channel_stim}.txt"),
                           points_nodes)

        # save electrode_pos in ElectrodeArray instances
        for i_channel_stim in range(self.n_channel_stim):
            for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
                _electrode_array.electrode_pos = electrode_pos[i_channel_stim][i_array]

        # estimate optimal electrode currents based on previous simulations
        # current_pos = 0
        # current_neg = 0
        # for i_channel_stim in range(self.n_channel_stim):
        #     for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
        #         for _electrode in _electrode_array.electrodes:
        #             _electrode.estimate_currents(electrode_pos[i_channel_stim])
        #
        #             # sum up current for scaling
        #             # estimator does not know about real total current and has to be corrected
        #             if _electrode.current < 0:
        #                 current_neg += _electrode.current
        #             else:
        #                 current_pos += _electrode.current
        #
        # scale current to match total current
        # current_temp = []
        # for i_channel_stim in range(self.n_channel_stim):
        #     for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
        #         for _electrode in _electrode_array.electrodes:
        #             if _electrode.current < 0:
        #                 _electrode.current = _electrode.current/np.abs(current_neg) * self.electrode[i_channel_stim].current_total
        #             else:
        #                 _electrode.current = _electrode.current/np.abs(current_pos) * self.electrode[i_channel_stim].current_total
        #             current_temp.append(_electrode.current)

        # self.logger.log(20, f'Estimating currents: { *current_temp, }')

        return node_idx_dict

    def update_field(self, electrode_pos, plot=False):
        """
        Calculate the E field for given electrode positions.

        Parameters
        ----------
        electrode_pos : list of list of np.ndarray of float [3] of length [n_channel_stim][n_ele_free]
            Spherical coordinates (beta, lambda) and orientation angle (alpha) for each electrode array.
                      electrode array 1                        electrode array 2
            [ np.array([beta_1, lambda_1, alpha_1]),   np.array([beta_2, lambda_2, alpha_2]) ]

        Returns
        -------
        e : list of list of np.ndarray [n_channel_stim][n_roi]
            Electric field for different stimulations in ROI(s).
        """

        e = [[0 for _ in range(self.n_roi)] for _ in range(self.n_channel_stim)]

        # assign surface nodes to electrode positions and estimate optimal currents
        # start = time.time()
        node_idx_dict = self.get_nodes_electrode(electrode_pos=electrode_pos, plot=plot)
        # stop = time.time()
        # print(f"Time: get_nodes_electrode: {stop-start}")

        # perform one electric field calculation for every stimulation condition (one at a time is on)
        for i_channel_stim in range(self.n_channel_stim):
            if type(node_idx_dict[i_channel_stim]) is str:
                self.logger.log(20, node_idx_dict)
                return None
            self.logger.log(20, "Electrode position: valid")

            # set RHS
            b = self.ofem.set_rhs(electrode=self.electrode[i_channel_stim])

            # solve system
            if self.dirichlet_correction[i_channel_stim]:
                v = self.ofem.solve_dirichlet_correction(b=b, electrode=self.electrode[i_channel_stim])
            else:
                v = self.ofem.solve(b)

            # Determine electric field in ROIs
            #start = time.time()
            for i_roi, r in enumerate(self.roi):
                if v is None:
                    e[i_channel_stim][i_roi] = None
                    self.logger.log(20, "Warning! Simulation failed! Returning e-field: None!")
                else:
                    e[i_channel_stim][i_roi] = r.calc_fields(v)
            #stop = time.time()
            #print(f"Time: calc fields: {stop - start}")

            # plot field
            if plot:
                for j in range(self.n_channel_stim):
                    for i, _e in enumerate(e[j]):
                        if _e is not None:
                            import pynibs

                            pynibs.write_geo_hdf5_surf(out_fn=os.path.join(self.plot_folder, f"e_roi_{i}_geo.hdf5"),
                                                       points=self.roi[i].points,
                                                       con=self.roi[i].con,
                                                       replace=True,
                                                       hdf5_path='/mesh')

                            pynibs.write_data_hdf5_surf(data=[np.mean(_e.flatten()[self.roi[i].con], axis=1)],
                                                        data_names=["E_mag"],
                                                        data_hdf_fn_out=os.path.join(self.plot_folder, f"e_stim_{j}_roi_{i}_data.hdf5"),
                                                        geo_hdf_fn=os.path.join(self.plot_folder, f"e_roi_{i}_geo.hdf5"),
                                                        replace=True)

        return e

        # determine QOIs in ROIs

        # if normalize:
        #     v_elec = [v[self.array_layout_node_idx[k] - 1] for k in range(len(I))]
        #
        #     v_mean = [np.mean(v_elec[k]) for k in range(len(I))]
        #     for k in range(len(I)):
        #         vn = (v_elec[k] - np.mean(v_elec[k])) / np.std(v_mean)
        #         In = (I[k] - np.mean(I[k])) / np.mean(I[k])
        #
        #         v_norm[k] = np.append(v_norm[k], vn.reshape(1, len(vn)), axis=0)
        #         I_norm[k] = np.append(I_norm[k], In.reshape(1, len(In)), axis=0)

    # def solve(self, node_idx, I, v_norm, I_norm):
    #     """
    #
    #     :param node_idx:
    #     :param I:
    #     :param v_norm:
    #     :param I_norm:
    #     :return:
    #     """
    #     # set RHS (in fem.py, check for speed)
    #     b = self.fem.assemble_tdcs_neumann_rhs(electrodes=node_idx, currents=TODO, input_type='node', areas=TODO)
    #     # solve
    #     v = self.fem.solve(b)
    #
    #     v_elec = [v[self.array_layout_node_idx[k] - 1] for k in range(len(I))]
    #
    #     v_mean = [np.mean(v_elec[k]) for k in range(len(I))]
    #     for k in range(len(I)):
    #         vn = (v_elec[k] - np.mean(v_elec[k])) / np.std(v_mean)
    #         In = (I[k] - np.mean(I[k])) / np.mean(I[k])
    #
    #         v_norm[k] = np.append(v_norm[k], vn.reshape(1, len(vn)), axis=0)
    #         I_norm[k] = np.append(I_norm[k], In.reshape(1, len(In)), axis=0)
    #
    #     return v, v_elec, v_norm, I_norm

    def run(self, parameters):
        """
        Run function for optimization algorithms.

        Parameters
        ----------
        parameters : np.ndarray of float [n_channel_stim * n_free_arrays * 3]
            Electrodes positions in spherical coordinates (theta, phi, alpha) for each freely movable array.
            e.g.: np.array([theta_stim_1_1, phi_stim_1_1, alpha_stim_1_1, theta_stim_1_2, phi_stim_1_2, alpha_2, ...])

        Returns
        -------
        y : float
            Goal function value
        """
        parameters_str = f"Parameters: {parameters}"
        self.logger.log(20, parameters_str)

        # reformat parameters
        parameters_array = np.reshape(parameters, (np.sum(self.n_ele_free), 3))
        self.electrode_pos = [[] for _ in range(self.n_channel_stim)]

        i_para = 0
        for i_channel_stim in range(self.n_channel_stim):
            for i_ele_free in range(self.n_ele_free[i_channel_stim]):
                self.electrode_pos[i_channel_stim].append(parameters_array[i_para, :])
                i_para += 1

        # update field, returns list of list e[n_channel_stim][n_roi] (None if position is not applicable)
        e = self.update_field(electrode_pos=self.electrode_pos, plot=False)

        # calculate goal function value for every ROI
        y = np.zeros((self.n_channel_stim, self.n_roi))     # shape: [n_channel_stim x n_roi]

        if e is None:
            y = np.ones(self.n_roi)
        else:
            if self.goal == "mean":
                for i_channel_stim in range(self.n_channel_stim):
                    for i_roi in range(self.n_roi):
                        if e[i_channel_stim] is None:
                            y[i_channel_stim, i_roi] = 1
                        else:
                            y[i_channel_stim, i_roi] = -np.mean(e[i_channel_stim][i_roi])
            elif self.goal == "max":
                for i_channel_stim in range(self.n_channel_stim):
                    for i_roi in range(self.n_roi):
                        if e[i_channel_stim] is None:
                            y[i_channel_stim, i_roi] = 1
                        else:
                            y[i_channel_stim, i_roi] = -np.percentile(e[i_channel_stim][i_roi], 99.9)
            else:
                raise NotImplementedError(f"Specified goal: '{self.goal}' not implemented as goal function.")

        # weight and sum the goal function values of the ROIs
        y_weighted_sum = np.sum(y * self.weights)

        self.logger.log(20, f"Goal ({self.goal}): {y_weighted_sum:.3f}")
        self.logger.log(20, "-"*len(parameters_str))

        return y_weighted_sum

    def optimize(self):
        """
        Runs the optimization problem

        Returns
        -------
        parameters : list of np.ndarray [3]
            List containing the optimal parameters for each freely movable electrode array
            [np.array([beta_1, lambda_1, alpha_1]), np.array([beta_2, lambda_2, alpha_2]), ...]
        """
        # run optimization
        ################################################################################################################
        start = time.time()
        if self.optimizer == "direct":
            result = direct(self.run,
                            bounds=self.optimizer_options["bounds"],
                            vol_tol=self.optimizer_options["vol_tol"],
                            len_tol=self.optimizer_options["len_tol"],
                            f_min_rtol=self.optimizer_options["f_min_rtol"],
                            maxiter=self.optimizer_options["maxiter"])

            # reformat parameters
            parameters_array = np.reshape(result.x, (np.sum(self.n_ele_free), 3))
            self.electrode_pos_opt = [[] for _ in range(self.n_channel_stim)]

            i_para = 0
            for i_channel_stim in range(self.n_channel_stim):
                for i_ele_free in range(self.n_ele_free[i_channel_stim]):
                    self.electrode_pos_opt[i_channel_stim].append(parameters_array[i_para, :])
                    i_para += 1

            fopt = result.fun
            nfev = result.nfev
            self.logger.log(20, f"Optimization finished! Best electrode position: {self.electrode_pos_opt}")
        else:
            raise NotImplementedError(f"Specified optimization method: '{self.optimizer}' not implemented.")

        stop = time.time()
        t_optimize = stop - start

        # plot final solution and electrode position
        ################################################################################################################
        # compute best field again, plot field and electrode position
        e = self.update_field(electrode_pos=self.electrode_pos_opt, plot=True)

        # print optimization summary
        save_optimization_results(fname=os.path.join(self.output_folder, "summary"),
                                  optimizer=self.optimizer,
                                  optimizer_options=self.optimizer_options,
                                  fopt=fopt,
                                  popt=self.electrode_pos_opt,
                                  nfev=nfev,
                                  e=e,
                                  time=t_optimize,
                                  msh=self.mesh,
                                  electrode=self.electrode,
                                  goal=self.goal)

        # plot skin surface
        import pynibs
        pynibs.write_geo_hdf5_surf(out_fn=os.path.join(self.plot_folder, f"skin_surface_geo.hdf5"),
                                   points=self.skin_surface.nodes,
                                   con=self.skin_surface.tr_nodes,
                                   replace=True,
                                   hdf5_path='/mesh')

        pynibs.write_data_hdf5_surf(data=[np.zeros(self.skin_surface.tr_nodes.shape[0])],
                                    data_names=["domain"],
                                    data_hdf_fn_out=os.path.join(self.plot_folder, f"skin_surface_data.hdf5"),
                                    geo_hdf_fn=os.path.join(self.plot_folder, f"skin_surface_geo.hdf5"),
                                    replace=True)

        # save fitted ellipsoid
        beta = np.linspace(-np.pi / 2, np.pi / 2, 180)
        lam = np.linspace(0, 2 * np.pi, 360)
        coords_sphere_jac = np.array(np.meshgrid(beta, lam)).T.reshape(-1, 2)
        eli_coords_jac = self.ellipsoid.jacobi2cartesian(coords=coords_sphere_jac, return_normal=False)
        np.savetxt(os.path.join(self.output_folder, "plots", "fitted_ellipsoid.txt"), eli_coords_jac)


def save_optimization_results(fname, optimizer, optimizer_options, fopt, popt, nfev, e, time, msh, electrode, goal):
    """
    Saves optimization settings and results in an <fname>.hdf5 file and prints a summary in a <fname>.txt file.

    Parameters
    ----------
    fname : str
        Filename of the summary.txt file.
    optimizer : str
        Name of optimization method.
    optimizer_options : dict
        Dictionary containing the optimization setting.
    fopt : float
        Objective function value in optimum.
    popt : list of list of np.ndarray of float [n_channel_stim][n_free_electrodes]
        List containing the optimal parameters for each freely movable electrode array
        [np.array([beta_1, lambda_1, alpha_1]), np.array([beta_2, lambda_2, alpha_2]), ...]
    nfev : int
        Number of function evaluations during optimization
    e : list of np.ndarray [n_roi]
        List of containing np.ndarrays of the electric fields in the ROIs
    time : float
        Runtime of optimization in s
    electrode : list of ElectrodeArray objects [n_channel_stim]
        List of ElectrodeArray objects for every stimulation
    goal : str

    """

    def sep(x):
        if x >= 0:
            sep = " "
        else:
            sep = ""

        return sep

    fname_txt = fname + ".txt"
    fname_hdf5 = fname + ".hdf5"

    # print summary in <fname>.txt file
    ####################################################################################################################
    d = datetime.datetime.now()

    with open(fname_txt, 'w') as f:
        f.write(f"Optimization summary:\n")
        f.write(f"===================================================================\n")
        f.write(f"Date: {d.year}-{d.month}-{d.day}, {d.hour}:{d.minute}:{d.second}\n")
        f.write(f"Simulation time: {str(datetime.timedelta(seconds=time))[:-7]}\n")
        f.write(f"headmodel: {msh.fn}\n")
        f.write(f"Goal: {goal}\n")
        f.write(f"Number of Channels: {len(e)}\n")
        f.write(f"Number of ROIs: {len(e[0])}\n")
        f.write(f"\n")
        f.write(f"Electrode coordinates:\n")
        f.write(f"===================================================================\n")
        f.write(f"Ellipsoid space (Jacobian coordinates):\n")
        f.write(f"---------------------------------------\n")

        for i_stim in range(len(popt)):
            f.write(f"Stimulation {i_stim}:\n")
            for i, p in enumerate(popt[i_stim]):
                f.write(f"Array {i}:\n")
                f.write(f"\tbeta:   {sep(p[0])}{p[0]:.3f}\n")
                f.write(f"\tlambda: {sep(p[1])}{p[1]:.3f}\n")
                f.write(f"\talpha:  {sep(p[2])}{p[2]:.3f}\n")
        f.write(f"\n")
        f.write(f"Subject space (Cartesian coordinates):\n")
        f.write(f"--------------------------------------\n")
        for i_stim in range(len(popt)):
            f.write(f"Stimulation {i_stim}:\n")
            for i_array, _electrode_array in enumerate(electrode[i_stim].electrode_arrays):
                f.write(f"Array {i_array}:\n")
                for i_electrode, _electrode in enumerate(_electrode_array.electrodes):
                    f.write(f"\tElectrode {i_electrode} ({_electrode.type}):\n")
                    for i_row in range(4):
                        f.write("\t\t" + sep(_electrode.posmat[i_row, 0]) + f"{_electrode.posmat[i_row, 0]:.3f}, " +
                                sep(_electrode.posmat[i_row, 1]) + f"{_electrode.posmat[i_row, 1]:.3f}, " +
                                sep(_electrode.posmat[i_row, 2]) + f"{_electrode.posmat[i_row, 2]:.3f}, " +
                                sep(_electrode.posmat[i_row, 3]) + f"{_electrode.posmat[i_row, 3]:.3f}\n")
        f.write("\n")
        f.write("Optimization method:\n")
        f.write("===================================================================\n")
        f.write(f"Optimizer: {optimizer}\n")
        f.write(f"Settings:\n")
        if optimizer_options is not None:
            for key in optimizer_options:
                if type(optimizer_options[key]) is Bounds:
                    f.write(f"\tlb: {optimizer_options[key].lb}\n")
                    f.write(f"\tub: {optimizer_options[key].ub}\n")
                else:
                    f.write(f"\t{key}: {optimizer_options[key]}\n")
        else:
            f.write("None\n")
        f.write(f"nfev: {nfev}\n")
        f.write(f"fopt: {fopt}\n")

    # save optimization settings and results in <fname>.hdf5 file
    ####################################################################################################################
    with h5py.File(fname_hdf5, "w") as f:
        # general info
        f.create_dataset(data=msh.fn, name="fnamehead")
        f.create_dataset(data=f"{d.year}-{d.month}-{d.day}, {d.hour}:{d.minute}:{d.second}", name="date")
        f.create_dataset(data=len(e), name="n_channel")
        f.create_dataset(data=len(e[0]), name="n_roi")

        # optimizer
        f.create_dataset(data=optimizer, name="optimizer/optimizer")
        f.create_dataset(data=fopt, name="optimizer/fopt")
        f.create_dataset(data=nfev, name="optimizer/nfev")
        f.create_dataset(data=time, name="optimizer/time")
        f.create_dataset(data=goal, name="optimizer/goal")

        for key in optimizer_options:
            if type(optimizer_options[key]) is Bounds:
                f.create_dataset(data=optimizer_options[key].lb, name=f"optimizer/optimizer_options/lb")
                f.create_dataset(data=optimizer_options[key].ub, name=f"optimizer/optimizer_options/ub")
            else:
                f.create_dataset(data=optimizer_options[key], name=f"optimizer/optimizer_options/{key}")

        # electrodes
        for i_stim in range(len(electrode)):
            f.create_dataset(data=electrode[i_stim].center, name=f"electrode/channel_{i_stim}/center")
            f.create_dataset(data=electrode[i_stim].radius, name=f"electrode/channel_{i_stim}/radius")
            f.create_dataset(data=electrode[i_stim].length_x, name=f"electrode/channel_{i_stim}/length_x")
            f.create_dataset(data=electrode[i_stim].length_y, name=f"electrode/channel_{i_stim}/length_y")
            f.create_dataset(data=electrode[i_stim].current, name=f"electrode/channel_{i_stim}/current")

            for i_array, _electrode_array in enumerate(electrode[i_stim].electrode_arrays):
                f.create_dataset(data=popt[i_stim][i_array][0], name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/beta")
                f.create_dataset(data=popt[i_stim][i_array][1], name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/lambda")
                f.create_dataset(data=popt[i_stim][i_array][2], name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/alpha")

                for i_electrode, _electrode in enumerate(_electrode_array.electrodes):
                    f.create_dataset(data=_electrode.posmat, name=f"electrode/channel_{i_stim}/posmat/electrode_array_{i_array}/electrode_{i_electrode}/posmat")

        # electric field in ROIs
        for i_stim in range(len(e)):
            for i_roi in range(len(e[i_stim])):
                f.create_dataset(data=e[i_stim][i_roi].flatten(), name=f"e/channel_{i_stim}/e_roi_{i_roi}")


def relabel_internal_air(m, subpath, label_skin=1005, label_new=1099, label_internal_air=501):
    """
    Relabels skin in internal air cavities to something else; relevant for charm meshes

    Parameters
    ----------
    m : Msh object
        Mesh object with internal air
    subpath :
        Path to subject m2m folder
    label_skin : int
        Original skin label
    label_new : int
        New skin label
    label_internal_air : int
        New label of internal air

    Returns
    -------
    m : Msh object
        Mesh with relabeled internal air
    """
    subject_files = SubjectFiles(subpath=subpath)

    # relabel internal skin to some other label
    label_nifti = nib.load(subject_files.labeling)
    label_affine = label_nifti.affine
    label_img = label_nifti.get_fdata().astype(int)
    label_img = label_img == label_internal_air
    label_img = mrph.binary_dilation(label_img, iterations=2)

    m = copy.copy(m)
    ed = mesh_io.ElementData.from_data_grid(m, label_img, label_affine, order=0)
    idx = ed.value * (m.elm.tag1 == label_skin)
    m.elm.tag1[idx] = label_new
    m.elm.tag2[:] = m.elm.tag1

    return m


def get_array_direction(electrode_pos, ellipsoid):
    """
    Determine electrode array direction given in ellipsoidal coordinates [theta, phi, alpha] and return direction
    angle alpha in Jacobian coordinates.

    Parameters
    ----------
    electrode_pos : np.ndarray of float [3]
        Spherical coordinates (theta, phi) and orientation angle (alpha) for electrode array (in this order)
    ellipsoid : Ellipsoid class instance
        Ellipsoid

    Returns
    -------
    alpha_jac : float
        Angle of array with respect to constant lambda (Jacobian coordinates)
    """
    # create cartesian vector in direction of electrode array in ellipsoidal coordinates (c1 -> c4)
    c1 = ellipsoid.ellipsoid2cartesian(coords=electrode_pos[:2])
    c2 = ellipsoid.ellipsoid2cartesian(coords=np.array([electrode_pos[0] - 1e-3, electrode_pos[1]]))
    c3 = ellipsoid.ellipsoid2cartesian(coords=np.array([electrode_pos[0], electrode_pos[1] - 1e-3]))
    c4 = c1 + (1e-3 * ((c3 - c1) * np.sin(electrode_pos[2]) + (c2 - c1) * np.cos(electrode_pos[2])))

    # create vector in direction of constant lambda
    l1 = ellipsoid.cartesian2jacobi(coords=c1)
    l2 = ellipsoid.jacobi2cartesian(coords=np.array([l1[0, 0] + 1e-3, l1[0, 1]]))

    l2c1 = ((l2 - c1) / np.linalg.norm(l2 - c1)).flatten()
    c4c1 = ((c4 - c1) / np.linalg.norm(c4 - c1)).flatten()
    angle_jac = np.arccos(np.dot((c4c1), (l2c1)))

    return angle_jac


def create_new_connectivity_list_point_mask(points, con, point_mask):
    """
    Creates a new point and connectivity list when applying a point mask (changes indices of points)

    Parameters
    ----------
    points : np.ndarray of float [n_points x 3]
        Point coordinates
    con : np.ndarray of float [n_tri x 3]
        Connectivity of triangles
    point_mask : nparray of bool [n_points]
        Mask of (True/False) which points are kept in the mesh

    Returns
    -------
    points_new : np.ndarray of float [n_points_new x 3]
        New point array containing the remaining points after applying the mask
    con_new : np.ndarray of float [n_tri_new x 3]
        New connectivity list containing the remaining points (includes reindexing)
    """
    con_global = con[point_mask[con].all(axis=1), :]
    unique_points = np.unique(con_global)
    points_new = points[unique_points, :]

    con_new = np.zeros(con_global.shape).astype(int)

    for i, idx in enumerate(unique_points):
        idx_where = np.where(con_global == idx)
        con_new[idx_where[0], idx_where[1]] = i

    return points_new, con_new


def setup_logger(logname, filemode='w', format='[ %(name)s ] %(levelname)s: %(message)s',
                 datefmt='%H:%M:%S'):
    """
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
                        level=logging.DEBUG)

    logger = logging.getLogger("simnibs")

    return logger

# def get_element_intersect_line_surface(p, w, points, con, triangle_center=None, triangle_normals=None):
#     """
#     Get a element indices of where a given line (point p, direction w) intersects with a surface
#     given by points (points) and connectivity list (con).
#
#     Parameters
#     ----------
#     p : np.ndarray of float [3]
#         Base point of ray
#     w : np.ndarray of float [3]
#         Direction of ray
#     points : np.ndarray of float [n_points x 3]
#         Surface points
#     con : np.ndarray of int [n_tri x 3]
#         Connectivity list of triangles
#     triangle_center : np.ndarray of float [n_tri x 3]
#         Center of triangles
#     triangle_normals : np.ndarray of float [n_tri x 3]
#         Normals of triangles
#
#     Returns
#     -------
#     ele_idx : np.ndarray of int [n_tri_intersect]
#         Index of triangles where ray intersects surface
#     dist : np.ndarray of float [n_tri_intersect]
#         Distances between source point p and intersection
#     """
#
#     p = np.tile(p, (con.shape[0], 1))
#     w = np.tile(w, (con.shape[0], 1))
#
#     p1 = points[con[:, 0], :]
#     p2 = points[con[:, 1], :]
#     p3 = points[con[:, 2], :]
#
#     if triangle_center is None:
#         triangle_center = 1 / 3. * (p1 + p2 + p3)
#
#     if triangle_normals is None:
#         triangle_normals = np.cross(p2 - p1, p3 - p1)
#
#     pi = p + w * (-np.sum((p - triangle_center) * triangle_normals, axis=1) / np.sum(w * triangle_normals, axis=1))[:, np.newaxis]
#
#     mask = np.ones(con.shape[0]).astype(bool)
#
#     p1p2 = p1 - p2
#     p3p2 = p3 - p2
#     pp2 = pi - p2
#     p1p3 = p1 - p3
#     p2p3 = p2 - p3
#     pp3 = pi - p3
#
#     l1 = np.sum(p1p2 * pp2, axis=1)
#     l2 = np.sum(p3p2 * pp2, axis=1)
#     l3 = np.sum(p1p3 * pp3, axis=1)
#     l4 = np.sum(p2p3 * pp3, axis=1)
#
#     mask *= (l1 >= 0) * (l2 >= 0) * (l3 >= 0) * (l4 >= 0)
#     ele_idx = np.where(mask)[0]
#     dist = np.linalg.norm(pi[mask] - triangle_center[mask], axis=1)
#
#     return ele_idx, dist

    # def flatten_surface(self):
    #     """
    #     Flatten surface
    #     :return:
    #     """
    #     import flatsurf.halfedge as halfedge
    #     import flatsurf.mccartney1999 as mc1999
    #
    #     # create surface from valid skin region points
    #
    #     triangulation = halfedge.HalfEdgeTriangulation.from_coords_and_simplices(self.skin_surface.nodes_valid, self.skin_surface.tr_nodes_valid.tolist())
    #     flat_trig = mc1999.McCartney1999Flattening(triangulation, Et=1e-4, delta=1e-5)
    #     flat_trig.flatten()
    #     coords_2d = [np.concatenate([vertex.coord_2d, [0.0]]) for vertex in triangulation.vertices]
    #
    #     import pynibs
    #     # write hdf5 _geo file
    #     pynibs.write_geo_hdf5_surf(out_fn="/home/kporzig/tmp/test_geo_head.hdf5",
    #                                points=self.skin_surface.nodes_valid,
    #                                con=self.skin_surface.tr_nodes_valid,
    #                                replace=True,
    #                                hdf5_path='/mesh')
    #
    #     pynibs.write_data_hdf5_surf(data=[np.zeros(self.skin_surface.tr_nodes_valid.shape[0])],
    #                                 data_names=["test"],
    #                                 data_hdf_fn_out="/home/kporzig/tmp/test_data_head.hdf5",
    #                                 geo_hdf_fn="/home/kporzig/tmp/test_geo_head.hdf5",
    #                                 replace=True)




# theta = np.pi / 2# np.linspace(0, np.pi, 20)
#         phi = np.pi / 2 # np.linspace(0, 2*np.pi, 20)
#         # coords_sphere = np.hstack((theta[:, np.newaxis], phi[:, np.newaxis]))
#         coords_sphere = np.array([theta, phi])[np.newaxis, :]
#         index, point = self.ellipsoid2subject(coords_sphere)
#
#         import pynibs
#         pynibs.write_geo_hdf5_surf(out_fn="/home/kporzig/tmp/test_geo_skin_surface_valid_project.hdf5",
#                                    points=self.skin_surface.nodes_valid,
#                                    con=self.skin_surface.tr_nodes_valid,
#                                    replace=True,
#                                    hdf5_path='/mesh')
#
#         data = np.zeros(self.skin_surface.tr_nodes_valid.shape[0])
#         data[index] = 1
#
#         pynibs.write_data_hdf5_surf(data=[data],
#                                     data_names=["intersect"],
#                                     data_hdf_fn_out="/home/kporzig/tmp/test_data_skin_surface_valid_project.hdf5",
#                                     geo_hdf_fn="/home/kporzig/tmp/test_geo_skin_surface_valid_project.hdf5",
#                                     replace=True)


# import pynibs
# # write hdf5 _geo files for visualization in paraview
# pynibs.write_geo_hdf5_surf(out_fn="/data/pt_01756/studies/ttf/upper_head_region_geo.hdf5",
#                            points=self.skin_surface.nodes,
#                            con=self.skin_surface.tr_nodes,
#                            replace=True,
#                            hdf5_path='/mesh')
#
# pynibs.write_data_hdf5_surf(data=[np.zeros(self.skin_surface.tr_nodes.shape[0])],
#                             data_names=["domain"],
#                             data_hdf_fn_out="/data/pt_01756/studies/ttf/upper_head_region_data.hdf5",
#                             geo_hdf_fn=os.path.join(self.output_folder, "plots", "upper_head_region_geo.hdf5"),
#                             replace=True)
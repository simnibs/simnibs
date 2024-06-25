import copy
import csv
import os
import time
import h5py
import types
import datetime
import logging
import sys

import numpy as np
import scipy.spatial
import nibabel
from scipy.optimize import (
    direct,
    Bounds,
    minimize,
    differential_evolution,
    shgo,
    basinhopping,
)

from ...simulation.sim_struct import ELECTRODE
from ...simulation.fem import get_dirichlet_node_index_cog
from ...utils.region_of_interest import RegionOfInterest
from .electrode_layout import (
    ElectrodeArray,
    ElectrodeArrayPair,
    ElectrodeLayout,
    create_tdcs_session_from_array,
    CircularArray,
)
from ...simulation.onlinefem import FemTargetPointCloud, OnlineFEM, postprocess_e
from ...mesh_tools import mesh_io, surface
from ...utils.file_finder import SubjectFiles, Templates
from ...utils.transformations import (
    subject2mni_coords,
    create_new_connectivity_list_point_mask,
)
from .ellipsoid import Ellipsoid, subject2ellipsoid, ellipsoid2subject
from ...utils.TI_utils import get_maxTI, get_dirTI
from .measures import AUC, integral_focality, ROC
from simnibs.utils import simnibs_logger
from simnibs.utils.simnibs_logger import logger
from simnibs import __version__

class TesFlexOptimization:
    """
    Defines a TES optimization problem using node-wise current sources.

    Parameters
    --------------
    electrode : Electrode Object
        Electrode object containing ElectrodeArray instances
        (see /simulation/array_layout.py for pre-implemented examples)
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
    disable_SPR_for_volume_roi : bool
        Weather to use SPR interpolation for volume rois, defaul: True
    anisotropy_type : str
        Specify type of anisotropy for simulation ('scalar', 'vn' or 'mc')
    weights : np.array of float [n_roi]
        Weights for optimizer for ROI specific goal function weighting
    min_electrode_distance : float, optional, default: None
        Minimum electrode distance to ensure during optimization (in mm).
    constrain_electrode_locations : bool, optional, default: False
        Constrains the possible locations of freely movable electrode arrays. Recommended for TTF optimizations,
        where two pairs of large electrode arrays are optimized. If True, parameter bounds for the optimization
        will be specified restricting the array locations to be frontal, parietal or occipital.
    overlap_factor : float, optional, default: 1
        Factor of overlap of allowed lambda regions to place electrodes. (1 corresponds to neatless regions,
        <1 regions have a gap in between, >1 regions are overlapping)
    plot : bool, optional, default: False
        Plot configurations in output folder for visualization and control
    polish : bool, optional, default: True
        If True (default), then scipy.optimize.minimize with the L-BFGS-B method is used to polish the best
        population member at the end, which can improve the minimization.
    optimize_init_vals : bool, optional, default: True
        If True, find initial values for optimization, guaranteeing a valid solution. If False, initial values
        are the center between bounds.
    e_postproc : str, optional, default: "norm"
        Specifies how the raw e-field in the ROI (Ex, Ey, Ez) is post-processed.
        - "norm": electric field magnitude (default)
        - "normal": determine normal component (required surface normals in dirvec)
        - "tangential": determine tangential component (required surface normals in dirvec)
        - "max_TI": maximum envelope for temporal interference fields
        - "dir_TI": directional sensitive maximum envelope for temporal interference fields
    optimizer : str, optional, default: "differential_evolution"
        Optimization algorithm
    goal : list of str [n_roi], or FunctionType, optional, default: ["mean"]
        Implemented or user provided goal functions:
        - "mean": maximize mean e-field in ROI
        - "max": maximize 99.9 percentile of electric field in ROI
        - "focality": Maximize focality  (goal: sensitivity = specificity = 1)
        - "focality_inv": Maximize inverse focality (goal: sensitivity(ROI) = 1, sensitivity(nonROI) = 1)
        - user provided function taking e-field as an input which is  a list of list of np.ndarrays of float [n_channel_stim][n_roi] containing np.array with e-field
    track_focality : bool, optional, default: False
        Tracks focality for each goal function value (requires ROI and non-ROI definition)
    run_final_electrode_simulation : bool, optional, default: True
        Runs final simulation with optimized parameters using real electrode model including remeshing.

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

    def __init__(self, settings_dict=None):
        """Initialized TESoptimize class instance"""
        # folders and I/O
        self.date = time.strftime("%Y-%m-%d %H:%M:%S")
        self.time_str = time.strftime("%Y%m%d-%H%M%S")
        self.output_folder = None
        self.plot_folder = None
        self.plot = False
        self.fn_final_sim = []
        self.fn_results_hdf5 = None
        self._prepared = False

        # headmodel
        self.fn_eeg_cap = "EEG10-20_Okamoto_2004.csv"
        self.fn_mesh = None
        self.subpath = None
        self._mesh = None
        self._mesh_relabel = None
        self._mesh_nodes_areas = None
        self._ff_templates = Templates()
        self._ff_subject = None
        self._skin_surface = None
        self._node_idx_msh = None
        self._ellipsoid = Ellipsoid()
        self._fn_electrode_mask = self._ff_templates.mni_volume_upper_head_mask

        # roi
        self.roi : list[RegionOfInterest] = []
        self._n_roi = None
        self._con = []
        self._nodes = []
        self._vol = []
        self.disable_SPR_for_volume_roi = True

        # electrode
        self.electrode: list[ElectrodeArray] = []
        self.electrode_pos = None
        self.electrode_pos_opt = None
        self.min_electrode_distance = 1e-3
        self.n_channel_stim = None
        self.n_iter_dirichlet_correction = None
        self.n_ele_free = None
        self.init_pos = None
        self.init_pos_subject_coords = None

        # goal function
        self.goal = None
        self._goal_dir = None
        self.e_postproc = None
        self.threshold = None
        self.optimizer = "differential_evolution"
        self.weights = None
        self.track_focality = False
        self.constrain_electrode_locations = False
        self.overlap_factor = None
        self.polish = False
        self.n_test = 0  # number of tries to place the electrodes
        self.n_sim = 0  # number of final simulations carried out (only valid electrode positions)
        self.optimize_init_vals = True
        self._bounds = None
        self.x0 = None

        # track goal fun value (in ROI 0) and focality measures for later analysis
        self.goal_fun_value = None
        self.AUC = None
        self.integral_focality = None

        # set default options for optimizer
        self.optimizer_options = None  # passed by user
        self._optimizer_options_std = {
            "len_tol": 1.0 / 3600000000.0,
            "f_min_rtol": 1e-12,
            "maxiter": 1000,
            "polish": True,
            "disp": True,
            "recombination": 0.7,  # differential evolution
            "mutation": [0.01, 0.5],  # differential evolution
            "popsize": 13,  # differential evolution
            "tol": 0.1,  # differential evolution
            "locally_biased": False,
        }

        # FEM
        self.run_final_electrode_simulation = True
        self.dirichlet_node = None
        self.dataType = None
        self.anisotropy_type = "scalar"
        self.solver_options = "pardiso"
        self._ofem = None

        # number of CPU cores (set by run(cpus=NN) )
        self._n_cpu = None
        self._log_handlers = []
        
        if settings_dict:
            self.from_dict(settings_dict)

    def _prepare(self):
        """
        Prepares TESoptimize
        
        """
        # check variable assignments
        ################################################################################################################
        if self.output_folder is None:
            raise ValueError("Please define TESoptimize.output_folder !")

        if self.subpath is None:
            raise ValueError(
                "Please define TESoptime.subpath: m2m_* folder containing the headmodel."
            )

        if self.roi is None:
            raise ValueError(
                "Please define TESoptime.roi using the simnibs.RegionOfInterest class."
            )

        if self.electrode is None:
            raise ValueError(
                "Please define TESoptime.electrode using the simnibs.ElectrodeArrayPair or "
                "simnibs.CircularArray class."
            )

        if self.goal is None:
            raise ValueError(
                "Please define type of goal function in TESoptimize.goal"
                " ('mean', 'focality', 'focality_inv')"
            )

        # setup output folders, logging and IO
        ################################################################################################################
        self.plot_folder = os.path.join(self.output_folder, "plots")
        self.fn_results_hdf5 = os.path.join(self.output_folder, "opt.hdf5")

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        if not os.path.exists(self.plot_folder):
            os.makedirs(self.plot_folder)

        logger.info( "Setting up output folders, logging and IO ...")

        # setup headmodel
        ################################################################################################################
        logger.info( "Setting up headmodel ...")

        # get subject specific filenames
        self._ff_subject = SubjectFiles(subpath=self.subpath)

        # read mesh or store in self
        self.fn_mesh = self._ff_subject.fnamehead
        self._mesh = mesh_io.read_msh(self.fn_mesh)

        # Calculate node areas for whole mesh
        self._mesh_nodes_areas = self._mesh.nodes_areas()

        # relabel internal air
        self._mesh_relabel = self._mesh.relabel_internal_air()

        # make final skin surface including some additional distance
        self._skin_surface = surface.Surface(mesh=self._mesh_relabel, labels=1005)
        self._skin_surface = valid_skin_region(
            skin_surface=self._skin_surface,
            fn_electrode_mask=self._fn_electrode_mask,
            mesh=self._mesh_relabel,
            additional_distance=0,
        )

        # get mapping between skin_surface node indices and global mesh nodes
        self._node_idx_msh = np.where(
            np.isin(self._mesh.nodes.node_coord, self._skin_surface.nodes).all(axis=1)
        )[0]

        # fit optimal ellipsoid to valid skin points
        self._ellipsoid.fit(points=self._skin_surface.nodes)

        # plot skin surface and ellipsoid
        if self.plot:
            # save skin surface
            np.savetxt(
                os.path.join(self.plot_folder, "skin_surface_nodes.txt"),
                self._skin_surface.nodes,
            )
            np.savetxt(
                os.path.join(self.plot_folder, "skin_surface_con.txt"),
                self._skin_surface.tr_nodes,
            )

            # save fitted ellipsoid
            beta = np.linspace(-np.pi / 2, np.pi / 2, 180)
            lam = np.linspace(0, 2 * np.pi, 360)
            coords_sphere_jac = np.array(np.meshgrid(beta, lam)).T.reshape(-1, 2)
            eli_coords_jac = self._ellipsoid.jacobi2cartesian(
                coords=coords_sphere_jac, return_normal=False
            )
            np.savetxt(
                os.path.join(self.plot_folder, "fitted_ellipsoid.txt"), eli_coords_jac
            )

        # setup ROI
        ################################################################################################################
        logger.info( "Setting up ROI ...")

        if type(self.roi) is not list:
            self.roi = [self.roi]

        if type(self._con) is not list:
            self._con = [self._con]

        if type(self._nodes) is not list:
            self._nodes = [self._nodes]

        # initialize ROIs if not done already
        self._roi = []
        for i in range(len(self.roi)):
            if self.subpath is not None and self.roi[i].subpath is None and self.roi[i].mesh is None:
                    self.roi[i].subpath = self.subpath
            elif self.roi[i].subpath is None and self.roi.mesh is None:
                self.roi[i].mesh = self._mesh
            self._roi.append(FemTargetPointCloud(self._mesh, self.roi[i].get_nodes(), nearest_neighbor=((self.roi[i].method == "volume" or self.roi[i].method == "volume_from_surface") and self.disable_SPR_for_volume_roi)))

        self._n_roi = len(self._roi)

        # setup electrode
        ################################################################################################################
        logger.info( "Setting up electrodes ...")

        if type(self.electrode) is not list:
            self.electrode = [self.electrode]

        # initialize electrodes if not done already
        for i in range(len(self.electrode)):
            self.electrode[i]._prepare()

        # number of independent stimulation channels
        self.n_channel_stim = len(self.electrode)

        # initialize lists with number of dirichlet correction iterations for convergence analysis
        self.n_iter_dirichlet_correction = [[] for _ in range(self.n_channel_stim)]

        # list containing the number of freely movable arrays for each channel [i_channel_stim]
        self.n_ele_free = [len(ele._electrode_arrays) for ele in self.electrode]

        # list containing beta, lambda, alpha for each freely movable array and for each stimulation channel
        self.electrode_pos = [
            [0 for _ in range(n_ele_free)] for n_ele_free in self.n_ele_free
        ]

        for i_channel_stim in range(self.n_channel_stim):
            for i_array, _electrode_array in enumerate(
                self.electrode[i_channel_stim]._electrode_arrays
            ):
                if _electrode_array.optimize_alpha:
                    self.electrode_pos[i_channel_stim][i_array] = np.zeros(3)
                else:
                    self.electrode_pos[i_channel_stim][i_array] = np.zeros(2)

        # parameter bounds for optimizer (constrain if desired)
        self._bounds = self.get_bounds(
            constrain_electrode_locations=self.constrain_electrode_locations,
            overlap_factor=self.overlap_factor,
        )

        # determine initial values
        self.x0 = self.get_init_vals(
            bounds=self._bounds, optimize=self.optimize_init_vals
        )

        # compile node arrays
        for _electrode in self.electrode:
            _electrode.compile_node_arrays()

        # plot electrodes
        if self.plot:
            for i_channel_stim in range(self.n_channel_stim):
                for i_array, _electrode_array in enumerate(
                    self.electrode[i_channel_stim]._electrode_arrays
                ):
                    _electrode_array.plot(
                        show=False,
                        fn_plot=os.path.join(
                            self.plot_folder,
                            f"electrode_stim_{i_channel_stim}_array_{i_array}.png",
                        ),
                    )

        # setup optimization
        ################################################################################################################
        logger.info( "Setting up optimization algorithm ...")

        # equal ROI weighting if None is provided
        if type(self.weights) is float or type(self.weights) is int:
            self.weights = None

        if self.weights is None:
            self.weights = np.ones(len(self._roi)) / len(self._roi)
        elif type(self.weights) is list:
            self.weights = np.array(self.weights)

        if type(self.goal) is not list:
            self.goal = [self.goal]

        if "focality" in self.goal and len(self.goal) != len(self._roi):
            self.goal = ["focality"] * len(self._roi)

        self._goal_dir = []

        for i_roi in range(len(self._roi)):
            nodes = None if i_roi >= len(self._nodes) else self._nodes[i_roi]
            con = None if i_roi >= len(self._con) else self._con[i_roi]
            vol, triangles_normals = get_element_properties(self.roi[i_roi] , nodes, con, len(self._roi[i_roi].center))
            self._vol.append(vol)
            if "normal" in self.e_postproc or "tangential" in self.e_postproc:
                self._goal_dir.append(triangles_normals)
            else:
                self._goal_dir.append(None)

        if (
            "focality" in self.goal or "focality_inv" in self.goal
        ) and self.threshold is None:
            raise ValueError(
                "Please define TESoptimze.threshold for focality optimization!"
            )

        if (
            type(self.threshold) != list and type(self.threshold) != np.ndarray
        ) and self.threshold is not None:
            self.threshold = [self.threshold]

        if type(self.e_postproc) is not list:
            self.e_postproc = [self.e_postproc]

        if len(self.e_postproc) != len(self._roi):
            self.e_postproc = self.e_postproc * len(self._roi)

        if self.track_focality and len(self._roi) != 2:
            raise ValueError(
                "Focality can not be computed and tracked with only one ROI (non-ROI missing)."
            )

        if not (
            (isinstance(self.goal[0], types.FunctionType) and len(self.goal) == 1)
            or "focality" in self.goal
            or "focality_inv" in self.goal
        ):
            assert len(self.goal) == len(
                self._roi
            ), "Please provide a goal function for each ROI."
            assert len(self.weights) == len(
                self._roi
            ), "Number of weights has to match the number ROIs"

        if "focality" in self.goal and len(self._roi) != 2:
            raise ValueError(
                "For focality optimization please provide ROI and non ROI region (in this order)."
            )

        # track goal fun value (in ROI 0) and focality measures for later analysis
        self.goal_fun_value = [[] for _ in range(self.n_channel_stim)]
        self.AUC = [[] for _ in range(self.n_channel_stim)]
        self.integral_focality = [[] for _ in range(self.n_channel_stim)]

        # direct and shgo optimizer do not take init vals
        if self.optimizer in ["direct", "shgo"]:
            self.optimize_init_vals = False

        # define gpc parameters for current estimator
        if self.electrode[0]._current_estimator is not None:
            if self.electrode[0]._current_estimator.method == "gpc":
                min_idx = 0
                max_idx = 0

                for i_channel_stim in range(self.n_channel_stim):
                    for _electrode_array in self.electrode[
                        i_channel_stim
                    ]._electrode_arrays:
                        if _electrode_array.optimize_alpha:
                            max_idx += 3
                        else:
                            max_idx += 2

                    self.electrode[i_channel_stim]._current_estimator.set_gpc_parameters(
                        lb=self._bounds.lb[min_idx:max_idx],
                        ub=self._bounds.ub[min_idx:max_idx],
                    )

                    min_idx = max_idx

        # set default options for optimizer
        self._optimizer_options_std["bounds"] = self._bounds
        self._optimizer_options_std["init_vals"] = self.x0
        self._optimizer_options_std["vol_tol"] = (
            1.0 / 3600000000.0 * 3 * np.sum(self.n_ele_free)
        )

        # insert user specific options
        if self.optimizer_options is not None:
            for key in self.optimizer_options:
                self._optimizer_options_std[key] = self.optimizer_options[key]

        # setup FEM
        ################################################################################################################
        # set dirichlet node to closest node of center of gravity of head model (indexing starting with 1)
        self.dirichlet_node = get_dirichlet_node_index_cog(mesh=self._mesh, roi=self._roi)

        # always compute e-field components (Ex, Ey, Ez), it will be postprocessed later according to self.e_postproc
        self.dataType = [1] * len(self._roi)

        # prepare FEM
        self._ofem = OnlineFEM(
            mesh=self._mesh,
            electrode=self.electrode,
            method="TES",
            roi=self._roi,
            anisotropy_type=self.anisotropy_type,
            solver_options=self.solver_options,
            fn_results=self.fn_results_hdf5,
            useElements=True,
            dataType=self.dataType,
            dirichlet_node=self.dirichlet_node,
        )

        # log summary
        ################################################################################################################
        logger.log(25, "=" * 100)
        logger.log(25, f"headmodel:                        {self._mesh.fn}")
        logger.log(25, f"n_roi:                            {self._n_roi}")
        logger.log(
            25, f"anisotropy type:                  {self._ofem.anisotropy_type}"
        )
        logger.log(25, f"n_channel_stim:                   {self.n_channel_stim}")
        logger.log(25, f"fn_eeg_cap:                       {self.fn_eeg_cap}")
        logger.log(
            25, f"fn_electrode_mask:                {self._fn_electrode_mask}"
        )
        logger.log(
            25, f"FEM solver options:               {self._ofem.solver_options}"
        )
        logger.log(
            25,
            f"dirichlet_correction:             {self.electrode[0].dirichlet_correction}",
        )
        logger.log(
            25,
            f"dirichlet_correction_detailed:    {self.electrode[0].dirichlet_correction_detailed}",
        )
        logger.log(
            25,
            f"current_outlier_correction:       {self.electrode[0].current_outlier_correction}",
        )
        logger.log(25, f"optimizer:                        {self.optimizer}")
        logger.log(25, f"goal:                             {self.goal}")
        logger.log(25, f"e_postproc:                       {self.e_postproc}")
        logger.log(25, f"threshold:                        {self.threshold}")
        logger.log(25, f"weights:                          {self.weights}")
        logger.log(25, f"output_folder:                    {self.output_folder}")
        logger.log(25, f"fn_results_hdf5:                  {self.fn_results_hdf5}")

        if self.optimizer_options is not None:
            for key in self.optimizer_options:
                if key != "bounds":
                    logger.log(
                        25, f"{key}:                {self.optimizer_options[key]}"
                    )

        for i_channel_stim in range(self.n_channel_stim):
            logger.log(
                25,
                f"Stimulation: {i_channel_stim} (n_ele_free: {self.n_ele_free[i_channel_stim]})",
            )
            logger.log(25, "---------------------------------------------")

            for i_array, _electrode_array in enumerate(
                self.electrode[i_channel_stim]._electrode_arrays
            ):
                logger.log(
                    25,
                    f"Electrode array [i_channel_stim][i_array]: [{i_channel_stim}][{i_array}]",
                )
                logger.log(25, f"\tn_ele: {(_electrode_array.n_ele)}")
                logger.log(25, f"\tcenter: {(_electrode_array.center)}")
                logger.log(25, f"\tradius: {(_electrode_array.radius)}")
                logger.log(25, f"\tlength_x: {(_electrode_array.length_x)}")
                logger.log(25, f"\tlength_y: {(_electrode_array.length_y)}")

        logger.log(25, "=" * 100)

        self._prepared = True

    def _set_logger(self, fname_prefix='simnibs_simulation', summary=True):
        """
        Set-up loggger to write to a file

        Parameters
        ----------
        fname_prefix: str, optional
            Prefix of log-file. Defaults to 'simnibs_simulation'.
        summary: bool, optional
            Create summary file 'fields_summary.txt'. Default: True.
        """
        if not os.path.isdir(self.output_folder):
            os.mkdir(self.output_folder)
        log_fn = os.path.join(
            self.output_folder,
            fname_prefix + '_{0}.log'.format(self.time_str))
        fh = logging.FileHandler(log_fn, mode='w')
        formatter = logging.Formatter(
            f'[ %(name)s {__version__} - %(asctime)s - %(process)d ]%(levelname)s: %(message)s')
        fh.setFormatter(formatter)
        fh.setLevel(logging.DEBUG)
        logger = logging.getLogger("simnibs")
        logger.addHandler(fh)
        self._log_handlers += [fh]

        if summary:
            fn_summary = os.path.join(self.output_folder, 'summary_AUTO.txt')
            fh_s = logging.FileHandler(fn_summary, mode='w')
            fh_s.setFormatter(logging.Formatter('%(message)s'))
            fh_s.setLevel(25)
            logger.addHandler(fh_s)
            self._log_handlers += [fh_s]
        simnibs_logger.register_excepthook(logger)
    
    def _finish_logger(self):
        logger = logging.getLogger("simnibs")
        [logger.removeHandler(lh) for lh in self._log_handlers]
        self._log_handlers = []
        simnibs_logger.unregister_excepthook()
        
    def add_electrode_layout(self, electrode_type, electrode=None):
        """
        Adds an electrode to the current TESoptimize

        Parameters
        ----------
        electrode: ElectrodeArrayPair or CircularArray class instance, optional, default=None
            Electrode structure.

        Returns
        -------
        electrode: ElectrodeLayout
            ElectrodeLayout structure added to TESoptimize
        """
        if electrode is None:
            if electrode_type == "CircularArray":
                electrode = CircularArray()
            elif electrode_type == "ElectrodeArrayPair":
                electrode = ElectrodeArrayPair()
        self.electrode.append(electrode)
        return electrode

    def add_roi(self, roi=None):
        """
        Adds an ROI to the current TESoptimize

        Parameters
        ----------
        roi: RegionOfInterest object, optional, default=None
            ROI structure.

        Returns
        -------
        roi: RegionOfInterestInitializer
            RegionOfInterestInitializer structure added to TESoptimize
        """
        if roi is None:
            roi = RegionOfInterest()
        self.roi.append(roi)
        return roi
    
    def to_dict(self) -> dict:
        """ Makes a dictionary storing all settings as key value pairs

        Returns
        --------------------
        dict
            Dictionary containing settings as key value pairs
        """
        # Generate dict from instance variables (excluding variables starting with _ or __)
        settings = {
            key:value for key, value in self.__dict__.items() 
            if not key.startswith('__')
            and not key.startswith('_')
            and not callable(value) 
            and not callable(getattr(value, "__get__", None))
            and value is not None
        }

        # Add class name as type (type is protected in python so it cannot be a instance variable)
        settings['type'] = 'TesFlexOptimization'

        roi_dicts = []
        for roi_class in self.roi:
            roi_dicts.append(roi_class.to_dict())
        settings["roi"] = roi_dicts

        electrode_dicts = []
        for electrode_class in self.electrode:
            electrode_dicts.append(electrode_class.to_dict())
        settings["electrode"] = electrode_dicts    

        return settings

    def from_dict(self, settings: dict) -> "TesFlexOptimization":
        """ Reads parameters from a dict

        Parameters
        ----------
        settings: dict
            Dictionary containing parameter as key value pairs
        """

        for key, value in self.__dict__.items():
            if key.startswith('__') or key.startswith('_') or callable(value) or callable(getattr(value, "__get__", None)):
                continue
            setattr(self, key, settings.get(key, value))

        self.roi = []
        if 'roi' in settings:
            roi_dicts = settings["roi"]
            if isinstance(roi_dicts, dict):
                roi_dicts = [roi_dicts]
            for roi_dict in roi_dicts:
                self.roi.append(RegionOfInterest(roi_dict))

        self.electrode = []
        if 'electrode' in settings:
            electrode_dicts = settings["electrode"]
            if isinstance(electrode_dicts, dict):
                electrode_dicts = [electrode_dicts]
            for elec_dict in electrode_dicts:
                if elec_dict['type'] == 'ElectrodeArrayPair':
                    self.electrode.append(ElectrodeArrayPair().from_dict(elec_dict))
                elif elec_dict['type'] == 'CircularArray':
                    self.electrode.append(CircularArray().from_dict(elec_dict))

        return self

    def run(self, cpus=None, allow_multiple_runs=False, save_mat=True, return_n_max=1):
        """
        Runs the tes optimization

        Parameters
        ----------
        cpus : int, optional, default: None
            Number of CPU cores to use (so far used only during ellipsoid-fitting; 
                                        still ignored during FEM)
        allow_multiple_runs: bool, optional, default: False
            Whether to allow multiple runs in one folder. (not implemented yet)
        save_mat: bool, optional, default: True
            Save the ".mat" file of this structure
        return_n_max: int, optional, default: 1
            Return n-th best solutions. (not implemented yet)

        Returns
        --------
        <files>: Results files (.hdf5) in self.output_folder.
        """
        self._set_logger()
        self._n_cpu = cpus
        
        # prepare optimization
        if not self._prepared:
            self._prepare()

        # save structure in .mat format
        if save_mat:
            mat = self.to_dict()
            scipy.io.savemat(
                os.path.join(
                    self.output_folder,
                    "simnibs_simulation_{0}.mat".format(self.time_str),
                ),
                mat,
            )

        # run optimization
        ################################################################################################################
        start = time.time()
        if self.optimizer == "direct":
            result = direct(
                self.goal_fun,
                bounds=self._optimizer_options_std["bounds"],
                vol_tol=self._optimizer_options_std["vol_tol"],
                len_tol=self._optimizer_options_std["len_tol"],
                f_min_rtol=self._optimizer_options_std["f_min_rtol"],
                maxiter=self._optimizer_options_std["maxiter"],
                locally_biased=self._optimizer_options_std["locally_biased"],
            )

        elif self.optimizer == "Nelder-Mead":
            result = minimize(
                self.goal_fun,
                self._optimizer_options_std["init_vals"],
                method="Nelder-Mead",
                bounds=self._optimizer_options_std["bounds"],
                options={"disp": self._optimizer_options_std["disp"]},
            )

        elif self.optimizer == "differential_evolution":
            result = differential_evolution(
                self.goal_fun,
                x0=self._optimizer_options_std["init_vals"],
                strategy="best1bin",
                recombination=self._optimizer_options_std["recombination"],
                mutation=tuple(self._optimizer_options_std["mutation"]),
                tol=self._optimizer_options_std["tol"],
                maxiter=self._optimizer_options_std["maxiter"],
                popsize=self._optimizer_options_std["popsize"],
                bounds=self._optimizer_options_std["bounds"],
                disp=self._optimizer_options_std["disp"],
                polish=False,
            )  # we will decide if to polish afterwards

        elif self.optimizer == "shgo":
            result = shgo(
                self.goal_fun,
                bounds=self._optimizer_options_std["bounds"],
                options={"disp": self._optimizer_options_std["disp"]},
            )

        elif self.optimizer == "basinhopping":
            result = basinhopping(
                self.goal_fun,
                x0=self._optimizer_options_std["init_vals"],
                disp=self._optimizer_options_std["disp"],
            )
        else:
            raise NotImplementedError(
                f"Specified optimization method: '{self.optimizer}' not implemented."
            )

        logger.info(f"Optimization finished! Best electrode position: {result.x}")
        fopt_before_polish = result.fun
        stop = time.time()
        t_optimize = stop - start

        # polish optimization
        ################################################################################################################
        if self.polish:
            logger.info( "Polishing optimization results!")
            result = minimize(
                self.goal_fun,
                x0=result.x,
                method="L-BFGS-B",
                bounds=self._optimizer_options_std["bounds"],
                jac="2-point",
                options={"finite_diff_rel_step": 0.01},
            )
            logger.info(f"Optimization finished! Best electrode position: {result.x}")

        # transform electrode pos from array to list of list
        self.electrode_pos_opt = self.get_electrode_pos_from_array(result.x)

        fopt = result.fun
        nfev = result.nfev

        # plot final solution and electrode position (with node-wise dirichlet correction)
        ################################################################################################################
        for _electrode in self.electrode:
            _electrode.dirichlet_correction = True
            _electrode.dirichlet_correction_detailed = True

        # compute best e-field again, plot field and electrode position
        e = self.update_field(electrode_pos=self.electrode_pos_opt, plot=True)

        # postprocess e-field
        e_pp = [[0 for _ in range(self._n_roi)] for _ in range(self.n_channel_stim)]
        e_plot = [[] for _ in range(self._n_roi)]
        e_plot_label = [[] for _ in range(self._n_roi)]

        if np.array(["TI" in _t for _t in self.e_postproc]).any():
            for i_roi in range(self._n_roi):
                e_pp[0][i_roi] = postprocess_e(
                    e=e[0][i_roi],
                    e2=e[1][i_roi],
                    dirvec=self._goal_dir[i_roi],
                    type=self.e_postproc[i_roi],
                )
                e_pp[1][i_roi] = e_pp[0][i_roi]

                e_plot[i_roi].append(e[0][i_roi])
                e_plot[i_roi].append(e[0][i_roi])
                e_plot[i_roi].append(e_pp[0][i_roi])
                e_plot_label[i_roi].append("e_stim_0")
                e_plot_label[i_roi].append("e_stim_1")
                e_plot_label[i_roi].append("e_pp")

                # plot field
                if self.plot:
                    fn_out = os.path.join(self.plot_folder, f"e_roi_{i_roi}")
                    plot_roi_field(
                        e=e_plot[i_roi],
                        roi=self.roi[i_roi],
                        e_label=e_plot_label[i_roi],
                        fn_out=fn_out,
                    )
        else:
            for i_roi in range(self._n_roi):
                for i_channel_stim in range(self.n_channel_stim):
                    e_pp[i_channel_stim][i_roi] = postprocess_e(
                        e=e[i_channel_stim][i_roi],
                        e2=None,
                        dirvec=self._goal_dir[i_roi],
                        type=self.e_postproc[i_roi],
                    )
                    e_plot[i_roi].append(e[i_channel_stim][i_roi])
                    e_plot[i_roi].append(e_pp[i_channel_stim][i_roi])
                    e_plot_label[i_roi].append(f"e_stim_{i_channel_stim}")
                    e_plot_label[i_roi].append(f"e_pp_stim_{i_channel_stim}")

                if self.plot:
                    fn_out = os.path.join(self.plot_folder, f"e_roi_{i_roi}")
                    plot_roi_field(
                        e=e_plot[i_roi],
                        roi=self.roi[i_roi],
                        e_label=e_plot_label[i_roi],
                        fn_out=fn_out,
                    )

        # run final simulation with real electrode including remeshing
        ################################################################################################################
        if self.run_final_electrode_simulation:
            for i_channel_stim in range(self.n_channel_stim):
                s = create_tdcs_session_from_array(
                    electrode_array=self.electrode[i_channel_stim],
                    fnamehead=self._mesh.fn,
                    pathfem=os.path.join(
                        self.output_folder, f"final_sim_{i_channel_stim}"
                    ),
                )
                self.fn_final_sim.append(s.run()[0])

        # print optimization summary
        save_optimization_results(
            fname=os.path.join(self.output_folder, "summary"),
            optimizer=self.optimizer,
            optimizer_options=self._optimizer_options_std,
            fopt=fopt,
            fopt_before_polish=fopt_before_polish,
            popt=self.electrode_pos_opt,
            nfev=nfev,
            e=e,
            e_pp=e_pp,
            time=t_optimize,
            msh=self._mesh,
            electrode=self.electrode,
            goal=self.goal,
            n_test=self.n_test,
            n_sim=self.n_sim,
            n_iter_dirichlet_correction=self.n_iter_dirichlet_correction,
            goal_fun_value=self.goal_fun_value,
            AUC=self.AUC,
            integral_focality=self.integral_focality,
        )
        self._finish_logger()

    def goal_fun(self, parameters):
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
        self.n_test += 1
        parameters_str = f"Parameters: {parameters}"
        logger.info( parameters_str)

        # transform electrode pos from array to list of list
        self.electrode_pos = self.get_electrode_pos_from_array(parameters)

        # update field, returns list of list e[n_channel_stim][n_roi] (None if position is not applicable)
        e = self.update_field(electrode_pos=self.electrode_pos, plot=False)

        # post-process raw electric field (components Ex, Ey, Ez)
        if e is None:
            e_pp = None
        else:
            e_pp = [[0 for _ in range(self._n_roi)] for _ in range(self.n_channel_stim)]
            if np.array(["TI" in _t for _t in self.e_postproc]).any():
                for i_roi in range(self._n_roi):
                    e_pp[0][i_roi] = postprocess_e(
                        e=e[0][i_roi],
                        e2=e[1][i_roi],
                        dirvec=self._goal_dir[i_roi],
                        type=self.e_postproc[i_roi],
                    )
                    e_pp[1][i_roi] = e_pp[0][i_roi]
            else:
                for i_channel_stim in range(self.n_channel_stim):
                    for i_roi in range(self._n_roi):
                        e_pp[i_channel_stim][i_roi] = postprocess_e(
                            e=e[i_channel_stim][i_roi],
                            e2=None,
                            dirvec=self._goal_dir[i_roi],
                            type=self.e_postproc[i_roi],
                        )

        # compute goal function value
        goal_fun_value = self.compute_goal(e_pp)

        logger.info(
            f"Goal ({self.goal}): {goal_fun_value:.3f} (n_sim: {self.n_sim}, n_test: {self.n_test})",
        )
        logger.info( "-" * len(parameters_str))

        return goal_fun_value

    def compute_goal(self, e):
        """
        Computes goal function value from electric field.

        Parameters
        ----------
        e: list of list of np.ndarrays of float [n_channel_stim][n_roi][n_roi_nodes]
            Post-processed electric fields from simulated simulation conditions and ROIs.

        Returns
        -------
        goal_fun_value : float
            Accumulated goal function value. The average is taken over all stimulation conditions and the weighted
            average is taken according to self.weights over the different goal functions of the ROIs.
        """
        # calculate goal function value for every ROI
        y = np.zeros(
            (self.n_channel_stim, self._n_roi)
        )  # shape: [n_channel_stim x n_roi]

        if e is None:
            logger.info( f"Goal ({self.goal}): 2.0")
            return 2.0

        # user provided goal function
        ################################################################################################################
        elif isinstance(self.goal[0], types.FunctionType):
            y_weighted_sum = self.goal[0](e)

        # implemented goal functions
        else:
            # focality based goal functions
            ############################################################################################################
            if "focality" in self.goal:
                # TI focality (total field was previously calculated by the 2 channels, no loop over channel_stim here)
                if np.array(["TI" in _t for _t in self.e_postproc]).any():
                    y[:, :] = -100 * (
                        np.sqrt(2)
                        - ROC(
                            e1=e[0][0],  # e-field in ROI
                            e2=e[0][1],  # e-field in non-ROI
                            threshold=self.threshold,
                            focal=True,
                        )
                    )

                # General focality (can be different for each channel for e.g. TTF)
                else:
                    for i_channel_stim in range(self.n_channel_stim):
                        y[i_channel_stim, :] = -100 * (
                            np.sqrt(2)
                            - ROC(
                                e1=e[i_channel_stim][0],  # e-field in ROI
                                e2=e[i_channel_stim][1],  # e-field in non-ROI
                                threshold=self.threshold,
                                focal=True,
                            )
                        )

            elif "focality_inv" in self.goal:
                # TI focality (total field was previously calculated by the 2 channels, no loop over channel_stim here)
                if np.array(["TI" in _t for _t in self.e_postproc]).any():
                    y[:, :] = -100 * ROC(
                        e1=e[0][0],  # e-field in ROI
                        e2=e[0][1],  # e-field in non-ROI
                        threshold=self.threshold,
                        focal=False,
                    )

                # General focality (can be different for each channel for e.g. TTF)
                else:
                    for i_channel_stim in range(self.n_channel_stim):
                        y[i_channel_stim, :] = -100 * ROC(
                            e1=e[i_channel_stim][0],  # e-field in ROI
                            e2=e[i_channel_stim][1],  # e-field in non-ROI
                            threshold=self.threshold,
                            focal=False,
                        )

            # Mean/Max/ etc. based goal functions
            ############################################################################################################
            else:
                for i_roi in range(self._n_roi):
                    for i_channel_stim in range(self.n_channel_stim):
                        if e[i_channel_stim][i_roi] is None:
                            logger.info(
                                f"Goal ({self.goal}): 2.0 (one e-field was None)"
                            )
                            return 2.0
                        else:
                            # mean electric field in the roi
                            if self.goal[i_roi] == "mean":
                                y[i_channel_stim, i_roi] = -np.mean(
                                    e[i_channel_stim][i_roi]
                                )
                            # negative mean electric field in the roi (for e.g. normal component)
                            elif self.goal[i_roi] == "neg_mean":
                                y[i_channel_stim, i_roi] = np.mean(
                                    e[i_channel_stim][i_roi]
                                )
                            # mean of absolute value of electric field in the roi (for e.g. normal component)
                            elif self.goal[i_roi] == "mean_abs":
                                y[i_channel_stim, i_roi] = -np.mean(
                                    np.abs(e[i_channel_stim][i_roi])
                                )
                            # max electric field in the roi (percentile)
                            elif self.goal[i_roi] == "max":
                                y[i_channel_stim, i_roi] = -np.percentile(
                                    e[i_channel_stim][i_roi], 99.9
                                )
                            # negative max electric field in the roi (percentile)
                            elif self.goal[i_roi] == "neg_max":
                                y[i_channel_stim, i_roi] = np.percentile(
                                    e[i_channel_stim][i_roi], 99.9
                                )
                            # max of absolute value of electric field in the roi (percentile)
                            elif self.goal[i_roi] == "max_abs":
                                y[i_channel_stim, i_roi] = -np.percentile(
                                    np.abs(e[i_channel_stim][i_roi]), 99.9
                                )

            # if desired, track focality measures and goal function values
            for i_channel_stim in range(self.n_channel_stim):
                if self.track_focality:
                    if np.array(["max_TI" in _t for _t in self.e_postproc]).any():
                        e1 = get_maxTI(E1_org=e[0][0], E2_org=e[1][0])
                        e2 = get_maxTI(E1_org=e[0][1], E2_org=e[1][1])
                    elif np.array(["dir_TI" in _t for _t in self.e_postproc]).any():
                        e1 = get_dirTI(E1=e[0][0], E2=e[1][0], dirvec_org=self._goal_dir)
                        e2 = get_dirTI(E1=e[0][1], E2=e[1][1], dirvec_org=self._goal_dir)
                    else:
                        e1 = e[i_channel_stim][0]
                        e2 = e[i_channel_stim][1]

                    # compute integral focality
                    self.integral_focality[i_channel_stim].append(
                        integral_focality(
                            e1=e1, e2=e2, v1=self._vol[0], v2=self._vol[1]
                        )
                    )

                    # compute auc
                    self.AUC[i_channel_stim].append(AUC(e1=e1, e2=e2))
                else:
                    self.AUC[i_channel_stim].append(0)
                    self.integral_focality[i_channel_stim].append(0)

                # goal fun value in roi 0
                self.goal_fun_value[i_channel_stim].append(
                    np.mean(y[i_channel_stim, 0])
                )

            # average over all stimulations (channel_stim)
            ############################################################################################################
            y = np.mean(y, axis=0)

            # weight and sum the goal function values of the ROIs
            ############################################################################################################
            y_weighted_sum = np.sum(y * self.weights)

        self.n_sim += 1

        return y_weighted_sum

    def get_bounds(self, constrain_electrode_locations, overlap_factor=1.0):
        """
        Get boundaries of freely movable electrode arrays for optimizer.

        Parameters
        ----------
        constrain_electrode_locations : bool
            Constrains the possible locations of freely movable electrode arrays. Recommended for TTF optimizations,
            where two pairs of large electrode arrays are optimized. If True, parameter bounds for the optimization
            will be specified restricting the array locations to be frontal, parietal or occipital.
        overlap_factor : float, optional, default: 1.
            Factor of overlap of allowed lambda regions to place electrodes. (1 corresponds to neatless regions,
            <1 regions have a gap in between, >1 regions are overlapping)

        Returns
        -------
        bounds : scipy.optimize.Bounds instance [n_ele_free]
            Boundaries of freely movable electrode arrays tuple of length [n_ele_free] with lower bounds (lb) and
            upper bounds (ub) of beta, lambda, alpha.
        """

        if constrain_electrode_locations:
            # read fiducials from subject data
            fn_fiducials = os.path.join(self._ff_subject.eeg_cap_folder, "Fiducials.csv")
            with open(fn_fiducials, newline="") as f:
                reader = csv.reader(f)
                for row in reader:
                    if "Fiducial" in row and "Nz" in row:
                        # Nz (Nasion), Iz (Inion), LPHA (left ear), RPA (right ear)
                        Nz = np.array([row[1], row[2], row[3]]).astype(float)

            # project nasion to skin surface and determine normal vector
            con_skin = (
                self._mesh.elm.node_number_list[self._mesh.elm.tag1 == 1005,][:, :3] - 1
            )
            tri_skin_center = np.mean(self._mesh.nodes.node_coord[con_skin,], axis=1)
            idx_min = np.argmin(np.linalg.norm(tri_skin_center - Nz, axis=1))
            p1_tri = self._mesh.nodes.node_coord[con_skin[:, 0], :]
            p2_tri = self._mesh.nodes.node_coord[con_skin[:, 1], :]
            p3_tri = self._mesh.nodes.node_coord[con_skin[:, 2], :]
            tri_normal = np.cross(p2_tri - p1_tri, p3_tri - p1_tri)
            tri_normal /= np.linalg.norm(tri_normal, axis=1)[:, np.newaxis]
            Nz_normal = tri_normal[idx_min]

            # project nasion from skin surface to ellipsoid
            Nz_eli = subject2ellipsoid(
                coords=Nz, normals=Nz_normal, ellipsoid=self._ellipsoid
            )
            Nz_cart = self._ellipsoid.ellipsoid2cartesian(coords=Nz_eli, norm=False)
            Nz_jacobi = self._ellipsoid.cartesian2jacobi(coords=Nz_cart, norm=False)

            # go a small step into positive z-direction
            Nz_cart_test = self._ellipsoid.jacobi2cartesian(
                coords=Nz_jacobi + np.array([0, 1e-2]), norm=False
            )

            if Nz_cart_test[0, 2] > Nz_cart[0, 2]:
                lambda_sign = 1
            else:
                lambda_sign = -1

            # if Nz_jacobi[0, 1] + np.pi / 2 > np.pi:
            #     lambda_sign = -1
            # else:
            #     lambda_sign = 1
            beta_min = np.array([-np.pi / 4, -np.pi / 4, -np.pi / 2, np.pi / 4])
            beta_max = np.array([np.pi / 4, np.pi / 4, -np.pi / 4, np.pi / 2])
            lambda_min = np.array(
                [
                    Nz_jacobi[0, 1],
                    Nz_jacobi[0, 1] + lambda_sign * (np.pi / 2 + np.pi / 16),
                    -np.pi,
                    -np.pi,
                ]
            )
            lambda_max = np.array(
                [
                    Nz_jacobi[0, 1] + lambda_sign * (np.pi / 2 + np.pi / 16),
                    Nz_jacobi[0, 1] + lambda_sign * 3 * np.pi / 2 - np.pi / 8,
                    np.pi,
                    np.pi,
                ]
            )
            alpha_min = -np.pi * np.ones(np.sum(self.n_ele_free))
            alpha_max = np.pi * np.ones(np.sum(self.n_ele_free))

            # rearrange reg_lam_center such that electrode array pairs for stim on opposite sites are coming next
            # to each other in order
            lb = np.ravel([beta_min, lambda_min, alpha_min], "F")
            ub = np.ravel([beta_max, lambda_max, alpha_max], "F")

            # sort bounds (min, max)
            lbub = np.sort(np.vstack((lb, ub)), axis=0)
            lb = lbub[0, :]
            ub = lbub[1, :]

            # write txt files with some points for visualization
            beta_region_0 = np.linspace(lb[0], ub[0], 100)
            beta_region_1 = np.linspace(lb[3], ub[3], 100)
            beta_region_2 = np.linspace(lb[6], ub[6], 100)
            beta_region_3 = np.linspace(lb[9], ub[9], 100)
            lam_region_0 = np.linspace(lb[1], ub[1], 100)
            lam_region_1 = np.linspace(lb[4], ub[4], 100)
            lam_region_2 = np.linspace(lb[7], ub[7], 100)
            lam_region_3 = np.linspace(lb[10], ub[10], 100)

            coords_region_0_jac = np.array(
                np.meshgrid(beta_region_0, lam_region_0)
            ).T.reshape(-1, 2)
            coords_region_1_jac = np.array(
                np.meshgrid(beta_region_1, lam_region_1)
            ).T.reshape(-1, 2)
            coords_region_2_jac = np.array(
                np.meshgrid(beta_region_2, lam_region_2)
            ).T.reshape(-1, 2)
            coords_region_3_jac = np.array(
                np.meshgrid(beta_region_3, lam_region_3)
            ).T.reshape(-1, 2)

            coords_region_0 = self._ellipsoid.jacobi2cartesian(
                coords=coords_region_0_jac, return_normal=False
            )
            coords_region_1 = self._ellipsoid.jacobi2cartesian(
                coords=coords_region_1_jac, return_normal=False
            )
            coords_region_2 = self._ellipsoid.jacobi2cartesian(
                coords=coords_region_2_jac, return_normal=False
            )
            coords_region_3 = self._ellipsoid.jacobi2cartesian(
                coords=coords_region_3_jac, return_normal=False
            )

            np.savetxt(
                os.path.join(self.plot_folder, "coords_region_0.txt"), coords_region_0
            )
            np.savetxt(
                os.path.join(self.plot_folder, "coords_region_1.txt"), coords_region_1
            )
            np.savetxt(
                os.path.join(self.plot_folder, "coords_region_2.txt"), coords_region_2
            )
            np.savetxt(
                os.path.join(self.plot_folder, "coords_region_3.txt"), coords_region_3
            )

        else:
            # beta, lambda, alpha
            # order: [stim_1_array_1, stim_1_array_2, stim_2_array_1, stim_2_array_2, ... ]
            lb = np.array([-np.pi / 2, -np.pi, -np.pi] * np.sum(self.n_ele_free))
            ub = np.array([np.pi / 2, np.pi, np.pi] * np.sum(self.n_ele_free))

        # check if optimize_alpha is set for movable electrode_arrays, if not, remove alpha from the bounds
        i_para = 2
        idx_alpha_remove = []
        for i_channel_stim in range(self.n_channel_stim):
            for _electrode_array in self.electrode[i_channel_stim]._electrode_arrays:
                if not _electrode_array.optimize_alpha:
                    idx_alpha_remove.append(i_para)
                i_para += 3
        lb = np.delete(lb, idx_alpha_remove)
        ub = np.delete(ub, idx_alpha_remove)

        # TODO: think this works only for one channel_stim right now (HDTES), test with 2 channel stim and adapt
        # add bounds of geometry parameters of electrode (if any)
        for i_channel_stim in range(self.n_channel_stim):
            if self.electrode[i_channel_stim]._any_free_geometry:
                lb = np.append(
                    lb,
                    self.electrode[i_channel_stim]._geo_para_bounds[
                        self.electrode[i_channel_stim]._free_geometry, 0
                    ],
                )
                ub = np.append(
                    ub,
                    self.electrode[i_channel_stim]._geo_para_bounds[
                        self.electrode[i_channel_stim]._free_geometry, 1
                    ],
                )

        bounds = Bounds(lb=lb, ub=ub)

        return bounds

    def get_electrode_pos_from_array(self, electrode_pos_array):
        """
        Transforms an electrode_pos in array (1D) format to a list [n_channel_stim] of list [n_electrode_arrays] of
        numpy arrays (2 or 3, i.e. with or without alpha optimization).
        Changes electrode geometry in place in self.electrode in case of geometrical optimization.

        Parameters:
        -----------
        electrode_pos_array : np.array of float
            Electrode position in array format

        Returns
        -------
        electrode_pos : list of list of np.ndarray of float [2 or 3] of length [n_channel_stim][n_ele_free]
            Spherical coordinates (beta, lambda) and orientation angle (alpha) for each electrode array.
                      electrode array 1                        electrode array 2
            [ np.array([beta_1, lambda_1, alpha_1]),   np.array([beta_2, lambda_2, alpha_2]) ]
        """

        # extract electrode positions from optimal parameters
        electrode_pos = [[] for _ in range(self.n_channel_stim)]

        i_para = 0
        for i_channel_stim in range(self.n_channel_stim):
            for i_ele_free in range(self.n_ele_free[i_channel_stim]):
                if (
                    self.electrode[i_channel_stim]
                    ._electrode_arrays[i_ele_free]
                    .optimize_alpha
                ):
                    i_para_increment = 3
                else:
                    i_para_increment = 2
                electrode_pos[i_channel_stim].append(
                    electrode_pos_array[i_para : (i_para + i_para_increment)]
                )
                i_para += i_para_increment

        # TODO: same here I think it only works for one channel_stim
        # extract geometrical electrode parameters from optimal parameters and update electrode
        for i_channel_stim in range(self.n_channel_stim):
            if self.electrode[i_channel_stim]._any_free_geometry:
                n_free_parameters = np.sum(self.electrode[i_channel_stim]._free_geometry)
                self.electrode[i_channel_stim].set_geometrical_parameters_optimization(
                    electrode_pos_array[i_para : (i_para + n_free_parameters)]
                )
                i_para += n_free_parameters

        return electrode_pos

    def get_init_vals(self, bounds, optimize=True):
        """
        Determine initial values for optimization, guaranteeing a valid electrode position.

        Parameters
        ----------
        bounds : Bounds instance
            Lower and upper bounds of optimization problem
        optimize : bool, optional, default: True
            If True, find initial values for optimization, guaranteeing a valid solution. If False, initial values
            are the center between bounds.

        Returns
        -------
        x0 : ndarray of float [n_para]
            Initial values
        """

        if optimize:
            logger.info(
                "Finding valid initial values for electrode position for optimization.",
            )
            n_max = 5000
            n_para = len(bounds.lb)

            # make a list of possible parameter combinations within bounds
            para_test_grid = np.random.rand(n_max, n_para)
            para_test_grid = para_test_grid * (bounds.ub - bounds.lb) + bounds.lb
            para_test_grid_orig = copy.deepcopy(para_test_grid)

            for i in range(n_max):
                # transform electrode pos from array to list of list
                electrode_pos = self.get_electrode_pos_from_array(para_test_grid[i, :])

                # test position
                node_idx_dict = self.get_nodes_electrode(electrode_pos=electrode_pos)

                if type(node_idx_dict[0]) is str:
                    valid = False
                    electrode_pos_valid = node_idx_dict[1]

                    # write valid electrode positions into test grid to ease future iterations
                    i_para = 0
                    for i_channel_stim in range(self.n_channel_stim):
                        for i_ele_free in range(self.n_ele_free[i_channel_stim]):
                            if (
                                self.electrode[i_channel_stim]
                                ._electrode_arrays[i_ele_free]
                                .optimize_alpha
                            ):
                                i_para_increment = 3
                            else:
                                i_para_increment = 2

                            if (
                                electrode_pos_valid[i_channel_stim][i_ele_free]
                                is not None
                            ):
                                para_test_grid[
                                    i:, i_para : (i_para + i_para_increment)
                                ] = electrode_pos_valid[i_channel_stim][i_ele_free]
                            else:
                                para_test_grid[
                                    :, i_para : (i_para + i_para_increment)
                                ] = para_test_grid_orig[
                                    :, i_para : (i_para + i_para_increment)
                                ]
                            i_para += i_para_increment
                else:
                    valid = True

                logger.info(f"Testing position #{i + 1}: {para_test_grid[i, :]} -> {valid}")

                if not valid:
                    logger.info( f"> electrode_pos_valid: {node_idx_dict[1]}")

                if valid:
                    return para_test_grid[i, :]

        else:
            if self.n_channel_stim > 1 and type(self.init_pos) is str:
                raise AssertionError(
                    "Please provide a list of initial positions for each stimulation channel"
                    "containing a list of initial positions for each freely movable array."
                )

            if self.n_channel_stim == 1 and type(self.init_pos) is str:
                self.init_pos = [[self.init_pos]]

            if self.n_channel_stim == 1 and type(self.init_pos) is np.ndarray:
                self.init_pos = [[self.init_pos]]

            self.init_pos_subject_coords = [[] for _ in range(self.n_channel_stim)]

            init_pos_list = [
                ["Fz", "Pz", "P3", "F4"],  # init defaults for first i_channel_stim
                ["C3", "C4", "F3", "P4"],
            ]  # init defaults for second i_channel_stim

            # set initial positions of electrodes if nothing is provided
            assert self.n_channel_stim <= len(
                init_pos_list
            ), "Please provide initial electrode positions."

            if self.init_pos is None:
                self.init_pos = [0 for _ in range(self.n_channel_stim)]
                for i_channel_stim in range(self.n_channel_stim):
                    if self.n_ele_free[i_channel_stim] > len(
                        init_pos_list[i_channel_stim]
                    ):
                        raise NotImplementedError(
                            "Please specify initial coordinates or EEG electrode positions for each"
                            "freely movable electrode array (init_pos)!"
                        )
                    self.init_pos[i_channel_stim] = [
                        init_pos_list[i_channel_stim][i]
                        for i in range(self.n_ele_free[i_channel_stim])
                    ]

            # get subject coordinates of initial positions
            x0 = np.array([])
            for i_channel_stim in range(self.n_channel_stim):
                # user provided EEG electrode position as str (e.g. "C3", ...)
                if type(self.init_pos[i_channel_stim][0]) is str:
                    for eeg_pos in self.init_pos[i_channel_stim]:
                        tmp = ELECTRODE()
                        tmp.centre = eeg_pos
                        tmp.substitute_positions_from_cap(
                            cap=self._ff_subject.get_eeg_cap(cap_name=self.fn_eeg_cap)
                        )
                        self.init_pos_subject_coords[i_channel_stim].append(tmp.centre)
                # user provided coordinates in subject space as np.array
                else:
                    self.init_pos_subject_coords[i_channel_stim] = self.init_pos[
                        i_channel_stim
                    ]

                # transform initial positions from subject to ellipsoid space
                for i_ele_free, coords in enumerate(
                    self.init_pos_subject_coords[i_channel_stim]
                ):
                    # get closest point idx on subject surface
                    point_idx = np.argmin(
                        np.linalg.norm(coords - self._skin_surface.nodes, axis=1)
                    )

                    # electrode positon in ellipsoid space (jacobi coordinates)
                    self.electrode_pos[i_channel_stim][i_ele_free][:2] = (
                        self._ellipsoid.cartesian2jacobi(
                            coords=self._ellipsoid.ellipsoid2cartesian(
                                coords=subject2ellipsoid(
                                    coords=self._skin_surface.nodes[point_idx, :],
                                    normals=self._skin_surface.nodes_normals[
                                        point_idx, :
                                    ],
                                    ellipsoid=self._ellipsoid,
                                )
                            )
                        )
                    )

                    # set initial orientation alpha to zero
                    if len(self.electrode_pos[i_channel_stim][i_ele_free]) > 2:
                        self.electrode_pos[i_channel_stim][i_ele_free][2] = 0.0

                    # append position to initial values
                    x0 = np.append(x0, self.electrode_pos[i_channel_stim][i_ele_free])

                    # add geometric parameters if applicable
                    x0 = np.append(
                        x0,
                        self.electrode[i_channel_stim]._geo_para_mean[
                            self.electrode[i_channel_stim]._free_geometry
                        ],
                    )

        return x0

    def get_nodes_electrode(self, electrode_pos):
        """
        Assigns the skin points of the electrodes in electrode array and writes the points in
        electrode_array.electrodes[i].nodes and electrode_array.electrodes[i].node_area.
        Estimate optimal electrode currents based on previous simulations.

        Parameters
        ----------
        electrode_pos : list of list of np.ndarray of float [2 or 3] of length [n_channel_stim][n_ele_free]
            Spherical coordinates (beta, lambda) and orientation angle (alpha) for each electrode array.
                      electrode array 1                        electrode array 2
            [ np.array([beta_1, lambda_1, alpha_1]),   np.array([beta_2, lambda_2, alpha_2]) ]

        Returns
        -------
        node_idx_dict : list of dict
            List [n_channel_stim] containing dicts with electrode channel IDs as keys and node indices.
        """

        electrode_coords_subject = [0 for _ in range(self.n_channel_stim)]
        electrode_pos_valid = [
            [None for _ in range(len(self.electrode[i_channel_stim]._electrode_arrays))]
            for i_channel_stim in range(self.n_channel_stim)
        ]
        node_idx_dict = [dict() for _ in range(self.n_channel_stim)]
        node_coords_list = [[] for _ in range(np.sum(self.n_ele_free))]
        i_array_global = 0

        n = []
        cx = []
        cy = []

        for i_channel_stim in range(self.n_channel_stim):
            # collect all parameters
            start = np.zeros((self.n_ele_free[i_channel_stim], 3))
            a = np.zeros((self.n_ele_free[i_channel_stim], 3))
            b = np.zeros((self.n_ele_free[i_channel_stim], 3))
            cx.append(np.zeros((self.n_ele_free[i_channel_stim], 3)))
            cy.append(np.zeros((self.n_ele_free[i_channel_stim], 3)))
            n_tmp = np.zeros((self.n_ele_free[i_channel_stim], 3))
            start_shifted_ = np.zeros((len(electrode_pos[i_channel_stim]), 3))
            distance = []
            alpha = []

            for i_array, _electrode_array in enumerate(
                self.electrode[i_channel_stim]._electrode_arrays
            ):
                start[i_array, :] = self._ellipsoid.jacobi2cartesian(
                    coords=electrode_pos[i_channel_stim][i_array][:2]
                )

                c0, n_tmp[i_array, :] = self._ellipsoid.jacobi2cartesian(
                    coords=electrode_pos[i_channel_stim][i_array][:2],
                    return_normal=True,
                )
                a[i_array, :] = (
                    self._ellipsoid.jacobi2cartesian(
                        coords=np.array(
                            [
                                electrode_pos[i_channel_stim][i_array][0] - 1e-2,
                                electrode_pos[i_channel_stim][i_array][1],
                            ]
                        )
                    )
                    - c0
                )
                b[i_array, :] = (
                    self._ellipsoid.jacobi2cartesian(
                        coords=np.array(
                            [
                                electrode_pos[i_channel_stim][i_array][0],
                                electrode_pos[i_channel_stim][i_array][1] - 1e-2,
                            ]
                        )
                    )
                    - c0
                )
                a[i_array, :] /= np.linalg.norm(a[i_array, :])
                b[i_array, :] /= np.linalg.norm(b[i_array, :])

                if len(electrode_pos[i_channel_stim][i_array]) > 2:
                    start_shifted_[i_array, :] = c0 + (
                        1e-3
                        * (
                            (a[i_array, :])
                            * np.cos(electrode_pos[i_channel_stim][i_array][2])
                            + (b[i_array, :])
                            * np.sin(electrode_pos[i_channel_stim][i_array][2])
                        )
                    )
                else:
                    start_shifted_[i_array, :] = c0 + 1e-3 * a[i_array, :]

                cy[i_channel_stim][i_array, :] = (
                    start_shifted_[i_array, :] - start[i_array, :]
                )
                cy[i_channel_stim][i_array, :] /= np.linalg.norm(
                    cy[i_channel_stim][i_array, :]
                )
                cx[i_channel_stim][i_array, :] = np.cross(
                    cy[i_channel_stim][i_array, :], -n_tmp[i_array, :]
                )
                cx[i_channel_stim][i_array, :] /= np.linalg.norm(
                    cx[i_channel_stim][i_array, :]
                )

                distance.append(_electrode_array.distance)
                if _electrode_array.optimize_alpha:
                    alpha.append(
                        electrode_pos[i_channel_stim][i_array][2]
                        + _electrode_array.angle
                    )
                else:
                    alpha.append(_electrode_array.angle)

            distance = np.array(distance).flatten()
            alpha = np.array(alpha).flatten()

            for i_a, _alpha in enumerate(alpha):
                if _alpha > np.pi:
                    alpha[i_a] = _alpha - 2 * np.pi
                elif _alpha < -np.pi:
                    alpha[i_a] = _alpha + 2 * np.pi

            start = np.vstack(
                [
                    np.tile(start[i_array, :], (_electrode_array.n_ele, 1))
                    for i_array, _electrode_array in enumerate(
                        self.electrode[i_channel_stim]._electrode_arrays
                    )
                ]
            )

            electrode_array_idx = np.hstack(
                [
                    i_array * np.ones(_electrode_array.n_ele)
                    for i_array, _electrode_array in enumerate(
                        self.electrode[i_channel_stim]._electrode_arrays
                    )
                ]
            )

            # determine electrode center on ellipsoid
            if not (distance == 0.0).all():
                if sys.platform == 'win32' and self._n_cpu != 1:
                    n_cpu = 1
                    logger.info("Restricting geodesic destination calculatoins on Windows to one CPU core.")
                else:
                    n_cpu = self._n_cpu
                electrode_coords_eli_cart = self._ellipsoid.get_geodesic_destination(
                    start=start, distance=distance, alpha=alpha, n_steps=400, n_cpu=n_cpu
                )
            else:
                electrode_coords_eli_cart = start

            n.append(self._ellipsoid.get_normal(coords=electrode_coords_eli_cart))

            # transform to ellipsoidal coordinates
            electrode_coords_eli_eli = self._ellipsoid.cartesian2ellipsoid(
                coords=electrode_coords_eli_cart
            )

            # project coordinates to subject
            tmp_arrays = []
            i_ele = 0
            for i_array, _electrode_array in enumerate(
                self.electrode[i_channel_stim]._electrode_arrays
            ):
                ele_idx, tmp = ellipsoid2subject(
                    coords=electrode_coords_eli_eli[electrode_array_idx == i_array, :],
                    ellipsoid=self._ellipsoid,
                    surface=self._skin_surface,
                )
                tmp_arrays.append(tmp)

                if len(ele_idx) != len(alpha[electrode_array_idx == i_array]):
                    return (
                        "Electrode position: invalid (not all electrodes in valid skin region)",
                        electrode_pos_valid,
                    )
                else:
                    electrode_pos_valid[i_channel_stim][i_array] = electrode_pos[
                        i_channel_stim
                    ][i_array]
                    # print("Electrode position: invalid (not all electrodes in valid skin region)")

                electrode_coords_subject[i_channel_stim] = np.vstack(tmp_arrays)

                # loop over electrodes and determine node indices
                for _electrode in _electrode_array.electrodes:
                    if _electrode.type == "spherical":
                        # mask with a sphere
                        mask = (
                            np.linalg.norm(
                                self._skin_surface.nodes
                                - electrode_coords_subject[i_channel_stim][i_ele, :],
                                axis=1,
                            )
                            < _electrode.radius
                        )

                        # save position of electrode in subject space to posmat field
                        _electrode.posmat[:3, 3] = electrode_coords_subject[
                            i_channel_stim
                        ][i_ele, :]

                    elif _electrode.type == "rectangular":
                        cx_local = np.cross(
                            n[i_channel_stim][i_ele, :], cy[i_channel_stim][i_array, :]
                        )

                        # rotate skin nodes to normalized electrode space
                        rotmat = np.array(
                            [
                                [
                                    cx_local[0],
                                    cy[i_channel_stim][i_array, 0],
                                    n[i_channel_stim][i_ele, 0],
                                ],
                                [
                                    cx_local[1],
                                    cy[i_channel_stim][i_array, 1],
                                    n[i_channel_stim][i_ele, 1],
                                ],
                                [
                                    cx_local[2],
                                    cy[i_channel_stim][i_array, 2],
                                    n[i_channel_stim][i_ele, 2],
                                ],
                            ]
                        )
                        center = np.array(
                            [
                                electrode_coords_subject[i_channel_stim][i_ele, 0],
                                electrode_coords_subject[i_channel_stim][i_ele, 1],
                                electrode_coords_subject[i_channel_stim][i_ele, 2],
                            ]
                        )

                        # save position of electrode in subject space to posmat field
                        _electrode.posmat = np.vstack(
                            (
                                np.hstack((rotmat, center[:, np.newaxis])),
                                np.array([0, 0, 0, 1]),
                            )
                        )

                        skin_nodes_rotated = (self._skin_surface.nodes - center) @ rotmat

                        # mask with a box
                        mask_x = np.logical_and(
                            skin_nodes_rotated[:, 0] > -_electrode.length_x / 2,
                            skin_nodes_rotated[:, 0] < +_electrode.length_x / 2,
                        )
                        mask_y = np.logical_and(
                            skin_nodes_rotated[:, 1] > -_electrode.length_y / 2,
                            skin_nodes_rotated[:, 1] < +_electrode.length_y / 2,
                        )
                        mask_z = np.logical_and(
                            skin_nodes_rotated[:, 2] > -30,
                            skin_nodes_rotated[:, 2] < +30,
                        )
                        mask = np.logical_and(np.logical_and(mask_x, mask_y), mask_z)
                    else:
                        raise AssertionError(
                            "Electrodes have to be either 'spherical' or 'rectangular'"
                        )

                    # node areas
                    _electrode.node_area = self._skin_surface.nodes_areas[mask]

                    # total effective area of all nodes
                    _electrode.area_skin = _electrode.node_area_total

                    # electrode position is invalid if it overlaps with invalid skin region and area is not "complete"
                    if _electrode.area_skin < 0.90 * _electrode.area:
                        electrode_pos_valid[i_channel_stim][i_array] = None
                        # print("Electrode position: invalid (partly overlaps with invalid skin region)")
                        return (
                            "Electrode position: invalid (partly overlaps with invalid skin region)",
                            electrode_pos_valid,
                        )

                    # save node indices (referring to global mesh)
                    _electrode.node_idx = self._node_idx_msh[mask]

                    # save node coords (refering to global mesh)
                    _electrode.node_coords = self._skin_surface.nodes[mask]

                    node_coords_list[i_array_global].append(_electrode.node_coords)

                    # group node indices of same channel IDs
                    if _electrode.channel_id in node_idx_dict[i_channel_stim].keys():
                        node_idx_dict[i_channel_stim][_electrode.channel_id] = (
                            np.append(
                                node_idx_dict[i_channel_stim][_electrode.channel_id],
                                _electrode.node_idx,
                            )
                        )
                    else:
                        node_idx_dict[i_channel_stim][
                            _electrode.channel_id
                        ] = _electrode.node_idx

                    i_ele += 1

                electrode_pos_valid[i_channel_stim][i_array] = electrode_pos[
                    i_channel_stim
                ][i_array]

                # gather all electrode node coords of freely movable arrays
                node_coords_list[i_array_global] = np.vstack(
                    node_coords_list[i_array_global]
                )
                i_array_global += 1

        # check if electrode distance is sufficient
        invalid = False
        i_array_global_lst = np.hstack(
            [
                np.arange(self.n_ele_free[i_channel_stim])
                for i_channel_stim in range(self.n_channel_stim)
            ]
        ).astype(int)
        i_channel_stim_global_lst = np.hstack(
            [
                i_channel_stim * np.ones(self.n_ele_free[i_channel_stim])
                for i_channel_stim in range(self.n_channel_stim)
            ]
        ).astype(int)
        if self.min_electrode_distance is not None and self.min_electrode_distance >= 0:
            i_array_test_start = 1
            # start with first array and test if all node coords are too close to other arrays
            for i_array_global in range(np.sum(self.n_ele_free)):
                for node_coord in node_coords_list[i_array_global]:
                    for i_array_test in range(
                        i_array_test_start, np.sum(self.n_ele_free)
                    ):
                        # calculate euclidean distance between node coords
                        min_dist = np.min(
                            np.linalg.norm(
                                node_coords_list[i_array_test] - node_coord, axis=1
                            )
                        )
                        # stop testing if an electrode is too close
                        if min_dist <= self.min_electrode_distance:
                            # remove tested array from valid list
                            electrode_pos_valid[
                                i_channel_stim_global_lst[i_array_test]
                            ][i_array_global_lst[i_array_test]] = None
                            # print("Electrode position: invalid (minimal distance between electrodes too small)")
                            invalid = True

                i_array_test_start += 1

        if invalid:
            return (
                "Electrode position: invalid (minimal distance between electrodes too small)",
                electrode_pos_valid,
            )

        # save electrode_pos in ElectrodeArray instances
        for i_channel_stim in range(self.n_channel_stim):
            for i_array, _electrode_array in enumerate(
                self.electrode[i_channel_stim]._electrode_arrays
            ):
                _electrode_array.electrode_pos = electrode_pos[i_channel_stim][i_array]

        # estimate optimal electrode currents from previous simulations
        for i_channel_stim in range(self.n_channel_stim):
            if self.electrode[i_channel_stim]._current_estimator is not None:

                # estimate optimal currents electrode wise
                currents_estimate = self.electrode[i_channel_stim].estimate_currents(
                    electrode_pos[i_channel_stim]
                )

                # write currents in electrodes
                if currents_estimate is not None:
                    currents_estimate = currents_estimate.flatten()
                    for _electrode_array in self.electrode[
                        i_channel_stim
                    ]._electrode_arrays:
                        for _ele in _electrode_array.electrodes:
                            mask_estimator = (
                                self.electrode[i_channel_stim]._current_estimator.ele_id
                                == _ele.ele_id
                            ) * (
                                self.electrode[
                                    i_channel_stim
                                ]._current_estimator.channel_id
                                == _ele.channel_id
                            )
                            _ele.ele_current = currents_estimate[mask_estimator]
            else:
                # reset to original currents
                for _electrode_array in self.electrode[i_channel_stim]._electrode_arrays:
                    for _ele in _electrode_array.electrodes:
                        _ele.ele_current = _ele.ele_current_init

                self.electrode[i_channel_stim].compile_node_arrays()

        # compile node arrays
        for i_channel_stim in range(self.n_channel_stim):
            self.electrode[i_channel_stim].compile_node_arrays()

        return node_idx_dict

    def update_field(self, electrode_pos, plot=False):
        """
        Calculate the E field for given electrode positions.

        Parameters
        ----------
        electrode_pos : list of list of np.ndarray of float [3] of length [n_channel_stim][n_ele_free][3]
            Spherical coordinates (beta, lambda) and orientation angle (alpha) for each electrode array.
                      electrode array 1                        electrode array 2
            [ np.array([beta_1, lambda_1, alpha_1]),   np.array([beta_2, lambda_2, alpha_2]) ]
        plot : bool, optional, default: False
            Save data to plot e-field and electrode positions

        Returns
        -------
        e : list of list of np.ndarray [n_channel_stim][n_roi]
            Electric field for different stimulations in ROI(s).
        """

        e = [[0 for _ in range(self._n_roi)] for _ in range(self.n_channel_stim)]

        # assign surface nodes to electrode positions and estimate optimal currents
        # start = time.time()
        node_idx_dict = self.get_nodes_electrode(electrode_pos=electrode_pos)
        # stop = time.time()
        # print(f"Time: get_nodes_electrode: {stop-start}")

        # perform one electric field calculation for every stimulation condition (one at a time is on)
        for i_channel_stim in range(self.n_channel_stim):
            if type(node_idx_dict[0]) is str:
                logger.info( node_idx_dict[0])
                return None
            logger.info(f"Electrode position for stimulation {i_channel_stim}: valid")

            # set RHS
            b = self._ofem.set_rhs(electrode=self.electrode[i_channel_stim])

            # solve system
            if self.electrode[i_channel_stim].dirichlet_correction:
                if plot:
                    fn_electrode_txt = os.path.join(
                        self.plot_folder,
                        f"electrode_coords_nodes_subject_{i_channel_stim}.txt",
                    )
                else:
                    fn_electrode_txt = None

                v = self._ofem.solve_dirichlet_correction(
                    electrode=self.electrode[i_channel_stim],
                    fn_electrode_txt=fn_electrode_txt,
                )

                # store number of dirichlet iterations for convergence analysis
                self.n_iter_dirichlet_correction[i_channel_stim].append(
                    self._ofem.n_iter_dirichlet_correction
                )
            else:
                v = self._ofem.solve(b)

                if plot:
                    fn_electrode_txt = os.path.join(
                        self.plot_folder,
                        f"electrode_coords_nodes_subject_{i_channel_stim}.txt",
                    )
                    np.savetxt(
                        fn_electrode_txt,
                        np.hstack(
                            (
                                self.electrode[i_channel_stim]._node_coords,
                                self.electrode[i_channel_stim]._node_current[
                                    :, np.newaxis
                                ],
                            )
                        ),
                    )

            # Determine electric field in ROIs
            # start = time.time()
            for i_roi, r in enumerate(self._roi):
                if v is None:
                    e[i_channel_stim][i_roi] = None
                    logger.info("Warning! Simulation failed! Returning e-field: None!")
                else:
                    e[i_channel_stim][i_roi] = r.calc_fields(
                        v, dataType=self.dataType[i_roi]
                    )
            # stop = time.time()
            # print(f"Time: calc fields: {stop - start}")

        return e


def save_optimization_results(
    fname,
    optimizer,
    optimizer_options,
    fopt,
    fopt_before_polish,
    popt,
    nfev,
    e,
    e_pp,
    time,
    msh,
    electrode: list[ElectrodeLayout],
    goal,
    n_test=None,
    n_sim=None,
    n_iter_dirichlet_correction=None,
    goal_fun_value=None,
    AUC=None,
    integral_focality=None,
):
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
    e : list of np.ndarray [n_channel_stim][n_roi]
        List of list containing np.ndarrays of the (raw) electric fields in the ROIs (Ex, Ey, Ez)
    e_pp : list of np.ndarray [n_channel_stim][n_roi]
        List of list containing np.ndarrays of the postprocessed electric fields in the ROIs (norm or normal etc...)
    time : float
        Runtime of optimization in s
    electrode : list of ElectrodeArray objects [n_channel_stim]
        List of ElectrodeArray objects for every stimulation
    goal : str
        Goal function definition
    n_test : int, optional, default: None
        Number of runs to place the electrodes
    n_sim : int, optional, default: None
        Number of actual FEM simulations (for valid electrode placements only)
    n_iter_dirichlet_correction : list of int [n_channel_stim][n_iter], optional, default: None
        Number of iterations required to determine optimal currents in case of dirichlet correction
        for each call of "solve_dirichlet_correction"
    goal_fun_value : list of list of float [n_channel_stim][n_opt_runs]
        Goal function values of all stimulation conditions during optimization of ROI 0
    AUC : list of list of float [n_channel_stim][n_opt_runs]
        Area under curve focality measure for all stimulation conditions.
    integral_focality : list of list of float [n_channel_stim][n_opt_runs]
        Integral focality measure for all stimulation conditions.
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

    with open(fname_txt, "w") as f:
        f.write("Optimization summary:\n")
        f.write(
            "===================================================================\n"
        )
        f.write(f"Date: {d.year}-{d.month}-{d.day}, {d.hour}:{d.minute}:{d.second}\n")
        f.write(f"Simulation time: {str(datetime.timedelta(seconds=time))[:-7]}\n")
        f.write(f"headmodel: {msh.fn}\n")
        f.write(f"Goal: {goal}\n")
        f.write(f"Number of Channels: {len(e)}\n")
        f.write(f"Number of ROIs: {len(e[0])}\n")
        f.write("\n")
        f.write("Electrode coordinates:\n")
        f.write(
            "===================================================================\n"
        )
        f.write("Ellipsoid space (Jacobian coordinates):\n")
        f.write("---------------------------------------\n")

        for i_stim in range(len(popt)):
            f.write(f"Stimulation {i_stim}:\n")
            for i, p in enumerate(popt[i_stim]):
                f.write(f"Array {i}:\n")
                f.write(f"\tbeta:   {sep(p[0])}{p[0]:.3f}\n")
                f.write(f"\tlambda: {sep(p[1])}{p[1]:.3f}\n")
                if len(p) > 2:
                    f.write(f"\talpha:  {sep(p[2])}{p[2]:.3f}\n")
                else:
                    f.write(f"\talpha:  {sep(0.000)}{0.000}\n")
                
                electrode_array = electrode[i_stim]
                if isinstance(electrode_array, CircularArray):
                    f.write(
                        f"\tradius_inner:  {sep(electrode_array.radius_inner)}{electrode_array.radius_inner}\n"
                    )
                    f.write(
                        f"\tradius_outer:  {sep(electrode_array.radius_outer)}{electrode_array.radius_outer}\n"
                    )
                    f.write(
                        f"\tdistance:  {sep(electrode_array.distance)}{electrode_array.distance}\n"
                    )
                    f.write(
                        f"\tn_outer:  {sep(electrode_array.n_outer)}{electrode_array.n_outer}\n"
                    )

        f.write("\n")
        f.write("Subject space (Cartesian coordinates):\n")
        f.write("--------------------------------------\n")
        for i_stim in range(len(popt)):
            f.write(f"Stimulation {i_stim}:\n")
            for i_array, _electrode_array in enumerate(
                electrode[i_stim]._electrode_arrays
            ):
                f.write(f"Array {i_array}:\n")
                for i_electrode, _electrode in enumerate(_electrode_array.electrodes):
                    f.write(f"\tElectrode {i_electrode} ({_electrode.type}):\n")
                    for i_row in range(4):
                        f.write(
                            "\t\t"
                            + sep(_electrode.posmat[i_row, 0])
                            + f"{_electrode.posmat[i_row, 0]:.3f}, "
                            + sep(_electrode.posmat[i_row, 1])
                            + f"{_electrode.posmat[i_row, 1]:.3f}, "
                            + sep(_electrode.posmat[i_row, 2])
                            + f"{_electrode.posmat[i_row, 2]:.3f}, "
                            + sep(_electrode.posmat[i_row, 3])
                            + f"{_electrode.posmat[i_row, 3]:.3f}\n"
                        )
        f.write("\n")
        f.write("Optimization method:\n")
        f.write("===================================================================\n")
        f.write(f"Optimizer: {optimizer}\n")
        f.write("Settings:\n")
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
        f.write(f"fopt_before_polish: {fopt_before_polish}\n")

        if n_sim is not None:
            f.write(f"n_sim: {n_sim}\n")

        if n_test is not None:
            f.write(f"n_test: {n_test}\n")

    # save optimization settings and results in <fname>.hdf5 file
    ####################################################################################################################
    with h5py.File(fname_hdf5, "w") as f:
        # general info
        f.create_dataset(data=msh.fn, name="fnamehead")
        f.create_dataset(
            data=f"{d.year}-{d.month}-{d.day}, {d.hour}:{d.minute}:{d.second}",
            name="date",
        )
        f.create_dataset(data=len(e), name="n_channel")
        f.create_dataset(data=len(e[0]), name="n_roi")

        # optimizer
        f.create_dataset(data=optimizer, name="optimizer/optimizer")
        f.create_dataset(data=fopt, name="optimizer/fopt")
        f.create_dataset(data=fopt_before_polish, name="optimizer/fopt_before_polish")
        f.create_dataset(data=nfev, name="optimizer/nfev")
        f.create_dataset(data=time, name="optimizer/time")
        f.create_dataset(data=goal, name="optimizer/goal")

        if goal_fun_value is not None:
            f.create_dataset(
                data=np.array(goal_fun_value), name="optimizer/goal_fun_value"
            )

        if AUC is not None:
            f.create_dataset(data=np.array(AUC), name="optimizer/AUC")

        if integral_focality is not None:
            f.create_dataset(
                data=np.array(integral_focality), name="optimizer/integral_focality"
            )

        if n_iter_dirichlet_correction is not None:
            for i_channel_stim, n_iter in enumerate(n_iter_dirichlet_correction):
                f.create_dataset(
                    data=n_iter,
                    name=f"optimizer/n_iter_dirichlet_correction/channel_{i_channel_stim}",
                )

        if n_sim is not None:
            f.create_dataset(data=n_sim, name="optimizer/n_sim")

        if n_test is not None:
            f.create_dataset(data=n_test, name="optimizer/n_test")

        for key in optimizer_options:
            if type(optimizer_options[key]) is Bounds:
                f.create_dataset(
                    data=optimizer_options[key].lb,
                    name="optimizer/optimizer_options/lb",
                )
                f.create_dataset(
                    data=optimizer_options[key].ub,
                    name="optimizer/optimizer_options/ub",
                )
            else:
                f.create_dataset(
                    data=optimizer_options[key],
                    name=f"optimizer/optimizer_options/{key}",
                )

        # electrodes
        for i_stim, elec in enumerate(electrode):
            f.create_dataset(
                data=elec.center, name=f"electrode/channel_{i_stim}/center"
            )
            f.create_dataset(
                data=elec.length_x,
                name=f"electrode/channel_{i_stim}/length_x",
            )
            f.create_dataset(
                data=elec.length_y,
                name=f"electrode/channel_{i_stim}/length_y",
            )
            f.create_dataset(
                data=elec.current,
                name=f"electrode/channel_{i_stim}/current",
            )

            if isinstance(elec, CircularArray):
                f.create_dataset(
                    data=elec.radius_inner,
                    name=f"electrode/channel_{i_stim}/radius_inner",
                )
                f.create_dataset(
                    data=elec.radius_outer,
                    name=f"electrode/channel_{i_stim}/radius_outer",
                )
                f.create_dataset(
                    data=elec.distance,
                    name=f"electrode/channel_{i_stim}/distance",
                )
                f.create_dataset(
                    data=elec.n_outer,
                    name=f"electrode/channel_{i_stim}/n_outer",
                )
            
            if isinstance(elec, ElectrodeArrayPair):
                f.create_dataset(
                    data=elec.radius, name=f"electrode/channel_{i_stim}/radius"
                )
                

            for i_array, _electrode_array in enumerate(
                electrode[i_stim]._electrode_arrays
            ):
                f.create_dataset(
                    data=popt[i_stim][i_array][0],
                    name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/beta",
                )
                f.create_dataset(
                    data=popt[i_stim][i_array][1],
                    name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/lambda",
                )
                if len(popt[i_stim][i_array]) > 2:
                    f.create_dataset(
                        data=popt[i_stim][i_array][2],
                        name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/alpha",
                    )
                else:
                    f.create_dataset(
                        data=0.0,
                        name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/alpha",
                    )

                for i_electrode, _electrode in enumerate(_electrode_array.electrodes):
                    f.create_dataset(
                        data=_electrode.posmat,
                        name=f"electrode/channel_{i_stim}/posmat/electrode_array_{i_array}/electrode_{i_electrode}/posmat",
                    )

        # electric field in ROIs
        for i_stim in range(len(e)):
            for i_roi in range(len(e[i_stim])):
                f.create_dataset(
                    data=e[i_stim][i_roi], name=f"e/channel_{i_stim}/e_roi_{i_roi}"
                )
                f.create_dataset(
                    data=e_pp[i_stim][i_roi],
                    name=f"e_pp/channel_{i_stim}/e_roi_{i_roi}",
                )


def valid_skin_region(skin_surface, mesh, fn_electrode_mask, additional_distance=0.0):
    """
    Determine the nodes of the scalp surface where the electrode can be applied (not ears and face etc.)

    Parameters
    ----------
    skin_surface : Surface object
        Surface of the mesh (mesh_tools/surface.py)
    mesh : Msh object
        Mesh object created by SimNIBS (mesh_tools/mesh_io.py)
    additional_distance : float, optional, default: 0
        Additional distance in anterior part to put between original MNI template registration
    """
    nodes_all = copy.deepcopy(skin_surface.nodes)
    tr_nodes_all = copy.deepcopy(skin_surface.tr_nodes)
    # load mask of valid electrode positions (in MNI space)
    mask_img = nibabel.load(fn_electrode_mask)
    mask_img_data = mask_img.get_fdata()

    # add a certain distance to mask out closer to the eyes
    for i_border in range(mask_img_data.shape[0]):
        if mask_img_data[:, :, i_border].all():
            break

    mask_img_data[:, :, (i_border - additional_distance) : i_border] = 1

    # transform skin surface points to MNI space
    skin_nodes_mni_ras = subject2mni_coords(
        coordinates=skin_surface.nodes,
        m2m_folder=os.path.split(mesh.fn)[0],
        transformation_type="nonl",
    )

    # transform coordinates to voxel space
    skin_nodes_mni_voxel = (
        np.floor(
            np.linalg.inv(mask_img.affine)
            @ np.hstack(
                (
                    skin_nodes_mni_ras,
                    np.ones(skin_nodes_mni_ras.shape[0])[:, np.newaxis],
                )
            ).transpose()
        )[:3, :]
        .transpose()
        .astype(int)
    )
    skin_nodes_mni_voxel[skin_nodes_mni_voxel[:, 0] >= mask_img.shape[0], 0] = (
        mask_img.shape[0] - 1
    )
    skin_nodes_mni_voxel[skin_nodes_mni_voxel[:, 1] >= mask_img.shape[1], 1] = (
        mask_img.shape[1] - 1
    )
    skin_nodes_mni_voxel[skin_nodes_mni_voxel[:, 2] >= mask_img.shape[2], 2] = (
        mask_img.shape[2] - 1
    )

    # get boolean mask of valid skin points
    skin_surface.mask_valid_nodes = mask_img_data[
        skin_nodes_mni_voxel[:, 0],
        skin_nodes_mni_voxel[:, 1],
        skin_nodes_mni_voxel[:, 2],
    ].astype(bool)

    # remove points outside of MNI space (lower neck)
    skin_surface.mask_valid_nodes[(skin_nodes_mni_voxel < 0).any(axis=1)] = False

    skin_surface.mask_valid_tr = np.zeros(skin_surface.tr_centers.shape).astype(bool)

    unique_points = np.unique(
        skin_surface.tr_nodes[
            skin_surface.mask_valid_nodes[skin_surface.tr_nodes].all(axis=1), :
        ]
    )
    for point in unique_points:
        idx_where = np.where(skin_surface.tr_nodes == point)
        skin_surface.mask_valid_tr[idx_where[0], idx_where[1]] = True
    skin_surface.mask_valid_tr = skin_surface.mask_valid_tr.all(axis=1)

    # determine connectivity list of valid skin region (creates new node and connectivity list)
    skin_surface.nodes, skin_surface.tr_nodes = create_new_connectivity_list_point_mask(
        points=skin_surface.nodes,
        con=skin_surface.tr_nodes,
        point_mask=skin_surface.mask_valid_nodes,
    )

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
            nodes_idx_of_domain = np.unique(
                np.append(
                    nodes_idx_of_domain, skin_surface.tr_nodes[tri_idx_of_domain, :]
                )
            ).astype(int)
            tri_idx_of_domain = np.isin(skin_surface.tr_nodes, nodes_idx_of_domain).any(
                axis=1
            )
            n_current = np.sum(tri_idx_of_domain)
            # print(f"domain: {domain}, n_current: {n_current}")

        tri_domain[tri_idx_of_domain] = domain
        point_domain[nodes_idx_of_domain] = domain
        domain += 1

    domain_idx_main = np.argmax([np.sum(point_domain == d) for d in range(domain)])

    skin_surface.nodes, skin_surface.tr_nodes = create_new_connectivity_list_point_mask(
        points=skin_surface.nodes,
        con=skin_surface.tr_nodes,
        point_mask=point_domain == domain_idx_main,
    )

    # update masks
    skin_surface.mask_valid_nodes = np.zeros(nodes_all.shape[0]).astype(bool)
    for i_p, p in enumerate(nodes_all):
        if p in skin_surface.nodes:
            skin_surface.mask_valid_nodes[i_p] = True

    skin_surface.mask_valid_tr = np.zeros(tr_nodes_all.shape[0]).astype(bool)
    for i_t, t in enumerate(tr_nodes_all):
        if t in skin_surface.tr_nodes:
            skin_surface.mask_valid_tr[i_t] = True

    skin_surface.nodes_areas = skin_surface.nodes_areas[skin_surface.mask_valid_nodes]
    skin_surface.nodes_normals = skin_surface.nodes_normals[
        skin_surface.mask_valid_nodes, :
    ]
    skin_surface.surf2msh_nodes = skin_surface.surf2msh_nodes[
        skin_surface.mask_valid_nodes
    ]
    skin_surface.surf2msh_triangles = skin_surface.surf2msh_triangles[
        skin_surface.mask_valid_tr
    ]
    skin_surface.tr_areas = skin_surface.tr_areas[skin_surface.mask_valid_tr]
    skin_surface.tr_centers = skin_surface.tr_centers[skin_surface.mask_valid_tr, :]
    skin_surface.tr_normals = skin_surface.tr_normals[skin_surface.mask_valid_tr, :]

    return skin_surface


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
    c2 = ellipsoid.ellipsoid2cartesian(
        coords=np.array([electrode_pos[0] - 1e-3, electrode_pos[1]])
    )
    c3 = ellipsoid.ellipsoid2cartesian(
        coords=np.array([electrode_pos[0], electrode_pos[1] - 1e-3])
    )
    c4 = c1 + (
        1e-3
        * ((c3 - c1) * np.sin(electrode_pos[2]) + (c2 - c1) * np.cos(electrode_pos[2]))
    )

    # create vector in direction of constant lambda
    l1 = ellipsoid.cartesian2jacobi(coords=c1)
    l2 = ellipsoid.jacobi2cartesian(coords=np.array([l1[0, 0] + 1e-3, l1[0, 1]]))

    l2c1 = ((l2 - c1) / np.linalg.norm(l2 - c1)).flatten()
    c4c1 = ((c4 - c1) / np.linalg.norm(c4 - c1)).flatten()
    angle_jac = np.arccos(np.dot((c4c1), (l2c1)))

    return angle_jac

def plot_roi_field(e, roi, fn_out, e_label=None):
    """
    Exports data to .txt files for plotting.

    Parameters
    ----------
    e : np.ndarray of float [n_roi_ele, ] or list of np.ndarray of float [n_e-fields]
        Electric fields to visualize. Multiple fields can be passed in a list. Each field has to match the ROI.
    roi : RegionOfInterest class instance
        RegionOfInterest Object the data is associated with.
    fn_out : str
        Prefix of output file name, will append *_geo.hdf5, *_data.hdf5, and *_data.xdmf
    e_label : str or list of str [n_e-fields]
        Data names
    """
    if type(e) is not list:
        e = [e]

    if e_label is None:
        e_label = [f"data_{str(i)}" for i in range(len(e))]

    if type(e_label) is not list:
        e_label = [e_label]

    np.savetxt(fn_out + "_nodes.txt", roi.nodes)
    np.savetxt(fn_out + "_data.txt", np.hstack(e))
    np.savetxt(fn_out + "_data_label.txt", e_label)

    # if we have a connectivity (surface):
    if roi.con is not None:
        np.savetxt(fn_out + "_con.txt", roi.con)

def get_element_properties(roi, nodes, con, n_center):

    roi_mesh = roi.get_roi_mesh()
    if roi_mesh.elm.nr > 0:
        if nodes is None:
            nodes = roi_mesh.nodes.node_coord
        
        if con is None:
            con = roi_mesh.elm.node_number_list - 1
    # determine element properties
    triangles_normals = None
    if nodes is not None and con is not None:
        # surface ROI
        if np.all(con[:, 3] < 0):
            p1 = nodes[con[:, 0], :]
            p2 = nodes[con[:, 1], :]
            p3 = nodes[con[:, 2], :]
            vol = 0.5 * np.linalg.norm(np.cross(p2 - p1, p3 - p1), axis=1)
            triangles_normals = np.cross(p2 - p1, p3 - p1)
            triangles_normals /= np.linalg.norm(triangles_normals, axis=1)[:, np.newaxis]
        # volume ROI
        else:
            p1 = nodes[con[:, 0], :]
            p2 = nodes[con[:, 1], :]
            p3 = nodes[con[:, 2], :]
            p4 = nodes[con[:, 3], :]
            vol = 1.0 / 6 * np.sum(np.multiply(np.cross(p2 - p1, p3 - p1), p4 - p1), axis=1)

    elif nodes is None or con is None:
        vol = np.ones(n_center)

    return vol, triangles_normals
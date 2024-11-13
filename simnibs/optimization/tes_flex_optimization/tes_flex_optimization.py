import copy
import csv
import os
import time
import h5py
import types
import logging
import sys
import nibabel

import numpy as np
import scipy.spatial
from scipy.optimize import (
    direct,
    Bounds,
    minimize,
    differential_evolution,
)

from simnibs import __version__
from simnibs.mesh_tools import mesh_io, surface
from simnibs.simulation.fem import get_dirichlet_node_index_cog
from simnibs.simulation.onlinefem import FemTargetPointCloud, OnlineFEM
from simnibs.utils import simnibs_logger
from simnibs.utils.simnibs_logger import logger
from simnibs.utils.region_of_interest import RegionOfInterest
from simnibs.utils.roi_result_visualization import RoiResultVisualization
from simnibs.utils.TI_utils import get_maxTI, get_dirTI
from simnibs.utils.file_finder import SubjectFiles, Templates
from simnibs.utils.transformations import (
    subject2mni_coords,
    create_new_connectivity_list_point_mask,
)
from simnibs.utils.mesh_element_properties import ElementTags, tissue_names

from .ellipsoid import Ellipsoid, subject2ellipsoid, ellipsoid2subject
from .measures import AUC, integral_focality, ROC
from .electrode_layout import (
    ElectrodeArray,
    ElectrodeArrayPair,
    create_tdcs_session_from_array,
    CircularArray,
)


class TesFlexOptimization:
    """
    Defines a TES optimization problem using node-wise current sources.
    
    Parameters (general)
    --------------------
    electrode : Electrode Object
        Electrode object containing ElectrodeArray instances
        (pre-implemented layouts: ElectrodeArrayPair, CircularArray)

    roi : list of RegionOfInterest class instances
        Region of interest(s) the field is evaluated in.
    
    goal : list of str [n_roi], or FunctionType, optional, default: ["mean"]
        Implemented or user provided goal functions:
        
        - "mean": maximize mean e-field in ROI
        - "max": maximize 99.9 percentile of electric field in ROI
        - "focality": Maximize focality  (goal: sensitivity = specificity = 1). NOTE: "focality" needs excatly two ROIs: The first will be treated as ROI, the second as non-ROI.
        - "focality_inv": Maximize inverse focality (goal: sensitivity(ROI) = 1, sensitivity(nonROI) = 1). NOTE: "focality" needs excatly two ROIs: The first will be treated as ROI, the second as non-ROI.
        - user provided function taking e-field as an input which is  a list of list of np.ndarrays of float [n_channel_stim][n_roi] containing np.array with e-field
    
    e_postproc : str, optional, default: "norm"
        Specifies how the raw e-field in the ROI (Ex, Ey, Ez) is post-processed:
        
        - "norm": electric field magnitude (default)
        - "normal": determine normal component (requires surface ROIS)
        - "tangential": determine tangential component (requires surface ROIS)
        - "max_TI": maximum envelope for temporal interference fields
        - "dir_TI": directional sensitive maximum envelope for temporal interference fields (requires surface ROIS)    
        
    min_electrode_distance : float, optional, default: 5.0
        Minimally ensured distance between electrodes of different arrays (in mm).    
    
    constrain_electrode_locations : bool, optional, default: False
        Constrains the possible locations of freely movable electrode arrays. Recommended for TTF optimizations,
        where two pairs of large electrode arrays are optimized. If True, parameter bounds for the optimization
        will be specified restricting the array locations to be frontal, parietal and occipital.
        
    overlap_factor : float, optional, default: 1
        Factor of overlap of allowed lambda regions to place electrodes. (1 corresponds to neatless regions,
        <1 regions have a gap in between, >1 regions are overlapping)
    
    weights : np.array of float [n_roi], optional, default: equal weights for each ROI
            Weights for optimizer for ROI specific goal function weighting 
            NOTE: Will be ignored for "focality" and "focality_inv" goals (see below),
            where ROI and non-ROI are combined into a single goal function value
        
    Parameters (optimizer + FEM)
    ------------------------------------
    optimizer : str, optional, default: "differential_evolution"
        Gobal optimization algorithm
    polish : bool, optional, default: False
            If True, then scipy.optimize.minimize with the L-BFGS-B method is used to polish the best
            population member at the end, which can improve the minimization.
    run_final_electrode_simulation : bool, optional, default: True
           Runs final simulation with optimized parameters using real electrode model including remeshing.
           Note: This is required to get final e-fields for visualization
    anisotropy_type : str, optional, default: 'scalar'
        Specify type of anisotropy for simulation ('scalar', 'vn' or 'mc')
    disable_SPR_for_volume_roi : bool, optional, default: True
            Whether to use SPR interpolation for volume rois
        
    Parameters (debugging)
    ----------------------
    initial_x0: numpy.array, optional, default: None
        starting values for optimization (will be automatically determined per default)
    detailed_results : bool, optional, default: False
        write detailed results into subfolder of output folder for visualization and control
    track_focality : bool, optional, default: False
        Tracks focality for each goal function value (requires ROI and non-ROI definition)

    """

    def __init__(self, settings_dict=None):
        """Initialized TESoptimize class instance"""
        # folders and I/O
        self.date = time.strftime("%Y-%m-%d %H:%M:%S")
        self.time_str = time.strftime("%Y%m%d-%H%M%S")
        self.output_folder = None
        self.run_final_electrode_simulation = True
        self.open_in_gmsh = True
        self.detailed_results = False
        self._detailed_results_folder = None
        self.fn_final_sim = []
        self._prepared = False

        # headmodel
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
        self._vol = []
        self.disable_SPR_for_volume_roi = True

        # electrode
        self.electrode: list[ElectrodeArray] = []
        self.electrode_pos = None
        self.electrode_pos_opt = None
        self.min_electrode_distance = 5.0
        self.n_channel_stim = None
        self.n_iter_dirichlet_correction = None
        self.n_ele_free = None

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
        self.initial_x0 = None
        self.x0 = None

        # track goal fun value (in ROI 0) and focality measures for later analysis
        self.optim_funvalue = None
        self.goal_fun_value = None
        self.AUC = None
        self.integral_focality = None

        # set default options for optimizer
        self.seed = None
        self.optimizer_options = None  # passed by user
        self._optimizer_options_std = {
            "len_tol": 1.0 / 3600000000.0,
            "f_min_rtol": 1e-12,
            "maxiter": 1000,
            "disp": True,
            "recombination": 0.7,  # differential evolution
            "mutation": [0.01, 0.5],  # differential evolution
            "popsize": 13,  # differential evolution
            "tol": 0.1,  # differential evolution
            "locally_biased": False,
            "seed": self.seed
        }

        # FEM
        self.dirichlet_node = None
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
        ####################################################################################################
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
        ####################################################################################################
        logger.info( "Setting up output folders, logging and IO ...")
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        
        if self.detailed_results:
            self._detailed_results_folder = os.path.join(self.output_folder, "detailed_results")
            if not os.path.exists(self._detailed_results_folder):
                os.makedirs(self._detailed_results_folder)
        
        # setup headmodel
        ####################################################################################################
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

        # setup ROI
        ####################################################################################################
        logger.info( "Setting up ROI ...")

        if type(self.roi) is not list:
            self.roi = [self.roi]

        # initialize ROIs if not done already
        self._roi = []
        for i in range(len(self.roi)):
            if self.subpath is not None and self.roi[i].subpath is None and self.roi[i].mesh is None:
                    self.roi[i].subpath = self.subpath
            elif self.roi[i].subpath is None and self.roi[i].mesh is None:
                self.roi[i].mesh = self._mesh
            self._roi.append(FemTargetPointCloud(self._mesh, self.roi[i].get_nodes(), nearest_neighbor=((self.roi[i].method == "volume" or self.roi[i].method == "volume_from_surface") and self.disable_SPR_for_volume_roi)))

        self._n_roi = len(self._roi)

        # setup electrode
        ####################################################################################################
        logger.info( "Setting up electrodes ...")

        if type(self.electrode) is not list:
            self.electrode = [self.electrode]

        # initialize electrodes if not done already
        for i in range(len(self.electrode)):
            self.electrode[i]._prepare()

        # number of independent stimulation channels
        self.n_channel_stim = len(self.electrode)

        # plot electrodes (has to be called here because initial values may be different)
        for i_channel_stim in range(self.n_channel_stim):
            el_layout = self.electrode[i_channel_stim]
            usesDirichlet = el_layout.dirichlet_correction or el_layout.dirichlet_correction_detailed

            for i_array, _electrode_array in enumerate(
                    el_layout._electrode_arrays
            ):
                _electrode_array.plot(
                    show=False,
                    fn_plot=os.path.join(
                        self.output_folder,
                        f"electrode_channel_{i_channel_stim}_array_{i_array}.png",
                    ),
                    usesDirichlet=usesDirichlet
                )

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
        if self.initial_x0 is None:
            self.x0 = self.get_init_vals(bounds=self._bounds)
        else:
            self.x0 = self.initial_x0

        # compile node arrays
        for _electrode in self.electrode:
            _electrode.compile_node_arrays()

        # setup optimization
        ####################################################################################################
        logger.info( "Setting up optimization algorithm ...")

        if type(self.goal) is not list:
            self.goal = [self.goal]

        self._goal_dir = []

        # equal ROI weighting if None is provided
        if type(self.weights) is float or type(self.weights) is int:
            self.weights = None

        if (self.weights is None) and ("focality" not in self.goal) and ("focality_inv" not in self.goal):
            self.weights = np.ones(len(self._roi)) / len(self._roi)
        elif type(self.weights) is list:
            self.weights = np.array(self.weights)

        for i_roi in range(len(self._roi)):
            vol, node_normals = get_element_properties(self.roi[i_roi])
            self._vol.append(vol)
            if "normal" in self.e_postproc or "tangential" in self.e_postproc or "dir_TI" in self.e_postproc:
                self._goal_dir.append(node_normals)
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

        # track goal fun value (in ROI 0) and focality measures for later analysi
        n_effective_channels = self.n_channel_stim
        if np.array(["TI" in _t for _t in self.e_postproc]).any():
            n_effective_channels = 1 # for TI, fields of the channels are combined to a common metric
            
        self.goal_fun_value = [[] for _ in range(n_effective_channels)]
        self.AUC = [[] for _ in range(n_effective_channels)]
        self.integral_focality = [[] for _ in range(n_effective_channels)]

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
        self._optimizer_options_std["seed"] = self.seed
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
        ####################################################################################################
        # set dirichlet node to closest node of center of gravity of head model (indexing starting with 1)
        self.dirichlet_node = get_dirichlet_node_index_cog(mesh=self._mesh, roi=self._roi)

        # prepare FEM
        self._ofem = OnlineFEM(
            mesh=self._mesh,
            electrode=self.electrode,
            method="TES",
            roi=self._roi,
            anisotropy_type=self.anisotropy_type,
            solver_options=self.solver_options,
            fn_logger=False, # set to True to get more FEM details in log file
            useElements=True,
            dataType=[1] * len(self._roi),
            dirichlet_node=self.dirichlet_node,
            cpus=self._n_cpu
        )
        self._prepared = True
        
    def _set_logger(self, fname_prefix='simnibs_optimization', summary=True):
        """
        Set-up logger to write to a file

        Parameters
        ----------
        fname_prefix: str, optional
            Prefix of log-file. Defaults to 'simnibs_optimization'.
        summary: bool, optional
            Create summary file 'fields_summary.txt'. Default: True.
        """
        if not os.path.isdir(self.output_folder):
            os.makedirs(self.output_folder)
        log_fn = os.path.join(
            self.output_folder, f"{fname_prefix}_{self.time_str}.log"
        )
        fh = logging.FileHandler(log_fn, mode='w')
        formatter = logging.Formatter(
            f'[ %(name)s {__version__} - %(asctime)s - %(process)d ]%(levelname)s: %(message)s')
        fh.setFormatter(formatter)
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)
        self._log_handlers += [fh]

        if summary:
            fn_summary = os.path.join(self.output_folder, 'summary.txt')
            fh_s = logging.FileHandler(fn_summary, mode='w')
            fh_s.setFormatter(logging.Formatter('%(message)s'))
            fh_s.setLevel(26) # 25 is used by normal FEM for summary; using 26 here to not induced those summaries
            logger.addHandler(fh_s)
            self._log_handlers += [fh_s]
        simnibs_logger.register_excepthook(logger)
    
    def _finish_logger(self):
        [logger.removeHandler(lh) for lh in self._log_handlers]
        self._log_handlers = []
        simnibs_logger.unregister_excepthook()
        
    def _log_summary_preopt(self):
        ''' summary of optimization setup'''
        logger.log(26, "=" * 100)
        logger.log(26, "Optimization summary:")
        logger.log(26, "=" * 100)
        logger.log(26, f"Date: {self.date}")
        logger.log(26, f"Headmodel:                        {self._mesh.fn}")
        logger.log(26, f"Electrode_mask:                   {self._fn_electrode_mask}")
        logger.log(26, f"Conductivity anisotropy type:     {self._ofem.anisotropy_type}")
        logger.log(26, f"Output_folder:                    {self.output_folder}")
        logger.log(26, " ")
        
        logger.log(26, "Optimization and FEM settings")
        logger.log(26, "-" * 100)
        logger.log(26, f"Optimizer:                        {self.optimizer}")
        logger.log(26, f"Goal:                             {self.goal}")
        logger.log(26, f"Postprocessing:                   {self.e_postproc}")
        logger.log(26, f"Threshold:                        {self.threshold}")
        logger.log(26, f"Number of ROIs:                   {self._n_roi}")
        logger.log(26, f"Number of Channels:               {self.n_channel_stim}")
        logger.log(26, f"ROI weights:                      {self.weights}")
        logger.log(26, f"Constrain electrode locations:    {self.constrain_electrode_locations}")
        logger.log(26, f"Polish (local optimization):      {self.polish}")        
        logger.log(26, "Optimizer settings:")
        if self._optimizer_options_std is not None:
            for key in self._optimizer_options_std:
                if type(self._optimizer_options_std[key]) is Bounds:
                    logger.log(26, f"\tlb: {self._optimizer_options_std[key].lb}")
                    logger.log(26, f"\tub: {self._optimizer_options_std[key].ub}")
                else:
                    logger.log(26, f"\t{key}: {self._optimizer_options_std[key]}")
        else:
            logger.log(26, "None")
        logger.log(26, " ")
        
        logger.log(26, f"FEM solver options:               {self._ofem.solver_options}")
        logger.log(26, f"Dirichlet correction:             {self.electrode[0].dirichlet_correction}")
        logger.log(26, f"Dirichlet correction (detailed):  {self.electrode[0].dirichlet_correction_detailed}")
        logger.log(26, f"Current outlier correction:       {self.electrode[0].current_outlier_correction}")
        logger.log(26, " ")
        
        logger.log(26, "=" * 100)
        logger.log(26, "Optimization results:")
        logger.log(26, "=" * 100)

                
    def _log_summary_postopt(self):
        ''' summary of optimization results'''
        
        def sep(x):
            if x >= 0:
                return " "
            return ""        
        
        logger.log(26, "Electrode array infos:")
        logger.log(26, "=" * 100)
        for i_stim in range(len(self.electrode_pos_opt)):
           electrode_array = self.electrode[i_stim]
           logger.log(26, f"Stimulation channel {i_stim}:")
           logger.log(26, f"\tcurrents: {electrode_array.current}")
           if isinstance(electrode_array, CircularArray):
                logger.log(26, f"\tradius_inner:  {sep(electrode_array.radius_inner)}{electrode_array.radius_inner}")
                logger.log(26, f"\tradius_outer:  {sep(electrode_array.radius_outer)}{electrode_array.radius_outer}")
                logger.log(26, f"\tdistance:  {sep(electrode_array.distance)}{electrode_array.distance}")
                logger.log(26, f"\tn_outer:  {sep(electrode_array.n_outer)}{electrode_array.n_outer}")
        logger.log(26, " ")
        
        logger.log(26, "Electrode coordinates (Cartesian space):")
        logger.log(26, "-" * 100)
        for i_stim in range(len(self.electrode_pos_opt)):
           logger.log(26, f"Stimulation channel {i_stim}:")
           for i_array, _electrode_array in enumerate(self.electrode[i_stim]._electrode_arrays):
               logger.log(26, f"Array {i_array}:")
               for i_electrode, _electrode in enumerate(_electrode_array.electrodes):
                   logger.log(26, f"\tElectrode {i_electrode} ({_electrode.type}):")
                   for i_row in range(4):
                       logger.log(26, 
                           "\t\t"
                           + sep(_electrode.posmat[i_row, 0])
                           + f"{_electrode.posmat[i_row, 0]:.3f}, "
                           + sep(_electrode.posmat[i_row, 1])
                           + f"{_electrode.posmat[i_row, 1]:.3f}, "
                           + sep(_electrode.posmat[i_row, 2])
                           + f"{_electrode.posmat[i_row, 2]:.3f}, "
                           + sep(_electrode.posmat[i_row, 3])
                           + f"{_electrode.posmat[i_row, 3]:.3f}"
                       )

    def _write_detailed_results_preopt(self):
        ''' write out some more results prior to optimization '''

        # save skin surface
        np.savetxt(
            os.path.join(self._detailed_results_folder, "skin_surface_nodes.txt"),
            self._skin_surface.nodes,
        )
        np.savetxt(
            os.path.join(self._detailed_results_folder, "skin_surface_con.txt"),
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
            os.path.join(self._detailed_results_folder, "fitted_ellipsoid.txt"), eli_coords_jac
        )
                
        # plot valid skin regions
        if self.constrain_electrode_locations:
            beta = []        
            lam = []
            min_idx = 0
            max_idx = 0
            for i_channel_stim in range(self.n_channel_stim):
                for _electrode_array in self.electrode[
                    i_channel_stim
                ]._electrode_arrays:
                    print(min_idx)
                    if _electrode_array.optimize_alpha:
                        max_idx += 3
                    else:
                        max_idx += 2    
                    beta.append([self._bounds.lb[min_idx], self._bounds.ub[min_idx]])
                    lam.append([self._bounds.lb[min_idx+1], self._bounds.ub[min_idx+1]])
                    min_idx = max_idx
            
            for i in range(len(beta)):
               beta_region = np.linspace(beta[i][0], beta[i][1], 100)
               lam_region = np.linspace(lam[i][0], lam[i][1], 100)

               coords_region_jac = np.array(
                    np.meshgrid(beta_region, lam_region)
                ).T.reshape(-1, 2)
 
               coords_region = self._ellipsoid.jacobi2cartesian(
                    coords=coords_region_jac, return_normal=False
                )

               np.savetxt(
                    os.path.join(self._detailed_results_folder, f"coords_region_{i}.txt"), coords_region
                )
  
    def _write_detailed_results_postopt(self, fopt=None):
        ''' write out some more results after to optimization
        
            The fields are re-evaluated using the online FEM with dirichlet corrections,
            and saved for visualization as text files.
            
            In addition, the final fields and the optimization settings and (intermediate)
            results are saved in a hdf5-file which contains the following:
                  optimizer : str
                      Name of optimization method.
                  optimizer_options : dict
                      Dictionary containing the optimization setting.
                  fopt : float
                      Objective function value in optimum.
                  popt : list of list of np.ndarray of float [n_channel_stim][n_free_electrodes]
                      List containing the optimal parameters for each freely movable electrode array
                      [np.array([beta_1, lambda_1, alpha_1]), np.array([beta_2, lambda_2, alpha_2]), ...]
                  e : list of np.ndarray [n_channel_stim][n_roi]
                      List of list containing np.ndarrays of the (raw) electric fields in the ROIs (Ex, Ey, Ez)
                  e_pp : list of np.ndarray [n_channel_stim][n_roi]
                      List of list containing np.ndarrays of the postprocessed electric fields in the ROIs (norm or normal etc...)
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
        '''
    
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

            for i in range(len(e)):
                if e[i].ndim == 1:
                    e[i] = e[i].reshape((-1,1))
            np.savetxt(fn_out + "_nodes.txt", roi.get_nodes())
            np.savetxt(fn_out + "_data.txt", np.hstack(e))
            np.savetxt(fn_out + "_data_label.txt", e_label, fmt='%s')
                
        # get final solution and electrode position with node-wise dirichlet correction for external plotting
        for _electrode in self.electrode:
             _electrode.dirichlet_correction = True
             _electrode.dirichlet_correction_detailed = True

        # compute best e-field again, plot field and electrode position
        e = self.update_field(electrode_pos=self.electrode_pos_opt, plot=True)

        # postprocess e-field
        e_plot = [[] for _ in range(self._n_roi)]
        e_plot_label = [[] for _ in range(self._n_roi)]

        if np.array(["TI" in _t for _t in self.e_postproc]).any():
            
            e_pp = [[0 for _ in range(self._n_roi)]]
            
            for i_roi in range(self._n_roi):
                 e_pp[0][i_roi] = postprocess_e(
                     e=e[0][i_roi],
                     e2=e[1][i_roi],
                     dirvec=self._goal_dir[i_roi],
                     type=self.e_postproc[i_roi],
                 )
                 e_plot[i_roi].append(e[0][i_roi])
                 e_plot[i_roi].append(e[1][i_roi])
                 e_plot[i_roi].append(e_pp[0][i_roi])
                 e_plot_label[i_roi].append("e_stim_0_Ex")
                 e_plot_label[i_roi].append("e_stim_0_Ey")
                 e_plot_label[i_roi].append("e_stim_0_Ez")
                 e_plot_label[i_roi].append("e_stim_1_Ex")
                 e_plot_label[i_roi].append("e_stim_1_Ey")
                 e_plot_label[i_roi].append("e_stim_1_Ez")
                 e_plot_label[i_roi].append("e_pp")
                 
                 fn_out = os.path.join(self._detailed_results_folder, f"e_roi_{i_roi}")
                 plot_roi_field(
                   e=e_plot[i_roi],
                   roi=self.roi[i_roi],
                   e_label=e_plot_label[i_roi],
                   fn_out=fn_out,
                 )
        else:
            
             e_pp = [[0 for _ in range(self._n_roi)] for _ in range(self.n_channel_stim)]
             
             for i_roi in range(self._n_roi):
                for i_channel_stim in range(self.n_channel_stim):
                     e_pp[i_channel_stim][i_roi] = postprocess_e(
                         e=e[i_channel_stim][i_roi],
                         e2=None,
                         dirvec=self._goal_dir[i_roi],
                         type=self.e_postproc[i_roi],
                     ).reshape((-1, 1))
                     e_plot[i_roi].append(e[i_channel_stim][i_roi])
                     e_plot[i_roi].append(e_pp[i_channel_stim][i_roi])
                     e_plot_label[i_roi].append(f"e_stim_{i_channel_stim}_Ex")
                     e_plot_label[i_roi].append(f"e_stim_{i_channel_stim}_Ey")
                     e_plot_label[i_roi].append(f"e_stim_{i_channel_stim}_Ez")
                     e_plot_label[i_roi].append(f"e_pp_stim_{i_channel_stim}")
                
                fn_out = os.path.join(self._detailed_results_folder, f"e_roi_{i_roi}")
                plot_roi_field(
                    e=e_plot[i_roi],
                    roi=self.roi[i_roi],
                    e_label=e_plot_label[i_roi],
                    fn_out=fn_out,
                )
                
        # save optimization settings and results in <fname>.hdf5 file
        fname_hdf5 = os.path.join(self._detailed_results_folder, "summary.hdf5")
        with h5py.File(fname_hdf5, "w") as f:
            # general info
            f.create_dataset(data=self._mesh.fn, name="fnamehead")
            f.create_dataset(
                data=f"{self.date}",
                name="date",
            )
            f.create_dataset(data=len(e), name="n_channel")
            f.create_dataset(data=len(e[0]), name="n_roi")

            # optimizer
            f.create_dataset(data=self.optimizer, name="optimizer/optimizer")
            f.create_dataset(data=fopt, name="optimizer/fopt")
            f.create_dataset(data=self.goal, name="optimizer/goal")

            if self.goal_fun_value is not None:
                f.create_dataset(
                    data=np.array(self.goal_fun_value), name="optimizer/goal_fun_value"
                )

            if self.AUC is not None:
                f.create_dataset(data=np.array(self.AUC), name="optimizer/AUC")

            if self.integral_focality is not None:
                f.create_dataset(
                    data=np.array(self.integral_focality), name="optimizer/integral_focality"
                )

            if self.n_iter_dirichlet_correction is not None:
                for i_channel_stim, n_iter in enumerate(self.n_iter_dirichlet_correction):
                    f.create_dataset(
                        data=n_iter,
                        name=f"optimizer/n_iter_dirichlet_correction/channel_{i_channel_stim}",
                    )

            if self.n_test is not None:
                f.create_dataset(data=self.n_test, name="optimizer/n_test")

            if self.n_sim is not None:
                f.create_dataset(data=self.n_sim, name="optimizer/n_sim")

            for key in self._optimizer_options_std:
                if type(self._optimizer_options_std[key]) is Bounds:
                    f.create_dataset(
                        data=self._optimizer_options_std[key].lb,
                        name="optimizer/optimizer_options/lb",
                    )
                    f.create_dataset(
                        data=self._optimizer_options_std[key].ub,
                        name="optimizer/optimizer_options/ub",
                    )
                else:
                    if self._optimizer_options_std[key] is None:
                        data = "None"
                    else:
                        data = self._optimizer_options_std[key]

                    f.create_dataset(
                        data=data,
                        name=f"optimizer/optimizer_options/{key}"
                    )

            # electrodes
            for i_stim, elec in enumerate(self.electrode):
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
                    self.electrode[i_stim]._electrode_arrays
                ):
                    f.create_dataset(
                        data=self.electrode_pos_opt[i_stim][i_array][0],
                        name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/beta",
                    )
                    f.create_dataset(
                        data=self.electrode_pos_opt[i_stim][i_array][1],
                        name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/lambda",
                    )
                    if len(self.electrode_pos_opt[i_stim][i_array]) > 2:
                        f.create_dataset(
                            data=self.electrode_pos_opt[i_stim][i_array][2],
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
            for i_stim in range(len(e_pp)):
                for i_roi in range(len(e_pp[i_stim])):     
                    f.create_dataset(
                        data=e_pp[i_stim][i_roi],
                        name=f"e_pp/channel_{i_stim}/e_roi_{i_roi}",
                    )

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

    def run(self, cpus=None, save_mat=True):
        """
        Runs the tes optimization

        Parameters
        ----------
        cpus : int, optional, default: None
            Number of CPU cores to use (so far used only during ellipsoid-fitting; 
                                        still ignored during FEM)
        save_mat: bool, optional, default: True
            Save the ".mat" file of this structure

        Returns
        --------
        <files>: Results files (.hdf5) in self.output_folder.
        """

        start = time.time()
        self._set_logger()
        self._n_cpu = cpus

        if not cpus is None:
            from numba import set_num_threads
            set_num_threads(int(cpus))
            from numba import get_num_threads
            logger.info(f"Numba reports {get_num_threads()} threads available")
        
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
            
        # log settings in summary text
        self._log_summary_preopt()
        
        # write settings and preparation details (skin surface, fitted ellipsoid, electrodes)
        if self.detailed_results:
            self._write_detailed_results_preopt()
           
        # run global optimization
        ######################################################################################################
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
                seed=self.seed
            )  # we will decide if to polish afterwards

        else:
            raise NotImplementedError(
                f"Specified optimization method: '{self.optimizer}' not implemented."
            )
        
        self.optim_funvalue = result.fun
        optim_x = result.x
        
        logger.info(f"Global optimization finished! Best electrode position: {optim_x}")
        logger.log(26, f"Number of function evaluations (global optimization):   {self.n_test}")
        logger.log(26, f"Number of FEM evaluations (global optimization):        {self.n_sim}")
        logger.log(26, f"Goal function value (global optimization):              {self.optim_funvalue}")
        
        # run local optimization to polish results
        #######################################################################################################
        if self.polish:
            logger.info("Polishing optimization results!")
            result = minimize(
                self.goal_fun,
                x0=optim_x,
                method="L-BFGS-B",
                bounds=self._optimizer_options_std["bounds"],
                jac="2-point",
                options={"finite_diff_rel_step": 0.01},
            )
            logger.info(f"Local optimization finished! Best electrode position: {result.x}")
    
            if self.optim_funvalue <= result.fun:
                logger.info("Local optimization did not improve the results, proceeding with global optimization results.")
            else:
                optim_x = result.x
                self.optim_funvalue = result.fun
                
        # transform optimal electrode pos from array to list of list
        self.electrode_pos_opt = self.get_electrode_pos_from_array(optim_x)
        
        # internally update electrodes to correspond to optimal electrode pos
        self.get_nodes_electrode(electrode_pos=self.electrode_pos_opt)
        
        logger.log(26, f"Total number of function evaluations:                   {self.n_test}")
        logger.log(26, f"Total number of FEM evaluations:                        {self.n_sim}")
        logger.log(26, f"Final goal function value:                              {self.optim_funvalue}")
        logger.log(26, f"Duration (setup and optimization):                      {time.time() - start}")
                
        # run final simulation with real electrode including remeshing
        #########################################################################################################
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
                
            # extract e-fields from FEM simulations, add extra data and show results
            base_file_name = os.path.splitext(os.path.basename(self._ff_subject.fnamehead))[0]
            base_file_name += '_tes_flex_opt'
            
            fn_vis, m_head, m_surf = write_visualization(self.output_folder,
                                                         base_file_name,
                                                         self.roi, 
                                                         self.fn_final_sim,
                                                         self.e_postproc,
                                                         self.goal)
            
            # extract key metrics from m_head, m_surf and add to summary log
            logger.log(26, make_summary_text(m_surf, m_head))
            
            if self.open_in_gmsh:
                for i in fn_vis:
                    mesh_io.open_in_gmsh(i, True)
                        
        # append optimization results to summary
        self._log_summary_postopt()
        
        # write results details (final fields via onlineFEM with dirichlet_corrections as txt and hdf5)
        if self.detailed_results:
            self._write_detailed_results_postopt(self.optim_funvalue)
        
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

        if e is None:
            # overlap --> skip further steps
            logger.info( f"Goal ({self.goal}): 2.0")
            logger.info( "-" * len(parameters_str))
            return 2.0
        
        self.n_sim += 1
        
        # post-process raw electric field (components Ex, Ey, Ez)
        if np.array(["TI" in _t for _t in self.e_postproc]).any():
            e_pp = [[0 for _ in range(self._n_roi)]]
            for i_roi in range(self._n_roi):
                e_pp[0][i_roi] = postprocess_e(
                    e=e[0][i_roi],
                    e2=e[1][i_roi],
                    dirvec=self._goal_dir[i_roi],
                    type=self.e_postproc[i_roi],
                )
        else:
            e_pp = [[0 for _ in range(self._n_roi)] for _ in range(self.n_channel_stim)]
            for i_channel_stim in range(self.n_channel_stim):
                for i_roi in range(self._n_roi):
                    e_pp[i_channel_stim][i_roi] = postprocess_e(
                        e=e[i_channel_stim][i_roi],
                        e2=None,
                        dirvec=self._goal_dir[i_roi],
                        type=self.e_postproc[i_roi],
                    )

        # compute goal function value
        if isinstance(self.goal[0], types.FunctionType):
            goal_fun_value = self.goal[0](e_pp) # user provided goal function
        else:
            goal_fun_value = self.compute_goal(e_pp)

        logger.info(
            f"Goal ({self.goal}): {goal_fun_value:.3f} (n_sim: {self.n_sim}, n_test: {self.n_test})",
        )
        logger.info( "-" * len(parameters_str))

        return goal_fun_value

    def compute_goal(self, e):
        """
        Computes goal function value from postprocessed electric field

        Parameters
        ----------
        e: list of list of np.ndarrays of float 
            e-field magnitude, normal component, tangential component: [n_channel_stim][n_roi][n_roi_nodes]
            dirTI, maxTI: [1][n_roi][n_roi_nodes] 
           Post-processed electric fields from simulated simulation conditions and ROIs.

        Returns
        -------
        goal_fun_value : float
            Accumulated goal function value. The average is taken over all stimulation conditions and the weighted
            average is taken according to self.weights over the different goal functions of the ROIs.
        """
        n_effective_channels = len(e) # for TI, fields of the channels were previously combined --> only one "effective channel"
                    
        # focality based goal functions
        if "focality" in self.goal or "focality_inv" in self.goal:
            
            y = np.zeros((n_effective_channels))  # one function value per channel

            for i_channel_stim in range(n_effective_channels):
                ROCval = ROC(
                    e1=e[i_channel_stim][0],  # e-field in ROI
                    e2=e[i_channel_stim][1],  # e-field in non-ROI
                    threshold=self.threshold,
                    focal="focality" in self.goal,  # True for "focality", False for "focality_inv"
                )
                if "focality" in self.goal:
                    y[i_channel_stim] = -100 * ( np.sqrt(2) - ROCval )
                else:
                    y[i_channel_stim] = -100 * ROCval
                                
            # average over all effective channels   
            y_weighted_sum = np.mean(y)

        # Mean/Max/ etc. based goal functions
        else:

            y = np.zeros((n_effective_channels, self._n_roi))  # one goal function value per channel and roi
            
            for i_roi in range(self._n_roi):
                for i_channel_stim in range(n_effective_channels):
                    
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

            # average over all effective channels
            y_weighted_sum = np.mean(y, axis=0)
    
            # weight and sum the goal function values of the ROIs
            y_weighted_sum = np.sum(y_weighted_sum * self.weights)
        
        # if desired, track focality measures and goal function values
        for i_channel_stim in range(n_effective_channels):
            # goal fun value in rois
            self.goal_fun_value[i_channel_stim].append( y[i_channel_stim] )
                        
            if self.track_focality:
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

        # TODO: think this works only for one channel_stim right now (HDTES), test with 2 channel stim and adapt (KW)
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

        # TODO: same here I think it only works for one channel_stim (KW)
        # extract geometrical electrode parameters from optimal parameters and update electrode
        for i_channel_stim in range(self.n_channel_stim):
            if self.electrode[i_channel_stim]._any_free_geometry:
                n_free_parameters = np.sum(self.electrode[i_channel_stim]._free_geometry)
                self.electrode[i_channel_stim].set_geometrical_parameters_optimization(
                    electrode_pos_array[i_para : (i_para + n_free_parameters)]
                )
                i_para += n_free_parameters

        return electrode_pos

    def get_init_vals(self, bounds):
        """
        Determine initial values for optimization, guaranteeing a valid electrode position.

        Parameters
        ----------
        bounds : Bounds instance
            Lower and upper bounds of optimization problem

        Returns
        -------
        x0 : ndarray of float [n_para]
            Initial values
        """

        if self.seed is not None:
            np.random.seed(self.seed)

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

            if valid:
                return para_test_grid[i, :]
            
            logger.info( f"> electrode_pos_valid: {node_idx_dict[1]}")

        raise RuntimeError("failed to find valid electrode position.")

    def get_nodes_electrode(self, electrode_pos, plot=False):
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

        # determine node coords on skin surface of electrodes
        # afterwards, this is a good entry point to determine the distance to the outer boundary of the valid
        # skin region and between the electrodes in order to quantify the "distance" to the constrained region
        node_coords_list, node_idx_dict, electrode_pos_valid = self.get_node_coords_from_electrode_pos(
            electrode_pos=electrode_pos
        )

        if type(node_coords_list) is str:
            return node_coords_list, electrode_pos_valid

        # check if distance between electrodes is sufficient
        valid_str, electrode_pos_valid = check_electrode_distance(
            node_coords_list=node_coords_list,
            electrode_pos_valid=electrode_pos_valid,
            n_ele_free=self.n_ele_free,
            n_channel_stim=self.n_channel_stim,
            min_electrode_distance=self.min_electrode_distance
        )

        if type(valid_str) is str:
            return valid_str, electrode_pos_valid

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

        if plot:
            fn_electrode_txt = os.path.join(
                self._detailed_results_folder,
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

        return node_idx_dict

    def get_node_coords_from_electrode_pos(self, electrode_pos):
        """
        Determine coordinates of skin surface nodes given the electrode positions

        Parameter
        ---------
        electrode_pos : list of list of np.ndarray of float [2 or 3] of length [n_channel_stim][n_ele_free]
            Spherical coordinates (beta, lambda) and orientation angle (alpha) for each electrode array.
                      electrode array 1                        electrode array 2
            [ np.array([beta_1, lambda_1, alpha_1]),   np.array([beta_2, lambda_2, alpha_2]) ]

        Returns
        -------
        node_coords_list : list of ndarray [n_array_global_total][n_nodes x 3]
            List containing arrays of coordinates of nodes for each electrode array
        node_idx_dict : list of dict
            List [n_channel_stim] containing dicts with electrode channel IDs as keys and node indices.
        electrode_pos_valid : list of list containing ndarray of float of len (2) or (3) [n_channel][n_array]
            List of list containing the ellipsoid position and orientation parameters of the electrode arrays
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

        # determine electrode coords on ellipsoid
        ################################################################################################################
        # loop over channels
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

            # loop over arrays in channel
            for i_array, _electrode_array in enumerate(
                    self.electrode[i_channel_stim]._electrode_arrays
            ):
                # electrode center on ellipsoid (start point for geodesic)
                start[i_array, :], n_tmp[i_array, :] = self._ellipsoid.jacobi2cartesian(
                    coords=electrode_pos[i_channel_stim][i_array][:2],
                    return_normal=True,
                )

                # vector in direction of angle beta on ellipsoid (constant lambda)
                a[i_array, :] = (
                        self._ellipsoid.jacobi2cartesian(
                            coords=np.array(
                                [
                                    electrode_pos[i_channel_stim][i_array][0] - 1e-2,
                                    electrode_pos[i_channel_stim][i_array][1],
                                ]
                            )
                        )
                        - start[i_array, :]
                )

                # vector in direction of angle lambda on ellipsoid (constant beta)
                b[i_array, :] = (
                        self._ellipsoid.jacobi2cartesian(
                            coords=np.array(
                                [
                                    electrode_pos[i_channel_stim][i_array][0],
                                    electrode_pos[i_channel_stim][i_array][1] - 1e-2,
                                ]
                            )
                        )
                        - start[i_array, :]
                )
                a[i_array, :] /= np.linalg.norm(a[i_array, :])
                b[i_array, :] /= np.linalg.norm(b[i_array, :])

                # determine target point from start point (electrode center) into y-direction of electrode
                # according to selected orientation of array with respect to vector of constant lambda
                if len(electrode_pos[i_channel_stim][i_array]) > 2:
                    start_shifted_[i_array, :] = start[i_array, :] + (
                            1e-3
                            * (
                                    (a[i_array, :])
                                    * np.cos(electrode_pos[i_channel_stim][i_array][2])
                                    + (b[i_array, :])
                                    * np.sin(electrode_pos[i_channel_stim][i_array][2])
                            )
                    )
                else:
                    start_shifted_[i_array, :] = start[i_array, :] + 1e-3 * a[i_array, :]

                # vector of array y-direction
                cy[i_channel_stim][i_array, :] = (
                        start_shifted_[i_array, :] - start[i_array, :]
                )
                cy[i_channel_stim][i_array, :] /= np.linalg.norm(
                    cy[i_channel_stim][i_array, :]
                )

                # vector of array x-direction
                cx[i_channel_stim][i_array, :] = np.cross(
                    cy[i_channel_stim][i_array, :], -n_tmp[i_array, :]
                )
                cx[i_channel_stim][i_array, :] /= np.linalg.norm(
                    cx[i_channel_stim][i_array, :]
                )

                # distances of electrodes in array wrt to center of array for geodesic
                distance.append(_electrode_array.distance)

                # direction angles of single electrodes wrt array orientation for geodesic
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

            # determine electrode center on ellipsoid (run geodesics for each electrode)
            if not (distance == 0.0).all():
                if sys.platform == 'win32' and self._n_cpu != 1:
                    n_cpu = 1
                    logger.debug("Restricting geodesic destination calculations on Windows to one CPU core.")
                else:
                    n_cpu = self._n_cpu
                electrode_coords_eli_cart = self._ellipsoid.get_geodesic_destination(
                    start=start, distance=distance, alpha=alpha, n_steps=400, n_cpu=n_cpu
                )
            else:
                electrode_coords_eli_cart = start

            # normal vector at destinations (electrode locations on ellipsoid)
            n.append(self._ellipsoid.get_normal(coords=electrode_coords_eli_cart))

            # transform to ellipsoidal coordinates
            electrode_coords_eli_eli = self._ellipsoid.cartesian2ellipsoid(
                coords=electrode_coords_eli_cart
            )

            # project electrode coordinates from ellipsoid to subject and stencil out the nodes on the skin surface
            ############################################################################################################
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
                        node_idx_dict,
                        electrode_pos_valid,
                    )
                else:
                    electrode_pos_valid[i_channel_stim][i_array] = electrode_pos[
                        i_channel_stim
                    ][i_array]
                    # print("Electrode position: invalid (not all electrodes in valid skin region)")

                electrode_coords_subject[i_channel_stim] = np.vstack(tmp_arrays)

                # loop over electrodes and determine skin node indices of electrodes
                for _electrode in _electrode_array.electrodes:
                    if _electrode.type == "spherical":
                        mask, posmat = get_node_mask_spherical_electrode(
                            node_coords_skin=self._skin_surface.nodes,
                            ele_coords_center=electrode_coords_subject[i_channel_stim][i_ele, :],
                            ele_radius=_electrode.radius
                        )

                    elif _electrode.type == "rectangular":
                        mask, posmat = get_node_mask_rectangular_electrode(
                            node_coords_skin=self._skin_surface.nodes,
                            ele_coords_center=electrode_coords_subject[i_channel_stim][i_ele, :],
                            ele_dir_y=cy[i_channel_stim][i_array, :],
                            ele_dir_n=n[i_channel_stim][i_ele, :],
                            ele_len_x=_electrode.length_x,
                            ele_len_y=_electrode.length_y
                        )

                    else:
                        raise AssertionError(
                            "Electrodes have to be either 'spherical' or 'rectangular'"
                        )

                    # save position of electrode in subject space to posmat field
                    _electrode.posmat = posmat

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
                            node_idx_dict,
                            electrode_pos_valid

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

        return node_coords_list, node_idx_dict, electrode_pos_valid

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

        # assign surface nodes to electrode positions and estimate optimal currents
        node_idx_dict = self.get_nodes_electrode(electrode_pos=electrode_pos, plot=plot)

        if type(node_idx_dict[0]) is str:
            # some electrodes overlap --> return
            logger.info( node_idx_dict[0] )
            return None
            
        # perform one electric field calculation for every stimulation condition (one at a time is on)
        logger.info("Electrode positions valid")
        e = [[] for _ in range(self.n_channel_stim)]
        for i_channel_stim in range(self.n_channel_stim):
            
            if plot:
                # when using Dirichlet correction, online FEM will save 
                # electrode node coords and the optimized currents to this file
                fn_electrode_txt = os.path.join(
                    self._detailed_results_folder,
                    f"electrode_coords_nodes_subject_dirichlet{i_channel_stim}.txt",
                )
            else:
                fn_electrode_txt = None
                
            e[i_channel_stim] = self._ofem.update_field(electrode=self.electrode[i_channel_stim],
                                                         dirichlet_correction=self.electrode[i_channel_stim].dirichlet_correction,
                                                         fn_electrode_txt=fn_electrode_txt)[0]
            
            # store number of dirichlet iterations for convergence analysis
            if self.electrode[i_channel_stim].dirichlet_correction:
                self.n_iter_dirichlet_correction[i_channel_stim].append(
                    self._ofem.n_iter_dirichlet_correction
                )
                
        return e


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


def get_element_properties(roi):
    
    vol = None
    node_normals = None
    if roi.method == 'surface':
        if np.any(roi._mesh.elm.tetrahedra):
            raise ValueError("surface ROI must not contain tetrahedra")
        # return surface node areas and surface node normals
        node_normals = roi._mesh.nodes_normals().value[roi._mask]
        vol = np.mean(roi._mesh.nodes_volumes_or_areas().value[roi._mask])
    elif roi.method in ["volume", "volume_from_surface"]:
        # return element volumes
        vol = np.mean(roi._mesh.elements_volumes_and_areas().value[roi._mask])
        
    return vol, node_normals


def postprocess_e(e, e2=None, dirvec=None, type="magn"):
    """
    Post-processing electric field according to specified type.

    Parameters
    ----------
    e : np.ndarray of float [n_nodes x 3]
        Electric field components in query points (Ex, Ey, Ez)
    e2 : np.ndarray of float [n_nodes x 3]
        Electric field components in query points of second channel (Ex, Ey, Ez) for TI fields.
    dirvec : np.ndarray of float [n_nodes x 3]
        Normal vectors for normal and tangential e-field component calculation or general direction vectors the
        directional TI fields are calculated for. Can be either a single vector (1 x 3)
        that is applied to all positions or one vector per position (N x 3).
    type : str, optional, default: "magn"
        Type of postprocessing to apply:
        - "magn": electric field magnitude (default)
        - "normal": determine normal component (required surface normals in dirvec)
        - "tangential": determine tangential component (required surface normals in dirvec)
        - "max_TI": maximum envelope for TI fields
        - "dir_TI": directional sensitive maximum envelope for TI fields

    Returns
    -------
    e_pp: np.ndarray of float [n_nodes, ]
        Post-processed electric field in query points
    """
    assert e.shape[1] == 3, "Shape of electric field does not match the requirement [n_roi x 3]. " \
                            "Electric field components needed (Ex, Ey, Ez)!"

    if dirvec is not None:
        if dirvec.shape[0] == 1:
            dirvec = np.repeat(dirvec, e.shape[0], axis=0)

    if type in ["max_TI", "dir_TI"] and e2 is None:
        raise ValueError("Please provide second e-field to calculate TI field!")

    if type == "magn":
        e_pp = np.linalg.norm(e, axis=1)

    elif type == "normal":
        e_pp = -np.sum(e * dirvec, axis=1)

    elif type == "tangential":
        e_pp = np.sqrt(np.linalg.norm(e, axis=1) ** 2 - np.sum(e * dirvec, axis=1) ** 2)

    elif type == "max_TI":
        e_pp = get_maxTI(E1_org=e, E2_org=e2)

    elif type == "dir_TI":
        e_pp = get_dirTI(E1=e, E2=e2, dirvec_org=dirvec)

    else:
        raise NotImplementedError(f"Specified type for e-field post-processing '{type}' not implemented.")

    return e_pp


def write_visualization(folder_path, base_file_name, roi_list, results_list, e_postproc, goal_list):
    ''' uses the final FEM simulation results to create
        surface and/or volume meshes with results for each metric in e_postproc
        
        writes the meshes including geo and opt files to disc
        returns volume and surface mesh
    '''
    # e_postproc has a postprocessing option for each roi --> keep only unique set
    e_postproc = list(np.unique(e_postproc))
    
    for i in e_postproc:
        if i not in ["max_TI", "dir_TI", "magn", "normal", "tangential"]:
            raise ValueError(f"postprocessing option {i} unknown")
    
    n_roi = len(roi_list)
    n_results = len(results_list)
    
    # 1) initialize visualization
    roi_txt = ['ROI']
    if n_roi > 1:
        if "focality" in goal_list or "focality_inv" in goal_list:
            roi_txt = ['ROI','non-ROI']
        else:
            roi_txt = ['ROI_'+str(x) for x in range(n_roi)]
    
    results_txt = ['']
    if n_results > 1:
        results_txt = ['channel_'+str(x) for x in range(n_results)]
    
    roi_result_vis = RoiResultVisualization(
        roi_list,
        results_list,
        folder_path,
        base_file_name,
        roi_txt,
        results_txt
        )
    
    # 2a) append head mesh data to visualization
    headmesh_newdata = []
    m_head = None
    if roi_result_vis.has_head_mesh():
        m_head = roi_result_vis.get_head_mesh()
    
        if 'magn' in e_postproc and n_results > 1:
            # append magnE averaged across channels
            fieldnames = [results_txt[i]+'__magnE' for i in range(len(results_txt))]
            idx = [i for i, data in enumerate(m_head.elmdata) if data.field_name in fieldnames]
    
            data = np.zeros_like(m_head.elmdata[idx[0]].value)
            for i in idx:
                data += m_head.elmdata[i].value
            data /= n_results
            headmesh_newdata.append(m_head.add_element_field(data, 'average__magnE'))
    
        if 'max_TI' in e_postproc and n_results == 2:
            # append maxTI
            fieldnames = [results_txt[i]+'__E' for i in range(len(results_txt))]
            idx = [i for i, data in enumerate(m_head.elmdata) if data.field_name in fieldnames]
            if len(idx) != 2:
                raise ValueError('Exact two E fields needed to calculate maxTI')
    
            data = postprocess_e(m_head.elmdata[idx[0]].value, 
                                 e2=m_head.elmdata[idx[1]].value, 
                                 dirvec=None, 
                                 type='max_TI'
                                 )
            headmesh_newdata.append(m_head.add_element_field(data, 'max_TI'))
    
    # 2b) append surface data to visualization
    surfacemesh_newdata = []
    m_surf = None
    if roi_result_vis.has_surface_mesh():
        m_surf = roi_result_vis.get_surface_mesh()
        dirvec = m_surf.nodes_normals().value
        
        if 'magn' in e_postproc and n_results > 1:
            # append magnE averaged across channels
            fieldnames = [results_txt[i]+'__magnE' for i in range(len(results_txt))]
            idx = [i for i, data in enumerate(m_surf.nodedata) if data.field_name in fieldnames]
    
            data = np.zeros_like(m_surf.nodedata[idx[0]].value)
            for i in idx:
                data += m_surf.nodedata[i].value
            data /= n_results
            surfacemesh_newdata.append(m_surf.add_node_field(data, 'average__magnE'))
        
        for metric in ["normal", "tangential"]:
            if metric in e_postproc:
                # append normal and tangential components (and average across channels if n_results > 1)
                if n_results == 1:
                    fieldnames = ['E']
                else:
                    fieldnames = [results_txt[i]+'__E' for i in range(len(results_txt))]
                idx = [i for i, data in enumerate(m_surf.nodedata) if data.field_name in fieldnames]
                
                for idx_txt, i_nodedata in enumerate(idx):
                    data = postprocess_e(m_surf.nodedata[i_nodedata].value, 
                                         e2=None, 
                                         dirvec=dirvec, 
                                         type=metric
                                         )
                    surfacemesh_newdata.append(m_surf.add_node_field(data, results_txt[idx_txt]+'__'+metric)) 
                    if idx_txt == 0:
                        d_avg = np.copy(data)
                    else:
                        d_avg += data
        
                if len(idx) > 1:
                    d_avg /= len(idx)
                    surfacemesh_newdata.append(m_surf.add_node_field(d_avg, 'average__'+metric))
                    
        for metric in ['max_TI', 'dir_TI']:
            if metric in e_postproc and n_results == 2:
                # append maxTI and dirTI
                fieldnames = [results_txt[i]+'__E' for i in range(len(results_txt))]
                idx = [i for i, data in enumerate(m_surf.nodedata) if data.field_name in fieldnames]
                if len(idx) != 2:
                    raise ValueError('Exact two E fields needed to calculate maxTI and dirTI')
        
                data = postprocess_e(m_surf.nodedata[idx[0]].value, 
                                     e2=m_surf.nodedata[idx[1]].value, 
                                     dirvec=dirvec, 
                                     type=metric
                                     )
                surfacemesh_newdata.append(m_surf.add_node_field(data, metric))   
            
    
    # 3) create new meshes including geo and opt file data
    roi_result_vis.create_visualization()
    
    # 4) delete raw E field data to save some disc space
    fieldnames = [results_txt[i]+'__E' for i in range(len(results_txt))]
    fieldnames.append('E')
    for i in fieldnames:
        if roi_result_vis.has_head_mesh():
            roi_result_vis.remove_field_from_head_mesh(i)
        if roi_result_vis.has_surface_mesh():    
            roi_result_vis.remove_field_from_surface_mesh(i)
        
    # 5) update visualization settings in opt-files
    for _, view in roi_result_vis.head_mesh_data_name_to_gmsh_view.items():
        view.Visible = 0 # disable visiblity of all data fields
    
    if roi_result_vis.has_head_mesh():
        # get tissue list from volume ROIs
        tissues = []
        for i in roi_list:
            if i.method in ["volume", "volume_from_surface"] and i.tissues is not None:
                tissues += i.tissues
    
        for i in headmesh_newdata:
            view = roi_result_vis.head_mesh_data_name_to_gmsh_view[i.field_name]
            view2 = i.view_options(visible_tags=tissues)
            
            view.Visible = 1
            view.RangeType = view2.RangeType
            view.CenterGlyphs = view2.CenterGlyphs
            view.GlyphLocation = view2.GlyphLocation
            view.VectorType = view2.VectorType
            view.CustomMax = view2.CustomMax
            view.CustomMin = view2.CustomMin
            view.SaturateValues = view2.SaturateValues 
            view.RangeType = view2.RangeType 
            view.ShowScale = view2.ShowScale
            view.ColormapNumber = view2.ColormapNumber
            view.ColorTable = view2.ColorTable
            
        if len(headmesh_newdata) == 0:
            view = roi_result_vis.head_mesh_data_name_to_gmsh_view[m_head.elmdata[-1].field_name]
            view.Visible = 1
    
    for _, view in roi_result_vis.surface_mesh_data_name_to_gmsh_view.items():
        view.Visible = 0 # disable visiblity of all data fields    
        
    if roi_result_vis.has_surface_mesh():
        for i in surfacemesh_newdata:
            view = roi_result_vis.surface_mesh_data_name_to_gmsh_view[i.field_name]
            view2 = i.view_options()
            
            view.Visible = 1
            view.RangeType = view2.RangeType
            view.CenterGlyphs = view2.CenterGlyphs
            view.GlyphLocation = view2.GlyphLocation
            view.VectorType = view2.VectorType
            view.CustomMax = view2.CustomMax
            view.CustomMin = view2.CustomMin
            view.SaturateValues = view2.SaturateValues 
            view.RangeType = view2.RangeType 
            view.ShowScale = view2.ShowScale
            view.ColormapNumber = view2.ColormapNumber
            view.ColorTable = view2.ColorTable
            
        if len(surfacemesh_newdata) == 0:
            view = roi_result_vis.surface_mesh_data_name_to_gmsh_view[m_surf.nodedata[-1].field_name]
            view.Visible = 1
                            
    # 6) write out visualization meshes together with their geo and opt files
    fn_vis = roi_result_vis.write_visualization() 
    
    return fn_vis, m_head, m_surf


def make_summary_text(m_surf, m_head, tissues_m_head = [ElementTags.GM]):
    """generates a summary text with peak fields, focality, median fields in ROIs

    Args:
        m_surf (mesh_io.Msh): surface mesh created by write_visualization
        m_head (mesh_io.Msh): volume mesh created by write_visualization
        tissues_m_head (list of ElementTags): optional, standard: [ElementTags.GM]
                        tissue types included in evaluation of summary metrics for m_head
        
    Returns:
        summary text (str)
    """
    def summary_for_mesh(m):
        # get fields with final results
        result_field_names = {'max_TI', 'dir_TI', 'magnE', 'E__normal', 'E__tangential', 
                            'average__magnE', 'average__normal', 'average__tangential'}
        result_fields = m.field.keys() & result_field_names

        # get ROI fields
        roi_field_names = {'ROI', 'non-ROI'}
        roi_fields = m.field.keys() & roi_field_names
        for key in m.field:
            if key.startswith('ROI_'):
                roi_fields.add(key)

        # get medians per ROI
        arr_medians=np.ndarray((len(result_fields)+1,len(roi_fields)+1),dtype=object)
        arr_medians[0,0] = ''
        for idx_r, r in enumerate(roi_fields):
            arr_medians[0,idx_r+1] = r
            idx=np.argwhere(m.field[r].value>0)+1
            for idx_f, f in enumerate(result_fields):
                arr_medians[idx_f+1,0] = f
                arr_medians[idx_f+1,idx_r+1] = f'{m.field[f].get_percentiles(percentile=[50], roi=idx)[0]: .2e}' 

        # assemble summary text
        summary_text = '======================\n'
        summary_text += m.fields_summary(fields=result_fields, percentiles=[99.9, 50], focality_cutoffs=[75, 50])
        summary_text += '\nMedian fields per ROI\n----------------------\n'
        summary_text += mesh_io._format_table(arr_medians)
        return summary_text

    summary_text = '\n'
    hlpStr = 'head_mesh (included tissues: '+' '.join(tissue_names[k] for k in tissues_m_head)+')'
    m_hlp = None
    if m_head is not None:
        m_hlp = m_head.crop_mesh(tags=tissues_m_head)
    for m, name in [[m_surf, 'surface_mesh'],
                    [m_hlp, hlpStr]]:
        if m is not None:
            summary_text += f'Results for {name}:\n'
            summary_text += summary_for_mesh(m) + '\n\n'
    return summary_text


def get_node_mask_spherical_electrode(node_coords_skin, ele_coords_center, ele_radius):
    """
    Return mask of skin nodes included in the spherical electrode.

    Parameters
    ----------
    node_coords_skin : ndarray of float (n_node_coords_skin, 3)
        Coordinates of skin nodes
    ele_coords_center : ndarray of float (3)
        Coordinates (x,y,z) of electrode center
    ele_radius : float
        Radius of electrode

    Returns
    -------
    mask : ndarray of bool (node_coords_skin.shape[0]])
        Node mask of spherical electrode for node_coords_skin
    posmat : ndarray of float (4, 4)
        Matrix containing the position and orientation of the electrode
    """

    # mask with a sphere
    mask = (
            np.linalg.norm(
                node_coords_skin - ele_coords_center,
                axis=1,
            )
            < ele_radius
    )

    posmat = np.eye(4)
    posmat[:3, 3] = ele_coords_center

    return mask, posmat


def get_node_mask_rectangular_electrode(node_coords_skin, ele_coords_center,
                                        ele_dir_y, ele_dir_n,
                                        ele_len_x, ele_len_y
                                        ):
    """
    Return mask of skin nodes included in a rectangular electrode.

    Parameters
    ----------
    node_coords_skin : ndarray of float (n_node_coords_skin, 3)
        Coordinates of skin nodes
    ele_coords_center : ndarray of float (3)
        Coordinates (x,y,z) of electrode center
    ele_dir_y : ndarray of float (3)
        Direction vector of electrode y-direction
    ele_dir_n : ndarray of float (3)
        Direction vector of electrode pointing into normal direction wrt skin surface
    ele_len_x : float
        Dimension of electrode in x-direction
    ele_len_y : float
        Dimension of electrode in y-direction

    Returns
    -------
    mask : np.array of bool [node_coords_skin.shape[0]]
        Node mask of spherical electrode for node_coords_skin
    posmat :
    """
    n = ele_dir_n
    cy = ele_dir_y
    center = ele_coords_center

    cx_local = np.cross(n, cy)

    # rotate skin nodes to normalized electrode space
    rotmat = np.array(
        [
            [
                cx_local[0],
                cy[0],
                n[0],
            ],
            [
                cx_local[1],
                cy[1],
                n[1],
            ],
            [
                cx_local[2],
                cy[2],
                n[2],
            ],
        ]
    )

    # save position of electrode in subject space to posmat field
    posmat = np.vstack(
        (
            np.hstack((rotmat, center[:, np.newaxis])),
            np.array([0, 0, 0, 1]),
        )
    )

    skin_nodes_rotated = (node_coords_skin - center) @ rotmat

    # mask with a box
    mask_x = np.logical_and(
        skin_nodes_rotated[:, 0] > -ele_len_x / 2,
        skin_nodes_rotated[:, 0] < +ele_len_x / 2,
    )
    mask_y = np.logical_and(
        skin_nodes_rotated[:, 1] > -ele_len_y / 2,
        skin_nodes_rotated[:, 1] < +ele_len_y / 2,
    )
    mask_z = np.logical_and(
        skin_nodes_rotated[:, 2] > -30,
        skin_nodes_rotated[:, 2] < +30,
    )
    mask = np.logical_and(np.logical_and(mask_x, mask_y), mask_z)

    return mask, posmat


def check_electrode_distance(node_coords_list, electrode_pos_valid, n_ele_free, n_channel_stim,
                             min_electrode_distance):
    """
    Checks if the distance between electrode positions is sufficient.

    Parameter
    ---------
    node_coords_list : list of ndarray [n_array_global_total][n_nodes x 3]
        List containing arrays of coordinates of nodes for each electrode array
    electrode_pos_valid : list of list containing ndarray of float of len (2) or (3) [n_channel][n_array]
        List of list containing the ellipsoid position and orientation parameters of the electrode arrays
    n_ele_free : list of int [n_channel_stim]
        Number of free moveable electrode arrays per channel
    n_channel_stim : int
        Number of stimulation channels
    min_electrode_distance : float
        Minimum allowable distance between electrodes

    Returns
    -------
    str or bool:
        Error message in case of invalid electrode position
    electrode:pos_valid : list of list of ndarrays of float of len (2) or (3) [n_channel][n_array]
        Updated electrode_pos_valid list. Contains None at identified invalid electrode positions.
    """
    invalid = False
    i_array_global_lst = np.hstack(
        [
            np.arange(n_ele_free[i_channel_stim])
            for i_channel_stim in range(n_channel_stim)
        ]
    ).astype(int)
    i_channel_stim_global_lst = np.hstack(
        [
            i_channel_stim * np.ones(n_ele_free[i_channel_stim])
            for i_channel_stim in range(n_channel_stim)
        ]
    ).astype(int)
    if min_electrode_distance is not None and min_electrode_distance >= 0:
        i_array_test_start = 1
        # start with first array and test if all node coords are too close to other arrays
        for i_array_global in range(np.sum(n_ele_free)):
            for node_coord in node_coords_list[i_array_global]:
                for i_array_test in range(
                        i_array_test_start, np.sum(n_ele_free)
                ):
                    # calculate euclidean distance between node coords
                    min_dist = np.min(
                        np.linalg.norm(
                            node_coords_list[i_array_test] - node_coord, axis=1
                        )
                    )
                    # stop testing if an electrode is too close
                    if min_dist <= min_electrode_distance:
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
    else:
        return True, electrode_pos_valid

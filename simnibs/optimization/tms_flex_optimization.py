import itertools
import logging
import os
import subprocess
import tempfile
import time
import glob
from typing import Callable, Optional, get_type_hints
import numpy.typing as npt

import scipy.optimize as opt
import scipy
import numpy as np

from simnibs.mesh_tools import mesh_io
from simnibs.simulation.onlinefem import FemTargetPointCloud
from simnibs.utils.region_of_interest import (
    RegionOfInterest,
)
from simnibs.simulation.tms_coil.tms_coil import TmsCoil
from simnibs.simulation.tms_coil.tms_coil_deformation import (
    TmsCoilDeformation,
    TmsCoilDeformationRange,
    TmsCoilRotation,
    TmsCoilTranslation,
)
from simnibs.simulation.tms_coil.tms_coil_element import DipoleElements, TmsCoilElements
from simnibs.utils.matlab_read import (
    matlab_field_to_list,
    matlab_struct_to_dict,
    matlab_sub_struct_to_matlab_struct,
    remove_None,
    try_to_read_matlab_field,
)

from ..simulation.sim_struct import POSITION
from ..utils.simnibs_logger import logger
from ..utils.file_finder import SubjectFiles
from ..utils.mesh_element_properties import ElementTags
from simnibs import Msh

from simnibs.utils import file_finder

from simnibs.utils import simnibs_logger
from .. import __version__
from simnibs.simulation import sim_struct


class TmsFlexOptimization:
    """Class that defines a flexible TMS optimization

    Parameters
    ------------------------
    matlab_struct: (optional) scipy.io.loadmat()
        matlab structure
    """

    date: str
    """Date when the optimization struct was initiated"""
    time_str: str
    """Time when the optimization struct was initiated"""

    fnamecoil: str
    """Path to the coil file (example: "path/to/coil")"""
    pos: POSITION
    """The initial coil position from where the optimization is started"""

    fnamehead: str
    """Path to the head mesh (example: "path/to/msh")"""
    subpath: str
    """Path to the m2m folder of the subject (example: "path/to/m2m")"""
    path_optimization: str
    """Path to the output folder where the result of the optimization is saved (example: "path/to/output")"""
    eeg_cap: str
    """Path to a eeg cap file (example: "path/to/csv")"""
    run_simulation: bool
    """Weather to run a full simulation at the optimized position after the optimization"""

    method: str
    """The method of optimization {"distance", "emag"}"""
    roi: RegionOfInterest
    """The region of interest in which the e field is simulated, method = "emag" """
    global_translation_ranges: list
    """Ranges that describe how far the coil is allowed to move in x,y,z direction [[min(x), max(x)],[min(y), max(y)], [min(z), max(z)]] (example: [[-1, 1], [-2, 2], [-10, 10]] | [-1, 1] [0, 0], [0, 0]])"""
    global_rotation_ranges: list
    """Ranges that describe how far the coil is allowed to rotate around the x,y,z axis [[min(x), max(x)],[min(y), max(y)], [min(z), max(z)]] (example: [[-1, 1], [-2, 2], [-10, 10]] | [-1, 1] [0, 0], [0, 0]])"""

    dither_skip: int
    """How many voxel positions should be skipped when creating the coil volume representation.
        Used to speed up the optimization. When set to 0, no dithering will be applied. (example: 2 | 0), default 6"""
    fem_evaluation_cutoff: float
    """If the penalty from the intersection and self intersection is greater than this cutoff value, the fem will not be evaluated to save time.
        Set to np.inf to always evaluate the fem (example: 500 | np.inf), default 1000"""

    run_global_optimization: bool
    """Weather to run the global optimization, the global optimization will always run first"""
    run_local_optimization: bool
    """Weather to run the local optimization, the local optimization will always run second"""

    direct_args: dict
    """Settings for the scipy direct global optimization"""
    l_bfgs_b_args: dict
    """Settings for the scipy L-BFGS-B local optimization"""

    def __init__(self, matlab_struct=None):
        # Date when the session was initiated
        self.date: str = time.strftime("%Y-%m-%d %H:%M:%S")
        self.time_str: str = time.strftime("%Y%m%d-%H%M%S")

        self.fnamecoil: str = None
        self.pos: POSITION = None

        self.fnamehead: str = None
        self.subpath: str = None
        self.path_optimization: str = None
        self.eeg_cap: str = None
        self.run_simulation: bool = True

        self.method: str = None
        self.roi: RegionOfInterest = None
        self.global_translation_ranges: list = None
        self.global_rotation_ranges: list = None

        self.dither_skip: int = 6
        self.fem_evaluation_cutoff: float = 1000

        self.run_global_optimization = True
        self.run_local_optimization = True

        self.direct_args = {}
        self.l_bfgs_b_args = {}

        self._prepared = False
        self._log_handlers = []

        if matlab_struct:
            self.read_mat_struct(matlab_struct)

    def add_position(self, position=None):
        """Adds the position to the current optimization

        Parameters
        -----
        position: POSITION (Optional)
            Position structure defining the coil position

        Returns
        ------
        position: POSITION
            POSITION structure defining the coil position
        """
        if position is None:
            position = POSITION()

        self.pos = position
        return position

    def add_region_of_interest(self, roi=None):
        """Adds the region of interest to the current optimization

        Parameters
        -----
        roi: RegionOfInterest (Optional)
            Region of interest structure defining the region of interest

        Returns
        ------
        roi: RegionOfInterest
            RegionOfInterest structure defining the region of interest
        """
        if roi is None:
            roi = RegionOfInterest()

        self.roi = roi
        return roi

    def _prepare(self):
        """Prepares session for simulations
        relative paths are made absolute,
        empty fields are set to default values,
        check if required fields exist
        """

        if self.fnamehead is not None:
            self.fnamehead = os.path.abspath(os.path.expanduser(self.fnamehead))
            if not os.path.isfile(self.fnamehead):
                raise IOError("Cannot locate head mesh file: %s" % self.fnamehead)
        else:
            self.subpath = os.path.abspath(os.path.expanduser(self.subpath))
            if not os.path.isdir(self.subpath):
                raise IOError("Cannot locate subjects m2m folder: %s" % self.subpath)

        sub_files = SubjectFiles(self.fnamehead, self.subpath)
        self.fnamehead = sub_files.fnamehead
        self.subpath = sub_files.subpath

        if not os.path.isdir(self.subpath):
            logger.warning("Cannot locate subjects m2m folder")
            logger.warning("some postprocessing options might fail")
            self.subpath = None

        if self.eeg_cap is None:
            self.eeg_cap = sub_files.eeg_cap_1010

        logger.info(f"Head Mesh: {self.fnamehead}")
        logger.info(f"Subject Path: {self.subpath}")
        self.path_optimization = os.path.abspath(
            os.path.expanduser(self.path_optimization)
        )
        logger.info(f"Optimization Folder: {self.path_optimization}")

        if self.pos is None:
            raise ValueError("No initial position specified")

        self.pos.eeg_cap = self.eeg_cap
        self.pos._prepare()

        try:
            fnamecoil = os.path.expanduser(self.fnamecoil)
        except TypeError:
            raise TypeError(f"fnamecoil ({self.fnamecoil}) is not a string")
        if os.path.isfile(fnamecoil):
            self.fnamecoil = fnamecoil
        else:
            fnamecoil = os.path.join(file_finder.coil_models, self.fnamecoil)
            if os.path.isfile(fnamecoil):
                self.fnamecoil = fnamecoil
            else:
                raise IOError(f"Could not find coil file: {self.fnamecoil}")
        logger.info(f"Coil Path: {self.fnamecoil}")
        self._coil = TmsCoil.from_file(fnamecoil)

        if self.method is None:
            raise ValueError("No method specified")

        self._mesh = Msh(fn=self.fnamehead)

        if self.method == "emag":
            if self.roi is not None:
                if (
                    self.subpath is not None
                    and self.roi.subpath is None
                    and self.roi.mesh is None
                ):
                    self.roi.subpath = self.subpath
                elif self.roi.subpath is None and self.roi.mesh is None:
                    self.roi.mesh = self._mesh
                self._roi = FemTargetPointCloud(self._mesh, self.roi.get_nodes())
            else:
                raise ValueError("No ROI specified")

        if self.global_translation_ranges is not None:
            self._global_translation_ranges = np.array(self.global_translation_ranges)

        if self.global_rotation_ranges is not None:
            self._global_rotation_ranges = np.array(self.global_rotation_ranges)

        self.pos.calc_matsimnibs(self._mesh)

        if "vol_tol" not in self.direct_args:
            if self.method == "distance":
                self.direct_args["vol_tol"] = 1e-19
            else:
                self.direct_args["vol_tol"] = 1e-16

        self._prepared = True

    @property
    def type(self):
        return self.__class__.__name__

    def _set_logger(self, fname_prefix="simnibs_optimization"):
        """
        Set-up logger to write to a file

        Parameters
        ----------
        fname_prefix: str, optional
            Prefix of log-file
        """
        if not os.path.isdir(self.path_optimization):
            os.mkdir(self.path_optimization)
        log_fn = os.path.join(
            self.path_optimization, f"{fname_prefix}_{self.time_str}.log"
        )
        fh = logging.FileHandler(log_fn, mode="w")
        formatter = logging.Formatter(
            f"[ %(name)s {__version__} - %(asctime)s - %(process)d ]%(levelname)s: %(message)s"
        )
        fh.setFormatter(formatter)
        fh.setLevel(logging.DEBUG)
        logger = logging.getLogger("simnibs")
        logger.addHandler(fh)
        self._log_handlers += [fh]
        simnibs_logger.register_excepthook(logger)

    def _finish_logger(self):
        logger = logging.getLogger("simnibs")
        [logger.removeHandler(lh) for lh in self._log_handlers]
        self._log_handlers = []
        simnibs_logger.unregister_excepthook()

    def run(self):
        """Runs the tms flex optimization"""
        self._set_logger()
        dir_name = os.path.abspath(os.path.expanduser(self.path_optimization))
        if os.path.isdir(dir_name):
            g = glob.glob(os.path.join(dir_name, "optimized_simnibs_simulation*.mat"))
            if len(g) > 0:
                raise IOError(
                    f"{os.linesep}Found already existing simulation results in directory."
                    f"{os.linesep}Please run the simulation in a new directory or delete"
                    f" the simnibs_simulation*.mat files from the folder : {dir_name}"
                )
            logger.info(f"Running optimization in the directory: {dir_name}")
        else:
            logger.info(f"Running optimization on new directory: {dir_name}")
            os.makedirs(dir_name)

        self._prepare()

        if self.method == "emag" and self._roi.n_center == 0:
            raise ValueError("The region of interest contains no positions")

        logger.info(f"Optimization options:{os.linesep}{self.to_str_formatted()}")

        logger.info(f"Running optimization ({self.method})")
        # Run simulations
        if self.method == "distance":
            initial_cost, optimized_cost, opt_matsimnibs, direct = optimize_distance(
                self._coil,
                self._mesh,
                self.pos.matsimnibs,
                self._global_translation_ranges,
                self._global_rotation_ranges,
                self.dither_skip,
                self.run_global_optimization,
                self.run_local_optimization,
                self.direct_args,
                self.l_bfgs_b_args,
            )
        elif self.method == "emag":
            initial_cost, optimized_cost, opt_matsimnibs, optimized_e_mag, direct = (
                optimize_e_mag(
                    self._coil,
                    self._mesh,
                    self._roi,
                    self.pos.matsimnibs,
                    self._global_translation_ranges,
                    self._global_rotation_ranges,
                    self.dither_skip,
                    self.fem_evaluation_cutoff,
                    self.run_global_optimization,
                    self.run_local_optimization,
                    self.direct_args,
                    self.l_bfgs_b_args,
                )
            )
        else:
            raise ValueError("method should be 'distance' or 'emag'")

        logger.info(
            f"Optimization result:{os.linesep}{direct}{os.linesep}"
            f"Initial cost: {initial_cost}{os.linesep}"
            f"Optimized cost: {optimized_cost}"
        )
        if self.method == "emag":
            logger.info(f"Optimized mean E-field magnitude in ROI: {optimized_e_mag}")

        logger.info(f"Matsimnibs result:{os.linesep}{opt_matsimnibs}")

        logger.info(f"Creating Visualizations")

        skin_mesh = self._mesh.crop_mesh(tags=[ElementTags.SCALP_TH_SURFACE])
        mesh_name = os.path.splitext(os.path.basename(self._mesh.fn))[0]
        coil_name = os.path.splitext(os.path.basename(self.fnamecoil))[0]
        fn_geo = os.path.join(
            self.path_optimization, f"{mesh_name}_{coil_name}_optimization.geo"
        )
        fn_out = os.path.join(
            self.path_optimization, f"{mesh_name}_{coil_name}_optimization.msh"
        )

        skin_mesh.write(fn_out)
        v = self._mesh.view(visible_tags=[ElementTags.SCALP_TH_SURFACE.value])
        self._coil.append_simulation_visualization(
            v, fn_geo, skin_mesh, self.pos.matsimnibs, visibility=0, infix="-initial"
        )
        self._coil.append_simulation_visualization(
            v, fn_geo, skin_mesh, opt_matsimnibs, infix="-optimized"
        )

        fn_optimized_coil = os.path.join(
            self.path_optimization, f"{mesh_name}_{coil_name}_optimized.tcd"
        )
        self._coil.freeze_deformations().write(fn_optimized_coil)

        v.add_merge(fn_geo)
        v.write_opt(fn_out)

        logger.info(f"Creating Simulation")
        S = sim_struct.SESSION()
        S.subpath = self.subpath
        S.fnamehead = self.fnamehead

        S.pathfem = os.path.join(self.path_optimization, "tms_simulation")

        ## Define the TMS simulation
        tms = S.add_tmslist()
        tms.fnamecoil = fn_optimized_coil

        # Define the coil position
        pos = tms.add_position(self.pos)
        pos.matsimnibs = opt_matsimnibs

        fn_sim = os.path.join(
            self.path_optimization, f"optimized_simnibs_simulation_{self.time_str}.mat"
        )
        sim_struct.save_matlab_sim_struct(S, fn_sim)

        if self.run_simulation:
            S.run()

    def to_mat(self):
        """Makes a dictionary for saving a matlab structure with scipy.io.savemat()

        Returns
        --------------------
        dict
            Dictionary for usage with scipy.io.savemat
        """
        # Generate dict from instance variables (excluding variables starting with _ or __)
        mat = {
            key: remove_None(value)
            for key, value in self.__dict__.items()
            if not key.startswith("__")
            and not key.startswith("_")
            and not callable(value)
            and not callable(getattr(value, "__get__", None))
        }

        # Add class name as type (type is protected in python so it cannot be a instance variable)
        mat["type"] = "TmsFlexOptimization"

        # Add all instance variables that are classes
        # Manually or by calling their to_mat function
        if self.roi is not None:
            mat["roi"] = self.roi.to_mat()

        if self.pos is not None:
            pos_dt = np.dtype(
                [
                    ("type", "O"),
                    ("name", "O"),
                    ("date", "O"),
                    ("matsimnibs", "O"),
                    ("didt", "O"),
                    ("fnamefem", "O"),
                    ("centre", "O"),
                    ("pos_ydir", "O"),
                    ("distance", "O"),
                ]
            )

            pos_array = np.array(
                [
                    (
                        "POSITION",
                        remove_None(self.pos.name),
                        remove_None(self.pos.date),
                        remove_None(self.pos.matsimnibs),
                        remove_None(self.pos.didt),
                        remove_None(self.pos.fnamefem),
                        remove_None(self.pos.centre),
                        remove_None(self.pos.pos_ydir),
                        remove_None(self.pos.distance),
                    )
                ],
                dtype=pos_dt,
            )
            mat["pos"] = pos_array

        return mat

    @classmethod
    def read_mat_struct(self, mat):
        """Reads parameters from matlab structure

        Parameters
        ----------
        mat: scipy.io.loadmat
            Loaded matlab structure
        """
        # Load all instance variables from the mat file
        # Datatypes come from the type hints of the instance variables
        types = get_type_hints(TmsFlexOptimization)
        for key, value in self.__dict__.items():
            if key in mat:
                setattr(
                    self, key, try_to_read_matlab_field(mat, key, types[key], value)
                )

        if "direct_args" in mat:
            self.direct_args = matlab_struct_to_dict(mat["direct_args"])

        if "l_bfgs_b_args" in mat:
            self.l_bfgs_b_args = matlab_struct_to_dict(mat["l_bfgs_b_args"])

        # Load all 2d arrays manually (try_to_read_matlab_field will load only array_name[0])
        self.global_translation_ranges = matlab_field_to_list(
            mat, "global_translation_ranges", 2
        )
        self.global_rotation_ranges = matlab_field_to_list(
            mat, "global_rotation_ranges", 2
        )

        # Load all instance variables that are classes
        if "pos" in mat and mat["pos"].size > 0:
            self.pos = POSITION(mat["pos"][0][0])

        if "roi" in mat and mat["roi"].size > 0:
            self.roi = RegionOfInterest(matlab_sub_struct_to_matlab_struct(mat["roi"]))

        return self

    def to_str_formatted(self):
        import pprint

        # Generate dict from instance variables (excluding variables starting with _ or __ or None or len == 0)
        tms_flex_dict = {
            key: value
            for key, value in self.__dict__.items()
            if not key.startswith("__")
            and not key.startswith("_")
            and not callable(value)
            and not callable(getattr(value, "__get__", None))
            and value is not None
            and (not hasattr(value, "__len__") or len(value) > 0)
        }

        # Add all instance variables that are classes manually
        tms_flex_dict["pos"] = {
            key: value
            for key, value in self.pos.__dict__.items()
            if not key.startswith("__")
            and not key.startswith("_")
            and not callable(value)
            and not callable(getattr(value, "__get__", None))
            and value is not None
            and (not hasattr(value, "__len__") or len(value) > 0)
        }

        return pprint.pformat(tms_flex_dict, indent=4)


def _get_fast_distance_score(
    distance_function: Callable,
    elements: list[TmsCoilElements],
    affine: npt.NDArray[np.float_],
) -> float:
    """Calculates the mean absolute distance from the min distance points of the coil using the distance function.
    If no min distance points are present the node positions of the casings will be used.

    Parameters
    ----------
    distance_function : Callable
        A distance function calculating the distance between input points and a target
    elements : list[TmsCoilElements]
        The coil elements to be used
    affine : npt.NDArray[np.float_]
        The affine transformation from coil to world space

    Returns
    -------
    float
        The mean absolute distance from the min distance points (coil casing nodes) using the distance function
    """
    casing_points = []
    min_distance_points = []

    for coil_element in elements:
        if coil_element.casing is not None:
            element_casing_points = coil_element.get_casing_coordinates(affine, True)
            if len(element_casing_points[0]) > 0:
                casing_points.append(element_casing_points[0])
            if len(element_casing_points[1]) > 0:
                min_distance_points.append(element_casing_points[1])

    if len(casing_points) > 0:
        casing_points = np.concatenate(casing_points, axis=0)
    if len(min_distance_points) > 0:
        min_distance_points = np.concatenate(min_distance_points, axis=0)

    min_distance_points = (
        min_distance_points if len(min_distance_points) > 0 else casing_points
    )

    return np.mean(np.abs(distance_function(min_distance_points)))


def _get_fast_intersection_penalty(
    element_voxel_volumes: dict[TmsCoilElements, npt.NDArray[np.bool_]],
    element_voxel_indexes: dict[TmsCoilElements, npt.NDArray[np.int_]],
    element_voxel_dither_factors: dict[TmsCoilElements, npt.NDArray[np.float_]],
    element_voxel_affines: dict[TmsCoilElements, npt.NDArray[np.float_]],
    target_voxel_distance: npt.NDArray[np.float_],
    target_voxel_affine: npt.NDArray[np.float_],
    self_intersection_elements: list[TmsCoilElements],
    affine: npt.NDArray[np.float_],
    order: int = 3,
) -> tuple[float, float]:
    """Evaluates how far the element voxel volumes intersect with the target voxel volume measured in mm^3 (volume) * mm (depth).
    Evaluates the amount of self intersection based on the self intersection element groups measured in mm^3.


    Parameters
    ----------
    element_voxel_volume : dict[TmsCoilElements, npt.NDArray[np.bool_]]
        The voxel volume (True inside, False outside) of each coil element
    element_voxel_indexes : dict[TmsCoilElements, npt.NDArray[np.int_]]
        The indexes of the inside voxels for each coil element
    element_voxel_dither_factors : dict[TmsCoilElements, npt.NDArray[np.float_]]
        The dither scaling factor for each voxel for each coil element
    element_voxel_affine : dict[TmsCoilElements, npt.NDArray[np.float_]]
        The affine transformations from world to voxel space for each coil element
    target_voxel_distance : npt.NDArray[np.float_]
        The voxel distance field of the target
    target_voxel_affine : npt.NDArray[np.float_]
        The affine transformations from voxel to world space for the target
    self_intersection_elements : list[TmsCoilElements]
        The groups of coil elements that need to be checked for self intersection with each other
    affine : npt.NDArray[np.float_]
        The affine transformation from coil to world space
    order : int
        The order of interpolation on the target_voxel_distance

    Returns
    -------
    weighted_target_intersection_quibic_mm : float
        Depth weighted volume intersection between the coil volume and the target volume in mm^3 (volume) * mm (depth)
    self_intersection_quibic_mm : float
        Sum of the volume intersection between the coil element inside the self intersection groups in mm^3
    """
    element_affines = {
        element: element.get_combined_transformation(affine)
        for element in element_voxel_volumes
    }
    element_inv_affines = {
        element: np.linalg.inv(element_affines[element])
        for element in element_voxel_volumes
    }

    target_voxel_inv_affine = np.linalg.inv(target_voxel_affine)
    weighted_target_intersection_quibic_mm = 0
    for element in element_voxel_volumes:
        dither_factors = element_voxel_dither_factors[element]
        indexes_in_vox1 = element_voxel_indexes[element]
        element_voxel_affine = element_voxel_affines[element]
        vox_to_vox_affine = (
            target_voxel_inv_affine @ element_affines[element] @ element_voxel_affine
        )
        indexes_in_vox2 = (
            vox_to_vox_affine[:3, :3] @ indexes_in_vox1.T
            + vox_to_vox_affine[:3, 3, None]
        )
        intersections = scipy.ndimage.map_coordinates(
            target_voxel_distance, indexes_in_vox2, order=order, prefilter=False
        )
        weighted_target_intersection_quibic_mm += np.sum(intersections * dither_factors)

    self_intersection_quibic_mm = 0
    for intersection_group in self_intersection_elements:
        for intersection_pair in itertools.combinations(intersection_group, 2):
            dither_factors = element_voxel_dither_factors[intersection_pair[0]]
            vox_to_vox_affine = (
                np.linalg.inv(element_voxel_affines[intersection_pair[1]])
                @ element_inv_affines[intersection_pair[1]]
                @ element_affines[intersection_pair[0]]
                @ element_voxel_affines[intersection_pair[0]]
            )
            indexes_in_vox1 = element_voxel_indexes[intersection_pair[0]]
            indexes_in_vox2 = (
                vox_to_vox_affine[:3, :3] @ indexes_in_vox1.T
                + vox_to_vox_affine[:3, 3, None]
            )
            voxel_volume1 = element_voxel_volumes[intersection_pair[1]]
            intersections = scipy.ndimage.map_coordinates(
                voxel_volume1, indexes_in_vox2, order=order, prefilter=False
            )
            self_intersection_quibic_mm += np.sum(intersections * dither_factors)

    return (
        weighted_target_intersection_quibic_mm,
        self_intersection_quibic_mm,
    )


def _prepare_skin_surface(mesh: Msh) -> Msh:
    """Preparing the skin surface by removing the internal air and closing the surface

    Parameters
    ----------
    mesh : Msh
        The head mesh where the skin surface should be created for

    Returns
    -------
    Msh
        The closed skin surface without internal air

    Raises
    ------
    subprocess.CalledProcessError
        If the meshfix call returned an error
    """
    try:
        skin = mesh.relabel_internal_air().crop_mesh(ElementTags.SCALP_TH_SURFACE)
        with tempfile.TemporaryDirectory() as tmp_dirname:
            mesh_io.write_off(skin, os.path.join(tmp_dirname, "mymesh.off"))
            cmd = [
                file_finder.path2bin("meshfix"),
                os.path.join(tmp_dirname, "mymesh.off"),
                "-a",
                "2.0",
                "-o",
                os.path.join(tmp_dirname, "mymesh.off"),
            ]
            p = subprocess.Popen(
                cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
            p.wait()
            cmd = [
                file_finder.path2bin("meshfix"),
                os.path.join(tmp_dirname, "mymesh.off"),
                "-a",
                "2.0",
                "-u",
                "1",
                "--vertexDensity",
                "0.2",
                "-o",
                os.path.join(tmp_dirname, "mymesh.off"),
            ]
            p = subprocess.Popen(
                cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
            )
            p.wait()
            if p.returncode != 0:
                raise subprocess.CalledProcessError(p.returncode, cmd)
            skin = mesh_io.read_off(os.path.join(tmp_dirname, "mymesh.off"))
    except Exception:
        skin = mesh.crop_mesh(ElementTags.SCALP_TH_SURFACE)

    return skin


def optimize_distance(
    coil: TmsCoil,
    head_mesh: Msh,
    affine: npt.NDArray[np.float_],
    coil_translation_ranges: Optional[npt.NDArray[np.float_]] = None,
    coil_rotation_ranges: Optional[npt.NDArray[np.float_]] = None,
    dither_skip: int = 0,
    global_optimization: bool = True,
    local_optimization: bool = True,
    direct_args: dict | None = None,
    l_bfgs_b_args: dict | None = None,
) -> tuple[float, float, npt.NDArray[np.float_], list]:
    """Optimizes the deformations of the coil elements as well as the global transformation to minimize the distance between the optimization_surface
    and the min distance points (if not present, the coil casing points) while preventing intersections of the
    optimization_surface and the intersect points (if not present, the coil casing points)

    Parameters
    ----------
    coil : TmsCoil
        The coil used in the optimization
    head_mesh : Msh
        The head mesh used in the TMS simulation and the head mesh where the scalp surface is used for coil head intersection
    affine : npt.NDArray[np.float_]
        The affine transformation that is applied to the coil
    coil_translation_ranges : Optional[npt.NDArray[np.float_]], optional
        If the global coil position is supposed to be optimized as well, these ranges in the format
        [[min(x), max(x)],[min(y), max(y)], [min(z), max(z)]] are used
        and the updated affine coil transformation is returned, by default None
    coil_translation_ranges : Optional[npt.NDArray[np.float_]], optional
        If the global coil rotation is supposed to be optimized as well, these ranges in the format
        [[min(x), max(x)],[min(y), max(y)], [min(z), max(z)]] are used
        and the updated affine coil transformation is returned, by default None
    dither_skip : int, optional
        How many voxel positions should be skipped when creating the coil volume representation.
        Used to speed up the optimization. When set to 0, no dithering will be applied, by default 0

    Returns
    -------
    initial_cost
        The initial cost
    optimized_cost
        The cost after optimization
    result_affine
        The affine matrix. If coil_translation_ranges is None than it's the input affine,
        otherwise it is the optimized affine.
    opt_results
        The results of the optimizations

    Raises
    ------
    ValueError
        If the coil has no deformations to optimize
    ValueError
        If the coil has no coil casing or no min distance points
    """

    coil_deformation_ranges = coil.get_deformation_ranges()

    if (
        len(coil_deformation_ranges) == 0
        and coil_translation_ranges is None
        and coil_rotation_ranges is None
    ):
        raise ValueError(
            "The coil has no deformations to optimize the coil element positions with."
        )

    if not np.any([np.any(arr) for arr in coil.get_casing_coordinates()]):
        raise ValueError("The coil has no coil casing or min_distance points.")

    if direct_args is None:
        direct_args = {}
    direct_args.setdefault("eps", 1e-4)
    direct_args.setdefault("maxiter", 1000)
    direct_args.setdefault("locally_biased", False)
    direct_args.setdefault("vol_tol", 1e-19)
    direct_args.setdefault("len_tol", 1e-6)

    if l_bfgs_b_args is None:
        l_bfgs_b_args = {}
    if "options" not in l_bfgs_b_args:
        l_bfgs_b_args["options"] = {}
    l_bfgs_b_args["options"].setdefault("maxls", 100)

    global_deformations = add_global_deformations(
        coil, coil_rotation_ranges, coil_translation_ranges
    )

    optimization_surface = _prepare_skin_surface(head_mesh)

    (
        element_voxel_volume,
        element_voxel_indexes,
        element_voxel_dither_factors,
        element_voxel_affine,
        self_intersection_elements,
    ) = get_voxel_volume(coil, global_deformations, dither_skip=dither_skip)
    (
        target_distance_function,
        target_voxel_distance,
        target_voxel_affine,
        cost_surface_tree,
    ) = optimization_surface.get_min_distance_on_grid()
    target_voxel_distance_inside = np.minimum(target_voxel_distance, 0) * -1

    coil_deformation_ranges = coil.get_deformation_ranges()
    initial_deformation_settings = np.array(
        [coil_deformation.current for coil_deformation in coil_deformation_ranges]
    )

    best_f = np.inf
    best_x = None

    def cost_f_x0_w(x):
        for coil_deformation, deformation_setting in zip(coil_deformation_ranges, x):
            coil_deformation.current = deformation_setting

        (
            intersection_penalty,
            self_intersection_penalty,
        ) = _get_fast_intersection_penalty(
            element_voxel_volume,
            element_voxel_indexes,
            element_voxel_dither_factors,
            element_voxel_affine,
            target_voxel_distance_inside,
            target_voxel_affine,
            self_intersection_elements,
            affine,
        )
        distance_penalty = _get_fast_distance_score(
            target_distance_function, element_voxel_volume.keys(), affine
        )

        f = distance_penalty + (intersection_penalty + self_intersection_penalty)

        nonlocal best_f
        if f < best_f:
            nonlocal best_x
            best_x = x.copy()
            best_f = f

        return f

    initial_cost = cost_f_x0_w(initial_deformation_settings)

    optimization_ranges = [deform.range for deform in coil_deformation_ranges]

    opt_results = []
    if global_optimization:
        direct_args.setdefault("maxfun", len(optimization_ranges) * 2000)
        direct = opt.direct(cost_f_x0_w, bounds=optimization_ranges, **direct_args)
        if direct.fun > best_f:
            direct.message += f""" Using better solution encountered during optimization.
            opt.fun: {direct.fun}, best_fun: {best_f}
            opt.x: {direct.x}, best_x: {best_x}"""
            direct.x = best_x
            direct.fun = best_f

        best_deformation_settings = direct.x

        for coil_deformation, deformation_setting in zip(
            coil_deformation_ranges, best_deformation_settings
        ):
            coil_deformation.current = deformation_setting
        opt_results.append(direct)

    if local_optimization:
        initial_deformation_settings = np.array(
            [coil_deformation.current for coil_deformation in coil_deformation_ranges]
        )
        local_opt = opt.minimize(
            cost_f_x0_w,
            x0=initial_deformation_settings,
            bounds=[deform.range for deform in coil_deformation_ranges],
            method="L-BFGS-B",
            **l_bfgs_b_args,
        )

        if local_opt.fun > best_f:
            local_opt.message += f""" Using better solution encountered during optimization.
            opt.fun: {local_opt.fun}, best_fun: {best_f}
            opt.x: {local_opt.x}, best_x: {best_x}"""
            local_opt.x = best_x
            local_opt.fun = best_f

        best_deformation_settings = local_opt.x

        for coil_deformation, deformation_setting in zip(
            coil_deformation_ranges, best_deformation_settings
        ):
            coil_deformation.current = deformation_setting

        opt_results.append(local_opt)

    optimized_cost = cost_f_x0_w(best_deformation_settings)

    result_affine = np.eye(4)
    if len(global_deformations) > 0:
        for global_deformation in global_deformations:
            for coil_element in coil.elements:
                coil_element.deformations.remove(global_deformation)
            result_affine = global_deformation.as_matrix() @ result_affine
    result_affine = affine.astype(float) @ result_affine

    return initial_cost, optimized_cost, result_affine, opt_results


def get_voxel_volume(
    coil: TmsCoil, global_deformations: list[TmsCoilDeformation], dither_skip: int = 0
) -> tuple[
    dict[TmsCoilElements, npt.NDArray[np.bool_]],
    dict[TmsCoilElements, npt.NDArray[np.int_]],
    dict[TmsCoilElements, npt.NDArray[np.float_]],
    dict[TmsCoilElements, npt.NDArray[np.float_]],
    list[TmsCoilElements],
]:
    """Generates voxel volume information about the coil.

    Parameters
    ----------
    global_deformations : list[TmsCoilDeformation]
        Global deformations used in the coil
    dither_skip : int, optional
        How many voxel positions should be skipped when creating the coil volume representation.
        Used to speed up the optimization. When set to 0, no dithering will be applied, by default 0

    Returns
    -------
    element_voxel_volume : dict[TmsCoilElements, npt.NDArray[np.bool_]]
        The voxel volume (True inside, False outside) of each coil element
    element_voxel_indexes : dict[TmsCoilElements, npt.NDArray[np.int_]]
        The indexes of the inside voxels for each coil element (interior voxel coordinates first, then edge coordinates)
    element_voxel_dither_factors : dict[TmsCoilElements, npt.NDArray[float]]
        The spacial factor that describes how many voxel are described by each interior dithered voxel for each voxel inside the volume for each coil element
    element_voxel_affine : dict[TmsCoilElements, npt.NDArray[np.float_]]
        The affine transformations from world to voxel space for each coil element
    self_intersection_elements : list[TmsCoilElements]
        The groups of coil elements that need to be checked for self intersection with each other
    """
    element_voxel_volume = {}
    element_voxel_indexes = {}
    element_voxel_dither_factors = {}
    element_voxel_affine = {}
    if coil.casing is not None:
        base_element = DipoleElements(
            None,
            np.zeros((1, 3)),
            np.zeros((1, 3)),
            casing=coil.casing,
            deformations=global_deformations,
        )
        (
            element_voxel_volume[base_element],
            element_voxel_affine[base_element],
            element_voxel_indexes[base_element],
            element_voxel_dither_factors[base_element],
            _,
        ) = coil.casing.mesh.get_voxel_volume(dither_skip=dither_skip)
    self_intersection_elements = []
    for self_intersection_group in coil.self_intersection_test:
        self_intersection_elements.append([])
        for self_intersection_index in self_intersection_group:
            if self_intersection_index == 0:
                self_intersection_elements[-1].append(base_element)
            else:
                self_intersection_elements[-1].append(
                    coil.elements[self_intersection_index - 1]
                )
    for element in coil.elements:
        if element.casing is not None:
            (
                element_voxel_volume[element],
                element_voxel_affine[element],
                element_voxel_indexes[element],
                element_voxel_dither_factors[element],
                _,
            ) = element.casing.mesh.get_voxel_volume(dither_skip=dither_skip)

    return (
        element_voxel_volume,
        element_voxel_indexes,
        element_voxel_dither_factors,
        element_voxel_affine,
        self_intersection_elements,
    )


def add_global_deformations(
    coil,
    coil_rotation_ranges: Optional[npt.NDArray[np.float_]] = None,
    coil_translation_ranges: Optional[npt.NDArray[np.float_]] = None,
) -> list[TmsCoilDeformation]:
    """Adds deformations to the coil. The deformations are added to all coil elements so that they are global.

    Parameters
    ----------
    coil_rotation_ranges : Optional[npt.NDArray[np.float_]], optional
        Adds global rotations to the coil, the format is [[min(x), max(x)],[min(y), max(y)], [min(z), max(z)]], by default None
    coil_translation_ranges : Optional[npt.NDArray[np.float_]], optional
        Adds global deformations to the coil, the format is [[min(x), max(x)],[min(y), max(y)], [min(z), max(z)]], by default None

    Returns
    -------
    global_deformations : list[TmsCoilDeformation]
        Returns a list of the added global deformations
    """
    global_deformations = []
    if coil_rotation_ranges is not None:
        if coil_rotation_ranges[0, 0] != coil_rotation_ranges[0, 1]:
            global_deformations.append(
                TmsCoilRotation(
                    TmsCoilDeformationRange(
                        0, (coil_rotation_ranges[0, 0], coil_rotation_ranges[0, 1])
                    ),
                    [0, 0, 0],
                    [1, 0, 0],
                )
            )

        if coil_rotation_ranges[1, 0] != coil_rotation_ranges[1, 1]:
            global_deformations.append(
                TmsCoilRotation(
                    TmsCoilDeformationRange(
                        0, (coil_rotation_ranges[1, 0], coil_rotation_ranges[1, 1])
                    ),
                    [0, 0, 0],
                    [0, 1, 0],
                )
            )

        if coil_rotation_ranges[2, 0] != coil_rotation_ranges[2, 1]:
            global_deformations.append(
                TmsCoilRotation(
                    TmsCoilDeformationRange(
                        0, (coil_rotation_ranges[2, 0], coil_rotation_ranges[2, 1])
                    ),
                    [0, 0, 0],
                    [0, 0, 1],
                )
            )
    if coil_translation_ranges is not None:
        if coil_translation_ranges[0, 0] != coil_translation_ranges[0, 1]:
            global_deformations.append(
                TmsCoilTranslation(
                    TmsCoilDeformationRange(
                        0,
                        (
                            coil_translation_ranges[0, 0],
                            coil_translation_ranges[0, 1],
                        ),
                    ),
                    0,
                )
            )
        if coil_translation_ranges[1, 0] != coil_translation_ranges[1, 1]:
            global_deformations.append(
                TmsCoilTranslation(
                    TmsCoilDeformationRange(
                        0,
                        (
                            coil_translation_ranges[1, 0],
                            coil_translation_ranges[1, 1],
                        ),
                    ),
                    1,
                )
            )

        if coil_translation_ranges[2, 0] != coil_translation_ranges[2, 1]:
            global_deformations.append(
                TmsCoilTranslation(
                    TmsCoilDeformationRange(
                        0,
                        (
                            coil_translation_ranges[2, 0],
                            coil_translation_ranges[2, 1],
                        ),
                    ),
                    2,
                )
            )
    for global_deformation in global_deformations:
        for coil_element in coil.elements:
            coil_element.deformations.append(global_deformation)
    return global_deformations


def optimize_e_mag(
    coil: TmsCoil,
    head_mesh: Msh,
    roi: FemTargetPointCloud,
    affine: npt.NDArray[np.float_],
    coil_translation_ranges: Optional[npt.NDArray[np.float_]] = None,
    coil_rotation_ranges: Optional[npt.NDArray[np.float_]] = None,
    dither_skip: int = 0,
    fem_evaluation_cutoff: float = 1000,
    global_optimization: bool = True,
    local_optimization: bool = True,
    direct_args: dict | None = None,
    l_bfgs_b_args: dict | None = None,
    debug: bool = False,
) -> tuple[float, float, npt.NDArray[np.float_], npt.NDArray[np.float_], list]:
    """Optimizes the deformations of the coil elements as well as the global transformation to maximize the mean e-field magnitude in the ROI while preventing intersections of the
    scalp surface and the coil casing

    Parameters
    ----------
    coil : TmsCoil
        The coil used in the optimization
    head_mesh : Msh
        The head mesh used in the TMS simulation and the head mesh where the scalp surface is used for coil head intersection
    roi: RegionOfInterest
        Region of interest for the calculation of the e field
    affine : npt.NDArray[np.float_]
        The affine transformation that is applied to the coil
    coil_translation_ranges : Optional[npt.NDArray[np.float_]], optional
        If the global coil position is supposed to be optimized as well, these ranges in the format
        [[min(x), max(x)],[min(y), max(y)], [min(z), max(z)]] are used
        and the updated affine coil transformation is returned, by default None
    coil_rotation_ranges : Optional[npt.NDArray[np.float_]], optional
        If the global coil rotation is supposed to be optimized as well, these ranges in the format
        [[min(x), max(x)],[min(y), max(y)], [min(z), max(z)]] are used
        and the updated affine coil transformation is returned, by default None
    dither_skip : int, optional
        How many voxel positions should be skipped when creating the coil volume representation.
        Used to speed up the optimization. When set to 0, no dithering will be applied, by default 0
    fem_evaluation_cutoff : float, optional
        If the penalty from the intersection and self intersection is greater than this cutoff value, the fem will not be evaluated to save time.
        Set to np.inf to always evaluate the fem, by default 1000

    Returns
    -------
    initial_cost : float
        The initial cost
    optimized_cost : float
        The cost after the optimization
    result_affine : npt.NDArray[np.float_]
        The optimized affine matrix
    optimized_e_mag : npt.NDArray[np.float_]
        The e field magnitude in the roi elements after the optimization
    opt_results : list
        The results of the optimizations

    Raises
    ------
    ValueError
        If the coil has no deformations to optimize
    ValueError
        If the coil has no casing
    """
    from simnibs.simulation.onlinefem import OnlineFEM

    coil_sampled = coil.as_sampled()
    coil_deformation_ranges = coil_sampled.get_deformation_ranges()

    if (
        len(coil_deformation_ranges) == 0
        and coil_translation_ranges is None
        and coil_rotation_ranges is None
    ):
        raise ValueError(
            "The coil has no deformations to optimize the coil element positions with."
        )

    if direct_args is None:
        direct_args = {}
    direct_args.setdefault("eps", 1e-4)
    direct_args.setdefault("maxiter", 1000)
    direct_args.setdefault("locally_biased", True)
    direct_args.setdefault("vol_tol", 1e-16)
    direct_args.setdefault("len_tol", 1e-6)

    if l_bfgs_b_args is None:
        l_bfgs_b_args = {}
    if "options" not in l_bfgs_b_args:
        l_bfgs_b_args["options"] = {}
    l_bfgs_b_args["options"].setdefault("maxls", 100)

    global_deformations = add_global_deformations(
        coil_sampled, coil_rotation_ranges, coil_translation_ranges
    )
    optimization_surface = _prepare_skin_surface(head_mesh)

    (
        element_voxel_volume,
        element_voxel_indexes,
        element_voxel_dither_factors,
        element_voxel_affine,
        self_intersection_elements,
    ) = get_voxel_volume(coil_sampled, global_deformations, dither_skip=dither_skip)
    (
        target_distance_function,
        target_voxel_distance,
        target_voxel_affine,
        cost_surface_tree,
    ) = optimization_surface.get_min_distance_on_grid()
    target_voxel_distance_inside = np.minimum(target_voxel_distance, 0) * -1

    coil_deformation_ranges = coil_sampled.get_deformation_ranges()

    if len(element_voxel_volume) == 0:
        raise ValueError(
            "The coil has no coil casing to be used for coil head intersection tests."
        )

    fem = OnlineFEM(
        head_mesh, "TMS", roi, coil=coil_sampled, dataType=[0], useElements=False
    )

    initial_deformation_settings = np.array(
        [coil_deformation.current for coil_deformation in coil_deformation_ranges]
    )
    if debug:
        tracking_deformations = []
        fs = []

    best_f = np.inf
    best_x = None

    def cost_f_x0_w(x):
        if debug:
            deformations = []

        for coil_deformation, deformation_setting in zip(coil_deformation_ranges, x):
            if debug:
                deformations.append(deformation_setting)
            coil_deformation.current = deformation_setting
        (
            intersection_penalty,
            self_intersection_penalty,
        ) = _get_fast_intersection_penalty(
            element_voxel_volume,
            element_voxel_indexes,
            element_voxel_dither_factors,
            element_voxel_affine,
            target_voxel_distance_inside,
            target_voxel_affine,
            self_intersection_elements,
            affine,
        )

        penalty = intersection_penalty + self_intersection_penalty
        if penalty > fem_evaluation_cutoff:
            return penalty

        roi_e_field = fem.update_field(matsimnibs=affine)

        f = penalty - 100 * np.mean(roi_e_field)

        nonlocal best_f
        if f < best_f:
            nonlocal best_x
            best_x = x.copy()
            best_f = f

        if debug:
            tracking_deformations.append(deformations)
            fs.append(f)
        return f

    initial_cost = cost_f_x0_w(initial_deformation_settings)

    opt_results = []
    if global_optimization:
        direct = opt.direct(
            cost_f_x0_w,
            bounds=[deform.range for deform in coil_deformation_ranges],
            **direct_args,
        )
        if direct.fun > best_f:
            direct.message += f""" Using better solution encountered during optimization.
            opt.fun: {direct.fun}, best_fun: {best_f}
            opt.x: {direct.x}, best_x: {best_x}"""
            direct.x = best_x
            direct.fun = best_f
        best_deformation_settings = direct.x

        for sampled_coil_deformation, deformation_setting in zip(
            coil_deformation_ranges, best_deformation_settings
        ):
            sampled_coil_deformation.current = deformation_setting
        opt_results.append(direct)

    if local_optimization:
        initial_deformation_settings = np.array(
            [coil_deformation.current for coil_deformation in coil_deformation_ranges]
        )
        local_opt = opt.minimize(
            cost_f_x0_w,
            x0=initial_deformation_settings,
            bounds=[deform.range for deform in coil_deformation_ranges],
            method="L-BFGS-B",
            **l_bfgs_b_args,
        )

        if local_opt.fun > best_f:
            local_opt.message += f""" Using better solution encountered during optimization.
            opt.fun: {local_opt.fun}, best_fun: {best_f}
            opt.x: {local_opt.x}, best_x: {best_x}"""
            local_opt.x = best_x
            local_opt.fun = best_f

        best_deformation_settings = local_opt.x

        for coil_deformation, deformation_setting in zip(
            coil_deformation_ranges, best_deformation_settings
        ):
            coil_deformation.current = deformation_setting

        opt_results.append(local_opt)

    optimized_cost = cost_f_x0_w(best_deformation_settings)
    optimized_e_mag = np.ravel(fem.update_field(matsimnibs=affine))

    result_affine = np.eye(4)
    if len(global_deformations) > 0:
        for global_deformation in global_deformations:
            for coil_element in coil_sampled.elements:
                coil_element.deformations.remove(global_deformation)
            result_affine = global_deformation.as_matrix() @ result_affine
    result_affine = affine.astype(float) @ result_affine

    for sampled_coil_deformation, coil_deformation in zip(
        coil_sampled.get_deformation_ranges(), coil.get_deformation_ranges()
    ):
        coil_deformation.current = sampled_coil_deformation.current
    if debug:
        return (
            initial_cost,
            optimized_cost,
            result_affine,
            optimized_e_mag,
            opt_results,
            tracking_deformations,
            fs,
        )
    else:
        return initial_cost, optimized_cost, result_affine, optimized_e_mag, opt_results

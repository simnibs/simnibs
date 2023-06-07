import json
from typing import Optional
import jsonschema
import numpy as np
import numpy.typing as npt
import nibabel as nib
import scipy.optimize as opt

from simnibs.mesh_tools.mesh_io import Elements, Msh, NodeData, Nodes
from simnibs.simulation.coil.coil_constants import CoilElementTag
from simnibs.simulation.coil.coil_deformation import (
    CoilDeformation,
    CoilRotation,
    CoilTranslation,
)
from simnibs.simulation.coil.tcd_element import TcdElement
from simnibs.simulation.coil.tms_stimulator import TMSStimulator, TMSWaveform
from simnibs.utils import file_finder

from .coil_model import CoilModel
from .coil_element import (
    CoilDipoles,
    CoilElement,
    CoilLineElements,
    CoilLinePoints,
    CoilSampledGridElements,
    DirectionalCoilElement,
    PositionalCoilElement,
)
from ... import __version__


class Coil(TcdElement):
    def __init__(
        self,
        name: Optional[str],
        brand: Optional[str],
        version: Optional[str],
        limits: Optional[npt.NDArray[np.float_]],
        resolution: Optional[npt.NDArray[np.float_]],
        coil_casing: Optional[CoilModel],
        coil_elements: list[CoilElement],
    ):
        self.name = name
        self.brand = brand
        self.version = version
        self.limits = limits
        self.resolution = resolution
        self.coil_casing = coil_casing
        self.coil_elements = coil_elements

        self.deformations: list[CoilDeformation] = []
        for coil_element in self.coil_elements:
            for coil_deformation in coil_element.deformations:
                if coil_deformation not in self.deformations:
                    self.deformations.append(coil_deformation)

    def get_da_dt(
        self,
        msh: Msh,
        coil_affine: npt.NDArray[np.float_],
        di_dt: float,
        eps: float = 1e-3,
    ):
        target_positions = msh.nodes.node_coord
        A = np.zeros_like(target_positions)
        for coil_element in self.coil_elements:
            A += coil_element.get_da_dt(target_positions, coil_affine, di_dt, eps)

        return NodeData(A)

    def get_a_field(
        self,
        points: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ):
        a_field = np.zeros_like(points)
        for coil_element in self.coil_elements:
            a_field += coil_element.get_a_field(points, coil_affine, eps)

        return a_field

    def get_coil_mesh(
        self,
        coil_affine: Optional[npt.NDArray[np.float_]] = None,
        apply_deformation: bool = True,
        include_casing: bool = True,
        include_optimization_points: bool = True,
        include_coil_elements: bool = True,
    ) -> Msh:
        if coil_affine is None:
            coil_affine = np.eye(4)

        coil_msh = Msh()
        if self.coil_casing is not None and include_casing:
            coil_msh = coil_msh.join_mesh(
                self.coil_casing.get_transformed_mesh(
                    coil_affine, include_optimization_points, 0
                )
            )
        for i, coil_element in enumerate(self.coil_elements):
            coil_msh = coil_msh.join_mesh(
                coil_element.get_mesh(
                    coil_affine,
                    apply_deformation,
                    include_casing,
                    include_optimization_points,
                    include_coil_elements,
                    (i + 1) * 100,
                )
            )
        return coil_msh

    def get_deformed_casing_coordinates(
        self, affine: Optional[npt.NDArray[np.float_]] = None
    ) -> tuple[npt.NDArray[np.float_], npt.NDArray[np.float_], npt.NDArray[np.float_]]:
        if affine is None:
            affine = np.eye(4)

        casing_points = (
            [self.coil_casing.get_points(affine)]
            if self.coil_casing is not None
            else []
        )
        min_distance_points = (
            [self.coil_casing.get_min_distance_points(affine)]
            if self.coil_casing is not None
            else []
        )
        intersect_points = (
            [self.coil_casing.get_intersect_points(affine)]
            if self.coil_casing is not None
            else []
        )

        for coil_element in self.coil_elements:
            if coil_element.element_casing is not None:
                element_casing_points = coil_element.get_deformed_casing_coordinates(
                    affine
                )
                casing_points.append(element_casing_points[0])
                min_distance_points.append(element_casing_points[1])
                intersect_points.append(element_casing_points[2])

        casing_points = np.concatenate(casing_points, axis=0)
        min_distance_points = np.concatenate(min_distance_points, axis=0)
        intersect_points = np.concatenate(intersect_points, axis=0)

        return casing_points, min_distance_points, intersect_points

    @staticmethod
    def _add_logo(mesh: Msh):
        """adds the simnibs logo to the coil surface"""

        msh_logo = Msh(fn=file_finder.templates.simnibs_logo)

        # 'simnibs' has tag 1, '3' has tag 2, '4' has tag 3
        # renumber tags, because they will be converted to color:
        # 0 gray, 1 red, 2 lightblue, 3 blue
        major_version = __version__.split(".")[0]
        if major_version == "3":
            msh_logo = msh_logo.crop_mesh(tags=[1, 2])
            msh_logo.elm.tag1[msh_logo.elm.tag1 == 2] = 3  # version in blue
        elif major_version == "4":
            msh_logo = msh_logo.crop_mesh(tags=[1, 3])
        else:
            msh_logo = msh_logo.crop_mesh(tags=1)
        msh_logo.elm.tag1[msh_logo.elm.tag1 == 1] = 2  # 'simnibs' in light blue

        # center logo in xy-plane, mirror at yz-plane and scale
        bbox_coil = np.vstack([np.min(mesh.nodes[:], 0), np.max(mesh.nodes[:], 0)])
        bbox_logo = np.vstack(
            [np.min(msh_logo.nodes[:], 0), np.max(msh_logo.nodes[:], 0)]
        )
        bbox_ratio = np.squeeze(np.diff(bbox_logo, axis=0) / np.diff(bbox_coil, axis=0))
        bbox_ratio = max(bbox_ratio[0:2])  # maximal size ratio in xy plane

        msh_logo.nodes.node_coord[:, 0:2] -= np.mean(bbox_logo[:, 0:2], axis=0)
        msh_logo.nodes.node_coord[:, 0] = -msh_logo.nodes.node_coord[:, 0]
        msh_logo.nodes.node_coord[:, 0:2] *= 1 / (4 * bbox_ratio)

        # shift logo along negative z to the top side of coil
        msh_logo.nodes.node_coord[:, 2] += bbox_coil[0, 2] - bbox_logo[0, 2] - 5

        mesh = mesh.join_mesh(msh_logo)
        return mesh

    @classmethod
    def from_file(
        cls,
        fn: str,
        name: Optional[str] = None,
        brand: Optional[str] = None,
        version: Optional[str] = None,
    ):
        if fn.endswith(".tcd"):
            return Coil.from_tcd_file(fn)
        elif fn.endswith(".ccd"):
            return Coil.from_ccd(fn, version)
        elif fn.endswith(".nii.gz") or fn.endswith(".nii"):
            return Coil.from_nifti(fn, name, brand, version)
        else:
            return Coil.from_tcd(fn)

    def write(self, fn: str):
        self.write_tcd(fn)

    @classmethod
    def from_ccd(
        cls,
        fn: str,
        version: Optional[str] = None,
    ):
        with open(fn, "r") as f:
            header = f.readline()

        pairs = header.replace("\n", "").split(";")[1:]

        header_dict = {}
        for pair in pairs:
            key, value = pair.split("=")

            if value == "none":
                value = None
            elif "." in value:
                value = float(value)
            elif value.isdigit():
                value = int(value)
            elif "," in value:
                value = np.fromstring(value, dtype=int, sep=",")

            header_dict[key] = value

        bb = []
        for dim in ("x", "y", "z"):
            a = header_dict.get(dim)
            if a is None:
                bb.append(None)
            else:
                if len(a) < 2:
                    bb.append((-np.abs(a[0]), np.abs(a[0])))
                else:
                    bb.append(a)

        res = []
        a = header_dict.get("resolution")
        if a is None:
            res.append(None)
        else:
            a = np.atleast_1d(a)
            if len(a) < 3:
                for i in range(len(a), 3):
                    a = np.concatenate((a, (a[i - 1],)))
            res = a

        ccd_file = np.atleast_2d(np.loadtxt(fn, skiprows=2))

        dipole_positions = ccd_file[:, 0:3]
        dipole_moments = ccd_file[:, 3:]

        if "dIdtmax" in header_dict.keys():
            stimulator = TMSStimulator(
                header_dict.get("stimulator"), None, header_dict["dIdtmax"], None
            )
        else:
            stimulator = None

        coil_elements = [
            CoilDipoles(None, None, [], dipole_positions, dipole_moments, stimulator)
        ]

        return cls(
            header_dict.get("coilname"),
            header_dict.get("brand"),
            version,
            np.array(bb),
            np.array(res),
            None,
            coil_elements,
        )
    
    def to_tcd(self) -> dict:
        tcd_coil_models = []
        coil_models = []
        if self.coil_casing is not None:
            tcd_coil_models.append(self.coil_casing.to_tcd())
            coil_models.append(self.coil_casing)

        tcd_deforms = []

        tcd_stimulators = []
        stimulators = []

        for deformation in self.deformations:
            tcd_deforms.append(deformation.to_tcd())

        tcd_coil_elements = []
        for coil_element in self.coil_elements:
            if (
                coil_element.element_casing not in coil_models
                and coil_element.element_casing is not None
            ):
                coil_models.append(coil_element.element_casing)
                tcd_coil_models.append(coil_element.element_casing.to_tcd())

            if (
                coil_element.stimulator not in stimulators
                and coil_element.stimulator is not None
            ):
                stimulators.append(coil_element.stimulator)
                tcd_stimulators.append(coil_element.stimulator.to_tcd())

            tcd_coil_elements.append(
                coil_element.to_tcd(stimulators, coil_models, self.deformations)
            )

        tcd_coil = {}
        if self.name is not None:
            tcd_coil["name"] = self.name
        if self.brand is not None:
            tcd_coil["brand"] = self.brand
        if self.version is not None:
            tcd_coil["version"] = self.version
        if self.limits is not None:
            tcd_coil["limits"] = self.limits.tolist()
        if self.resolution is not None:
            tcd_coil["resolution"] = self.resolution.tolist()
        if self.coil_casing is not None:
            tcd_coil["coilCasing"] = coil_models.index(self.coil_casing)
        tcd_coil["coilElementList"] = tcd_coil_elements
        if len(tcd_stimulators) > 0:
            tcd_coil["stimulatorList"] = tcd_stimulators
        if len(tcd_deforms) > 0:
            tcd_coil["deformList"] = tcd_deforms
        if len(tcd_coil_models) > 0:
            tcd_coil["coilModels"] = tcd_coil_models

        return tcd_coil

    @classmethod
    def from_tcd(cls, coil: dict, validate=True):
        if validate:
            with open(file_finder.templates.tcd_json_schema, "r") as fid:
                tcd_schema = json.loads(fid.read())

            try:
                jsonschema.validate(coil, tcd_schema)
            except jsonschema.ValidationError as e:
                instance = str(e.instance)
                e.instance = (
                    instance
                    if len(instance) < 900
                    else f"{instance[:400]} ... {instance[-400:]}"
                )
                raise e

        coil_models = []
        for coil_model in coil.get("coilModels", []):
            coil_models.append(CoilModel.from_tcd(coil_model))

        deformations = []
        for deform in coil.get("deformList", []):
            deformations.append(CoilDeformation.from_tcd(deform))

        stimulators = []
        for stimulator in coil.get("stimulatorList", []):
            stimulators.append(TMSStimulator.from_tcd(stimulator))

        coil_elements = []
        for coil_element in coil["coilElementList"]:
            coil_elements.append(
                CoilElement.from_tcd(
                    coil_element, stimulators, coil_models, deformations
                )
            )

        coil_casing = (
            None if coil.get("coilCasing") is None else coil_models[coil["coilCasing"]]
        )

        return cls(
            coil.get("name"),
            coil.get("brand"),
            coil.get("version"),
            None if coil.get("limits") is None else np.array(coil["limits"]),
            None if coil.get("resolution") is None else np.array(coil["resolution"]),
            coil_casing,
            coil_elements,
        )

    @classmethod
    def from_tcd_file(cls, fn: str, validate=True):
        with open(fn, "r") as fid:
            coil = json.loads(fid.read())

        return cls.from_tcd(coil, validate)

    def write_tcd(self, fn: str):
        with open(fn, "w") as json_file:
            json.dump(self.to_tcd(), json_file, indent=4)

    @classmethod
    def from_nifti(
        cls,
        fn: str,
        name: Optional[str] = None,
        brand: Optional[str] = None,
        version: Optional[str] = None,
    ):
        nifti = nib.load(fn)
        data = nifti.get_fdata()
        affine = nifti.affine

        resolution = np.array(
            [
                affine[0][0],
                affine[1][1],
                affine[2][2],
            ]
        )

        limits = np.array(
            [
                [affine[0][3], data.shape[0] * resolution[0] + affine[0][3]],
                [affine[1][3], data.shape[1] * resolution[1] + affine[1][3]],
                [affine[2][3], data.shape[2] * resolution[2] + affine[2][3]],
            ]
        )

        coil_elements = [CoilSampledGridElements(name, None, [], data, affine, None)]

        return cls(name, brand, version, limits, resolution, None, coil_elements)

    def write_nifti(
        self,
        fn: str,
        limits: Optional[npt.NDArray[np.float_]] = None,
        resolution: Optional[npt.NDArray[np.float_]] = None,
    ):
        limits = limits or self.limits
        if limits is None:
            raise ValueError("Limits needs to be set")
        resolution = resolution or self.resolution
        if resolution is None:
            raise ValueError("resolution needs to be set")

        dims = [
            int((max_ - min_) // res) for [min_, max_], res in zip(limits, resolution)
        ]

        dx = np.spacing(1e4)
        x = np.linspace(limits[0][0], limits[0][1] - resolution[0] + dx, dims[0])
        y = np.linspace(limits[1][0], limits[1][1] - resolution[0] + dx, dims[1])
        z = np.linspace(limits[2][0], limits[2][1] - resolution[0] + dx, dims[2])
        points = np.array(np.meshgrid(x, y, z, indexing="ij"))
        points = points.reshape((3, -1)).T

        data = self.get_a_field(points, np.eye(4)).reshape((len(x), len(y), len(z), 3))

        affine = np.array(
            [
                [resolution[0], 0, 0, limits[0][0]],
                [0, resolution[1], 0, limits[1][0]],
                [0, 0, resolution[2], limits[2][0]],
                [0, 0, 0, 1],
            ]
        )

        nib.save(nib.Nifti1Image(data, affine), fn)

    def optimize_deformations(
        self, optimization_surface: Msh, affine: npt.NDArray[np.float_]
    ):
        # TODO Test for any casings and raise exception
        coil_deformations = self.deformations
        if len(coil_deformations) == 0:
            raise ValueError(
                "The coil has no deformations to optimize the coil element positions with."
            )

        cost_surface_tree = optimization_surface.get_AABBTree()
        deformation_ranges = np.array([deform.range for deform in coil_deformations])

        intersecting, min_found_distance = self._get_current_deformation_scores(
            cost_surface_tree, affine
        )

        if intersecting:
            raise ValueError("Initial intersection detected.")

        initial_abs_mean_dist = np.abs(
            self._get_current_deformation_scores(cost_surface_tree, affine)[1]
        )
        initial_deformation_settings = np.array(
            [coil_deformation.current for coil_deformation in coil_deformations]
        )
        best_deformation_settings = np.copy(initial_deformation_settings)

        def cost_f_x0(x, x0):
            for coil_deformation, deformation_setting in zip(coil_deformations, x0 + x):
                coil_deformation.current = deformation_setting
            intersecting, distance = self._get_current_deformation_scores(
                cost_surface_tree, affine
            )
            f = initial_abs_mean_dist * intersecting + distance
            if not intersecting:
                nonlocal min_found_distance
                if f < min_found_distance:
                    nonlocal best_deformation_settings
                    best_deformation_settings = x0 + x
                    min_found_distance = f
            return f

        cost_f = lambda x: cost_f_x0(x, initial_deformation_settings)
        min_found_distance = cost_f(np.zeros_like(initial_deformation_settings))

        opt.direct(
            cost_f,
            bounds=list(
                deformation_ranges - initial_deformation_settings[:, np.newaxis]
            ),
        )

        cost_f = lambda x: cost_f_x0(x, 0)
        opt.minimize(
            cost_f,
            x0=np.copy(best_deformation_settings),
            method="L-BFGS-B",
            options={"eps": 0.001, "maxls": 100},
            bounds=deformation_ranges,
        )

        intermediate_best_deformation_settings = np.copy(best_deformation_settings)

        # refine univariate
        for i in range(len(intermediate_best_deformation_settings)):
            cost1 = lambda xx: cost_f(
                np.concatenate(
                    (
                        intermediate_best_deformation_settings[:i],
                        [xx],
                        intermediate_best_deformation_settings[i + 1 :],
                    ),
                    axis=None,
                )
            )
            opt.minimize(
                cost1,
                x0=intermediate_best_deformation_settings[i],
                method="L-BFGS-B",
                options={"eps": 0.001, "maxls": 100},
                bounds=[deformation_ranges[i]],
            )

        for coil_deformation, deformation_setting in zip(
            coil_deformations, best_deformation_settings
        ):
            coil_deformation.current = deformation_setting
        return initial_abs_mean_dist, min_found_distance

    def _get_current_deformation_scores(
        self, cost_surface_tree, affine: npt.NDArray[np.float_]
    ) -> tuple[bool, float]:
        (
            casing_points,
            min_distance_points,
            intersect_points,
        ) = self.get_deformed_casing_coordinates(affine)

        min_distance_points = (
            min_distance_points if len(min_distance_points) > 0 else casing_points
        )
        intersect_points = (
            intersect_points if len(intersect_points) > 0 else casing_points
        )

        return cost_surface_tree.any_point_inside(intersect_points), np.mean(
            np.sqrt(cost_surface_tree.min_sqdist(min_distance_points))
        )

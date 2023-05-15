import json
from typing import Optional
import numpy as np
import numpy.typing as npt
import nibabel as nib
import re

from simnibs.mesh_tools.mesh_io import Elements, Msh, NodeData, Nodes
from simnibs.simulation.coil.coil_deformation import CoilRotation, CoilTranslation
from simnibs.simulation.coil.tms_stimulator import TMSStimulator

from .coil_model import CoilModel
from .coil_element import (
    CoilDipoles,
    CoilElement,
    CoilLineElements,
    CoilLinePoints,
    CoilSampledElements,
    DirectionalCoilElement,
)


class Coil:
    def __init__(
        self,
        name: str,
        brand: str,
        version: str,
        limits: npt.NDArray[np.float_],
        resolution: npt.NDArray[np.float_],
        coil_elements: list[CoilElement],
        coil_models: list[CoilModel],
    ):
        self.name = name
        self.brand = brand
        self.version = version
        self.limits = limits
        self.resolution = resolution
        self.coil_elements = coil_elements
        self.coil_models = coil_models

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

    @classmethod
    def from_file(
        cls,
        fn: str,
        name: str = "Unknown",
        brand: str = "Unknown",
        version: str = "Unknown",
    ):
        if fn.endswith(".tcd"):
            return Coil.from_tcd(fn)
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
        version: str = "Unknown",
    ):
        # TODO CCL file needed?

        with open(fn, "r") as f:
            header = f.readline()

        pairs = header.replace("\n", "").split(';')[1:]

        header_dict = {}
        for pair in pairs:
            key, value = pair.split('=')

            if value == 'none':
                value = None
            elif '.' in value:
                value = float(value)
            elif value.isdigit():
                value = int(value)
            elif ',' in value:
                value = np.fromstring(value, dtype=int, sep=',')

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

        # parse resolution
        res = []
        a = header_dict.get('resolution')
        if a is None:
            res.append(None)
        else:
            a = np.atleast_1d(a)
            if len(a) < 3:
                for i in range(len(a), 3):
                    a = np.concatenate((a, (a[i - 1],)))
            res = a

        ccd_file = np.loadtxt(fn, skiprows=2)

        # if there is only 1 dipole, loadtxt return as array of the wrong shape
        if len(np.shape(ccd_file)) == 1:
            a = np.zeros([1, 6])
            a[0, 0:3] = ccd_file[0:3]
            a[0, 3:] = ccd_file[3:]
            ccd_file = a

        dipole_positions = ccd_file[:, 0:3]
        dipole_moments = ccd_file[:, 3:]
        stimulator = TMSStimulator(header_dict.get("stimulator", 'Unknown'), 'Unknown', header_dict.get("dIdtmax", 0))

        coil_elements = [CoilDipoles('Unknown', [], dipole_positions, dipole_moments, stimulator)]

        return cls(
            header_dict.get("coilname", 'Unknown'),
            header_dict.get("brand", 'Unknown'),
            version,
            np.array(bb),
            np.array(res),
            coil_elements,
            [],
        )

    @classmethod
    def from_tcd(cls, fn: str):
        with open(fn, "r") as fid:
            coil = json.loads(fid.read())

        coil_models = []
        for coil_model in coil["coilModels"]:
            points = np.array(coil_model["points"])
            faces = np.array(coil_model["faces"]) + 1
            min_distance_points = np.array(coil_model["minDistancePoints"])
            intersect_points = np.array(coil_model["intersectPoints"])

            mesh = Msh(Nodes(points), Elements(triangles=faces))

            coil_models.append(CoilModel(mesh, min_distance_points, intersect_points))

        deformations = []
        for deform in coil["deformList"]:
            initial = deform.get("initial", 0.0)
            deform_range = deform.get("range")
            if deform["type"] == "x":
                deformations.append(CoilTranslation(initial, deform_range, 0))
            elif deform["type"] == "y":
                deformations.append(CoilTranslation(initial, deform_range, 1))
            elif deform["type"] == "z":
                deformations.append(CoilTranslation(initial, deform_range, 2))
            elif deform["type"] == "rot2p":
                point_1 = np.array(deform["point1"])
                point_2 = np.array(deform["point2"])
                deformations.append(
                    CoilRotation(initial, deform_range, point_1, point_2)
                )
            else:
                raise ValueError(f"Invalid deformation type: {deform['type']}")

        stimulators = []
        for stimulator in coil["stimulatorList"]:
            name = stimulator["name"]
            brand = stimulator["brand"]
            max_di_dt = stimulator["maxdIdt"]
            stimulators.append(TMSStimulator(name, brand, max_di_dt))

        coil_elements = []
        for coil_element in coil["coilElementList"]:
            name = coil_element.get("name", "Unknown")
            stimulator = stimulators[coil_element["stimulator"]]
            element_deformations = [
                deformations[i] for i in coil_element["deformations"]
            ]
            points = np.array(coil_element["points"])
            weights = coil_element.get("weights", None)

            if coil_element["type"] == 1:
                values = np.array(coil_element["values"])
                coil_elements.append(
                    CoilDipoles(
                        name, element_deformations, points, values, stimulator, weights
                    )
                )
            elif coil_element["type"] == 2:
                values = np.array(coil_element["values"])
                coil_elements.append(
                    CoilLineElements(
                        name, element_deformations, points, values, stimulator, weights
                    )
                )
            elif coil_element["type"] == 3:
                coil_elements.append(
                    CoilLinePoints(
                        name, element_deformations, points, stimulator, weights
                    )
                )
            elif coil_element["type"] == 4:
                values = np.array(coil_element["values"])
                coil_elements.append(
                    CoilSampledElements(
                        name, element_deformations, points, values, stimulator, weights
                    )
                )

        return cls(
            coil["name"],
            coil["brand"],
            coil["version"],
            coil["limits"],
            coil["resolution"],
            coil_elements,
            coil_models,
        )

    def write_tcd(self, fn: str):
        json_coil_models = []
        for coil_model in self.coil_models:
            json_coil_model = {}
            json_coil_model["points"] = coil_model.mesh.nodes.node_coord.tolist()
            json_coil_model["faces"] = (
                coil_model.mesh.elm.node_number_list - 1
            ).tolist()
            json_coil_model[
                "minDistancePoints"
            ] = coil_model.min_distance_points.tolist()
            json_coil_model["intersectPoints"] = coil_model.intersect_points.tolist()
            json_coil_models.append(json_coil_model)

        json_deforms = []
        deformations = []

        json_stimulators = []
        stimulators = []

        json_coil_elements = []
        for coil_element in self.coil_elements:
            for deformation in coil_element.deformations:
                if deformation in deformations:
                    continue

                deformations.append(deformation)
                json_deformation = {}
                json_deformation["initial"] = deformation.initial
                json_deformation["range"] = list(deformation.range)

                if isinstance(deformation, CoilTranslation):
                    if deformation.axis == 0:
                        json_deformation["type"] = "x"
                    elif deformation.axis == 1:
                        json_deformation["type"] = "y"
                    elif deformation.axis == 2:
                        json_deformation["type"] = "z"
                    else:
                        raise ValueError(
                            f"Translation axis ({deformation.axis}) out of range (0-2)"
                        )
                elif isinstance(deformation, CoilRotation):
                    json_deformation["type"] = "rot2p"
                    json_deformation["point1"] = deformation.point_1.tolist()
                    json_deformation["point2"] = deformation.point_2.tolist()

                json_deforms.append(json_deformation)

            if coil_element.stimulator not in stimulators:
                stimulators.append(coil_element.stimulator)
                json_stimulator = {}
                json_stimulator["name"] = coil_element.stimulator.name
                json_stimulator["brand"] = coil_element.stimulator.brand
                json_stimulator["maxdIdt"] = coil_element.stimulator.max_di_dt
                json_stimulators.append(json_stimulator)

            json_coil_element = {}
            json_coil_element["name"] = coil_element.name
            json_coil_element["stimulator"] = stimulators.index(coil_element.stimulator)
            if isinstance(coil_element, CoilDipoles):
                json_coil_element["type"] = 1
            elif isinstance(coil_element, CoilLineElements):
                json_coil_element["type"] = 2
            elif isinstance(coil_element, CoilLinePoints):
                json_coil_element["type"] = 3
            elif isinstance(coil_element, CoilSampledElements):
                json_coil_element["type"] = 4
            json_coil_element["deformations"] = [
                deformations.index(x) for x in coil_element.deformations
            ]
            json_coil_element["points"] = coil_element.points.tolist()
            if isinstance(coil_element, DirectionalCoilElement):
                json_coil_element["values"] = coil_element.values.tolist()
            if coil_element.weights is not None:
                json_coil_element["weights"] = coil_element.weights.tolist()

            json_coil_elements.append(json_coil_element)

        coil = {}
        coil["name"] = self.name
        coil["brand"] = self.brand
        coil["version"] = self.version
        coil["limits"] = self.limits.tolist()
        coil["resolution"] = self.resolution.tolist()
        coil["coilElementList"] = json_coil_elements
        coil["stimulatorList"] = json_stimulators
        coil["deformList"] = json_deforms
        coil["coilModels"] = json_coil_models

        with open(fn, "w") as json_file:
            json.dump(coil, json_file)

    @classmethod
    def from_nifti(
        cls,
        fn: str,
        name: str = "Unknown",
        brand: str = "Unknown",
        version: str = "Unknown",
    ):
        nifti = nib.load(fn)
        data = nifti.get_fdata()
        affine = nifti.affine

        coords = np.array(
            np.meshgrid(
                np.arange(data.shape[0]),
                np.arange(data.shape[1]),
                np.arange(data.shape[2]),
                indexing="ij",
            )
        )
        coords = coords.reshape((3, -1)).T
        coords = coords @ affine[:3, :3].T + affine[None, :3, 3]
        data = data.reshape((3, -1)).T

        limits = np.array(
            [
                [np.min(coords[:, 0]), np.max(coords[:, 0])],
                [np.min(coords[:, 1]), np.max(coords[:, 1])],
                [np.min(coords[:, 2]), np.max(coords[:, 2])],
            ]
        )
        resolution = np.array(
            [
                affine[0][0],
                affine[1][1],
                affine[2][2],
            ]
        )

        print(limits)
        print(resolution)

        # TODO Figure out the Stimulator
        coil_elements = [CoilSampledElements(name, [], coords, data, None)]

        return cls(
            name,
            brand,
            version,
            limits,
            resolution,
            coil_elements,
            [],
        )

    def write_nifti(
        self,
        fn: str,
        limits: Optional[npt.NDArray[np.float_]] = None,
        resolution: Optional[npt.NDArray[np.float_]] = None,
    ):
        limits = limits or self.limits
        resolution = resolution or self.resolution

        dims = [int((max_ - min_) // res) for [min_, max_], res in zip(limits, resolution)]

        dx = np.spacing(1e4)
        x = np.linspace(limits[0][0], limits[0][1] + dx, dims[0])
        y = np.linspace(limits[1][0], limits[1][1] + dx, dims[1])
        z = np.linspace(limits[2][0], limits[2][1] + dx, dims[2])
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

        print(data)

        img = nib.Nifti1Image(data, affine)

        # Save the image
        nib.save(img, fn)

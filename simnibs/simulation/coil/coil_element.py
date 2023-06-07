from abc import ABC, abstractmethod
from copy import deepcopy
from typing import Optional

import fmm3dpy
import numpy as np
import numpy.typing as npt
from scipy import ndimage

from simnibs.mesh_tools.mesh_io import Elements, Msh, Nodes
from simnibs.simulation.coil.coil_constants import CoilElementTag
from simnibs.simulation.coil.coil_model import CoilModel
from simnibs.simulation.coil.tcd_element import TcdElement

from .coil_deformation import CoilDeformation
from .tms_stimulator import TMSStimulator


class CoilElement(ABC, TcdElement):
    def __init__(
        self,
        name: Optional[str],
        element_casing: Optional[CoilModel],
        deformations: Optional[list[CoilDeformation]],
        stimulator: Optional[TMSStimulator],
        weights: Optional[npt.NDArray[np.float_]] = None,
    ):
        self.name = name
        self.element_casing = element_casing
        self.deformations = deformations if deformations is not None else []
        self.stimulator = stimulator
        self.weights = weights if weights is not None else np.array([])

    @abstractmethod
    def get_a_field(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        pass

    def get_da_dt(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        di_dt: float,
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        return di_dt * self.get_a_field(target_positions, coil_affine, eps)

    def get_combined_transformation(
        self, affine_matrix: Optional[npt.NDArray[np.float_]] = None
    ) -> npt.NDArray[np.float_]:
        if affine_matrix is None:
            affine_matrix = np.eye(4)
        for deformation in self.deformations:
            affine_matrix = affine_matrix @ deformation.as_matrix()

        return affine_matrix

    def get_deformed_casing_coordinates(
        self, affine_matrix: Optional[npt.NDArray[np.float_]] = None
    ) -> tuple[npt.NDArray[np.float_], npt.NDArray[np.float_], npt.NDArray[np.float_]]:
        affine_matrix = self.get_combined_transformation(affine_matrix)
        transformed_coordinates = [np.array([]), np.array([]), np.array([])]
        if self.element_casing is not None:
            transformed_coordinates[0] = self.element_casing.get_points(
                affine_matrix
            )
            transformed_coordinates[
                1
            ] = self.element_casing.get_min_distance_points(affine_matrix)
            transformed_coordinates[
                2
            ] = self.element_casing.get_intersect_points(affine_matrix)

        return tuple(transformed_coordinates)

    def get_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        apply_deformation: bool = True,
        include_element_casing: bool = True,
        include_optimization_points: bool = True,
        include_coil_element: bool = True,
        element_tag: int = 0,
    ) -> Msh:
        element_mesh = Msh()
        if self.element_casing is not None and include_element_casing:
            if apply_deformation:
                element_mesh = element_mesh.join_mesh(
                    self.element_casing.get_transformed_mesh(
                        self.get_combined_transformation(affine_matrix),
                        include_optimization_points,
                        element_tag,
                    )
                )
            else:
                element_mesh = element_mesh.join_mesh(
                    self.element_casing.get_transformed_mesh(
                        affine_matrix, include_optimization_points, element_tag
                    )
                )
        return element_mesh
    
    def to_tcd(self, stimulators: list[TMSStimulator], coil_models : list[CoilModel], deformations: list[CoilDeformation]) -> dict:
        tcd_coil_element = {}
        if self.name is not None:
            tcd_coil_element["name"] = self.name

        if self.stimulator is not None:
            tcd_coil_element["stimulator"] = stimulators.index(
                self.stimulator
            )

        if self.element_casing is not None:
            tcd_coil_element["elementCasing"] = coil_models.index(
                self.element_casing
            )    
        
        if len(self.deformations) > 0:
            tcd_coil_element["deformations"] = [
                deformations.index(x) for x in self.deformations
            ]

        if len(self.weights) > 0:
            tcd_coil_element["weights"] = self.weights.tolist()

        return tcd_coil_element
    
    @classmethod
    def from_tcd(cls, tcd_coil_element: dict, stimulators: list[TMSStimulator], coil_models: list[CoilModel], deformations: list[CoilDeformation]):
        name = tcd_coil_element.get("name")
        stimulator = (
            None
            if tcd_coil_element.get("stimulator") is None
            else stimulators[tcd_coil_element["stimulator"]]
        )
        element_casing = (
            None
            if tcd_coil_element.get("elementCasing") is None
            else coil_models[tcd_coil_element["elementCasing"]]
        )
        element_deformations = (
            None
            if tcd_coil_element.get("deformations") is None
            else [deformations[i] for i in tcd_coil_element["deformations"]]
        )
        points = np.array(tcd_coil_element["points"])
        weights = (
            None
            if tcd_coil_element.get("weights") is None
            else np.array(tcd_coil_element["weights"])
        )

        if tcd_coil_element["type"] == 1:
            values = np.array(tcd_coil_element["values"])
            return CoilDipoles(
                name,
                element_casing,
                element_deformations,
                points,
                values,
                stimulator,
                weights,
            )
            
        elif tcd_coil_element["type"] == 2:
            values = np.array(tcd_coil_element["values"])
            return CoilLineElements(
                name,
                element_casing,
                element_deformations,
                points,
                values,
                stimulator,
                weights,
            )
        elif tcd_coil_element["type"] == 3:
            return CoilLinePoints(
                name,
                element_casing,
                element_deformations,
                points,
                stimulator,
                weights,
            )
        elif tcd_coil_element["type"] == 4:
            data = np.array(tcd_coil_element["data"])
            affine = np.array(tcd_coil_element["affine"])
            return CoilSampledGridElements(
                    name,
                    element_casing,
                    element_deformations,
                    data,
                    affine,
                    stimulator,
                    weights,
                )
        else:
             raise ValueError(f"Invalid coil element type: {tcd_coil_element['type']}")


class PositionalCoilElement(CoilElement, ABC):
    def __init__(
        self,
        name: Optional[str],
        element_casing: Optional[CoilModel],
        deformations: Optional[list[CoilDeformation]],
        points: npt.NDArray[np.float_],
        stimulator: Optional[TMSStimulator],
        weights: Optional[npt.NDArray[np.float_]] = None,
    ):
        super().__init__(name, element_casing, deformations, stimulator, weights)
        self.points = points

    @abstractmethod
    def get_a_field(self, msh: Msh, coil_affine: npt.NDArray[np.float_]):
        pass

    def get_transformed_points(
        self,
        affine_matrix: Optional[npt.NDArray[np.float_]] = None,
        apply_deformation: bool = True,
    ):
        if affine_matrix is None:
            affine_matrix = np.eye(4)
        if apply_deformation:
            affine_matrix = self.get_combined_transformation(affine_matrix)
        return self.points @ affine_matrix[:3, :3].T + affine_matrix[None, :3, 3]


class DirectionalCoilElement(PositionalCoilElement, ABC):
    def __init__(
        self,
        name: Optional[str],
        element_casing: Optional[CoilModel],
        deformations: Optional[list[CoilDeformation]],
        points: npt.NDArray[np.float_],
        values: npt.NDArray[np.float_],
        stimulator: Optional[TMSStimulator],
        weights: Optional[npt.NDArray[np.float_]] = None,
    ):
        super().__init__(
            name, element_casing, deformations, points, stimulator, weights
        )
        self.values = values

    @abstractmethod
    def get_a_field(self, msh: Msh, coil_affine: npt.NDArray[np.float_]):
        pass

    def get_transformed_values(
        self,
        affine_matrix: Optional[npt.NDArray[np.float_]] = None,
        apply_deformation: bool = True,
    ):
        if affine_matrix is None:
            affine_matrix = np.eye(4)
        if apply_deformation:
            affine_matrix = self.get_combined_transformation(affine_matrix)
        return self.values @ affine_matrix[:3, :3].T


class CoilDipoles(DirectionalCoilElement):
    def get_a_field(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        dipole_moment = self.get_transformed_values(coil_affine)
        dipole_position_m = self.get_transformed_points(coil_affine) * 1e-3
        target_positions_m = target_positions * 1e-3
        if dipole_moment.shape[0] < 300:
            out = fmm3dpy.l3ddir(
                charges=dipole_moment.T,
                sources=dipole_position_m.T,
                targets=target_positions_m.T,
                nd=3,
                pgt=2,
            )
        else:
            out = fmm3dpy.lfmm3d(
                charges=dipole_moment.T,
                sources=dipole_position_m.T,
                targets=target_positions_m.T,
                eps=eps,
                nd=3,
                pgt=2,
            )

        A = np.empty((target_positions_m.shape[0], 3), dtype=float)

        A[:, 0] = out.gradtarg[1][2] - out.gradtarg[2][1]
        A[:, 1] = out.gradtarg[2][0] - out.gradtarg[0][2]
        A[:, 2] = out.gradtarg[0][1] - out.gradtarg[1][0]

        A *= -1e-7

        return A

    def get_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        apply_deformation: bool = True,
        include_element_casing: bool = True,
        include_optimization_points: bool = True,
        include_coil_element: bool = True,
        element_tag: int = 0,
    ) -> Msh:
        element_mesh = super().get_mesh(
            affine_matrix,
            apply_deformation,
            include_element_casing,
            include_optimization_points,
            include_coil_element,
            element_tag,
        )
        if not include_coil_element:
            return element_mesh

        transformed_points = self.get_transformed_points(affine_matrix)
        transformed_values = transformed_points + self.get_transformed_values(
            affine_matrix
        )
        point_mesh = Msh(
            Nodes(np.concatenate((transformed_points, transformed_values))),
            Elements(
                lines=np.column_stack(
                    (
                        np.arange(len(transformed_points)),
                        np.arange(len(transformed_points)) + len(transformed_values),
                    )
                )
                + 1
            ),
        )
        point_mesh.elm.tag1[:] = element_tag + CoilElementTag.DIPOLES
        point_mesh.elm.tag2[:] = element_tag + CoilElementTag.DIPOLES

        return element_mesh.join_mesh(point_mesh)
    
    def to_tcd(self, stimulators: list[TMSStimulator], coil_models : list[CoilModel], deformations: list[CoilDeformation]) -> dict:
        tcd_coil_element = super().to_tcd(stimulators, coil_models, deformations)
        tcd_coil_element["type"] = 1
        tcd_coil_element["points"] = self.points.tolist()
        tcd_coil_element["values"] = self.values.tolist()
        return tcd_coil_element


class CoilLineElements(DirectionalCoilElement):
    def get_a_field(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        directions_m = self.get_transformed_values(coil_affine) * 1e-3
        segment_position_m = self.get_transformed_points(coil_affine) * 1e-3
        target_positions_m = target_positions * 1e-3

        if directions_m.shape[0] >= 300:
            A = fmm3dpy.l3ddir(
                sources=segment_position_m.T,
                charges=directions_m.T,
                targets=target_positions_m.T,
                nd=3,
                pgt=1,
            )
        else:
            A = fmm3dpy.lfmm3d(
                sources=segment_position_m.T,
                charges=directions_m.T,
                targets=target_positions_m.T,
                nd=3,
                eps=eps,
                pgt=1,
            )

        A = 1e-7 * A.pottarg.T
        return A

    def get_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        apply_deformation: bool = True,
        include_element_casing: bool = True,
        include_optimization_points: bool = True,
        include_coil_element: bool = True,
        element_tag: int = 0,
    ) -> Msh:
        element_mesh = super().get_mesh(
            affine_matrix,
            apply_deformation,
            include_element_casing,
            include_optimization_points,
            include_coil_element,
            element_tag,
        )
        if not include_coil_element:
            return element_mesh

        transformed_points = self.get_transformed_points(affine_matrix)
        transformed_values = transformed_points + self.get_transformed_values(
            affine_matrix
        )
        point_mesh = Msh(
            Nodes(np.concatenate((transformed_points, transformed_values))),
            Elements(
                lines=np.column_stack(
                    (
                        np.arange(len(transformed_points)),
                        np.arange(len(transformed_points)) + len(transformed_values),
                    )
                )
                + 1
            ),
        )
        point_mesh.elm.tag1[:] = element_tag + CoilElementTag.LINE_ELEMENTS
        point_mesh.elm.tag2[:] = element_tag + CoilElementTag.LINE_ELEMENTS

        return element_mesh.join_mesh(point_mesh)
    
    def to_tcd(self, stimulators: list[TMSStimulator], coil_models : list[CoilModel], deformations: list[CoilDeformation]) -> dict:
        tcd_coil_element = super().to_tcd(stimulators, coil_models, deformations)
        tcd_coil_element["type"] = 2
        tcd_coil_element["points"] = self.points.tolist()
        tcd_coil_element["values"] = self.values.tolist()
        return tcd_coil_element


class CoilLinePoints(PositionalCoilElement):
    def get_a_field(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        segment_position_m = self.get_transformed_points(coil_affine) * 1e-3
        target_positions_m = target_positions * 1e-3

        directions_m = np.zeros(segment_position_m.shape)
        directions_m[:, :-1] = np.diff(segment_position_m, axis=1)
        directions_m[:, -1] = segment_position_m[:, 0] - segment_position_m[:, -1]

        if directions_m.shape[0] >= 300:
            A = fmm3dpy.l3ddir(
                sources=segment_position_m.T,
                charges=directions_m.T,
                targets=target_positions_m.T,
                nd=3,
                pgt=1,
            )
        else:
            A = fmm3dpy.lfmm3d(
                sources=segment_position_m.T,
                charges=directions_m.T,
                targets=target_positions_m.T,
                nd=3,
                eps=eps,
                pgt=1,
            )
        A = 1e-7 * A.pottarg.T
        return A

    def get_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        apply_deformation: bool = True,
        include_element_casing: bool = True,
        include_optimization_points: bool = True,
        include_coil_element: bool = True,
        element_tag: int = 0,
    ) -> Msh:
        element_mesh = super().get_mesh(
            affine_matrix,
            apply_deformation,
            include_element_casing,
            include_optimization_points,
            include_coil_element,
            element_tag,
        )
        if not include_coil_element:
            return element_mesh

        transformed_points = self.get_transformed_points(
            affine_matrix, apply_deformation
        )
        point_mesh = Msh(
            Nodes(transformed_points),
            Elements(points=np.arange(len(transformed_points)) + 1),
        )
        point_mesh.elm.tag1[:] = element_tag + CoilElementTag.LINE_POINTS
        point_mesh.elm.tag2[:] = element_tag + CoilElementTag.LINE_POINTS

        return element_mesh.join_mesh(point_mesh)
    
    def to_tcd(self, stimulators: list[TMSStimulator], coil_models : list[CoilModel], deformations: list[CoilDeformation]) -> dict:
        tcd_coil_element = super().to_tcd(stimulators, coil_models, deformations)
        tcd_coil_element["type"] = 3
        tcd_coil_element["points"] = self.points.tolist()
        return tcd_coil_element


class CoilSampledGridElements(CoilElement):
    def __init__(
        self,
        name: Optional[str],
        element_casing: Optional[CoilModel],
        deformations: Optional[list[CoilDeformation]],
        data: npt.NDArray[np.float_],
        affine: npt.NDArray[np.float_],
        stimulator: Optional[TMSStimulator],
        weights: Optional[npt.NDArray[np.float_]] = None,
    ):
        super().__init__(name, element_casing, deformations, stimulator, weights)
        self.data = data
        self.affine = affine

    def get_a_field(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        combined_affine = self.get_combined_transformation(coil_affine)
        iM = np.linalg.pinv(self.affine) @ np.linalg.pinv(combined_affine)

        target_voxle_coordinates = iM[:3, :3] @ target_positions.T + iM[:3, 3][:, np.newaxis]

        # Interpolates the values of the field in the given coordinates
        out = np.zeros((3, target_voxle_coordinates.shape[1]))
        for dim in range(3):
            out[dim] = ndimage.map_coordinates(
                np.asanyarray(self.data)[..., dim], target_voxle_coordinates, order=1
            )

        # Rotates the field
        return (coil_affine[:3, :3] @ out).T
    
    def get_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        apply_deformation: bool = True,
        include_element_casing: bool = True,
        include_optimization_points: bool = True,
        include_coil_element: bool = True,
        element_tag: int = 0,
    ) -> Msh:
        element_mesh = super().get_mesh(
            affine_matrix,
            apply_deformation,
            include_element_casing,
            include_optimization_points,
            include_coil_element,
            element_tag,
        )
        if not include_coil_element:
            return element_mesh
        
        combined_affine = self.affine @ self.get_combined_transformation(affine_matrix)

        voxel_coordinates = np.array(list(np.ndindex(self.data.shape[0],self.data.shape[1],self.data.shape[2])))

        points = voxel_coordinates @ combined_affine[:3, :3].T + combined_affine[None, :3, 3]
        targets = points + self.data.reshape(-1, 3)

        point_mesh = Msh(
            Nodes(np.concatenate((points, targets))),
            Elements(
                lines=np.column_stack(
                    (
                        np.arange(len(points)),
                        np.arange(len(targets)) + len(points),
                    )
                )
                + 1
            ),
        )
        point_mesh.elm.tag1[:] = element_tag + CoilElementTag.SAMPLED_GRID_ELEMENTS
        point_mesh.elm.tag2[:] = element_tag + CoilElementTag.SAMPLED_GRID_ELEMENTS

        return element_mesh.join_mesh(point_mesh)
    
    def to_tcd(self, stimulators: list[TMSStimulator], coil_models : list[CoilModel], deformations: list[CoilDeformation]) -> dict:
        tcd_coil_element = super().to_tcd(stimulators, coil_models, deformations)
        tcd_coil_element["type"] = 4
        tcd_coil_element["data"] = self.data.tolist()
        tcd_coil_element["affine"] = self.affine.tolist()
        return tcd_coil_element

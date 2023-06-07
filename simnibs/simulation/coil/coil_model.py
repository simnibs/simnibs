from copy import deepcopy
from typing import Optional

from simnibs.simulation.coil.coil_constants import CoilElementTag
from simnibs.simulation.coil.tcd_element import TcdElement

from ...mesh_tools.mesh_io import Elements, Msh, Nodes

import numpy as np
import numpy.typing as npt


class CoilModel(TcdElement):
    def __init__(
        self,
        mesh: Msh,
        min_distance_points: Optional[npt.NDArray[np.float_]],
        intersect_points: Optional[npt.NDArray[np.float_]],
    ):
        self.mesh = mesh
        self.min_distance_points = (
            np.array([]) if min_distance_points is None else min_distance_points
        )
        self.intersect_points = (
            np.array([]) if intersect_points is None else intersect_points
        )

    def get_transformed_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        include_optimization_points: bool = True,
        model_tag: int = 0,
    ) -> Msh:
        transformed_mesh = deepcopy(self.mesh)
        transformed_mesh.nodes.node_coord = self.get_points(affine_matrix)
        transformed_mesh.elm.tag1[:] = model_tag + CoilElementTag.COIL_CASING
        transformed_mesh.elm.tag2[:] = model_tag + CoilElementTag.COIL_CASING
        if not include_optimization_points:
            return transformed_mesh

        if len(self.min_distance_points) > 0:
            transformed_min_distance_points = self.get_min_distance_points(
                affine_matrix
            )
            point_mesh = Msh(
                Nodes(transformed_min_distance_points),
                Elements(points=np.arange(len(transformed_min_distance_points)) + 1),
            )
            point_mesh.elm.tag1[:] = (
                model_tag + CoilElementTag.COIL_CASING_MIN_DISTANCE_POINTS
            )
            point_mesh.elm.tag2[:] = (
                model_tag + CoilElementTag.COIL_CASING_MIN_DISTANCE_POINTS
            )
            transformed_mesh = transformed_mesh.join_mesh(point_mesh)
        if len(self.intersect_points) > 0:
            transformed_intersect_points = self.get_intersect_points(affine_matrix)
            point_mesh = Msh(
                Nodes(transformed_intersect_points),
                Elements(points=np.arange(len(transformed_intersect_points)) + 1),
            )
            point_mesh.elm.tag1[:] = (
                model_tag + CoilElementTag.COIL_CASING_INTERSECT_POINTS
            )
            point_mesh.elm.tag2[:] = (
                model_tag + CoilElementTag.COIL_CASING_INTERSECT_POINTS
            )
            transformed_mesh = transformed_mesh.join_mesh(point_mesh)
        return transformed_mesh

    def get_points(
        self, affine_matrix: Optional[npt.NDArray[np.float_]] = None
    ) -> npt.NDArray[np.float_]:
        if affine_matrix is None:
            affine_matrix = np.eye(4)
        if len(self.mesh.nodes.node_coord) == 0:
            return np.array([])
        return (
            self.mesh.nodes.node_coord @ affine_matrix[:3, :3].T
            + affine_matrix[None, :3, 3]
        )

    def get_min_distance_points(
        self, affine_matrix: Optional[npt.NDArray[np.float_]] = None
    ) -> npt.NDArray[np.float_]:
        if affine_matrix is None:
            affine_matrix = np.eye(4)
        if len(self.min_distance_points) == 0:
            return np.array([])
        return (
            self.min_distance_points @ affine_matrix[:3, :3].T
            + affine_matrix[None, :3, 3]
        )

    def get_intersect_points(
        self, affine_matrix: Optional[npt.NDArray[np.float_]] = None
    ) -> npt.NDArray[np.float_]:
        if affine_matrix is None:
            affine_matrix = np.eye(4)
        if len(self.intersect_points) == 0:
            return np.array([])
        return (
            self.intersect_points @ affine_matrix[:3, :3].T + affine_matrix[None, :3, 3]
        )

    def to_tcd(self) -> dict:
        tcd_coil_model = {}
        tcd_coil_model["points"] = self.mesh.nodes.node_coord.tolist()
        tcd_coil_model["faces"] = (self.mesh.elm.node_number_list[:, :3] - 1).tolist()
        if len(self.min_distance_points) > 0:
            tcd_coil_model["minDistancePoints"] = self.min_distance_points.tolist()
        if len(self.intersect_points) > 0:
            tcd_coil_model["intersectPoints"] = self.intersect_points.tolist()

        return tcd_coil_model

    @classmethod
    def from_tcd(cls, tcd_coil_model: dict):
        points = np.array(tcd_coil_model["points"])
        faces = np.array(tcd_coil_model["faces"]) + 1
        min_distance_points = (
            None
            if tcd_coil_model.get("minDistancePoints") is None
            else np.array(tcd_coil_model["minDistancePoints"])
        )
        intersect_points = (
            None
            if tcd_coil_model.get("intersectPoints") is None
            else np.array(tcd_coil_model["intersectPoints"])
        )

        mesh = Msh(Nodes(points), Elements(triangles=faces))

        return cls(mesh, min_distance_points, intersect_points)

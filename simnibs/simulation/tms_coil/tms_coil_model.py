from copy import deepcopy
from typing import Optional

import numpy as np
import numpy.typing as npt

from simnibs.simulation.tms_coil.tcd_element import TcdElement
from simnibs.simulation.tms_coil.tms_coil_constants import TmsCoilElementTag

from ...mesh_tools.mesh_io import Elements, Msh, Nodes


class TmsCoilModel(TcdElement):
    """A representation of a coil model including the casing and optimization points

    Parameters
    ----------
    mesh : Msh
        The coil model
    min_distance_points : Optional[npt.NDArray[np.float_]]
        Min distance points used for optimization of coil deformation
    intersect_points : Optional[npt.NDArray[np.float_]]
        Intersection points used for optimization of coil deformation

    Attributes
    ----------------------
    mesh : Msh
        The coil model
    min_distance_points : Optional[npt.NDArray[np.float_]]
        Min distance points used for optimization of coil deformation
    intersect_points : Optional[npt.NDArray[np.float_]]
        Intersection points used for optimization of coil deformation

    """
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

    def get_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        include_casing: bool = True,
        include_optimization_points: bool = True,
        model_index: int = 0,
    ) -> Msh:
        """Returns the casing as a mesh, optionally including the min distance points and the intersection points

        Parameters
        ----------
        affine_matrix : npt.NDArray[np.float_]
            The affine transformation that is applied to the coil model
        include_casing : bool, optional
            Whether or not to include the casing mesh, by default True
        include_optimization_points : bool, optional
            Whether or not to include the min distance and intersection points, by default True
        model_index : int, optional
            The index of this coil model, by default 0

        Returns
        -------
        Msh
            The coil casing as a mesh
        """
        mesh = Msh()
        model_base_tag = model_index * TmsCoilElementTag.INDEX_OFFSET

        if include_casing:
            transformed_mesh = deepcopy(self.mesh)
            transformed_mesh.nodes.node_coord = self.get_points(affine_matrix)
            transformed_mesh.elm.tag1[:] = model_base_tag + TmsCoilElementTag.COIL_CASING
            transformed_mesh.elm.tag2[:] = model_base_tag + TmsCoilElementTag.COIL_CASING
            mesh = mesh.join_mesh(transformed_mesh)

        if not include_optimization_points:
            return mesh

        if len(self.min_distance_points) > 0:
            transformed_min_distance_points = self.get_min_distance_points(
                affine_matrix
            )
            point_mesh = Msh(
                Nodes(transformed_min_distance_points),
                Elements(points=np.arange(len(transformed_min_distance_points)) + 1),
            )
            point_mesh.elm.tag1[:] = (
                model_base_tag + TmsCoilElementTag.COIL_CASING_MIN_DISTANCE_POINTS
            )
            point_mesh.elm.tag2[:] = (
                model_base_tag + TmsCoilElementTag.COIL_CASING_MIN_DISTANCE_POINTS
            )
            mesh = mesh.join_mesh(point_mesh)
        if len(self.intersect_points) > 0:
            transformed_intersect_points = self.get_intersect_points(affine_matrix)
            point_mesh = Msh(
                Nodes(transformed_intersect_points),
                Elements(points=np.arange(len(transformed_intersect_points)) + 1),
            )
            point_mesh.elm.tag1[:] = (
                model_base_tag + TmsCoilElementTag.COIL_CASING_INTERSECT_POINTS
            )
            point_mesh.elm.tag2[:] = (
                model_base_tag + TmsCoilElementTag.COIL_CASING_INTERSECT_POINTS
            )
            mesh = mesh.join_mesh(point_mesh)
        return mesh

    def get_points(
        self, affine_matrix: Optional[npt.NDArray[np.float_]] = None
    ) -> npt.NDArray[np.float_]:
        """Returns the coil model points, optionally transformed by an affine transformation

        Parameters
        ----------
        affine_matrix : Optional[npt.NDArray[np.float_]], optional
            The affine transformation that is applied to the coil model, by default None

        Returns
        -------
        npt.NDArray[np.float_]
            The coil model points
        """
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
        """Returns the coil model min distance points, optionally transformed by an affine transformation

        Parameters
        ----------
        affine_matrix : Optional[npt.NDArray[np.float_]], optional
            The affine transformation that is applied to the coil model, by default None

        Returns
        -------
        npt.NDArray[np.float_]
            The coil model min distance points
        """
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
        """Returns the coil model intersection points, optionally transformed by an affine transformation

        Parameters
        ----------
        affine_matrix : Optional[npt.NDArray[np.float_]], optional
            The affine transformation that is applied to the coil model, by default None

        Returns
        -------
        npt.NDArray[np.float_]
            The coil model coil model intersection points
        """
        if affine_matrix is None:
            affine_matrix = np.eye(4)
        if len(self.intersect_points) == 0:
            return np.array([])
        return (
            self.intersect_points @ affine_matrix[:3, :3].T + affine_matrix[None, :3, 3]
        )
    
    def merge(self, other_coil_model: "TmsCoilModel") -> "TmsCoilModel":
        """Merges two coil models into one and returning the combination as a new coil model

        Parameters
        ----------
        other_coil_model : TmsCoilModel
            The coil model that is supposed to be merged with this coil model

        Returns
        -------
        TmsCoilModel
            The combined coil model
        """
        merged_coil_model = deepcopy(self)
        merged_coil_model.mesh = merged_coil_model.mesh.join_mesh(other_coil_model.mesh)
        merged_coil_model.min_distance_points = np.concatenate((merged_coil_model.min_distance_points, other_coil_model.min_distance_points), axis=0)
        merged_coil_model.intersect_points = np.concatenate((merged_coil_model.intersect_points, other_coil_model.intersect_points), axis=0)

        return merged_coil_model


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
    def from_tcd_dict(cls, tcd_coil_model: dict):
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

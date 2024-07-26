import base64
from copy import deepcopy
from typing import Optional
import zlib

import numpy as np
import numpy.typing as npt
from scipy import ndimage
from scipy.spatial import KDTree

from simnibs.segmentation.marching_cube import marching_cubes_lewiner

from simnibs.simulation.tms_coil.tcd_element import TcdElement
from simnibs.simulation.tms_coil.tms_coil_constants import TmsCoilElementTag

from ...mesh_tools.mesh_io import Elements, Msh, Nodes


class TmsCoilModel(TcdElement):
    """A representation of a coil model including the casing and optimization points

    Parameters
    ----------
    mesh : Msh
        The coil model
    min_distance_points : Optional[npt.ArrayLike], optional
        Min distance points used for optimization of coil deformation, by default None

    Attributes
    ----------------------
    mesh : Msh
        The coil model
    min_distance_points : npt.NDArray[np.float64]
        Min distance points used for optimization of coil deformation

    """

    def __init__(
        self,
        mesh: Msh,
        min_distance_points: Optional[npt.ArrayLike] = None,
    ):
        self.mesh = mesh
        self.mesh.nodes.node_coord = np.array(
            self.mesh.nodes.node_coord, dtype=np.float64
        )
        self.mesh.elm.node_number_list = np.array(
            self.mesh.elm.node_number_list, dtype=np.int64
        )

        self.min_distance_points = (
            np.array([], dtype=np.float64)
            if min_distance_points is None
            else np.array(min_distance_points, dtype=np.float64)
        )

        if min_distance_points is not None and (
            self.min_distance_points.ndim != 2 or self.min_distance_points.shape[1] != 3
        ):
            raise ValueError(
                f"Expected 'min_distance_points' to have the shape (N, 3) but shape was {self.min_distance_points.shape}"
            )

    def get_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        include_casing: bool = True,
        include_optimization_points: bool = True,
        model_index: int = 0,
    ) -> Msh:
        """Returns the casing as a mesh, optionally including the min distance points

        Parameters
        ----------
        affine_matrix : npt.NDArray[np.float_]
            The affine transformation that is applied to the coil model
        include_casing : bool, optional
            Whether or not to include the casing mesh, by default True
        include_optimization_points : bool, optional
            Whether or not to include the min distance, by default True
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
            transformed_mesh.elm.tag1[:] = (
                model_base_tag + TmsCoilElementTag.COIL_CASING
            )
            transformed_mesh.elm.tag2[:] = (
                model_base_tag + TmsCoilElementTag.COIL_CASING
            )
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
        return mesh

    def get_points(
        self, affine_matrix: Optional[npt.NDArray[np.float_]] = None
    ) -> npt.NDArray[np.float64]:
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
    ) -> npt.NDArray[np.float64]:
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
        merged_coil_model.min_distance_points = np.concatenate(
            (
                merged_coil_model.min_distance_points,
                other_coil_model.min_distance_points,
            ),
            axis=0,
        )

        return merged_coil_model

    def apply_deformations(
        self, affine_matrix: npt.NDArray[np.float_]
    ) -> "TmsCoilModel":
        deformed_min_distance_points = self.get_min_distance_points(affine_matrix)
        return TmsCoilModel(
            Msh(Nodes(self.get_points(affine_matrix)), self.mesh.elm),
            deformed_min_distance_points
            if len(deformed_min_distance_points) > 0
            else None,
        )

    def to_tcd(self, ascii_mode: bool = False) -> dict:
        tcd_coil_model = {}

        if ascii_mode:
            tcd_coil_model["points"] = self.mesh.nodes.node_coord.tolist()
            tcd_coil_model["faces"] = (
                self.mesh.elm.node_number_list[:, :3] - 1
            ).tolist()
            if len(self.min_distance_points) > 0:
                tcd_coil_model["minDistancePoints"] = self.min_distance_points.tolist()
        else:
            tcd_coil_model["points"] = base64.b64encode(
                zlib.compress(self.mesh.nodes.node_coord.tobytes())
            ).decode("ascii")
            tcd_coil_model["faces"] = base64.b64encode(
                zlib.compress((self.mesh.elm.node_number_list[:, :3] - 1).tobytes())
            ).decode("ascii")
            if len(self.min_distance_points) > 0:
                tcd_coil_model["minDistancePoints"] = base64.b64encode(
                    zlib.compress(self.min_distance_points.tobytes())
                ).decode("ascii")

        return tcd_coil_model

    @classmethod
    def from_points(cls, points: npt.ArrayLike, distance: float, grid_spacing: float):
        """Creates a coil model by generating a casing based on the distance from the points.

        Parameters
        ----------
        points : npt.ArrayLike (N, 3)
            The points the casing is generated from
        distance : float
            The distance of the casings to the points
        grid_spacing : float
            The spacing of the discretization grid used to generate the casing.
            Should be smaller than distance. Determines the resolution of the casing.
        """
        points = np.array(points, dtype=np.float64)
        limits = [
            [
                np.min(points[:, 0]) - 3 * distance - grid_spacing,
                np.max(points[:, 0]) + 3 * distance + grid_spacing,
            ],
            [
                np.min(points[:, 1]) - 3 * distance - grid_spacing,
                np.max(points[:, 1]) + 3 * distance + grid_spacing,
            ],
            [
                np.min(points[:, 2]) - 3 * distance - grid_spacing,
                np.max(points[:, 2]) + 3 * distance + grid_spacing,
            ],
        ]

        grid_x = np.arange(limits[0][0], limits[0][1], grid_spacing)
        grid_y = np.arange(limits[1][0], limits[1][1], grid_spacing)
        grid_z = np.arange(limits[2][0], limits[2][1], grid_spacing)
        grid_points = np.stack(
            np.meshgrid(grid_x, grid_y, grid_z, indexing="ij"), axis=-1
        ).reshape(-1, 3)

        tree = KDTree(points)
        distances, _ = tree.query(grid_points, k=1)
        volume_data = distances.reshape((len(grid_x), len(grid_y), len(grid_z)))

        volume_data = ndimage.gaussian_filter(volume_data, sigma=1)

        vertices, faces, _, _ = marching_cubes_lewiner(
            volume_data,
            level=distance,
            spacing=tuple(
                np.array([grid_spacing, grid_spacing, grid_spacing], dtype="float32")
            ),
            step_size=1,
            allow_degenerate=False,
        )
        vertices += [limits[0][0], limits[1][0], limits[2][0]]

        casing = Msh(Nodes(vertices), Elements(faces + 1))
        casing.smooth_surfaces(40)

        return cls(casing)

    @classmethod
    def from_tcd_dict(cls, tcd_coil_model: dict):
        if isinstance(tcd_coil_model["points"], str):
            points = np.frombuffer(
                zlib.decompress(base64.b64decode(tcd_coil_model["points"])), dtype=np.float64
            ).reshape(-1, 3)
        else:
            points = np.array(tcd_coil_model["points"])

        if isinstance(tcd_coil_model["faces"], str):
            faces = (
                np.frombuffer(
                    zlib.decompress(base64.b64decode(tcd_coil_model["faces"])), dtype=np.int64
                ).reshape(-1, 3)
                + 1
            )
        else:
            faces = np.array(tcd_coil_model["faces"]) + 1

        if tcd_coil_model.get("minDistancePoints") is None:
            min_distance_points = None
        elif isinstance(tcd_coil_model["minDistancePoints"], str):
            min_distance_points = np.frombuffer(
                zlib.decompress(base64.b64decode(tcd_coil_model["minDistancePoints"])), dtype=np.float64
            ).reshape(-1, 3)
        else:
            min_distance_points = np.array(tcd_coil_model["minDistancePoints"])

        mesh = Msh(Nodes(points), Elements(triangles=faces))

        return cls(mesh, min_distance_points)

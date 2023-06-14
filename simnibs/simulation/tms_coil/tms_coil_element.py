from abc import ABC, abstractmethod
from typing import Optional

import fmm3dpy
import numpy as np
import numpy.typing as npt
from scipy import ndimage

from simnibs.mesh_tools.mesh_io import Elements, Msh, Nodes
from simnibs.simulation.tms_coil.tcd_element import TcdElement
from simnibs.simulation.tms_coil.tms_coil_constants import TmsCoilElementTag
from simnibs.simulation.tms_coil.tms_coil_model import TmsCoilModel

from .tms_coil_deformation import TmsCoilDeformation
from .tms_stimulator import TmsStimulator


class TmsCoilElements(ABC, TcdElement):
    """A representation of a stimulating element of a TMS coil

    Parameters
    ----------
    name : Optional[str]
        The name of the element
    casing : Optional[TmsCoilModel]
        The casing of the element
    deformations : Optional[list[TmsCoilDeformation]]
        A list of all deformations of the element
    stimulator : Optional[TmsStimulator]
        The stimulator used for this element
    weights : Optional[npt.NDArray[np.float_]], optional
        An additional weighting array with one scalar per element, by default None

    Attributes
    ----------------------
    name : Optional[str]
        The name of the element
    casing : Optional[TmsCoilModel]
        The casing of the element
    deformations : Optional[list[TmsCoilDeformation]]
        A list of all deformations of the element
    stimulator : Optional[TmsStimulator]
        The stimulator used for this element
    weights : Optional[npt.NDArray[np.float_]], optional
        An additional weighting array with one scalar per element, by default None
    """

    def __init__(
        self,
        name: Optional[str],
        casing: Optional[TmsCoilModel],
        deformations: Optional[list[TmsCoilDeformation]],
        stimulator: Optional[TmsStimulator],
        weights: Optional[npt.NDArray[np.float_]] = None,
    ):
        self.name = name
        self.casing = casing
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
        """Calculates the A field applied by the coil element at each target position.

        Parameters
        ----------
        target_positions : npt.NDArray[np.float_] (N x 3)
            The points at which the A field should be calculated (in mm)
        coil_affine : npt.NDArray[np.float_] (4 x 4)
            The affine transformation that is applied to the coil
        eps : float, optional
            The requested precision, by default 1e-3

        Returns
        -------
        npt.NDArray[np.float_] (N x 3)
            The A field at every target positions in Tesla*meter
        """
        pass

    def get_da_dt(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        di_dt: float,
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        """Calculate the dA/dt field applied by the coil element at each target point

        Parameters
        ----------
        target_positions : npt.NDArray[np.float_]
            The target positions at which the dA/dt field should be calculated
        coil_affine : npt.NDArray[np.float_]
            The affine transformation that is applied to the coil
        di_dt : float
            dI/dt in A/s
        eps : float, optional
            The requested precision, by default 1e-3

        Returns
        -------
        npt.NDArray[np.float_]
            The dA/dt field at every target position
        """
        return di_dt * self.get_a_field(target_positions, coil_affine, eps)

    def get_combined_transformation(
        self, affine_matrix: Optional[npt.NDArray[np.float_]] = None
    ) -> npt.NDArray[np.float_]:
        """Returns the affine matrix that combines the deformations and the input affine matrix into one

        Parameters
        ----------
        affine_matrix : Optional[npt.NDArray[np.float_]], optional
            The affine transformation that is applied to the coil element, by default None

        Returns
        -------
        npt.NDArray[np.float_]
            The affine matrix that combines the deformations and the input affine matrix into one
        """
        if affine_matrix is None:
            affine_matrix = np.eye(4)
        for deformation in self.deformations:
            affine_matrix = affine_matrix @ deformation.as_matrix()

        return affine_matrix

    def get_casing_coordinates(
        self,
        affine_matrix: Optional[npt.NDArray[np.float_]] = None,
        apply_deformation: bool = True,
    ) -> tuple[npt.NDArray[np.float_], npt.NDArray[np.float_], npt.NDArray[np.float_]]:
        """Returns the casing points, min distance points and intersection points,
        optionally transformed by the affine matrix and deformed by the element deformation

        Parameters
        ----------
        affine_matrix : Optional[npt.NDArray[np.float_]], optional
            The affine transformation that is applied to the coil element, by default None
        apply_deformation : bool, optional
            Whether or not to apply the current deformations, by default True

        Returns
        -------
        tuple[npt.NDArray[np.float_], npt.NDArray[np.float_], npt.NDArray[np.float_]]
            The casing points, min distance points and intersection points
        """
        if apply_deformation:
            affine_matrix = self.get_combined_transformation(affine_matrix)
        else:
            affine_matrix = affine_matrix if affine_matrix is not None else np.eye(4)

        transformed_coordinates = [np.array([]), np.array([]), np.array([])]
        if self.casing is not None:
            transformed_coordinates[0] = self.casing.get_points(affine_matrix)
            transformed_coordinates[1] = self.casing.get_min_distance_points(
                affine_matrix
            )
            transformed_coordinates[2] = self.casing.get_intersect_points(affine_matrix)

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
        """Generates a mesh of the coil element, optionally transformed by the affine matrix, deformed by the element deformation,
        including the element casing, including the min distance and intersection points and including the coil element

        Parameters
        ----------
        affine_matrix : npt.NDArray[np.float_]
            The affine transformation that is applied to the coil element
        apply_deformation : bool, optional
            Whether or not to apply the current coil element deformations, by default True
        include_element_casing : bool, optional
            Whether or not to include the casing mesh, by default True
        include_optimization_points : bool, optional
            Whether or not to include the min distance and intersection points, by default True
        include_coil_element : bool, optional
            Whether or not to include the stimulating elements in the mesh, by default True
        element_tag : int, optional
            The base value used as a tag for this coil element, by default 0

        Returns
        -------
        Msh
            The generated mesh
        """
        element_mesh = Msh()
        if self.casing is not None:
            if apply_deformation:
                element_mesh = element_mesh.join_mesh(
                    self.casing.get_mesh(
                        self.get_combined_transformation(affine_matrix),
                        include_element_casing,
                        include_optimization_points,
                        element_tag,
                    )
                )
            else:
                element_mesh = element_mesh.join_mesh(
                    self.casing.get_mesh(
                        affine_matrix,
                        include_element_casing,
                        include_optimization_points,
                        element_tag,
                    )
                )
        if include_coil_element:
            element_mesh = element_mesh.join_mesh(
                self.generate_element_mesh(
                    affine_matrix, apply_deformation, element_tag
                )
            )
        return element_mesh

    @abstractmethod
    def generate_element_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        apply_deformation: bool = True,
        element_tag: int = 0,
    ) -> Msh:
        """Generate a visualization of the coil element as a mesh

        Parameters
        ----------
        affine_matrix : npt.NDArray[np.float_]
            The affine transformation that is applied to the coil element
        apply_deformation : bool, optional
            Whether or not to apply the current coil element deformations, by default True
        element_tag : int, optional
            The base value used as a tag for this coil element, by default 0

        Returns
        -------
        Msh
            The generated mesh representing the coil element
        """
        pass

    def to_tcd(
        self,
        stimulators: list[TmsStimulator],
        coil_models: list[TmsCoilModel],
        deformations: list[TmsCoilDeformation],
    ) -> dict:
        tcd_coil_element = {}
        if self.name is not None:
            tcd_coil_element["name"] = self.name

        if self.stimulator is not None:
            tcd_coil_element["stimulator"] = stimulators.index(self.stimulator)

        if self.casing is not None:
            tcd_coil_element["elementCasing"] = coil_models.index(self.casing)

        if len(self.deformations) > 0:
            tcd_coil_element["deformations"] = [
                deformations.index(x) for x in self.deformations
            ]

        if len(self.weights) > 0:
            tcd_coil_element["weights"] = self.weights.tolist()

        return tcd_coil_element

    @classmethod
    def from_tcd_dict(
        cls,
        tcd_coil_element: dict,
        stimulators: list[TmsStimulator],
        coil_models: list[TmsCoilModel],
        deformations: list[TmsCoilDeformation],
    ):
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
        weights = (
            None
            if tcd_coil_element.get("weights") is None
            else np.array(tcd_coil_element["weights"])
        )

        if tcd_coil_element["type"] == 1:
            points = np.array(tcd_coil_element["points"])
            values = np.array(tcd_coil_element["values"])
            return DipoleElements(
                name,
                element_casing,
                element_deformations,
                points,
                values,
                stimulator,
                weights,
            )

        elif tcd_coil_element["type"] == 2:
            points = np.array(tcd_coil_element["points"])
            values = np.array(tcd_coil_element["values"])
            return LineSegmentElements(
                name,
                element_casing,
                element_deformations,
                points,
                values,
                stimulator,
                weights,
            )
        elif tcd_coil_element["type"] == 3:
            points = np.array(tcd_coil_element["points"])
            return LinePointElements(
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
            return SampledGridPointElements(
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


class PositionalTmsCoilElements(TmsCoilElements, ABC):
    """A representation of a positional stimulating element of a TMS coil

    Parameters
    ----------
    name : Optional[str]
        The name of the element
    casing : Optional[TmsCoilModel]
        The casing of the element
    deformations : Optional[list[TmsCoilDeformation]]
        A list of all deformations of the element
    points : npt.NDArray[np.float_]
        The positions of the stimulation elements
    stimulator : Optional[TmsStimulator]
        The stimulator used for this element
    weights : Optional[npt.NDArray[np.float_]], optional
        An additional weighting array with one scalar per element, by default None
    """

    def __init__(
        self,
        name: Optional[str],
        casing: Optional[TmsCoilModel],
        deformations: Optional[list[TmsCoilDeformation]],
        points: npt.NDArray[np.float_],
        stimulator: Optional[TmsStimulator],
        weights: Optional[npt.NDArray[np.float_]] = None,
    ):
        super().__init__(name, casing, deformations, stimulator, weights)
        self.points = points

    def get_points(
        self,
        affine_matrix: Optional[npt.NDArray[np.float_]] = None,
        apply_deformation: bool = True,
    ) -> npt.NDArray[np.float_]:
        """Returns the positions of the stimulation elements,
        optionally transformed by the affine matrix and deformed by the element deformation

        Parameters
        ----------
        affine_matrix : Optional[npt.NDArray[np.float_]], optional
            The affine transformation that is applied to the coil element, by default None
        apply_deformation : bool, optional
            Whether or not to apply the current coil element deformations, by default True

        Returns
        -------
        npt.NDArray[np.float_]
        The positions of the stimulation elements, optionally transformed by the affine matrix and deformed by the element deformation
        """
        if affine_matrix is None:
            affine_matrix = np.eye(4)
        if apply_deformation:
            affine_matrix = self.get_combined_transformation(affine_matrix)
        return self.points @ affine_matrix[:3, :3].T + affine_matrix[None, :3, 3]


class DirectionalTmsCoilElements(PositionalTmsCoilElements, ABC):
    """A representation of a directional stimulating element of a TMS coil

    Parameters
    ----------
    name : Optional[str]
        The name of the element
    casing : Optional[TmsCoilModel]
        The casing of the element
    deformations : Optional[list[TmsCoilDeformation]]
        A list of all deformations of the element
    points : npt.NDArray[np.float_]
        The positions of the stimulation elements
    values: npt.NDArray[np.float_]
        The values of the stimulation elements
    stimulator : Optional[TmsStimulator]
        The stimulator used for this element
    weights : Optional[npt.NDArray[np.float_]], optional
        An additional weighting array with one scalar per element, by default None
    """

    def __init__(
        self,
        name: Optional[str],
        casing: Optional[TmsCoilModel],
        deformations: Optional[list[TmsCoilDeformation]],
        points: npt.NDArray[np.float_],
        values: npt.NDArray[np.float_],
        stimulator: Optional[TmsStimulator],
        weights: Optional[npt.NDArray[np.float_]] = None,
    ):
        super().__init__(name, casing, deformations, points, stimulator, weights)
        self.values = values

    def get_values(
        self,
        affine_matrix: Optional[npt.NDArray[np.float_]] = None,
        apply_deformation: bool = True,
    ) -> npt.NDArray[np.float_]:
        """Returns the values of the stimulation elements,
        optionally transformed by the affine matrix and deformed by the element deformation

        Parameters
        ----------
        affine_matrix : Optional[npt.NDArray[np.float_]], optional
            The affine transformation that is applied to the coil element, by default None
        apply_deformation : bool, optional
            Whether or not to apply the current coil element deformations, by default True

        Returns
        -------
        npt.NDArray[np.float_]
            the values of the stimulation elements, optionally transformed by the affine matrix and deformed by the element deformation
        """
        if affine_matrix is None:
            affine_matrix = np.eye(4)
        if apply_deformation:
            affine_matrix = self.get_combined_transformation(affine_matrix)
        return self.values @ affine_matrix[:3, :3].T


class DipoleElements(DirectionalTmsCoilElements):
    def get_a_field(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        """Calculates the A field applied by the dipole elements at each target positions.

        Parameters
        ----------
        target_positions : npt.NDArray[np.float_] (N x 3)
            The points at which the A field should be calculated (in mm)
        coil_affine : npt.NDArray[np.float_] (4 x 4)
            The affine transformation that is applied to the coil
        eps : float, optional
            The requested precision, by default 1e-3

        Returns
        -------
        npt.NDArray[np.float_] (N x 3)
            The A field at every target positions in Tesla*meter
        """
        dipole_moment = self.get_values(coil_affine)
        dipole_position_m = self.get_points(coil_affine) * 1e-3
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

    def generate_element_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        apply_deformation: bool = True,
        element_tag: int = 0,
    ) -> Msh:
        """Generates a visualization of the sampled grid point elements visualized as lines as a mesh

        Parameters
        ----------
        affine_matrix : npt.NDArray[np.float_]
            The affine transformation that is applied to the coil element
        apply_deformation : bool, optional
            Whether or not to apply the current coil element deformations, by default True
        element_tag : int, optional
            The base value used as a tag for this coil element, by default 0

        Returns
        -------
        Msh
            The generated mesh representing the coil element
        """

        transformed_points = self.get_points(affine_matrix, apply_deformation)
        transformed_values = (
            transformed_points + self.get_values(affine_matrix, apply_deformation) * 1e3
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
        point_mesh.elm.tag1[:] = element_tag + TmsCoilElementTag.DIPOLES
        point_mesh.elm.tag2[:] = element_tag + TmsCoilElementTag.DIPOLES

        return point_mesh

    def to_tcd(
        self,
        stimulators: list[TmsStimulator],
        coil_models: list[TmsCoilModel],
        deformations: list[TmsCoilDeformation],
    ) -> dict:
        tcd_coil_element = super().to_tcd(stimulators, coil_models, deformations)
        tcd_coil_element["type"] = 1
        tcd_coil_element["points"] = self.points.tolist()
        tcd_coil_element["values"] = self.values.tolist()
        return tcd_coil_element


class LineSegmentElements(DirectionalTmsCoilElements):
    def get_a_field(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        """Calculates the A field applied by the line segment elements at each target positions.

        Parameters
        ----------
        target_positions : npt.NDArray[np.float_] (N x 3)
            The points at which the A field should be calculated (in mm)
        coil_affine : npt.NDArray[np.float_] (4 x 4)
            The affine transformation that is applied to the coil
        eps : float, optional
            The requested precision, by default 1e-3

        Returns
        -------
        npt.NDArray[np.float_] (N x 3)
            The A field at every target positions in Tesla*meter
        """
        directions_m = self.get_values(coil_affine) * 1e-3
        segment_position_m = self.get_points(coil_affine) * 1e-3
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

    def generate_element_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        apply_deformation: bool = True,
        element_tag: int = 0,
    ) -> Msh:
        """Generates a visualization of the line segment elements visualized as lines as a mesh

        Parameters
        ----------
        affine_matrix : npt.NDArray[np.float_]
            The affine transformation that is applied to the coil element
        apply_deformation : bool, optional
            Whether or not to apply the current coil element deformations, by default True
        element_tag : int, optional
            The base value used as a tag for this coil element, by default 0

        Returns
        -------
        Msh
            The generated mesh representing the coil element
        """
        transformed_points = self.get_points(affine_matrix, apply_deformation)
        transformed_values = transformed_points + self.get_values(
            affine_matrix, apply_deformation
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
        point_mesh.elm.tag1[:] = element_tag + TmsCoilElementTag.LINE_ELEMENTS
        point_mesh.elm.tag2[:] = element_tag + TmsCoilElementTag.LINE_ELEMENTS

        return point_mesh

    def to_tcd(
        self,
        stimulators: list[TmsStimulator],
        coil_models: list[TmsCoilModel],
        deformations: list[TmsCoilDeformation],
    ) -> dict:
        tcd_coil_element = super().to_tcd(stimulators, coil_models, deformations)
        tcd_coil_element["type"] = 2
        tcd_coil_element["points"] = self.points.tolist()
        tcd_coil_element["values"] = self.values.tolist()
        return tcd_coil_element


class LinePointElements(PositionalTmsCoilElements):
    def get_a_field(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        """Calculates the A field applied by the line point elements at each target positions.

        Parameters
        ----------
        target_positions : npt.NDArray[np.float_] (N x 3)
            The points at which the A field should be calculated (in mm)
        coil_affine : npt.NDArray[np.float_] (4 x 4)
            The affine transformation that is applied to the coil
        eps : float, optional
            The requested precision, by default 1e-3

        Returns
        -------
        npt.NDArray[np.float_] (N x 3)
            The A field at every target positions in Tesla*meter
        """
        segment_position_m = self.get_points(coil_affine) * 1e-3
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

    def generate_element_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        apply_deformation: bool = True,
        element_tag: int = 0,
    ) -> Msh:
        """Generates a visualization of the line point elements visualized as nodes as a mesh

        Parameters
        ----------
        affine_matrix : npt.NDArray[np.float_]
            The affine transformation that is applied to the coil element
        apply_deformation : bool, optional
            Whether or not to apply the current coil element deformations, by default True
        element_tag : int, optional
            The base value used as a tag for this coil element, by default 0

        Returns
        -------
        Msh
            The generated mesh representing the coil element
        """
        transformed_points = self.get_points(affine_matrix, apply_deformation)
        point_mesh = Msh(
            Nodes(transformed_points),
            Elements(points=np.arange(len(transformed_points)) + 1),
        )
        point_mesh.elm.tag1[:] = element_tag + TmsCoilElementTag.LINE_POINTS
        point_mesh.elm.tag2[:] = element_tag + TmsCoilElementTag.LINE_POINTS

        return point_mesh

    def to_tcd(
        self,
        stimulators: list[TmsStimulator],
        coil_models: list[TmsCoilModel],
        deformations: list[TmsCoilDeformation],
    ) -> dict:
        tcd_coil_element = super().to_tcd(stimulators, coil_models, deformations)
        tcd_coil_element["type"] = 3
        tcd_coil_element["points"] = self.points.tolist()
        return tcd_coil_element


class SampledGridPointElements(TmsCoilElements):
    def __init__(
        self,
        name: Optional[str],
        casing: Optional[TmsCoilModel],
        deformations: Optional[list[TmsCoilDeformation]],
        data: npt.NDArray[np.float_],
        affine: npt.NDArray[np.float_],
        stimulator: Optional[TmsStimulator],
        weights: Optional[npt.NDArray[np.float_]] = None,
    ):
        super().__init__(name, casing, deformations, stimulator, weights)
        self.data = data
        self.affine = affine

    def get_a_field(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        """Calculates the A field interpolated from the sampled grid point elements at each target positions.

        Parameters
        ----------
        target_positions : npt.NDArray[np.float_] (N x 3)
            The points at which the A field should be calculated (in mm)
        coil_affine : npt.NDArray[np.float_] (4 x 4)
            The affine transformation that is applied to the coil
        eps : float, optional
            The requested precision, by default 1e-3

        Returns
        -------
        npt.NDArray[np.float_] (N x 3)
            The A field at every target positions in Tesla*meter
        """
        combined_affine = self.get_combined_transformation(coil_affine)
        iM = np.linalg.pinv(self.affine) @ np.linalg.pinv(combined_affine)

        target_voxle_coordinates = (
            iM[:3, :3] @ target_positions.T + iM[:3, 3][:, np.newaxis]
        )

        # Interpolates the values of the field in the given coordinates
        out = np.zeros((3, target_voxle_coordinates.shape[1]))
        for dim in range(3):
            out[dim] = ndimage.map_coordinates(
                np.asanyarray(self.data)[..., dim], target_voxle_coordinates, order=1
            )

        # Rotates the field
        return (coil_affine[:3, :3] @ out).T

    def generate_element_mesh(
        self,
        affine_matrix: npt.NDArray[np.float_],
        apply_deformation: bool = True,
        element_tag: int = 0,
    ) -> Msh:
        """Generates a visualization of the sampled grid point elements visualized as lines as a mesh

        Parameters
        ----------
        affine_matrix : npt.NDArray[np.float_]
            The affine transformation that is applied to the coil element
        apply_deformation : bool, optional
            Whether or not to apply the current coil element deformations, by default True
        element_tag : int, optional
            The base value used as a tag for this coil element, by default 0

        Returns
        -------
        Msh
            The generated mesh representing the coil element
        """
        combined_affine = self.affine
        if apply_deformation:
            combined_affine = self.get_combined_transformation(
                affine_matrix) @ combined_affine 
        else:
            combined_affine = affine_matrix @ combined_affine

        voxel_coordinates = np.array(
            list(np.ndindex(self.data.shape[0], self.data.shape[1], self.data.shape[2]))
        )

        points = (
            voxel_coordinates @ combined_affine[:3, :3].T + combined_affine[None, :3, 3]
        )
        targets = points + self.data.reshape(-1, 3)

        print(targets)

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
        point_mesh.elm.tag1[:] = element_tag + TmsCoilElementTag.SAMPLED_GRID_ELEMENTS
        point_mesh.elm.tag2[:] = element_tag + TmsCoilElementTag.SAMPLED_GRID_ELEMENTS

        return point_mesh

    def to_tcd(
        self,
        stimulators: list[TmsStimulator],
        coil_models: list[TmsCoilModel],
        deformations: list[TmsCoilDeformation],
    ) -> dict:
        tcd_coil_element = super().to_tcd(stimulators, coil_models, deformations)
        tcd_coil_element["type"] = 4
        tcd_coil_element["data"] = self.data.tolist()
        tcd_coil_element["affine"] = self.affine.tolist()
        return tcd_coil_element

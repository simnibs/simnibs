from abc import ABC, abstractmethod
from typing import Optional

import fmm3dpy
import numpy as np
import numpy.typing as npt
from scipy.interpolate import LinearNDInterpolator

from simnibs.mesh_tools.mesh_io import Msh

from .coil_deformation import CoilDeformation
from .tms_stimulator import TMSStimulator


class CoilElement(ABC):
    def __init__(
        self,
        name: str,
        deformations: list[CoilDeformation],
        points: npt.NDArray[np.float_],
        stimulator: TMSStimulator,
        weights: Optional[npt.NDArray[np.float_]] = None,
    ):
        self.name = name
        self.deformations = deformations
        self.points = points
        self.stimulator = stimulator
        self.weights = weights

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
    ):
        if affine_matrix is None:
            affine_matrix = np.eye(4)
        for deformation in self.deformations:
            affine_matrix = deformation.as_matrix() @ affine_matrix

        return affine_matrix

    def get_transformed_points(
        self, affine_matrix: Optional[npt.NDArray[np.float_]] = None
    ):
        affine_matrix = self.get_combined_transformation(affine_matrix)
        return self.points @ affine_matrix[:3, :3].T + affine_matrix[None, :3, 3]


class DirectionalCoilElement(CoilElement, ABC):
    def __init__(
        self,
        name: str,
        deformations: list[CoilDeformation],
        points: npt.NDArray[np.float_],
        values: npt.NDArray[np.float_],
        stimulator: TMSStimulator,
        weights: Optional[npt.NDArray[np.float_]] = None,
    ):
        super().__init__(name, deformations, points, stimulator, weights)
        self.values = values

    @abstractmethod
    def get_a_field(self, msh: Msh, coil_affine: npt.NDArray[np.float_]):
        pass

    def get_transformed_values(
        self, affine_matrix: Optional[npt.NDArray[np.float_]] = None
    ):
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


class CoilLinePoints(CoilElement):
    def get_a_field(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        segment_position_m = self.get_transformed_points(coil_affine) * 1e-3
        target_positions_m = target_positions * 1e-3

        directions_m = np.zeros(segment_position_m.shape)
        directions_m[:,:-1] = np.diff(segment_position_m,axis=1)
        directions_m[:,-1] = segment_position_m[:,0]-segment_position_m[:,-1]

        if directions_m.shape[0] >= 300:
            A = fmm3dpy.l3ddir(
                sources=segment_position_m.T,
                charges=directions_m.T,
                targets=target_positions_m.T,
                nd=3,
                pgt=1,
            ).T
        else:
            A = fmm3dpy.lfmm3d(
                sources=segment_position_m.T,
                charges=directions_m.T,
                targets=target_positions_m.T,
                nd=3,
                eps=eps,
                pgt=1,
            ).T
        A = 1e-7 * A.pottarg
        return A


class CoilSampledElements(DirectionalCoilElement):
    def get_a_field(
        self,
        target_positions: npt.NDArray[np.float_],
        coil_affine: npt.NDArray[np.float_],
        eps: float = 1e-3,
    ) -> npt.NDArray[np.float_]:
        values = self.get_transformed_values(coil_affine)
        positions = self.get_transformed_points(coil_affine)

        interpolation = LinearNDInterpolator(positions, values, fill_value=0) 

        return interpolation(target_positions)

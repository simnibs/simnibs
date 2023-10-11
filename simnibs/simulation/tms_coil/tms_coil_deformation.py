from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt
from scipy.spatial.transform import Rotation

from simnibs.simulation.tms_coil.tcd_element import TcdElement


class TmsCoilDeformationRange(TcdElement):
    """Represents a deformation range

    Parameters
    ----------
    initial : float
        The initial value of the deformation (e.g angle of a rotation deformation)
    range : npt.ArrayLike
        The allowed range of the deformation ([min, max])

    Attributes
    ----------------------
    initial : float
        The initial value of the deformation (e.g angle of a rotation deformation)
    range : npt.NDArray[np.float64]
        The allowed range of the deformation ([min, max])
    current : float
        The current value of the deformation
    """

    def __init__(self, initial: float, range: npt.ArrayLike):
        self.initial = initial
        self._current = initial
        self.range = np.array(range, dtype=np.float64)

        if self.range.ndim != 1 or self.range.shape[0] != 2:
            raise ValueError(
                f"Expected 'range' to be in the format [min, max] ({self.range})"
            )
        elif self.range[0] > self.range[1]:
            raise ValueError(
                f"Expected 'range' to be in the format [min, max] but min was greater than max ({self.range})"
            )

        if self.initial < self.range[0] or self.initial > self.range[1]:
            raise ValueError(
                f"Expected 'initial' to be in the range of 'range' ({initial}, {self.range})"
            )

    def reset(self):
        self.current = self.initial

    @property
    def current(self) -> float:
        return self._current

    @current.setter
    def current(self, value):
        if value < self.range[0] or value > self.range[1]:
            raise ValueError(f"Value must be within the range ({value}, {self.range})")
        else:
            self._current = value

    def to_tcd(self, ascii_mode: bool = False) -> dict:
        tcd_deformation = {}
        tcd_deformation["initial"] = self.initial
        tcd_deformation["range"] = list(self.range)
        return tcd_deformation

    @classmethod
    def from_tcd(cls, tcd_deformation: dict):
        initial = tcd_deformation["initial"]
        deform_range = tcd_deformation["range"]

        return TmsCoilDeformationRange(initial, deform_range)


class TmsCoilDeformation(ABC, TcdElement):
    """Represents a deformation

    Parameters
    ----------
    deformation_range : TmsCoilDeformationRange
            The range, initial and current value of the deformation

    Attributes
    ----------------------
    deformation_range : TmsCoilDeformationRange
            The range, initial and current value of the deformation
    """

    def __init__(self, deformation_range: TmsCoilDeformationRange):
        self.deformation_range = deformation_range

    @abstractmethod
    def apply(self, points: npt.NDArray[np.float_]) -> npt.NDArray[np.float_]:
        """Applies the deformation to every point

        Parameters
        ----------
        points : npt.NDArray[np.float_] (Nx3)
            The points that should be deformed


        Returns
        -------
        npt.NDArray[np.float_] (Nx3)
            The deformed points
        """
        pass

    @abstractmethod
    def as_matrix(self) -> npt.NDArray[np.float_]:
        """An affine matrix representation of the deformation

        Returns
        -------
        npt.NDArray[np.float_] (4x4)
            The affine matrix representing the deformation
        """
        pass

    def to_tcd(
        self,
        deformation_ranges: list[TmsCoilDeformationRange],
        ascii_mode: bool = False,
    ) -> dict:
        tcd_deformation = {}
        tcd_deformation["deformRange"] = deformation_ranges.index(
            self.deformation_range
        )
        return tcd_deformation

    @classmethod
    def from_tcd(
        cls, tcd_deformation: dict, deformation_ranges: list[TmsCoilDeformationRange]
    ):
        deformation_range = deformation_ranges[tcd_deformation["deformRange"]]
        if tcd_deformation["type"] == "x":
            return TmsCoilTranslation(deformation_range, 0)
        elif tcd_deformation["type"] == "y":
            return TmsCoilTranslation(deformation_range, 1)
        elif tcd_deformation["type"] == "z":
            return TmsCoilTranslation(deformation_range, 2)
        elif tcd_deformation["type"] == "rot2p":
            point_1 = np.array(tcd_deformation["point1"])
            point_2 = np.array(tcd_deformation["point2"])
            return TmsCoilRotation(deformation_range, point_1, point_2)
        else:
            raise ValueError(f"Invalid deformation type: {tcd_deformation['type']}")


class TmsCoilTranslation(TmsCoilDeformation):
    """Represents a translation

    Parameters
    ----------
    deformation_range : TmsCoilDeformationRange
            The range, initial and current value of the deformation
    axis : int
        The axis to be used for the translation

    Attributes
    ----------------------
    deformation_range : TmsCoilDeformationRange
            The range, initial and current value of the deformation
    current : float
        The current value of the deformation
    axis : int
        The axis to be used for the translation
    """

    def __init__(
        self,
        deformation_range: TmsCoilDeformationRange,
        axis: int,
    ):
        super().__init__(deformation_range)
        self.axis = axis

        if self.axis < 0 or self.axis > 2:
            raise ValueError(
                f"Expected 'axis' to be between 0 and 2 but was {self.axis}"
            )

    def apply(self, points: npt.NDArray[np.float_]):
        return points + self.get_translation()

    def get_translation(self) -> npt.NDArray[np.float_]:
        """Returns the translation vector for this coil translation

        Returns
        -------
        npt.NDArray[np.float_] (3)
            The translation vector for this coil translation
        """
        translation = np.zeros(3)
        translation[self.axis] = self.deformation_range.current
        return translation

    def as_matrix(self) -> npt.NDArray[np.float_]:
        affine_matrix = np.eye(4)
        affine_matrix[:3, 3] = self.get_translation()
        return affine_matrix

    def to_tcd(
        self,
        deformation_ranges: list[TmsCoilDeformationRange],
        ascii_mode: bool = False,
    ) -> dict:
        tcd_deformation = super().to_tcd(deformation_ranges, ascii_mode)
        if self.axis == 0:
            tcd_deformation["type"] = "x"
        elif self.axis == 1:
            tcd_deformation["type"] = "y"
        elif self.axis == 2:
            tcd_deformation["type"] = "z"
        else:
            raise ValueError(f"Translation axis ({self.axis}) out of range (0-2)")
        return tcd_deformation


class TmsCoilRotation(TmsCoilDeformation):
    """Represents a rotation around an axis defined by two points

    Parameters
    ----------
    deformation_range : TmsCoilDeformationRange
            The range, initial and current value of the deformation
    point_1: npt.ArrayLike
        The first point of the rotation axis
    point_2: npt.ArrayLike
        The second point of the rotation axis

    Attributes
    ----------------------
    deformation_range : TmsCoilDeformationRange
            The range, initial and current value of the deformation
    current : float
        The current value of the deformation
    point_1: npt.NDArray[np.float64]
        The first point of the rotation axis
    point_2: npt.NDArray[np.float64]
        The second point of the rotation axis
    """

    def __init__(
        self,
        deformation_range: TmsCoilDeformationRange,
        point_1: npt.ArrayLike,
        point_2: npt.ArrayLike,
    ):
        super().__init__(deformation_range)
        self.point_1: npt.NDArray[np.float64] = np.array(point_1, dtype=np.float64)
        self.point_2: npt.NDArray[np.float64] = np.array(point_2, dtype=np.float64)

        if self.point_1.ndim != 1 or self.point_1.shape[0] != 3:
            raise ValueError(
                f"Expected 'point_1' to be in the format [x, y, z] but shape was {self.point_1.shape}"
            )

        if self.point_2.ndim != 1 or self.point_2.shape[0] != 3:
            raise ValueError(
                f"Expected 'point_2' to be in the format [x, y, z] but shape was {self.point_2.shape}"
            )

    def get_rotation(self) -> npt.NDArray[np.float_]:
        """Returns the affine matrix of the rotation around the axis defined by the two points

        Returns
        -------
        npt.NDArray[np.float_] (4x4)
            The affine matrix representing the rotation around the axis defined by the two points
        """
        v = (self.point_2 - self.point_1) / np.linalg.norm(self.point_2 - self.point_1)
        T = np.identity(4)
        T[:3, 3] = self.point_1
        iT = np.identity(4)
        iT[:3, 3] = -self.point_1
        R = np.identity(4)
        R[:3, :3] = Rotation.from_rotvec(
            v * self.deformation_range.current, degrees=True
        ).as_matrix()
        Q = T @ R @ iT
        return Q

    def apply(self, points: npt.NDArray[np.float_]):
        rotation_matrix = self.get_rotation()
        points = points @ rotation_matrix[:3, :3].T + rotation_matrix[None, :3, 3]
        return points

    def as_matrix(self) -> npt.NDArray[np.float_]:
        return self.get_rotation()

    def to_tcd(
        self,
        deformation_ranges: list[TmsCoilDeformationRange],
        ascii_mode: bool = False,
    ) -> dict:
        tcd_deformation = super().to_tcd(deformation_ranges, ascii_mode)
        tcd_deformation["type"] = "rot2p"
        tcd_deformation["point1"] = self.point_1.tolist()
        tcd_deformation["point2"] = self.point_2.tolist()
        return tcd_deformation

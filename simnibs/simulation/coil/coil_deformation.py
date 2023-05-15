from abc import ABC, abstractmethod

import numpy as np
import numpy.typing as npt
from scipy.spatial.transform import Rotation

class CoilDeformation(ABC):
    def __init__(self, initial: float, range: tuple[float, float]):
        self.initial = initial
        self._current = initial
        self.range = range

    def reset(self):
        self.current = self.initial

    @property
    def current(self):
        return self._current  
    
    @current.setter
    def current(self, value):
        if self.current < self.range[0] or self.current > self.range[1]:
            raise ValueError(f"Value must be within the range ({self.range[0]}, {self.range[1]})")
        else:
            self._current = value  

    @abstractmethod
    def apply(self, points:npt.NDArray[np.float_]):
        pass

    @abstractmethod
    def as_matrix(self) -> npt.NDArray[np.float_]:
        pass


class CoilTranslation(CoilDeformation):
    def __init__(
        self,
        initial: float,
        range: tuple[float, float],
        axis: int,
    ):
        super().__init__(initial, range)
        self.axis = axis

    def apply(self, points:npt.NDArray[np.float_]):
        return points + self.get_translation()
    
    def get_translation(self) -> npt.NDArray[np.float_]:
        translation = np.zeros(3)
        translation[self.axis] = self.current
        return translation
    
    def as_matrix(self) -> npt.NDArray[np.float_]:
        affine_matrix = np.eye(4)
        affine_matrix[:3, 3] = self.get_translation()
        return affine_matrix


class CoilRotation(CoilDeformation):
    def __init__(
        self,
        initial: float,
        range: tuple[float, float],
        point_1: npt.NDArray[np.float_],
        point_2: npt.NDArray[np.float_],
    ):
        super().__init__(initial, range)
        self.point_1 = point_1
        self.point_2 = point_2

    def get_rotation(self) -> Rotation:
        v = (self.point_2 - self.point_1) / np.linalg.norm(self.point_2 - self.point_1)
        rotation = Rotation.from_rotvec(self.current * v, degrees=True)
        #TODO Check correctness 
        return rotation

    def apply(self, points:npt.NDArray[np.float_]):
        return self.get_rotation().apply(points)

    def as_matrix(self) -> npt.NDArray[np.float_]:
        rotation_matrix = self.get_rotation().as_matrix()
        affine_matrix = np.eye(4)
        affine_matrix[:3, :3] = rotation_matrix
        return affine_matrix
from ...mesh_tools.mesh_io import Msh

import numpy as np
import numpy.typing as npt

class CoilModel:
    def __init__(
        self, mesh : Msh, 
        min_distance_points : npt.NDArray[np.float_], 
        intersect_points : npt.NDArray[np.float_]
    ):
        self.mesh = mesh
        self.min_distance_points = min_distance_points
        self.intersect_points = intersect_points
       
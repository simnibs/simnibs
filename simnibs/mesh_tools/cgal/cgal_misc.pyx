# distutils: language = c++
# cython: language_level=3
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.math cimport abs
import cython
import numpy as np
cimport numpy as np


cdef extern from "_cgal_intersect.cpp" nogil:
    pair[vector[int], vector[float]] _segment_triangle_intersection(
        float* vertices, int n_vertices, int* tris, int n_faces,
        float* segment_start, float* segment_end, int n_segments)


def segment_triangle_intersection(vertices, faces, segment_start, segment_end):
    ''' Calculates the intersection between a triangular mesh and line segments
    '''
    cdef np.ndarray[float] vert = np.ascontiguousarray(vertices, dtype=np.float32).reshape(-1)
    cdef np.ndarray[int] fac = np.ascontiguousarray(faces, dtype=np.int32).reshape(-1)
    cdef np.ndarray[float] ss = np.ascontiguousarray(segment_start, dtype=np.float32).reshape(-1)
    cdef np.ndarray[float] se = np.ascontiguousarray(segment_end, dtype=np.float32).reshape(-1)
    cdef pair[vector[int], vector[float]] out

    out = _segment_triangle_intersection(
        &vert[0], len(vertices), &fac[0], len(faces),
        &ss[0], &se[0], len(segment_start)
    )

    pairs = np.array(out.first, dtype=int).reshape(-1, 2)
    positions = np.array(out.second, dtype=float).reshape(-1, 3)

    return pairs, positions

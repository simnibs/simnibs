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
    cdef cppclass TreeC:
        TreeC(float* vertices, int n_vertices, int* tris, int n_faces) except+
        TreeC() except+
        pair[vector[int], vector[float]] _intersections(
            float* segment_start, float* segment_end, int n_segments)
        bool _any_intersections(
            float* segment_start, float* segment_end, int n_segments)
        bool _any_point_inside(
            float* pts, int n_points)
        vector[int] _points_inside(
            float* pts, int n_points)

cdef class pyAABBTree:
    cdef TreeC *thisptr
    def ___cinit___(self):
        self.thisptr = new TreeC()

    def __dealloc___(self):
        del self.thisptr

    def __del__(self):
        self.__dealloc___()

    def set_data(self,vertices, faces):
        cdef np.ndarray[float] vert = np.ascontiguousarray(vertices, dtype=np.float32).reshape(-1)
        cdef np.ndarray[int] fac = np.ascontiguousarray(faces, dtype=np.int32).reshape(-1)
        del self.thisptr 
        self.thisptr = new TreeC(&vert[0], len(vertices), &fac[0], len(faces))

    def intersection(self,segment_start, segment_end):
        cdef np.ndarray[float] ss = np.ascontiguousarray(segment_start, dtype=np.float32).reshape(-1)
        cdef np.ndarray[float] se = np.ascontiguousarray(segment_end, dtype=np.float32).reshape(-1)
        cdef pair[vector[int], vector[float]] out
        out = self.thisptr._intersections(&ss[0], &se[0], len(segment_start))
        pairs = np.array(out.first, dtype=int).reshape(-1, 2)
        positions = np.array(out.second, dtype=float).reshape(-1, 3)
        return pairs, positions
    
    def any_intersection(self, segment_start, segment_end):
        cdef np.ndarray[float] ss = np.ascontiguousarray(segment_start, dtype=np.float32).reshape(-1)
        cdef np.ndarray[float] se = np.ascontiguousarray(segment_end, dtype=np.float32).reshape(-1)
        out = self.thisptr._any_intersections(&ss[0], &se[0], len(segment_start))
        return out

    def any_point_inside(self, points):
        cdef np.ndarray[float] pts = np.ascontiguousarray(points, dtype=np.float32).reshape(-1)
        out = self.thisptr._any_point_inside(&pts[0], len(points))
        return out
    
    def points_inside(self, points):
        cdef np.ndarray[float] pts = np.ascontiguousarray(points, dtype=np.float32).reshape(-1)
        cdef vector[int] out
        out = self.thisptr._points_inside(&pts[0], len(points))
        return out


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

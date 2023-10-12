# distutils: language = c++
# cython: language_level=3
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.math cimport abs
import cython
import numpy as np
cimport numpy as np


cdef extern from "_mesh_volumes.cpp" nogil:
    int _mesh_image(
        char *fn_image,
        char *fn_out,
        float facet_angle,
        float facet_size,
        float facet_distance,
        float cell_radius_edge_ratio,
        float cell_size,
        int num_threads,
        bool do_perturb,
        bool do_exude,
        bool do_lloyd
    )
    int _mesh_image_sizing_field(
        char *fn_image,
        char *fn_out,
        float facet_angle,
        float *facet_size,
        float *facet_distance,
        float cell_radius_edge_ratio,
        float *cell_size,
        int num_threads,
        bool do_perturb,
        bool do_exude,
        bool do_lloyd
    )


def mesh_image(
        fn_image,
        fn_out,
        facet_angle: float,
        facet_size: float,
        facet_distance: float,
        cell_radius_edge_ratio: float,
        cell_size: float,
        num_threads: int,
        do_perturb: bool = False,
        do_exude: bool = False,
        do_lloyd: bool = False,
    ):
    return _mesh_image(
        fn_image, fn_out, facet_angle,
        facet_size, facet_distance,
        cell_radius_edge_ratio, cell_size,
        num_threads, do_perturb, do_exude, do_lloyd,
    )


def mesh_image_sizing_field(
        fn_image,
        fn_out,
        facet_angle: float,
        facet_size,
        facet_distance,
        cell_radius_edge_ratio: float,
        cell_size,
        num_threads: int,
        do_perturb: bool = False,
        do_exude: bool = False,
        do_lloyd: bool = False,
    ):

    cdef np.ndarray[float, ndim=3] sf_facet_size = np.array(
        facet_size, dtype=np.float32, order='F', copy=False
    )
    cdef np.ndarray[float, ndim=3] sf_cell_size = np.array(
        cell_size, dtype=np.float32, order='F', copy=False
    )
    cdef np.ndarray[float, ndim=3] sf_facet_distance = np.array(
        facet_distance, dtype=np.float32, order='F', copy=False
    )

    return _mesh_image_sizing_field(
        fn_image, fn_out, facet_angle,
        &sf_facet_size[0, 0, 0],
        &sf_facet_distance[0, 0, 0],
        cell_radius_edge_ratio,
        &sf_cell_size[0, 0, 0],
        num_threads, do_perturb, do_exude, do_lloyd,
    )



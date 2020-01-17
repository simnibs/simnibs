# distutils: language = c++
# cython: language_level=3

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.math cimport abs
import cython
import numpy as np
cimport numpy as np


cdef extern from "_mesh.cpp" nogil:
    int _mesh_image(
        char *fn_image, char *fn_out, float facet_angle,
        float facet_size, float facet_distance,
        float cell_radius_edge_ratio, float cell_size,
        bool optimize)
    int _mesh_surfaces(
        vector[char *]filenames, vector[pair[int, int]] incident_subdomains,
        char *fn_out,
        float facet_angle, float facet_size, float facet_distance,
        float cell_radius_edge_ratio, float cell_size,
        bool optimize)

def mesh_image(fn_image, fn_out, float facet_angle, float facet_size,
               float facet_distance, float cell_radius_edge_ratio, float cell_size,
               bool optimize):
    ret =  _mesh_image(
        fn_image, fn_out, facet_angle,
        facet_size, facet_distance,
        cell_radius_edge_ratio, cell_size,
        optimize
    )
    return ret

def mesh_surfaces(fn_surfaces, incident_subdomains, fn_out,
                  float facet_angle, float facet_size,
                  float facet_distance, float cell_radius_edge_ratio, float cell_size,
                  bool optimize):
    cdef vector[pair[int, int]] subdomains_pairs = incident_subdomains
    cdef vector[char *] surfaces = fn_surfaces
    ret  = _mesh_surfaces(
      surfaces, subdomains_pairs, fn_out,
      facet_angle, facet_size, facet_distance,
      cell_radius_edge_ratio, cell_size,
      optimize
    )
    return ret


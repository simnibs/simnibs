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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def gauss_smooth(
    np.uint_t[::1] surf_nodes,
    np.float_t[:, ::1] nodes_pos,
    np.uint_t[:, ::1] tetrahedra,
    np.uint_t[::1] adj_indices,
    np.uint_t[::1] adj_indptr,
    np.uint_t[::1] adj_th_indices,
    np.uint_t[::1] adj_th_indptr,
    float factor):
    ''' Gaussian smoothing in-place
    checks for violations if the volume of the node stays the same before and after the
    move
    '''
    cdef int count_cancelled = 0
    cdef np.uint_t n
    cdef np.float_t[3] bar
    cdef np.float_t[::1] pos_before
    cdef Py_ssize_t i, j, k

    # for speed, I dont use numpy functions
    for i in range(len(surf_nodes)):
        n = surf_nodes[i]
        for k in range(3):
            bar[k] = 0.

        for j in range(adj_indptr[n], adj_indptr[n+1]):
            for k in range(3):
                bar[k] = bar[k] + nodes_pos[adj_indices[j], k]

        for k in range(3):
            bar[k] = bar[k]/float(adj_indptr[n+1] - adj_indptr[n])

        pos_before = nodes_pos[n].copy()

        for k in range(3):
            nodes_pos[n, k] += factor * (bar[k] - nodes_pos[n, k])

        for j in range(adj_th_indptr[n], adj_th_indptr[n+1]):
            if not _test_sign(nodes_pos, tetrahedra[adj_th_indices[j]]):
                nodes_pos[n] = pos_before
                count_cancelled += 1
                break

    print('{0} moves cancelled'.format(count_cancelled))



@cython.boundscheck(False)
@cython.wraparound(False)
cdef bool _test_sign(
    np.float_t[:, ::1] nodes_pos,
    np.uint_t[:] th):
    cdef np.float_t[3][3] M
    cdef Py_ssize_t i, j, k
    cdef np.float_t det
    for i in range(3):
        for j in range(3):
            M[i][j] = nodes_pos[th[i], j] - nodes_pos[th[3], j]
    
    det = M[0][0]*(M[1][1]*M[2][2] - M[2][1]*M[1][2]) - \
          M[0][1]*(M[1][0]*M[2][2] - M[2][0]*M[1][2]) + \
          M[0][2]*(M[1][0]*M[2][1] - M[2][0]*M[1][1])
    return det < 0

# Only for testing
def test_sign(
    np.float_t[:, ::1] nodes_pos,
    np.uint_t[:] th):
    return _test_sign(nodes_pos, th)


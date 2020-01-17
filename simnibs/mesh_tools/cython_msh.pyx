# cython: language_level=2

from __future__ import division
import numpy as np
import cython
import scipy
cimport numpy as np
from libcpp cimport bool
from libc.math cimport exp
from libc.math cimport abs
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

@cython.boundscheck(False)
@cython.wraparound(False)
def interp_grid(np.ndarray[np.int_t, ndim=1] n_voxels,
                np.ndarray[np.float_t, ndim=2] field,
                np.ndarray[np.float_t, ndim=2] nd,
                np.ndarray[np.int_t, ndim=2] tetrahedra):
    # image
    cdef np.ndarray[np.float_t, ndim=4] image = np.zeros((n_voxels[0], n_voxels[1],
                                                      n_voxels[2], field.shape[1]), np.double)

    cdef np.int_t nr_components = field.shape[1]
    cdef np.int_t node_data = field.shape[0] == nd.shape[0]
    ## Create bounding box with each tetrahedra
    cdef np.ndarray[np.float_t, ndim=3] th_coords = nd[tetrahedra]

    cdef np.ndarray[np.float_t, ndim=3] invM = \
            np.linalg.inv(
                np.transpose(th_coords[:, :3, :3] - th_coords[:, 3, None, :], (0, 2, 1)))

    cdef np.ndarray[np.int_t, ndim=2] th_boxes_min = np.rint(
        np.min(th_coords, axis=1)).astype(np.int)

    cdef np.ndarray[np.int_t, ndim=2] th_boxes_max = np.rint(
        np.max(th_coords, axis=1)).astype(np.int)

    cdef np.ndarray[np.int_t, ndim=1] in_roi = np.where(
        np.all((th_boxes_min <= n_voxels) * (th_boxes_max >= 0), axis=1))[0].astype(np.int)

    th_boxes_max = np.minimum(th_boxes_max, np.array(n_voxels) - 1)
    th_boxes_min = np.maximum(th_boxes_min, 0)
    # pre-calculate the inverse of M (this is faster than solving every M[j])
    #invM[in_roi] = np.linalg.inv(M[in_roi])
 
    cdef int i, j, k, x, y, z, info
    cdef np.ndarray[np.float_t, ndim=1] b = np.zeros((4, ), dtype=np.float)
    cdef np.float_t xc, yc, zc
    cdef np.float_t eps = 1e-5

    for j in in_roi:
        for x in range(th_boxes_min[j, 0], th_boxes_max[j, 0] + 1):
            xc = x - th_coords[j, 3, 0]
            for y in range(th_boxes_min[j, 1], th_boxes_max[j, 1] + 1):
                yc = y - th_coords[j, 3, 1]
                for z in range(th_boxes_min[j, 2], th_boxes_max[j, 2] + 1):
                    zc = z - th_coords[j, 3, 2]
                    b[0] = invM[j, 0, 0] * xc + \
                           invM[j, 0, 1] * yc + \
                           invM[j, 0, 2] * zc
                    if b[0] > -eps and b[0] < 1. + eps:
                        b[1] = invM[j, 1, 0] * xc + \
                               invM[j, 1, 1] * yc + \
                               invM[j, 1, 2] * zc
                        if b[1] > -eps and b[0] + b[1] < 1. + eps:
                            b[2] = invM[j, 2, 0] * xc + \
                                   invM[j, 2, 1] * yc + \
                                   invM[j, 2, 2] * zc
                            b[3] = 1. - b[0] - b[1] - b[2]
                            if b[2] > -eps and b[3] > -eps:
                                for k in range(nr_components):
                                    if node_data:
                                        image[x, y, z, k] = 0.
                                        for i in range(4):
                                            image[x, y, z, k] += b[i] * field[tetrahedra[j, i], k]
                                    else:
                                        image[x, y, z, k] = field[j, k]
                   
    del invM
    del th_boxes_min
    del th_boxes_max
    del in_roi
    del th_coords
    return image

@cython.boundscheck(False)
@cython.wraparound(False)
def find_tetrahedron_with_points(np.ndarray[np.float_t, ndim=2] points,
                                 np.ndarray[np.float_t, ndim=3] th_nodes,
                                 np.ndarray[np.int_t, ndim=1] starting_th,
                                 np.ndarray[np.int_t, ndim=2] th_faces,
                                 np.ndarray[np.int_t, ndim=2] adjacency_list):

    cdef np.ndarray[np.int_t, ndim=2] face_points = np.array(
        [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]], np.int)
    cdef np.ndarray[np.int_t, ndim=1] th_with_points = -np.ones(points.shape[0],
                                                                dtype=np.int)
    # We can now start the walking algorithm
    cdef np.ndarray[np.int_t, ndim=1] face_order = np.arange(4, dtype=np.int)
    cdef np.ndarray[np.float_t, ndim=1] p
    cdef np.int_t previous_t, adjacent, nr_cycles
    cdef np.uint_t end, outside, t, j, face
    cdef np.int_t pt = len(th_faces) + 1

    for i, p, t in zip(range(points.shape[0]), points, starting_th):
        previous_t = pt
        end = 0
        nr_cycles = 0
        while not end:
            #randomize face order
            np.random.shuffle(face_order)
            end = 1
            outside = 0
            # For each faces
            for j in range(4):
                f = face_order[j]
                face = th_faces[t, f]
                # Calculate where point lies in relation to face
                # If it is in the other side of the faces
                if orientation(p, th_nodes[t], face_points[f]) < 0.:
                    # See which triangle is adjacent to the current triangle throught
                    # that face
                    adjacent = adjacency_list[face, 0]
                    if adjacent == t:
                        adjacent = adjacency_list[face, 1]
                    # if the face has no adjacent tetrahedra
                    if adjacent == -1:
                        outside = 1
                        # We will only say that it is really outisde when we find a
                        # tetrahedron were p is in the same side of 3 faces, but in the
                        # other side of another face, wich points outwards
                    # If this is not the triangle we just came from
                    elif adjacent != previous_t:
                        # Move to the adjacent triangle
                        previous_t = t
                        t = adjacent
                        end = 0
                        break
            # this is only here for ensure that it will not loop forever
            nr_cycles += 1
            if nr_cycles >= 1000:
                outside = 1
                break
        if not outside:
            th_with_points[i] = t
        if outside:
            th_with_points[i] = -1

    return th_with_points

'''
@cython.boundscheck(False)
@cython.wraparound(False)
def orientation(np.ndarray[np.float_t, ndim=1] p, int th, int face,
                np.ndarray[np.float_t, ndim=3] th_nodes,
                np.ndarray[np.int_t, ndim=2] face_points):
    cdef np.ndarray[np.float_t, ndim=2] d = np.empty((3, 3), dtype=np.float)
    for i in range(3):
        for j in range(3):
            d[i, j] = th_nodes[th, face_points[face][j], i] - p[i]
    cdef np.float_t det = 0.
    det += d[0, 0] * d[1, 1] * d[2, 2]
    det += d[0, 1] * d[1, 2] * d[2, 0]
    det += d[0, 2] * d[1, 0] * d[2, 1]
    det -= d[2, 0] * d[1, 1] * d[0, 2]
    det -= d[2, 1] * d[1, 2] * d[0, 0]
    det -= d[2, 2] * d[1, 0] * d[0, 1]
    return np.sign(det)
'''
@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.float_t orientation(np.float_t[:] p,
                            np.float_t[:, :] th_nodes,
                            np.int_t[:] face_points):
    #cdef np.ndarray[np.float_t, ndim=2] d = np.empty((3, 3), dtype=np.float)
    cdef np.float_t[3][3] d
    for i in range(3):
        for j in range(3):
            d[i][j] = th_nodes[face_points[j], i] - p[i]
            #d[i, j] = th_nodes[face_points[j], i] - p[i]
    cdef np.float_t det = 0.
    det += d[0][0] * d[1][1] * d[2][2]
    det += d[0][1] * d[1][2] * d[2][0]
    det += d[0][2] * d[1][0] * d[2][1]
    det -= d[2][0] * d[1][1] * d[0][2]
    det -= d[2][1] * d[1][2] * d[0][0]
    det -= d[2][2] * d[1][0] * d[0][1]
    return (det > 0) - (det < 0)#np.sign(det)


@cython.boundscheck(False)
@cython.wraparound(False)
def test_point_in_triangle(np.ndarray[np.float_t, ndim=2] points,
                           np.ndarray[np.float_t, ndim=2] triangle0,
                           np.ndarray[np.float_t, ndim=2] edge0,
                           np.ndarray[np.float_t, ndim=2] edge1,
                           np.ndarray[np.float_t, ndim=2] bb_max,
                           np.ndarray[np.float_t, ndim=2] bb_min,
                           np.ndarray[np.float_t, ndim=1] inv_det,
                           np.ndarray[np.uint8_t, ndim=1] det_zero,
                           np.ndarray[np.float_t, ndim=2] P):

    cdef float eps = 1e-5
    cdef np.ndarray[np.float_t, ndim=1] ray_direction = \
            np.array([1, 0, 0], dtype=np.float)
    cdef np.ndarray[np.uint8_t, ndim=1] inside = np.zeros(len(points), dtype=np.uint8)
    cdef np.ndarray[np.uint8_t, ndim=1] can_cross
    cdef np.ndarray[np.float_t, ndim=2] T, Q
    cdef np.ndarray[np.float_t, ndim=1] p, u, v, q, t
    cdef int n_crosses
    for i in range(len(points)):
        p = points[i]
        can_cross = np.ones(len(triangle0), dtype=np.uint8)
        can_cross[det_zero.view(np.bool)] = 0
        can_cross[np.any(bb_max[:, 1:] < p[1:], axis=1)] = 0
        can_cross[np.any(bb_min[:, 1:] > p[1:], axis=1)] = 0
        can_cross[bb_max[:, 0] < p[0]] = 0
        T = p - triangle0[can_cross.view(np.bool)]
        u = (T * P[can_cross.view(np.bool)]).sum(axis=1) * \
            inv_det[can_cross.view(np.bool)]
        Q = np.cross(T, edge0[can_cross.view(np.bool)])
        v = np.dot(ray_direction, Q.T) * inv_det[can_cross.view(np.bool)]
        q = 1 - v - u
        t = (edge1[can_cross.view(np.bool)] * Q).sum(axis=1) * \
            inv_det[can_cross.view(np.bool)]
        n_crosses = np.sum(
            (u > 0.0 - eps) * (u < 1.0 + eps) * \
            (v > 0.0 - eps) * (q > 0.0 - eps) * \
            (t > eps))
        if n_crosses % 2 == 1:
            inside[i] = 1

    return inside.view(np.bool)


def calc_quantities_for_test_point_in_triangle(triangles):
    cdef float eps = 1e-5
    edge0 = triangles[:, 1, :] - triangles[:, 0, :]
    edge1 = triangles[:, 2, :] - triangles[:, 0, :]
    # Ray direcion is x
    ray_direction = np.array([1., 0., 0.])
    P = np.cross(ray_direction, edge1)
    det = (edge0 * P).sum(axis=1)
    det_zero = np.abs(det) < eps
    bb_min = np.min(triangles, axis=1)
    bb_max = np.max(triangles, axis=1)
    inv_det = np.ones(len(triangles), dtype=float)
    inv_det[~det_zero] = 1.0 / det[~det_zero]
    
    triangle0 = np.array(triangles[:, 0, :], dtype=np.float)
    edge0 = np.array(edge0, dtype=np.float)
    edge1 = np.array(edge1, dtype=np.float)
    bb_max = np.array(bb_max, dtype=np.float)
    bb_min = np.array(bb_min, dtype=np.float)
    inv_det = np.array(inv_det, dtype=np.float)
    det_zero = np.array(det_zero, dtype=np.uint8)
    P = np.array(P, np.float)

    return triangle0, edge0, edge1, bb_max, bb_min, inv_det, det_zero, P

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
    float factor,
    np.uint_t[::1] nodes_mask
    ):
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
        if not nodes_mask[i]:
            # do not smooth
            continue 
        n = surf_nodes[i]
        for k in range(3):
            bar[k] = 0.

        # move
        for j in range(adj_indptr[n], adj_indptr[n+1]):
            for k in range(3):
                bar[k] = bar[k] + nodes_pos[adj_indices[j], k]

        for k in range(3):
            bar[k] = bar[k]/float(adj_indptr[n+1] - adj_indptr[n])

        pos_before = nodes_pos[n].copy()

        for k in range(3):
            nodes_pos[n, k] += factor * (bar[k] - nodes_pos[n, k])

        # check
        for j in range(adj_th_indptr[n], adj_th_indptr[n+1]):
            if not _test_sign(nodes_pos, tetrahedra[adj_th_indices[j]]):
                nodes_pos[n] = pos_before
                count_cancelled += 1
                break

    return count_cancelled


@cython.boundscheck(False)
@cython.wraparound(False)
cdef int _test_sign(
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


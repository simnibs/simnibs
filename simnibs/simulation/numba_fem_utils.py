# -*- coding: utf-8 -*-

"""
Numba implementation of scipy.ndimage.map_coordinates, numpy.sum and numpy.bincount.
This file is a part of the real time SimNIBS package.

Author: Kristoffer H Madsen
Email: khma@dtu.dk
License: Please contact Axel Thielscher (axthi@dtu.dk)

"""

from numba import prange
from numba import njit
from numpy import empty, zeros, negative, sqrt
import numpy as np

@njit(fastmath=True,parallel=True,nogil=True,cache=True)
def map_coord_lin(x, coords):
    ijk = np.floor(coords).astype(np.int16) #add np.floor for better behavior
    fijk = (coords - ijk)
    n = ijk.shape[0]
    out = empty((n,3), dtype=np.float64)
    bx, by, bz = x.shape[1]-1, x.shape[2]-1, x.shape[3]-1
    # i = np.clip(ijk[:,0],1,bx)
    # j = np.clip(ijk[:,1],1,by)
    # k = np.clip(ijk[:,2],1,bz)
    for l in prange(n):
        i0 = ijk[l,0]
        j0 = ijk[l,1]
        k0 = ijk[l,2]
        i1, j1, k1 = i0+1, j0+1, k0+1
        fi1 = fijk[l, 0]
        fj1 = fijk[l, 1]
        fk1 = fijk[l, 2]

        fi0, fj0, fk0 = 1-fi1, 1-fj1, 1-fk1
        
        if i0<0 or j0<0 or k0<0 or i1>bx or j1>by or k1>bz:
            out[l,0] = out[l,1] = out[l,2] = 0.0
        else:
            out[l, 0] = (
                fi0 * fj0 * fk0 * x[0, i0, j0, k0] +
                fi0 * fj0 * fk1 * x[0, i0, j0, k1] +
                fi0 * fj1 * fk0 * x[0, i0, j1, k0] +
                fi0 * fj1 * fk1 * x[0, i0, j1, k1] +
                fi1 * fj0 * fk0 * x[0, i1, j0, k0] +
                fi1 * fj0 * fk1 * x[0, i1, j0, k1] +
                fi1 * fj1 * fk0 * x[0, i1, j1, k0] +
                fi1 * fj1 * fk1 * x[0, i1, j1, k1]
                )
            out[l, 1] = (
                fi0 * fj0 * fk0 * x[1, i0, j0, k0] +
                fi0 * fj0 * fk1 * x[1, i0, j0, k1] +
                fi0 * fj1 * fk0 * x[1, i0, j1, k0] +
                fi0 * fj1 * fk1 * x[1, i0, j1, k1] +
                fi1 * fj0 * fk0 * x[1, i1, j0, k0] +
                fi1 * fj0 * fk1 * x[1, i1, j0, k1] +
                fi1 * fj1 * fk0 * x[1, i1, j1, k0] +
                fi1 * fj1 * fk1 * x[1, i1, j1, k1]
                )
            out[l, 2] = (
                fi0 * fj0 * fk0 * x[2, i0, j0, k0] +
                fi0 * fj0 * fk1 * x[2, i0, j0, k1] +
                fi0 * fj1 * fk0 * x[2, i0, j1, k0] +
                fi0 * fj1 * fk1 * x[2, i0, j1, k1] +
                fi1 * fj0 * fk0 * x[2, i1, j0, k0] +
                fi1 * fj0 * fk1 * x[2, i1, j0, k1] +
                fi1 * fj1 * fk0 * x[2, i1, j1, k0] +
                fi1 * fj1 * fk1 * x[2, i1, j1, k1]
                )
    return out.T

@njit(fastmath=True,parallel=True,nogil=True,cache=True)
def map_coord_lin_trans(x, coords, M1, t, M2):
    coords = coords.T@M1.T
    coords[:,0] += t[0]
    coords[:,1] += t[1]
    coords[:,2] += t[2]
    ijk = np.floor(coords).astype(np.int16) #can be speeded up by avoiding np.floor
    fijk = (coords - ijk)
    n = ijk.shape[0]
    out = empty((n,3), dtype=np.float64)
    bx, by, bz = x.shape[1]-1, x.shape[2]-1, x.shape[3]-1
    # i = np.clip(ijk[:,0],1,bx)
    # j = np.clip(ijk[:,1],1,by)
    # k = np.clip(ijk[:,2],1,bz)
    for l in prange(n):
        i0 = ijk[l,0]
        j0 = ijk[l,1]
        k0 = ijk[l,2]
        i1, j1, k1 = i0+1, j0+1, k0+1
        fi1 = fijk[l, 0]
        fj1 = fijk[l, 1]
        fk1 = fijk[l, 2]

        fi0, fj0, fk0 = 1-fi1, 1-fj1, 1-fk1
        
        if i0<0 or j0<0 or k0<0 or i1>bx or j1>by or k1>bz:
            out[l,0] = out[l,1]= out[l,2] = 0.0
        else:
            out[l, 0] = (
                fi0 * fj0 * fk0 * x[0, i0, j0, k0] +
                fi0 * fj0 * fk1 * x[0, i0, j0, k1] +
                fi0 * fj1 * fk0 * x[0, i0, j1, k0] +
                fi0 * fj1 * fk1 * x[0, i0, j1, k1] +
                fi1 * fj0 * fk0 * x[0, i1, j0, k0] +
                fi1 * fj0 * fk1 * x[0, i1, j0, k1] +
                fi1 * fj1 * fk0 * x[0, i1, j1, k0] +
                fi1 * fj1 * fk1 * x[0, i1, j1, k1]
                )
            out[l, 1] = (
                fi0 * fj0 * fk0 * x[1, i0, j0, k0] +
                fi0 * fj0 * fk1 * x[1, i0, j0, k1] +
                fi0 * fj1 * fk0 * x[1, i0, j1, k0] +
                fi0 * fj1 * fk1 * x[1, i0, j1, k1] +
                fi1 * fj0 * fk0 * x[1, i1, j0, k0] +
                fi1 * fj0 * fk1 * x[1, i1, j0, k1] +
                fi1 * fj1 * fk0 * x[1, i1, j1, k0] +
                fi1 * fj1 * fk1 * x[1, i1, j1, k1]
                )
            out[l, 2] = (
                fi0 * fj0 * fk0 * x[2, i0, j0, k0] +
                fi0 * fj0 * fk1 * x[2, i0, j0, k1] +
                fi0 * fj1 * fk0 * x[2, i0, j1, k0] +
                fi0 * fj1 * fk1 * x[2, i0, j1, k1] +
                fi1 * fj0 * fk0 * x[2, i1, j0, k0] +
                fi1 * fj0 * fk1 * x[2, i1, j0, k1] +
                fi1 * fj1 * fk0 * x[2, i1, j1, k0] +
                fi1 * fj1 * fk1 * x[2, i1, j1, k1]
                )
    return (M2@out.T)

@njit(parallel=True,fastmath=True,nogil=True,cache=True)
def sumf(x, y, z):
    for j in prange(x.shape[1]):
        for i in prange(x.shape[0]):
            z[i,j] = x[i,j,0]*y[j,0]
            for k in prange(1,x.shape[2]):
                z[i,j] += x[i,j,k]*y[j,k]

@njit(parallel=True,fastmath=True,nogil=True,cache=True)
def sumf2(x, y, w):
    z = empty(x.shape[:2])
    for j in prange(x.shape[1]):
        for i in prange(x.shape[0]):
            z[i,j] = x[i,j,0]*y[j,0]
            for k in prange(1,x.shape[2]):
                z[i,j] += x[i,j,k]*y[j,k]
    return np.bincount(w, z.ravel())

@njit(parallel=True,fastmath=True,nogil=True,cache=True)
def sumf3(v, dadt, nn, g, c):
    sigma_dadt = zeros((v.shape[0], 3))
    elm_node_integral = zeros((v.shape[0], 4))
    for j in prange(v.shape[0]):
        for k in range(3):
            for n in range(3):
                sigma_dadt[j, n] += c[j, n, k] * dadt[j, k]
        for i in range(4):
            for k in range(3):
                elm_node_integral[j,i] -= v[j]*sigma_dadt[j, k]*g[j,i,k]
    return np.bincount(nn, elm_node_integral.ravel()*1e-6)


@njit(parallel=True,fastmath=True,nogil=True,cache=True)
def postp(g, v, dadt, idx1, idx2):#, A, iA, jA, out2):
    out = negative(dadt[idx2])
    for i in prange(g.shape[0]):
        for j in range(g.shape[1]):
            vv = -1000.0 * v[idx1[i,j]]
            for k in range(g.shape[2]):
                out[i,k] += g[i,j,k] * vv
    return out

@njit(parallel=True,fastmath=True,nogil=True,cache=True)
def postp_mag(g, v, dadt, idx1, idx2):#, A, iA, jA, out2):
    out = negative(dadt[idx2])
    out2 = zeros((out.shape[0],1), dtype=out.dtype)
    for i in prange(g.shape[0]):
        for j in range(g.shape[1]):
            vv = -1000.0 * v[idx1[i,j]]
            for k in range(g.shape[2]):
                out[i,k] += g[i,j,k] * vv
        for k in range(g.shape[2]):
            out2[i,0] += out[i,k]**2
        out2[i,0] = sqrt(out2[i,0])
    return out2

# parallel sparse matrix multiplication, note that the array is incremented
# multiple runs will there sum up not reseting the array in-between
@njit(parallel=True,fastmath=True,nogil=True,cache=True)
def spmatmul(A, iA, jA, B, out):
    for i in prange(out.shape[0]):
        for j in range(iA[i],iA[i+1]):
            for k in range(out.shape[1]):
                out[i,k] += A[j]*B[jA[j],k]

# goes from node to elements with flattened index array             
@njit(parallel=True,fastmath=True,nogil=True,cache=True)
def node2elmf(x, nodecoords):
    N = nodecoords.shape[0]//4
    out = empty((N, x.shape[0]), dtype=x.dtype)
    for i in prange(N):
        j = 4*i
        for k in range(out.shape[1]):
            out[i,k] = x[k,nodecoords[j]]
            out[i,k] += x[k,nodecoords[j+1]]
            out[i,k] += x[k,nodecoords[j+2]]
            out[i,k] += x[k,nodecoords[j+3]]
            out[i,k] *= 0.25
    return out

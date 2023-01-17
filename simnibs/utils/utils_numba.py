# -*- coding: utf-8 -*-

"""
Numba implementation of scipy.ndimage.map_coordinates, numpy.sum and numpy.bincount.
This file is a part of the real time SimNIBS package.

Author: Kristoffer H Madsen
Email: khma@dtu.dk
License: Please contact Axel Thielscher (axthi@dtu.dk)

"""

import numba as nb
import numpy as np

@nb.njit(parallel=True,fastmath=True)
def map_coord_lin(x, p, fout):
    N = x.shape[0]
    M = p.shape[1]
    xout=p[0]
    yout=p[1]
    zout=p[2]
    #if fout is None:
    #    fout = np.empty((N,M), dtype=np.float64)
    for mi in nb.prange(M):
        xx = xout[mi]
        yy = yout[mi]
        zz = zout[mi]
        if xx<0 or yy<0 or zz<0 or xx>x.shape[1]-1 or yy>x.shape[2]-1 or zz>x.shape[3]-1:
            for n in range(N):
                fout[n,mi] = 0.0
        else:
            ix = int(xx)
            iy = int(yy)
            iz = int(zz)
            ratx = xx - (ix+0.5)
            raty = yy - (iy+0.5)
            ratz = zz - (iz+0.5)
            asx = np.empty(2)
            asy = np.empty(2)
            asz = np.empty(2)
            asx[0] = 0.5 - ratx
            asx[1] = 0.5 + ratx
            asy[0] = 0.5 - raty
            asy[1] = 0.5 + raty
            asz[0] = 0.5 - ratz
            asz[1] = 0.5 + ratz
            for n in range(N):
                fout[n,mi] = 0.0
                for i in range(2):
                    ixi =  ix + i
                    for j in range(2):
                        iyj = iy + j
                        for k in range(2):
                            izk = iz + k
                            fout[n,mi] += x[n,ixi,iyj,izk]*asx[i]*asy[j]*asz[k]

@nb.njit(parallel=True,fastmath=True)
def map_coord_nn(x, p, fout):
    N = x.shape[0]
    M = p.shape[1]
    #pp = np.round(p).astype(np.int)
    xout=p[0]
    yout=p[1]
    zout=p[2]
    #if fout is None:
    #    fout = np.empty((N,M), dtype=np.float64)
    for mi in nb.prange(M):
        xx = xout[mi]
        yy = yout[mi]
        zz = zout[mi]
        if xx<0 or yy<0 or zz<0 or xx>x.shape[1]-1 or yy>x.shape[2]-1 or zz>x.shape[3]-1:
            for n in range(N):
                fout[n,mi] = 0.0
        else:
            for n in range(N):
                fout[n,mi] = x[n,round(xx),round(yy),round(zz)]


@nb.njit(parallel=True,fastmath=True)
def sumf(x, y, z):
    for i in nb.prange(x.shape[0]):
        for j in nb.prange(x.shape[1]):
            z[i,j] = x[i,j,0]*y[j,0]
            for k in nb.prange(1,x.shape[2]):
                z[i,j] += x[i,j,k]*y[j,k]

@nb.njit(parallel=True,fastmath=True)
def bincount_nb(x, w):
    # the x and w arrays should be in 'C' order.
    # Assumes indexing from 0 and no holes.
    o = np.zeros(int(np.max(x)+1),dtype=np.float64)
    for i in nb.prange(x.shape[0]):
        o[x[i]] += w[i]
    return o

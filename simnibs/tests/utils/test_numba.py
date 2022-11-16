import numpy as np
import pytest
import numba
import time

@numba.njit
def bincount_numba(x,w,N):
    out = np.zeros(N,dtype='float64')
    for i in range(x.shape[0]):
        out[x[i]] += w[i]
    return out

def bincount_python(x,w,N):
    out = np.zeros(N,dtype='float64')
    for i in range(x.shape[0]):
        out[x[i]] += w[i]
    return out

@numba.njit(parallel=True, fastmath=True)
def sum_numba(y):
    s = 0
    for i in numba.prange(y.shape[0]):
        s += y[i]
    return s

class TestNumba:
    def test_numba(self):
        #just a test to see if numba will work after being packaged 
        x = np.random.randint(0,5000,size=5000000)
        w = np.random.randn(x.shape[0])
        t0=time.time()
        bc = bincount_python(x, w, 5000)
        print('time python: %.0fms'%(1000.*(time.time()-t0)))
        t0=time.time()
        bc_np = np.bincount(x, w, 5000)
        print('time numpy: %.0fms'%(1000.*(time.time()-t0)))
        t0=time.time()
        bc_nb = bincount_numba(x, w, 5000)
        print('time numba: %.0fms'%(1000.*(time.time()-t0)))
        assert np.allclose(bc_np, bc)
        assert np.allclose(bc_np, bc_nb)
        assert np.allclose(sum_numba(w), w.sum())

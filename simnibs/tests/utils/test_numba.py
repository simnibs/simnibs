import numpy as np
import pytest
import numba


@numba.njit
def bincount_numba(x,w,N):
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
        x = np.random.randint(0,5000,size=1000000)
        w = np.random.randn(x.shape[0])
        bc_np = np.bincount(x, w, 5000)
        bc_nb = bincount_numba(x, w, 5000)
        assert np.allclose(bc_np, bc_nb)
        assert np.allclose(sum_numba(w), w.sum())

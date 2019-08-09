''' Interface to the MKL PARDISO solver

This code is a modified version of the PyPardiso project

https://github.com/haasad/PyPardisoProject
Copyright (c) 2016, Adrian Haas and ETH Zürich
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this 
list of conditions and the following disclaimer. Redistributions in binary 
form must reproduce the above copyright notice, this list of conditions and the 
following disclaimer in the documentation and/or other materials provided 
with the distribution.
Neither the name of ETH Zürich nor the names of its contributors may be used 
to endorse or promote products derived from this software without specific 
prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Modifications done by Guilherme Saturnino, 2019
'''

import sys
import ctypes
import warnings
import time
import numpy as np
import scipy.sparse as sp
from scipy.sparse import SparseEfficiencyWarning

from simnibs.utils.simnibs_logger import logger


def test_mkl():
    try:
        get_libmkl()
    except OSError:
        return False
    else:
        return True


def get_libmkl():
    if sys.platform == 'darwin':
        return ctypes.CDLL('libmkl_rt.dylib')
    elif sys.platform == 'win32':
        return ctypes.CDLL('mkl_rt.dll')
    else:
        return ctypes.CDLL('libmkl_rt.so')


class Solver:
    """
    Python interface to the Intel MKL PARDISO library for solving large sparse linear systems of equations Ax=b.

    Pardiso documentation: https://software.intel.com/en-us/node/470282

    Parameters
    ------------
    A: scipy.sparse csr or csc
        Sparse square left-hand side matrix
    mtype (optional): int
        Type of matrix. Please see
        https://software.intel.com/en-us/mkl-developer-reference-fortran-pardiso
    """
    # I get segfaults if mtype=2, so usng mtype=11
    def __init__(self, A, mtype=11):
        self._libmkl = get_libmkl()

        self._mkl_pardiso = self._libmkl.pardiso

        # determine 32bit or 64bit architecture
        if ctypes.sizeof(ctypes.c_void_p) == 8:
            self._pt_type = (ctypes.c_int64, np.int64)
        else:
            self._pt_type = (ctypes.c_int32, np.int32)

        self._mkl_pardiso.argtypes = [
            ctypes.POINTER(self._pt_type[0]),    # pt
            ctypes.POINTER(ctypes.c_int32),      # maxfct
            ctypes.POINTER(ctypes.c_int32),      # mnum
            ctypes.POINTER(ctypes.c_int32),      # mtype
            ctypes.POINTER(ctypes.c_int32),      # phase
            ctypes.POINTER(ctypes.c_int32),      # n
            ctypes.POINTER(None),                # a
            ctypes.POINTER(ctypes.c_int32),      # ia
            ctypes.POINTER(ctypes.c_int32),      # ja
            ctypes.POINTER(ctypes.c_int32),      # perm
            ctypes.POINTER(ctypes.c_int32),      # nrhs
            ctypes.POINTER(ctypes.c_int32),      # iparm
            ctypes.POINTER(ctypes.c_int32),      # msglvl
            ctypes.POINTER(None),                # b
            ctypes.POINTER(None),                # x
            ctypes.POINTER(ctypes.c_int32)]      # error

        self._mkl_pardiso.restype = None

        self._pt = np.zeros(64, dtype=self._pt_type[1])
        self._iparm = np.zeros(64, dtype=np.int32)
        self._perm = np.zeros(0, dtype=np.int32)

        self._mtype = mtype
        self._msglvl = False

        self._solve_transposed = False
        self._factorize(A)


    def _factorize(self, A):
        """
        Factorize the matrix A, the factorization will automatically be used if the same
        matrix A is passed to the solve method.

        Parameters
        -------------
        A: scipy.sparse csr or csc
            Spase square matrix
        """

        self._check_A(A)
        self._A = A.copy()
        logger.info('Factorizing FEM matrix unsing MKL PARDISO')
        start = time.time()
        b = np.zeros((A.shape[0],1))
        self._call_pardiso(b, 12)     
        logger.info(f'{time.time()-start:.2f} seconds to factorize matrix')

    def solve(self, b):
        """ solve Ax=b for x

        Parameters
        ----------
        b: numpy ndarray
           right-hand side(s), b.shape[0] needs to be the same as A.shape[0]

        Returns
        --------
        x: numpy ndarray
           solution of the system of linear equations, same shape as input b
        """

        logger.info('Solving system using MKL PARDISO')
        start = time.time()
        b = self._check_b(b)
        x = self._call_pardiso(b, 33)
        logger.info(f'{time.time()-start:.2f} seconds to solve system')
        return x


    def _check_A(self, A):
        if A.shape[0] != A.shape[1]:
            raise ValueError('Matrix A needs to be square, but has shape: {}'.format(A.shape))

        if sp.isspmatrix_csr(A):
            self._solve_transposed = False
            self._iparm[11] = 0
        elif sp.isspmatrix_csc(A):
            self._solve_transposed = True
            self._iparm[11] = 1
        else:
            msg = 'Pardiso requires matrix A to be in CSR or CSC format,' \
                  ' but matrix A is: {}'.format(type(A))
            raise TypeError(msg)

        # scipy allows unsorted csr-indices, which lead to completely wrong pardiso results
        if not A.has_sorted_indices:
            A.sort_indices()

        # scipy allows csr matrices with empty rows. a square matrix with an empty row is singular. calling 
        # pardiso with a matrix A that contains empty rows leads to a segfault, same applies for csc with 
        # empty columns
        if not np.diff(A.indptr).all():
            row_col = 'column' if self._solve_transposed else 'row'
            raise ValueError('Matrix A is singular, because it contains empty'
                             ' {}(s)'.format(row_col))

        if A.dtype != np.float64:
            raise TypeError('Pardiso currently only supports float64, '
                            'but matrix A has dtype: {}'.format(A.dtype))

    def _check_b(self, b):
        if sp.isspmatrix(b):
            warnings.warn('Pardiso requires the right-hand side b'
                          'to be a dense array for maximum efficiency',
                          SparseEfficiencyWarning)
            b = b.todense()

        # pardiso expects fortran (column-major) order if b is a matrix
        if b.ndim == 2:
            b = np.asfortranarray(b)

        if b.shape[0] != self._A.shape[0]:
            raise ValueError("Dimension mismatch: Matrix A {} and array b "
                             "{}".format(self._A.shape, b.shape))

        if b.dtype != np.float64:
            if b.dtype in [np.float16, np.float32, np.int16, np.int32, np.int64]:
                warnings.warn("Array b's data type was converted from "
                              "{} to float64".format(str(b.dtype)), 
                              PardisoWarning)
                b = b.astype(np.float64)
            else:
                raise TypeError('Dtype {} for array b is '
                                'not supported'.format(str(b.dtype)))
        
        return b
        
    def _call_pardiso(self, b, phase):
        x = np.zeros_like(b)
        pardiso_error = ctypes.c_int32(0)
        c_int32_p = ctypes.POINTER(ctypes.c_int32)
        c_float64_p = ctypes.POINTER(ctypes.c_double)

        # 1-based indexing
        ia = self._A.indptr + 1
        ja = self._A.indices + 1

        self._mkl_pardiso(
            self._pt.ctypes.data_as(ctypes.POINTER(self._pt_type[0])), # pt
            ctypes.byref(ctypes.c_int32(1)), # maxfct
            ctypes.byref(ctypes.c_int32(1)), # mnum
            ctypes.byref(ctypes.c_int32(self._mtype)), # mtype -> 11 for real-nonsymetric
            ctypes.byref(ctypes.c_int32(phase)), # phase
            ctypes.byref(ctypes.c_int32(self._A.shape[0])), #N -> number of equations/size of matrix
            self._A.data.ctypes.data_as(c_float64_p), # A -> non-zero entries in matrix
            ia.ctypes.data_as(c_int32_p), # ia -> csr-indptr
            ja.ctypes.data_as(c_int32_p), # ja -> csr-indices
            self._perm.ctypes.data_as(c_int32_p), # perm -> empty
            ctypes.byref(ctypes.c_int32(1 if b.ndim == 1 else b.shape[1])), # nrhs
            self._iparm.ctypes.data_as(c_int32_p), # iparm-array
            ctypes.byref(ctypes.c_int32(self._msglvl)), # msg-level -> 1: statistical info is printed
            b.ctypes.data_as(c_float64_p), # b -> right-hand side vector/matrix
            x.ctypes.data_as(c_float64_p), # x -> output
            ctypes.byref(pardiso_error) # pardiso error
        )
        if pardiso_error.value != 0:
            raise PardisoError(pardiso_error.value)
        else:
            return np.ascontiguousarray(x) # change memory-layout back from fortran to c order

    def __del__(self):
        self._A = sp.csr_matrix((0,0))
        b = np.zeros(0)
        self._call_pardiso(b, 0)



class PardisoWarning(UserWarning):
    pass


class PardisoError(Exception):
    def __init__(self, value):
        self.value = value
    
    def __str__(self):
        return ('The Pardiso solver failed with error code {}. '
                'See Pardiso documentation for details.'.format(self.value))

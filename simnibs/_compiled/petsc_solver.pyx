# cython: language_level=2

import time
import cython
import numpy as np
cimport numpy as np
import tempfile
import os
from ..utils.simnibs_logger import logger

cdef extern from "stdlib.h":
   void free(void* ptr)
   void* malloc(size_t size)

cdef extern from "stdio.h" nogil:
    ctypedef struct FILE
    FILE *fopen(
        const char *filename, const char *opentype)
    int fclose(FILE *stream)

cdef extern from "petscksp.h" nogil:
    ctypedef int PetscInt
    ctypedef float PetscScalar
    ctypedef int PetscErrorCode
    struct _p_KSP
    ctypedef _p_KSP* PetscKSP "KSP"



cdef extern from "_solver.c" nogil:
    PetscErrorCode _petsc_prepare_ksp(
        int argc,char **args, PetscInt N,
        PetscInt row_indices[], PetscInt column_indices[],
        PetscScalar matrix_values[], FILE *stream, PetscKSP *ksp)
    PetscErrorCode _petsc_solve_with_ksp(
        PetscKSP ksp, PetscInt N, PetscScalar rhs[],
        FILE *stream, PetscScalar solution[])
    PetscErrorCode _print_ksp_info(PetscKSP ksp, FILE *stream)
    PetscErrorCode _dealloc(PetscKSP KSP)
    PetscErrorCode _petsc_initialize()
    PetscErrorCode _petsc_finalize()


def petsc_initialize():
    err = _petsc_initialize()
    if err:
        raise SolverError(
            'There was an error initializing PETSc.\n'
            'PETSc returned error code: {}'.format(err))



def petsc_finalize():
    err = _petsc_finalize()
    if err:
        raise SolverError(
            'There was an error finalizing PETSc.\n'
            'PETSc returned error code: {}'.format(err))


class SolverError(Exception):
    pass


cdef class Solver:
    ''' Solver for FEM equations 
    Parameters
    --------------
    petsc_options: str
        Options to be passed on to PETSc, with information like preconditiorer type
    A: scipy.sparse.csr
        csr sparse matrix
    log_level: int (optional)
        Logger level. Default: 20 (Info)
    Attributes
    ------------
    ksp: PETSc KSP object
        Petsc solver object
    N: int
        number of DOF in the system
    log_level: int
        Logger level
    Notes
    ----------
    Having 2 solver objects existting within the same process with different PETSc options might lead to errors!
    '''

    cdef PetscKSP ksp
    cdef int log_level
    cdef int _N
    cdef str _log_file
    cdef FILE * _log_stream

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __cinit__(self, petsc_options, A, log_level=20):
        ''' Set-up
        '''
        petsc_options = petsc_options.split(' ')
        petsc_options = ['dummy'] + petsc_options
        A.sort_indices()
        self._set_log_stream()
        self.log_level = log_level

        cdef PetscKSP new_ksp
        cdef int argc = len(petsc_options)
        cdef char **args  = <char **> malloc(argc * sizeof(char *))
        cdef int N = A.shape[0]
        cdef PetscErrorCode err
        cdef np.ndarray[PetscInt, ndim=1] indptr = A.indptr
        cdef np.ndarray[PetscInt, ndim=1] indices = A.indices
        cdef np.ndarray[PetscScalar, ndim=1] data = A.data
        logger.log(self.log_level, 'Preparing the KSP')
        # Prepare the arguments
        petsc_options = [o.encode() for o in petsc_options]
        for i, o in enumerate(petsc_options):
            args[i] = o
        start = time.time()
        with nogil:
            err = _petsc_prepare_ksp(argc, args, N, &indptr[0], &indices[0],
                                     &data[0], self._log_stream, &new_ksp)

        if err:
            self.log_level = 50
            raise SolverError('There was an error setting up the solver.\n'
                              'PETSc returned error code: {}'.format(err))
        else:
            self._log_record()
            self.ksp = new_ksp
            self._N = N
            end = time.time()
            logger.log(self.log_level, 'Time to prepare the KSP: {0:.2f}s'.format(end-start))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def solve(self, b):
        ''' Solve a series of equations
        Parameters
        -----------
        b: (N x n_sims) np.ndarray 
            Numpy array with "n_sim" right hand sides

        Returns
        ------------
        x: (N x n_sims) np.ndarray
            Solutions
        '''
        if b.ndim == 1:
            b = b[:, None]
        cdef int n_sims = b.shape[1]
        cdef int i
        cdef PetscErrorCode err
        cdef np.ndarray[PetscScalar, ndim=1] rhs = b.T.reshape(-1) # get the right-hand side in C order
        cdef np.ndarray[PetscScalar, ndim=1] solution = np.zeros(b.shape, dtype=float).reshape(-1)
        for i in range(n_sims):
            # Print KSP info
            logger.log(self.log_level, 'Solving system {0} of {1}'.format(i+1, n_sims))
            err = _print_ksp_info(self.ksp, self._log_stream)
            start = time.time()
            if err:
                self.log_level = 50
                raise SolverError('There was an error during the solve.\n'
                                  'PETSc returned error code: {}'.format(err))
            self._log_record()
            # Solve
            with nogil:
                err = _petsc_solve_with_ksp(
                    self.ksp, self._N, &rhs[self._N*i], self._log_stream, &solution[self._N*i])

            if err:
                self.log_level = 50
                raise SolverError('There was an error during the solve.\n'
                                  'PETSc returned error code: {}'.format(err))
            self._log_record()
            end = time.time()
            logger.log(self.log_level, 'Time to solve system: {0:.2f}s'.format(end-start))

        return solution.reshape(n_sims, -1).T

    def _set_log_stream(self):
        with tempfile.NamedTemporaryFile(suffix='.log', delete=False) as tmpfile:
            self._log_file = tmpfile.name
        self._log_stream = fopen(self._log_file.encode(), "w")

    def __dealloc__(self):
        try:
            os.remove(self._log_file)
        except:
            pass
        _dealloc(self.ksp)

    def _log_record(self):
        if self._log_stream:
            fclose(self._log_stream)
        if self._log_file:
            try:
                with open(self._log_file, 'rb') as f:
                    for line in f:
                        line = line.decode().strip()
                        if len(line) > 0:
                            logger.log(self.log_level, line)
                try:
                    os.remove(self._log_file)
                except:
                    pass
            except IOError:
                logger.log(self.log_level, 'Could not log PETSc output')
            except OSError, ValueError:
                pass
        self._set_log_stream()

def petsc_solve(petsc_options, A, b):
    ''' Solves equations of the form Ax = b using PETSc

    Parameters
    --------------
    petsc_options: str
        Options to be passed on to PETSc, with information like preconditiorer type
    A: scipy.sparse.csr
        csr sparse matrix
    b: (N x n_sims) np.ndarray 
        Numpy array with "n_sim" right hand sides

    Returns
    ------------
    x: (N x n_sims) np.ndarray
        Solutions
    '''
    solver = Solver(petsc_options, A)
    return solver.solve(b)

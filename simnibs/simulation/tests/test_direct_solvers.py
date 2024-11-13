import sys
import scipy.sparse as sp
import numpy as np
import pytest

from simnibs.simulation.fem import KSPSolver, MUMPS_Solver, pardiso

def create_matrix(dim, alpha=0.95, smallest_coef=0.1, largest_coef=.9):
    ''' Based o scikit-learn make_sparse_spd_matrix'''
    rng = np.random.default_rng(seed=42)
    chol = -np.eye(dim)
    aux = rng.random((dim, dim))
    aux[aux < alpha] = 0
    aux[aux > alpha] = (smallest_coef
                        + (largest_coef - smallest_coef)
                        * rng.random(np.sum(aux > alpha)))
    aux = np.tril(aux, k=-1)

    # Permute the lines: we don't want to have asymmetries in the final
    # SPD matrix
    permutation = rng.permutation(dim)
    aux = aux[permutation].T[permutation]
    chol += aux
    A = sp.csc_matrix(np.dot(chol.T, chol))

    x = rng.random(dim)
    b = A.dot(x)

    return A,b,x


def test_create_matrix():
    A, b, x = create_matrix(100, .95)

    assert np.allclose(A.dot(x), b)

    Ad = A.todense()
    assert np.all(np.linalg.eigvals(Ad) > 0)


@pytest.mark.parametrize(
    ["ksp_type", "pc_type","factor_solver_type"],
    [
        pytest.param(
            "preonly",
            "lu",
            "mkl_pardiso",
            marks=pytest.mark.skipif(
                sys.platform == "darwin",
                reason="PETSc is not built with Intel MKL on macos."
            ),
        ),
        pytest.param(
            "preonly",
            "cholesky",
            "mkl_pardiso",
            marks=pytest.mark.skipif(
                sys.platform == "darwin",
                reason="PETSc is not built with Intel MKL on macos."
            ),
        ),
        pytest.param(
            "preonly",
            "lu",
            "mumps",
            marks=pytest.mark.skipif(
                sys.platform != "darwin",
                reason="PETSc only built with MUMPS on macos."
            ),
        ),
        pytest.param(
            "preonly",
            "cholesky",
            "mumps",
            marks=pytest.mark.skipif(
                sys.platform != "darwin",
                reason="PETSc only built with MUMPS on macos."
            ),
        ),
    ]
)
class TestPETScFactorSolver:
    def test_solve(self, ksp_type, pc_type, factor_solver_type):
        A, b, x = create_matrix(1000, .99)
        solver = KSPSolver(A, ksp_type, pc_type, factor_solver_type)
        x_pd = solver.solve(b)
        np.testing.assert_allclose(x, x_pd, atol=1e-12)

    def test_solve_csr(self, ksp_type, pc_type, factor_solver_type):
        A, b, x = create_matrix(1000, .99)
        A = A.tocsr()
        solver = KSPSolver(A, ksp_type, pc_type, factor_solver_type)
        x_pd = solver.solve(b)
        np.testing.assert_allclose(x, x_pd, atol=1e-12)

    def test_love_many_rhs(self, ksp_type, pc_type, factor_solver_type):
        A, _, _ = create_matrix(1000, .99)
        x = np.random.random((1000, 3))
        b = A.dot(x)
        solver = KSPSolver(A, ksp_type, pc_type, factor_solver_type)
        x_pd = solver.solve(b)
        np.testing.assert_allclose(x, x_pd, atol=1e-12)


class TestPythonMUMPS:
    def test_solve(self):
        A, b, x = create_matrix(1000, .99)
        solver = MUMPS_Solver(A)
        x_pd = solver.solve(b)
        np.testing.assert_allclose(x, x_pd, atol=1e-12)

    def test_solve_csr(self):
        A, b, x = create_matrix(1000, .99)
        A = A.tocsr()
        solver = MUMPS_Solver(A)
        x_pd = solver.solve(b)
        np.testing.assert_allclose(x, x_pd, atol=1e-12)

    def test_love_many_rhs(self):
        A, _, _ = create_matrix(1000, .99)
        x = np.random.random((1000, 3))
        b = A.dot(x)
        solver = MUMPS_Solver(A)
        x_pd = solver.solve(b)
        np.testing.assert_allclose(x, x_pd, atol=1e-12)

class TestPythonPardiso:
    def test_solve(self):
        A, b, x = create_matrix(1000, .99)
        solver = pardiso.Solver(A)
        x_pd = solver.solve(b)
        np.testing.assert_allclose(x, x_pd, atol=1e-12)

    def test_solve_csr(self):
        A, b, x = create_matrix(1000, .99)
        A = A.tocsr()
        solver = pardiso.Solver(A)
        x_pd = solver.solve(b)
        np.testing.assert_allclose(x, x_pd, atol=1e-12)

    def test_love_many_rhs(self):
        A, _, _ = create_matrix(1000, .99)
        x = np.random.random((1000, 3))
        b = A.dot(x)
        solver = pardiso.Solver(A)
        x_pd = solver.solve(b)
        np.testing.assert_allclose(x, x_pd, atol=1e-12)
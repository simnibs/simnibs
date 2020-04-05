import scipy.sparse as sp
import numpy as np

from .. import pardiso

def create_matrix(dim, alpha=0.95, smallest_coef=0.1, largest_coef=.9):
    ''' Based o scikit-learn make_sparse_spd_matrix'''
    chol = -np.eye(dim)
    aux = np.random.rand(dim, dim)
    aux[aux < alpha] = 0
    aux[aux > alpha] = (smallest_coef
                        + (largest_coef - smallest_coef)
                        * np.random.rand(np.sum(aux > alpha)))
    aux = np.tril(aux, k=-1)

    # Permute the lines: we don't want to have asymmetries in the final
    # SPD matrix
    permutation = np.random.permutation(dim)
    aux = aux[permutation].T[permutation]
    chol += aux
    A = sp.csc_matrix(np.dot(chol.T, chol))

    x = np.random.rand(dim)
    b = A.dot(x)

    return A,b,x



def test_create_matrix():
    A, b, x = create_matrix(100, .95)

    assert np.allclose(A.dot(x), b)

    Ad = A.todense()
    assert np.all(np.linalg.eigvals(Ad) > 0)


class TestPardiso():
    def test_solve(self):
        A, b, x = create_matrix(1000, .99)
        solver = pardiso.Solver(A)
        x_pd = solver.solve(b)
        assert np.allclose(x, x_pd)

    def test_solve_csr(self):
        A, b, x = create_matrix(1000, .99)
        A = A.tocsr()
        solver = pardiso.Solver(A)
        x_pd = solver.solve(b)
        assert np.allclose(x, x_pd)

    def test_love_many_rhs(self):
        A, _, _ = create_matrix(1000, .99)
        x = np.random.random((1000, 3))
        b = A.dot(x)
        solver = pardiso.Solver(A)
        x_pd = solver.solve(b)
        assert np.allclose(x, x_pd)

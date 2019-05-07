import pytest
import numpy as np
import functools
from simnibs.pygpc import adaptive, grid

@pytest.fixture
def uniform_dist():
    pdftype = ['beta', 'beta']
    gridshape = [[1, 1], [1, 1]]
    limits = [[-1, -1], [1, 1]]
    return pdftype, gridshape, limits

class TestExpandPolynomial:
    def test_expand_poly(self):
        active_set = [(1, 0), (0, 1)]
        old_set = [(0, 0)]
        to_expand = (1, 0)
        order_max = 10
        interaction_max = 10
        active_set, old_set, expand = adaptive._expand_polynomial(
            active_set, old_set, to_expand, order_max, interaction_max)
        assert active_set == [(0, 1), (2, 0)]
        assert old_set == [(0, 0), (1, 0)]
        assert expand == [(2, 0)]

        to_expand = (0, 1)
        active_set, old_set, expand = adaptive._expand_polynomial(
            active_set, old_set, to_expand, order_max, interaction_max)
        assert active_set == [(2, 0), (1, 1), (0, 2)]
        assert old_set == [(0, 0), (1, 0), (0, 1)]
        assert expand == [(1, 1), (0, 2)]



class TestRelativeError:
    def test_2d(self):
        a = np.ones((10, 2))
        b = 2 * a
        assert np.isclose(adaptive._relative_error(a, b), 1)
    def test_1d(self):
        a = np.ones(2)
        b = 2 * a
        assert np.isclose(adaptive._relative_error(a, b), 1)

class TestTikhonov:
    def test_1d_zero_alpha(self):
        A = np.random.rand(100, 20)
        x = np.random.rand(20)
        b = A.dot(x)
        alpha = 0.
        x_solve = adaptive._tikhonov(A, b, alpha)
        assert np.allclose(x, x_solve)

    def test_1d_alpha(self):
        A = np.random.rand(100, 20)
        b = np.random.rand(100)
        alpha = 5.
        x = adaptive._tikhonov(A, b, alpha)
        ''' Test solution using the KKT conditions '''
        kkt = A.T.dot(A).dot(x) - A.T.dot(b) + alpha*x
        ''' Figure out how to write this test '''
        assert np.allclose(kkt, 0)

    def test_2d_zero_alpha(self):
        A = np.random.rand(100, 20)
        x = np.random.rand(20, 2)
        b = A.dot(x)
        alpha = 0.
        x_solve = adaptive._tikhonov(A, b, alpha)
        assert np.allclose(x, x_solve)

    def test_2d_alpha(self):
        A = np.random.rand(100, 20)
        b = np.random.rand(100, 2)
        alpha = 5.
        x = adaptive._tikhonov(A, b, alpha)
        ''' Test solution using the KKT conditions '''
        kkt = A.T.dot(A).dot(x) - A.T.dot(b) + alpha*x
        ''' Figure out how to write this test '''
        assert np.allclose(kkt, 0)

class TestCV:
    def test_k_fold_cv_regression(self):
        np.random.seed(1)
        A = np.eye(50)
        x = np.ones(50)
        r = functools.partial(
            adaptive._tikhonov,
            alpha=0.)
        b = A.dot(x)
        noise = .5 * np.ones(50)
        b_noise = b + noise
        regression = lambda A, b: x
        error_eq = lambda data, model: np.linalg.norm(data - model)/data.shape[0]
        eps = adaptive._k_fold_cv_regression(
            A, b_noise,
            regression,
            error_eq=error_eq,
            k=10)
        assert np.isclose(eps, .5)


class TestRegularizedRegression:
    def test_expand(self, uniform_dist):
        pdfshape, gridshape, limits = uniform_dist
        # Define grid
        g = grid.randomgrid(pdfshape, gridshape, limits, 50)
        # Shifted Legendre Polynomials
        data = np.vstack(
            [2 * g.coords[:, 0],
             3 * g.coords[:, 1],
             2 * g.coords[:, 0] * g.coords[:, 1]]).T
        reg = adaptive.RegularizedRegression(
            pdfshape, gridshape, limits, [1, 1], 10, 2, g)
        coeffs, err = reg.expand(data)
        d_approx = reg.evaluate(coeffs, g.coords_norm)
        assert np.allclose(
            np.linalg.norm(data - d_approx, axis=1) / np.linalg.norm(data, axis=1),
            0, atol=1e-3)
        assert np.allclose(err, 0, atol=1e-3)

    def test_add_sampling_points(self, uniform_dist):
        pdfshape, gridshape, limits = uniform_dist
        # Define grid
        g = grid.randomgrid(pdfshape, gridshape, limits, 1)
        reg = adaptive.RegularizedRegression(
            pdfshape, gridshape, limits, [1, 1], 10, 1, g)
        new_points = reg.add_n_sampling_points(50)
        assert new_points.shape == (50, 2)
        assert reg.construct_gpc_matrix().shape == (51, 3)


class TestAdaptive:
    def test_adaptive_expand(self,  uniform_dist):
        np.random.seed(1)
        pdftype, pdfshape, limits = uniform_dist
        function = lambda x: np.array((3 * x[0], 1 * x[1], 2 * x[0]**2))
        reg, _ = adaptive.run_reg_adaptive_grid(
                pdftype, pdfshape,
                limits, function,
                data_poly_ratio=2)
        assert reg.construct_gpc_matrix().shape == (8,4)


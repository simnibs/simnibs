import itertools
from unittest import mock
import numpy as np
import scipy.optimize
import pytest
import warnings

from .. import optimization_methods
from ...simulation.analytical_solutions import fibonacci_sphere

@pytest.fixture()
def optimization_variables():
    l = np.array([4, 1, 1, -5], dtype=float)
    A = np.array([[1, 2, 3, 4],
                  [3, 1, 2, 2],
                  [2, 3, 1, 1],
                  [1, 3, 5, 1]], dtype=float)
    Q = A.T.dot(A)
    P = np.array([[-1, -1, -1, -1],
                  [1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])
    return l, Q, P

def objective(l, Q, lam, x):
    P = np.vstack([-np.ones(len(l)), np.eye(len(l))])
    return l.dot(x) + 0.5 * x.T.dot(Q.dot(x)) + lam * np.linalg.norm(P.dot(x), 1)

@pytest.fixture()
def optimization_variables_avg():
    l = np.array([4, 1, 1, -5, 3], dtype=float)
    A = np.array([[1, 2, 3, 4],
                  [3, 1, 2, 2],
                  [2, 3, 1, 1],
                  [1, 3, 5, 1]], dtype=float)
    A = np.vstack([-np.mean(A, axis=1), A - np.mean(A, axis=1)])
    Q = A.dot(A.T)
    return l, Q, np.ones(5)

@pytest.fixture()
def optimization_variables_avg_QCQP():
    l = np.array([4, 1, 1, -5, 3], dtype=float)
    A = np.array([[1, 2, 3, 4],
                  [3, 1, 2, 2],
                  [2, 3, 1, 1],
                  [1, 3, 5, 1]], dtype=float)
    A = np.vstack([-np.mean(A, axis=1), A - np.mean(A, axis=1)])
    Q = A.dot(A.T)
    A_in = np.array([[1, 5, 2, 4],
                     [5, 2, 1, 2],
                     [0, 1, 1, 1],
                     [2, 4, 3, 5]], dtype=float)
    A_in = np.vstack([-np.mean(A_in, axis=1), A_in - np.mean(A_in, axis=1)])
    Q_in = A_in.dot(A_in.T)
    Q_in += np.outer(l, l)
    return l, Q, np.ones(5), Q_in


def optimize_scipy(l, Q, C, d, A=None, b=None, atol=1e-5):
    objective = lambda x: l.dot(x) + .5 * x.dot(Q.dot(x))
    x0 = np.zeros_like(l)
    jac = lambda x: l + x.dot(Q)
    #hess = lambda x: Q
    iq_constraint = {}
    iq_constraint['type'] = 'ineq'
    iq_constraint['fun'] = lambda x: d - C.dot(x)
    iq_constraint['jac'] = lambda x: -C
    constraints = (iq_constraint, )
    if A is not None:
        eq_constraint = {}
        eq_constraint['type'] = 'eq'
        eq_constraint['fun'] = lambda x: A.dot(x) - b
        eq_constraint['jac'] = lambda x: A
        constraints = (iq_constraint, eq_constraint)

    res = scipy.optimize.minimize(
        objective, x0, jac=jac,
        constraints=constraints)
    return res.x

def optimize_comp(l, A, max_el_current=None, max_total_current=None):
    objective = lambda x: -l.dot(x)
    x0 = np.zeros_like(l)
    jac = lambda x: -l

    eq_constraint = {}
    eq_constraint['type'] = 'eq'
    eq_constraint['fun'] = lambda x: A.dot(x)
    eq_constraint['jac'] = lambda x: A
    constraints = (eq_constraint, )

    if max_el_current is not None:
        bounds = scipy.optimize.Bounds(-max_el_current, max_el_current)
    else:
        bounds = None

    if max_total_current is not None:
        iq_constraint = {}
        iq_constraint['type'] = 'ineq'
        iq_constraint['fun'] = lambda x: 2 * max_total_current - np.linalg.norm(x, 1)
        constraints += (iq_constraint, )

    res = scipy.optimize.minimize(
        objective, x0, jac=jac,
        bounds=bounds,
        constraints=constraints)
    return res.x


def optimize_focality(l, Q, b,
                      max_el_current=None,
                      max_total_current=None,
                      Qin=None,
                      max_angle=None):
    objective = lambda x: .5 * x.dot(Q.dot(x))
    x0 = np.zeros(Q.shape[0])
    jac = lambda x: x.dot(Q)
    #hess = lambda x: Q
    m_constraint = {}
    m_constraint['type'] = 'eq'
    m_constraint['fun'] = lambda x: l.dot(x) - b
    m_constraint['jac'] = lambda x: l

    e_constraint = {}
    e_constraint['type'] = 'eq'
    e_constraint['fun'] = lambda x: np.sum(x)
    e_constraint['jac'] = lambda x: np.ones_like(x)
    constraints = (m_constraint, e_constraint)

    if max_el_current is not None:
        bounds = scipy.optimize.Bounds(-max_el_current, max_el_current)
    else:
        bounds = None
    if max_total_current is not None:
        iq_constraint = {}
        iq_constraint['type'] = 'ineq'
        iq_constraint['fun'] = lambda x: 2 * max_total_current - np.linalg.norm(x, 1)
        constraints += (iq_constraint, )
    if Qin is not None:
        a_constraint = {}
        a_constraint['type'] = 'ineq'
        a_constraint['fun'] = lambda x: (b/np.cos(np.deg2rad(max_angle)))**2 - x.dot(Qin.dot(x))
        constraints += (a_constraint, )

    res = scipy.optimize.minimize(
        objective, x0, jac=jac,
        bounds=bounds,
        constraints=constraints)
    return res.x

def optimize_lstsq(A, b,
                   max_el_current=None,
                   max_total_current=None,
                   atol=1e-5):
    objective = lambda x: np.linalg.norm(A.dot(x) - b)**2
    x0 = np.zeros(A.shape[1])
    #jac = lambda x: b.dot(A) + x.dot(Q)
    jac = None

    e_constraint = {}
    e_constraint['type'] = 'eq'
    e_constraint['fun'] = lambda x: np.sum(x)
    e_constraint['jac'] = lambda x: np.ones_like(x)
    constraints = (e_constraint, )

    if max_el_current is not None:
        bounds = scipy.optimize.Bounds(-max_el_current, max_el_current)
    else:
        bounds = None
    if max_total_current is not None:
        iq_constraint = {}
        iq_constraint['type'] = 'ineq'
        iq_constraint['fun'] = lambda x: 2 * max_total_current - np.linalg.norm(x, 1)
        constraints += (iq_constraint, )

    warnings.filterwarnings("ignore", message="Values in x were outside bounds during a minimize step, clipping to bounds")

    res = scipy.optimize.minimize(
        objective, x0, jac=jac,
        bounds=bounds,
        constraints=constraints)

    return res.x

class TestActiveSetQp:
    def test_inactive_constraints(self, optimization_variables):
        l, Q, P = optimization_variables
        C = np.eye(len(l))
        x0 = np.zeros(len(l))
        d = 1e3 * np.ones(C.shape[0])
        x = optimization_methods._active_set_QP(l, Q, C, d, x0)
        x_sp = optimize_scipy(l, Q, C, d)
        assert np.linalg.norm(C.dot(x) <= d)
        assert np.isclose(objective(l, Q, 0, x), objective(l, Q, 0, x_sp),
                          rtol=1e-4, atol=1e-4)

    def test_active_bound_constraints(self, optimization_variables):
        l, Q, P = optimization_variables
        C = np.vstack([np.eye(len(l)), -np.eye(len(l))])
        x0 = np.zeros(len(l))
        d = 0.3 * np.ones(C.shape[0])
        x = optimization_methods._active_set_QP(l, Q, C, d, x0)

        x_sp = optimize_scipy(l, Q, C, d)

        assert np.linalg.norm(C.dot(x) <= d + 1e-4)
        assert np.isclose(objective(l, Q, 0, x), objective(l, Q, 0, x_sp),
                          rtol=1e-4, atol=1e-4)

    def test_active_other_constraints(self, optimization_variables):
        l, Q, P = optimization_variables
        C = np.vstack([np.eye(len(l)), np.ones(len(l)), -np.eye(len(l)), -np.ones(len(l))])
        x0 = np.zeros(len(l))
        d =  0.3 * np.ones(C.shape[0])
        x = optimization_methods._active_set_QP(l, Q, C, d, x0)

        x_sp = optimize_scipy(l, Q, C, d)
        assert np.linalg.norm(C.dot(x) <= d)
        assert np.isclose(objective(l, Q, 0, x), objective(l, Q, 0, x_sp),
                          rtol=1e-4, atol=1e-4)

    def test_equality_constraints(self, optimization_variables):
        l, Q, P = optimization_variables
        C = np.vstack([np.eye(len(l)), -np.eye(len(l))])
        A = np.ones((1, len(l)))
        b = np.zeros(1)
        x0 = np.zeros(len(l))
        d = 0.3 * np.ones(C.shape[0])
        x = optimization_methods._active_set_QP(l, Q, C, d, x0, A=A, b=b)
        x_sp = optimize_scipy(l, Q, C, d, A, b)
        assert np.isclose(A.dot(x), b)
        assert np.linalg.norm(C.dot(x) <= d + 1e-4)
        assert np.isclose(objective(l, Q, 0, x), objective(l, Q, 0, x_sp),
                          rtol=1e-4, atol=1e-4)

    def test_remove_constraint(self, optimization_variables):
        optimum = np.array([2.5, 3, 1])
        Q = np.array([[8, 0, 0], [0, 4, 0], [0, 0, 1]])
        l = Q.dot(optimum)
        C = np.vstack([np.eye(len(l)), -np.eye(len(l))])
        x0 = 2 * np.ones(len(l))
        d = 2 * np.ones(C.shape[0])
        x = optimization_methods._active_set_QP(l, Q, C, d, x0)

        x_sp = optimize_scipy(l, Q, C, d)
        assert np.linalg.norm(C.dot(x) <= d + 1e-4)
        assert np.isclose(objective(l, Q, 0, x), objective(l, Q, 0, x_sp),
                          rtol=1e-4, atol=1e-4)


class TestTESOptimizationProblem():
    def test_calc_Q(self):
        A = np.random.random((2, 5, 3))
        volumes = np.array([1, 2, 2, 2, 4])
        tes_opt = optimization_methods.TESOptimizationProblem(
            A, 1e3, 1e3, volumes
        )
        currents = np.array([1, 0])
        field = np.vstack([A[..., i].T.dot(currents) for i in range(3)])
        energy = np.sum(np.linalg.norm(field, axis=0)**2 * volumes)
        energy /= np.sum(volumes)
        currents = np.array([-1, 1, 0])
        assert np.allclose(currents.dot(tes_opt.Q.dot(currents)), energy)

    def test_bound_constraints(self):
        A = np.random.random((2, 5, 3))
        tes_opt = optimization_methods.TESOptimizationProblem(
            A, 1e4, 1
        )
        C, d = tes_opt._bound_contraints()
        currents = tes_opt.extend_currents(
            np.array([0.5, 2, -3])
        )
        assert np.all((C.dot(currents) < d)[:6] ==
                       [True, False, True, True, True, False])
        assert np.all((C.dot(currents) <= d)[6:])

    def test_l1_constraint(self):
        A = np.random.random((2, 5, 3))
        tes_opt = optimization_methods.TESOptimizationProblem(
            A, 1, 1e4
        )
        C, d = tes_opt._l1_constraint()
        currents = tes_opt.extend_currents(
            np.array([0.3, -0.9, 0.6])
        )
        assert np.all(C.dot(currents) < d)
        assert np.isclose(C.dot(currents), 1.8)

    def test_l1_constraint_fail(self):
        A = np.random.random((2, 5, 3))
        tes_opt = optimization_methods.TESOptimizationProblem(
            A, 1, 1e4
        )
        C, d = tes_opt._l1_constraint()
        currents = tes_opt.extend_currents(
            np.array([1.2, -1.2, 0])
        )
        assert ~np.all(C.dot(currents) < d)
        assert np.isclose(C.dot(currents), 2.4)

class TestCalcl:
    def test_calc_l(self):
        # test if l.dot(currents) really is the average current
        A = np.random.random((2, 5, 3))
        targets = [0, 1]
        target_direction = [[1, 0, 0], [1, 0, 0]]
        volumes = np.array([1, 2, 2, 2, 2])
        l = optimization_methods._calc_l(A, targets, target_direction, volumes)
        currents = [1, 1]
        field_x = A[..., 0].T.dot(currents)
        avg = np.average(field_x[[0, 1]], weights=[1, 2])
        assert np.allclose(avg, l.dot([-2, 1, 1]))

    def test_calc_different_directions(self):
        A = np.random.random((2, 5, 3))
        targets = [0, 1]
        target_direction = np.array([[1, 0, 0], [0, 1, 0]])
        volumes = np.array([1, 2, 2, 2, 2])
        l = optimization_methods._calc_l(
            A, targets, target_direction, volumes
        )
        currents = [1, 1]
        field_x = A[..., 0].T.dot(currents)
        field_y = A[..., 1].T.dot(currents)
        avg = (field_x[0]*1 + field_y[1]*2)/3
        assert np.allclose(avg, l.dot([-2, 1, 1]))


class TestTESLinearConstrained():
    def test_calc_l(self):
        A = np.random.random((2, 5, 3))
        targets = [0, 1]
        target_direction = [[1, 0, 0], [1, 0, 0]]
        volumes = np.array([1, 2, 2, 2, 2])
        tes_problem = optimization_methods.TESLinearConstrained(
            A, weights=volumes
        )
        tes_problem.add_linear_constraint(targets, target_direction, 0.2)
        currents = [-2, 1, 1]
        field_x = A[..., 0].T.dot(currents[1:])
        avg = np.average(field_x[[0, 1]], weights=[1, 2])
        assert np.allclose(avg, tes_problem.l.dot(currents))
        assert np.allclose(0.2, tes_problem.target_means)

    def test_calc_different_directions(self):
        A = np.random.random((2, 5, 3))
        targets = [0, 1]
        target_direction = np.array([[1, 0, 0], [0, 1, 0]])
        volumes = np.array([1, 2, 2, 2, 2])
        tes_problem = optimization_methods.TESLinearConstrained(
            A, weights=volumes
        )
        tes_problem.add_linear_constraint(targets, target_direction, 0.2)
        currents = [1, 1]
        field_x = A[..., 0].T.dot(currents)
        field_y = A[..., 1].T.dot(currents)
        avg = (field_x[0]*1 + field_y[1]*2)/3
        assert np.allclose(avg, tes_problem.l.dot([-2, 1, 1]))

    def test_calc_many_targets(self):
        A = np.random.random((2, 5, 3))
        volumes = np.array([1, 2, 2, 2, 2])
        tes_problem = optimization_methods.TESLinearConstrained(
            A, weights=volumes
        )
        tes_problem.add_linear_constraint([0], [1, 0, 0], 0.2)
        tes_problem.add_linear_constraint([1], [1, 0, 0], 0.4)
        currents = [-2, 1, 1]
        field_x = A[..., 0].T.dot(currents[1:])
        assert np.allclose(field_x[[0, 1]], tes_problem.l.dot(currents))
        assert np.allclose([0.2, 0.4], tes_problem.target_means)

    @pytest.mark.parametrize('target_mean', [1e-4, 5e-2])
    def test_feasible(self, target_mean):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)
        targets = [0, 1]
        target_direction = [[1, 0, 0], [1, 0, 0]]

        max_total_current = 0.2
        max_el_current = 0.1

        tes_problem = optimization_methods.TESLinearConstrained(
            leadfield, max_total_current, max_el_current
        )
        tes_problem.add_linear_constraint(targets, target_direction, target_mean)

        x = tes_problem.solve()

        x_sp = optimize_focality(
            tes_problem.l, tes_problem.Q,
            target_mean, max_el_current,
            max_total_current
        )

        assert np.all(np.abs(x) <= max_el_current)
        assert np.all(np.linalg.norm(x) <= 2*max_total_current)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(tes_problem.l.dot(x), tes_problem.target_means)
        assert np.allclose(
            x.dot(tes_problem.Q.dot(x)),
            x_sp.dot(tes_problem.Q.dot(x_sp)),
            rtol=1e-4, atol=1e-5
        )

    def test_infeasible(self):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)
        targets = [0, 1]
        target_direction = [[1, 0, 0], [1, 0, 0]]

        target_mean = 1e4
        max_total_current = 0.2
        max_el_current = 0.1

        tes_problem = optimization_methods.TESLinearConstrained(
            leadfield, max_total_current, max_el_current
        )
        tes_problem.add_linear_constraint(targets, target_direction, target_mean)

        x = tes_problem.solve()

        x_sp = optimize_comp(
            np.squeeze(tes_problem.l),
            np.ones((1, 6)),
            max_el_current,
            max_total_current
        )

        assert np.linalg.norm(x, 1) <= 2 * max_total_current
        assert np.all(np.abs(x) <= max_el_current)
        assert np.isclose(np.sum(x), 0)
        assert np.allclose(tes_problem.l.dot(x), tes_problem.l.dot(x_sp),
                           rtol=1e-4, atol=1e-5)

    def test_multi_target_feasible(self):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)

        max_total_current = 0.2
        max_el_current = 0.1

        tes_problem = optimization_methods.TESLinearConstrained(
            leadfield, max_total_current, max_el_current
        )
        tes_problem.add_linear_constraint([0], [1, 0, 0], 1e-3)
        tes_problem.add_linear_constraint([1], [1, 0, 0], 2e-3)

        x = tes_problem.solve()

        assert np.all(np.abs(x) <= max_el_current)
        assert np.all(np.linalg.norm(x) <= 2*max_total_current)
        assert np.isclose(np.sum(x), 0)
        assert np.allclose(tes_problem.l.dot(x), tes_problem.target_means)


class TestCalcQnorm:
    def test_calc_Qnorm(self):
        A = np.random.random((2, 5, 3))
        targets = [0, 1]
        volumes = np.array([1, 2, 2, 2, 4])
        Qin = optimization_methods._calc_Qnorm(
            A, targets, volumes,
        )
        currents = np.array([1, 0])
        field = np.einsum('ijk, i->jk', A, currents)
        avg_in = np.sum(
            np.linalg.norm(field[targets], axis=1)**2 * volumes[targets])
        avg_in /= np.sum(volumes[targets])
        currents = np.array([-1, 1, 0])
        assert np.allclose(currents.dot(Qin.dot(currents)), avg_in)

class TestTESLinearAngleConstrained:
    def test_limit_angle_inactive(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        Q_in = 5e1 * np.eye(len(l))
        max_angle = 45
        max_el_current = 2e-3
        max_total_current = 4e-3
        f = .02
        x = optimization_methods._linear_angle_constrained_tes_opt(
            np.atleast_2d(l), np.atleast_1d(f), Q,
            max_el_current, max_total_current, Q_in, max_angle
        )
        x_sp = optimize_focality(
            np.atleast_2d(l), Q,
            np.atleast_1d(f),
            max_el_current, max_total_current,
            Qin=Q_in, max_angle=max_angle
        )

        assert np.linalg.norm(x, 1) <= 2 * max_total_current + 1e-4
        assert np.all(np.abs(x) <= max_el_current + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), f)
        assert np.arccos(l.dot(x) / np.sqrt(x.dot(Q_in).dot(x))) <= np.deg2rad(max_angle)
        assert np.allclose(x.dot(Q.dot(x)), x_sp.dot(Q.dot(x_sp)), rtol=1e-4, atol=1e-5)

    def test_limit_angle_reduce_target(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        Q_in = 5e1 * np.eye(len(l))
        max_angle = 20
        max_el_current = 2e-3
        max_total_current = 4e-3
        f = .02
        x = optimization_methods._linear_angle_constrained_tes_opt(
            np.atleast_2d(l), np.atleast_1d(f), Q,
            max_el_current, max_total_current, Q_in, max_angle
        )
        assert np.linalg.norm(x, 1) <= 2 * max_total_current + 1e-4
        assert np.all(np.abs(x) <= max_el_current + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.arccos(l.dot(x) / np.sqrt(x.dot(Q_in).dot(x))) <= np.deg2rad(max_angle)

    def test_both_limit_angle_Q_iteration(self, optimization_variables_avg_QCQP):
        l, Q, A, Q_in = optimization_variables_avg_QCQP
        max_angle = 10
        max_el_current = 2e-3
        max_total_current = 4e-3
        f = .01
        x = optimization_methods._linear_angle_constrained_tes_opt(
            np.atleast_2d(l), np.atleast_1d(f), Q,
            max_el_current, max_total_current, Q_in, max_angle
        )

        x_sp = optimize_focality(
            l, Q, f, max_el_current=max_el_current,
            max_total_current=max_total_current,
            Qin=Q_in, max_angle=max_angle
        )
        assert np.linalg.norm(x, 1) <= 2 * max_total_current + 1e-4
        assert np.all(np.abs(x) <= max_el_current + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), f)
        assert np.arccos(l.dot(x) / np.sqrt(x.dot(Q_in).dot(x))) <= np.deg2rad(max_angle)
        assert np.allclose(x.dot(Q.dot(x)), x_sp.dot(Q.dot(x_sp)), rtol=1e-4, atol=1e-5)

    @pytest.mark.parametrize('target_intensity', [1e-2, 1e-1])
    def test_linear_angle_constrained_Q(self, target_intensity):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)

        max_total_current = 0.2
        max_el_current = 0.1
        max_angle = 20

        tes_problem = optimization_methods.TESLinearAngleConstrained(
            [0], [1, 0, 0], target_intensity, max_angle,
            leadfield, max_total_current,
            max_el_current
        )
        x = tes_problem.solve()
        field = np.array([x[1:].dot(leadfield[..., i]) for i in range(3)]).T
        angle = np.rad2deg(np.arccos(field[0, 0]/np.linalg.norm(field[0])))
        assert np.all(np.abs(x) <= max_el_current)
        assert np.all(np.linalg.norm(x) <= 2*max_total_current)
        assert np.isclose(np.sum(x), 0)
        assert angle <= max_angle


class TestBBAlgorithm:
    def test_bb_node(self):
        def bounds_func(s):
            a = np.array([0, 4, 1, 3, 2, 6, 7, 8, 9, 5])
            if len(s.active) > 4:
                return np.inf
            to_consider = s.active + s.unassigned
            ub = np.sum(a[to_consider])
            lb = np.sum(np.sort(a[to_consider])[:4])
            for i in np.argsort(a[to_consider]):
                el = to_consider[i]
                if el in s.unassigned:
                    split1 = s.inactivate(el)
                    split2 = s.activate(el)
                    break
            if len(s.unassigned) == 0:
                split1 = None
                split2 = None
            return ub, lb, split1, split2

        init = optimization_methods.bb_state([], [], list(range(10)))
        bb_node = optimization_methods.bb_node(init, bounds_func)
        assert bb_node.lb_val == 0 + 1 + 2 + 3
        assert bb_node.ub_val == sum(range(10))
        assert np.all(bb_node.child1.active == [])
        assert np.all(bb_node.child1.inactive == [0])
        assert np.all(bb_node.child2.active == [0])
        assert np.all(bb_node.child2.inactive == [])
        c1, c2 = bb_node.split()
        assert c1.lb_val == 1 + 2 + 3 + 4
        assert c1.ub_val == sum(range(1, 10))

        assert c2.lb_val == 0 + 1 + 2 + 3
        assert c2.ub_val == sum(range(0, 10))

    def test_bb_algorithm(self):
        # Function with gives out the bounds as well as the new sets to split
        def bounds_func(s):
            a = np.array([4, 0, 1, 3, 2, 6, 7, 8, 9, 5]) * .1

            if len(s.active) > 4:
                return np.inf, np.inf, None, None

            to_consider = s.active + s.unassigned
            v_in_consideration = a[to_consider]
            ub = np.sum(v_in_consideration)
            lb = np.sum(np.sort(v_in_consideration)[:4])
            for i in np.argsort(v_in_consideration):
                el = to_consider[i]
                if el in s.unassigned:
                    split1 = s.inactivate(el)
                    split2 = s.activate(el)
                    break
            if len(s.unassigned) == 0:
                split1 = None
                split2 = None

            return ub, lb, split1, split2

        init = optimization_methods.bb_state([], [], list(range(10)))
        eps = 1e-1
        final_state = optimization_methods._branch_and_bound(init, bounds_func, eps, 100)
        assert np.all(final_state.active == [1, 2, 4, 3])

    def test_bb_bounds_tes_problem(self):
        state = optimization_methods.bb_state([1], [0], [2, 3, 4]) # 1 active, 0 inactive
        linear = [np.arange(5)[None, :]]
        quadratic = [np.diag(np.arange(5))]
        max_l0 = 3
        max_el_current = 4

        mock_fun = mock.Mock(
            side_effect=[
                (np.array([4, 3, 2, 1]), 2), # First call (lower bound)
                (np.array([4, 3, 2, 1]), None), # second call (preparation upper bound)
                (np.array([4, 3, 2]), 4), # Third call (calculation upper bound)
            ]
        )
        ub, lb, child1, child2 = optimization_methods._bb_bounds_tes_problem(
            state, max_l0, linear, quadratic,
            max_el_current, mock_fun,
        )
        assert ub == 4
        assert lb == 2
        assert child2.inactive == [0, 2]
        assert child1.active == [1, 2]
        # Assertions realted to first call
        assert np.allclose(mock_fun.call_args_list[0][0][0][0], np.arange(1, 5))
        assert np.allclose(mock_fun.call_args_list[0][0][1][0], np.diag(np.arange(1, 5)))
        assert np.allclose(mock_fun.call_args_list[0][1]['extra_ineq'][0],
                           np.tile(np.array([0., 0.25, 0.25, 0.25]), 2))
        assert np.allclose(mock_fun.call_args_list[0][1]['extra_ineq'][1], 2)
        # Second call
        assert np.allclose(mock_fun.call_args_list[1][0][0][0], np.arange(1, 5))
        assert np.allclose(mock_fun.call_args_list[1][0][1][0], np.diag(np.arange(1, 5)))
        # Third call
        assert np.allclose(mock_fun.call_args_list[2][0][0][0], np.arange(1, 4))
        assert np.allclose(mock_fun.call_args_list[2][0][1][0], np.diag(np.arange(1, 4)))


class TestLinearElecConstrained:
    @pytest.mark.parametrize('init_startegy', ['compact', 'full'])
    def test_solve_feasible(self, init_startegy):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)
        targets = [0, 1]
        target_direction = [[1, 0, 0], [1, 0, 0]]

        max_total_current = 0.2
        max_el_current = 0.1
        target_mean = 5e-2
        n_elec = 4

        tes_problem = optimization_methods.TESLinearElecConstrained(
            n_elec, leadfield, max_total_current, max_el_current
        )
        tes_problem.add_linear_constraint(targets, target_direction, target_mean)

        x = tes_problem.solve(init_startegy=init_startegy)

        x_bf = None
        obj_bf = np.inf
        for c in itertools.combinations(range(6),  n_elec):
            l_ = tes_problem.l[:, c]
            Q_ = tes_problem.Q[np.ix_(c, c)]
            x_ = optimization_methods._linear_constrained_tes_opt(
                l_, tes_problem.target_means,
                Q_, max_el_current,
                max_total_current, log_level=0
            )
            obj = x_.dot(Q_).dot(x_)
            if np.all(l_.dot(x_) > target_mean*0.99) and obj < obj_bf:
                obj_bf = obj
                x_bf = np.zeros(6)
                x_bf[list(c)] = x_

        assert np.linalg.norm(x, 1) <= 2 * max_total_current + 1e-4
        assert np.linalg.norm(x, 0) <= n_elec
        assert np.all(np.abs(x) <= max_el_current + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(tes_problem.l.dot(x), target_mean)
        assert np.allclose(x_bf, x)

    def test_solve_infeasible(self):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)
        targets = [0, 1]
        target_direction = [[1, 0, 0], [1, 0, 0]]

        max_total_current = 0.2
        max_el_current = 0.1
        target_mean = 20
        n_elec = 4

        tes_problem = optimization_methods.TESLinearElecConstrained(
            n_elec, leadfield, max_total_current, max_el_current
        )
        tes_problem.add_linear_constraint(targets, target_direction, target_mean)

        x = tes_problem.solve()

        x_sp = optimize_comp(
            np.squeeze(tes_problem.l),
            np.ones((1, 6)),
            max_el_current,
            max_total_current
        )

        assert np.linalg.norm(x, 1) <= 2 * max_total_current + 1e-4
        assert np.linalg.norm(x, 0) <= n_elec
        assert np.all(np.abs(x) <= max_el_current + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.allclose(x_sp, x, atol=1e-4)

class TestLinearAngleElecConstrained:
    @pytest.mark.parametrize('init_startegy', ['compact', 'full'])
    def test_solve_feasible(self, init_startegy):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)

        max_total_current = 0.2
        max_el_current = 0.1
        target_mean = 1e-3
        max_angle = 40
        n_elec = 4

        tes_problem = optimization_methods.TESLinearAngleElecConstrained(
            n_elec, [0], [1, 0, 0], target_mean, max_angle,
            leadfield, max_total_current, max_el_current
        )
        x = tes_problem.solve()
        field = np.array([x[1:].dot(leadfield[..., i]) for i in range(3)]).T
        angle = np.rad2deg(np.arccos(field[0, 0]/np.linalg.norm(field[0])))

        obj_bf = np.inf
        for c in itertools.combinations(range(6),  n_elec):
            l_ = tes_problem.l[:, c]
            Q_ = tes_problem.Q[np.ix_(c, c)]
            Qnorm_ = tes_problem.Qnorm[np.ix_(c, c)]
            x_ = optimization_methods._linear_angle_constrained_tes_opt(
                l_, tes_problem.target_mean,
                Q_, max_el_current, max_total_current,
                Qnorm_, max_angle,
                log_level=0
            )
            obj = x_.dot(Q_).dot(x_)
            if np.all(l_.dot(x_) > target_mean*0.99) and obj < obj_bf:
                obj_bf = obj


        assert np.all(np.abs(x) <= max_el_current)
        assert np.all(np.linalg.norm(x) <= 2*max_total_current)
        assert np.isclose(np.sum(x), 0)
        assert angle <= max_angle
        assert np.isclose(tes_problem.l.dot(x), target_mean)
        assert x.dot(tes_problem.Q).dot(x) <= obj_bf + 1.0e-14

    def test_solve_angle_elec_infeasible(self):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)

        max_total_current = 0.2
        max_el_current = 0.1
        target_mean = 1.
        max_angle = 40
        n_elec = 4

        tes_problem = optimization_methods.TESLinearAngleElecConstrained(
            n_elec, [0], [1, 0, 0], target_mean, max_angle,
            leadfield, max_total_current, max_el_current
        )
        x = tes_problem.solve()
        field = np.array([x[1:].dot(leadfield[..., i]) for i in range(3)]).T
        angle = np.rad2deg(np.arccos(field[0, 0]/np.linalg.norm(field[0])))

        obj_bf = -np.inf
        for c in itertools.combinations(range(6),  n_elec):
            l_ = tes_problem.l[:, c]
            Q_ = tes_problem.Q[np.ix_(c, c)]
            Qnorm_ = tes_problem.Qnorm[np.ix_(c, c)]
            x_ = optimization_methods._linear_angle_constrained_tes_opt(
                l_, tes_problem.target_mean,
                Q_, max_el_current, max_total_current,
                Qnorm_, max_angle,
                log_level=10
            )
            obj = l_.dot(x_)
            if np.all(obj > obj_bf):
                obj_bf = obj

        assert np.all(np.abs(x) <= max_el_current)
        assert np.all(np.linalg.norm(x) <= 2*max_total_current)
        assert np.isclose(np.sum(x), 0)
        assert angle <= max_angle
        res = tes_problem.l.dot(x)
        assert np.allclose(res, obj_bf) or np.all(res > obj_bf)


class TestEqConstrained:
    def test_equality_constraints(self, optimization_variables):
        l, Q, P = optimization_variables
        A = np.ones((1, len(l)))
        b = np.zeros(1)
        x = optimization_methods._eq_constrained_QP(l, Q, A, b)

        ob = lambda x: l.dot(x) + .5 * x.dot(Q.dot(x))
        x0 = np.zeros(Q.shape[0])
        jac = lambda x: l + x.dot(Q)
        constraint = {}
        constraint['type'] = 'eq'
        constraint['fun'] = lambda x: A.dot(x) - b
        constraint['jac'] = lambda x: A
        constraints = (constraint, )

        res = scipy.optimize.minimize(
            ob, x0, jac=jac,
            constraints=constraints)
        x_sp = res.x

        assert np.isclose(A.dot(x), b)
        assert np.isclose(objective(l, Q, 0, x), objective(l, Q, 0, x_sp),
                          rtol=1e-4, atol=1e-4)
        assert np.allclose(x, x_sp, rtol=1e-4)


class TestConstrainedEigenvalue:
    def test_constrained_eigv(self):
        A = np.random.rand(5, 100)
        Q = A.dot(A.T)
        eigvals, eigvec = optimization_methods._constrained_eigenvalue(Q)
        assert np.allclose(np.sum(eigvec, axis=0), 0)
        assert np.allclose(eigvec.T.dot(eigvec)[1:, 1:], np.eye(4))
        assert np.allclose(eigvals, np.diagonal(eigvec.T.dot(Q).dot(eigvec)))


class TestTESNormContrained:
    def test_magn_constrained_tes_opt(self):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)

        tes_opt = optimization_methods.TESOptimizationProblem(leadfield)
        Q = tes_opt.Q
        Qnorm = optimization_methods._calc_Qnorm(leadfield, 0, np.ones(10))

        max_total_current = 0.2
        max_el_current = 0.1
        target_magn = 0.1

        x = optimization_methods._norm_constrained_tes_opt(
            Qnorm, target_magn, Q,
            max_el_current, max_total_current
        )

        energy_bf = np.inf
        directions = fibonacci_sphere(51)
        directions = directions[directions[:, 2] > 0]
        for d in directions:
            l = optimization_methods._calc_l(
                leadfield, 0, d, np.ones(10)
            )
            x_ = optimization_methods._linear_constrained_tes_opt(
                l[None, :], np.atleast_1d(target_magn),
                Q, max_el_current,
                max_total_current, log_level=0
            )
            norm_ = np.sqrt(x_.dot(Qnorm).dot(x_))
            energy_ = x_.dot(Q).dot(x_)
            if np.isclose(norm_, target_magn, rtol=1e-1) and energy_ < energy_bf:
                energy_bf = energy_

        assert np.all(np.abs(x) <= max_el_current)
        assert np.all(np.linalg.norm(x) <= 2*max_total_current)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(np.sqrt(x.dot(Qnorm).dot(x)), target_magn, rtol=1e-1)
        assert x.dot(Q).dot(x) < energy_bf

    def test_magn_constrained_tes_opt_infeasible(self):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)

        tes_opt = optimization_methods.TESOptimizationProblem(leadfield)
        Q = tes_opt.Q
        Qnorm = optimization_methods._calc_Qnorm(leadfield, 0, np.ones(10))

        max_total_current = 0.2
        max_el_current = 0.1
        target_magn = 20

        x = optimization_methods._norm_constrained_tes_opt(
            Qnorm[None, ...], target_magn, Q,
            max_el_current, max_total_current
        )

        norm_bf = -np.inf
        directions = fibonacci_sphere(51)
        directions = directions[directions[:, 2] > 0]
        for d in directions:
            l = optimization_methods._calc_l(
                leadfield, 0, d, np.ones(10)
            )
            x_ = optimization_methods._linear_constrained_tes_opt(
                l[None, :], np.atleast_1d(target_magn),
                Q, max_el_current,
                max_total_current, log_level=0
            )
            norm_ = np.sqrt(x_.dot(Qnorm).dot(x_))
            if norm_ > norm_bf:
                norm_bf = norm_

        assert np.all(np.abs(x) <= max_el_current)
        assert np.all(np.linalg.norm(x) <= 2*max_total_current)
        assert np.isclose(np.sum(x), 0)
        assert np.sqrt(x.dot(Qnorm).dot(x)) > norm_bf*0.999


    def test_magn_constrained_tes_opt_multi(self):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)

        tes_opt = optimization_methods.TESOptimizationProblem(leadfield)
        Q = tes_opt.Q
        Qnorm1 = optimization_methods._calc_Qnorm(leadfield, 0, np.ones(10))
        Qnorm2 = optimization_methods._calc_Qnorm(leadfield, 1, np.ones(10))
        Qnorm = np.stack([Qnorm1, Qnorm2])

        max_total_current = 0.2
        max_el_current = 0.1
        target_magn = np.array([0.1, 0.1])

        x = optimization_methods._norm_constrained_tes_opt(
            Qnorm, target_magn, Q,
            max_el_current, max_total_current
        )

        assert np.all(np.abs(x) <= max_el_current)
        assert np.all(np.linalg.norm(x) <= 2*max_total_current)
        assert np.isclose(np.sum(x), 0)
        assert np.allclose(np.sqrt(x.dot(Qnorm).dot(x)), target_magn, rtol=1e-1)


    def test_magn_constrained_tes_opt_multi_infeasible(self):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)

        tes_opt = optimization_methods.TESOptimizationProblem(leadfield)
        Q = tes_opt.Q
        Qnorm1 = optimization_methods._calc_Qnorm(leadfield, 0, np.ones(10))
        Qnorm2 = optimization_methods._calc_Qnorm(leadfield, 1, np.ones(10))
        Qnorm = np.stack([Qnorm1, Qnorm2])

        max_total_current = 0.2
        max_el_current = 0.1
        target_magn = np.array([0.1, 0.2])

        x = optimization_methods._norm_constrained_tes_opt(
            Qnorm, target_magn, Q,
            max_el_current, max_total_current
        )

        assert np.all(np.abs(x) <= max_el_current*1.001)
        assert np.all(np.linalg.norm(x) <= 2*max_total_current*1.001)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(np.sqrt(x.dot(Qnorm).dot(x))[0], target_magn[0], rtol=1e-1)

    def test_add_magn_constraint(self):
        A = np.random.random((2, 5, 3))
        targets = [0, 1]
        volumes = np.array([1, 2, 2, 2, 2])
        tes_problem = optimization_methods.TESNormConstrained(
            A, weights=volumes
        )
        tes_problem.add_norm_constraint(targets, 0.2)
        currents = np.array([-2, 1, 1])
        field_norm = np.linalg.norm(A.T.dot(currents[1:]), axis=0)
        avg_sq = np.average(field_norm[[0, 1]]**2, weights=[1, 2])
        assert np.allclose(avg_sq, currents.T @ tes_problem.Qnorm @ currents)
        assert np.allclose(0.2, tes_problem.target_means)

    def test_solve_magn_constraint(self):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)

        max_total_current = 0.2
        max_el_current = 0.1
        target_magn = np.array([0.1, 0.1])

        tes_opt = optimization_methods.TESNormConstrained(
            leadfield, max_total_current, max_el_current
        )
        tes_opt.add_norm_constraint(0, target_magn[0])
        tes_opt.add_norm_constraint(1, target_magn[1])
        x = tes_opt.solve()

        assert np.all(np.abs(x) <= max_el_current)
        assert np.all(np.linalg.norm(x) <= 2*max_total_current)
        assert np.isclose(np.sum(x), 0)
        assert np.allclose(np.sqrt(x.dot(tes_opt.Qnorm).dot(x)), target_magn, rtol=1e-1)


class TestMagnElecConstrained:
    @pytest.mark.parametrize('init_startegy', ['compact', 'full'])
    def test_solve_feasible(self, init_startegy):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)
        targets = [0, 1]

        max_total_current = 0.2
        max_el_current = 0.1
        target_mean = 5e-2
        n_elec = 4

        tes_problem = optimization_methods.TESNormElecConstrained(
            n_elec, leadfield, max_total_current, max_el_current
        )
        tes_problem.add_norm_constraint(targets, target_mean)

        x = tes_problem.solve(init_startegy=init_startegy)

        x_bf = None
        obj_bf = np.inf
        for c in itertools.combinations(range(6),  n_elec):
            Qnorm_ = tes_problem.Qnorm[0][np.ix_(c, c)]
            Q_ = tes_problem.Q[np.ix_(c, c)]
            x_ = optimization_methods._norm_constrained_tes_opt(
                Qnorm_, tes_problem.target_means,
                Q_, max_el_current,
                max_total_current, log_level=0
            )
            obj = x_.dot(Q_).dot(x_)
            if np.all(np.sqrt(x_.dot(Qnorm_).dot(x_)) > target_mean*0.99) and obj < obj_bf:
                obj_bf = obj
                x_bf = np.zeros(6)
                x_bf[list(c)] = x_

        assert np.linalg.norm(x, 1) <= 2 * max_total_current + 1e-4
        assert np.linalg.norm(x, 0) <= n_elec
        assert np.all(np.abs(x) <= max_el_current + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(np.sqrt(x @ tes_problem.Qnorm @ x), target_mean, rtol=1e-2)
        assert np.allclose(x_bf, x) or np.allclose(x_bf, -x)

    def test_solve_infeasible(self):
        np.random.seed(1)
        leadfield = np.random.random((5, 10, 3))
        np.random.seed(None)
        targets = [0, 1]

        max_total_current = 0.2
        max_el_current = 0.1
        target_mean = 20
        n_elec = 4

        tes_problem = optimization_methods.TESNormElecConstrained(
            n_elec, leadfield, max_total_current, max_el_current
        )
        tes_problem.add_norm_constraint(targets, target_mean)

        x = tes_problem.solve()

        assert np.linalg.norm(x, 1) <= 2 * max_total_current + 1e-4
        assert np.linalg.norm(x, 0) <= n_elec
        assert np.all(np.abs(x) <= max_el_current + 1e-4)
        assert np.isclose(np.sum(x), 0)


class TestDistributed:
    @pytest.mark.parametrize('max_el_current', [1e5, 1e-4])
    @pytest.mark.parametrize('max_total_current', [1e5, 1e-4])
    def test_least_squares_opt(self, max_el_current, max_total_current):
        np.random.seed(1)
        A = np.random.rand(30, 5)
        b = np.random.rand(30)
        np.random.seed(None)

        Q = A.T.dot(A)
        l = -2*b.dot(A)

        x = optimization_methods._least_squares_tes_opt(
            l, Q,
            max_el_current,
            max_total_current,
        )
        x_sp = optimize_lstsq(A, b, max_el_current, max_total_current)
        objective = lambda x: np.linalg.norm(A.dot(x) - b)**2
        assert np.linalg.norm(x, 1) <= 2 * max_total_current + 1e-4
        assert np.all(np.abs(x) <= max_el_current + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert objective(x) < objective(x_sp) *1.1

    def test_calc_l_Q(self):
        leadfield = np.zeros((5, 30, 3))
        for i in range(5):
            leadfield[i, i, :] = 1

        target_field = np.zeros((30, 3))
        target_field[0, 0] = 1

        tes_problem = optimization_methods.TESDistributed(
            leadfield, target_field
        )
        x = optimization_methods._eq_constrained_QP(
            tes_problem.l, 2*tes_problem.Q, np.ones((1, 6)), [0]
        )
        field = leadfield.T.dot(x[1:]).T
        assert np.allclose(field[0], 1/3)
        assert np.allclose(field[1:], 0)

    def test_calc_l_Q_weighted(self):
        leadfield = np.zeros((5, 30, 3))
        for i in range(5):
            leadfield[i, i, :] = 1

        target_field = np.zeros((30, 3))
        target_field[0, 0] = 1

        weights = np.ones((30, 3))
        weights[0] = [-2, -1, -1]

        tes_problem = optimization_methods.TESDistributed(
            leadfield, target_field, weights
        )
        x = optimization_methods._eq_constrained_QP(
            tes_problem.l, 2*tes_problem.Q, np.ones((1, 6)), [0]
        )
        field = leadfield.T.dot(x[1:]).T
        assert np.allclose(field[0], 1/2)
        assert np.allclose(field[1:], 0)

    @pytest.mark.parametrize('max_el_current', [1e5, 1e-2])
    @pytest.mark.parametrize('max_total_current', [1e5, 1e-2])
    def test_solve_scalar_w(self, max_el_current, max_total_current):
        np.random.seed(1)
        leadfield = np.random.rand(5, 30, 3)
        target_field = np.random.rand(30, 3)
        np.random.seed(None)

        weights = np.ones(30)
        weights[20:] = 0

        tes_problem = optimization_methods.TESDistributed(
            leadfield, target_field, weights, max_total_current,
            max_el_current
        )

        x = tes_problem.solve()

        P = np.linalg.pinv(np.vstack([-np.ones(5), np.eye(5)]))
        A = leadfield.reshape(5, -1).T.dot(P)
        A[60:] = 0
        b = target_field.reshape(-1)
        b[60:] = 0

        x_sp = optimize_lstsq(
            A, b,
            max_el_current, max_total_current
        )

        objective = lambda x: np.linalg.norm(A.dot(x) - b)**2
        assert np.linalg.norm(x, 1) <= 2 * max_total_current + 1e-4
        assert np.all(np.abs(x) <= max_el_current + 1e-4)
        assert np.isclose(np.sum(x), 0, atol=1e-6)
        assert objective(x) < objective(x_sp) *1.1


    @pytest.mark.parametrize('max_el_current', [1e5, 1e-2])
    @pytest.mark.parametrize('max_total_current', [1e5, 1e-2])
    def test_solve_vector_w(self, max_el_current, max_total_current):
        np.random.seed(1)
        leadfield = np.random.rand(5, 30, 3)
        target_field = np.random.rand(30, 3)
        np.random.seed(None)

        weights = np.zeros((30, 3))
        weights[:, 1] = 1.

        tes_problem = optimization_methods.TESDistributed(
            leadfield, target_field, weights, max_total_current,
            max_el_current
        )

        x = tes_problem.solve()

        P = np.linalg.pinv(np.vstack([-np.ones(5), np.eye(5)]))
        A = leadfield[..., 1].T.dot(P)
        b = target_field[..., 1]

        x_sp = optimize_lstsq(
            A, b,
            max_el_current, max_total_current
        )

        objective = lambda x: np.linalg.norm(A.dot(x) - b)**2
        assert np.linalg.norm(x, 1) <= 2 * max_total_current + 1e-4
        assert np.all(np.abs(x) <= max_el_current + 1e-4)
        assert np.isclose(np.sum(x), 0, atol=1e-6)
        assert objective(x) < objective(x_sp) *1.1


class TestDistributedElec:
    @pytest.mark.parametrize('max_el_current', [1e5, 1e-2])
    @pytest.mark.parametrize('max_total_current', [1e5, 1e-2])
    def test_solve(self, max_el_current, max_total_current):
        np.random.seed(1)
        leadfield = np.random.rand(5, 30, 3)
        target_field = np.random.rand(30, 3)
        np.random.seed(None)

        weights = np.zeros((30, 3))
        weights[:, 1] = 1.

        tes_problem = optimization_methods.TESDistributedElecConstrained(
            4, leadfield, target_field,
            weights, max_total_current,
            max_el_current
        )
        x = tes_problem.solve()

        assert np.linalg.norm(x, 1) <= 2 * max_total_current + 1e-4
        assert np.linalg.norm(x, 0) <= 4
        assert np.all(np.abs(x) <= max_el_current + 1e-4)

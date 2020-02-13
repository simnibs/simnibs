import functools
import numpy as np
import scipy.optimize
import pytest
from simnibs.optimization import optimization_methods

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

class TestLinearContrainedBB:
    def test_bb_lower_bound(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 4e-3
        f = 1e-2
        max_active = 2
        state = optimization_methods.bb_state([1], [2], [3, 0, 4])
        x, val = optimization_methods._linear_constrained_bb_lower_bound(
            np.atleast_2d(l), Q, np.atleast_1d(f), m1, m,
            max_active, state)

        s = np.ones(len(l)) / m
        s[state.active] = 0
        select = state.active + state.unassigned
        l_ = l[select]
        Q_ = Q[np.ix_(select, select)]
        A_ = A[select]
        s_ = np.diag(s[select])

        ob = lambda x: .5 * x.dot(Q_.dot(x))
        x0 = np.zeros(Q_.shape[0])
        jac = lambda x: x.dot(Q_)
        #hess = lambda x: Q
        m_constraint = {}
        m_constraint['type'] = 'eq'
        m_constraint['fun'] = lambda x: l_.dot(x) - f
        m_constraint['jac'] = lambda x: l_

        e_constraint = {}
        e_constraint['type'] = 'eq'
        e_constraint['fun'] = lambda x: A_.dot(x)
        e_constraint['jac'] = lambda x: A_

        bounds = scipy.optimize.Bounds(-m, m)

        iq_constraint = {}
        iq_constraint['type'] = 'ineq'
        iq_constraint['fun'] = lambda x: 2 * m1 - np.linalg.norm(x, 1)

        r_constraint = {}
        r_constraint['type'] = 'ineq'
        r_constraint['fun'] = lambda x: 1 - np.linalg.norm(s_*x, 1)

        constraints = (m_constraint, e_constraint, iq_constraint, r_constraint)
        res = scipy.optimize.minimize(
            ob, x0, jac=jac,
            bounds=bounds,
            constraints=constraints)
        x_sp = np.zeros(len(l))
        x_sp[select] = res.x

        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(x[state.inactive], 0)
        assert np.isclose(l.dot(x), f)
        assert np.allclose(x.dot(Q.dot(x)), x_sp.dot(Q.dot(x_sp)), rtol=1e-4, atol=1e-5)

    def test_bb_upper_bound(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 4e-3
        f = 1e-3
        max_active = 2
        state = optimization_methods.bb_state([1], [2], [3, 0, 4])
        x_ub, ub = optimization_methods._linear_constrained_bb_upper_bound(
            np.atleast_2d(l), Q, np.atleast_1d(f), m1, m, max_active, state
        )
        assert np.linalg.norm(x_ub, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x_ub) <= m + 1e-4)
        assert np.isclose(np.sum(x_ub), 0)
        assert np.isclose(x_ub[state.inactive], 0)
        assert np.isclose(l.dot(x_ub), f)


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


class TestOptimizeTargetIntensity:
    def test_bound_constraints(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        x = optimization_methods.optimize_field_component(l, max_el_current=m)
        x_sp = optimize_comp(l, A, m)
        assert np.all(np.abs(x) <= m+1e-6)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), l.dot(x_sp),
                          rtol=1e-4, atol=1e-4)

    def test_l1_constraint(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        x = optimization_methods.optimize_field_component(l, max_total_current=m)
        x_sp = optimize_comp(l, A, max_total_current=m)
        assert np.linalg.norm(x, 1) <= 2 *m + 1e-4
        assert np.isclose(l.dot(x), l.dot(x_sp),
                          rtol=1e-4, atol=1e-4)
        assert np.isclose(np.sum(x), 0)

    def test_both_constraints(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 3e-3
        x = optimization_methods.optimize_field_component(l, max_el_current=m,
                                                          max_total_current=m1)
        x_sp = optimize_comp(l, A, m, m1)
        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m+1e-6)
        assert np.isclose(l.dot(x), l.dot(x_sp),
                          rtol=1e-4, atol=1e-4)
        assert np.isclose(np.sum(x), 0)


class TestOptimizeFocality:
    def test_bound_constraints_feasible(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        f = 2e-2
        x = optimization_methods.optimize_focality(l, Q, f, max_el_current=m)
        x_sp = optimize_focality(l, Q, f, m)
        assert np.all(np.abs(x) <= m)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), f)
        assert np.allclose(x.dot(Q.dot(x)), x_sp.dot(Q.dot(x_sp)), rtol=1e-4, atol=1e-5)

    def test_bound_constraints_feasible(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 3e-3
        f = 2e-2
        x = optimization_methods.optimize_focality(l, Q, f, max_total_current=m)
        x_sp = optimize_focality(l, Q, f, max_total_current=m)
        assert np.linalg.norm(x, 1) <= 2 * m + 1e-4
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), f)
        assert np.allclose(x.dot(Q.dot(x)), x_sp.dot(Q.dot(x_sp)), rtol=1e-4, atol=1e-5)

    def test_both_constraints_feasible(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 4e-3
        f = 2e-2
        x = optimization_methods.optimize_focality(l, Q, f, max_el_current=m,
                                                   max_total_current=m1)

        x_sp = optimize_focality(
            l, Q, f, max_el_current=m, max_total_current=m1)

        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), f)
        assert np.allclose(x.dot(Q.dot(x)), x_sp.dot(Q.dot(x_sp)), rtol=1e-4, atol=1e-5)

    def test_both_constraints_infeasible(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 4e-3
        f = 4e-2
        x = optimization_methods.optimize_focality(l, Q, f,
                                                   max_el_current=m,
                                                   max_total_current=m1)

        x_sp = optimize_comp(l, np.ones_like(l), max_el_current=m, max_total_current=m1)
        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.allclose(l.dot(x), l.dot(x_sp), rtol=1e-4, atol=1e-5)

    def test_both_limit_angle_inactive(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        Q_in = 5e1 * np.eye(len(l))
        max_angle = 45
        m = 2e-3
        m1 = 4e-3
        f = .02
        x = optimization_methods.optimize_focality(l, Q, f,
                                                   max_el_current=m,
                                                   max_total_current=m1,
                                                   Qin=Q_in,
                                                   max_angle=max_angle)
        x_sp = optimize_focality(
            l, Q, f, max_el_current=m, max_total_current=m1,
            Qin=Q_in, max_angle=max_angle)

        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), f)
        assert np.arccos(l.dot(x) / np.sqrt(x.dot(Q_in).dot(x))) <= np.deg2rad(max_angle)
        assert np.allclose(x.dot(Q.dot(x)), x_sp.dot(Q.dot(x_sp)), rtol=1e-4, atol=1e-5)

    def test_both_limit_angle_limit(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        Q_in = 5e1 * np.eye(len(l))
        max_angle = 25
        m = 2e-3
        m1 = 4e-3
        f = .02
        x = optimization_methods.optimize_focality(l, Q, f,
                                                   max_el_current=m,
                                                   max_total_current=m1,
                                                   Qin=Q_in,
                                                   max_angle=max_angle)
        x_sp = optimize_focality(
            l, Q, f, max_el_current=m, max_total_current=m1,
            Qin=Q_in, max_angle=max_angle)
        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), f)
        assert np.arccos(l.dot(x) / np.sqrt(x.dot(Q_in).dot(x))) <= np.deg2rad(max_angle)
        assert np.allclose(x.dot(Q.dot(x)), x_sp.dot(Q.dot(x_sp)), rtol=1e-4, atol=1e-5)

    def test_both_limit_angle_reduce_target(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        Q_in = 5e1 * np.eye(len(l))
        max_angle = 20
        m = 2e-3
        m1 = 4e-3
        f = .02
        x = optimization_methods.optimize_focality(l, Q, f,
                                                   max_el_current=m,
                                                   max_total_current=m1,
                                                   Qin=Q_in,
                                                   max_angle=max_angle)

        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.arccos(l.dot(x) / np.sqrt(x.dot(Q_in).dot(x))) <= np.deg2rad(max_angle)

    def test_both_limit_angle_Q_iteration(self, optimization_variables_avg_QCQP):
        l, Q, A, Q_in = optimization_variables_avg_QCQP
        # l, Q, A = optimization_variables_avg
        # Q_in = 6 * np.eye(len(l)) + np.outer(l, l)
        max_angle = 20
        m = 2e-3
        m1 = 4e-3
        f = .01
        x = optimization_methods.optimize_focality(l, Q, f,
                                                   max_el_current=m,
                                                   max_total_current=m1,
                                                   Qin=Q_in,
                                                   max_angle=max_angle)
        x_sp = optimize_focality(
            l, Q, f, max_el_current=m, max_total_current=m1,
            Qin=Q_in, max_angle=max_angle)
        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), f)
        assert np.arccos(l.dot(x) / np.sqrt(x.dot(Q_in).dot(x))) <= np.deg2rad(max_angle)
        assert np.allclose(x.dot(Q.dot(x)), x_sp.dot(Q.dot(x_sp)), rtol=1e-4, atol=1e-5)

    def test_both_limit_angle_infeasible_field(self, optimization_variables_avg_QCQP):
        l, Q, A, Q_in = optimization_variables_avg_QCQP
        max_angle = 15
        m = 2e-3
        m1 = 4e-3
        f = 2
        x = optimization_methods.optimize_focality(l, Q, f,
                                                   max_el_current=m,
                                                   max_total_current=m1,
                                                   Qin=Q_in,
                                                   max_angle=max_angle)

        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.arccos(l.dot(x) / np.sqrt(x.dot(Q_in).dot(x))) <= np.deg2rad(max_angle)


    def test_limit_nr_electrodes_infeasible_comp(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 4e-3
        f = 2e3
        n = 3
        x = optimization_methods.optimize_focality(l, Q, f, max_el_current=m,
                                                   max_total_current=m1,
                                                   max_active_electrodes=n)

        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.linalg.norm(x, 0) <= n
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)

    def test_limit_nr_electrodes_feasible_comp(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 4e-3
        f = 1e-2
        n = 3
        x = optimization_methods.optimize_focality(l, Q, f, max_el_current=m,
                                                   max_total_current=m1,
                                                   max_active_electrodes=n)

        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.linalg.norm(x, 0) <= n
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), f)

    def test_limit_nr_angle(seld, optimization_variables_avg_QCQP):
        l, Q, A, Q_in = optimization_variables_avg_QCQP
        max_angle = 15
        m = 2e-3
        m1 = 4e-3
        f = 2
        n = 4
        x = optimization_methods.optimize_focality(l, Q, f,
                                                   max_el_current=m,
                                                   max_total_current=m1,
                                                   Qin=Q_in,
                                                   max_angle=max_angle,
                                                   max_active_electrodes=n)

        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.linalg.norm(x, 0) <= n
        assert np.arccos(l.dot(x) / np.sqrt(x.dot(Q_in).dot(x))) <= np.deg2rad(max_angle)

    def test_limit_nr_angle_change_Q(seld, optimization_variables_avg_QCQP):
        l, Q, A, Q_in = optimization_variables_avg_QCQP
        max_angle = 15
        m = 2e-3
        m1 = 4e-3
        f = .01
        n = 4
        x = optimization_methods.optimize_focality(l, Q, f,
                                                   max_el_current=m,
                                                   max_total_current=m1,
                                                   Qin=Q_in,
                                                   max_angle=max_angle,
                                                   max_active_electrodes=n)

        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.linalg.norm(x, 0) <= n
        assert np.arccos(l.dot(x) / np.sqrt(x.dot(Q_in).dot(x))) <= np.deg2rad(max_angle)

    def test_2_targets_field_component(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        l2 = l[::-1]
        l = np.vstack([l ,l2])
        m = 2e-3
        m1 = 4e-3
        x = optimization_methods.optimize_field_component(l, max_el_current=m,
                                                          max_total_current=m1)

        l_avg = np.average(l, axis=0)
        x_sp = optimize_comp(l_avg, np.ones_like(l2), max_el_current=m, max_total_current=m1)

        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-6)
        assert np.isclose(l_avg.dot(x), l_avg.dot(x_sp),
                          rtol=1e-4, atol=1e-4)
        assert np.isclose(np.sum(x), 0)

    def test_2_targets_feasible(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        l2 = l[::-1]
        l = np.vstack([l, l2])
        m = 2e-3
        m1 = 4e-3
        f = [1e-2, 1e-2]
        x = optimization_methods.optimize_focality(l, Q, f, max_el_current=m,
                                                   max_total_current=m1)
        x_sp = optimize_focality(
            l, Q, f, max_el_current=m, max_total_current=m1)
        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.allclose(l.dot(x), f)
        assert np.allclose(x.dot(Q.dot(x)), x_sp.dot(Q.dot(x_sp)), rtol=1e-4, atol=1e-5)


    def test_2_targets_limit_nr_electrodes(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        l2 = l[::-1]
        l = np.vstack([l, l2])
        m = 2e-3
        m1 = 4e-3
        f = np.array([1e-2, 1e-2])
        n = 3
        x = optimization_methods.optimize_focality(
            l, Q, f, max_el_current=m,
            max_total_current=m1,
            max_active_electrodes=n)
        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.linalg.norm(x, 0) <= n
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.allclose(l.dot(x), f)


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


class TestBB:
    def test_bb_node(self):
        # Function with gives out the bounds as well as the new sets to split
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

    def test_bb_lower_bound(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 4e-3
        f = 1e-2
        max_active = 2
        state = optimization_methods.bb_state([1], [2], [3, 0, 4])
        x, val = optimization_methods._bb_lower_bound(
            np.atleast_2d(l), Q, f, m1, m,
            max_active, state)

        s = np.ones(len(l)) / m
        s[state.active] = 0
        select = state.active + state.unassigned
        l_ = l[select]
        Q_ = Q[np.ix_(select, select)]
        A_ = A[select]
        s_ = np.diag(s[select])

        ob = lambda x: .5 * x.dot(Q_.dot(x))
        x0 = np.zeros(Q_.shape[0])
        jac = lambda x: x.dot(Q_)
        #hess = lambda x: Q
        m_constraint = {}
        m_constraint['type'] = 'eq'
        m_constraint['fun'] = lambda x: l_.dot(x) - f
        m_constraint['jac'] = lambda x: l_

        e_constraint = {}
        e_constraint['type'] = 'eq'
        e_constraint['fun'] = lambda x: A_.dot(x)
        e_constraint['jac'] = lambda x: A_

        bounds = scipy.optimize.Bounds(-m, m)

        iq_constraint = {}
        iq_constraint['type'] = 'ineq'
        iq_constraint['fun'] = lambda x: 2 * m1 - np.linalg.norm(x, 1)

        r_constraint = {}
        r_constraint['type'] = 'ineq'
        r_constraint['fun'] = lambda x: 1 - np.linalg.norm(s_*x, 1)

        constraints = (m_constraint, e_constraint, iq_constraint, r_constraint)
        res = scipy.optimize.minimize(
            ob, x0, jac=jac,
            bounds=bounds,
            constraints=constraints)
        x_sp = np.zeros(len(l))
        x_sp[select] = res.x

        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(x[state.inactive], 0)
        assert np.isclose(l.dot(x), f)
        assert np.allclose(x.dot(Q.dot(x)), x_sp.dot(Q.dot(x_sp)), rtol=1e-4, atol=1e-5)

    def test_bb_upper_bound(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 4e-3
        f = 1e-3
        max_active = 2
        state = optimization_methods.bb_state([1], [2], [3, 0, 4])
        x_ub, ub = optimization_methods._bb_upper_bound(
            np.atleast_2d(l), Q, f, m1, m, max_active, state)

        assert np.linalg.norm(x_ub, 1) <= 2 * m1 + 1e-4
        assert np.all(np.abs(x_ub) <= m + 1e-4)
        assert np.isclose(np.sum(x_ub), 0)
        assert np.isclose(x_ub[state.inactive], 0)
        assert np.isclose(l.dot(x_ub), f)

    def test_bb_split(self):
        state = optimization_methods.bb_state([1], [2], [3, 0, 4])
        x = np.array([4, 3, 0, 1, -5])
        child1, child2 = optimization_methods._bb_split(x, state)
        assert child1.active == [1, 4]
        assert child1.inactive == [2]
        assert child2.active == [1]
        assert child2.inactive == [2, 4]

    def test_bb_full(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 4e-3
        f = 1e-2
        max_active = 2
        init = optimization_methods.bb_state([], [], list(range(len(l))))
        bf = functools.partial(
            optimization_methods._bb_bounds_function,
            np.atleast_2d(l), Q, f, m1, m, max_active)
        final_state = optimization_methods._branch_and_bound(
            init, bf, 1e-2, 100)
        x = final_state.x_lb
        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.linalg.norm(x, 0) <= max_active
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), f)

    def test_constrained_l0_branch_and_bound(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 4e-3
        f = 1e-2
        max_active = 2
        x = optimization_methods._constrained_l0_branch_and_bound(
            l, Q, f, m1, m, max_active)
        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.linalg.norm(x, 0) <= max_active
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)
        assert np.isclose(l.dot(x), f)

    def test_constrained_l0_branch_and_bound_ac(self, optimization_variables_avg):
        l, Q, A = optimization_variables_avg
        m = 2e-3
        m1 = 4e-3
        f = 1e-2
        max_active = 3
        x = optimization_methods._constrained_l0_branch_and_bound(
            l, Q, f, m1, m, max_active, start_inactive=[2], start_active=[1])
        assert np.linalg.norm(x, 1) <= 2 * m1 + 1e-4
        assert np.linalg.norm(x, 0) <= max_active
        assert np.isclose(x[2], 0)
        assert np.all(np.abs(x) <= m + 1e-4)
        assert np.isclose(np.sum(x), 0)

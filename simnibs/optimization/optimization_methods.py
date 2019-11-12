# -*- coding: utf-8 -*-\
import copy
import functools
import numpy as np
import scipy.optimize
import scipy.linalg

from simnibs.utils.simnibs_logger import logger

def _active_set_QP(l, Q, C, d, x0, eps=1e-5, A=None, b=None):
    ''' Solves the problem
    minimize l^T x + 1/2 x^T Q x
    subject to  Cx <= d
                Ax = b
    '''
    # get the active set:
    x = np.copy(x0)
    n = len(x0)
    n_iter = 0
    active = np.abs(C.dot(x) - d) < eps
    if np.any(np.abs(C.dot(x) - d) < -eps):
        raise ValueError('Infeasible Start!')
    if A is not None and not np.allclose(A.dot(x0), b, atol=eps):
        raise ValueError('infeasible start!')

    if A is None:
        A = np.empty((0, len(x0)))
    else:
        A = np.atleast_2d(A)
        if A.shape[0] > A.shape[1]:
            A = A.T

    max_iter = max(200, 2 * len(x0))
    added_iq_constraint = -1
    duals_inequality = np.zeros(C.shape[0])
    Active = np.vstack([A, C[active, :]])
    n_active = Active.shape[0]
    if n_active > 0:
        ZY, R = np.linalg.qr(Active.T, 'complete')

    else:
        ZY = np.eye(n)
        R = np.zeros((n, n))

    Y = ZY[:, :n_active]
    Z = ZY[:, n_active:]
    '''
    if n_active >= n:
        L = None
    else:
        L = np.linalg.cholesky(Z.T.dot(Q).dot(Z))
    '''
    while n_iter <= max_iter:
        l_i = l + Q.dot(x)
        Y = ZY[:, :n_active]
        Z = ZY[:, n_active:]
        if n_active >= n:
            p = np.zeros(n)
        else:
            try:
                w_z = np.linalg.solve(Z.T.dot(Q).dot(Z), - Z.T.dot(l_i))
            except np.linalg.LinAlgError:
                w_z = np.linalg.lstsq(Z.T.dot(Q).dot(Z), - Z.T.dot(l_i), rcond=None)[0]
            p = Z.dot(w_z)
            p = np.squeeze(p.T)
        # If no update, check the duals and try to terminate
        if np.all(np.abs(p) < eps):
            # calculate duals
            if Active.shape[0] > 0:
                try:
                    duals = scipy.linalg.solve_triangular(R[:n_active, :], -Y.T.dot(l_i))
                except (scipy.linalg.LinAlgError, ValueError, TypeError):
                    duals = np.linalg.lstsq(Active.T, -l_i, rcond=None)[0]

                duals = np.atleast_1d(np.squeeze(duals.T))
                duals_inequality = duals[A.shape[0]:]

            if np.all(duals_inequality > -eps):
                break

            min_C_dual = np.argmin(duals_inequality) + A.shape[0]
            # anti-cycling
            if n_iter > .5 * max_iter and min_C_dual == Active.shape[0] - 1:
                break

            Active = np.delete(Active, min_C_dual, axis=0)
            ZY, R = scipy.linalg.qr_delete(ZY, R,
                                           min_C_dual,
                                           which='col')
            n_active -= 1

        else:
            den = C.dot(p)
            s = den > 1e-9
            alpha_violation = (d[s] - C[s, :].dot(x)) / den[s]
            # if no constraint can be violated, step size is one
            if len(alpha_violation) == 0 or np.min(alpha_violation >= 1):
                alpha = 1.

            # if one constraint can be violated, calculate the maximum step size
            else:
                alpha = np.min(alpha_violation)
                went_to_bounds = np.argmin(alpha_violation)
                indices = np.where(s)[0]
                added_iq_constraint = indices[went_to_bounds]
                Active = np.vstack([Active, C[added_iq_constraint, :]])
                # Update the QR decomposition
                if n_active > 0:
                    ZY, R = scipy.linalg.qr_insert(ZY, R,
                                                   C[added_iq_constraint, :], n_active,
                                                   which='col')
                    n_active += 1
                else:
                    ZY, R = np.linalg.qr(np.atleast_2d(C[added_iq_constraint, :]).T, 'complete')
                    n_active += 1
            x = x + alpha * p

        n_iter += 1
    if n_iter >= max_iter:
        #raise ValueError('Maximal number of iterations reached!')
        pass

    return x


def _constrained_angle(l, Q, target_mean,
                       max_total_current,
                       max_el_current,
                       Qin, max_angle,
                       max_active_electrodes=None,
                       eps=1e-5,
                       eps_angle=1e-1,
                       log_level=20):
    logger.log(log_level, 'Running optimization with angle constraint')
    max_iter = 20
    # Try to find where we can find values bellow and above the target
    it = 0
    above = 0

    def angle(x):
        tan = np.sqrt(np.abs(x.dot(Qin).dot(x) - l.dot(x) ** 2))
        # the abs is here for numerical stability
        return np.abs(np.arctan2(tan, l.dot(x)))

    # Check if the maximal focality solution alreay fulfills the constraint
    x = optimize_focality(l, Q,
                          target_mean, max_total_current,
                          max_el_current,
                          max_active_electrodes=max_active_electrodes,
                          log_level=10)

    if angle(x) <= max_angle:
        logger.log(log_level, 'Max focality field fullfills angle constraint')
        return x

    x_above = np.copy(x)
    target_field = l.dot(x)

    # calculate the smallest angle, given a fixed intensity
    def calc_smallest_angle(alpha):
        x_l = optimize_focality(l, Qin,
                                alpha * target_field, max_total_current,
                                max_el_current,
                                max_active_electrodes=max_active_electrodes,
                                log_level=10)
        return x_l

    # if the target field is not reached, we don't need to calculate the smallest angle,
    # as the set of feasible solutions only contains one element
    if target_field >= target_mean * (1 - eps):
        x = calc_smallest_angle(1.)

    f = angle(x)

    # if we cannot reduce the angle to the target while keeting l^t x at the target
    # intensity, reduce the target intensity untill it's achievable.
    if f > max_angle:
        logger.log(log_level, "Target intensity can't be reached, reducing it")
        above = 1.
        f_above = f
        bellow = None
        f_bellow = None
        # lowers the target mean untill the desired angle can be achieved
        alpha = 1.
        it = 0
        # find lower bound
        while not (f > max_angle * (1 - eps_angle) and f < max_angle):
            if bellow is not None:
                alpha = above + \
                    (max_angle * (1 - eps_angle * .5) - f_above) * \
                    (bellow - above)/(f_bellow - f_above)
            else:
                alpha *= .8
            x = calc_smallest_angle(alpha)
            f = angle(x)
            logger.log(log_level,
                       '{0} alpha: {1}, angle: {2}, max_angle: {3}'.format(
                        it, alpha, f, max_angle))
            if f < max_angle:
                bellow = alpha
                f_bellow = f
                x_bellow = np.copy(x)
            else:
                above = alpha
                f_above = f
            it += 1
            if it > max_iter:
                if bellow is None:
                    return x
                else:
                    return x_bellow
        return x

    # In this case, we know that by minimizing x^t Qin x while keeping l^t x = t, we can
    # achieve the bound

    # find a combination between Q and Qin that maximizes focality while keeping Qin in
    # the bound
    else:
        logger.log(log_level, "Target intensity reached, optimizing focality with angle constraint")
        f_bellow = f
        bellow = 1.
        x_bellow = np.copy(x)
        f_above = angle(x_above)
        above = 0

        # Start bisection
        it = 0
        alpha = 1
        while not (f > max_angle * (1 - eps_angle) and f < max_angle):
            alpha = above + \
                (max_angle * (1 - eps_angle * .5) - f_above) * \
                (bellow - above)/(f_bellow - f_above)
            x = optimize_focality(l, (1 - alpha) * Q + alpha * Qin,
                                  target_field, max_total_current,
                                  max_el_current,
                                  max_active_electrodes=max_active_electrodes,
                                  log_level=10)
            f = angle(x)
            logger.log(log_level,
                       '{0} alpha: {1}, angle: {2}, max_angle: {3}'.format(
                        it, alpha, f, max_angle))
            if f > max_angle:
                above = alpha
                f_above = f
            else:
                bellow = alpha
                f_bellow = f
                x_bellow = np.copy(x)

            it += 1
            if it > max_iter:
                return x_bellow

        return x

def _constrained_l0(l, Q, target_mean,
                    max_total_current,
                    max_el_current,
                    max_l0, eps=1e-5,
                    method='bb_compact',
                    log_level=20):
    logger.log(log_level, 'Running optimization with constrained number of electrodes')
    max_l0 = int(max_l0)
    if max_total_current is None and max_el_current is None:
        raise ValueError('at least one of max_total_current or max_el_current must not be None')

    if max_el_current is not None:
        if max_total_current is not None:
            max_total_current = min(max_total_current, max_l0 * max_el_current / 2.)
        else:
            max_total_current = max_l0 * max_el_current
    else:
        max_el_current = max_total_current

    if method == 'proj':
        x = optimize_focality(l, Q,
                              target_mean, max_total_current,
                              max_el_current, log_level=10)
        active = np.argsort(np.abs(x))[-max_l0:]
        x_active = optimize_focality(l[:, active], Q[np.ix_(active, active)],
                                     target_mean, max_total_current,
                                     max_el_current, log_level=10)
        x = np.zeros(l.shape[1])
        x[active] = x_active
        return x

    elif method == 'bb_full':
        x = _constrained_l0_branch_and_bound(
            l, Q, target_mean, max_total_current, max_el_current,
            max_l0, eps_bb=1e-1, max_bb_iter=200, log_level=log_level)

    elif method == 'bb_compact':
        logger.log(log_level, 'Using the compact Branch-and-bound method')
        x = optimize_focality(l, Q,
                              target_mean, 3 * max_total_current,
                              max_el_current, log_level=10)
        active = np.abs(x) > 1e-3 * max_el_current
        x_active = _constrained_l0_branch_and_bound(
            l[:, active], Q[np.ix_(active, active)],
            target_mean, max_total_current, max_el_current,
            max_l0, eps_bb=1e-1, max_bb_iter=200, log_level=log_level)

        x = np.zeros(l.shape[1])
        x[active] = x_active

    else:
        raise Exception('Uknown method')

    return x


def target_matrices(leadfield, target_indices, target_direction, weights, avg_reference=True):
    ''' Calculate a matrices l and Qin such that
    l.dot(x) = a, x.dot(Qin).dot(x) = "b"
    "x" are the electrode currents
    "a" is the average electric in the target_direction across the target_indices.
    "b" is the average squared electric field norm in the target_indices

    Parameters
    ------------
    leadfield: ndarray (n_electrodes x n_points x 3)
        Leadfield
    target_indices: list of ints
        Indices of targets. If a 2 dimensional list the output "l" is also
        2-dimensional
    target_direction: ndarray
        The electric field direction to be optimized for each target
    weights: array
        Weight for each element (each leadfield column)
    avg_reference: bool (optional)
        Wether or not to re-reference to an average frame. Default: True

    Returns
    ----------
    l: np.ndarray
        Matrix calculating the average electric field
    Q_in: np.ndarray
        Matrix calculating the average squared electric field norm
    '''
    lf = leadfield
    target_indices = np.atleast_1d(target_indices)
    target_direction = np.atleast_2d(target_direction)
    if target_direction.shape[1] != 3:
        target_direction = target_direction.T
    if len(target_indices) != len(target_direction):
        raise ValueError('Please define one direction per target')
    if target_direction.shape[1] != 3:
        raise ValueError('A direction must have 3 dimentions')
    target_direction = target_direction/np.linalg.norm(target_direction, axis=1)[:, None]

    lf_t = lf[:, target_indices]
    w_idx = weights[target_indices]
    Q_in = sum(
        lf_t[..., i].dot((lf_t[..., i]*w_idx[None, :]).T) for i in range(3))
    Q_in /= np.sum(w_idx)
    lf_t *= w_idx[None, :, None]
    lf_t /= np.sum(w_idx)
    l = np.einsum('ijk, jk -> i', lf_t, target_direction)

    if avg_reference:
        P = np.linalg.pinv(
            np.vstack([-np.ones(len(l)), np.eye(len(l))]))
        l = l.dot(P)
        Q_in = P.T.dot(Q_in).dot(P)

    return l, Q_in


def energy_matrix(leadfield, weights, avg_reference=True):
    ''' Calculate the energy matrix for optimization
    x.dot(Q.dot(x)) = e
    "x" are the electrode currents
    "e" is the average squared electric field norm

    Parameters
    -----------
    leadfield: ndarray (n_electrodes x n_points x 3)
        Leadfield
    weights: array
        weith for each element (each leadfield column), for example volume or area
    avg_reference: bool (optional)
        Wether or not to re-reference to an average frame. Default: True


    Returns:
    ----------
    l: np.ndarray
        linear part of objective, the mean field in the target indices and target
        direction
    Q_out: np.ndarray
        Quadratic part of objective, the mean energy in the head
    Q_in: np.ndarray
        The squared average of the norm of the field in the target area in the target area
    '''

    lf = leadfield
    Q = sum(lf[..., i].dot((lf[..., i] * weights).T) for i in range(3))
    Q /= np.sum(weights)

    if avg_reference:
        P = np.linalg.pinv(
            np.vstack([-np.ones(Q.shape[0]), np.eye(Q.shape[0])]))
        Q = P.T.dot(Q).dot(P)

    return Q


def optimize_field_component(l, max_total_current=None, max_el_current=None):
    ''' Optimize the intensity of the field in the given direction without regard to
    focality
    
    Parameters
    ------------
        l: np.ndarray
            Linear objective, obtained from target_matrices.
            If "l" is muti-targeted, optimizes the average of the targets
        max_total_current: float (optional)
            Maximal current flow though all electrodes. Default: No maximum
        max_el_current: float (optional)
            Maximal current flow though each electrode. Default: No maximum

    Returns
    --------
        x: np.ndarray
            Optimal electrode currents
    '''
    if max_total_current is None and max_el_current is None:
        raise ValueError('Please define a maximal total current or maximal electrode ' +
                         'current')

    l = np.average(np.atleast_2d(l), axis=0)
    n = l.shape[0]

    if max_total_current is not None:
        A_ub = np.ones((1, 2 * n))
        b_ub = np.array([2 * max_total_current])
    else:
        A_ub = None
        b_ub = None

    l_ = np.hstack([l, -l])
    A_eq = np.hstack([np.ones((1, n)), -np.ones((1, n))])
    b_eq = np.array([0.])
    sol = scipy.optimize.linprog(-l_, A_ub, b_ub, A_eq, b_eq,
                                 bounds=(0, max_el_current))
    x_ = sol.x

    return x_[:n] - x_[n:]

def _lp_variables(l, target_mean, max_total_current, max_el_current):
        n = l.shape[1]
        if max_el_current is None and max_total_current is None:
            raise ValueError(
                'max_el_current and max_total_current can be simultaneously None')
        if max_total_current is not None:
            A_ub = [np.ones((1, 2 * n))]
            b_ub = [2 * max_total_current]
        else:
            A_ub = []
            b_ub = []
        #Constraint on target intensity
        l_ = np.hstack([l, -l])
        # the LP will maximize the average of all targets, and limit the electric field
        # at each individual target
        l_avg = np.average(l_, axis=0)
        A_ub = np.vstack(A_ub + [l_])
        b_ub = np.hstack(b_ub + [target_mean])
        A_eq = np.hstack([np.ones((1, n)), -np.ones((1, n))])
        b_eq = np.array([0.])
        bounds = (0, max_el_current)
        return l_avg, A_ub, b_ub, A_eq, b_eq, bounds

def _qp_variables(l, Q, target_mean, max_total_current, max_el_current):
        n = l.shape[1]
        # Bound constraints
        C = np.vstack([np.eye(n), -np.eye(n)])
        # Equality constraints
        A = np.ones(n)
        if max_el_current is not None:
            d = max_el_current * np.ones(C.shape[0])
            eps = 1e-3 * max_el_current
        else:
            d = 1e10 * np.ones(C.shape[0])
            max_el_current = 1e10

        if max_total_current is not None:
            max_l1 = 2 * max_total_current
            eps = 1e-6 * max_total_current
        else:
            max_l1 = 1e10

        Q_ = np.vstack([np.hstack([Q, -Q]), np.hstack([-Q, Q])])
        A_ = np.vstack([np.hstack([A, -A]), np.hstack([l, -l])])
        b_ = np.hstack([0, target_mean])
        C_ = np.vstack([np.eye(2 * n), - np.eye(2 * n), np.ones(2 * n)])
        d_ = np.hstack([d, np.zeros(2 * n), max_l1])
        return np.zeros(2 * n), Q_, C_, d_, A_, b_, eps


def optimize_focality(l, Q, target_mean, max_total_current=None,
                      max_el_current=None, Qin=None, max_angle=None,
                      max_active_electrodes=None, log_level=20):
    ''' Optimizes the focality, with the specifield mean field in the target area

    Parameters
    -----------
        l: np.ndarray
            Linear objective, obtained from target_matrices
        Q: np.ndarray
            Quadratic objective, obtained from energy_matrix
        target_mean: float or Nx1 array of floats
            Mean field component that we will try to reach (l^t currents = target_mean)
        max_total_current: float (optional)
            Maximal current flow though all electrodes. Default: No maximum
        max_el_current: float (optional)
            Maximal current flow though each electrode. Default: No maximum
        Qin: np.ndarray (optional)
            Q_in matrix from target_matrices. Needed when constraining angle
        max_angle: float (optional)
            Maximum angle between the average taget component and the average field in
            the region, in degrees. Default: No Maximum
        max_active_electrodes: int (optional)
            Maximum number of active electrodes
    Returns
    ----------
        x: np.ndarray
            Optimal electrode currents
    '''
    if max_total_current is None and max_el_current is None:
        raise ValueError('Please define a maximal total current or maximal electrode ' +
                         'current')

    if max_angle is not None and Qin is None:
        raise ValueError('When setting a max angle, a Qin must also be given')

    if max_el_current is not None:
        eps = 1e-3 * max_el_current

    elif max_total_current is not None:
        eps = 1e-3 * max_total_current

    else:
        eps = 1e-5

    l = np.copy(np.atleast_2d(np.array(l)))
    assert np.all(np.array(target_mean) > 0), \
        'The target mean values have to be positive'

    logger.log(log_level, 'Began to run optimization')
    if max_angle is None and max_active_electrodes is None:
        target_mean = np.atleast_1d(np.array(target_mean))
        assert l.shape[0] == target_mean.shape[0], \
            "Please specify one target mean per target"
        n = l.shape[1]
        l_avg, A_ub, b_ub, A_eq, b_eq, bounds = \
            _lp_variables(l, target_mean, max_total_current, max_el_current)
        sol = scipy.optimize.linprog(-l_avg, A_ub, b_ub, A_eq, b_eq,
                                     bounds=bounds)
        x_ = sol.x
        f = l.dot(x_[:n] - x_[n:])
        if np.any(np.abs(f - target_mean) >= np.abs(1e-2 * target_mean)):
            logger.log(log_level, 'Could not reach target intensities')
            return x_[:n] - x_[n:]

        else:
            logger.log(log_level, 'Target intensity reached, optimizing focality')
            # Do the QP
            a_, Q_, C_, d_, A_, b_, eps = \
                _qp_variables(l, Q, f, max_total_current, max_el_current)
            x_ = _active_set_QP(a_,
                                Q_, C_, d_, x_,
                                A=A_, b=b_, eps=eps)
            x = x_[:n] - x_[n:]
            return x

    elif max_angle is not None:
        assert l.shape[0] == 1,\
                'Angle constraint only avaliable for single targets'
        max_angle = np.deg2rad(max_angle)
        x = _constrained_angle(l, Q, target_mean,
                               max_total_current,
                               max_el_current,
                               Qin, max_angle,
                               max_active_electrodes=max_active_electrodes,
                               eps=eps,
                               log_level=log_level)
        return x

    elif max_active_electrodes is not None:
        target_mean = np.atleast_1d(np.array(target_mean))
        assert l.shape[0] == target_mean.shape[0], \
            "Please specify one target mean per target"
        n = l.shape[1]
        x = _constrained_l0(
            l, Q, target_mean,
            max_total_current,
            max_el_current,
            max_active_electrodes,
            eps=eps,
            log_level=log_level)
        return x

    else:
        raise NotImplementedError(
            'Cant handle angle constraints and maximum number of '
            'electrodes simultaneouly')


def optimize_focality_cvxpy(l, Q, target_mean, max_total_current=None,
                    max_el_current=None, Qin=None, max_angle=None, return_duals=False,
                    none_on_infeasibility=False):
    import cvxpy

    if max_total_current is None and max_el_current is None:
        raise ValueError('Please define a maximal total current or maximal electrode ' +
                         'current')
    n = l.shape[0]
    C = np.vstack([np.eye(n), -np.eye(n)])
    A = np.ones(n)

    if max_el_current is not None:
        d = max_el_current * np.ones(C.shape[0])
        eps = 1e-3 * max_el_current
    else:
        d = 1e10 * np.ones(C.shape[0])

    if max_total_current is not None:
        max_l1 = 2 * max_total_current
    else:
        max_l1 = None

    C = np.vstack([C, l])
    d = np.hstack([d, target_mean])
    x = cvxpy.Variable(n)
    p = cvxpy.Problem(cvxpy.Maximize(l * x))
    p.constraints = [C * x <= d,
                     A * x == 0]
    if max_total_current is not None:
        p.constraints.append(cvxpy.norm(x, 1) <= max_l1)
    p.solve(solver=cvxpy.SCS)
    v = np.squeeze(np.array(x.value)).T

    field_component = l.dot(v)
    if field_component * (1 + 1e-4) < target_mean:
        return v
    else:
        target_field = target_mean
        # Solve the QP
        eq_constraints = np.vstack([A, l])
        b = np.array([0, target_field])
        C = C[:-1]
        d = d[:-1]
        p = cvxpy.Problem(cvxpy.Minimize(cvxpy.quad_form(x, Q)))
        p.constraints = [C * x <= d,
                         eq_constraints * x == b]
        if max_total_current is not None:
            p.constraints.append(cvxpy.norm(x, 1) <= max_l1)
        p.solve(solver=cvxpy.SCS)
        v = np.squeeze(np.array(x.value)).T

        return v


def _eq_constrained_QP(l, Q, A, b):
    Q_qr, _ = np.linalg.qr(A.T, 'complete')
    m, n = A.shape
    Y = Q_qr[:, np.arange(m)]
    Z = Q_qr[:, np.arange(m, n)]
    w_y = np.linalg.solve(A.dot(Y), b)
    p = Y.dot(w_y)
    w_z = np.linalg.solve(Z.T.dot(Q).dot(Z), -Z.T.dot(l + Q.dot(p)))
    x = p + Z.dot(w_z)
    x = np.squeeze(x.T)
    return x


class bb_state(object):
    ''' State for branch and bound algorithm. Contains a list of active, inactive and
    unassigned electrodes '''
    def __init__(self, active, inactive, unassigned):
        self.active = active
        self.inactive = inactive
        self.unassigned = unassigned
        self.x_lb = None
        self.x_ub = None

    def inactivate(self, i):
        if i not in self.unassigned:
            raise ValueError('Can only inactivate unassigned element')

        active = copy.copy(self.active)
        inactive = copy.copy(self.inactive)
        unassigned = copy.copy(self.unassigned)

        unassigned.remove(i)
        inactive.append(i)
        return bb_state(active, inactive, unassigned)

    def activate(self, i):
        if i not in self.unassigned:
            raise ValueError('Can only activate unassigned element')

        active = copy.copy(self.active)
        inactive = copy.copy(self.inactive)
        unassigned = copy.copy(self.unassigned)

        unassigned.remove(i)
        active.append(i)
        return bb_state(active, inactive, unassigned)

class bb_node(object):
    ''' Node for branch and bound algorithm.
    Contains the current state
    bounds_funct is a funtiom wich takes in a state and return the upper bound, lower
    bound, children1 and children2 '''
    def __init__(self, state, bounds_func):
        self.state = state
        self.bounds_func = bounds_func
        self.ub_val, self.lb_val, self.child1, self.child2 = self.bounds_func(self.state)

    def split(self):
        ''' Returns 2 child nodes '''
        return bb_node(self.child1, self.bounds_func), bb_node(self.child2, self.bounds_func)

def _branch_and_bound(init, function, eps, max_k, log_level=20):
    '''Brach and Bound Algorithm
    Parameters:
    --------
    init: object
        initial state
    function: func
        Function which takes the state and returns an upper bound, lower bound and 2
        child states
    eps: float
        Tolerance between upper and lower bound
    max_k: int
        Maximum depth
    '''
    active_nodes = [bb_node(init, function)]
    k = 0
    return_val = None
    while True:
        lb = np.array([n.lb_val for n in active_nodes])
        ub = np.array([n.ub_val for n in active_nodes])
        # Prune
        keep = lb <= ub.min()
        keep[ub.argmin()] = True
        lb = lb[keep]
        ub = ub[keep]
        active_nodes = [n for i, n in enumerate(active_nodes) if keep[i]]
        logger.log(log_level,
                   "{0} Upper Bound: {1}, Lower Bound: {2}".format(
                    k, ub.min(), lb.min()))
        if ub.min() - lb.min() <= eps * np.abs(lb.min()) or k >= max_k:
            if ub.min() - lb.min() <= eps * np.abs(lb.min()):
                logger.log(log_level, 'Tolerance reached, returning')
            else:
                logger.log(log_level, 'Maximum number of iterations reached, retunning')
            return_val = active_nodes[ub.argmin()].state
            break
        q = active_nodes.pop(lb.argmin())
        c1, c2 = q.split()
        active_nodes.append(c1)
        active_nodes.append(c2)
        k += 1

    return return_val

def _bb_lower_bound(l, Q, target_mean, max_total_current, max_el_current,
                    max_l0, state):

    n = l.shape[1]
    if len(state.active) > max_l0:
        return np.zeros(n), 1e20
    # Do the optimization without the inactive electrodes

    l_avg, A_ub, b_ub, A_eq, b_eq, bounds = \
        _lp_variables(l, target_mean, max_total_current, max_el_current)

    # Append the relaxed l0 constraint if at least one electrode is active
    # We don't do it if there's no active electrodes, because that is the same as the l1
    # constraint, and having both of them can cause numerical instabilities
    if len(state.active) > 0 and len(state.active) < max_l0:
        v = np.ones(n) / max_el_current
        v[state.active] = 0
        A_ub = np.vstack([A_ub, np.tile(v, 2)])
        b_ub = np.hstack([b_ub, max_l0 - len(state.active)])

    # de-select inactive electrodes
    ac = np.zeros(n, dtype=bool)
    if len(state.active) < max_l0:
        ac[state.active + state.unassigned] = True
        n_ac = len(state.active) + len(state.unassigned)
    elif len(state.active) == max_l0:
        ac[state.active] = True
        n_ac = len(state.active)

    ac_ = np.tile(ac, 2)
    l_avg = l_avg[ac_]
    A_ub = A_ub[:, ac_]
    A_eq = A_eq[:, ac_]

    sol = scipy.optimize.linprog(-l_avg, A_ub, b_ub, A_eq, b_eq,
                                 bounds=bounds)
    x_active = sol.x[:n_ac] - sol.x[n_ac:]
    x = np.zeros(n)
    x[ac] = x_active
    f = l.dot(x)
    if np.any(np.abs(f - target_mean) >= np.abs(1e-2 * target_mean)):
        return x, 1e10

    else:
        # Do the QP
        a_, Q_, C_, d_, A_, b_, eps = \
            _qp_variables(l, Q, f, max_total_current, max_el_current)

        # Add relaxed l0 constraint
        if len(state.active) > 0 and len(state.active) < max_l0:
            v = np.ones(n) / max_el_current
            v[state.active] = 0
            C_ = np.vstack([C_, np.tile(v, 2)])
            d_ = np.hstack([d_, max_l0 - len(state.active)])

        # de-select inactive electrodes
        a_ = a_[ac_]
        Q_ = Q_[np.ix_(ac_, ac_)]
        C_ = C_[:, ac_]
        # Remove the rows of zeros
        rows_to_keep = np.linalg.norm(C_, axis=1) > eps
        C_ = C_[rows_to_keep]
        d_ = d_[rows_to_keep]
        A_ = A_[:, ac_]
        x_ = sol.x

        x_ = _active_set_QP(a_,
                            Q_, C_, d_, x_,
                            A=A_, b=b_, eps=eps)
        x = np.zeros(n)
        x[ac] = x_[:n_ac] - x_[n_ac:]

        return x, x.dot(Q).dot(x)

def _bb_upper_bound(l, Q, target_mean, max_total_current, max_el_current,
                    max_l0, state):
    n = l.shape[1]
    if len(state.active) > max_l0:
        return np.zeros(n), 1e20
    # we can't use the lower bound x because it is biased
    ac = state.active + state.unassigned
    x_ac = optimize_focality(l[:, ac], Q[np.ix_(ac, ac)],
                             target_mean, max_total_current,
                             max_el_current, log_level=10)
    x = np.zeros(n)
    x[ac] = x_ac
    # Select the "l0 - active" largest unassigned electrodes
    order_unasigned = np.argsort(-np.abs(x[state.unassigned]))
    selected_unasigned = \
        [state.unassigned[i] for i in order_unasigned[:max_l0-len(state.active)]]
    # Select the active electrodes plus the largest unassigned electrodes
    s = state.active + selected_unasigned
    x_ = optimize_focality(l[:, s], Q[np.ix_(s, s)],
                           target_mean, max_total_current,
                           max_el_current, log_level=10)
    x_ub = np.zeros(n)
    x_ub[s] = x_

    if np.any(np.abs(l.dot(x_ub) - target_mean) >= np.abs(1e-2 * target_mean)):
        return x_ub, 1e10

    return x_ub, x_ub.dot(Q).dot(x_ub)

def _bb_split(x, state):
    ''' Spliting rule for branch and bound '''
    # Split by activating / deactivating the unassigned electrode with the most current
    split_var = state.unassigned[np.argmax(np.abs(x[state.unassigned]))]
    child1 = state.activate(split_var)
    child2 = state.inactivate(split_var)
    return child1, child2

def _bb_bounds_function(l, Q, target_mean,
                        max_total_current,
                        max_el_current,
                        max_l0,
                        state):

    x_lb, lb = \
        _bb_lower_bound(l, Q, target_mean, max_total_current, max_el_current,
                        max_l0, state)
    state.x_lb = x_lb
    x_ub, ub = \
        _bb_upper_bound(l, Q, target_mean, max_total_current, max_el_current,
                        max_l0, state)
    state.x_ub = x_ub
    child1, child2 = _bb_split(x_ub, state)
    return ub, lb, child1, child2

def _constrained_l0_branch_and_bound(
        l, Q, target_mean, max_total_current,
        max_el_current, max_l0, eps_bb=1e-1,
        max_bb_iter=100, start_inactive=[], start_active=[],
        log_level=20):
    ''' Solves the constrained L0 problem using a branch-and-bound algorithm '''
    logger.log(log_level, "Starting BB")
    max_l0 = int(max_l0)
    l = copy.copy(np.atleast_2d(l))
    if max_total_current is None and max_el_current is None:
        raise ValueError('at least one of max_l1 or max_el_current must not be None')
    if max_el_current is not None:
        if max_total_current is not None:
            max_total_current = min(max_total_current, max_l0 * max_el_current / 2.)
        else:
            max_total_current = max_l0 * max_el_current
    else:
        max_el_current = max_total_current

    n = l.shape[1]
    # Define the initial state
    if len(start_inactive) > 0:
        assert min(start_inactive) >= 0 and max(start_inactive) < n, \
            'Invalid electrode number in start_inactive'
    if len(start_active) > 0:
        assert min(start_active) >= 0 and max(start_active) < n, \
            'Invalid electrode number in start_active'
    unassigned = [i for i in range(n) if i not in start_inactive + start_active]
    init = bb_state(start_active, start_inactive, unassigned)
    # partilize _bb_bounds_function
    bf = functools.partial(
        _bb_bounds_function, l, Q, target_mean, max_total_current, max_el_current, max_l0)
    # do the branch and bound
    final_state = _branch_and_bound(init, bf, eps_bb, max_bb_iter, log_level=log_level)
    x = final_state.x_ub
    return x


def _constrained_eigenvalue(Q, c):
    '''
    Finds the stationary values for
    maximize x^TQx
    such that x^Tx = 1
              c^Tx = 0
    Q is real and symmetric

    Golub, Gene H. "Some modified matrix eigenvalue problems." Siam Review 15.2 (1973) 318-334.
    https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.454.9868&rep=rep1&type=pdf
    '''
    n = Q.shape[0]
    c = c/np.linalg.norm(c)
    P = np.eye(n) - np.outer(c, c)
    K = P.dot(Q).dot(P)
    eigvals, z = np.linalg.eigh(K)
    eigvec = P.dot(z)
    return eigvals, eigvec



def _maximize_quadratic(Q, bounds=None, A_ub=None, b_ub=None, A_eq=None, b_eq=None,
                        n_x0=None, rel_tol=1e-1):

    ''' Maximizes the quadratic functions
    
    maximize  x^T Q x
    such that A_ub <= b_ub
              A_eq == b_eq

    Kind of assumes that A_eq == c^T and b_eq == 0
    uses Convex-Concave Programing. Will use n_x0 initial points, drawn from the
    eigenvalues of Q

    Returns
    ---------
    maxima_vals: np.ndarray of floats
        Values of the objective at each initial point, sortex
    maxima_x: np.ndarray of size (n_x0, len(Q))
        The coresponding x value for each maxima,

    Lipp, Thomas, and Stephen Boyd. "Variations and extension of the convex–concave
    procedure." Optimization and Engineering 17.2 (2016): 263-287.
    '''
    # Calculate contrained eigenvatues to use as initial guesses for the jacobian
    if A_eq is not None:
        _, eigvec = _constrained_eigenvalue(Q, A_eq[0])
    else:
        _, eigvec = np.linalg.eigh(Q)
    n = Q.shape[0]

    # Number of starting points
    if n_x0 is None:
        n_x0 = n-1


    maxima_x = []
    maxima_vals = []
    for i in range(n_x0):
        l = Q.dot(eigvec[:, -1-i])
        last_min = 1e20
        # Run the CCP
        while True:
            sol = scipy.optimize.linprog(
                -l, A_ub, b_ub, A_eq, b_eq,
                bounds=bounds
            )
            local_val = sol.x.dot(Q.dot(sol.x))
            if (last_min - local_val)/local_val < rel_tol:
                break
            else:
                last_min = local_val
                l = sol.x.dot(Q)

        maxima_vals.append(local_val)
        maxima_x.append(sol.x)

    maxima_vals = np.array(maxima_vals)
    maxima_x = np.array(maxima_x)

    order = maxima_vals.argsort()[::-1]
    maxima_vals = maxima_vals[order]
    maxima_x = maxima_x[order]

    return maxima_vals, maxima_x


def optimize_norm(Qin, Q, target_mean, max_total_current=None,
                  max_el_current=None, max_active_electrodes=None,
                  n_init=None, rel_tol=1e-1, log_level=20):
    # TODO: First try to maximize the quadratic to evaluate feasibility and get initial
    # points

    # Store the current best values
    global_x = None
    global_min = np.inf
    for i in range(n_init):
        cur_x = eigvec[:, -1-i]
        last_min = 1e20
        # Run the CCP
        while True:
            # Project the current solution to the boundary
            cur_energy_target = cur_x.T.dot(Qin.dot(cur_x))
            cur_x *= np.sqrt(target_mean**2/cur_energy_target)
            l = cur_x.T.dot(Qin)
            # Not sure if this will work, might need to add as inequality constraint
            x = optimize_focality(
                l, Q, target_mean**2, max_total_current=max_total_current,
                max_el_current=max_el_current,
                max_active_electrodes=max_active_electrodes,
                log_level=20
            )
            energy_total = x.T.dot(Q.dot(x))
            energy_target = x.T.dot(Qin.dot(x))
            # infeasible
            if energy_target < target_mean**2:
                # Not sure if I should stop here or if it can become feasible after some
                # iterations
                break
            # convergence
            elif (last_min - energy_total)/energy_total < rel_tol:
                break
            else:
                last_min = energy_total
                l = x.dot(Qin)

        if energy_total < global_min:
            global_min = energy_total
            global_x = x

    return global_x


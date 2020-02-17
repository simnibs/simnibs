# -*- coding: utf-8 -*-\
import copy
import functools
import numpy as np
import scipy.optimize
import scipy.linalg

from simnibs.utils.simnibs_logger import logger


class TESConstraints:
    def __init__(self, n, max_total_current, max_el_current):
        self.max_total_current = max_total_current
        self.max_el_current = max_el_current
        self.n = n

    def _bound_contraints(self):
        # The bound constraints are defined in the extended system
        # format C_.dot(x_) < d_
        C_ = np.vstack([np.eye(2 * self.n), -np.eye(2 * self.n)])
        d_ = np.hstack([self.max_el_current * np.ones(2 * self.n),
                        np.zeros(2*self.n)])
        return C_, d_

    def _l1_constraint(self):
        # The L1 constraints are defined in the extended system
        # format C_.dot(x_) < d_
        C_ = np.ones((1, 2 * self.n))
        d_ = np.array(2 * self.max_total_current)
        return C_, d_

    def _kirschoff_constraint(self):
        # Kirschoff law, defined in the extended system
        # A_.dot(x_) = b_
        A_ = np.hstack([np.ones((1, self.n)), -np.ones((1, self.n))])
        b_ = np.array(0.)
        return A_, b_

class TESOptimizationProblem(TESConstraints):
    ''' Base class for TES Optimization problems. Provide basic functionality for
    manipulating leadfields and the quadratic portion of the objective

    Parameters
    -------------
    leadfield: N_elec x N_roi x N_comp ndarray
        Leadfield

    max_total_current: float
        Maximum total current flow through all electrodes

    max_el_current: float
        Maximum current flow through each electrode
    '''
    def __init__(self, leadfield, max_total_current=1e4, max_el_current=1e4, weights=None):
        super().__init__(leadfield.shape[0] + 1, max_total_current, max_el_current)
        self.leadfield = leadfield

        if weights is None:
            self.weights = np.ones(leadfield.shape[1])
        else:
            self.weights = weights

        if len(self.weights) != leadfield.shape[1]:
            raise ValueError('Define a weight per leadfield element')

        self.Q = self._quadratic_component()

    def _quadratic_component(self):
        ''' Calculate the energy matrix for optimization
        x.dot(Q.dot(x)) = e
        "x" are the electrode currents
        "e" is the average squared electric field norm

        Returns
        ----------
        Q: np.ndarray
            Quadratic component
        '''
        lf = self.leadfield
        Q = sum(lf[..., i].dot((lf[..., i] * self.weights).T)
                for i in range(lf.shape[2]))
        Q /= np.sum(self.weights)

        P = np.linalg.pinv(np.vstack([-np.ones(Q.shape[0]), np.eye(Q.shape[0])]))

        Q = P.T.dot(Q).dot(P)
        return Q

    def extend_currents(self, x):
        '''
         Returns the extended version of the input currents x
         x_ext = [x, -x]
         with x_ext > 0
         This representation in very usefull when dealing with L1 constraints

        Parameters
        ------------
        x: ndarray
            Variable to be extended (normally currents)

        Returns
        ----------
        x_: ndarray
            Extended version
        '''

        x_ = np.hstack([x, -x])
        x_[x_ < 0] = 0
        return x_


    def solve(self):
        raise NotImplementedError()



class TESLinearConstrained(TESOptimizationProblem):
    ''' Class for solving the TES Problem with linear constraints

    This corresponds to Problem 8 in Saturnino et al., 2019
    '''
    def __init__(self, leadfield, max_total_current=1e5,
                 max_el_current=1e5, weights=None):

        super().__init__(leadfield, max_total_current, max_el_current, weights)
        self.l = np.empty((0, self.n), dtype=float)
        self.target_means = np.empty(0, dtype=float)


    def add_linear_constraint(self, target_indices, target_direction, target_mean,
                              target_weights=None):
        ''' Add a linear constrait to the problem

        Parameters
        ------------
        target_indices: list of ints
            Indices of targets.
        target_direction: ndarray
            The electric field direction to be optimized for each target position
        target_mean: float
            Target mean electric field in region
        target_weights: ndarray (optional)
            Weights (such are areas/volumes) for calculating the mean. Defined for every
            index
        '''
        if target_weights is None:
            target_weights = self.weights
        l = _calc_l(self.leadfield, target_indices, target_direction, target_weights)
        l *= np.sign(target_mean)
        self.l = np.vstack([self.l, l])
        self.target_means = np.hstack([self.target_means, np.abs(target_mean)])


    def solve(self, log_level=20):
        ''' Solves the optimization problem

        Returns
        ----------
        x: np.array
            Optimized currents
        '''
        return _linear_constrained_tes_opt(
            self.l, self.target_means, self.Q,
            self.max_el_current, self.max_total_current,
            log_level=log_level
        )

    #TODO: Add a "solve_nr_contraint" method to solve with l0 contraint


class TESLinearAngleConstrained(TESOptimizationProblem):
    ''' Class for solving the TES Problem with linear and angle contraints

    This corresponds to Problem 6 in Saturnino et al., 2019
    '''
    def __init__(self, target_indices, target_direction, target_mean, max_angle,
                 leadfield, max_total_current=1e5,
                 max_el_current=1e5, weights=None, target_weights=None):

        super().__init__(leadfield, max_total_current, max_el_current, weights)
        if target_weights is None:
            target_weights = self.weights

        self.l = np.atleast_2d(
            _calc_l(leadfield, target_indices, target_direction, target_weights)
        )
        self.Qnorm = _calc_Qnorm(leadfield, target_indices, self.weights)
        self.target_mean = np.atleast_1d(target_mean)
        self.max_angle = max_angle

    def solve(self, log_level=20):

        return _linear_angle_constrained_tes_opt(
            self.l, self.target_mean, self.Q,
            self.max_total_current,
            self.max_el_current,
            self.Qnorm, self.max_angle,
            log_level=log_level
        )


class TESLinearElecConstrained(TESLinearConstrained):
    ''' Class for solving the TES Problem with linear and number of electrodes
    constraints

    This corresponds to Problem 10 in Saturnino et al., 2019
    '''
    def __init__(self, n_elec, leadfield,
                 max_total_current=1e5,
                 max_el_current=1e5, weights=None):

        super().__init__(leadfield, max_total_current, max_el_current, weights)
        self.n_elec = n_elec

    def _solve_reduced(self, linear, quadratic, extra_ineq=None):
        l = linear[0]
        Q = quadratic[0]
        x = _linear_constrained_tes_opt(
            l, self.target_means, Q,
            self.max_el_current, self.max_total_current,
            extra_ineq=extra_ineq,
            log_level=10
        )
        if np.any(l.dot(x) < self.target_means*0.99):
            return x, 1e20
        else:
            return x, x.dot(Q).dot(x)

    def solve(self, log_level=20, eps_bb=1e-1, max_bb_iter=100, init_startegy='compact'):
        # Heuristically eliminate electrodes
        max_el_current = min(self.max_el_current, self.max_total_current)
        el = np.arange(self.n)
        if init_startegy == 'compact':
            x = _linear_constrained_tes_opt(
                self.l, self.target_means, self.Q,
                self.max_el_current, self.max_total_current,
                log_level=10
            )
            active = np.abs(x) > 1e-3 * max_el_current
            init = bb_state([], el[~active].tolist(), el[active].tolist())

        elif init_startegy == 'full':
            init = bb_state([], [], el.tolist())

        else:
            raise ValueError('Invalid initialization strategy')

        bounds_function = functools.partial(
             _bb_bounds_tes_problem,
              max_l0=self.n_elec,
             linear=[self.l],
             quadratic=[self.Q],
             max_el_current=max_el_current,
             func=self._solve_reduced
         )

        final_state = _branch_and_bound(
            init, bounds_function,
            eps_bb, max_bb_iter,
            log_level=log_level
        )

        return final_state.x_ub


class TESLinearAngleElecConstrained(TESLinearAngleConstrained):
    ''' Class for solving the TES Problem with linear, angle and number of electrodes
    constraints

    This corresponds to Problem 7 in Saturnino et al., 2019
    '''
    def __init__(self, n_elec,
                 target_indices, target_direction,
                 target_mean, max_angle,
                 leadfield, max_total_current=1e5,
                 max_el_current=1e5, weights=None,
                 target_weights=None):

        super().__init__(
            target_indices, target_direction,
            target_mean, max_angle,
            leadfield, max_total_current,
            max_el_current, weights, target_weights)

        self.n_elec = n_elec
        self._feasible = True

    def _solve_reduced(self, linear, quadratic, extra_ineq=None, feasible=False):
        l = linear[0]
        Q = quadratic[0]
        Qnorm = quadratic[1]

        x = _linear_angle_constrained_tes_opt(
            l, self.target_mean, Q,
            self.max_el_current,
            self.max_total_current,
            Qnorm, self.max_angle,
            extra_ineq=extra_ineq,
            log_level=10,
        )

        field = l.dot(x)
        if not self._feasible:
            return x, -field

        elif np.any(field < self.target_mean*0.99):
            return x, 1e20

        else:
            return x, x.dot(Q).dot(x)

    def solve(self, log_level=20, eps_bb=1e-1, max_bb_iter=100, init_startegy='compact'):
        # Heuristically eliminate electrodes
        max_el_current = min(self.max_el_current, self.max_total_current)
        el = np.arange(self.n)
        #  first determine if ploblem is feasible
        x = _linear_angle_constrained_tes_opt(
            self.l, self.target_mean, self.Q,
            self.max_el_current, self.max_total_current,
            self.Qnorm, self.max_angle,
            log_level=10
        )

        feasible = np.allclose(self.l.dot(x), self.target_mean, rtol=1e-2)
        if not feasible:
            init_startegy = 'full'

        if init_startegy == 'compact':
            active = np.abs(x) > 1e-3 * max_el_current
            init = bb_state([], el[~active].tolist(), el[active].tolist())

        elif init_startegy == 'full':
            init = bb_state([], [], el.tolist())

        else:
            raise ValueError('Invalid initialization strategy')

        bounds_function = functools.partial(
             _bb_bounds_tes_problem,
             max_l0=self.n_elec,
             linear=[self.l],
             quadratic=[self.Q, self.Qnorm],
             max_el_current=max_el_current,
             func=self._solve_reduced
         )
        self._feasible = feasible
        final_state = _branch_and_bound(
            init, bounds_function,
            eps_bb, max_bb_iter,
            log_level=log_level
        )

        return final_state.x_ub

def _calc_l(leadfield, target_indices, target_direction, weights):
    ''' Calculates the matrix "l" (eq. 14 in Saturnino et al. 2019)
    '''
    target_indices = np.atleast_1d(target_indices)
    target_direction = np.atleast_2d(target_direction)
    if target_direction.shape[1] != 3:
        target_direction = target_direction.T
    if len(target_indices) != len(target_direction):
        raise ValueError('Please define one direction per target')
    if target_direction.shape[1] != 3:
        raise ValueError('A direction must have 3 dimentions')

    target_direction = target_direction/\
        np.linalg.norm(target_direction, axis=1)[:, None]

    lf_t = leadfield[:, target_indices]
    w_idx = weights[target_indices]
    lf_t *= w_idx[None, :, None]
    lf_t /= np.sum(w_idx)
    l = np.einsum('ijk, jk -> i', lf_t, target_direction)

    P = np.linalg.pinv(
        np.vstack([-np.ones(len(l)), np.eye(len(l))])
    )

    return l.dot(P)


def _calc_Qnorm(leadfield, target_indices, weights):
    ''' Calculates the matrix "Qnorm" (like eq. 21 in Saturnino et al. 2019,
    but for all field components and not just)
    '''
    n = leadfield.shape[0]
    target_indices = np.atleast_1d(target_indices)
    lf_t = leadfield[:, target_indices]
    w_idx = weights[target_indices]
    Q_in = sum(
        lf_t[..., i].dot((lf_t[..., i]*w_idx[None, :]).T) for i in range(3))
    Q_in /= np.sum(w_idx)

    P = np.linalg.pinv(
        np.vstack([-np.ones(n), np.eye(n)])
    )

    Qnorm = P.T.dot(Q_in).dot(P)
    return Qnorm


def _linear_constrained_tes_opt(l, target_mean, Q,
                                max_el_current, max_total_current,
                                extra_ineq=None, extra_eq=None,
                                log_level=10):

        assert l.shape[0] == target_mean.shape[0], \
            "Please specify one target mean per target"
        assert l.shape[1] == Q.shape[0]


        n = l.shape[1]
        tes_constraints = TESConstraints(n, max_total_current, max_el_current)
        # First solve an LP to get a feasible starting point
        l_ = np.hstack([l, -l])

        # L1 constaints
        C_, d_ = tes_constraints._l1_constraint()
        # Kirschoffs law
        A_, b_ = tes_constraints._kirschoff_constraint()

        # Inequality constraints
        if extra_ineq is not None:
            C_ = np.vstack([C_, extra_ineq[0]])
            d_ = np.hstack([d_, extra_ineq[1]])

        if extra_eq is not None:
            A_ = np.vstack([A_, extra_eq[0]])
            b_ = np.hstack([b_, extra_eq[1]])

        sol = scipy.optimize.linprog(
            -np.average(l_, axis=0),
            np.vstack([C_, l_]), np.hstack([d_, target_mean]),
            A_, b_,
            bounds=(0, max_el_current)
        )
        x_ = sol.x

        # Test if the objective can be reached
        f = l.dot(x_[:n] - x_[n:])

        if np.any(np.abs(f - target_mean) >= np.abs(1e-2 * target_mean)):
            logger.log(log_level, 'Could not reach target intensities')
            return x_[:n] - x_[n:]

        logger.log(log_level, 'Target intensity reached, optimizing focality')

        # Do the QP
        eps = 1e-3*np.max(np.abs(x_))
        Q_ = np.vstack([np.hstack([Q, -Q]),
                        np.hstack([-Q, Q])])
        C_b, d_b = tes_constraints._bound_contraints()

        x_ = _active_set_QP(
            np.zeros(2*n), Q_,
            np.vstack([C_b, C_]), np.hstack([d_b, d_]),
            x_, eps,
            np.vstack([A_, l_]), np.hstack([b_, f])  # I use "f"
        )

        x = x_[:n] - x_[n:]
        return x


def _calc_angle(x, Qin, l):
    tan = np.sqrt(np.abs(x.dot(Qin).dot(x) - l.dot(x) ** 2))
    return np.abs(np.arctan2(tan, l.dot(x)))[0]


def _linear_angle_constrained_tes_opt(
    l, target_mean, Q,
    max_el_current, max_total_current,
    Qin, max_angle,
    extra_ineq=None,
    extra_eq=None,
    eps_linear=1e-5,
    eps_angle=1e-1, log_level=20):

    max_angle = np.deg2rad(max_angle)
    logger.log(log_level, 'Running optimization with angle constraint')
    max_iter = 20
    # Try to find where we can find values below and above the target
    it = 0
    above = 0

    # Check if the maximal focality solution alreay fulfills the constraint
    x = _linear_constrained_tes_opt(
        l, target_mean, Q,
        max_el_current, max_total_current,
        extra_ineq=extra_ineq, extra_eq=extra_eq,
        log_level=log_level-10
    )

    if _calc_angle(x, Qin, l) <= max_angle:
        logger.log(log_level, 'Max focality solution fullfills angle constraint')
        return x

    x_above = np.copy(x)

    # calculate the smallest angle, given a fixed intensity
    def _minimize_angle(alpha):
        x_l = _linear_constrained_tes_opt(
            l, alpha*target_mean, Qin,
            max_el_current, max_total_current,
            extra_ineq=extra_ineq, extra_eq=extra_eq,
            log_level=log_level-10
        )
        return x_l

    x = _minimize_angle(1.)
    angle = _calc_angle(x, Qin, l)

    # if we cannot reduce the angle to the target while keeting l^t x at the target
    # intensity, reduce the target intensity untill it's achievable.
    if angle > max_angle:
        logger.log(log_level, "Target intensity can't be reached, reducing it")
        above = 1.
        angle_above = angle
        below = 0
        angle_below = 0
        it = 0
        # Use the secant method
        while not (angle > max_angle * (1 - eps_angle) and angle < max_angle):
            alpha = above + \
                (max_angle * (1 - eps_angle * .5) - angle_above) * \
                (below - above)/(angle_below - angle_above)
            x = _minimize_angle(alpha)
            angle = _calc_angle(x, Qin, l)
            logger.log(log_level,
                       '{0} alpha: {1:.3e}, angle: {2:.2e}, max_angle: {3:.2e}'.format(
                        it, alpha, angle, max_angle))
            if angle < max_angle:
                below = alpha
                angle_below = angle
                x_below = np.copy(x)
            else:
                above = alpha
                angle_above = angle
            it += 1
            if it > max_iter:
                if below == 0:
                    return x
                else:
                    return x_below
        return x

    # In this case, we know that by minimizing x^t Qin x while keeping l^t x = t, we can
    # achieve the bound
    # find a combination between Q and Qin that maximizes focality while keeping Qin in
    # the bound
    else:
        logger.log(
            log_level,
           "Target intensity reached, optimizing focality with angle constraint"
        )
        angle_below = angle
        below = 1.
        x_below = np.copy(x)
        angle_above = _calc_angle(x_above, Qin, l)
        above = 0

        # Start the secant method
        it = 0
        alpha = 1
        while not (angle > max_angle * (1 - eps_angle) and angle < max_angle):
            alpha = above + \
                (max_angle * (1 - eps_angle * .5) - angle_above) * \
                (below - above)/(angle_below - angle_above)

            x = _linear_constrained_tes_opt(
                l, target_mean, (1 - alpha) * Q + alpha * Qin,
                max_el_current, max_total_current,
                extra_ineq=extra_ineq, extra_eq=extra_eq,
                log_level=log_level-10
            )
            angle = _calc_angle(x, Qin, l)
            logger.log(
                log_level,
                f'{it} alpha: {alpha:.2f}, angle: {angle:.2e}, max_angle: {max_angle:.2e}'
            )

            if angle > max_angle:
                above = alpha
                angle_above = angle
            else:
                below = alpha
                angle_below = angle
                x_below = np.copy(x)

            it += 1
            if it > max_iter:
                return x_below

        return x



def _bb_bounds_tes_problem(state, max_l0, linear, quadratic, max_el_current, func):
    ''' Returns upper bound, lower bound a child states for a TES optimization problem

    This is meant ro be a general interface for the BB algorithm. It requireas a function
    "func" with the following call
        x, objective = func(linear, quadratic, extra_ineq)
    "linear" is a list of linear factors, quadratic is a list of quadratic factors
    '''
    if len(state.active) > max_l0:
        return 1e20, 1e20, None, None

    n = linear[0].shape[1]

    # Create a list of active + unasigned electrodes
    ac = np.zeros(n, dtype=bool)
    if len(state.active) < max_l0:
        ac[state.active + state.unassigned] = True
    elif len(state.active) == max_l0:
        ac[state.active] = True

    ##  Lower bound calculation

    # Create the extra l1 constraint (PS1.6)
    if len(state.active) > 0 and len(state.active) < max_l0:
        v = np.ones(n) / max_el_current
        v[state.active] = 0
        v = np.tile(v[ac], 2)
        extra_ineq = (v, max_l0 - len(state.active))
    else:
        extra_ineq = None

    # Run problem in reduced system
    linear_ac = [l[:, ac] for l in linear]
    quadratic_ac = [Q[np.ix_(ac, ac)] for Q in quadratic]
    x_ac, objective_lb = func(
        linear_ac, quadratic_ac,
        extra_ineq=extra_ineq
    )
    x_lb = np.zeros(n)
    x_lb[ac] = x_ac
    ## Upper bound calculation
    # Solve PS2
    x_ac, _ = func(
        linear_ac, quadratic_ac
    )
    x_ub1 = np.zeros(n)
    x_ub1[ac] = x_ac

    # Solve problem PS3
    # Select the "l0 - active" largest unassigned electrodes
    order_unasigned = np.argsort(-np.abs(x_ub1[state.unassigned]))
    selected_unasigned = \
        [state.unassigned[i] for i in order_unasigned[:max_l0-len(state.active)]]

    # Select the active electrodes plus the largest unassigned electrodes
    s = state.active + selected_unasigned
    linear_s = [l[:, s] for l in linear]
    quadratic_s = [Q[np.ix_(s, s)] for Q in quadratic]
    x_s, objective_ub = func(
        linear_s, quadratic_s
    )
    x_ub = np.zeros(n)
    x_ub[s] = x_s

    # Split by activating / deactivating the unassigned electrode with the most current
    split_var = state.unassigned[np.argmax(np.abs(x_ub[state.unassigned]))]
    child1 = state.activate(split_var)
    child2 = state.inactivate(split_var)
    state.x_ub = x_ub
    state.x_lb = x_lb

    return objective_ub, objective_lb, child1, child2


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
                   f"{k} Upper Bound: {ub.min():.2e}, Lower Bound: {lb.min():.2e}")
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

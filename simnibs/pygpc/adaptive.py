import multiprocessing
import functools
import numpy as np
import scipy.linalg

from .grid import randomgrid
from .ni import reg


def _expand_polynomial(active_set, old_set, to_expand, order_max, interaction_max):
    ''' Algorithm by Gerstner and Griebel '''
    active_set = [tuple(a) for a in active_set]
    old_set = [tuple(o) for o in old_set]
    to_expand = tuple(to_expand)
    active_set.remove(to_expand)
    old_set += [to_expand]
    expand = []
    for e in range(len(to_expand)):
        forward = np.asarray(to_expand, dtype=int)
        forward[e] += 1
        has_predecessors = True
        for e2 in range(len(to_expand)):
            if forward[e2] > 0:
                predecessor = forward.copy()
                predecessor[e2] -= 1
                predecessor = tuple(predecessor)
                has_predecessors *= predecessor in old_set
        if has_predecessors and np.sum(np.abs(forward)) <= order_max \
           and np.sum(forward > 0) <= interaction_max:
            expand += [tuple(forward)]
            active_set += [tuple(forward)]

    return active_set, old_set, expand


def _choose_to_expand(poly_idx, active_set, old_set, coeffs, order_max, interaction_max):
    coeffs = np.linalg.norm(coeffs, axis=1)
    for idx in poly_idx[coeffs.argsort()[::-1]]:
        if tuple(idx) in active_set:
            _, _, expand = _expand_polynomial(
                active_set, old_set, tuple(idx), order_max, interaction_max)
            if len(expand) > 0:
                return tuple(idx)

    raise ValueError('Could not find a polynomial in the expansion and in the active set')


def _relative_error(data, approx):
    if data.ndim == 1:
        data = data[None, :]
    if approx.ndim == 1:
        approx = approx[None, :]
    err = np.average(
        np.linalg.norm(data - approx, axis=1) /
        np.linalg.norm(data, axis=1))
    return err


def _tikhonov(A, b, alpha):
    '''Solves inv(A.T A + alpha*Id) * A.T b
    based on the _solve_cholesky function from sklearn
    '''
    n_samples, n_features = A.shape
    AtA = A.T.dot(A)
    Atb = A.T.dot(b)
    AtA.flat[::n_features + 1] += alpha
    return scipy.linalg.solve(AtA, Atb, sym_pos=True,
                              overwrite_a=True)


def _k_fold_cv_regression(A, data, regression_function, error_eq=_relative_error, k=10):
    if data.ndim == 1:
        data = data[:, None]
    n_simulations = data.shape[0]
    if n_simulations <= k:
        k = n_simulations
    shuffle = np.arange(n_simulations, dtype=int)
    np.random.shuffle(shuffle)
    groups = np.zeros((k, n_simulations), dtype=bool)
    group_sizes = float(n_simulations) / k
    for i in range(k):
        f = int(i * group_sizes)
        t = int((i + 1) * group_sizes)
        if i == k - 1:
            t = n_simulations
        groups[i, shuffle[f:t]] = True
    eps = 0.0
    for g in groups:
        # determine regression coefficients
        x = regression_function(A[~g, :], data[~g, :])
        eps += sum(g) * error_eq(data[g, :], A[g, :].dot(x))
    eps /= n_simulations
    return eps


class RegularizedRegression(reg):
    ''' RegularizedRegression Regression
    is a subclass of reg
    '''
    def __init__(self, pdftype, pdfshape, limits, order, order_max, interaction_order,
                 grid, regularization_factors=np.logspace(-5, 3, 9)):
        super().__init__(pdftype, pdfshape, limits, order, order_max, interaction_order, grid)
        self.regularization_factors = regularization_factors

    def construct_gpc_matrix(self):
        """ construct the gpc matrix A [N_samples x N_poly] """
        A = np.zeros([self.N_grid, self.N_poly])

        for i_poly in range(self.N_poly):
            A1 = np.ones(self.N_grid)
            for i_DIM in range(self.DIM):
                A1 *= self.poly[self.poly_idx[i_poly][i_DIM]][i_DIM](self.grid.coords_norm[:,i_DIM])
            A[:, i_poly] = A1
        return A

    def expand(self, data, return_error=True, return_reg_factor=False):
        """ Determine the gPC coefficients by the regression method
        input:    data ... results from simulations with N_out output quantities,
        np.array() [N_samples x N_out]
        output: coeffs ... gPC coefficients, np.array() [N_coeffs x N_out]
        """
        if data.ndim == 1:
            data = data[:, None]

        errors = np.zeros_like(self.regularization_factors)
        A = self.construct_gpc_matrix()
        for i, reg_factor in enumerate(self.regularization_factors): 
            regression_fun = functools.partial(_tikhonov, alpha=reg_factor)
            errors[i] = _k_fold_cv_regression(A, data, regression_fun)

        min_error = np.argmin(errors)
        reg_error = errors[min_error]
        selected_reg = self.regularization_factors[min_error]
        coeff = _tikhonov(A, data, selected_reg)
        if not return_error and not return_reg_factor:
            return coeff
        else:
            ret = [coeff]
            if return_error:
                ret += [reg_error]
            if return_reg_factor:
                ret += [selected_reg]
            return tuple(ret)


    def add_n_sampling_points(self, n):
        ''' Adds a given number of sampling points

        Parameters
        -----------
        n: float
            Will add ceil(n) sampling points

        Returns
        --------
        sampling_points: ndarray
            The new points added
        '''
        newgridpoints = randomgrid(self.pdftype, self.pdfshape, self.limits,
                                   int(np.ceil(n)))
            
        # append points to existing grid
        self.grid.coords = np.vstack([self.grid.coords, newgridpoints.coords])
        self.grid.coords_norm = np.vstack([self.grid.coords_norm, newgridpoints.coords_norm])
        self.N_grid = self.grid.coords.shape[0]

        return newgridpoints.coords


def run_reg_adaptive_grid(pdftype, pdfshape, limits, func, args=(),
                          data_poly_ratio=2, max_iter=1000,
                          order_max=None, interaction_max=None,
                          eps=1E-3, n_cpus=1, print_function=None):
    """  
    Adaptive regression approach based on leave one out cross validation error
    estimation
    
    Parameters
    ----------

    pdftype : list
              Type of probability density functions of input parameters,
              i.e. ["beta", "norm",...]
    pdfshape : list of lists
               Shape parameters of probability density functions
               s1=[...] "beta": p, "norm": mean
               s2=[...] "beta": q, "norm": std
               pdfshape = [s1,s2]
    limits : list of lists
             Upper and lower bounds of random variables (only "beta")
             a=[...] "beta": lower bound, "norm": n/a define 0
             b=[...] "beta": upper bound, "norm": n/a define 0
             limits = [a,b]
    func : callable func(x,*args)
           The objective function to be minimized.
    args : tuple, optional
           Extra arguments passed to func, i.e. f(x,*args).         
    data_poly_ratio: float, optional
        Number of function evaluations per polynomial in the basis. Default: 1.5
    max_iter: int, optional
        Maximum number of iterations. Default: 1000
    order_max: int, optional
        Maximum L1 norm for polynomials. Default: no maximum
    interaction_max: int, optional
        Maximum interaction order. Default: no maximum
    eps: float, optional
        Relative mean error of leave one out cross validation
    n_cpus: int, optional
        Number of cpus to evaluate "func". Default: 1
    print_function : function
        A function to print convergence information. Default: Does not print

    Returns
    -------
    gobj: object
       gpc object
    res: ndarray
       Function values at grid points of the N_out output variables
       size: [N_samples x N_out]        
    """
    
    # initialize iterators
    eps_gpc = eps+1
    i_samples = 0
    i_iter = 0
    DIM = len(pdftype)
    order = 0
    active_set = [tuple(0 for d in range(DIM))]
    old_set = []
    to_expand = tuple(0 for d in range(DIM))
    if order_max is None:
        order_max = 200
    if interaction_max is None:
        interaction_max = DIM

    if n_cpus > 1:
        pool = multiprocessing.Pool(processes=n_cpus)

    while i_iter < max_iter:
        i_iter = i_iter + 1
        if print_function:
            print_function("Iteration #{}".format(i_iter))
            print_function("=============")

        if i_iter == 1:
            grid_init = randomgrid(pdftype, pdfshape, limits, 0)
            regobj = RegularizedRegression(
                pdftype, pdfshape, limits, order*np.ones(DIM),
                order_max=order, interaction_order=DIM, grid=grid_init)
            RES = np.empty((0, 0))

        else:
            active_set, old_set, expand = _expand_polynomial(
                active_set, old_set, to_expand, order_max, interaction_max)
            # add polynomials to gpc expansion
            regobj.enrich_polynomial_basis(
                np.asarray(expand, dtype=int), form_A=False)

        n = int(np.ceil(data_poly_ratio * len(regobj.poly_idx) - RES.shape[0]))
        # This here will ensure we always use n_cpus
        if n % n_cpus > 0:
            n += n_cpus - (n % n_cpus)
        regobj.add_n_sampling_points(n)

        # run repeated simulations
        # MP version
        if n_cpus > 1:
            processes = []
            for s in range(i_samples, regobj.grid.coords.shape[0]):
                if print_function:
                    print_function("Performing simulation #{}".format(s+1))
                x = regobj.grid.coords[s, :]
                processes.append(pool.apply_async(func, (x, ) + args))
            if i_samples == 0:
                RES = np.vstack([p.get() for p in processes])
            else:
                RES = np.vstack([RES] + [p.get() for p in processes])
        # Sequential version
        else:
            for s in range(i_samples, regobj.grid.coords.shape[0]):
                if print_function:
                    print_function("Performing simulation #{}".format(s+1))
                # read radom variables from grid
                x = regobj.grid.coords[s, :]
                # evaluate function at grid points
                res = func(x, *(args))
                # append result to solution matrix (RHS)
                if s == 0:
                    RES = res[None, :]
                else:
                    RES = np.vstack([RES, res])

        i_samples = s + 1
        # Expand
        coeffs, eps_gpc = regobj.expand(RES)
        if print_function:
            print_function("gPC Error Estimate: {0:1e}".format(eps_gpc))
        if (eps_gpc > eps):
            to_expand = _choose_to_expand(
                regobj.poly_idx, active_set, old_set, coeffs, order_max, interaction_max)
        else:
            break

    if i_iter >= max_iter:
        raise ValueError('Maximum number of iterations reached')

    if n_cpus > 1:
        pool.close()
        pool.join()

    return regobj, RES

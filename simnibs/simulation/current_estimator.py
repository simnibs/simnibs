import pygpc
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.distance import cdist
from sklearn.linear_model import LinearRegression
from collections import OrderedDict


class CurrentEstimator():
    """
    CurrentEstimator class to track TES electrode currents for faster estimation of fake Dirichlet boundary conditions
    """
    def __init__(self, electrode_pos=None, current=None, method="linear"):
        """
        Constructor for CurrentEstimator class instance

        Parameters
        ----------
        electrode : CircularArray or ElectrodeArrayPair instance
            TES/TTF/TI Electrode
        electrode_pos : list of np.ndarray of float [n_array_free][n_pos x 3]
            Positions and orientations of ElectrodeArrayPair or CircularArray
        current : list of np.ndarray of float [n_array_free][n_pos x 3]
            Optimal currents corresponding to electrode positions
        method : str, optional, default: "linear"
            Method to estimate the electrode currents:
            - "GP": Gaussian process
            - "linear": linear regression

        Attributes
        ----------
        self.electrode : CircularArray or ElectrodeArrayPair instance
            TES/TTF/TI Electrode
        self.electrode_pos : list of np.ndarray of float [n_array_free][n_pos x 3]
            Positions and orientations of ElectrodeArrayPair or CircularArray
        self.currents : list of np.ndarray of float [n_array_free][n_pos x 3]
            Optimal currents corresponding to electrode positions
        self.method : str
            Method to estimate the electrode currents:
            - "GP": Gaussian process
        """
        self.method = method   # "GP" (Gaussian process) "linear" (sklearn.linear_model.LinearRegression) "gpc" (pygpc)

        self.gpc_grid = None
        self.gpc_coeffs = None
        self.gpc_session = None
        self.gpc_algorithm = None
        self.gpc_order_list = None
        self.gpc_parameters = None
        self.gpc_n_coeffs_list = None

        # reshape electrode_pos from list to nparray [n_train x n_ele]
        if electrode_pos is not None:
            self.electrode_pos = np.hstack(electrode_pos)

            if self.electrode_pos.ndim == 1:
                self.electrode_pos = self.electrode_pos[np.newaxis, :]
        else:
            self.electrode_pos = None

        # reshape currents from list to nparray [n_train x n_ele]
        if current is not None:
            self.current = np.hstack(current)

            if self.current.ndim == 1:
                self.current = self.current[np.newaxis, :]
        else:
            self.current = None

    def add_training_data(self, electrode_pos, current):
        """
        Add training data to the model.

        Parameters
        ----------
        electrode_pos : list of np.ndarray of float [n_array_free][n_pos x 3] or np.ndarray of float [n_pos x (n_array_free * 3)]
            Positions and orientations of ElectrodeArrayPair or CircularArray
            either a list [n_array_free] of np.ndarray [n_pos x 3]
            [ np.array([[beta_1_pos_1, lambda_1_pos_1, alpha_1_pos_1],    np.array([[beta_2_pos_1, lambda_2_pos_1, alpha_2_pos_1],
                        [beta_1_pos_2, lambda_1_pos_2, alpha_1_pos_2],              [beta_2_pos_2, lambda_2_pos_2, alpha_2_pos_2],
                                          ...                       ])                                ...                       ]) ]
                                          array 1                                                   array 2

            or a np.ndarray of [n_pos x (n_array_free * 3)]
            np.array([[beta_1_pos_1, lambda_1_pos_1, alpha_1_pos_1, beta_2_pos_1, lambda_2_pos_1, alpha_2_pos_1],
                      [beta_1_pos_2, lambda_1_pos_2, alpha_1_pos_2, beta_2_pos_2, lambda_2_pos_2, alpha_2_pos_2],
                                      ...                                               ...                    ])
                                     array 1                                           array 2
        current : np.ndarray of float [n_pos] or float
            Optimal current corresponding to electrode positions
         """

        if type(electrode_pos) in [float, int, np.float64]:
            electrode_pos = np.array([electrode_pos])

        if type(current) in [float, int, np.float64]:
            current = np.array([current])

        # reshape and append passed electrode position to set of all electrode positions [n_train x n_ele]
        if self.electrode_pos is None:
            self.electrode_pos = np.hstack(electrode_pos)

            if self.electrode_pos.ndim == 1:
                self.electrode_pos = self.electrode_pos[np.newaxis, :]

        else:
            # do not add training data, which is already included
            if np.hstack(electrode_pos).tolist() in self.electrode_pos.tolist():
                return

            self.electrode_pos = np.vstack((self.electrode_pos, np.hstack(electrode_pos)))

        # reshape and append passed electrode currents to set of all electrode currents [n_train x n_ele]
        if self.current is None:
            self.current = np.hstack(current)

            if self.current.ndim == 1:
                self.current = self.current[np.newaxis, :]

        else:
            self.current = np.vstack((self.current, np.hstack(current)))

        # determine list of approximation orders and associated number of gpc coefficients
        if self.method == "gpc":
            if self.gpc_grid is None:
                self.gpc_grid = pygpc.RandomGrid(parameters_random=self.gpc_parameters, coords=self.electrode_pos)
            else:
                if electrode_pos.ndim == 1:
                    electrode_pos = electrode_pos[np.newaxis, :]
                self.gpc_grid.extend_random_grid(coords=electrode_pos)

            if self.gpc_n_coeffs_list is None:
                self.gpc_order_list = np.arange(11) + 1
                self.gpc_n_coeffs_list = np.zeros(len(self.gpc_order_list)).astype(int)
                dim = self.electrode_pos.shape[1]

                for i_order, order in enumerate(self.gpc_order_list):
                    self.gpc_n_coeffs_list[i_order] = pygpc.get_num_coeffs_sparse(order_dim_max=[order] * dim,
                                                                                  order_glob_max=order,
                                                                                  order_inter_max=2,
                                                                                  dim=dim,
                                                                                  order_inter_current=None,
                                                                                  order_glob_max_norm=1)

    def estimate_current(self, electrode_pos):
        """
        Provides an estimate of the optimal currents at a given electrode position.

        Parameters
        ----------
        electrode_pos : list of np.ndarray of float [n_array_free][3]
            Electrode position

        Returns
        -------
        currents : ndarray of float [n_ele]
            Estimates optimal currents corresponding to electrode position
        """
        if self.current is None:
            return None

        if electrode_pos.ndim == 1:
            electrode_pos = electrode_pos[np.newaxis, :]

        if self.method == "GP":
            current = get_estimate_gaussian_process(x=electrode_pos,
                                                    x_train=self.electrode_pos,
                                                    y_train=self.current)
        elif self.method == "linear":
            if self.electrode_pos.shape[0] > 1:
                current = get_estimate_linear_regression(x=electrode_pos,
                                                         x_train=self.electrode_pos,
                                                         y_train=self.current)
            else:
                return None

        elif self.method == "gpc":
            current = self.get_estimate_gpc(x=electrode_pos,
                                            x_train=self.electrode_pos,
                                            y_train=self.current)
            if current is None:
                return None

        else:
            raise NotImplementedError(f"Specified current estimation method '{self.method}' not implemented.")

        return current[0]

    def set_gpc_parameters(self, lb, ub):
        """
        Sets gpc parameters in self.gpc_parameters

        Parameters
        ----------
        lb : np.array of float [n_para]
            Lower bounds of parameters
        ub : np.array of float [n_para]
            Upper bounds of parameters
        """
        self.gpc_parameters = OrderedDict()
        for i in range(len(lb)):
            self.gpc_parameters[str(i)] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[lb[i], ub[i]])

    def get_estimate_gpc(self, x, x_train, y_train):
        """
        Determine coordinates at highest variance determined by gpc.

        Parameters
        ----------
        x : ndarray of float [n_para]
            Query point
        x_train : ndarray of float [n_train x n_para]
            Set of training points (parameters)
        y_train : ndarray of float [n_train x n_para]
            Set of training points (solutions)

        Returns
        -------
        y : ndarray of float
            Estimate
        """
        # select order based on number of samples
        idx_order = np.where(self.gpc_grid.n_grid >= (2*self.gpc_n_coeffs_list))[0]

        if len(idx_order) > 0:
            order = self.gpc_order_list[idx_order[-1]]
        else:
            return None

        options = dict()
        options["method"] = "reg"
        options["solver"] = "Moore-Penrose"
        options["settings"] = None
        options["order"] = [order] * self.electrode_pos.shape[1]
        options["order_max"] = order
        options["interaction_order"] = 2
        options["error_type"] = None
        options["n_samples_validation"] = None
        options["fn_results"] = None
        options["save_session_format"] = ".pkl"
        options["backend"] = "omp"
        options["verbose"] = False

        # define algorithm
        self.gpc_algorithm = pygpc.Static_IO(parameters=self.gpc_parameters,
                                             options=options,
                                             grid=self.gpc_grid,
                                             results=y_train)

        # initialize gPC Session
        self.gpc_session = pygpc.Session(algorithm=self.gpc_algorithm)

        # run gPC algorithm
        self.gpc_session, self.gpc_coeffs, _ = self.gpc_session.run()

        # approximate current
        x_norm = self.gpc_grid.get_normalized_coordinates(coords=x)
        current = self.gpc_session.gpc[0].get_approximation(coeffs=self.gpc_coeffs, x=x_norm)

        return current


def get_parameters_gaussian_process(Xtrain, ytrain):
    """
    Determine optimal hyperparameters for Gaussian Process Regression (lengthscale, variance), without noise.

    Parameters
    ----------
    Xtrain : np.ndarray of float [N_train x dim]
        Coordinates of the training data
    ytrain : np.ndarray of float [N_train x n_out]
        Function values at the training data points for several output variables

    Returns
    -------
    lengthscale : ndarray of float [n_out]
        Lengthscale parameter for individual output variables
    variance : ndarray of float [n_out]
        Output variance for individual output variables
    """
    if ytrain.ndim == 1:
        ytrain = ytrain[:, np.newaxis]

    n_y = ytrain.shape[1]

    lengthscale_init = .2
    kernel_variance_init = 1
    bounds = ((1e-3, 1e2), (1e-3, 1e2))
    initial_parameters = np.array([lengthscale_init, kernel_variance_init])

    lengthscale = np.zeros(n_y)
    variance = np.zeros(n_y)

    for i in range(n_y):
        result = minimize(compute_neg_loglik, initial_parameters, (Xtrain, ytrain[:, i]),
                          method='l-bfgs-b', bounds=bounds)
        lengthscale[i] = result.x[0]
        variance[i] = result.x[1]

    return lengthscale, variance


def compute_neg_loglik(parameters, Xtrain, ytrain):
    """
    Computes the negative log likelihood of the hyperparameters of the Gaussian Process Regression.

    Parameters
    ----------
    parameters : np.ndarray of float [2]
        Hyperparameters (lengthscale, variance)
    Xtrain : np.ndarray of float [N_train x dim]
        Coordinates of the training data
    ytrain : np.ndarray of float [N_train]
        Function values at the training data points

    Returns
    -------
    log_likelihood : float
        Negative log likelihood
    """
    if Xtrain.ndim == 1:
        Xtrain = Xtrain[np.newaxis, :]

    lengthscale, variance = parameters
    K = squared_exponential_kernel(Xtrain, Xtrain, lengthscale, variance)  # n_train x n_train

    try:
        L = np.linalg.cholesky(K)
    except np.linalg.LinAlgError:
        return 0

    alpha = np.linalg.solve(L.T, np.linalg.solve(L, ytrain))
    log_likelihood = - 0.5 * ytrain.T @ alpha - np.log(np.diag(L)).sum() - len(ytrain) / 2 * np.log(2 * np.pi)

    return - log_likelihood.squeeze()


def squared_exponential_kernel(x, y, lengthscale, variance):
    """
    Computes the squared exponential kernel for Gaussian Processes.

    Parameters
    ----------
    x : np.ndarray of float [N x dim]
        Input observation locations
    y : np.ndarray of float [M x dim]
        Output observation locations
    lengthscale : float
        Lengthscale parameter
    variance : float
        Output variance

    Returns
    -------
    k : np.ndarray of float [M x X]
        Kernel function values (covariance function or covariance matrix)
    """
    sqdist = cdist(x, y, 'sqeuclidean')
    k = variance * np.exp(-0.5 * sqdist * (1/lengthscale**2))
    return k





def get_estimate_linear_regression(x, x_train, y_train):
    """
    Determine coordinates at highest variance determined by linear regression.

    Parameters
    ----------
    x : ndarray of float [n_para]
        Query point
    x_train : ndarray of float [n_train x n_para]
        Set of training points (parameters)
    y_train : ndarray of float [n_train x n_para]
        Set of training points (solutions)

    Returns
    -------
    y : ndarray of float
        Estimate
    """

    reg = LinearRegression().fit(x_train, y_train)

    # reg.score(x_train, y_train)
    # reg.coef_
    # reg.intercept_

    return reg.predict(x)


def get_estimate_gaussian_process(x, x_train, y_train, lengthscale=None, variance=None):
    """
    Determine coordinates at highest variance determined by Gaussian Process Regression

    Parameters
    ----------
    x : ndarray of float [n_para]
        Query point
    x_train : ndarray of float [n_train x n_para]
        Set of training points (parameters)
    y_train : ndarray of float [n_train x n_para]
        Set of training points (solutions)
    lengthscale : float, optional, default: 1.
        Lengthscale parameter
    variance : float, optional, default: 1.
        Output variance

    Returns
    -------
    y : ndarray of float
        Estimate
    """

    if y_train.ndim == 1:
        y_train = y_train[:, np.newaxis]

    n_y = y_train.shape[1]

    # update lengthscale and variance estimates
    if (lengthscale is None) or (variance is None):
        lengthscale, variance = get_parameters_gaussian_process(Xtrain=x_train,
                                                                ytrain=y_train)

    # estimate y
    y = np.zeros(n_y)

    for i in range(n_y):
        # kernels
        K = squared_exponential_kernel(x=x_train, y=x_train, lengthscale=lengthscale[i], variance=variance[i])  # n_train x n_train
        Ks = squared_exponential_kernel(x=x_train, y=x, lengthscale=lengthscale[i], variance=variance[i])  # n_train x n_test
        # Kss = squared_exponential_kernel(x=x, y=x, lengthscale=lengthscale, variance=variance)  # n_test x n_test

        try:
            # cholesky decomposition
            L = np.linalg.cholesky(K)

            # compute v
            # v = np.linalg.solve(L, Ks)

            # alpha
            alpha = np.linalg.solve(L.T, np.linalg.solve(L, y_train[:, i]))

        except np.linalg.LinAlgError:
            # print("Warning: Cholesky decomposition of K* matrix did not converge ... using Moore-Penrose pseudo inverse.")
            # v = np.linalg.pinv(K) @ Ks

            alpha = np.linalg.pinv(K) @ y_train[:, i]

        # compute the mean function
        y[i] = Ks.T @ alpha

    # compute the covariance
    # covariance = Kss - (v.T @ v)

    return y

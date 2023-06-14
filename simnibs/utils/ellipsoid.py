import numpy as np
import multiprocessing.pool

from _functools import partial
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, minimize
import warnings
warnings.filterwarnings("ignore")
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

from ..mesh_tools import cgal
from ..mesh_tools.surface import Surface


class Ellipsoid():
    """
    Ellipsoid class containing transformations and geodesic distance functions

    Parameters
    ----------
    center : np.ndarray of float [3], optional, default: None
        Center of ellipsoid
    radii : np.ndarray of float [3], optional, default: None
        Length of semiaxis
    rotmat : np.ndarray of float [3 x 3], optional, default: None
        Rotation matrix of ellipsoid (direction of principal axes)

    Attributes
    ----------
    E_x : float
        Linear eccentricity between first and third axis
    E_y : float
        Linear eccentricity between second and third axis
    E_e : float
        Linear eccentricity between first and second axis
    E_x2 : float
        Linear eccentricity between first and third axis (squared)
    E_y2 : float
        Linear eccentricity between second and third axis (squared)
    E_e2 : float
        Linear eccentricity between first and second axis (squared)
    e_x2 : float
        Squared eccentricity between first and third axis
    e_e2 : float
        Squared eccentricity between first and second axis
    """
    def __init__(self, center=None, radii=None, rotmat=None):
        """
        Constructor of Ellipsoid class instance
        """
        self._radii = None

        # ellipsoid parameters
        if center is not None:
            self.center = center
        else:
            self.center = None

        if radii is not None:
            self.radii = np.flip(np.sort(radii))
        else:
            self.radii = None

        if rotmat is not None:
            self.rotmat = rotmat
        else:
            self.rotmat = None

    # in place functions
    def H(self,coords):
        return coords[0] ** 2 + coords[1] ** 2 / (1 - self.e_e2) ** 2 + coords[2] / (1 - self.e_x2) ** 2  # eq. (43)

    def G(self, beta, lam):
        return np.sqrt(1 - self.e_x2 * np.cos(beta) ** 2 - self.e_e2 * np.sin(beta) ** 2 * np.sin(lam) ** 2)  # eq. (44)'

        # H = lambda coords: coords[0] ** 2 + coords[1] ** 2 / (1 - self.e_e2) ** 2 + coords[2] / (
        #             1 - self.e_x2) ** 2  # eq. (43)
        # G = lambda beta, lam: np.sqrt(
        #     1 - self.e_x2 * np.cos(beta) ** 2 - self.e_e2 * np.sin(beta) ** 2 * np.sin(lam) ** 2)  # eq. (44)'

    @property
    def radii(self):
        return self._radii

    @radii.setter
    def radii(self, radii):
        if radii is not None:
            self._radii = radii

            # squared semi axis
            self.radii2 = self._radii ** 2

            # linear eccentricities
            self.E_x = np.sqrt(self.radii2[0] - self.radii2[2])
            self.E_y = np.sqrt(self.radii2[1] - self.radii2[2])
            self.E_e = np.sqrt(self.radii2[0] - self.radii2[1])

            # linear eccentricities (squared)
            self.E_x2 = self.radii2[0] - self.radii2[2]
            self.E_y2 = self.radii2[1] - self.radii2[2]
            self.E_e2 = self.radii2[0] - self.radii2[1]

            # squared eccentricities
            self.e_x2 = (self.radii2[0] - self.radii2[2]) / self.radii2[0]
            self.e_e2 = (self.radii2[0] - self.radii2[1]) / self.radii2[0]

        else:
            # squared semi axis
            self.radii2 = None

            # linear eccentricities
            self.E_x = None
            self.E_y = None
            self.E_e = None

            # linear eccentricities (squared)
            self.E_x2 = None
            self.E_y2 = None
            self.E_e2 = None

            # squared eccentricities
            self.e_x2 = None
            self.e_e2 = None

    def get_normal(self, coords: np.ndarray, norm=False, return_norm=False):
        """
        Returns the normal vector in cartesian coordinates of a given point in cartesian coordinates

        Parameters
        ----------
        coords : np.ndarray of float [n_points x 3]
            Cartesian coordinates of given points
        norm : bool, optional, default: False
            Cartesian coordinates are given in normalized ellipsoidal space.
        return_norm : bool, optional, default: False
            Return normals in normalized space.

        Returns
        -------
        normal : np.ndarray of float [n_points x 3]
            Surface normal of given points (pointing outwards)
        """
        if coords.ndim == 1:
            coords = coords[np.newaxis, :]

        # rotate to normalized space if not already provided in normalized space
        if not norm:
            coords = self.rot_org2norm(coords=coords)

        # calculate normal in normalized space
        normal = 2 * np.hstack(((coords[:, 0] / self.radii2[0])[:, np.newaxis],
                                (coords[:, 1] / self.radii2[1])[:, np.newaxis],
                                (coords[:, 2] / self.radii2[2])[:, np.newaxis]))
        normal = normal / np.linalg.norm(normal, axis=1)[:, np.newaxis]

        # rotate normal to original space if desired
        if not return_norm:
            normal = (self.rotmat @ normal.T).T

        return normal

    def ellipsoid2cartesian(self, coords: np.ndarray, norm=False, return_normal=False) -> (np.ndarray,  np.ndarray):
        """
        Transforms spherical coordinates on ellipsoid to cartesian coordinates and optionally returns surface normal.

        Parameters
        ----------
        coords : np.ndarray of float [n_points x 2]
            Spherical coordinates [theta, phi] in radiant.
        norm : bool, optional, default: False
            Return coordinates in normalized space (ellipsoid centered at (0,0,0) and not rotated).
        return_normal : bool, optional, default: False
            Additionally return surface normal (pointing outwards)

        Returns
        -------
        coords_cart : np.ndarray of float [n_points x 3]
            Cartesian coordinates of given points
        normal : np.ndarray of float [n_points x 3]
            Surface normal of given points (pointing outwards)
        """

        if coords.ndim == 1:
            coords = coords[np.newaxis, :]

        phi = coords[:, 1]
        theta = coords[:, 0]

        # Cartesian coordinates that correspond to the spherical angles
        x = self.radii[0] * np.cos(phi) * np.sin(theta)
        y = self.radii[1] * np.sin(phi) * np.sin(theta)
        z = self.radii[2] * np.cos(theta)

        x_flatten = x.flatten()
        y_flatten = y.flatten()
        z_flatten = z.flatten()

        coords_cart = np.hstack((x_flatten[:, np.newaxis], y_flatten[:, np.newaxis], z_flatten[:, np.newaxis]))

        if not norm:
            coords_cart = self.rot_norm2org(coords=coords_cart)
            # coords_cart = ((self.rotmat @ coords_cart) + self.center[:, np.newaxis]).T

        if return_normal:
            normal = self.get_normal(coords=coords_cart, norm=norm, return_norm=norm)

            return coords_cart, normal
        else:
            return coords_cart

    def cartesian2ellipsoid(self, coords: np.ndarray, norm=False, return_norm=False, return_normal=False) -> (np.ndarray,  np.ndarray):
        """
        Transforms cartesian coordinates on ellipsoid to spherical coordinates and optionally returns surface normals.

        see eq. (3) - (7) in:
        Panou, G., & Korakitis, R. (2019). Geodesic equations and their numerical solution in Cartesian
        coordinates on a triaxial ellipsoid. Journal of Geodetic Science, 9(1), 1-12.

        Parameters
        ----------
        coords : np.ndarray of float [n_points x 3]
            Cartesian coordinates [x, y, z].
        norm : bool, optional, default: False
            Cartesian coordinates are given in normalized ellipsoidal space.
        return_norm : bool, optional, default: False
            Return normals in normalized space.
        return_normal : bool, optional, default: False
            Additionally return surface normal (pointing outwards)

        Returns
        -------
        coords_sphere : np.ndarray of float [n_points x 3]
            Spherical coordinates of given points [theta, phi], in radiant.
        normal : np.ndarray of float [n_points x 3]
            Surface normal of given points (pointing outwards)
        """

        if coords.ndim == 1:
            coords = coords[np.newaxis, :]

        if not norm:
            # rotate cartesian coordinates back to center and align with coordinate axes
            coords = self.rot_org2norm(coords=coords)
            # coords_rot = (coords - self.center) @ self.rotmat

        # quadrant (sign of x and y axes)
        mask_q1 = np.logical_and(coords[:, 0] >= 0, coords[:, 1] >= 0)
        mask_q2 = np.logical_and(coords[:, 0] < 0,  coords[:, 1] > 0)
        mask_q3 = np.logical_and(coords[:, 0] <= 0, coords[:, 1] <= 0)
        mask_q4 = np.logical_and(coords[:, 0] > 0,  coords[:, 1] < 0)

        # Cartesian coordinates that correspond to the spherical angles
        theta = np.arccos(coords[:, 2] / self.radii[2])
        phi = np.zeros(theta.shape)
        phi[mask_q1] = np.arctan(coords[mask_q1, 1] / coords[mask_q1, 0] * (self.radii[0] / self.radii[1]))
        phi[mask_q2] = np.arctan(coords[mask_q2, 1] / coords[mask_q2, 0] * (self.radii[0] / self.radii[1])) + np.pi
        phi[mask_q3] = np.arctan(coords[mask_q3, 1] / coords[mask_q3, 0] * (self.radii[0] / self.radii[1])) + np.pi
        phi[mask_q4] = np.arctan(coords[mask_q4, 1] / coords[mask_q4, 0] * (self.radii[0] / self.radii[1])) + 2 * np.pi

        coords_eli = np.hstack((theta[:, np.newaxis], phi[:, np.newaxis]))

        if return_normal:
            normal = self.get_normal(coords=coords, norm=True, return_norm=return_norm)

            return coords_eli, normal

        return coords_eli

    def cartesian2jacobi_system(self, x, coords):
        """

        :param x:
        :return:
        """
        beta = x[0]
        lam = x[1]

        # B = self.E_x2 * np.cos(beta)**2 + self.E_e2 * np.sin(beta)**2
        # L = self.E_x2 - self.E_e2 * np.cos(lam)**2
        #
        # f1 = self.radii[0] / self.E_x * np.sqrt(B) * np.cos(lam) - coords[0]
        # f2 = self.radii[1] * np.cos(beta) * np.sin(lam) - coords[1]
        # f3 = self.radii[2] / self.E_x * np.sin(beta) * np.sqrt(L) - coords[2]

        v = self.radii[0] / np.sqrt(
            1 - self.e_x2 * np.sin(beta) ** 2 - self.e_e2 * np.cos(beta) ** 2 * np.sin(lam) ** 2)

        # cartesian coordinates in normalized space
        f1 = v * np.cos(beta) * np.cos(lam) - coords[0]
        f2 = v * (1 - self.e_e2) * np.cos(beta) * np.sin(lam) - coords[1]
        f3 = v * (1 - self.e_x2) * np.sin(beta) - coords[2]

        f = np.sqrt(f1**2 + f2**2 + f3**2)

        return f

    def cartesian2jacobi(self, coords: np.ndarray, norm=False,  return_norm=False, return_normal=False) -> (np.ndarray, np.ndarray):
        """
        Transforms cartesian coordinates on ellipsoid to spherical coordinates in Jacobi form and optionally
        returns the surface normals.

        see eq. (8) - (16) and ()-() in:
        Panou, G., & Korakitis, R. (2019). Geodesic equations and their numerical solution in Cartesian
        coordinates on a triaxial ellipsoid. Journal of Geodetic Science, 9(1), 1-12.

        Parameters
        ----------
        coords : np.ndarray of float [n_points x 3]
            Cartesian coordinates [x, y, z]
        norm : bool, optional, default: False
            Coordinates are given in normalized space (ellipsoid centered at (0,0,0) and not rotated).
        return_norm : bool, optional, default: False
            Return normals in normalized space.
        return_normal : bool, optional, default: False
            Additionally return surface normal (pointing outwards)

        Returns
        -------
        coords_sphere : np.ndarray of float [n_points x 3]
            Spherical coordinates of given points [beta, lambda], in radiant.
            (beta: [-np.pi/2, +np.pi/2], lambda: [-np.pi, np.pi]
        normal : np.ndarray of float [n_points x 3]
            Surface normal of given points (pointing outwards)
        """
        if coords.ndim == 1:
            coords = coords[np.newaxis, :]

        # if coords are provided in original space, normalize them
        if not norm:
            coords = self.rot_org2norm(coords=coords)

        # calculate beta and lambda from normalized coordinates
        # beta = np.arctan((1-self.e_e2)/(1-self.e_x2) * coords[:, 2] / np.sqrt((1-self.e_e2)**2 * coords[:, 0]**2 + coords[:, 1]**2))
        # lam = np.arctan(1/(1-self.e_e2) * coords[:, 1] / coords[:, 0])

        # c0 = self.radii2[0] * self.radii2[1] + self.radii2[0]*self.radii2[2] + self.radii2[1]*self.radii2[2] \
        #      - (self.radii2[1] + self.radii2[2]) * coords[:, 0] \
        #      - (self.radii2[0] + self.radii2[2]) * coords[:, 1] \
        #      - (self.radii2[0] + self.radii2[1]) * coords[:, 2]
        # c1 = coords[:, 0]**2 + coords[:, 1]**2 + coords[:, 2]**2 - (self.radii2[0] + self.radii2[1] + self.radii2[2])
        # if np.isnan(np.sqrt(c1**2 - 4 * c0)):
        #     t2 = -c1 / 2.
        # else:
        #     t2 = (-c1 + np.sqrt(c1**2 - 4 * c0)) / 2.
        # t1 = c0 / t2
        # beta = np.arctan(np.sqrt((t1-self.radii2[2])/(self.radii2[1]-t1)))
        # lam = np.arctan(np.sqrt((t2-self.radii2[1])/(self.radii2[0]-t2)))

        beta = np.zeros(coords.shape[0])
        lam = np.zeros(coords.shape[0])

        for i, _coords in enumerate(coords):
            diff = np.inf
            while diff > 1:
                beta0 = np.random.rand(1)*2 - 1
                lambda0 = np.random.rand(1)*2 - 1
                res = minimize(self.cartesian2jacobi_system, (beta0[0], lambda0[0]),
                               method='Nelder-Mead',
                               bounds=((-np.pi/2, np.pi/2), (-np.pi, np.pi)),
                               args=_coords)
                diff = res.fun

            beta[i] = res.x[0]
            lam[i] = res.x[1]

        coords_jacobi = np.hstack((beta[:, np.newaxis], lam[:, np.newaxis]))

        if return_normal:
            normal = self.get_normal(coords=coords, norm=True, return_norm=return_norm)

            return coords_jacobi, normal

        return coords_jacobi

    def jacobi2cartesian(self, coords: np.ndarray, norm=False, return_norm=False, return_normal=False) -> (np.ndarray,  np.ndarray):
        """
        Transforms spherical coordinates on ellipsoid in Jacobi form to cartesian coordinates and optionally
        returns surface normals.

        Parameters
        ----------
        coords : np.ndarray of float [n_points x 2]
            Spherical coordinates [beta, lambda] in radiant.
        norm : bool, optional, default: False
            Return cartesian coordinates in normalized space (ellipsoid centered at (0,0,0) and not rotated)
        return_norm : bool, optional, default: False
            Return normals in normalized space.
        return_normal : bool, optional, default: False
            Additionally return surface normal (pointing outwards)

        Returns
        -------
        coords_cart : np.ndarray of float [n_points x 3]
            Cartesian coordinates of given points
        normal : np.ndarray of float [n_points x 3]
            Surface normal of given points (pointing outwards)
        """
        if coords.ndim == 1:
            coords = coords[np.newaxis, :]

        coords_cart = np.zeros((coords.shape[0], 3))

        v = lambda beta, lam: self.radii[0] / np.sqrt(
            1 - self.e_x2 * np.sin(beta) ** 2 - self.e_e2 * np.cos(beta) ** 2 * np.sin(lam) ** 2)

        v_ = v(coords[:, 0], coords[:, 1])

        # cartesian coordinates in normalized space
        coords_cart[:, 0] = v_ * np.cos(coords[:, 0]) * np.cos(coords[:, 1])
        coords_cart[:, 1] = v_ * (1 - self.e_e2) * np.cos(coords[:, 0]) * np.sin(coords[:, 1])
        coords_cart[:, 2] = v_ * (1 - self.e_x2) * np.sin(coords[:, 0])

        # transform to original space
        if not norm:
            coords_cart = self.rot_norm2org(coords=coords_cart)

        if return_normal:
            normal = self.get_normal(coords=coords_cart, norm=norm, return_norm=return_norm)
            return coords_cart, normal

        return coords_cart

    def rot_org2norm(self, coords):
        """
        Rotate cartesian coordinates from original to normalized space (ellipsoid centered and alligned with axes)

        Parameters
        ----------
        coords : np.ndarray of float [n_points, 3]
            Coordinates in original space

        Returns
        -------
        coords_norm : np.ndarray of float [n_points, 3]
            Coordinates in normalized space
        """
        return (coords - self.center) @ self.rotmat

    def rot_norm2org(self, coords):
        """
        Rotate cartesian coordinates from normalized (ellipsoid centered and alligned with axes) to original space.

        Parameters
        ----------
        coords : np.ndarray of float [n_points, 3]
            Coordinates in normalized space

        Returns
        -------
        coords_norm : np.ndarray of float [n_points, 3]
            Coordinates in original space
        """
        return (coords @ self.rotmat.T) + self.center[np.newaxis, :]

    def fit(self, points: np.ndarray):
        """
        Determines best fitting ellipsoid to a given point cloud.

        Parameters
        ----------
        points : np.ndarray of float [n_points x 3]
            Scattered points an ellipsoid is fitted to (cartesian)
        """
        # get convex hull
        hullV = ConvexHull(points)
        lH = len(hullV.vertices)
        hull = np.zeros((lH, 3))
        for i in range(len(hullV.vertices)):
            hull[i] = points[hullV.vertices[i]]
        hull = np.transpose(hull)

        # fit ellipsoid on convex hull
        eansa = ls_ellipsoid(hull[0], hull[1], hull[2])  # get ellipsoid polynomial coefficients
        center, radii, rotmat = polyToParams3D(eansa)    # get ellipsoid 3D parameters

        # sort axis such that radii[0] > radii[1] > radii[2]
        sort_idx = np.flip(np.argsort(radii))

        if (sort_idx == np.array([0, 1, 2])).all():
            sign_flip = np.array([1, 1, 1])
        elif (sort_idx == np.array([0, 2, 1])).all():
            sign_flip = np.array([1, 1, -1])
        elif (sort_idx == np.array([1, 0, 2])).all():
            sign_flip = np.array([1, 1, -1])
        elif (sort_idx == np.array([1, 2, 0])).all():
            sign_flip = np.array([1, 1, 1])
        elif (sort_idx == np.array([2, 0, 1])).all():
            sign_flip = np.array([1, 1, 1])
        elif (sort_idx == np.array([2, 1, 0])).all():
            sign_flip = np.array([1, 1, -1])
        else:
            sign_flip = np.array([1, 1, 1])

        self.center = center
        self.radii = radii[sort_idx]
        self.rotmat = rotmat[:, sort_idx] * np.tile(sign_flip, (3, 1))

    def get_geodesic_destination(self, start, alpha, distance, n_steps=400, method="Euler", n_cpu=None):
        """
        Calculates the destination point (in cartesian coordinates) on the triaxial ellipsoid using the geodesic.
        Originating from a start position (x, y, z)_start we walk into the direction alpha (angle with respect to
        constant lambda) for a given distance and end at a destination point (x, y, z)_end.

        Parameters
        ----------
        start : np.ndarray of float [n_start x 3]
            Starting points on ellipsoid in cartesian coordinates
        alpha : np.ndarray of float [n_start]
            Directions of travel (angle with respect to line of constant lambda)
        distance : np.ndarray of float [n_start]
            Traveling distances for each start point
        n_steps: int, optional, default: 1000
            Number of steps to solve differential equation
        method : str, optional, default: "euler"
            Integration method ("Euler", "RK45")
        n_cpu : int, optional, default: None
            Number of parallel threads (default: all available)

        Returns
        -------
        end : np.ndarray of float [n_start x 3]
            End points in cartesian coordinates
        """
        n_points = start.shape[0]

        # prepare input data
        input_data = [np.hstack(([start[i, :].flatten(), alpha[i], distance[i]])) for i in range(n_points)]

        if n_cpu is None:
            n_cpu = multiprocessing.cpu_count()
        n_cpu = min(n_cpu, len(input_data))

        workhorse = partial(self._workhorse_geodesic_destination_wrapper, n_steps=n_steps, method=method)

        pool = multiprocessing.Pool(n_cpu)
        coords_destination = np.vstack((pool.map(workhorse, input_data)))

        pool.close()
        pool.join()

        return coords_destination

    def _workhorse_geodesic_destination_wrapper(self, input_data, n_steps, method):
        return self._get_geodesic_destination(start=input_data[:3], alpha=input_data[3], distance=input_data[4],
                                              n_steps=n_steps, method=method)

    def _get_geodesic_destination(self, start, alpha, distance, n_steps=400, method="Euler", norm=False):
        """
        Calculates the destination point (in cartesian coordinates) on the triaxial ellipsoid using the geodesic.
        Originating from a start position (x, y, z)_start we walk into the direction alpha (angle with respect to
        constant lambda) for a given distance and end at a destination point (x, y, z)_end.

        Parameters
        ----------
        start : np.ndarray of float [3]
            Starting point on ellipsoid in cartesian coordinates
        alpha : float
            Direction of travel (angle with respect to line of constant lambda)
        distance : float
            Traveling distance
        n_steps: int, optional, default: 1000
            Number of steps to solve differential equation
        method : str, optional, default: "euler"
            Integration method ("Euler", "DOP853", "RK45", ...)
        norm : bool, optional, default: False
            Coordinates are provided in normalized space (ellipsoid centered and rotated)

        Returns
        -------
        end : np.ndarray of float [3]
            End point in cartesian coordinates
        """
        # transform coordinates to normalized space (ellipsoid centered and axes oriented)
        start = self.rot_org2norm(coords=start)
        # start = (self.rotmat.T @ (start - self.center).T).T

        initial_conditions = self.get_initial_conditions(start=start, alpha=alpha, norm=True)

        # Solve
        ################################################################################################################
        s = np.linspace(0, 1.1*distance, n_steps)
        ds = 1.1*distance / n_steps
        distance_current = 0

        if method == "Euler":
            x1 = np.zeros(n_steps)
            x2 = np.zeros(n_steps)
            y1 = np.zeros(n_steps)
            y2 = np.zeros(n_steps)
            z1 = np.zeros(n_steps)
            z2 = np.zeros(n_steps)
            x1[0] = initial_conditions[0]
            x2[0] = initial_conditions[1]
            y1[0] = initial_conditions[2]
            y2[0] = initial_conditions[3]
            z1[0] = initial_conditions[4]
            z2[0] = initial_conditions[5]
            Cx = lambda t, x: self.radii2[0]*x/(t+self.radii2[0])
            Cy = lambda t, y: self.radii2[1]*y/(t+self.radii2[1])
            Cz = lambda t, z: self.radii2[2]*z/(t+self.radii2[2])

            t = -50
            # walk along path in steps of ds and update coordinates
            for i, s_ in enumerate(s[:-1]):
                # apply correction step after every 20 steps (project point to ellipsoid surface)
                # https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
                if not (i % 5):
                    t = fsolve(self.F_t, t,
                               args=np.array([x1[i], y1[i], z1[i]]).astype(np.float64),
                               xtol=1e-9,
                               maxfev=200)
                    x1[i] = Cx(t, x1[i])
                    y1[i] = Cy(t, y1[i])
                    z1[i] = Cz(t, z1[i])

                # constants
                h = x2[i] ** 2 + y2[i] ** 2 / (1 - self.e_e2) + z2[i] ** 2 / (1 - self.e_x2)  # eq. (44)
                H = x1[i] ** 2 + y1[i] ** 2 / (1 - self.e_e2) + z1[i] ** 2 / (1 - self.e_x2)  # eq. (43)
                h_H = h / H

                # update equations of system of 1st order ODE
                x1[i + 1] = x2[i] * ds + x1[i]
                x2[i + 1] = -h_H * x1[i] * ds + x2[i]
                y1[i + 1] = y2[i] * ds + y1[i]
                y2[i + 1] = -h_H * y1[i] / (1. - self.e_e2) * ds + y2[i]
                z1[i + 1] = z2[i] * ds + z1[i]
                z2[i + 1] = -h_H * z1[i] / (1. - self.e_x2) * ds + z2[i]

                distance_current += np.linalg.norm(np.array([x1[i+1]-x1[i], y1[i+1]-y1[i], z1[i+1]-z1[i]]))

                if distance_current > distance:
                    break

            coords_destination = np.array([x1[i+1], y1[i+1], z1[i+1]])
            #coords_path = np.hstack((x1[:, np.newaxis], y1[:, np.newaxis], z1[:, np.newaxis]))

        else:
            x_rk = solve_ivp(self.direct_geodesic_cartesian, [0, distance], initial_conditions,
                             method=method,
                             t_eval=np.linspace(0, distance, 1000),
                             atol=1e-12,
                             rtol=1e-9)
            x_rk = x_rk.y
            coords_destination = np.array([x_rk[0, -1], x_rk[2, -1], x_rk[4, -1]])
            #coords_path = np.hstack((x_rk[0, :][:, np.newaxis], x_rk[2, :][:, np.newaxis], x_rk[4, :][:, np.newaxis]))

        #distance_check = np.sum(np.linalg.norm(coords_path[1:(i+2), :] - coords_path[:(i+1), :], axis=1))
        #print(f"target distance: {distance}, actual distance: {distance_check}, distance_current: {distance_current}")

        # rotate coordinates back to original space
        coords_destination = self.rot_norm2org(coords=coords_destination)

        return coords_destination

        # import matplotlib
        # matplotlib.use('Qt5Agg')
        # import matplotlib.pyplot as plt
        # fig = plt.figure()
        # ax = fig.add_subplot(projection='3d')
        #
        # theta = np.linspace(0, np.pi, 30)
        # phi = np.linspace(0, 2 * np.pi, 50)
        # coords_sphere = np.array(np.meshgrid(theta, phi)).T.reshape(-1, 2)
        # # coords_sphere = np.hstack([theta[:, np.newaxis], phi[:, np.newaxis]])
        # eli_coords = self.ellipsoid2cartesian(coords=coords_sphere, return_normal=False).astype(np.float64)
        # eli_coords_rot = (self.rotmat.T @ (eli_coords - self.center).T).T.astype(np.float64)
        #
        # ax.scatter(x1.astype(np.float64), y1.astype(np.float64), z1.astype(np.float64))
        # ax.scatter(eli_coords_rot[:, 0], eli_coords_rot[:, 1], eli_coords_rot[:, 2])
        # # ax.scatter(x_rk[0, :], x_rk[2, :], x_rk[4, :])
        # ax.set_xlabel("x")
        # ax.set_ylabel("y")
        # ax.set_zlabel("z")
        # plt.savefig("/home/kporzig/tmp/plots/ellipse_geodesic_test.png", dpi=300)

        # x = odeint(self.direct_geodesic_cartesian, initial_conditions, s)
        # start_time = time.time()
        # x_rk = scipy.integrate.solve_ivp(self.direct_geodesic_cartesian, [0, distance], initial_conditions,
        #                                  method='RK45',
        #                                  t_eval=np.linspace(0,distance, 1000),
        #                                  atol=1e-12,
        #                                  rtol=1e-9)
        # x_rk = x_rk.y
        # end_time = time.time()
        # print(end_time - start_time)

    def get_initial_conditions(self, start, alpha, norm=False):
        """
        Determine initial conditions of geodesic walk algorithm

        Parameters
        ----------
        start : np.ndarray of float [3]
            Starting point on ellipsoid in cartesian coordinates
        alpha : float
            Direction of travel (angle with respect to line of constant lambda)
        norm : bool, optional, default: False
            Coordinates are given in normalized space (ellipsoid centered and axes oriented)

        Returns
        -------
        initial_conditions : np.ndarray of float [6]
            Initial conditions [x0, dxds_0, y0, dyds_0, z0, dzds_0]
        """
        if start.ndim == 2:
            start = start.flatten()

        # Parameter definition and initialization
        ################################################################################################################
        if not norm:
            start = self.rot_org2norm(coords=start)
            # start = (self.rotmat.T @ (start - self.center).T).T

        start_jacobi = self.cartesian2jacobi(coords=start, norm=norm).flatten()

        t1 = lambda beta: self.radii2[1]*np.sin(beta)**2 + self.radii2[2]*np.cos(beta)**2                   # eq. (13)
        t2 = lambda lam: self.radii2[0]*np.sin(lam)**2 + self.radii2[1]*np.cos(lam)**2                      # eq. (14)
        L = lambda lam: t2(lam) - self.radii2[2]                                                            # eq. (67)
        F = lambda lam, beta: t2(lam) - t1(beta)                                                            # eq. (68)
        B = lambda beta: self.E_x2/self.E_y2 * (self.radii2[1] - t1(beta)) + \
                         self.E_e2/self.E_y2 * (t1(beta) - self.radii2[2])                                  # eq. (66)

        # Initial conditions
        ################################################################################################################
        # Dirichlet boundary conditions
        x0 = start[0]
        y0 = start[1]
        z0 = start[2]
        beta0 = start_jacobi[0]
        lambda0 = start_jacobi[1]

        # Calculate Neumann boundary conditions
        H0 = self.H(start)
        H0_sqrt = np.sqrt(H0)
        G0 = self.G(beta0, lambda0)

        # normal vector to start point on ellipsoid, eq. (58)
        _, n = self.cartesian2ellipsoid(coords=start, norm=norm, return_norm=True, return_normal=True)
        n = n.flatten()
        # n = np.zeros(3)
        # n[0] = x0 / H0_sqrt
        # n[1] = y0 / ((1 - self.e_e2) * H0_sqrt)
        # n[2] = z0 / ((1 - self.e_x2) * H0_sqrt)

        # unit vector tangent to line of constant beta, eq. (63)-(65)
        # L0 = L(lambda0)
        # F0 = F(lambda0, beta0)
        # B0 = B(beta0)
        # t10 = t1(beta0)
        # t20 = t2(lambda0)
        # p = np.zeros(3)
        # p[0] = -np.sign(y0) * np.sqrt(L0/(F0*t20)) * self.radii[0]/(self.E_x*self.E_e) * \
        #        np.sqrt(B0) * np.sqrt(t20-self.radii2[1])
        # p[1] = np.sign(x0) * np.sqrt(L0/(F0*t20)) * self.radii[1]/(self.E_y*self.E_e) * \
        #        np.sqrt((self.radii2[1]-t10) * (self.radii2[0]-t20))
        # p[2] = np.sign(x0)*np.sign(y0)*np.sign(z0) * (1/np.sqrt(F0*t20)) * self.radii[2]/(self.E_x*self.E_y) * \
        #        np.sqrt((t10-self.radii2[2])*(t20-self.radii2[1])*(self.radii2[0]-t20))

        # unit vector tangent to line of constant lambda, eq. (71)-(73)
        # q = np.cross(n, p)
        q = np.zeros(3)
        q[0] = -1 / G0 * np.sin(beta0) * np.cos(lambda0)
        q[1] = -1 / G0 * np.sqrt(1 - self.e_e2) * np.sin(beta0) * np.sin(lambda0)
        q[2] = 1 / G0 * np.sqrt(1 - self.e_x2) * np.cos(beta0)

        # unit vector tangent to line of constant beta, eq. (63)-(65)
        p = np.cross(q, n)

        # unit vector tangent to geodesic, eq. (57)
        sigma = p * np.sin(alpha) + q * np.cos(alpha)

        # Neumann
        dxds_0 = sigma[0]
        dyds_0 = sigma[1]
        dzds_0 = sigma[2]

        # Initial conditions
        initial_conditions = np.array([x0, dxds_0, y0, dyds_0, z0, dzds_0])

        return initial_conditions

    def F_t(self, t, coords):
        """
        Function to find root for t to find closest point on ellipsoid

        Parameters
        ----------
        t : float
            t-value
        coords : np.ndarray of float [3]
            Coordinates of query point
        Returns
        -------
        y : float
            Funtion value of F(t)
        """
        y = (self.radii[0] * coords[0] / (t + self.radii2[0])) ** 2 + \
            (self.radii[1] * coords[1] / (t + self.radii2[1])) ** 2 + \
            (self.radii[2] * coords[2] / (t + self.radii2[2])) ** 2 - 1
        return y

    def coords_correction(self, coords, norm=False, return_norm=False):
        """
        Apply correction of coordinates on ellipsoid by finding closest point on ellipsoid.
        https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf

        Parameters
        ----------
        coords : np.ndarray of float [3]
            Coordinates of query point
        norm : bool, optional, default: False
            Coordinates are given in normalized space
        return_norm : bool, optional, default: False
            Return corrected coordinates in normalized space

        Returns
        -------
        coords : np.ndarray of float [3]
            Coordinates of query point on ellipsoid (closest point to ellipsoid surface)
        """
        # normalize coordinates if they are given in original space
        if not norm:
            coords = self.rot_org2norm(coords)

        t = fsolve(self.F_t, -50.,
                   args=np.array([coords[0], coords[1], coords[2]]).astype(np.float64),
                   xtol=1e-6,
                   maxfev=1000)

        coords[0] = self.radii2[0]*coords[0]/(t+self.radii2[0])
        coords[1] = self.radii2[1]*coords[1]/(t+self.radii2[1])
        coords[2] = self.radii2[2]*coords[2]/(t+self.radii2[2])

        if not return_norm:
            coords = self.rot_norm2org(coords=coords)

        return coords

    def direct_geodesic_cartesian(self, s, x):
        """
        System of differential equations to determine the (direct) geodesic in cartesian coordinates.
        It determines the end point given a starting point, a direction and a distance.

        see eq. (51)-(56) in:
        Panou, G., & Korakitis, R. (2019). Geodesic equations and their numerical solution in Cartesian
        coordinates on a triaxial ellipsoid. Journal of Geodetic Science, 9(1), 1-12.

        Parameters
        ----------
        x : np.ndarray of float [6]
            Coordinates and derivatives
        s : np.ndarray of float [n_steps]
            Distance array with discrete steps of length ds=distance/n_steps between [0, distance]

        Returns
        -------
        x : np.ndarray of float [6]
            Coordinates and derivatives
        """
        # apply coordinate correction
        x_corr = self.coords_correction(coords=np.array([x[0], x[2], x[4]]), norm=True, return_norm=True)
        x[0] = x_corr[0]
        x[2] = x_corr[1]
        x[4] = x_corr[2]

        # determine constants
        h = x[1]**2 + x[3]**2/(1-self.e_e2) + x[5]**2/(1-self.e_x2)     # eq. (44)
        H = x[0]**2 + x[2]**2/(1-self.e_e2) + x[4]**2/(1-self.e_x2)     # eq. (43)
        h_H = h / H

        # system of 1st order ODE
        dx1dt = x[1]
        dx2dt = -h_H * x[0]
        dy1dt = x[3]
        dy2dt = -h_H * x[2] / (1-self.e_e2)
        dz1dt = x[5]
        dz2dt = -h_H * x[4] / (1-self.e_x2)

        return [dx1dt, dx2dt, dy1dt, dy2dt, dz1dt, dz2dt]


def ls_ellipsoid(xx: np.ndarray, yy: np.ndarray, zz: np.ndarray) -> np.ndarray:
    """
    Finds best fit ellipsoid. (http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html)
    least squares fit to a 3D-ellipsoid:

    Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz  = 1

    Note that sometimes it is expressed as a solution to
    Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz  = 1
    where the last six terms have a factor of 2 in them
    This is in anticipation of forming a matrix with the polynomial coefficients.
    Those terms with factors of 2 are all off diagonal elements.  These contribute
    two terms when multiplied out (symmetric) so would need to be divided by two
    """
    x = xx[:, np.newaxis]
    y = yy[:, np.newaxis]
    z = zz[:, np.newaxis]

    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz = 1
    J = np.hstack((x * x, y * y, z * z, x * y, x * z, y * z, x, y, z))
    K = np.ones_like(x)
    JT = J.transpose()
    JTJ = np.dot(JT, J)
    InvJTJ = np.linalg.inv(JTJ)
    ABC = np.dot(InvJTJ, np.dot(JT, K))

    # Rearrange, move the 1 to the other side
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz - 1 = 0
    #    or
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz + J = 0
    #  where J = -1
    eansa = np.append(ABC, -1)

    return (eansa)


def polyToParams3D(vec: np.ndarray) -> (np.ndarray,  np.ndarray,  np.ndarray):
    """
    Gets 3D parameters of an ellipsoid. (http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html)
    Convert the polynomial form of a 3D-ellipsoid to parameters center, axes, and transformation matrix.
    Algebraic form: X.T * Amat * X --> polynomial form

    Parameters
    ----------
    vec : np.ndarray of float
        Vector whose elements are the polynomial coefficients A...J

    Returns
    -------
    center : np.ndarray of float [3]
        Center of ellipsoid
    radii : np.ndarray of float [3]
        Axes length
    rotmat : np.ndarray of float [3 x 3]
        Rotation matrix of ellipsoid (direction of principal axes)
    """

    # polynomial
    Amat = np.array(
        [
            [vec[0], vec[3] / 2.0, vec[4] / 2.0, vec[6] / 2.0],
            [vec[3] / 2.0, vec[1], vec[5] / 2.0, vec[7] / 2.0],
            [vec[4] / 2.0, vec[5] / 2.0, vec[2], vec[8] / 2.0],
            [vec[6] / 2.0, vec[7] / 2.0, vec[8] / 2.0, vec[9]]
        ])

    # Algebraic form of polynomial
    # See B.Bartoni, Preprint SMU-HEP-10-14 Multi-dimensional Ellipsoidal Fitting
    # eq. (20) for the following method for finding the center
    # https://www.physics.smu.edu/~scalise/SMUpreprints/SMU-HEP-10-14.pdf
    A3 = Amat[0:3, 0:3]
    A3inv = np.linalg.inv(A3)
    ofs = vec[6:9] / 2.0

    # center of ellipsoid
    center = -np.dot(A3inv, ofs)

    # Center the ellipsoid at the origin
    Tofs = np.eye(4)
    Tofs[3, 0:3] = center
    R = np.dot(Tofs, np.dot(Amat, Tofs.T))

    # Algebraic form translated to center
    R3 = R[0:3, 0:3]
    s1 = -R[3, 3]
    R3S = R3 / s1
    (el, ec) = np.linalg.eig(R3S)
    recip = 1.0 / np.abs(el)

    # radii
    radii = np.sqrt(recip)

    # rotation matrix
    rotmat = np.linalg.inv(ec)

    return (center, radii, rotmat)


def check_ellipsoid_fit(center, radii, rotmat, points):
    """
    Check solution of ellipsoidal fit
    Convert to unit sphere centered at origin
     1) Subtract off center
     2) Rotate points so bulges are aligned with axes (no xy,xz,yz terms)
     3) Scale the points by the inverse of the axes gains
     4) Back rotate
    Rotations and gains are collected into single transformation matrix M

    Parameters
    ----------
    center : np.ndarray of float [3]
        Center of ellipsoid
    radii : np.ndarray of float [3]
        Axes length
    rotmat : np.ndarray of float [3 x 3]
        Rotation matrix of ellipsoid (direction of principal axes)
    points : np.ndarray of float [n_points x 3]
        XYZ coordinates of points to fit.
    Returns
    -------
    average_radius : float
        Average radius (truth is 1.0)
    std_radius : float
        Standard deviation of radius
    """
    # subtract the offset so ellipsoid is centered at origin
    xc = points[:, 0] - center[0]
    yc = points[:, 1] - center[1]
    zc = points[:, 2] - center[2]

    # create transformation matrix
    L = np.diag([1 / radii[0], 1 / radii[1], 1 / radii[2]])
    M = np.dot(rotmat.T, np.dot(L, rotmat))

    # apply the transformation matrix
    [xm, ym, zm] = np.dot(M, [xc, yc, zc])
    # Calculate distance from origin for each point (ideal = 1.0)
    rm = np.sqrt(xm * xm + ym * ym + zm * zm)

    average_radius = np.mean(rm)
    std_radius = np.std(rm)

    return average_radius, std_radius


def subject2ellipsoid(coords: np.ndarray, normals: np.ndarray, ellipsoid: Ellipsoid) -> np.ndarray:
    """
    Transform coordinates from subject to ellipsoid space
    Line intersection of skin surface point along normal and ellipsoid.
    https://math.stackexchange.com/questions/3309397/line-ellipsoid-intersection

    Parameters
    ----------
    coords : np.ndarray of float [n_points x 3]
        Cartesian coordinates in subject space (x, y, z)
    normals : np.ndarray of float [n_points x 3]
        Surface normals of points
    ellipsoid : Ellipsoid class instance
        Ellipsoid the data

    Returns
    -------
    coords : np.ndarray of float [n_coords x 2]
        Spherical coordinates [theta, phi] in radiant.
    """

    if coords.ndim == 1:
        coords = coords[np.newaxis, :]

    if normals.ndim == 1:
        normals = normals[np.newaxis, :]

    R = ellipsoid.rotmat
    A = np.diag(np.array(1/(ellipsoid.radii**2)))
    v = coords - ellipsoid.center
    B = R @ A @ R.transpose()

    p = np.zeros(3)
    p[0] = normals @ B @ normals.transpose()
    p[1] = 2 * v @ B @ normals.transpose()
    p[2] = v @ B @ v.transpose() - 1

    lam = np.roots(p)
    lam = np.min(np.abs(lam))

    eli_intersect_cart = coords + lam * normals

    eli_intersect_sphere = ellipsoid.cartesian2ellipsoid(coords=eli_intersect_cart, return_normal=False)

    return eli_intersect_sphere


def ellipsoid2subject(coords: np.ndarray, ellipsoid: Ellipsoid, surface: Surface) -> (int, np.ndarray):
    """
    Transform coordinates from ellipsoid to given surface by projecting the points given in "coords" (on ellipsoid)
    in normal direction to the given surface. Returns the index and the coordinates on the surface of the triangle hit.
    (uses raytracer of cgal)

    Parameters
    ----------
    coords : np.ndarray of float [n_points x 2]
        Spherical coordinates [theta, phi] in radiant on the ellipsoid surface.
    ellipsoid : Ellipsoid class instance
        Ellipsoid
    surface : Surface class instance
        Surface of subject the coordinates are projected on

    Returns
    -------
    indices : int
        Triangle indices on surface mesh
    coords : np.ndarray of float [n_coords x 3]
        Intersection points (x, y, z) on surface.
    """
    if coords.ndim == 1:
        coords = coords[np.newaxis, :]

    # get cartesian coordinates and normal
    x0, n = ellipsoid.ellipsoid2cartesian(coords=coords, return_normal=True)

    # define near and far point
    near = (x0 - 100*n)
    far = (x0 + 100*n)

    indices, points = cgal.segment_triangle_intersection(
        surface.nodes,
        surface.tr_nodes,
        [near[i, :] for i in range(coords.shape[0])],
        [far[i, :] for i in range(coords.shape[0])]  # near[i, :], far[i, :]
    )

    return (indices, points)



# coords_jacobi = np.zeros((coords.shape[0], 2))
        # eps0 = 1e-15
        # eps = np.inf
        # beta = 1.
        # lam = 1.
        # J = np.zeros((3, 2))
        # L = lambda lam: self.E_x2 - self.E_e2 * np.cos(lam)**2
        # B = lambda beta: self.E_x2*np.cos(beta)**2 + self.E_e2*np.sin(beta)**2
        #
        # for i_point, point in enumerate(coords):
        #     while eps > eps0:
        #         try:
        #             B_sqrt = np.sqrt(B(beta))
        #             L_sqrt = np.sqrt(L(lam))
        #         except RuntimeWarning:
        #             beta = np.random.rand(1)
        #             lam = np.random.rand(1)
        #
        #         J[0, 0] = -self.radii[0]*self.E_y2/(2*self.E_x)*np.sin(2*beta)/B_sqrt*np.cos(lam)   # dx/dbeta
        #         J[1, 0] = -self.radii[1]*np.sin(beta)*np.sin(lam)                                   # dy/dbeta
        #         J[2, 0] = self.radii[2]/self.E_x*np.cos(beta)*L_sqrt                                # dz/dbeta
        #         J[0, 1] = -self.radii[0]/self.E_x*B_sqrt*np.sin(lam)                                # dx/dlam
        #         J[1, 1] = self.radii[1]*np.cos(beta)*np.cos(lam)                                    # dy/dlam
        #         J[2, 1] = self.radii[2]*self.E_e2/(2*self.E_x)*np.sin(beta)*np.sin(2*lam)/L_sqrt    # dz/dlam
        #
        #         N = J.T @ J
        #         N_inv = np.linalg.inv(N)
        #         # N_inv = 1/(N[0, 0]*N[1, 1]-N[0, 1]*N[1, 0]) * np.array([[N[1, 1], -N[0, 1]], [-N[1, 0],  N[0, 0]]])
        #         coords_current = self.jacobi2cartesian(coords=np.array([[beta, lam]]), norm=True)
        #         dI = point - coords_current
        #
        #         dbeta_dlam = N_inv @ J.T @ dI.T  # np.linalg.inv(N)
        #
        #         beta += dbeta_dlam[0]
        #         lam += dbeta_dlam[1]
        #
        #         eps = np.linalg.norm(dbeta_dlam)
        #
        #         print(coords, coords_current, eps) #coords_current, coords,
        #
        #     coords_jacobi[i_point, :] = np.hstack((beta, lam))
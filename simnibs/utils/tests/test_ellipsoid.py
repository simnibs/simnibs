import os

import numpy as np
import pytest

# from .. import SIMNIBSDIR
# from ...mesh_tools import mesh_io
# from .. import transformations

from simnibs.utils.ellipsoid import Ellipsoid

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt


class TestEllipsoid:

    def test_fit(self):
        """ Tests ellipsoid fit """
        # parameters of test ellipsoid
        center = np.array([10, 20, 30])
        radii = np.array([10, 20, 30])
        rotmat = np.array([[1, 0, 0],
                           [0, np.cos(30. / 180. * np.pi), np.sin(30. / 180. * np.pi)],
                           [0, -np.sin(30. / 180. * np.pi), np.cos(30. / 180. * np.pi)]])
        ellipsoid = Ellipsoid(center=center, radii=radii, rotmat=rotmat)

        # generate some samples on test ellipsoid
        theta = np.linspace(0, np.pi, 30)
        phi = np.linspace(0, 2 * np.pi, 20)
        coords_sphere = np.array(np.meshgrid(theta, phi)).T.reshape(-1, 2)
        eli_coords = ellipsoid.ellipsoid2cartesian(coords=coords_sphere, return_normal=False)

        # fit ellipsoid to test data
        ellipsoid_fit = Ellipsoid()
        ellipsoid_fit.fit(points=eli_coords)

        # compare ellipsoid parameters
        np.testing.assert_allclose(ellipsoid.center, ellipsoid_fit.center,
                                    err_msg="Ellipsoid fit failed: center does not match!")
        np.testing.assert_allclose(ellipsoid.radii, ellipsoid_fit.radii,
                                   err_msg="Ellipsoid fit failed: radii do not match!")
        np.testing.assert_allclose(ellipsoid.rotmat, ellipsoid_fit.rotmat,
                                   err_msg="Ellipsoid fit failed: rotmat does not match!")

    def test_transformations(self):
        """ Tests coordinates transformations (cartesian, spherical, jacobi) """
        ellipsoid = Ellipsoid(center=np.array([10, 20, 30]),
                              radii=np.array([10, 20, 30]),
                              rotmat=np.array([[1, 0, 0],
                                               [0,  np.cos(30./180.*np.pi), np.sin(30./180.*np.pi)],
                                               [0, -np.sin(30./180.*np.pi), np.cos(30./180.*np.pi)]]))

        # base coordinates [theta, phi], theta: 0:pi, phi: 0:2*pi
        coords_eli = np.array([[  np.pi/4,   2*np.pi/4],
                               [3*np.pi/4, 3*2*np.pi/4],
                               [  np.pi/4,   2*np.pi/4],
                               [3*np.pi/4, 3*2*np.pi/4]])

        for coords_test in coords_eli:
            # transform back and forth (ellipsoid <-> cartesian)
            coords_eli2cart = ellipsoid.ellipsoid2cartesian(coords=coords_test)
            coords_cart2eli = ellipsoid.cartesian2ellipsoid(coords=coords_eli2cart)

            assert np.isclose(coords_test, coords_cart2eli).all(), \
                "Cartesian coordinates of ellipsoidal transformation differ compared to reference"

            # transform back and forth (jacobian <-> cartesian)
            coords_cart2jac = ellipsoid.cartesian2jacobi(coords=coords_eli2cart)  # coords_eli2cart
            coords_jac2cart = ellipsoid.jacobi2cartesian(coords=coords_cart2jac)

            assert np.isclose(coords_eli2cart, coords_jac2cart, rtol=1e-4).all(), \
                "Cartesian coordinates of jacobian transformation differ compared to reference"

    def test_geodesic_distance(self):
        """ Tests geodesic distance calculation """
        # test ellipsoid
        center = np.array([2.00564077, -16.61394422, 13.9095578])
        radii = np.array([103.45868648, 91.52166666, 73.89797068])
        rotmat = np.array([[-0.09566444, -0.04310434, 0.99447993],
                           [-0.89824862, 0.43425991, -0.06758504],
                           [0.42894956, 0.89975571, 0.08026164]])
        ellipsoid = Ellipsoid(center=center, radii=radii, rotmat=rotmat)

        # start and end points in jacobi coordinates [beta, lambda], beta: -pi/2:pi/2, lambda: -pi:pi
        start_jacobi = np.array([[np.pi/4, np.pi/2],
                                 [-np.pi/4, np.pi/2],
                                 [-np.pi/4, -np.pi/2],
                                 [np.pi/4, -np.pi/2]])

        for _start_jacobi in start_jacobi:
            # directions
            alpha = np.linspace(-np.pi, np.pi, 25)

            # traveling distance
            distance = 20. * np.ones(len(alpha))

            # start point in cartesian coordinates
            start_cart = np.tile(ellipsoid.jacobi2cartesian(coords=_start_jacobi), (len(alpha), 1))

            # determine geodesic distances
            end = ellipsoid.get_geodesic_destination(start=start_cart,
                                                     alpha=alpha,
                                                     distance=distance,
                                                     n_steps=400,
                                                     method="Euler")

            # difference between the endpoints
            diff_between = np.linalg.norm(end[:-1]-end[1:], axis=1)
            diff_between_rstd = np.std(diff_between)/np.mean(diff_between)

            # difference between end and start point
            diff_maxrel = np.max(np.abs(np.linalg.norm(start_cart - end, axis=1) - distance))/distance[0]

            assert diff_maxrel < 0.02, \
                f"Wrong geodesic distance travelled (diff_maxrel = {diff_maxrel*100}% > 2%)"
            assert diff_between_rstd < 0.01, \
                f"End points not equally spaced (diff_between_rstd = {diff_between_rstd * 100}% > 1%)"




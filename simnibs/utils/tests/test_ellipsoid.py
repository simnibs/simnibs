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
    def test_init(self):
        pass

    def test_geodesic_distance(self):
        pass

    def test_fit(self):
        pass

    def test_transformations(self):
        ellipsoid = Ellipsoid(center=np.array([10, 20, 30]),
                              radii=np.array([10, 20, 30]),
                              rotmat=np.array([[1, 0, 0],
                                               [0,  np.cos(30./180.*np.pi), np.sin(30./180.*np.pi)],
                                               [0, -np.sin(30./180.*np.pi), np.cos(30./180.*np.pi)]]))

        # base coordinates
        coords_eli = np.array([0.2, 0.5])  # [theta, phi]

        # transform back and forth (ellipsoid <-> cartesian)
        coords_eli2cart = ellipsoid.ellipsoid2cartesian(coords=coords_eli)
        coords_cart2eli = ellipsoid.cartesian2ellipsoid(coords=coords_eli2cart)

        assert np.isclose(coords_eli, coords_cart2eli).all(), "Cartesian coordinates of ellipsoidal transformation " \
                                                              "differ compared to reference"

        # transform back and forth (jacobian <-> cartesian)
        coords_cart2jac = ellipsoid.cartesian2jacobi(coords=coords_eli2cart)  # coords_eli2cart
        coords_jac2cart = ellipsoid.jacobi2cartesian(coords=coords_cart2jac)

        assert np.isclose(coords_eli2cart, coords_jac2cart).all(), "Cartesian coordinates of jacobian transformation " \
                                                                   "differ compared to reference"

        # import matplotlib
        # matplotlib.use('Qt5Agg')
        # import matplotlib.pyplot as plt
        #
        # beta = np.linspace(-np.pi/2, np.pi/2, 10)
        # lam = np.linspace(0, 2 * np.pi, 30)
        #
        # # theta = np.linspace(0, np.pi, 180)
        # # phi = np.linspace(0, 2 * np.pi, 360)
        # theta = np.linspace(0, np.pi, 10)
        # phi = np.linspace(0, 2 * np.pi, 20)
        # coords_sphere = np.array(np.meshgrid(theta, phi)).T.reshape(-1, 2)
        # coords_sphere_jac = np.array(np.meshgrid(beta, lam)).T.reshape(-1, 2)
        # # coords_sphere = np.hstack([theta[:, np.newaxis], phi[:, np.newaxis]])
        # eli_coords = ellipsoid.ellipsoid2cartesian(coords=coords_sphere, return_normal=False)
        # eli_coords_jac = ellipsoid.jacobi2cartesian(coords=coords_sphere_jac, return_normal=False)
        # # eli_coords_rot = (ellipsoid.rotmat.T @ (eli_coords - ellipsoid.center).T).T
        # eli_coords_rot = (eli_coords - ellipsoid.center) @ ellipsoid.rotmat
        # eli_coords_backrot = (eli_coords_rot @ ellipsoid.rotmat.T) + ellipsoid.center[np.newaxis, :]
        # # coords_sphere_test = self.ellipsoid.cartesian2ellipsoid(coords=eli_coords)
        # # np.isclose(coords_sphere, coords_sphere_test)
        # fig = plt.figure()
        # ax = fig.add_subplot(projection='3d')
        # ax.scatter(eli_coords[:, 0], eli_coords[:, 1], eli_coords[:, 2])
        # ax.scatter(eli_coords_jac[:, 0], eli_coords_jac[:, 1], eli_coords_jac[:, 2])
        # ax.scatter(eli_coords_rot[:, 0], eli_coords_rot[:, 1], eli_coords_rot[:, 2])
        # ax.scatter(eli_coords_backrot[:, 0], eli_coords_backrot[:, 1], eli_coords_backrot[:, 2])

    def test_cartesian2ellipsoid(self):
        pass

    def test_ellipsoid2cartesian(self):
        pass

    def test_cartesian2jacobi(self):
        pass

    def test_jacobi2cartesian(self):
        pass



# # start = np.array([-56.8752571, -12.5540671, 73.08751634])
# # start = np.array([56.28987403, -14.76619973,  74.3554083])
# # start = np.array([-30.48893514, -32.2734708 , 100.19494273])
# start_jacobi = np.array([0.1, 0.5])
# alpha = np.pi/4 # 0. #np.pi/2*0.8
# distance = 100
#
# center = np.array([  1.9327103 , -16.91330219,  14.18419223])
# # radii = np.array([102.96377824,  91.28007274,  73.89982405])
# radii = np.array([102.96377824,  100.,  98])
# rotmat = np.array([[-0.09554968, -0.04621891,  0.99435108],
#                    [-0.90557844,  0.41876494, -0.06755448],
#                    [ 0.41327708,  0.90691771,  0.0818677 ]])
#
# ellipsoid = Ellipsoid(center=center, radii=radii, rotmat=rotmat)
# start_cart = ellipsoid.jacobi2cartesian(coords=start_jacobi)
# #start_jacobi = ellipsoid.cartesian2jacobi(coords=start, norm=False)
# #start_test = ellipsoid.jacobi2cartesian(coords=start_jacobi)
#
# ellipsoid.get_geodesic_destination(start=start_cart.flatten(), alpha=alpha, distance=distance,
#                                                           n_steps=1000, method="Euler")
#
# n_steps_list = [10,50, 100,200,300,400,500,1000,10000]
# time_euler = np.zeros(len(n_steps_list))
# time_rk45 = np.zeros(len(n_steps_list))
# end_euler = np.zeros((len(n_steps_list), 3))
# end_rk45 = np.zeros((len(n_steps_list), 3))
#
# for i_step, n_steps in enumerate(n_steps_list):
#     # Euler
#     start_time = time.time()
#     end_euler[i_step, :] = ellipsoid.get_geodesic_destination(start=start_cart.flatten(), alpha=alpha, distance=distance,
#                                                         n_steps=n_steps, method="Euler")
#     stop_time = time.time()
#     time_euler[i_step] = stop_time-start_time
#
#     # RK45
#     start_time = time.time()
#     end_rk45[i_step, :] = ellipsoid.get_geodesic_destination(start=start_cart.flatten(), alpha=alpha, distance=distance,
#                                                         n_steps=n_steps, method="RK45")
#     stop_time = time.time()
#     time_rk45[i_step] = stop_time - start_time
#
#
#     print(f"N={n_steps} (Euler:{time_euler[i_step]}, RK45:{time_rk45[i_step]})")
#
# nrmsd_euler = (np.linalg.norm(end_euler[:-1, :]-end_euler[-1, :],axis=1)) / np.linalg.norm(end_euler[-1, :])
# nrmsd_rk45 = (np.linalg.norm(end_rk45[:-1, :]-end_rk45[-1, :],axis=1)) / np.linalg.norm(end_rk45[-1, :])
# nrmsd_euler_rk45 = (np.linalg.norm(end_rk45-end_euler[-1, :],axis=1)) / np.linalg.norm(end_euler[-1, :])
#
# plt.figure()
# plt.semilogx(n_steps_list, end_euler[:, 0], 'ro-')
# plt.semilogx(n_steps_list, end_rk45[:, 0], 'ro--')
# plt.semilogx(n_steps_list, end_euler[:, 1], 'go-')
# plt.semilogx(n_steps_list, end_rk45[:, 1], 'go--')
# plt.semilogx(n_steps_list, end_euler[:, 2], 'bo-')
# plt.semilogx(n_steps_list, end_rk45[:, 2], 'bo--')
# plt.legend(["x (Euler)", "x (RK45)", "y (Euler)", "y (RK45)", "z (Euler)", "z (RK45)"])
# plt.grid()
#
# plt.figure()
# plt.loglog(np.array(n_steps_list[:-1]), nrmsd_euler, 'o-')
# plt.grid()
#
# plt.figure()
# plt.loglog(np.array(n_steps_list[:-1]), nrmsd_rk45, 'o-')
# plt.grid()
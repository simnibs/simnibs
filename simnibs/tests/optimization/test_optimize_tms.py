import os
import numpy as np
import pytest
import scipy.spatial
from simnibs.optimization import optimize_tms
from simnibs.msh import mesh_io


@pytest.fixture
def sphere3_msh():
    fn = os.path.join(os.path.dirname(os.path.realpath(
        __file__)), '..', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn)


class TestGetOptGrid:
    def test_get_pos(self, sphere3_msh):
        pos, normals = optimize_tms._create_grid(
            sphere3_msh,
            [90, 0, 0],
            0, 10, 2
        )
        R = np.linalg.norm(pos - [94.60, -2.81, 5.79], axis=1)
        assert np.all(R <= 10)
        assert np.isclose(np.max(R), 10, atol=1.5)
        assert np.allclose(
            pos/np.linalg.norm(pos, axis=1)[:, None],
            normals, rtol=1e-1, atol=1e-1
        )
        dists = scipy.spatial.distance_matrix(pos, pos)
        dists += 1e6*np.eye(len(pos))
        assert np.allclose(
            np.min(dists, axis=1), 2, atol=1e-1
        )
    def test_rotate_system(self):
        R = np.array(
            [[1, 0, 0],
             [0, 0, -1],
             [0, 1, 0]]
        )
        matrices = optimize_tms._rotate_system(
            R, [-90, 90], 90
        )
        assert len(matrices) == 3
        # fix these
        assert np.allclose(
            matrices[0],
            np.array(
                [[0, 1, 0],
                 [0, 0, -1],
                 [-1, 0, 0]]
            )
        )
        assert np.allclose(
            matrices[1], R)
        assert np.allclose(
            matrices[2],
            np.array(
                [[0, -1, 0],
                 [0, 0, -1],
                 [1, 0, 0]]
            )
        )

    def test_get_opt_grid(self, sphere3_msh):
        matrices = optimize_tms.get_opt_grid(
            sphere3_msh, [90, 0, 0], handle_direction_ref=None,
            distance=1., radius=10, resolution_pos=1,
            resolution_angle=20, angle_limits=None)
        M = np.array(matrices)
        c = M[:, :3, 3]
        R = np.linalg.norm(c - [94.60, -2.81, 5.79], axis=1)
        assert np.all(R <= 10)
        assert np.isclose(np.max(R), 10, atol=1.5)
        n_pos = len(np.unique(R))
        n_angles = 360//20
        assert len(matrices) == n_pos * n_angles
        for i in range(n_pos):
            A = M[i*n_angles:((i + 1)*n_angles - 1), :3, :3]
            B = M[(i*n_angles + 1):(i + 1)*n_angles, :3, :3]
            AdotB = np.sum(A*B, axis=1)
            assert np.allclose(AdotB[:, 2], np.cos(0))
            assert np.allclose(AdotB[:, :2], np.cos(np.deg2rad(20)))
            assert np.allclose(
                np.sum(A[0, :, :2]*B[-1, :, :2], axis=0),
                np.cos(np.deg2rad(20))
            )

    def test_get_opt_handle_dir(self, sphere3_msh):
        matrices = optimize_tms.get_opt_grid(
            sphere3_msh, [90, 0, 0], handle_direction_ref=[90, 0, 1],
            distance=1., radius=10, resolution_pos=1,
            resolution_angle=30, angle_limits=[-60, 60])
        M = np.array(matrices)
        c = M[:, :3, 3]
        R = np.linalg.norm(c - [94.60, -2.81, 5.79], axis=1)
        assert np.all(R <= 10)
        assert np.isclose(np.max(R), 10, atol=1.5)
        n_pos = len(np.unique(R))
        n_angles = 120//30 + 1
        assert len(matrices) == n_pos * n_angles
        for i in range(n_pos):
            A = M[i*n_angles:((i + 1)*n_angles - 1), :3, :3]
            B = M[(i*n_angles + 1):(i + 1)*n_angles, :3, :3]
            AdotB = np.sum(A*B, axis=1)
            assert np.allclose(AdotB[:, 2], np.cos(0))
            assert np.allclose(AdotB[:, :2], np.cos(np.deg2rad(30)))
            assert np.allclose(
                np.sum(A[0, :, :2]*B[-1, :, :2], axis=0),
                np.cos(np.deg2rad(120))
            )
            assert np.allclose(
                np.sum(A[0, :, 1]*[0, 0, 1], axis=0),
                np.cos(np.deg2rad(60)),
                atol=1e-1
            )
            assert np.allclose(
                np.sum(B[-1, :, 1]*[0, 0, 1], axis=0),
                np.cos(np.deg2rad(60)),
                atol=1e-1
            )

def test_define_target_region(sphere3_msh):
    c = [85, 0, 0]
    r = 5
    t = 3
    elm = optimize_tms.define_target_region(sphere3_msh, c, r, t)
    outside = sphere3_msh.elm.elm_number[
        ~np.isin(sphere3_msh.elm.elm_number, elm)
    ]
    bar = sphere3_msh.elements_baricenters()
    dist = np.linalg.norm(bar[:] - c, axis=1)
    assert np.all(dist[elm - 1] < r)
    assert np.all(sphere3_msh.elm.tag1[elm - 1] == 3)
    assert np.all(sphere3_msh.elm.elm_type[elm - 1] == 4)

def test_get_opt_grid_ADM(sphere3_msh):
    matrices, coil_dir = optimize_tms.get_opt_grid_ADM(
        sphere3_msh, [90, 0, 0], handle_direction_ref=[90, 0, 1],
        distance=1., radius=10, resolution_pos=1,
        resolution_angle=30, angle_limits=[-60, 60]
    )
    assert np.allclose(np.linalg.det(matrices.transpose(2, 1, 0)), 1)
    centers = matrices[:3, 3, :]
    assert np.all(np.linalg.norm(centers.T - [94.60, -2.81, 5.79], axis=1) <= 10)
    normal_component = np.sum(
        matrices[:3, 2, :] *
        centers/np.linalg.norm(centers, axis=0),
        axis=0
    )
    assert np.all(normal_component < -0.9)
    handle_dir = np.sum(
        matrices[:3, 1, :].T * [0, 0, 1], axis=1
    )
    assert np.all(handle_dir > 0.9)

    assert np.allclose(
        np.rad2deg(np.arctan2(coil_dir[0, :], coil_dir[1, :])),
        [-60., -30., 0., 30., 60.]
    )

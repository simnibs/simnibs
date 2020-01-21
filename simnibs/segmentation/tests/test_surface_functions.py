import os

import numpy as np
import scipy.io
import pytest


from simnibs import SIMNIBSDIR, mesh_io
from simnibs.segmentation import surface_functions

@pytest.fixture
def sphere_surf():
    fn = os.path.join(
            SIMNIBSDIR, 'resources', 'testing_files', 'sphere3.msh')
    mesh = mesh_io.read_msh(fn)
    return mesh.crop_mesh(1005)



class TestRayTriangleIntersect:
    def test_ray_triangle_intersect(self, sphere_surf):
        bar = sphere_surf.elements_baricenters()[:]
        normals = sphere_surf.triangle_normals()[:]
        intersect, int_points, int_dist = \
            surface_functions.ray_triangle_intersect(
                    sphere_surf.nodes[sphere_surf.elm[:, :3]],
                    (bar[0] + normals[0])[None, :],
                    -normals[None, 0, :],
                    return_int_points=True,
                    return_int_dists=True
                )
        assert intersect[0, 0]
        assert np.allclose(int_points[0, 0], bar[0])
        assert np.allclose(int_dist[0, 0], 1.)

        assert np.linalg.norm(int_points[intersect][1] + bar[0]) < 5
        assert np.allclose(int_dist[intersect][1], 190, rtol=1e-2)

    def test_ray_triangle_no_intersect(self, sphere_surf):
        intersect, int_points, int_dist = \
            surface_functions.ray_triangle_intersect(
                    sphere_surf.nodes[sphere_surf.elm[:, :3]],
                    np.array([[100, 0, 0]]),
                    np.array([[0, 1, 0]]),
                    return_int_points=True,
                    return_int_dists=True
                )

        assert not np.any(intersect)


class TestExpandCS:
    def test_expand(self, sphere_surf):
        vertices = sphere_surf.nodes[:]
        faces = sphere_surf.elm[:, :3] - 1
        vertices_e = surface_functions.expandCS(
            vertices, faces, 5*np.ones(len(vertices)))
        assert np.allclose(np.linalg.norm(vertices_e, axis=1), 100, atol=1)

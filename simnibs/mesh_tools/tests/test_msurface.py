import os

import pytest
import numpy as np

from ... import SIMNIBSDIR
from .. import mesh_io, surface


@pytest.fixture(scope='module')
def sphere3_msh():
    fn = os.path.join(
        SIMNIBSDIR, '_internal', 'testing_files', 'sphere3.msh')

    return mesh_io.read_msh(fn)


@pytest.fixture
def sphere3_surf(sphere3_msh):
    return surface.Surface(sphere3_msh, 1005)


class TestSurfaceSetup:
    def test_nodes(self, sphere3_surf, sphere3_msh):
        """Test if all nodes are in the surface """
        R = np.linalg.norm(sphere3_surf.nodes, axis=1)
        np.testing.assert_allclose(R, 95 * np.ones(R.shape), rtol=0.05)

    def test_triangle_baricenters(self, sphere3_surf, sphere3_msh):
        R = np.linalg.norm(sphere3_surf.tr_centers, axis=1)
        np.testing.assert_allclose(R, 95 * np.ones(R.shape), rtol=0.05)

    def test_triangle_normals(self, sphere3_surf, sphere3_msh):
        cos = np.einsum(
            'ij,ij->i', sphere3_surf.tr_centers / 95, sphere3_surf.tr_normals)
        np.testing.assert_allclose(cos, np.ones(cos.shape), rtol=0.05)

    def test_triangle_nodes(self, sphere3_surf, sphere3_msh):
        b = np.average(sphere3_surf.nodes[sphere3_surf.tr_nodes[0]], axis=0)
        np.testing.assert_allclose(b, sphere3_surf.tr_centers[0])

    def test_triangle_dict(self, sphere3_surf, sphere3_msh):
        assert np.all(sphere3_surf.surf2msh_triangles == np.where(
            sphere3_msh.elm.tag1 == 1005)[0])

    def test_tiangle_area(self, sphere3_surf):
        np.testing.assert_allclose(
            np.sum(sphere3_surf.tr_areas), 4 * np.pi * 95**2, rtol=0.1)

    def test_node_normal(self, sphere3_surf):
        cos = np.einsum(
            'ij,ij->i', sphere3_surf.nodes / 95, sphere3_surf.nodes_normals)
        np.testing.assert_allclose(cos, np.ones(cos.shape), rtol=0.05)

    def test_node_area(self, sphere3_surf):
        triangles_with_node = np.where(sphere3_surf.tr_nodes == 0)[0]

        area_of_triangles_w_node = np.sum(
            sphere3_surf.tr_areas[triangles_with_node])

        np.testing.assert_allclose(
            sphere3_surf.nodes_areas[0], 1. / 3. * area_of_triangles_w_node)


if __name__ == '__main__':
    msh = sphere3_msh()
    s = sphere3_surf(msh)
    import IPython
    IPython.embed()

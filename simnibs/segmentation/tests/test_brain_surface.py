import os

import numpy as np
import scipy.io
import pytest


from simnibs import SIMNIBSDIR, mesh_io
from simnibs.segmentation import brain_surface

@pytest.fixture
def sphere_surf():
    fn = os.path.join(
            SIMNIBSDIR, 'resources', 'testing_files', 'sphere3.msh')
    mesh = mesh_io.read_msh(fn)
    return mesh.crop_mesh(1005)

@pytest.fixture
def sphere_2_surfs():
    fn = os.path.join(
            SIMNIBSDIR, 'resources', 'testing_files', 'sphere3.msh')
    mesh = mesh_io.read_msh(fn)
    return mesh.crop_mesh([1004, 1005])


class TestSegmentTriangleIntersect:
    def test_segment_triangle_intersect(self, sphere_surf):
        bar = sphere_surf.elements_baricenters()[:]
        normals = sphere_surf.triangle_normals()[:]
        intersect, int_points = \
            brain_surface.segment_triangle_intersect(
                    sphere_surf.nodes[:], sphere_surf.elm[:, :3] - 1,
                    bar - 1e-1 * normals,
                    bar + 1e-1 * normals,
                )
        assert np.all(intersect[:, 0] == intersect[:, 1])
        assert np.allclose(int_points, bar)

    def test_segment_triangle_no_intersect(self, sphere_surf):
        bar = sphere_surf.elements_baricenters()[:]
        normals = sphere_surf.triangle_normals()[:]
        intersect, int_points = \
            brain_surface.segment_triangle_intersect(
                    sphere_surf.nodes[:], sphere_surf.elm[:, :3] - 1,
                    bar + 1e-1 * normals,
                    bar + 2e-1 * normals,
                )

        assert len(intersect) == 0


class TestExpandCS:
    def test_expand(self, sphere_surf):
        vertices = sphere_surf.nodes[:]
        faces = sphere_surf.elm[:, :3] - 1
        vertices_e = brain_surface.expandCS(
            vertices, faces, 5*np.ones(len(vertices)))
        assert np.allclose(np.linalg.norm(vertices_e, axis=1), 100, atol=1)

    def test_reduce(self, sphere_surf):
        vertices = sphere_surf.nodes[:]
        faces = sphere_surf.elm[:, :3] - 1
        vertices_e = brain_surface.expandCS(
            vertices, faces, 2*np.ones(len(vertices)), deform="shrink")
        assert np.allclose(np.linalg.norm(vertices_e, axis=1), 92, atol=3)

    def test_no_move(self, sphere_surf):
        vertices = sphere_surf.nodes[:]
        faces = sphere_surf.elm[:, :3] - 1
        vertices_e = brain_surface.expandCS(vertices, faces, np.zeros(len(vertices)))
        assert np.allclose(np.linalg.norm(vertices_e, axis=1), 95, atol=1)


    def test_2_surfs(self, sphere_2_surfs):
        vertices = sphere_2_surfs.nodes[:]
        faces = sphere_2_surfs.elm[:, :3] - 1
        nodes_surf1 = np.unique(
            sphere_2_surfs.elm[sphere_2_surfs.elm.tag1 == 1004, :3]) - 1
        nodes_surf2 = np.unique(
            sphere_2_surfs.elm[sphere_2_surfs.elm.tag1 == 1005, :3]) - 1
        shift = np.ones(len(vertices))
        shift[nodes_surf1] = 2.
        shift[nodes_surf2] = 4.
        vertices_e = brain_surface.expandCS(vertices, faces, shift)
        assert np.allclose(np.linalg.norm(vertices_e[nodes_surf1], axis=1), 92, atol=1)
        assert np.allclose(np.linalg.norm(vertices_e[nodes_surf2], axis=1), 99, atol=1)

    def test_2_surfs_colision(self, sphere_2_surfs):
        vertices = sphere_2_surfs.nodes[:]
        faces = sphere_2_surfs.elm[:, :3] - 1
        nodes_surf1 = np.unique(
            sphere_2_surfs.elm[sphere_2_surfs.elm.tag1 == 1004, :3]) - 1
        nodes_surf2 = np.unique(
            sphere_2_surfs.elm[sphere_2_surfs.elm.tag1 == 1005, :3]) - 1
        # rotate nodes for numerical reasons
        angle = np.pi/4
        rot = np.array([[np.cos(angle), -np.sin(angle), 0.],
                        [np.sin(angle), np.cos(angle), 0.],
                        [0, 0, 1]])
        vertices[nodes_surf2] = vertices[nodes_surf2].dot(rot)
        shift = np.ones(len(vertices))
        shift[nodes_surf1] = 10.
        shift[nodes_surf2] = 0.
        vertices_e = brain_surface.expandCS(
            vertices, faces, shift, ensure_distance=0.5, nsteps=5)
        sphere_2_surfs.nodes.node_coord = vertices_e
        assert np.allclose(np.linalg.norm(vertices_e[nodes_surf2], axis=1), 95, atol=1)
        assert np.allclose(np.linalg.norm(vertices_e[nodes_surf1], axis=1), 93, atol=2)

class TestCreateSurfaceMask:

    @pytest.mark.parametrize('axis', ['z', 'y', 'x'])
    def test_rasterize_surface(self, axis, sphere_surf):
        vertices = sphere_surf.nodes[:]
        faces = sphere_surf.elm[:, :3] - 1
        affine = np.eye(4)
        affine[:3, 3] = -100
        shape = [200, 200, 200]
        mask = brain_surface._rasterize_surface(vertices, faces, affine, shape, axis=axis)
        assert np.all(mask[100, 100, 5:195])
        assert np.all(mask[100, 5:195, 100])
        assert np.all(mask[5:195, 100, 100])
        assert not np.any(mask[100, 100, :5]) and not np.any(mask[100, 100, 195:])
        assert not np.any(mask[100, :5, 100]) and not np.any(mask[100, 195:, 100])
        assert not np.any(mask[:5, 100, 100]) and not np.any(mask[195:, 100, 100])

    @pytest.mark.parametrize('axis', ['z', 'y', 'x'])
    def test_rasterize_surface_concave(self, axis, sphere_2_surfs):
        vertices = sphere_2_surfs.nodes[:]
        faces = sphere_2_surfs.elm[:, :3] - 1
        affine = np.eye(4)
        affine[:3, 3] = -100
        shape = [200, 200, 200]
        mask = brain_surface._rasterize_surface(vertices, faces, affine, shape, axis=axis)
        assert np.all(mask[100, 100, 5:10]) and np.all(mask[100, 100, 190:195])
        assert np.all(mask[100, 5:10, 100]) and np.all(mask[100, 190:195, 100])
        assert np.all(mask[5:10, 100, 100]) and np.all(mask[190:195, 100, 100])
        assert not np.any(mask[100, 100, 10:190])
        assert not np.any(mask[100, 10:190, 100])
        assert not np.any(mask[10:190, 100, 100])

    @pytest.mark.parametrize('axis', ['z', 'y', 'x'])
    def test_rasterize_surface_aniso(self, axis, sphere_surf):
        vertices = sphere_surf.nodes[:]
        faces = sphere_surf.elm[:, :3] - 1
        affine = np.diag([1., 1., 2., 1])
        affine[:3, 3] = -100
        shape = [200, 200, 100]
        mask = brain_surface._rasterize_surface(vertices, faces, affine, shape, axis=axis)
        assert np.all(mask[100, 100, 3:97])
        assert np.all(mask[100, 5:195, 50])
        assert np.all(mask[5:195, 100, 50])
        assert not np.any(mask[100, 100, :3]) and not np.any(mask[100, 100, 97:])
        assert not np.any(mask[100, :5, 50]) and not np.any(mask[100, 195:, 50])
        assert not np.any(mask[:5, 100, 50]) and not np.any(mask[195:, 100, 50])


    @pytest.mark.parametrize('axis', ['z', 'y', 'x'])
    def test_rasterize_surface_slice(self, axis, sphere_surf):
        vertices = sphere_surf.nodes[:]
        faces = sphere_surf.elm[:, :3] - 1
        affine = np.eye(4)
        affine[:3, 3] = [-100, -100, 0]
        shape = [200, 200, 1]
        mask = brain_surface._rasterize_surface(vertices, faces, affine, shape, axis=axis)
        assert np.all(mask[100, 5:195, 0])
        assert np.all(mask[5:195, 100, 0])
        assert not np.any(mask[100, :5, 0]) and not np.any(mask[100, 195:, 0])
        assert not np.any(mask[:5, 100, 0]) and not np.any(mask[195:, 100, 0])

    def test_rasterize_surface_bad_topology(self, sphere_surf):
        vertices = sphere_surf.nodes[:]
        faces = sphere_surf.elm[:, :3] - 1
        # Clone part of the faces to create an "outer shell" in part of the sphere
        faces_bar = np.average(vertices[faces], axis=1)
        faces_to_clone = faces[faces_bar[:, 2] > 50]
        verts_to_clone, reverse = np.unique(faces_to_clone, return_inverse=True)
        verts_cloned = vertices[verts_to_clone] + [0, 0, 2]
        faces_cloned = reverse.reshape(-1, 3) + len(vertices)
        faces = np.vstack([faces, faces_cloned])
        vertices = np.vstack([vertices, verts_cloned])
        affine = np.eye(4)
        affine[:3, 3] = -100
        shape = [200, 200, 200]
        mask = brain_surface._rasterize_surface(vertices, faces, affine, shape)
        assert np.all(mask[100, 100, 5:195])
        assert np.all(mask[100, 5:195, 100])
        assert np.all(mask[5:195, 100, 100])
        assert not np.any(mask[100, 100, :5]) and not np.any(mask[100, 100, 195:])
        assert not np.any(mask[100, :5, 100]) and not np.any(mask[100, 195:, 100])
        assert not np.any(mask[:5, 100, 100]) and not np.any(mask[195:, 100, 100])

    def test_mask_from_surface(self, sphere_surf):
        vertices = sphere_surf.nodes[:]
        faces = sphere_surf.elm[:, :3] - 1
        affine = np.eye(4)
        affine[:3, 3] = -100
        shape = [200, 200, 200]
        mask = brain_surface.mask_from_surface(vertices, faces, affine, shape)
        assert np.all(mask[100, 100, 5:195])
        assert np.all(mask[100, 5:195, 100])
        assert np.all(mask[5:195, 100, 100])
        assert not np.any(mask[100, 100, :5]) and not np.any(mask[100, 100, 195:])
        assert not np.any(mask[100, :5, 100]) and not np.any(mask[100, 195:, 100])
        assert not np.any(mask[:5, 100, 100]) and not np.any(mask[195:, 100, 100])

    def test_mask_from_surface_concave(self, sphere_2_surfs):
        vertices = sphere_2_surfs.nodes[:]
        faces = sphere_2_surfs.elm[:, :3] - 1
        affine = np.eye(4)
        affine[:3, 3] = -100
        shape = [200, 200, 200]
        mask = brain_surface.mask_from_surface(vertices, faces, affine, shape)
        assert np.all(mask[100, 100, 5:10]) and np.all(mask[100, 100, 190:195])
        assert np.all(mask[100, 5:10, 100]) and np.all(mask[100, 190:195, 100])
        assert np.all(mask[5:10, 100, 100]) and np.all(mask[190:195, 100, 100])
        assert not np.any(mask[100, 100, 10:190])
        assert not np.any(mask[100, 10:190, 100])
        assert not np.any(mask[10:190, 100, 100])

    def test_mask_from_surface_bad_topology(self, sphere_surf):
        vertices = sphere_surf.nodes[:]
        faces = sphere_surf.elm[:, :3] - 1
        # Clone part of the faces to create an "outer shell" in part of the sphere
        faces_bar = np.average(vertices[faces], axis=1)
        faces_to_clone = faces[faces_bar[:, 2] > 93]
        verts_to_clone, reverse = np.unique(faces_to_clone, return_inverse=True)
        verts_cloned = vertices[verts_to_clone] + [0, 0, 3]
        faces_cloned = reverse.reshape(-1, 3) + len(vertices)
        faces = np.vstack([faces, faces_cloned])
        vertices = np.vstack([vertices, verts_cloned])
        affine = np.eye(4)
        affine[:3, 3] = -100
        shape = [200, 200, 200]
        mask = brain_surface.mask_from_surface(vertices, faces, affine, shape)
        assert np.all(mask[100, 100, 5:195])
        assert np.all(mask[100, 5:195, 100])
        assert np.all(mask[5:195, 100, 100])
        assert not np.any(mask[100, :5, 100]) and not np.any(mask[100, 195:, 100])
        assert not np.any(mask[:5, 100, 100]) and not np.any(mask[195:, 100, 100])

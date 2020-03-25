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


def one_sphere_image():
    """
    create image with 1 sphere in center (WM is 3, GM is 2, CSF is 1)

    Returns
    -------
    img : 3d image
    vox2mm : affine trafo
    radius_wm : radius of WM-GM boundary
    radius_gm : radius of GM-CSF boundary
    radius_csf : radius of CSF-background boundary

    """  

    imgsize=[161, 141, 141]
    voxsize=0.25
    
    radius_wm=8
    radius_gm=12
    radius_csf=15
    
    vox2mm=np.float32(([voxsize, 0, 0, -voxsize*(imgsize[0]-1)/2],
                       [0, voxsize, 0, -voxsize*(imgsize[1]-1)/2],
                       [0, 0, voxsize, -voxsize*(imgsize[2]-1)/2],
                       [0, 0, 0, 1]))
    
    G = np.meshgrid(np.arange(imgsize[0]), np.arange(imgsize[1]), 
                    np.arange(imgsize[2]), indexing='ij')
    
    xyz_mm=np.inner(vox2mm, np.vstack((G[0].flatten(),
                                       G[1].flatten(),
                                       G[2].flatten(),
                                       np.ones_like(G[0].flatten()))).T
                    )[0:3]
    
    R_sp = np.linalg.norm(xyz_mm, axis=0)
    
    img=np.zeros_like(R_sp)
    img[R_sp<radius_csf]=1
    img[R_sp<radius_gm]=2
    img[R_sp<radius_wm]=3
    img=img.reshape(imgsize)

    return img, vox2mm, radius_wm, radius_gm, radius_csf


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
        affine[:3, 3] = -99.9 # I need to have a broken number here because the mesh has
        # nodes aligned with the x, y, and z axis
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
        affine[:3, 3] = -99.9
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
        affine[:3, 3] = -99.9
        shape = [200, 200, 100]
        mask = brain_surface._rasterize_surface(vertices, faces, affine, shape, axis=axis)
        assert np.all(mask[100, 100, 3:97])
        assert np.all(mask[100, 5:195, 50])
        assert np.all(mask[5:195, 100, 50])
        assert not np.any(mask[100, 100, :3]) and not np.any(mask[100, 100, 98:])
        assert not np.any(mask[100, :5, 50]) and not np.any(mask[100, 195:, 50])
        assert not np.any(mask[:5, 100, 50]) and not np.any(mask[195:, 100, 50])


    @pytest.mark.parametrize('axis', ['z', 'y', 'x'])
    def test_rasterize_surface_slice(self, axis, sphere_surf):
        vertices = sphere_surf.nodes[:]
        faces = sphere_surf.elm[:, :3] - 1
        affine = np.eye(4)
        affine[:3, 3] = [-99.9, -99.9, 0.1]
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
        affine[:3, 3] = -99.9
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
        affine[:3, 3] = -99.9
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
        affine[:3, 3] = -99.9
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
        affine[:3, 3] = -99.9
        shape = [200, 200, 200]
        mask = brain_surface.mask_from_surface(vertices, faces, affine, shape)
        assert np.all(mask[100, 100, 5:195])
        assert np.all(mask[100, 5:195, 100])
        assert np.all(mask[5:195, 100, 100])
        assert not np.any(mask[100, :5, 100]) and not np.any(mask[100, 195:, 100])
        assert not np.any(mask[:5, 100, 100]) and not np.any(mask[195:, 100, 100])


class Test_cat_vol_pbt_AT:
    def test_with_spheres(self):
        testdata=one_sphere_image()
        img=testdata[0]
        vox2mm=testdata[1]
        
        Yth, Ypp = brain_surface.cat_vol_pbt_AT(img, 0.25, False)
    
        # test for 2.5% accuracy in volume-average of GM thickness and percent position profile
        assert abs(1-Yth[img==2].mean()/4) < 0.025
        assert abs(1-Ypp[img==2].mean()/0.4342) < 0.025 
        # volume average of profile varying lineary from 1 to 0 
        # should be 0.4342... for r_wm=12, r_gm=8
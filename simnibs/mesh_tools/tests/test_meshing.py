import copy
import os
import tempfile
import pytest
import numpy as np

from ... import SIMNIBSDIR
from .. import meshing, mesh_io


@pytest.fixture
def labeled_image():
    img = np.zeros((50, 50, 50), dtype=np.uint8)
    img[10:40, 15:35, 20:30] = 1
    img[20:30, 20:30, 20:30] = 2
    return img


@pytest.fixture
def cube_image():
    img = np.zeros((50, 50, 50), dtype=np.float)
    img[10:40, 10:40, 10:40] = 1
    return img


@pytest.fixture
def surface():
    fn = os.path.join(
            SIMNIBSDIR, '_internal', 'testing_files', 'cube.off')
    return mesh_io.read_off(fn)

@pytest.fixture
def sphere3():
    fn = os.path.join(
            SIMNIBSDIR, '_internal', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn)


def volumes(mesh):
    vols = mesh.elements_volumes_and_areas()
    vols = [
        np.sum(vols[(mesh.elm.tag1 == t) * (mesh.elm.elm_type == 4)])
        for t in np.unique(mesh.elm.tag1)
    ]
    return vols

def test_write_inr(labeled_image):
    fn = os.path.join(tempfile.gettempdir(), tempfile.gettempprefix() + '.inr')
    meshing._write_inr(labeled_image, [1, 1, 1], fn)
    assert os.path.exists(fn)
    os.remove(fn)

class Test_decompose_affine():
    def test_diagonal(self):
        R, scaling, shearing = meshing._decompose_affine(np.eye(4))
        assert np.allclose(R, np.eye(3))
        assert np.allclose(scaling, 1)
        assert np.allclose(shearing, np.eye(3))
        R, scaling,  shearing = meshing._decompose_affine(2 * np.eye(4))
        assert np.allclose(R, np.eye(3))
        assert np.allclose(scaling, 2)
        assert np.allclose(shearing, np.eye(3))

    def test_off_diagonal(self):
        aff = np.array([
            [0, 1, 0, 0],
            [2, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        R, scaling, shearing = meshing._decompose_affine(aff)
        assert np.allclose(
            R,
            np.array([
                [0, 1, 0],
                [1, 0, 0],
                [0, 0, 1],
            ]))
        assert np.allclose(scaling, [2, 1, 1])
        assert np.allclose(shearing, np.eye(3))

    def test_partial(self, labeled_image):
        theta = np.pi/3
        R = np.array([
            [np.cos(theta), -np.sin(theta), 0, 0],
            [np.sin(theta), np.cos(theta), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        scaling = np.diag([1, 2, 3, 1])
        affine = R.dot(scaling)
        rot, scaling, shearing = meshing._decompose_affine(affine)
        assert np.allclose(R[:3, :3], rot)
        assert np.allclose(scaling, [1, 2, 3])
        assert np.allclose(shearing, np.eye(3))

    def test_shearing(self, labeled_image):
        theta = np.pi/3
        R = np.array([
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta), np.cos(theta), 0],
            [0, 0, 1],
        ])
        Z = np.diag([1, 2, 3])
        S = np.array([
            [1, 1.2, 0],
            [0, 1, 0],
            [0, 0, 1],
        ])
        affine = R.dot(Z.dot(S))
        rot, scaling, shearing = meshing._decompose_affine(affine)
        assert np.allclose(R, rot)
        assert np.allclose(scaling, np.diagonal(Z))
        assert np.allclose(shearing, S)


class TestResample2Iso():
    def test_resample(self, cube_image):
        aff = np.eye(4)
        #nibabel.viewers.OrthoSlicer3D(cube_image, aff).show()
        iso, aff2 = meshing._resample2iso(cube_image, aff, .5, order=0)
        #nibabel.viewers.OrthoSlicer3D(iso, aff2).show()
        mask = np.zeros((100, 100, 100), dtype=bool)
        mask[19:79, 19:79, 19:79] = True
        assert iso.shape == (100, 100, 100)
        assert np.allclose(aff2, np.diag([.5, .5, .5, 1]))
        assert np.all(iso[mask] > 0)
        assert np.allclose(iso[~mask], 0)


    def test_diagonal(self, cube_image):
        aff = np.diag([1, 2, 3, 1])
        #nibabel.viewers.OrthoSlicer3D(cube_image, aff).show()
        iso, aff2 = meshing._resample2iso(cube_image, aff, 1, order=0)
        #nibabel.viewers.OrthoSlicer3D(iso, aff2).show()
        mask = np.zeros((50, 100, 150), dtype=bool)
        mask[10:40, 19:79, 29:119] = True
        assert iso.shape == (50, 100, 150)
        assert np.allclose(aff2, np.eye(4))
        assert np.all(iso[mask] > 0)
        assert np.allclose(iso[~mask], 0)

    def test_rotation(self, cube_image):
        aff = np.array([
            [0, 1, 0, 0],
            [2, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        iso, aff2 = meshing._resample2iso(cube_image, aff, 1, order=0)
        mask = np.zeros((50, 100, 50), dtype=bool)
        mask[10:40, 19:79, 10:40] = True
        assert iso.shape == (50, 100, 50)
        assert np.allclose(
            aff2,
            np.array([
                [0, 1, 0, 0],
                [1, 0, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
            ]))
        assert np.all(iso[mask] > 0)
        assert np.allclose(iso[~mask], 0)


class TestImage2mesh():
    def test_diagonal(self, labeled_image):
        affine = np.eye(4)
        mesh = meshing.image2mesh(
            labeled_image, affine, facet_distance=0.5, facet_size=5, cell_size=5
        )
        assert np.allclose(np.min(mesh.nodes[:], axis=0), [9.5, 14.5, 19.5], rtol=1e-2)
        assert np.allclose(np.max(mesh.nodes[:], axis=0), [39.5, 34.5, 29.5], rtol=1e-2)
        vol_1, vol_2 = volumes(mesh)
        assert np.isclose(vol_2, 1e3, rtol=1e-1)
        assert np.isclose(vol_1, (30 * 20 * 10) - 1e3, rtol=1e-1)

    def test_scaling(self, labeled_image):
        affine = 2*np.eye(4)
        mesh = meshing.image2mesh(
            labeled_image, affine, facet_distance=1, facet_size=5, cell_size=5
        )
        assert np.allclose(np.min(mesh.nodes[:], axis=0), [2*9.5, 2*14.5, 2*19.5], rtol=1e-2)
        assert np.allclose(np.max(mesh.nodes[:], axis=0), [2*39.5, 2*34.5, 2*29.5], rtol=1e-2)
        vol_1, vol_2 = volumes(mesh)
        assert np.isclose(vol_2, 8e3, rtol=1e-1)
        assert np.isclose(vol_1, 8*((30 * 20 * 10) - 1e3), rtol=1e-1)


    def test_translation(self, labeled_image):
        affine = np.eye(4)
        affine[:3, 3] = -25
        mesh = meshing.image2mesh(
            labeled_image, affine, facet_distance=0.5, facet_size=5, cell_size=5
        )
        assert np.allclose(np.min(mesh.nodes[:], axis=0), [9.5-25, 14.5-25, 19.5-25], rtol=1e-2)
        assert np.allclose(np.max(mesh.nodes[:], axis=0), [39.5-25, 34.5-25, 29.5-25], rtol=1e-2)
        vol_1, vol_2 = volumes(mesh)
        assert np.isclose(vol_2, 1e3, rtol=1e-1)
        assert np.isclose(vol_1, (30 * 20 * 10) - 1e3, rtol=1e-1)

    def test_rotation_90(self, labeled_image):
        affine = np.array([
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        mesh = meshing.image2mesh(
            labeled_image, affine, facet_distance=0.5, facet_size=5, cell_size=5
        )
        assert np.allclose(np.min(mesh.nodes[:], axis=0), [14.5, 9.5, 19.5], rtol=1e-2)
        assert np.allclose(np.max(mesh.nodes[:], axis=0), [34.5, 39.5, 29.5], rtol=1e-2)
        vol_1, vol_2 = volumes(mesh)
        assert np.isclose(vol_2, 1e3, rtol=1e-1)
        assert np.isclose(vol_1, (30 * 20 * 10) - 1e3, rtol=1e-1)

    def test_rotation(self, labeled_image):
        theta = np.pi/4
        R = np.array([
            [np.cos(theta), -np.sin(theta), 0, 0],
            [np.sin(theta), np.cos(theta), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        scaling = np.diag([1, 2, 3, 1])
        affine = R.dot(scaling)
        mesh = meshing.image2mesh(
            labeled_image, affine, facet_distance=0.5, facet_size=5, cell_size=5
        )
        # transform nodes back to voxel space
        nodes_t = np.linalg.inv(affine[:3, :3]).dot(mesh.nodes[:].T).T
        assert np.allclose(np.min(nodes_t, axis=0), [9.5, 14.5, 19.5], rtol=1e-2)
        assert np.allclose(np.max(nodes_t, axis=0), [39.5, 34.5, 29.5], rtol=1e-2)
        mesh.nodes.node_coord = nodes_t
        vol_1, vol_2 = volumes(mesh)
        assert np.isclose(vol_2, 1e3, rtol=1e-1)
        assert np.isclose(vol_1, (30 * 20 * 10) - 1e3, rtol=1e-1)

    @pytest.mark.parametrize('axis', [0, 1, 2])
    def test_sizing_field(self, axis):
        label_image = np.zeros((50, 50, 50), dtype=np.uint8)
        label_image[10:40, 10:40, 10:40] = 1
        affine = np.eye(4)
        sizing_field = np.ones_like(label_image) * 10
        if axis == 0:
            sizing_field[25:, ...] = 3
        elif axis == 1:
            sizing_field[:, 25:, :] = 3
        elif axis == 2:
            sizing_field[..., 25:] = 3
        mesh = meshing.image2mesh(
            label_image, affine, facet_distance=0.5, facet_size=10,
            cell_size=sizing_field
        )
        assert np.allclose(np.min(mesh.nodes[:], axis=0), 9.5, rtol=1e-2)
        assert np.allclose(np.max(mesh.nodes[:], axis=0), 39.5, rtol=1e-2)
        vols = mesh.elements_volumes_and_areas()[mesh.elm.tetrahedra]
        bar = mesh.elements_baricenters()[mesh.elm.tetrahedra]
        assert np.average(vols[bar[:, axis] < 25]) > np.average(vols[bar[:, axis] > 25])

    def test_sizing_field_facet_size(self):
        label_image = np.zeros((50, 50, 50), dtype=np.uint8)
        label_image[10:40, 10:40, 10:40] = 1
        affine = np.eye(4)
        sizing_field = np.ones_like(label_image) * 10
        sizing_field[25:, ...] = 1
        mesh = meshing.image2mesh(
            label_image, affine, facet_distance=0.5, facet_size=sizing_field, cell_size=10
        )
        assert np.allclose(np.min(mesh.nodes[:], axis=0), 9.5, rtol=1e-2)
        assert np.allclose(np.max(mesh.nodes[:], axis=0), 39.5, rtol=1e-2)
        areas = mesh.elements_volumes_and_areas()[mesh.elm.triangles]
        bar = mesh.elements_baricenters()[mesh.elm.triangles]
        assert np.average(areas[bar[:, 0] < 25]) > np.average(areas[bar[:, 0] > 25])

    def test_sizing_field_facet_distance(self):
        label_image = np.zeros((50, 50, 50), dtype=np.uint8)
        label_image[10:40, 10:40, 10:40] = 1
        affine = np.eye(4)
        sizing_field = np.ones_like(label_image, dtype=np.float32) * .5
        sizing_field[25:, ...] = .3
        mesh = meshing.image2mesh(
            label_image, affine, facet_size=10, cell_size=10,
            facet_distance=sizing_field
        )
        assert np.allclose(np.min(mesh.nodes[:], axis=0), 9.5, rtol=1e-2)
        assert np.allclose(np.max(mesh.nodes[:], axis=0), 39.5, rtol=1e-2)
        areas = mesh.elements_volumes_and_areas()[mesh.elm.triangles]
        bar = mesh.elements_baricenters()[mesh.elm.triangles]
        assert np.average(areas[bar[:, 0] < 25]) > np.average(areas[bar[:, 0] > 25])

    @pytest.mark.parametrize('axis', [0, 1, 2])
    def test_sizing_field_scaling(self, axis):
        label_image = np.zeros((50, 50, 50), dtype=np.uint8)
        label_image[10:40, 10:40, 10:40] = 1
        affine = np.eye(4)
        affine[axis, axis] = 2
        sizing_field = np.ones_like(label_image) * 10
        if axis == 0:
            sizing_field[25:, ...] = 3
        elif axis == 1:
            sizing_field[:, 25:, :] = 3
        elif axis == 2:
            sizing_field[..., 25:] = 3
        mesh = meshing.image2mesh(
            label_image, affine, facet_size=10, cell_size=sizing_field,
            facet_distance=0.5
        )
        vols = mesh.elements_volumes_and_areas()[mesh.elm.tetrahedra]
        bar = mesh.elements_baricenters()[mesh.elm.tetrahedra]
        assert np.average(vols[bar[:, axis] < 50]) > np.average(vols[bar[:, axis] > 50])


def test_mesh_surfaces(surface):
    surface2 = copy.deepcopy(surface)
    surface2.nodes.node_coord *= 2
    mesh = meshing._mesh_surfaces(
        [surface, surface2], [[1, 2], [2, 0]], 30,
        1, 0.1, 2, 1, optimize=True
    )
    assert np.allclose(np.min(mesh.nodes[:], axis=0), -2, rtol=1e-2)
    assert np.allclose(np.max(mesh.nodes[:], axis=0), 2, rtol=1e-2)
    vol_1, vol_2 = volumes(mesh)
    assert np.isclose(vol_2, 4**3 - 2**3, rtol=1e-1)
    assert np.isclose(vol_1, 2**3, rtol=1e-1)

def test_remesh(sphere3):
    mesh = meshing.remesh(sphere3, 10, 10, facet_distance=1, optimize=False)
    assert np.allclose(np.min(mesh.nodes[:], axis=0), -95, rtol=1e-2)
    assert np.allclose(np.max(mesh.nodes[:], axis=0), 95, rtol=1e-2)
    vols = volumes(mesh)
    assert np.isclose(vols[0], 4/3*np.pi*85**3, rtol=1e-1)
    assert np.isclose(vols[1], 4/3*np.pi*(90**3 - 85**3), rtol=1e-1)
    assert np.isclose(vols[2], 4/3*np.pi*(95**3 - 90**3), rtol=1e-1)

class TestRelabelSpiles:
    def test_relabel_spikes(self, sphere3):
        mesh = copy.deepcopy(sphere3)
        mesh.elm.tag1[9033] = 3
        mesh.elm.tag2[9033] = 3
        mesh.elm.tag1[8343] = 5
        mesh.elm.tag2[8343] = 5
        meshing.relabel_spikes(
            mesh, 3, 5, 4
        )
        assert np.all(mesh.elm.tag1 == sphere3.elm.tag1)
        assert np.all(mesh.elm.tag2 == sphere3.elm.tag2)

    def test_despike(self, sphere3):
        mesh = copy.deepcopy(sphere3)
        mesh.elm.tag1[9033] = 3
        mesh.elm.tag2[9033] = 3
        mesh.elm.tag1[8343] = 5
        mesh.elm.tag2[8343] = 5
        meshing.despike(mesh, adj_threshold=3)
        assert np.all(mesh.elm.tag1 == sphere3.elm.tag1)
        assert np.all(mesh.elm.tag2 == sphere3.elm.tag2)

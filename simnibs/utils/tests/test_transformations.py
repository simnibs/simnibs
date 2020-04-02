import os

import numpy as np
import pytest

from ... import SIMNIBSDIR
from ...mesh_tools import mesh_io
from .. import transformations


@pytest.fixture
def sphere3_msh():
    fn = os.path.join(SIMNIBSDIR, 'resources', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn)

@pytest.fixture
def image():
    x, y, z = np.meshgrid(np.arange(-5, 6),
                          np.arange(-5, 6),
                          np.arange(-5, 6),
                          indexing='ij')
    image = np.array((x, y, z)).transpose(1, 2, 3, 0).astype(float)
    #image = x[:, :, :, None]
    affine = np.array([[1, 0, 0, -5],
                       [0, 1, 0, -5],
                       [0, 0, 1, -5],
                       [0, 0, 0, 1]], dtype=float)
    return image, affine

@pytest.fixture
def nonl_transfomation_id():
    x, y, z = np.meshgrid(np.arange(-6, 7),
                          np.arange(-6, 7),
                          np.arange(-6, 7),
                          indexing='ij')
    nonl_transform = np.concatenate(
        (x[..., None], y[..., None], z[..., None]), axis=3).astype(float)
    affine = np.array([[1, 0, 0, -6],
                       [0, 1, 0, -6],
                       [0, 0, 1, -6],
                       [0, 0, 0, 1]], dtype=float)
    return nonl_transform, affine


@pytest.fixture
def nonl_transfomation_reflex():
    x, y, z = np.meshgrid(np.arange(-6, 7)[::-1],
                          np.arange(-6, 7),
                          np.arange(-6, 7),
                          indexing='ij')
    nonl_transform = np.concatenate(
        (x[..., None], y[..., None], z[..., None]), axis=3).astype(float)
    affine = np.array([[1, 0, 0, -6],
                       [0, 1, 0, -6],
                       [0, 0, 1, -6],
                       [0, 0, 0, 1]], dtype=float)
    return nonl_transform, affine


@pytest.fixture
def affine_transformation():
    affine = np.array([[1, 0, 0, -1],
                       [0, 1, 0, -1],
                       [0, 0, 1, -1],
                       [0, 0, 0, 1]], dtype=float)

    target_affine = np.array(
        [[1, 0, 0, -6],
         [0, 1, 0, -6],
         [0, 0, 1, -6],
         [0, 0, 0, 1]], dtype=float)
    dim = [13, 13, 13]
    return affine, target_affine, dim

@pytest.fixture
def nonl_transfomation_rotate():
    x, y, z = np.meshgrid(np.arange(-6, 7, 2),
                          np.arange(-6, 7, 2),
                          np.arange(-6, 7, 2),
                          indexing='ij')
    nonl_transform = np.concatenate(
        (-y[..., None], x[..., None], z[..., None]), axis=3).astype(float)
    affine = np.array([[2, 0, 0, -6],
                       [0, 2, 0, -6],
                       [0, 0, 2, -6],
                       [0, 0, 0, 1]], dtype=float)
    return nonl_transform, affine


@pytest.fixture
def nonl_transfomation_shear():
    x, y, z = np.meshgrid(np.arange(-6, 7),
                          np.arange(-12, 13, 2),
                          np.arange(-6, 7),
                          indexing='ij')
    nonl_transform = np.concatenate(
        (-y[..., None], x[..., None], z[..., None]), axis=3).astype(float)
    affine = np.array([[1, 0, 0, -6],
                       [0, 1, 0, -6],
                       [0, 0, 1, -6],
                       [0, 0, 0, 1]], dtype=float)
    return nonl_transform, affine



@pytest.fixture
def nonl_large_rotate():
    x, y, z = np.meshgrid(np.linspace(-120, 120, 10),
                          np.linspace(-120, 120, 10),
                          np.linspace(-120, 120, 10),
                          indexing='ij')
    nonl_transform = np.concatenate(
        (-y[..., None], x[..., None], z[..., None]), axis=3).astype(float)
    s = 240 / 9.
    affine = np.array([[s, 0, 0, -120],
                       [0, s, 0, -120],
                       [0, 0, s, -120],
                       [0, 0, 0, 1]], dtype=float)
    return nonl_transform, affine



@pytest.fixture
def affine_rotate():
    affine = np.array([[0, -1, 0, 1],
                       [1, 0, 0, -1],
                       [0, 0, 1, -1],
                       [0, 0, 0, 1]], dtype=float)

    target_affine = np.array(
        [[1, 0, 0, -6],
         [0, 1, 0, -6],
         [0, 0, 1, -6],
         [0, 0, 0, 1]], dtype=float)
    dim = [13, 13, 13]
    return affine, target_affine, dim

@pytest.fixture
def affine_reflex():
    affine = np.array(
        [[-1, 0, 0, 1],
         [0, 1, 0, -1],
         [0, 0, 1, -1],
         [0, 0, 0, 1]], dtype=float)

    target_affine = np.array(
        [[1, 0, 0, -6],
         [0, 1, 0, -6],
         [0, 0, 1, -6],
         [0, 0, 0, 1]], dtype=float)
    dim = [13, 13, 13]
    return affine, target_affine, dim

@pytest.fixture
def shear():
    affine = np.array([[0, -2, 0, 0],
                       [1, 0, 0,  0],
                       [0, 0, 1,  0],
                       [0, 0, 0, 1]], dtype=float)

    target_affine = np.array(
        [[1, 0, 0, -6],
         [0, 1, 0, -6],
         [0, 0, 1, -6],
         [0, 0, 0, 1]], dtype=float)
    dim = [13, 13, 13]
    return affine, target_affine, dim



class TestNonlTransformVolume:
    def test_nonl_identity(self, image, nonl_transfomation_id):
        transformed = transformations.volumetric_nonlinear(image, nonl_transfomation_id)
        assert np.allclose(image[0], transformed[1:-1, 1:-1, 1:-1])

    def test_nonl_ressample(self, image, nonl_transfomation_id):
        t_affine = nonl_transfomation_id[1]
        dimensions = (13, 13, 13)
        transformed = transformations.volumetric_nonlinear(image,
                                                           nonl_transfomation_id,
                                                           target_space_affine=t_affine,
                                                           target_dimensions=dimensions)
        assert np.allclose(image[0], transformed[1:-1, 1:-1, 1:-1])

    def test_nonl_ressample2(self, image, nonl_transfomation_id):
        t_affine = nonl_transfomation_id[1].copy()
        t_affine[0, 0] = 2
        dimensions = [7, 13, 13]
        transformed = transformations.volumetric_nonlinear(image,
                                                           nonl_transfomation_id,
                                                           target_space_affine=t_affine,
                                                           target_dimensions=dimensions)
        assert np.allclose(image[0][1::2, :, :], transformed[1:-1, 1:-1, 1:-1])

    def test_nonl_reflex(self, image, nonl_transfomation_reflex):
        transformed = transformations.volumetric_nonlinear(
            image, nonl_transfomation_reflex)
        assert np.allclose(image[0][..., 0], transformed[-2:0:-1, 1:-1, 1:-1, 0])
        assert np.allclose(image[0][..., 1], transformed[-2:0:-1, 1:-1, 1:-1, 1])
        assert np.allclose(image[0][..., 2], transformed[-2:0:-1, 1:-1, 1:-1, 2])

        transformed = transformations.volumetric_nonlinear(
            image, nonl_transfomation_reflex,
            inverse_deformation=nonl_transfomation_reflex)
        assert np.allclose(image[0][..., 0],-transformed[-2:0:-1, 1:-1, 1:-1, 0])
        assert np.allclose(image[0][..., 1], transformed[-2:0:-1, 1:-1, 1:-1, 1])
        assert np.allclose(image[0][..., 2], transformed[-2:0:-1, 1:-1, 1:-1, 2])

    def test_nonl_rotate(self, image, nonl_transfomation_rotate):
        transformed = transformations.volumetric_nonlinear(
            image, nonl_transfomation_rotate)
        assert np.allclose(image[0][..., 0][1::2, 1::2, 1::2],
                           transformed[..., 0].transpose(1, 0, 2)[-2:0:-1, 1:-1, 1:-1])
        assert np.allclose(image[0][..., 1][1::2, 1::2, 1::2],
                           transformed[..., 1].transpose(1, 0, 2)[-2:0:-1, 1:-1, 1:-1])
        assert np.allclose(image[0][..., 2][1::2, 1::2, 1::2],
                           transformed[..., 2].transpose(1, 0, 2)[-2:0:-1, 1:-1, 1:-1])

        inv = nonl_transfomation_rotate[0].copy()
        inv[..., 1] *= -1
        inv[..., 0] *= -1
        inv = (inv, nonl_transfomation_rotate[1])
        transformed = transformations.volumetric_nonlinear(
            image, nonl_transfomation_rotate,
            inverse_deformation=inv)
        assert np.allclose(-image[0][..., 0][1::2, 1::2, 1::2],
                           transformed[..., 1].transpose(1, 0, 2)[-2:0:-1, 1:-1, 1:-1])
        assert np.allclose(image[0][..., 1][1::2, 1::2, 1::2],
                           transformed[..., 0].transpose(1, 0, 2)[-2:0:-1, 1:-1, 1:-1])
        assert np.allclose(image[0][..., 2][1::2, 1::2, 1::2],
                           transformed[..., 2].transpose(1, 0, 2)[-2:0:-1, 1:-1, 1:-1])

    def test_nonl_shear(self, image, nonl_transfomation_shear):
        transformed = transformations.volumetric_nonlinear(
            image, nonl_transfomation_shear)
        assert np.allclose(image[0][..., 0][1::2, ...],
                           transformed[..., 0].transpose(1, 0, 2)[-5:3:-1, 1:-1, 1:-1])
        assert np.allclose(image[0][..., 1][1::2, ...],
                           transformed[..., 1].transpose(1, 0, 2)[-5:3:-1, 1:-1, 1:-1])
        assert np.allclose(image[0][..., 2][1::2, ...],
                           transformed[..., 2].transpose(1, 0, 2)[-5:3:-1, 1:-1, 1:-1])

        inv = nonl_transfomation_shear[0].copy()
        inv[..., 1] *= -.5
        inv[..., 0] *= -.5
        inv = (inv, nonl_transfomation_shear[1])
        transformed = transformations.volumetric_nonlinear(
            image, nonl_transfomation_shear,
            inverse_deformation=inv, keep_vector_length=False)
        assert np.allclose(-.5 * image[0][..., 0][1::2, ...],
                           transformed[..., 1].transpose(1, 0, 2)[-5:3:-1, 1:-1, 1:-1])
        assert np.allclose(image[0][..., 1][1::2, ...],
                           transformed[..., 0].transpose(1, 0, 2)[-5:3:-1, 1:-1, 1:-1])
        assert np.allclose(image[0][..., 2][1::2, ...],
                           transformed[..., 2].transpose(1, 0, 2)[-5:3:-1, 1:-1, 1:-1])

        transformed = transformations.volumetric_nonlinear(
            image, nonl_transfomation_shear,
            inverse_deformation=inv, keep_vector_length=True)
        assert np.allclose(np.linalg.norm(image[0][1::2, ...], axis=3),
                           np.linalg.norm(
                               transformed[1:-1, -5:3:-1, 1:-1], axis=3).transpose(1, 0, 2))


class TestAffineTransfom:
    def test_affine(self, image, affine_transformation):
        transformed = transformations.volumetric_affine(image,
                                                        affine_transformation[0],
                                                        affine_transformation[1],
                                                        affine_transformation[2])
        assert np.allclose(image[0], transformed[2:, 2:, 2:])

    def test_affine_reflex(self, image, affine_reflex):
        transformed = transformations.volumetric_affine(image,
                                                        affine_reflex[0],
                                                        affine_reflex[1],
                                                        affine_reflex[2])
        assert np.allclose(image[0][..., 0],-transformed[-1:1:-1, 2:, 2:, 0])
        assert np.allclose(image[0][..., 1], transformed[-1:1:-1, 2:, 2:, 1])
        assert np.allclose(image[0][..., 2], transformed[-1:1:-1, 2:, 2:, 2])

    def test_affine_shear(self, image, shear):
        im = (image[0][...,  0], image[1])
        transformed = transformations.volumetric_affine(im,
                                                        shear[0],
                                                        shear[1],
                                                        shear[2])
        assert np.allclose(image[0][-2:0:-2, :, :, 0].transpose(1, 0, 2), transformed[1:-1, 4:-4, 1:-1])

        transformed = transformations.volumetric_affine(image,
                                                        shear[0],
                                                        shear[1],
                                                        shear[2],
                                                        keep_vector_length=False)
        assert np.allclose(-.5 * image[0][-2:0:-2, :, :, 0].transpose(1, 0, 2),
                           transformed[1:-1, 4:-4, 1:-1, 1])
        assert np.allclose(image[0][-2:0:-2, :, :, 1].transpose(1, 0, 2),
                           transformed[1:-1, 4:-4, 1:-1, 0])
        assert np.allclose(image[0][-2:0:-2, :, :, 2].transpose(1, 0, 2),
                           transformed[1:-1, 4:-4, 1:-1, 2])

        transformed = transformations.volumetric_affine(image,
                                                        shear[0],
                                                        shear[1],
                                                        shear[2],
                                                        keep_vector_length=True)
        assert np.allclose(np.linalg.norm(image[0][-2:0:-2, :, :, :], axis=3).transpose(1, 0, 2),
                           np.linalg.norm(transformed[1:-1, 4:-4, 1:-1], axis=3))


class TestWarpCoordinates:
    def test_nonl_ident(self, nonl_transfomation_id):
        coords = np.random.rand(10, 3) * 12 - 6
        coords_transf = transformations.coordinates_nonlinear(
            coords, nonl_transfomation_id)
        assert np.allclose(coords, coords_transf)

    def test_nonl_reflex(self, nonl_transfomation_reflex):
        coords = np.random.rand(10, 3) * 12 - 6
        coords_transf = transformations.coordinates_nonlinear(
            coords, nonl_transfomation_reflex)
        coords[:, 0] *= -1
        assert np.allclose(coords, coords_transf)

    def test_nonl_vec(self, nonl_transfomation_rotate):
        coords = np.random.rand(4, 3) * 8 - 4
        vec = np.random.rand(4, 3)
        coords_transf, vec_transf = transformations.coordinates_nonlinear(
            coords, nonl_transfomation_rotate, vectors=vec, keep_length=True)
        assert np.allclose(coords_transf[:, 1], coords[:, 0])
        assert np.allclose(coords_transf[:, 0],-coords[:, 1])
        assert np.allclose(coords_transf[:, 2], coords[:, 2])

        assert np.allclose(vec_transf[:, 1], vec[:, 0], rtol=1e-2, atol=1e-2)
        assert np.allclose(vec_transf[:, 0],-vec[:, 1], rtol=1e-2, atol=1e-2)
        assert np.allclose(vec_transf[:, 2], vec[:, 2], rtol=1e-2, atol=1e-2)


    def test_affine(self, affine_transformation):
        coords = np.random.rand(10, 3) * 12 - 6
        coords_transf = transformations.coordinates_affine(
            coords, affine_transformation[0])
        assert np.allclose(coords - 1, coords_transf)

    def test_affine_vector(self, affine_rotate):
        coords = np.random.rand(10, 3) * 12 - 6

        coords_transf = transformations.coordinates_affine(
            coords, affine_rotate[0])

        assert np.allclose(coords_transf[:, 1], coords[:, 0] - 1)
        assert np.allclose(coords_transf[:, 0],-coords[:, 1] + 1)
        assert np.allclose(coords_transf[:, 2], coords[:, 2] - 1)

        vec = np.random.rand(10, 3)
        vec_transf = transformations.vectors_affine(
            vec, affine_rotate[0], keep_length=False)
        assert np.allclose(vec_transf[:, 1], vec[:, 0], rtol=1e-2, atol=1e-2)
        assert np.allclose(vec_transf[:, 0],-vec[:, 1], rtol=1e-2, atol=1e-2)
        assert np.allclose(vec_transf[:, 2], vec[:, 2], rtol=1e-2, atol=1e-2)

        R = np.eye(3, 3) * 2
        vec_transf = transformations.vectors_affine(
            vec, R, keep_length=True)
        assert np.allclose(vec_transf, vec)

class TestProjectOnScalp():
    def test_project_scalp_0_dist(self, sphere3_msh):
        coords = np.array([[100., 0., 0.], [0., -100., 0]])
        c = transformations.project_on_scalp(coords, sphere3_msh)
        assert np.allclose(c, np.array([[95., 0., 0.], [0., -95., 0.]]))

    def test_project_scalp_dist(self, sphere3_msh):
        coords = np.array([[100., 0., 0.], [0., -100., 0]])
        c = transformations.project_on_scalp(coords, sphere3_msh, distance=2)
        assert np.allclose(c, np.array([[97., 0., 0.], [0., -97., 0.]]),
                           atol=1e-1)

class TestTransformPositions:
    def test_transform_tdcs_position_linear(self, sphere3_msh, affine_rotate):
        coords = np.array([[100., 0., 0.], [0., -100., 0]])
        pos_y = np.array([[-100., 0., 0.], [0., 100., 0]])
        c, y = transformations.transform_tdcs_positions(
            coords, 'affine', affine_rotate[0], pos_y, sphere3_msh)
        assert np.allclose(c, [[0., 95, 0], [95, 0, 0]])
        assert np.allclose(y, [[0.,-95, 0], [-95, 0, 0]])

    def test_transform_tdcs_position_nonl(self, sphere3_msh, nonl_large_rotate):
        coords = np.array([[100., 0., 0.], [0., -100., 0]])
        pos_y = np.array([[-100., 0., 0.], [0., 100., 0]])
        c, y = transformations.transform_tdcs_positions(
            coords, 'nonl', nonl_large_rotate, pos_y, sphere3_msh)
        assert np.allclose(c, [[0., 95, 0], [95, 0, 0]])
        assert np.allclose(y, [[0.,-95, 0], [-95, 0, 0]])

    def test_transform_tms_position_linear(self, sphere3_msh, affine_rotate):
        coords = np.array([[100., 0., 0.], [0., -100., 0]])
        v_y = np.array([[-1. / np.sqrt(2), 1. / np.sqrt(2), 1.],
                        [ 1. / np.sqrt(2),-1. / np.sqrt(2), 1.]])
        v_z = np.array([[0, 0, 1.],
                        [0, 0, 2.]])
        distances = np.array([1, 2])
        c, y, z = transformations.transform_tms_positions(
            coords, v_y, v_z, 'affine', affine_rotate[0], sphere3_msh, distances)
        assert np.allclose(c, [[0., 96, 0], [97, 0, 0]], atol=1e-2)
        assert np.allclose(y, [[-1. / np.sqrt(2),-1. / np.sqrt(2), 0],
                               [1. / np.sqrt(2),1. / np.sqrt(2), 0]])
        assert np.allclose(z, [[0, 0, 1],
                               [0, 0, 1]])


    def test_transform_tms_position_nonl(self, sphere3_msh, nonl_large_rotate):
        coords = np.array([[100., 0., 0.], [0., -100., 0]])
        v_y = np.array([[-1. / np.sqrt(2), 1. / np.sqrt(2), 1.],
                        [ 1. / np.sqrt(2),-1. / np.sqrt(2), 1.]])
        v_z = np.array([[0, 0, 1.],
                        [0, 0, 2.]])
        distances = np.array([1, 2])
        c, y, z = transformations.transform_tms_positions(
            coords, v_y, v_z, 'nonl', nonl_large_rotate, sphere3_msh, distances)
        assert np.allclose(c, [[0., 96, 0], [97, 0, 0]], atol=1e-2)
        assert np.allclose(y, [[-1. / np.sqrt(2),-1. / np.sqrt(2), 0],
                               [1. / np.sqrt(2),1. / np.sqrt(2), 0]], atol=1e-2)
        assert np.allclose(z, [[0, 0, 1],
                               [0, 0, 1]])



class TestCropVol:
    '''
    def test_crop_vol_data(self):
        data_in = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_resize_AT',
            'before_reduceBrain.mat'),
            squeeze_me=False
        )
        data_out = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_resize_AT',
            'after_reduceBrain.mat'),
            squeeze_me=False
        )
        cropped, new_affine, BB = transformations.crop_vol(
            data_in['Ymfs'], data_in['V']['mat'][0, 0], data_in['mask'], 4
        )
        assert np.allclose(data_out['Ymfs'], cropped)
        assert np.allclose(data_out['BB']['BB'][0, 0], BB.reshape(-1, order='C') + 1)
    '''

    def test_crop_vol_simple(self):
        data = np.zeros((50, 100, 150))
        data[:] = np.arange(50)[:, None, None]

        mask = np.zeros(data.shape, dtype=bool)
        mask[10:20, 30:40, 50:60] = True

        affine = np.eye(4)
        affine[:3, :3] *= np.arange(1, 4)[::-1]

        cropped, new_affine, bb = transformations.crop_vol(data, affine, mask)
        assert cropped.shape == (10, 10, 10)
        assert np.allclose(cropped, np.arange(10, 20)[:, None, None])
        assert np.allclose(new_affine[:3, :3], affine[:3, :3])
        assert np.allclose(new_affine[:3, 3], [30, 60, 50])
        assert np.allclose(bb, [[10, 19], [30, 39], [50, 59]])

    def test_crop_vol_simple_boundary(self):
        data = np.zeros((50, 100, 200))
        data[:] = np.arange(50)[:, None, None]

        mask = np.zeros(data.shape, dtype=bool)
        mask[10:20, 30:40, 50:60] = True

        affine = np.eye(4)
        affine[:3, :3] *= [4, 2, 1]

        cropped, new_affine, bb = transformations.crop_vol(data, affine, mask, 4)

        assert cropped.shape == (12, 14, 18)
        assert np.allclose(cropped, np.arange(9, 21)[:, None, None])
        assert np.allclose(new_affine[:3, :3], affine[:3, :3])
        assert np.allclose(new_affine[:3, 3], [(10-1)*4, (30-2)*2, 50-4])
        assert np.allclose(bb, [[9, 20], [28, 41], [46, 63]])

class TestPasteVol:
    '''
    def test_paste_vol_data(self):
        data = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_resize_AT',
            'dereduce_brain.mat'),
            squeeze_me=False
        )
        pasted = transformations.paste_vol(
            data['Ymfs'],
            data['BB']['sizeT'][0, 0][0],
            data['BB']['BB'][0, 0][0].reshape(3, 2) - 1
        )

        assert np.allclose(pasted, data['dereduced'])
    '''
    def test_paste_simple(self):
        data = np.random.random((10, 20, 30))
        target_size = (100, 100, 100)
        bb = np.array([[10, 19], [15, 34], [20, 49]])

        pasted = transformations.paste_vol(data, target_size, bb)

        assert pasted.shape == target_size
        assert np.allclose(pasted[10:20, 15:35, 20:50], data)

class TestGetVoxSize:
    def test_get_vox_size(self):
        voxsize = transformations.get_vox_size(
            np.array([[0, 0, 4, 0],
                      [0, 2, 0, 1],
                      [1, 0, 0, 5],
                      [0, 0, 0, 1]])
        )
        assert np.allclose(voxsize, (1, 2, 4))


class TestResample:
    '''
    def test_interp_data(self):
        data_in = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_resize_AT',
            'before_interp.mat'),
            squeeze_me=False
        )
        data_out = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_resize_AT',
            'after_interp.mat'),
            squeeze_me=False
        )
        vol = data_in['Ymfs']
        vol[vol < 1] = 1
        affine = data_in['V']['mat'][0, 0]
        affine = np.eye(4)
        resampled, new_affine = transformations.resample_vol(
            vol, affine, data_in['interpV'],
            order=1
        )
        out_mat = data_out['resI']['hdrN'][0, 0]['mat'][0, 0]
        #assert np.allclose(resampled, data_out['Ymfs'])
    '''

    def test_resample_iso_0(self):
        vol = np.random.randint(0, 100, (50, 50, 50))
        vol = vol.astype(np.int32)
        affine = np.eye(4)
        target_res = 0.5

        resampled, new_affine, orig_res = transformations.resample_vol(
            vol, affine, target_res, order=0
        )
        # OPEN SEPARATELY
        # Overlaying does not work, I don't know why
        #vol_nii = nib.Nifti1Pair(vol, affine)
        #nib.save(vol_nii, 'in.nii.gz')
        #resampled_nii = nib.Nifti1Image(resampled, new_affine)
        #nib.save(resampled_nii, 'out.nii.gz')
        assert np.allclose(orig_res, 1.)
        assert np.all(resampled.shape == (100, 100, 100))
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    assert np.all(
                        resampled[i::2, j::2, k::2] ==
                        vol
                    )



    def test_resample_iso(self):
        vol = np.zeros((50, 50, 50))
        vol[:] = np.arange(50)[:, None, None]
        affine = np.eye(4)
        target_res = 0.5

        resampled, new_affine, orig_res = transformations.resample_vol(
            vol, affine, target_res, order=1
        )
        #nib.save(nib.Nifti1Image(vol, affine), 'in.nii.gz')
        #nib.save(nib.Nifti1Image(resampled, new_affine), 'out.nii.gz')
        assert np.allclose(orig_res, 1.)
        assert np.allclose(resampled.shape, (100, 100, 100))
        assert np.allclose(
            resampled[1:-1],
            np.arange(0.25, 49, step=0.5)[:, None, None]
        )
        assert np.allclose(new_affine[:3, :3], 0.5*np.eye(3))
        assert np.allclose(new_affine[:3, 3], -0.25)
        assert np.allclose(new_affine[3, :], [0, 0, 0, 1])


    def test_resample_aniso(self):
        vol = np.zeros((40, 20, 80))
        vol[:] = np.arange(80, step=4)[None, :, None]
        affine = np.eye(4)
        affine[:3, :3] = np.array(
            [[0, 4, 0],
             [2, 0, 0],
             [0, 0, 1]])
        target_res = np.array([1, 1, 1])

        resampled, new_affine, orig_res = transformations.resample_vol(
            vol, affine, target_res
        )

        assert np.allclose(orig_res, [2, 4, 1])
        assert np.allclose(resampled.shape, (80, 80, 80))
        assert np.allclose(
            resampled[:, 2:-2, :],
            np.arange(0.5, 76)[None, :, None]
        )
        assert np.allclose(new_affine[:3, :3], [[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        assert np.allclose(new_affine[:3, 3], [-1.5, -0.5, 0])

    def test_decrease_res(self):
        vol = np.zeros((50, 50, 50))
        vol[:] = np.arange(50)[:, None, None]
        affine = np.eye(4)
        target_res = 2

        resampled, new_affine, orig_res = transformations.resample_vol(
            vol, affine, target_res, order=1
        )
        #nib.save(nib.Nifti1Image(vol, affine), 'in.nii.gz')
        #nib.save(nib.Nifti1Image(resampled, new_affine), 'out.nii.gz')
        assert np.allclose(orig_res, [1, 1, 1])
        assert np.allclose(resampled.shape, (25, 25, 25))
        assert np.allclose(
            resampled,
            np.arange(0.5, 49, step=2)[:, None, None]
        )
        assert np.allclose(new_affine[:3, :3], 2*np.eye(3))
        assert np.allclose(new_affine[:3, 3], 0.5)
        assert np.allclose(new_affine[3, :], [0, 0, 0, 1])

def test_surf2surf(sphere3_msh):
    in_surf = sphere3_msh.crop_mesh(1005)
    in_nodes = in_surf.nodes.node_coord
    in_nodes /= np.average(np.linalg.norm(in_nodes, axis=1))
    field = in_nodes[:, 0]
    out_surf = sphere3_msh.crop_mesh(1004)
    out_nodes = out_surf.nodes.node_coord
    out_nodes /= np.average(np.linalg.norm(out_nodes, axis=1))
    out_field, _ = transformations._surf2surf(field, in_surf, out_surf)
    assert np.allclose(out_field, out_nodes[:, 0], atol=1e-1)

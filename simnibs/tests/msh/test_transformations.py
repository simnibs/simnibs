import os
import tempfile

import numpy as np
import pytest

import simnibs.msh.mesh_io as mesh_io
import simnibs.msh.transformations as transformations


@pytest.fixture
def sphere3_msh():
    fn = os.path.join(os.path.dirname(os.path.realpath(
        __file__)), '..', 'testing_files', 'sphere3.msh')
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





class TestCSV:
    def test_read_csv_generic(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'type,x,y,z,extra1,extra2\n')
        f.write(b'Generic,1,2,3,a1,a2\n')
        f.write(b'Generic,1.2,2.4,3.6,b1,b2\n')
        f.close()
        type_, coordinates, extra, name, extra_cols, header = transformations._read_csv(f.name)
        assert type_ == ['Generic', 'Generic']
        assert np.allclose(coordinates, [[1, 2, 3], [1.2, 2.4, 3.6]])
        assert header == ['type','x', 'y', 'z', 'extra1', 'extra2']
        assert name == ['a1', 'b1']
        assert extra_cols == [['a2'], ['b2']]

        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'type,x,y,z,extra1,extra2\n')
        f.write(b'Generic,1,2,3\n')
        f.write(b'Generic,1.2,2.4,3.6\n')
        f.close()
        type_, coordinates, extra, name, extra_cols, header = transformations._read_csv(f.name)
        assert type_ == ['Generic', 'Generic']
        assert np.allclose(coordinates, [[1, 2, 3], [1.2, 2.4, 3.6]])
        assert header == ['type','x', 'y', 'z', 'extra1', 'extra2']
        assert name == [None, None]
        assert extra_cols == [None, None]


    def test_read_csv_fiducial_electrode(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'Fiducial,1,2,3,fiducial,a2\n')
        f.write(b'Electrode,1.2,2.4,3.6,9,8,7,electrode,b2,b3\n')
        f.write(b'ReferenceElectrode,1.2,2.4,7.1\n')
        f.close()
        type_, coordinates, extra, name, extra_cols, header = transformations._read_csv(f.name)
        assert type_ == ['Fiducial', 'Electrode', 'ReferenceElectrode']
        assert np.allclose(coordinates, [[1, 2, 3],
                                         [1.2, 2.4, 3.6],
                                         [1.2, 2.4, 7.1]])
        assert header == []
        assert name == ['fiducial', 'electrode', None]
        assert extra_cols == [['a2'], ['b2', 'b3'], None]
        assert extra[0] is None
        assert np.allclose(extra[1], [9, 8, 7])
        assert extra[2] is None

    def test_read_csv_coil(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'Fiducial,1,2,3,fiducial,a2\n')
        f.write(b'ReferenceElectrode,1.2,2.4,7.1\n')
        f.write(b'CoilPos,1,2,3,4,5,6,7,8,9,10,coil,comment\n')
        f.close()
        type_, coordinates, extra, name, extra_cols, header = transformations._read_csv(f.name)
        assert type_ == ['Fiducial', 'ReferenceElectrode', 'CoilPos']
        assert np.allclose(coordinates, [[1, 2, 3],
                                         [1.2, 2.4, 7.1],
                                         [1, 2, 3]])
        assert header == []
        assert name == ['fiducial',  None, 'coil']
        assert extra_cols == [['a2'], None, ['comment']]
        assert extra[0] is None
        assert extra[1] is None
        assert np.allclose(extra[2], [4, 5, 6, 7, 8, 9, 10])

    def test_write_csv_generic(self):
        type_ = ['Generic', 'Generic']
        coordinates = np.array([[1, 2, 3], [1.2, 2.4, 3.6]])
        header = ['type', 'x', 'y', 'z', 'extra1', 'extra2']
        name = ['a1', 'b1']
        extra_cols = [['a2'], ['b2']]
        extra = [None, None]
        f = tempfile.NamedTemporaryFile(delete=False)
        transformations._write_csv(f.name, type_, coordinates, extra, name, extra_cols, header)
        with open(f.name, 'r') as f:
            assert f.readline().strip() == 'type,x,y,z,extra1,extra2'
            assert f.readline().strip() == 'Generic,1.0,2.0,3.0,a1,a2'
            assert f.readline().strip() == 'Generic,1.2,2.4,3.6,b1,b2'

    def test_write_csv_electrode(self):
        type_ = ['Fiducial', 'Electrode', 'ReferenceElectrode']
        coordinates = np.array([[1, 2, 3],
                                [1.2, 2.4, 3.6],
                                [1.2, 2.4, 7.1]])
        header = []
        name = ['fiducial', 'electrode', None]
        extra_cols = [['a2'], ['b2', 'b3'], None]
        extra = [None, np.array([9, 8, 7], dtype=float), None]

        f = tempfile.NamedTemporaryFile(delete=False)
        transformations._write_csv(f.name, type_, coordinates, extra, name, extra_cols, header)
        with open(f.name, 'r') as f:
            assert f.readline().strip() == 'Fiducial,1.0,2.0,3.0,fiducial,a2'
            assert f.readline().strip() == 'Electrode,1.2,2.4,3.6,9.0,8.0,7.0,electrode,b2,b3'
            assert f.readline().strip() == 'ReferenceElectrode,1.2,2.4,7.1'


    def test_write_csv_coil(self):
        type_ = ['Fiducial', 'ReferenceElectrode', 'CoilPos']
        coordinates = np.array([[1, 2, 3],
                                [1.2, 2.4, 7.1],
                                [1, 2, 3]])
        header = []
        name = ['fiducial',  None, 'coil']
        extra_cols = [['a2'], None, ['comment']]
        extra = [None, None, np.array([4., 5., 6., 7., 8., 9., 10], dtype=float)]

        f = tempfile.NamedTemporaryFile(delete=False)
        transformations._write_csv(f.name, type_, coordinates, extra, name, extra_cols, header)
        with open(f.name, 'r') as f:
            assert f.readline().strip() == 'Fiducial,1.0,2.0,3.0,fiducial,a2'
            assert f.readline().strip() == 'ReferenceElectrode,1.2,2.4,7.1'
            assert f.readline().strip() == \
                'CoilPos,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,coil,comment'

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

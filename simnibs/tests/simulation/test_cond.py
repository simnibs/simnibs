import copy
import os
import numpy as np
import pytest
import nibabel

import simnibs.simulation.cond as cond
import simnibs.msh.mesh_io as mesh_io

@pytest.fixture
def sphere3_msh():
    fn = os.path.join(os.path.dirname(os.path.realpath(
        __file__)), '..', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn)

@pytest.fixture
def tensor():
    v1 = np.random.rand(3)
    v1 /= np.linalg.norm(v1)
    v2 = np.random.rand(3)
    v2 = v2 - v1.dot(v2) * v1 / np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    v3 = np.random.rand(3)
    v3 = v3 - v1.dot(v3) * v1 / np.linalg.norm(v1)
    v3 = v3 - v2.dot(v3) * v2 / np.linalg.norm(v2)
    v3 /= np.linalg.norm(v3)
    t = np.outer(v1, v1) + 2 * np.outer(v2, v2) + 3 * np.outer(v3, v3)
    return t


class TestCond2Elmdata:
    def test_isotropic(self, sphere3_msh):
        cond_list = [None, None, 1, 2, 3]
        c = cond.cond2elmdata(sphere3_msh, cond_list).value
        assert np.all(c[sphere3_msh.elm.tag1==3] == 1)
        assert np.all(c[sphere3_msh.elm.tag1==1003] == 1)
        assert np.all(c[sphere3_msh.elm.tag1==4] == 2)
        assert np.all(c[sphere3_msh.elm.tag1==1004] == 2)
        assert np.all(c[sphere3_msh.elm.tag1==5] == 3)
        assert np.all(c[sphere3_msh.elm.tag1==1005] == 3)

    def test_anisotropic(self, sphere3_msh, tensor):
        cond_list = [None, None, 1, 2, 3]
        v = np.zeros((255, 255, 255, 6))
        for i, j in enumerate([0, 1, 2, 4, 5, 8]):
            v[:, :, :, i] = tensor.reshape(-1)[j]
        affine = np.array([[-1, 0, 0, 128],
                           [0, 1, 0, -128],
                           [0, 0, 1, -127],
                           [0, 0, 0, 1]])

        msh = sphere3_msh
        elmcond = cond.cond2elmdata(msh, cond_list, v, affine,
                                    aniso_tissues=3,
                                    max_cond=np.inf,
                                    max_ratio=np.inf,
                                    correct_intensity=False)
        t = affine[:3, :3].dot(tensor).dot(affine[:3, :3])
        assert np.allclose(elmcond.value[msh.elm.tag1 == 3],
                           t.reshape(-1))
        assert np.allclose(elmcond.value[msh.elm.tag1 == 4],
                           [2, 0, 0, 0, 2, 0, 0, 0, 2])
        assert np.allclose(elmcond.value[msh.elm.tag1 == 5],
                           [3, 0, 0, 0, 3, 0, 0, 0, 3])
        vol = np.linalg.eig(tensor)[0].prod() ** (1./3.)
        scaling = (vol * 1 + vol * 2) / (2 * (vol ** 2))
        elmcond = cond.cond2elmdata(msh, cond_list, v, affine,
                                    aniso_tissues=[3, 4],
                                    max_cond=np.inf,
                                    max_ratio=np.inf,
                                    correct_intensity=True)
        assert np.allclose(elmcond.value[msh.elm.tag1 == 3],
                           scaling * t.reshape(-1))
        assert np.allclose(elmcond.value[msh.elm.tag1 == 4],
                           scaling * t.reshape(-1))
        assert np.allclose(elmcond.value[msh.elm.tag1 == 5],
                           [3, 0, 0, 0, 3, 0, 0, 0, 3])

    def test_anisotropic_normalized(self, sphere3_msh, tensor):
        t = tensor
        cond_list = [None, None, 2, 7, 9]
        v = np.zeros((255, 255, 255, 6))
        for i, j in enumerate([0, 1, 2, 4, 5, 8]):
            v[..., i] = t.reshape(-1)[j]
        affine = np.array([[-1, 0, 0, 128],
                           [0, 1, 0, -128],
                           [0, 0, 1, -127],
                           [0, 0, 0, 1]])

        elmcond = cond.cond2elmdata(sphere3_msh, cond_list, v, affine, aniso_tissues=3,
                                    normalize=True)
        tensor = elmcond.value[sphere3_msh.elm.tag1 == 3].reshape(-1, 3, 3)
        t = affine[:3, :3].dot(t).dot(affine[:3, :3])
        assert np.allclose(np.linalg.det(tensor), 2 ** 3)
        assert np.allclose(np.linalg.eig(tensor)[1][:, 0], np.linalg.eig(t)[1][0])
        assert np.allclose(np.linalg.eig(tensor)[1][:, 1], np.linalg.eig(t)[1][1])
        assert np.allclose(np.linalg.eig(tensor)[1][:, 2], np.linalg.eig(t)[1][2])
        assert np.allclose(elmcond.value[sphere3_msh.elm.tag1 == 4],
                           [7, 0, 0, 0, 7, 0, 0, 0, 7])
        assert np.allclose(elmcond.value[sphere3_msh.elm.tag1 == 5],
                           [9, 0, 0, 0, 9, 0, 0, 0, 9])

    def test_anisotropic_zero_region(self, sphere3_msh, tensor):
        t = tensor
        cond_list = [None, None, 2, 7, 9]
        v = np.zeros((255, 255, 255, 6))
        affine = np.array([[-1, 0, 0, 128],
                           [0, 1, 0, -128],
                           [0, 0, 1, -127],
                           [0, 0, 0, 1]])

        elmcond = cond.cond2elmdata(sphere3_msh, cond_list, v, affine, aniso_tissues=3,
                                    normalize=True)
        assert np.allclose(elmcond.value[sphere3_msh.elm.tag1 == 3],
                           [2, 0, 0, 0, 2, 0, 0, 0, 2])
        '''
        elmcond = cond.cond2elmdata(sphere3_msh, cond_list, v, affine, aniso_tissues=3,
                                    normalize=False)
        assert np.allclose(elmcond.value[sphere3_msh.elm.tag1 == 3],
                           [2, 0, 0, 0, 2, 0, 0, 0, 2])
        '''

    def test_anisotropic_scale_excentricity(self, sphere3_msh, tensor):
        t = tensor
        cond_list = [None, None, 2, 7, 9]
        v = np.zeros((255, 255, 255, 6))
        for i, j in enumerate([0, 1, 2, 4, 5, 8]):
            v[..., i] = t.reshape(-1)[j]
        affine = np.array([[-1, 0, 0, 128],
                           [0, 1, 0, -128],
                           [0, 0, 1, -127],
                           [0, 0, 0, 1]])

        elmcond = cond.cond2elmdata(sphere3_msh, cond_list, v, affine, aniso_tissues=3,
                                    normalize=True, excentricity_scaling=0.)
        tensor = elmcond.value[sphere3_msh.elm.tag1 == 3].reshape(-1, 3, 3)
        t = affine[:3, :3].dot(t).dot(affine[:3, :3])
        assert np.allclose(np.linalg.det(tensor), 2 ** 3)
        assert np.allclose(tensor, 2 * np.eye(3))



class TestExcentricity:
    def test_excentricity_no_change(self, tensor):
        t = np.tile(np.array([10., 5, 1]), [10, 1])
        t[1] = [1, 1, 1]
        t2 = cond._adjust_excentricity(t, .5)
        assert np.allclose(t2, t)

    def test_excentricity_to_iso(self, tensor):
        t = np.tile(np.array([10., 5, 1]), [10, 1])
        t[1, :] = 50 ** (1./3.)
        t2 = cond._adjust_excentricity(t, 0)
        assert np.allclose(t2.prod(axis=1), t.prod(axis=1))
        assert np.allclose(t2, 50 ** (1./3.))

    def test_iso_to_excentric(self):
        t = np.tile(np.array([1., 1, 1]), [10, 1])
        t *= np.arange(1, 11)[:, None]
        t2 = cond._adjust_excentricity(t, .99)
        assert np.allclose(t2, t)


class TestFixTensors:
    def test_change_negative_eigv(self):
        eigv = np.array([[10, 1, -2], [20, 2, -1]])
        eigv = cond._fix_eigv(eigv, 100, 10, 1)
        assert np.allclose(eigv, [[10, 1, 1], [20, 2, 2]])

    def test_change_large_eigv(self):
        eigv = np.array([[100, 1, 1], [200, 2, 2]])
        eigv = cond._fix_eigv(eigv, 50, 100, 1)
        assert np.allclose(eigv, [[50, 1, 1], [50, 2, 2]])

    def test_change_small_eigv(self):
        eigv = np.array([[100, 1, 1], [200, 2, 2]])
        eigv = cond._fix_eigv(eigv, 500, 10, 1)
        assert np.allclose(eigv, [[100, 10, 10], [200, 20, 20]])

    def test_change_negative_semidef(self):
        eigv = np.array([[100, 1, 1], [-1, -2, -1]])
        eigv = cond._fix_eigv(eigv, 500, 10, 1)
        assert np.allclose(eigv, [[100, 10, 10], [1, 1, 1]])


class TestTensorVisualization:
    def test_visualization(self, tensor):
        eig_val, eig_vec = np.linalg.eig(tensor)
        max_eig = eig_val.argmax()
        min_eig = eig_val.argmin()
        t = np.tile(tensor, [10, 1, 1])
        factor = np.arange(1, 11)
        t *= factor[:, None, None]
        t = mesh_io.ElementData(t.reshape(-1, 9))
        fields = cond.TensorVisualization(t, None, all_compoents=True)
        assert np.allclose(fields[0].value / factor[:, None],
                           eig_val[max_eig] * eig_vec[:, max_eig])
        assert np.allclose(fields[2].value / factor[:, None],
                           eig_val[min_eig] * eig_vec[:, min_eig])
        assert np.allclose(fields[3].value / factor,
                           np.prod(eig_val) ** (1.0/3.0))
'''
        v1 = np.random.rand(3)
        v1 /= np.linalg.norm(v1)
        v2 = np.random.rand(3)
        v2 = v2 - v1.dot(v2) * v1 / np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)
        v3 = np.random.rand(3)
        v3 = v3 - v1.dot(v3) * v1 / np.linalg.norm(v1)
        v3 = v3 - v2.dot(v3) * v2 / np.linalg.norm(v2)
        v3 /= np.linalg.norm(v3)
        t = np.outer(v1, v1) + 2 * np.outer(v2, v2) + 3 * np.outer(v3, v3)
        return t
'''


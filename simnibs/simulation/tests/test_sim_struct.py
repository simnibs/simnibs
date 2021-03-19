import os
import copy
import shutil
import tempfile

import pytest
import h5py
from mock import Mock, patch, call
import scipy.io
import numpy as np
import nibabel

from ... import SIMNIBSDIR
from .. import sim_struct
from ...mesh_tools import mesh_io

@pytest.fixture(scope='module')
def sphere3_fn():
    return os.path.join(
        SIMNIBSDIR, '_internal_resources',
        'testing_files', 'sphere3.msh')

@pytest.fixture
def sphere3_msh():
    return mesh_io.read_msh(os.path.join(
        SIMNIBSDIR, '_internal_resources',
        'testing_files', 'sphere3.msh'))



@pytest.fixture(scope='module')
def mat_session():
    mat = scipy.io.loadmat(os.path.join(
        SIMNIBSDIR, '_internal_resources',
        'testing_files', 'session.mat'),
        struct_as_record=True, squeeze_me=False)

    return mat


@pytest.fixture(scope='module')
def sphere3_folders():
    """ Creates folder structure
    """
    path_to_dir = os.path.join(
        SIMNIBSDIR, '_internal_resources',
        'testing_files', 'd2c_sphere3')
    if not os.path.exists(path_to_dir):
        os.mkdir(path_to_dir)

    with open(os.path.join(path_to_dir, 'CTI_vn_tensor.nii.gz'), 'w') as f:
        f.write(b'')
    with open(os.path.join(path_to_dir, 'CTI_dir_tensor.nii.gz'), 'w') as f:
        f.write(b'')
    with open(os.path.join(path_to_dir, 'CTI_dir_mc.nii.gz'), 'w') as f:
        f.write(b'')

    def fin():
        shutil.rmtree(path_to_dir)

    return os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files')

class TestFiducials:
    def test_fiducial_from_csv(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'Fiducial,0,82.9,-43,Nz\n')
        f.write(b'Fiducial,0,-116.2,-30.5,Iz\n')
        f.write(b'Fiducial,-79.6,-20.6,-48.6,LPA\n')
        f.write(b'Fiducial,80.6,-20.6,-48.1,RPA\n')
        f.close()
        fiducials = sim_struct.FIDUCIALS()
        fiducials.from_csv(f.name)
        assert np.allclose(fiducials.Nz, [0, 82.9 ,-43])
        assert np.allclose(fiducials.Iz, [0,-116.2,-30.5])
        assert np.allclose(fiducials.LPA,[-79.6,-20.6,-48.6])
        assert np.allclose(fiducials.RPA,[80.6,-20.6,-48.1])

class TestList:
    def test_list_anisotropy(self):
        l = sim_struct.SimuList()
        with pytest.raises(ValueError):
            l.anisotropy_type = None

    def test_list_postprocess(self):
        l = sim_struct.SimuList()
        with pytest.raises(ValueError):
            l.postprocess = 'Y'

    def test_list_read_mat_cond(self, mat_session):
        l = sim_struct.SimuList()
        l.read_cond_mat_struct(mat_session['poslist'][0][0][0][0])
        assert l.cond[1].name == 'GM'

    def test_list_cond_mat_struct(self, mat_session):
        l = sim_struct.SimuList()
        l.cond[0].distribution_type = 'uniform'
        l.cond[0].distribution_parameters = [1.5, 3.0]
        mat = l.cond_mat_struct()
        scipy.io.savemat('tmp.mat', mat)

        mat = scipy.io.loadmat(
            'tmp.mat', struct_as_record=True, squeeze_me=False)
        l2 = sim_struct.SimuList()
        l2.read_cond_mat_struct(mat)
        os.remove('tmp.mat')
        assert l2.cond[0].distribution_type == 'uniform'
        assert np.all(l2.cond[0].distribution_parameters == [1.5, 3.0])

    def test_list_cond_array(self, sphere3_msh):
        l = sim_struct.SimuList()
        l.cond[0].value = None
        l.cond[1].value = None
        l.cond[2].value = 21
        l.cond[3].value = 31
        l.cond[4].value = 41
        l.mesh = sphere3_msh
        elmcond = l.cond2elmdata()
        assert np.all(elmcond.value[sphere3_msh.elm.tag1 == 3] == 21)
        assert np.all(elmcond.value[sphere3_msh.elm.tag1 == 4] == 31)
        assert np.all(elmcond.value[sphere3_msh.elm.tag1 == 5] == 41)

    def test_get_conductivity(self):
        l = sim_struct.SimuList()
        l.cond[0].value = 5
        assert l.conductivity[1].value == 5

    def test_set_conductivity(self):
        l = sim_struct.SimuList()
        l.conductivity[1].value = 5
        assert l.cond[0].value == 5

    def test_list_cond_array_aniso(self, sphere3_msh):
        l = sim_struct.SimuList()
        l.anisotropy_type = 'dir'
        l.cond[3].value = 7
        l.cond[4].value = 9
        with tempfile.NamedTemporaryFile(suffix='.nii.gz') as f:
            l.fn_tensor_nifti = f.name
        v = np.zeros((255, 255, 255, 6))
        for i in range(6):
            v[:, :, :, i] = i
        affine = np.array([[-1, 0, 0, 128],
                           [0, 1, 0, -128],
                           [0, 0, 1, -127],
                           [0, 0, 0, 1]])

        img = nibabel.Nifti1Image(v, affine)
        nibabel.save(img, l.fn_tensor_nifti)
        msh = copy.deepcopy(sphere3_msh)
        msh.elm.tag1[msh.elm.tag1 == 3] = 1

        l.mesh = msh
        elmcond = l.cond2elmdata()
        assert np.allclose(elmcond.value[sphere3_msh.elm.tag1 == 1],
                           [0, -1, -2, -1, 3, 4, -2, 4, 5])
        assert np.allclose(elmcond.value[sphere3_msh.elm.tag1 == 4],
                           [7, 0, 0, 0, 7, 0, 0, 0, 7])
        assert np.allclose(elmcond.value[sphere3_msh.elm.tag1 == 5],
                           [9, 0, 0, 0, 9, 0, 0, 0, 9])
        os.remove(l.fn_tensor_nifti)

    def test_list_cond_array_aniso_vn(self, sphere3_msh):
        l = sim_struct.SimuList()
        l.anisotropy_type = 'vn'
        l.conductivity[3].value = 2
        l.conductivity[4].value = 7
        l.conductivity[5].value = 9
        l.anisotropic_tissues = [3]
        l.fn_tensor_nifti = '/tmp/test_aniso.nii.gz'
        l.mesh = sphere3_msh
        v = np.zeros((255, 255, 255, 6))
        # set-up tensor
        v1 = np.random.rand(3)
        v1 /= np.linalg.norm(v1)
        v2 = np.random.rand(3)
        v2 = v2 - v1.dot(v2) * v1 / np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)
        v3 = np.random.rand(3)
        v3 = v3 - v1.dot(v3) * v1 / np.linalg.norm(v1)
        v3 = v3 - v2.dot(v3) * v2 / np.linalg.norm(v2)
        v3 /= np.linalg.norm(v3)
        t = np.outer(v1, v1) + 2*np.outer(v2, v2) + 3*np.outer(v3, v3)
        v[:, :, :] = t.reshape(-1)[[0, 1, 2, 4, 5, 8]]
        affine = np.array([[-1, 0, 0, 128],
                           [0, 1, 0, -128],
                           [0, 0, 1, -127],
                           [0, 0, 0, 1]])

        l.anisotropy_vol = v
        l.anisotropy_affine = affine

        elmcond = l.cond2elmdata()
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


    def test_write_to_hdf5(self):
        l = sim_struct.SimuList()
        l.anisotropy_type = 'vn'
        l.cond[3].value = 2
        l.cond[0].distribution_type = 'beta'
        l.anisotropic_tissues = [3]
        v = np.zeros((255, 255, 255, 6))
        for i in range(6):
            v[:, :, :, i] = i
        affine = np.array([[-1, 0, 0, 128],
                           [0, 1, 0, -128],
                           [0, 0, 1, -127],
                           [0, 0, 0, 1]])
        l.anisotropy_vol = v
        l.anisotropy_affine = affine
        with tempfile.NamedTemporaryFile(suffix='.hdf5') as f:
            fn_hdf5 = f.name
        l._write_conductivity_to_hdf5(fn_hdf5)

        with h5py.File(fn_hdf5, 'r') as f:
            assert np.isclose(f['cond/values'][3], 2)
            assert np.isnan(f['cond/values'][46])
            assert f['cond/names'][0] == b'WM'
            assert np.all(f['cond'].attrs['anisotropic_tissues'] == [3])
            assert f['cond'].attrs['anisotropy_type'] == 'vn'
            assert np.allclose(f['cond/anisotropy_affine'], affine)
            assert np.allclose(f['cond/anisotropy_vol'], v)
            assert f['cond/distribution_types'][0] == b'beta'

        os.remove(fn_hdf5)

    def test_read_from_hdf5(self):
        l = sim_struct.SimuList()
        l.anisotropy_type = 'vn'
        l.cond[3].value = 2
        l.anisotropic_tissues = [3]
        l.cond[0].distribution_type = 'beta'
        v = np.zeros((255, 255, 255, 6))
        for i in range(6):
            v[:, :, :, i] = i
        affine = np.array([[-1, 0, 0, 128],
                           [0, 1, 0, -128],
                           [0, 0, 1, -127],
                           [0, 0, 0, 1]])
        l.anisotropy_vol = v
        l.anisotropy_affine = affine

        with tempfile.NamedTemporaryFile(suffix='.hdf5') as f:
            fn_hdf5 = f.name
        l._write_conductivity_to_hdf5(fn_hdf5)
        l2 = sim_struct.SimuList()
        l2._get_conductivity_from_hdf5(fn_hdf5)

        assert np.isclose(l2.cond[3].value, 2)
        assert l2.cond[46].value is None
        assert l2.cond[0].name == 'WM'
        assert l2.cond[46].name is None
        assert np.all(l2.anisotropic_tissues == [3])
        assert l2.anisotropy_type == 'vn'
        assert np.allclose(l2.anisotropy_affine, affine)
        assert np.allclose(l2.anisotropy_vol, v)
        assert l2.cond[0].distribution_type == 'beta'
        assert l2.cond[1].distribution_type is None

        os.remove(fn_hdf5)





class TestPoslist:
    def test_poslist_csv(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'CoilPos,1,2,3,0,0,1,0,1,0,10,coil,comment\n')
        f.write(b'CoilPos,3,2,1,1,0,0,0,1,0,10,coil,comment\n')
        f.close()
        p = sim_struct.TMSLIST()
        p.add_positions_from_csv(f.name)
        assert np.allclose(p.pos[0].matsimnibs, [[1, 0, 0, 1],
                                                 [0, 1, 0, 2],
                                                 [0, 0, 1, 3],
                                                 [0, 0, 0, 1]])
        assert np.allclose(p.pos[1].matsimnibs, [[0, 0, 1, 3],
                                                 [0, 1, 0, 2],
                                                 [-1, 0, 0, 1],
                                                 [0, 0, 0, 1]])

class TestPosition:
    def test_position_read_mat_matsimnibs(self, mat_session):
        pos = sim_struct.POSITION(mat_session['poslist'][0][0][0]['pos'][0][0][0])
        assert pos.matsimnibs == [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]

    def test_position_read_didt(self, mat_session):
        pos = sim_struct.POSITION(mat_session['poslist'][0][0][0]['pos'][0][0][0])
        assert pos.didt == 123

class TestCond:
    def test_set_distribution(self):
        c = sim_struct.COND()
        c.distribution_type = 'uniform'
        assert True

    def test_set_distribution_fail(self):
        c = sim_struct.COND()
        with pytest.raises(ValueError):
            c.distribution_type = 'b'

class TestTryToReadMat:
    def test_try_to_read_string(self):
        mat = {}
        mat['string'] = ['abc']
        r = sim_struct.try_to_read_matlab_field(mat, 'string', str)
        assert r == 'abc'

    def test_try_to_read_wrong_key(self):
        mat = {}
        mat['s'] = ['abc']
        r = sim_struct.try_to_read_matlab_field(mat, 'st', str)
        assert r is None

class TestLeadfield:
    def test_prepare_find_m2m(self, sphere3_fn):
        l = sim_struct.LEADFIELD()
        l.pathfem = ''
        l.fnamehead = sphere3_fn
        path, n = os.path.split(sphere3_fn)
        dir_fn = os.path.join(path, 'm2m_'+n[:-4])
        os.mkdir(dir_fn)
        l._prepare()
        assert os.path.abspath(l.subpath) == os.path.abspath(dir_fn)
        os.rmdir(dir_fn)
        l = sim_struct.LEADFIELD()
        l.pathfem = ''
        l.fnamehead = sphere3_fn
        path, n = os.path.split(sphere3_fn)
        dir_fn = os.path.join(path, 'm2m_'+n[:-4])
        os.mkdir(dir_fn)
        l._prepare()
        assert os.path.abspath(l.subpath) == os.path.abspath(dir_fn)
        os.rmdir(dir_fn)


class TestTDCSLEADFIELD:
    def test_tdcslist_csv(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        f.write(b'Fiducial,1,2,3,fiducial,a2\n')
        f.write(b'Electrode,1.2,2.4,3.6,9,8,7,electrode,b2,b3\n')
        f.write(b'ReferenceElectrode,1.2,2.4,7.1\n')
        f.close()
        p = sim_struct.TDCSLEADFIELD()
        p.eeg_cap = f.name
        p._add_electrodes_from_cap()

        assert np.allclose(p.electrode[1].centre, [1.2,2.4,3.6])
        assert np.allclose(p.electrode[1].pos_ydir, [9, 8, 7])
        assert np.allclose(p.electrode[0].centre, [1.2,2.4,7.1])
        assert p.electrode[0].pos_ydir == []

        p = sim_struct.TDCSLEADFIELD()
        p.eeg_cap = f.name
        p.electrode.shape = 'ellipse'
        p._add_electrodes_from_cap()

        assert np.allclose(p.electrode[1].centre, [1.2, 2.4, 3.6])
        assert np.allclose(p.electrode[1].pos_ydir, [9, 8, 7])
        assert p.electrode[1].shape == 'ellipse'

        assert np.allclose(p.electrode[0].centre, [1.2,2.4,7.1])
        assert p.electrode[0].pos_ydir == []
        assert p.electrode[0].shape == 'ellipse'

        p = sim_struct.TDCSLEADFIELD()
        p.eeg_cap = f.name
        p.electrode = [sim_struct.ELECTRODE(), sim_struct.ELECTRODE()]
        p.electrode[0].shape = 'ellipse'
        p.electrode[1].shape = 'rect'
        p._add_electrodes_from_cap()

        assert np.allclose(p.electrode[1].centre, [1.2, 2.4, 3.6])
        assert np.allclose(p.electrode[1].pos_ydir, [9, 8, 7])
        assert p.electrode[1].shape == 'ellipse'

        assert np.allclose(p.electrode[0].centre, [1.2,2.4,7.1])
        assert p.electrode[0].pos_ydir == []
        assert p.electrode[0].shape == 'rect'


class TestGetSurroundPos:
    def test_get_surround_pos(self,sphere3_fn):
        P = sim_struct.get_surround_pos([0, 0, 95], sphere3_fn,
                                        radius_surround = 94.5*np.pi, N = 3,
                                        pos_dir_1stsurround = [10, 0, 95])
        P = np.asarray(P)
        assert np.max(np.abs(np.diff(P,axis=0))) < 0.1
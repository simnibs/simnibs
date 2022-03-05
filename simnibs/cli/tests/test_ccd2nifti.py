import os
import numpy as np
import nibabel as nib
import pytest

from ... import SIMNIBSDIR
from .. import ccd2nifti

@pytest.fixture(scope='module')
def testcoil_ccd():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'testcoil.ccd')
    return fn

@pytest.fixture(scope='module')
def testcoil_nii():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'testcoil.nii.gz')
    return nib.load(fn)

@pytest.fixture(scope='module')
def testcoil_niiB():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'testcoilB.nii.gz')
    return nib.load(fn)

class TestCcd2nifti:
    def test_ccd2nifti(self,testcoil_ccd,testcoil_nii,testcoil_niiB):
        nii = ccd2nifti.ccd2nifti(testcoil_ccd)
        assert np.allclose(nii.dataobj, testcoil_nii.dataobj)
        assert np.allclose(nii.affine[:3,:3],10*np.identity(3))
        assert np.allclose(nii.affine[:4,3],(-100,-100,-100,1))
        assert nii.header['descrip']==b'dIdtmax=100.0;coilname=Test coil;stimulator=None;brand=None'
        niiB = ccd2nifti.ccd2nifti(testcoil_ccd, Bfield=True)
        assert np.allclose(niiB.dataobj, testcoil_niiB.dataobj)

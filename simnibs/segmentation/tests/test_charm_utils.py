import numpy as np
import pytest
import nibabel as nib
import os
from ... import SIMNIBSDIR
from .. import charm_utils

@pytest.fixture(scope='module')
def testernie_nii():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'ernie_T1_ds5.nii.gz')
    return fn

@pytest.fixture(scope='module')
def testmni_nii():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'MNI_test_ds5.nii.gz')
    return fn


def generate_label_arr(ndim=2):
    """Generate an array of labels."""
    size = ndim * (10,)
    arr = np.zeros(size, dtype=np.uint16)
    arr[5:] = 1
    arr[0, :5] = 2
    arr[5, 5] = 3
    arr[7:, 7:] = 5
    return arr


# all cases invoke remapping!
test_label_unassigned_elements_inputs = (
    (2, [0, 1, 2, 3, 5], 3, None, 0),
    (3, None, 3, None, 1),
    (3, None, 3, [1], 0),
    (5, None, 3, None, 1),
    (5, None, 3, [1], 5),  # should show a warning
)


# @pytest.mark.filterwarnings("ignore:Some elements could not be labeled")
@pytest.mark.parametrize(
    "label_unassign, labels, window_size, ignore_labels, expected_label",
    test_label_unassigned_elements_inputs,
)
def test_label_unassigned_elements(
    label_unassign,
    labels,
    window_size,
    ignore_labels,
    expected_label,
):
    label_arr = generate_label_arr(3)
    expected_arr = label_arr.copy()
    expected_arr[label_arr == label_unassign] = expected_label
    np.testing.assert_allclose(
        expected_arr,
        charm_utils.label_unassigned_elements(
            label_arr, label_unassign, labels, window_size, ignore_labels
        ),
    )

def test_sanlm(tmpdir, testernie_nii):
    denoised_scan = tmpdir.mkdir("denoised").join("denoised.nii.gz")
    input_scan = nib.load(testernie_nii)
    input_data = input_scan.get_fdata()
    charm_utils._denoise_input_and_save(testernie_nii, denoised_scan)
    output_scan = nib.load(denoised_scan)
    output_data = output_scan.get_fdata()
    assert input_data.var() > output_data.var()

def test_mni_affine(tmpdir, testmni_nii):
    trans_scan_name = tmpdir.mkdir("shifted").join("shifted_MNI.nii.gz")
    input_scan = nib.load(testmni_nii)
    trans_mat = np.eye(4)
    trans_mat[:3, 3] = -10
    trans_affine = trans_mat@input_scan.affine
    trans_mni = nib.Nifti1Image(input_scan.get_fdata(), trans_affine)
    nib.save(trans_mni, trans_scan_name)
    trans_settings = {"translation_scale": -100,
                      "max_iter": 10,
                      "shrink_factors": [0],
                      "smoothing_factors": [4.0],
                      "center_of_mass": True,
                      "samp_factor": 1.0,
                      "bg_value": 0}
    RAS2LPS = np.diag([-1, -1, 1, 1])
    estimated_trans_mat = charm_utils._init_atlas_affine(str(trans_scan_name),
                                                         testmni_nii,
                                                         trans_settings)
    np.testing.assert_allclose(trans_mat,
                               RAS2LPS@estimated_trans_mat@RAS2LPS)

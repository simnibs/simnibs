import pytest
import os
import numpy as np
import tempfile

from ... import SIMNIBSDIR
from ...simulation.sim_struct import TMSLIST, POSITION
from ..nnav import localite, softaxic, brainsight


FIXTURE_DIR = os.path.join(
    SIMNIBSDIR, 
    "_internal_resources", 
    "testing_files", 
    "nnav_testdata"
)

## LOCALITE test files and data ##
@pytest.fixture()
def tm_fn():
    fn = os.path.join(FIXTURE_DIR,
                      "localite",
                      "Session_20120925103137328/TMSTrigger/"
                      "TriggerMarkers_Coil1_20210409170817799.xml".replace('/', os.sep))
    return fn


@pytest.fixture()
def tm_no_pos_fn():
    fn = os.path.join(FIXTURE_DIR,
                      "localite",
                      "Session_20120925103137328/TMSTrigger/"
                      "no_pos_data_TriggerMarkers_Coil0_20200220155529299.xml".replace('/', os.sep))
    return fn


@pytest.fixture()
def im_lps_fn1():
    """Localite data for LPS oriented (-> DICOM) T1."""
    fn = os.path.join(FIXTURE_DIR,
                      "localite",
                      "Session_20120925103137328/InstrumentMarkers/"
                      "LPS_InstrumentMarker20120925123537421.xml".replace('/', os.sep))
    return fn


@pytest.fixture()
def im_ras_fn1():
    """Localite data for RAS oriented T1."""
    fn = os.path.join(FIXTURE_DIR,
                      "localite",
                      "Session_20120925103137328/InstrumentMarkers/"
                      "RAS_InstrumentMarker20170316160859032.xml".replace('/', os.sep))
    return fn


@pytest.fixture()
def im_ras_fn2():
    """Localite data for RAS oriented T1."""
    fn = os.path.join(FIXTURE_DIR,
                      "localite",
                      "Session_20120925103137328/InstrumentMarkers/"
                      "RAS_InstrumentMarker20170403154022328.xml".replace('/', os.sep))
    return fn


@pytest.fixture
def arr_tm_0():
    return np.array([[-6.18652573e-01, 7.39047765e-01, 2.66519129e-01, -2.08740412e+01],
                     [7.54371396e-01, 6.53553556e-01, -6.12094487e-02, -2.76457908e+01],
                     [-2.19426163e-01, 1.63190693e-01, -9.61860701e-01, 1.02090554e+02],
                     [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00]])


@pytest.fixture
def arr_im_0():
    return np.array([[-4.54280857e-01, +5.15500953e-01, +7.26594267e-01, -4.18376747e+01],
                     [+7.82924568e-01, +6.20186401e-01, +4.94924946e-02, -3.41115599e+01],
                     [-4.25099496e-01, +5.91336751e-01, -6.85319504e-01, +8.86091674e+01],
                     [+0.00000000e+00, +0.00000000e+00, -0.00000000e+00, +1.00000000e+00]])


@pytest.fixture
def arr_im_lps_0():
    return np.array([[-7.47832761e-01, +5.64773443e-01, +3.49056141e-01, -3.62375709e+01],
                     [+6.05764663e-01, +7.95614077e-01, +1.05110853e-02, -3.56913917e+01],
                     [-2.71769027e-01, +2.19299494e-01, -9.37076517e-01, +1.01916610e+02],
                     [+0.00000000e+00, +0.00000000e+00, +0.00000000e+00, +1.00000000e+00]])


## SOFTAXIC test files and data ##
@pytest.fixture()
def softaxicTestFile():
    fn = os.path.join(FIXTURE_DIR,'softaxic','neuronav_coil_test.stmpx')
    return fn


@pytest.fixture
def pos0():
    pos=np.array([[ -0.335503 ,  -0.171625 ,   0.926273 , -51.3255   ],
           [ -0.234913 ,   0.967445 ,   0.0941655,   2.26082  ],
           [ -0.912279 ,  -0.186001 ,  -0.364897 ,  68.1565   ],
           [  0.       ,   0.       ,   0.       ,   1.       ]])
    return pos


## BRAINSIGHT test files and data ##
@pytest.fixture()
def dcm_brainsight_fn():
    # Brainsight export for DICOM file, brainsight world coordinate system
    dcm_fn = os.path.join(FIXTURE_DIR,
                          "brainsight",
                          "dcmimage_brainsightcoord.txt".replace('/', os.sep))
    return dcm_fn


@pytest.fixture()
def dcm_world_fn():
    # Brainsight export for DICOM file, dicom coordinate system
    dcm_fn = os.path.join(FIXTURE_DIR,
                          "brainsight",
                          "dcmimage_worldcoord.txt".replace('/', os.sep))
    return dcm_fn


@pytest.fixture()
def nii_brainsight_fn():
    # Brainsight export for nifti file, brainsight coordinate system
    nii_fn = os.path.join(FIXTURE_DIR,
                          "brainsight",
                          "niiImage_brainsightcoord.txt".replace('/', os.sep))
    return nii_fn


@pytest.fixture()
def nii_nifti_fn():
    # Brainsight export for nifti file, nifti coordinate system
    nii_fn = os.path.join(FIXTURE_DIR,
                          "brainsight",
                          "niiImage_nifticoord.txt".replace('/', os.sep))
    return nii_fn


@pytest.fixture()
def nii_nifti_empty_fn():
    # Brainsight export for nifti file, nifti coordinate system
    # no targets, not samples
    nii_fn = os.path.join(FIXTURE_DIR,
                          "brainsight",
                          "niiImage_nifticoord_empty.txt".replace('/', os.sep))
    return nii_fn


@pytest.fixture()
def nii_nifti_no_samples_fn():
    # Brainsight export for nifti file, nifti coordinate system
    # not samples
    nii_fn = os.path.join(FIXTURE_DIR,
                          "brainsight",
                          "niiImage_nifticoord_no_samples.txt".replace('/', os.sep))
    return nii_fn


@pytest.fixture()
def nii_nifti_no_targets_fn():
    # Brainsight export for nifti file, nifti coordinate system
    # no targets
    nii_fn = os.path.join(FIXTURE_DIR,
                          "brainsight",
                          "niiImage_nifticoord_no_targets.txt".replace('/', os.sep))
    return nii_fn


@pytest.fixture()
def arr_sample1():
    return np.array([[-0.534, 0.553, 0.64, -46.889],
                     [0.529, 0.808, -0.258, 44.464],
                     [-0.66, 0.201, -0.724, 73.812],
                     [0., 0., 0., 1.]])


@pytest.fixture()
def arr_sample2():
    return np.array([[-7.1200e-01, 3.2100e-01, 6.2400e-01, -4.4198e+01],
                     [-1.3000e-02, 8.8300e-01, -4.6900e-01, 7.5810e+01],
                     [-7.0200e-01, -3.4200e-01, -6.2400e-01, 6.5742e+01],
                     [0.0000e+00, 0.0000e+00, 0.0000e+00, 1.0000e+00]])



class TestLocalite:
    def test_read_tm(self, tm_fn, arr_tm_0):
        """Read TriggerMarker file"""
        poslist = localite().read(tm_fn)
        assert len(poslist.pos) == 8
        assert np.all(np.isclose(poslist.pos[0].matsimnibs, arr_tm_0))
        assert poslist.pos[7].didt == 62000000

    def test_read_tm_no_pos(self, tm_no_pos_fn, arr_tm_0):
        """Read TriggerMarker file that has some positions with untracked coil"""
        with pytest.warns(UserWarning):
            poslist = localite().read(tm_no_pos_fn)
        assert len(poslist.pos) == 1
        assert np.all(np.isclose(poslist.pos[0].matsimnibs, arr_tm_0))
        assert poslist.pos[0].didt == 62000000

    def test_read_im(self, im_ras_fn1, arr_im_0):
        """Read InstrumentMarker file"""
        poslist = localite().read(im_ras_fn1)

        assert len(poslist.pos) == 6
        assert poslist.pos[0].name == 'M0'
        assert poslist.pos[0].didt == 1e6
        assert np.all(np.isclose(poslist.pos[0].matsimnibs, arr_im_0))

    def test_read_im_lps(self, im_lps_fn1, arr_im_lps_0):
        """Read LPS file"""
        poslist = localite().read(im_lps_fn1)
        assert len(poslist.pos) == 4
        assert(poslist.pos[0].name == 'hotspot')
        assert np.all(np.isclose(poslist.pos[0].matsimnibs, arr_im_lps_0))

    def test_read_multiple_files(self, im_ras_fn1, im_ras_fn2):
        """Read InstrumentMarkers from multiple files"""
        poslist_list = localite().read([im_ras_fn1, im_ras_fn2])
        assert len(poslist_list) == 2

    def test_write_arr(self, arr_tm_0):
        """Write InstrumentMarker file from matsimnibs"""
        with tempfile.TemporaryDirectory() as folder:
            fn = os.path.join(folder, 'im_from_arr.xml')

            localite().write(arr_tm_0, fn)
            poslist_from_xml = localite().read(fn, markertype='InstrumentMarker')

            assert np.all(np.isclose(poslist_from_xml.pos[0].matsimnibs, arr_tm_0))

    def test_write_pos(self, arr_tm_0):
        """Write InstrumentMarker file from POSITION()"""
        pos = POSITION()
        pos.matsimnibs = arr_tm_0
        with tempfile.TemporaryDirectory() as folder:
            fn = os.path.join(folder, 'im_from_pos.xml')

            localite().write(pos, fn)
            poslist_from_xml = localite().read(fn, markertype='InstrumentMarker')

            assert np.all(np.isclose(poslist_from_xml.pos[0].matsimnibs, arr_tm_0))

    def test_write_tmslist(self, arr_tm_0, arr_im_0):
        """Write InstrumentMarker file from TMSLIST()"""
        tmslist = TMSLIST()
        for pos in [arr_tm_0, arr_im_0]:
            p = POSITION()
            p.matsimnibs = pos
            tmslist.add_position(p)

        with tempfile.TemporaryDirectory() as folder:
            fn = os.path.join(folder, 'im_from_tmslist.xml')

            localite().write(tmslist, fn)
            poslist_from_xml = localite().read(fn, markertype='InstrumentMarker')

            for pos_org, pos_xml in zip(tmslist.pos, poslist_from_xml.pos):
                assert np.all(np.isclose(pos_org.matsimnibs, pos_xml.matsimnibs))


class TestSoftaxic():
    def test_read(self, softaxicTestFile, pos0):
        poslist = softaxic().read(softaxicTestFile)
        assert np.allclose(poslist.pos[0].matsimnibs, pos0)
        

class TestBrainsight:
    def test_read_dcm_brainsight(self, dcm_brainsight_fn):
        # No support for brainsight NLR
        with pytest.raises(ValueError) as e:
            brainsight().read(dcm_brainsight_fn)

    def test_read_dcm_world(self, dcm_world_fn, arr_sample1):
        tmslist_targets, tms_list_samples = brainsight().read(dcm_world_fn)
        assert tmslist_targets.name == 'Targets'
        assert tms_list_samples.name == 'Samples'
        assert len(tmslist_targets.pos) == 1
        assert len(tms_list_samples.pos) == 18
        assert tms_list_samples.pos[0].name == 'Sample 1'
        np.all(np.isclose(tms_list_samples.pos[0].matsimnibs, arr_sample1))

    def test_read_nii_brainsight(self, nii_brainsight_fn):
        # No support for brainsight NLR
        with pytest.raises(ValueError) as e:
            brainsight().read(nii_brainsight_fn)

    def test_read_nii_nifti(self, nii_nifti_fn, arr_sample2):
        tmslist_targets, tms_list_samples = brainsight().read(nii_nifti_fn)
        assert tmslist_targets.name == 'Targets'
        assert tms_list_samples.name == 'Samples'
        assert len(tmslist_targets.pos) == 7
        assert len(tms_list_samples.pos) == 101
        assert tms_list_samples.pos[0].name == 'Sample 1'
        np.all(np.isclose(tms_list_samples.pos[0].matsimnibs, arr_sample2))

    def test_read_nii_nifti_empty(self, nii_nifti_empty_fn):
        with pytest.raises(ValueError) as e:
            tmslist_targets, tms_list_samples = brainsight().read(nii_nifti_empty_fn)

    def test_read_nii_nifti_no_samples(self, nii_nifti_no_samples_fn):
        tmslist_targets, tmslist_samples = brainsight().read(nii_nifti_no_samples_fn)
        assert len(tmslist_targets.pos) == 7
        assert len(tmslist_samples.pos) == 0

    def test_read_nii_nifti_no_targets(self, nii_nifti_no_targets_fn):
        tmslist_targets, tmslist_samples = brainsight().read(nii_nifti_no_targets_fn)
        assert len(tmslist_targets.pos) == 0
        assert len(tmslist_samples.pos) == 101

    def test_write_nii_space(self, arr_sample2):
        pos = POSITION()
        pos.matsimnibs = arr_sample2
        with tempfile.TemporaryDirectory() as folder:
            print(folder)
            fn = os.path.join(folder, 'export.txt')

            brainsight().write(pos, fn, out_coord_space='NIfTI:Scanner')
            tmslist_targets, tmslist_samples = brainsight().read(fn)
            assert len(tmslist_targets.pos) == 1
            assert len(tmslist_samples.pos) == 0
            assert np.all(np.isclose(tmslist_targets.pos[0].matsimnibs, arr_sample2))

    def test_write_dicom_space(self, arr_sample1):
        pos = POSITION()
        pos.matsimnibs = arr_sample1
        with tempfile.TemporaryDirectory() as folder:
            print(folder)
            fn = os.path.join(folder, 'export.txt')

            brainsight().write(pos, fn, out_coord_space='World')
            tmslist_targets, tmslist_samples = brainsight().read(fn)
            assert len(tmslist_targets.pos) == 1
            assert len(tmslist_samples.pos) == 0
            assert np.all(np.isclose(tmslist_targets.pos[0].matsimnibs, arr_sample1))
            
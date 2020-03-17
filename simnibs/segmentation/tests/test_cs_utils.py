import os

import numpy as np
import scipy.io
import pytest


from simnibs.segmentation import _cs_utils

FUNCTION_TESTS_FOLDER = os.path.join(
    os.path.dirname(__file__),
   '..', '..', 'function_tests',
)


@pytest.fixture
def cube_image():
    img = np.zeros((50, 60, 70), dtype=np.float)
    img[10:40, 10:40, 10:40] = 1
    return img


class TestSanlm:
    '''
    def test_data(self):
        data = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_sanlm',
            'testdata_cat_sanlm.mat'
            ))
        Ym_filtered = _cs_utils.sanlm(data['Ym'], 3, 1)
        assert np.allclose(Ym_filtered, data['Yms'])
    '''
    def test_cube(self, cube_image):
        img = cube_image
        noisy = cube_image + 0.1*np.random.normal(size=cube_image.shape)
        filtered = _cs_utils.sanlm(noisy, 3, 1)
        assert np.linalg.norm(filtered-img) < np.linalg.norm(noisy-img)

class TestMedian3:
    '''
    def test_data(self):
        data_before = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_median3',
            'before_cat_vol_median3.mat'
            )
        )
        data_after = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_median3',
            'after_cat_vol_median3.mat'
            ))
        filtered = _cs_utils.cat_vol_median3(
            data_before['YcsfdM'],
            data_before['YM'],
            data_before['YM'])
        assert np.allclose(data_after['YcsfdM'], filtered)
    '''

    def test_cube(self, cube_image):
        img = cube_image
        noisy = cube_image + 0.1*np.random.normal(size=cube_image.shape)
        filtered = _cs_utils.cat_vol_median3(noisy)
        assert np.linalg.norm(filtered-img) < np.linalg.norm(noisy-img)

    def test_cube_mask(self, cube_image):
        img = cube_image
        noisy = cube_image + 0.1*np.random.normal(size=cube_image.shape)
        filtered = _cs_utils.cat_vol_median3(
            noisy,
            np.zeros_like(img),
            np.zeros_like(img)
        )
        assert np.allclose(noisy, filtered)


class TestVolEidist:
    '''
    def test_data(self):
        data = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_eidist',
            'input_cat_vol_eidist.mat'
            )
        )
        dist, _ = _cs_utils.cat_vol_eidist(data['YM'], data['F'], [1, 1, 1], 1, 1, 0, 1)
        ref = data['Ycsfd']
        assert np.all(np.isnan(ref) == np.isnan(dist))
        assert np.all(np.isinf(ref) == np.isinf(dist))
        assert np.allclose(dist[~np.isnan(dist) * ~np.isinf(dist)],
                           ref[~np.isnan(ref) * ~np.isinf(ref)])
    '''
    def test_simple(self):
        img = np.zeros((50, 60, 70), dtype=np.float)
        img[:25] = np.linspace(0, 1, 25)[:, None, None]
        img[25:] *= np.nan
        F = np.ones_like(img)
        dist, _ = _cs_utils.cat_vol_eidist(img, F, [1, 1, 1], 1, 1, 0, 0)
        assert np.allclose(dist[:13], np.arange(13)[::-1][:, None, None])
        assert np.allclose(dist[12:25], 0)
        assert np.all(np.isnan(dist[25:]))

class TestVolLocalstat:
    '''
    def test_data(self):
        in_ = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_localstat',
            'input_cat_vol_localstat.mat'
            )
        )
        out = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_localstat',
            'output_cat_vol_localstat.mat'
            )
        )

        Ygmts, _, _ = _cs_utils.cat_vol_localstat(in_['Ygmts'], in_['Ygmt1'] > 0, 1, 1)
        assert np.allclose(Ygmts, out['Ygmts'])
    '''
    def test_simple(self):
        img = np.zeros((50, 60, 70), dtype=np.float)
        img[:] = np.arange(50)[:, None, None]
        mask = np.zeros_like(img, dtype=bool)
        mask[25:] = True
        mean, _, _ = _cs_utils.cat_vol_localstat(img, mask, 1, 1)
        assert np.allclose(mean[:25], 0)
        assert np.allclose(mean[25, 1:-1, 1], (25 * 5 + 26*1)/6)
        assert np.allclose(mean[26:-1, 1:-1, 1], np.arange(26, 49)[:, None])

class TestVolPbtp:
    '''
    def test_data(self):
        data = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_pbtp',
            'testdata_cat_vol_pbtp.mat'
            )
        )
        Ygmt1, _ = _cs_utils.cat_vol_pbtp(data['Ymf'], data['Ywmd'], data['Ycsfd'])
        # Sometimes this test fails... why?
        assert np.allclose(Ygmt1, data['Ygmt1'])
    '''
    def test_simple(self):
        SEG = np.zeros((60, 50, 70), dtype=np.float)
        SEG[:20] = 1    # CSF
        SEG[20:40] = 2  # GM
        SEG[40:] = 3    # WM
        WMD = np.zeros_like(SEG)
        WMD[:40] = np.arange(40)[::-1, None, None]

        CSFD = np.zeros_like(SEG)
        CSFD[20:] = np.arange(40)[:, None, None]

        dist, _ = _cs_utils.cat_vol_pbtp(SEG, WMD, CSFD)
        assert np.allclose(dist[np.isclose(SEG, 1)], 0)
        assert np.allclose(dist[np.isclose(SEG, 3)], 0)
        assert np.allclose(dist[np.isclose(SEG, 2)], 20)

class TestVolGenus0:
    '''
    def test_data(self):
        data = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vol_genus0',
            'testdata_cat_vol_genus0.mat'
            ),
        )
        tmp, faces, vertices = _cs_utils.cat_vol_genus0(
            data['Yppt'], data['th_initial'])
        assert np.allclose(tmp, data['tmp'])
        assert np.allclose(faces, data['CS']['faces'][0, 0])
        assert np.allclose(vertices, data['CS']['vertices'][0, 0])
    '''
    def test_simple(self):
        x = np.linspace(-25, 25, 50)
        y = np.linspace(-30, 30, 60)
        z = np.linspace(-35, 35, 70)
        G = np.meshgrid(x, y, z)
        R = np.linalg.norm(G, axis=0)
        tmp, faces, vertices = _cs_utils.cat_vol_genus0((R < 10).astype(float), 0.9)
        assert np.allclose(tmp, R < 10)
        vertices -= [25, 30, 35]
        assert np.allclose(np.linalg.norm(vertices, axis=1), 10, atol=2)

class TestVbdist:
    '''
    def test_data(self):
        data = scipy.io.loadmat(os.path.join(
            FUNCTION_TESTS_FOLDER, 'cat_vbdist',
            'testdata_cat_vbdist.mat'
            ),
            squeeze_me=True
        )
        Ycsfd, _, _ = _cs_utils.cat_vbdist(
            data['Ymf'] < 1.5,
            data['Ymf'] > 1,
            data['vx_vol']
        )
        assert np.allclose(Ycsfd, data['Ycsfd'])
    '''
    def test_simple(self):
        P = np.zeros((50, 60, 70), dtype=bool)
        P[:30] = True
        R = np.zeros((50, 60, 70), dtype=bool)
        R[20:] = True
        dist, _, _ = _cs_utils.cat_vbdist(P, R)
        assert np.allclose(dist[30:], np.arange(1, 21)[:, None, None])


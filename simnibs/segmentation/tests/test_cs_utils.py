import os

import numpy as np
import scipy.io
import pytest


from .. import _cs_utils

@pytest.fixture
def cube_image():
    img = np.zeros((50, 60, 70), dtype=np.float)
    img[10:40, 10:40, 10:40] = 1
    return img


class TestSanlm:
    def test_cube(self, cube_image):
        img = cube_image
        noisy = cube_image + 0.1*np.random.normal(size=cube_image.shape)
        filtered = _cs_utils.sanlm(noisy, 3, 1)
        assert np.linalg.norm(filtered-img) < np.linalg.norm(noisy-img)

class TestMedian3:
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
    def test_simple(self):
        P = np.zeros((50, 60, 70), dtype=bool)
        P[:30] = True
        R = np.zeros((50, 60, 70), dtype=bool)
        R[20:] = True
        dist, _, _ = _cs_utils.cat_vbdist(P, R)
        assert np.allclose(dist[30:], np.arange(1, 21)[:, None, None])


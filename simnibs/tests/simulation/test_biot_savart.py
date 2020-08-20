import os
import numpy as np
import pytest
from simnibs import msh
import simnibs.simulation.biot_savart as biot_savart


@pytest.fixture
def sphere3_msh():
    fn = os.path.join(os.path.dirname(os.path.realpath(
        __file__)), '..', 'testing_files', 'sphere3.msh')
    return msh.read_msh(fn)

def _smooth_field(r, x, y, z):
    '''Defines a smooth field, usefull for testing'''
    # Transform from mm to m
    r = r * 1e-3
    x = x * 1e-3
    y = y * 1e-3
    z = z * 1e-3
    mask = (x**2 + y**2 + z**2) < r**2
    exp = np.exp(-r**2/np.abs(r**2 - (x**2+y**2+z**2)))*mask
    R = (r**2 - (x**2+y**2+z**2))

    Bx = (4*x*y**2*exp)/(-R)**2\
        - 3*x*exp/r**2\
        + (2*x*z**2*exp)/(-R)**2

    By = 3*y*exp/r**2\
        - (4*x**2*y*exp)/(-R)**2\
        - (2*y*z**2*exp)/(-R)**2

    Bz = (2*y**2*z*exp)/(-R)**2\
        - (2*x**2*z*exp)/(-R)**2

    B = np.stack((Bx, By, Bz), axis=-1)

    Jx = (14*y*z*exp)/(-R)**2\
        - (8*y**3*z*exp)/(-R)**3\
        - (4*y*z**3*exp*r**2)/(-R)**4\
        - (4*y**3*z*exp*r**2)/(-R)**4\
        - (8*y*z**3*exp)/(-R)**3\
        - (8*x**2*y*z*exp)/(-R)**3\
        - (4*x**2*y*z*exp*r**2)/(-R)**4

    Jy = (14*x*z*exp)/(-R)**2\
        - (8*x**3*z*exp)/(-R)**3\
        - (4*x*z**3*exp*r**2)/(-R)**4\
        - (4*x**3*z*exp*r**2)/(-R)**4\
        - (8*x*z**3*exp)/(-R)**3\
        - (8*x*y**2*z*exp)/(-R)**3\
        - (4*x*y**2*z*exp*r**2)/(-R)**4

    Jz = (16*x*y**3*exp)/(-R)**3\
        + (16*x**3*y*exp)/(-R)**3\
        + (8*x*y**3*exp*r**2)/(-R)**4\
        + (8*x**3*y*exp*r**2)/(-R)**4\
        - (28*x*y*exp)/(-R)**2\
        + (16*x*y*z**2*exp)/(-R)**3\
        + (8*x*y*z**2*exp*r**2)/(-R)**4

    J = np.stack((Jx, Jy, Jz), axis=-1)/biot_savart._mu0

    return J, B


class TestCompDomain:
    def test_mesh_larger(self, sphere3_msh):
        nvox = np.array([3, 3, 3])
        affine = np.array([[1, 0, 0, 1],
                           [0, 1, 0, 1],
                           [0, 0, 1, 1],
                           [0, 0, 0, 1]])
        domain = biot_savart._comp_domain(sphere3_msh, nvox, affine)
        assert np.allclose(domain, [[-95, -95, -95], [95, 95, 95]])

    @pytest.mark.parametrize('diag_entry', [1, -1])
    def test_target_larger(self, diag_entry, sphere3_msh):
        nvox = np.array([200, 200, 200])
        affine = np.array([[diag_entry, 0, 0, -100*diag_entry],
                           [0, 1, 0, -100],
                           [0, 0, 1, -100],
                           [0, 0, 0, 1]])
        domain = biot_savart._comp_domain(sphere3_msh, nvox, affine)
        assert np.allclose(domain, [[-100, -100, -100], [100, 100, 100]])

    @pytest.mark.parametrize('d', [0, 1, 2])
    def test_target_larger_dims(self, d, sphere3_msh):
        nvox = np.array([100, 100, 100])
        affine = np.array([[1, 0, 0, -50],
                           [0, 1, 0, -50],
                           [0, 0, 1, -50],
                           [0, 0, 0, 1]])
        nvox[d] = 200
        affine[d, 3] = -100
        domain = biot_savart._comp_domain(sphere3_msh, nvox, affine)
        res = np.array([[-95, -95, -95], [95, 95, 95]])
        res[:, d] = [-100, 100]
        assert np.allclose(domain, res)

    @pytest.mark.parametrize('d', [0, 1, 2])
    def test_target_shifted(self, d, sphere3_msh):
        nvox = np.array([100, 100, 100])
        affine = np.array([[1, 0, 0, -50],
                           [0, 1, 0, -50],
                           [0, 0, 1, -50],
                           [0, 0, 0, 1]])
        affine[d, 3] = 0
        domain = biot_savart._comp_domain(sphere3_msh, nvox, affine)
        res = np.array([[-95, -95, -95], [95, 95, 95]])
        res[1, d] = 100
        assert np.allclose(domain, res)


class TestVoxelize:
    @pytest.mark.parametrize('d', [100, 120])
    @pytest.mark.parametrize('r', [1., 2.])
    def test_voxelize(self, d, r, sphere3_msh):
        J = sphere3_msh.add_element_field(np.ones((sphere3_msh.elm.nr, 3)), 'J')
        domain = np.array([[-d, -d, -d], [d, d, d]])
        J_v, _ = biot_savart._voxelize(J, domain, r)
        grid_pos = np.meshgrid(
            *[np.arange(domain[0, i], domain[1, i], r) for i in range(3)],
            indexing='ij', sparse=True
        )
        R = np.sqrt(grid_pos[0]**2 + grid_pos[1]**2 + grid_pos[2]**2)
        '''
        plt.subplot(131)
        plt.imshow(J_v[..., 50, 0].T, origin='lower')
        plt.subplot(132)
        plt.imshow(J_v[..., 50, 1].T, origin='lower')
        plt.subplot(133)
        plt.imshow(J_v[50,..., 2].T, origin='lower')
        plt.show()
        '''

        assert J_v.shape == (2*d//r, 2*d//r, 2*d//r, 3)
        # Some points can be exacly in the voxel
        assert np.linalg.norm(J_v[R < 94]-1)/np.linalg.norm(J_v[R<94]) < 1e-2
        assert np.allclose(J_v[R > 96], 0)

class TestBSVol:
    def test_bs(self):
        N = 50  # Number of voxels
        r = 100.  # Radius
        res = 4.  # Resolution
        x = np.arange((-N/2)*res, (N/2)*res, res)
        x += 1/2
        x, y, z = np.meshgrid(x, x, x, indexing='ij', sparse=True)

        J, B = _smooth_field(r, x, y, z)

        Bh = biot_savart._bs_ft(J, res)

        Erms = np.linalg.norm(B - Bh)/np.linalg.norm(B)
        '''
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(3, 2)
        for i in range(3):
            axes[i, 0].imshow(B[..., N//2, i])
            axes[i, 1].imshow(Bh[..., N//2, i])
        plt.show()
        '''
        assert Erms < 1e-2

class TestCurl:
    def test_curl(self):
        N = 100  # Number of voxels
        r = 100.  # Radius
        res = 2.  # Resolution
        x = np.arange((-N/2)*res, (N/2)*res, res)
        x += 1/2
        x, y, z = np.meshgrid(x, x, x, indexing='ij', sparse=True)

        J, B = _smooth_field(r, x, y, z)
        J_curl = biot_savart._curl(B, res*1e-3)/biot_savart._mu0

        Erms = np.linalg.norm(J - J_curl)/np.linalg.norm(J)

        '''
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(3, 2)
        for i in range(3):
            im = axes[i, 0].imshow(J[..., N//2, i], cmap="coolwarm")
            fig.colorbar(im, ax=axes[i,0])
            im = axes[i, 1].imshow(J_curl[..., N//2, i], cmap="coolwarm")
            fig.colorbar(im, ax=axes[i,1])
        plt.show()
        '''
        assert Erms < .1

    def test_back_and_fourth(self):
        N = 100  # Number of voxels
        r = 100.  # Radius
        res = 2.  # Resolution
        x = np.arange((-N/2)*res, (N/2)*res, res)
        x += 1/2
        x, y, z = np.meshgrid(x, x, x, indexing='ij', sparse=True)

        J, B = _smooth_field(r, x, y, z)

        Bh = biot_savart._bs_ft(J, res)
        Jh = biot_savart._curl(Bh, res*1e-3)/biot_savart._mu0
        Erms = np.linalg.norm(J - Jh)/np.linalg.norm(J)
        '''
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(3, 2)
        for i in range(3):
            im = axes[i, 0].imshow(J[..., N//2, i], cmap="coolwarm")
            fig.colorbar(im, ax=axes[i,0])
            im = axes[i, 1].imshow(Jh[..., N//2, i], cmap="coolwarm")
            fig.colorbar(im, ax=axes[i,1])
        plt.show()
        print(Erms)
        '''
        assert Erms < .1

class Test_calcB:
    def test_calc_on_mesh(self, sphere3_msh):
        bar = sphere3_msh.elements_baricenters()
        J, _ = _smooth_field(
            95.,
            bar[:, 0],
            bar[:, 1],
            bar[:, 2]
        )
        J_elmdata = sphere3_msh.add_element_field(J, 'J')
        n_voxels = [200, 200, 200]
        affine = np.eye(4)
        affine[:3, 3] = [-100, -100, -100]

        B_grid = biot_savart.calc_B(J_elmdata, n_voxels, affine)

        x = np.arange(-100, 100, dtype=float)
        x += .1
        x, y, z = np.meshgrid(x, x, x, indexing='ij', sparse=True)
        _, B_ref = _smooth_field(95., x, y, z)

        Erms = np.linalg.norm(B_grid - B_ref)/np.linalg.norm(B_ref)
        '''
        import nibabel
        img = nibabel.Nifti1Image(B_grid, affine)
        nibabel.save(img, '/tmp/B_grid.nii.gz')
        img = nibabel.Nifti1Image(B_ref, affine)
        nibabel.save(img, '/tmp/B_ref.nii.gz')
        '''
        # The RMS is quite large because the J varies quickly and the mesh is coarse
        assert Erms < .5

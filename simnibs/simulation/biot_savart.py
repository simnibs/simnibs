# -*- coding: utf-8 -*-\
'''
    Functions for calculating magnetic fields
'''

import numpy as np
from ..utils import transformations
from ..mesh_tools import mesh_io

'''
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2019 Hassan Yazdanian, Guilherme B Saturnino


    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

_mu0 = 4 * np.pi * 1e-7

def calc_B(J, n_voxels, affine, calc_res=2., mask=False):
    ''' Calculates the magnetic field caused by a current density distribution J

    Parameters
    ----------
    J: simnibs.msh.ElementData
        ElementData field with the current density field
    n_voxels: list or tuple
        number of voxels in x, y, and z directions for the output B
    affine: ndarray
        A 4x4 matrix specifying the transformation from voxels to xyz for the output B
    calc_res: float (optional)
        Resolution of the space used for the FFT calculations, in mm. Default: 2
    mask: bool (optional)
        Wether to mask the resul using the mesh. Default: False

    Returns
    --------
    B: ndarray
        numpy array of size n_voxesx3 with the magnetic field        

    Reference
    ----------
    Yazdanian, H., Saturnino, G. B., Thielscher, A., & Knudsen, K. (2020).
    Fast evaluation of the Biot-Savart integral using FFT for electrical conductivity imaging.
    Journal of Computational Physics, 109408.
    https://doi.org/10.1016/j.jcp.2020.109408

    Example
    ---------
    Load a simulation result with a 'J' field can calculate B in an ROI
    Defined by a nifti file
    >>> import simnibs
    >>> import nibabel
    >>> simulation_J = simnibs.read_msh('simulation_result_with_J.msh')
    >>> reference_nifti = nibabel.load('reference.nii)
    >>> B = simnibs.simulation.calc_B(simulation_J.field['J'], reference_nifti.shape, reference_nifti.affine)
    '''
    mesh = J.mesh
    # Computational Domain (before zero-padding)
    domain = _comp_domain(mesh, n_voxels, affine)
    # Voxelize J
    J_vol, affine_vol = _voxelize(J, domain, calc_res)
    # Calculate B in the whole domain
    B_vol = _bs_ft(J_vol, calc_res)
    # Interpolate B to the domain of interest
    B = transformations.volumetric_affine(
        (B_vol, affine_vol), np.eye(4), affine, n_voxels,
        intorder=1, keep_vector_length=False
    )
    if mask:
        msk = mesh_io.ElementData(np.ones(mesh.elm.nr), mesh=mesh)
        msk = msk.interpolate_to_grid(n_voxels, affine, method='assign')
        B *= msk[..., None]

    return B


def _comp_domain(mesh, n_voxels, affine):
    ''' Compute the domain (before zero-padding) '''
    lim_mesh = np.array([np.min(mesh.nodes[:], axis=0), np.max(mesh.nodes[:], axis=0)])
    lim_t = np.array([affine[:3, 3], affine[:3, :3].dot(n_voxels) + affine[:3, 3]])
    lim_t = np.sort(lim_t, axis=0)
    domain = np.zeros((2, 3))
    domain[0, :] = np.minimum(lim_mesh[0, :], lim_t[0, :])
    domain[1, :] = np.maximum(lim_mesh[1, :], lim_t[1, :])
    return domain

def _voxelize(J, domain, res):
    ''' Voxelize a field given a domain and a resolution'''
    domain_size = domain[1] - domain[0]
    affine = np.eye(4)
    affine[:3, :3] *= res
    affine[:3, 3] = domain[0]
    nvox = np.ceil(domain_size/res).astype(int)
    nvox += nvox % 2
    J_v = J.interpolate_to_grid(nvox, affine, method='linear', continuous=False)
    return J_v, affine


def _bs_ft(J, res, pad=None):
    '''Calculates the denominator the biot-savart kernel. Outputs a real matrix'''
    if pad is None:
        pad = [d for d in J.shape[:3]]

    #pad = (2, 2, 2)
    kx, ky, kz = np.meshgrid(
        *[2*np.pi*np.fft.fftfreq(d + p, res * 1e-3) for d, p in zip(J.shape[:3], pad)],
        indexing='ij', sparse=True)

    kx[0] = 1e-9
    K = _mu0/(kx**2+ky**2+kz**2)
    K[0, 0, 0] = 0
    kx[0] = 0

    Jx_fft = np.fft.fftn(np.pad(J[..., 0], [(p//2, ) for p in pad],
                                'constant', constant_values=0))
    By = 1j*K*kz*Jx_fft
    Bz = -1j*K*ky*Jx_fft
    del Jx_fft

    Jy_fft = np.fft.fftn(np.pad(J[..., 1], [(p//2, ) for p in pad],
                                'constant', constant_values=0))
    Bx = -1j*K*kz*Jy_fft
    Bz += 1j*K*kx*Jy_fft
    del Jy_fft

    Jz_fft = np.fft.fftn(np.pad(J[..., 2], [(p//2, ) for p in pad],
                                'constant', constant_values=0))
    Bx += 1j*K*ky*Jz_fft
    By += -1j*K*kx*Jz_fft
    del Jz_fft

    Bx = np.real(np.fft.ifftn(Bx))[pad[0]//2:-pad[0]//2, pad[1]//2:-pad[1]//2, pad[2]//2:-pad[2]//2]
    By = np.real(np.fft.ifftn(By))[pad[0]//2:-pad[0]//2, pad[1]//2:-pad[1]//2, pad[2]//2:-pad[2]//2]
    Bz = np.real(np.fft.ifftn(Bz))[pad[0]//2:-pad[0]//2, pad[1]//2:-pad[1]//2, pad[2]//2:-pad[2]//2]

    B = np.stack((Bx, By, Bz), axis=-1)

    del Bx, By, Bz

    return B

def _curl(M, res):
    _, dMxdy, dMxdz = np.gradient(M[..., 0], res)
    dMydx, _, dMydz = np.gradient(M[..., 1], res)
    dMzdx, dMzdy, _ = np.gradient(M[..., 2], res)
    curlM = np.zeros_like(M)
    curlM[..., 0] = dMzdy - dMydz
    curlM[..., 1] = dMxdz - dMzdx
    curlM[..., 2] = dMydx - dMxdy
    return curlM

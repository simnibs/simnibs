'''
    Manipulation of conductivity information for SimNIBS
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2018 Guilherme Saturnino, Andre Antunes, Axel Thielscher

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.



'''


from __future__ import division
import warnings
import numpy as np
from . import sim_struct
from ..msh import mesh_io
from ..utils.simnibs_logger import logger


def standard_cond():
    S = []
    for i in range(1000):
        S.append(sim_struct.COND())

    # WM
    S[0].name = 'WM'
    S[0].value = 0.126
    S[0].descrip = 'brain white matter (from Wagner 2004)'

    # GM
    S[1].name = 'GM'
    S[1].value = 0.275
    S[1].descrip = 'brain gray matter (from Wagner 2004)'

    # CSF
    S[2].name = 'CSF'
    S[2].value = 1.654
    S[2].descrip = 'cerebrospinal fluid (from Wagner 2004)'

    # Bone
    S[3].name = 'Bone'
    S[3].value = 0.010
    S[3].descrip = 'average bone (from Wagner 2004)'

    # Scalp
    S[4].name = 'Scalp'
    S[4].value = 0.465
    S[4].descrip = 'average scalp (from Wagner 2004)'

    # Eye balls (vitreous humour)
    S[5].name = 'Eye_balls'
    S[5].value = 0.5
    S[5].descrip = 'vitreous humour (from Opitz, Paulus, Thielscher, submitted)'

    # Compact bone
    S[6].name = 'Compact_bone'
    S[6].value = 0.008
    S[6].descrip = 'compact bone (from Opitz, Paulus, Thielscher, submitted)'

    # Spongy bone
    S[7].name = 'Spongy_bone'
    S[7].value = 0.025
    S[7].descrip = 'spongy bone (from Opitz, Paulus, Thielscher, submitted)'

    # Blood
    S[8].name = 'Blood'
    S[8].value = 0.6
    S[8].descrip = 'Blood (from Gabriel et al, 2009)'

    # Muscle
    S[9].name = 'Muscle'
    S[9].value = 0.16
    S[9].descrip = 'Muscle (from Gabriel et al, 2009)'

    # Rubber
    S[99].name = 'Electrode_rubber'
    S[99].value = 29.4
    S[99].descrip = 'for tDCS rubber electrodes'

    # Saline
    S[499].name = 'Saline'
    S[499].value = 1.
    S[499].descrip = 'for tDCS sponge electrodes'

    return S


def cond2elmdata(mesh, cond_list, anisotropy_volume=None, affine=None,
                 aniso_tissues=[1, 2], correct_FSL=True, normalize=False,
                 excentricity_scaling=None, max_ratio=10, max_cond=2, correct_intensity=True):
    ''' Define conductivity ElementData from a conductivity list or anisotropy
    information
    Parameters
    ------------
    mesh: simnibs.mesh_io.Msh
        Mesh with geometry information
    cond_list: list
        List of conductivity values in each volume tag of the mesh. The first element
        corresponds to the conductivity value at the tag == 1, and so on
    anisotropy_volume: np.ndarray
        4-dimensional array with anisotropy information
    affine: np.ndarray
        4x4 matrix defining an affine transformation from the grid to the mesh space.
        Mandatory if anisotropy_volue is not None
    aniso_tissues: list
        List of tissues where anisotrpic conductivities should be applied. Default: 1 and 2
    correct_FSL: bool
        Wether to correct the tensors, if the data was preprocessed with FSL.
    normalize: bool
        Wether or not to normalize the tensors and use the conductivity provided in
        cond_list. This is equivalend to using the 'vn' method. Default: False
    excentricity_scaling: float between 0 and 1
        Wether to scale the excentricities. 0 will turn all tensors isotropic, 1 will
        make the 2 smallest eigenvalues 0, .5 will not change the tensors. Default: Do
        not scale excentricity
    max_ratio: positive float >= 1
        Maximum value for V1/V2 and V1/V3 where V1, V2, and V3 are the largest, middle
        and smallest eigenvalues of each tensor
    max_cond: Positive float
        Maximum eigenvalue for each tensor
    correct_intensity: bool
        Wether or not to fit the tensor sizes according to the scalar values (See
        Rullmann et. al. 2009). This procedure scales the entire tensor field with a
        single scalar. Does not run if normalize==True

    Returns:
    ------
    cond: simnibs.mesh_io.ElementData
        ElementData field with conductivity information
    '''
    try:
        aniso_tissues[0]
    except TypeError:
        aniso_tissues = [aniso_tissues]

    vol_tags = np.unique(mesh.elm.tag1[mesh.elm.elm_type == 4])

    def test_numerical(v, i):
        try:
            c = float(v)
        except (TypeError, ValueError):
            raise TypeError('The value {0} in cond_list is not numerical'.format(i))
        return c

    # Isotropic
    if anisotropy_volume is None:
        cond = mesh_io.ElementData(np.zeros(mesh.elm.nr, dtype=float),
                                'conductivity', mesh=mesh)
        for t in vol_tags:
            if len(cond_list) < t:
                raise ValueError('The cond_list size is too small'
                                 ', should be at least of size {0}'
                                 ''.format(np.max(vol_tags)))
            c = test_numerical(cond_list[t-1], t-1)
            cond.value[(mesh.elm.tag1 == t) * (mesh.elm.elm_type == 4)] = c
    # Anisotropic
    else:
        # Interpolate conductivity values
        assert len(anisotropy_volume.shape) == 4
        assert anisotropy_volume.shape[-1] == 6
        assert anisotropy_volume.shape[-1] == 6
        if affine is None:
            raise ValueError('Please define a 4x4 affine from the grid to the mesh space')

        cond = mesh_io.ElementData.from_data_grid(
                    mesh, anisotropy_volume, affine, 'conductivities',
                    order=1, cval=0.0, prefilter=True)
        cond.value = cond.value[:, [0, 1, 2, 1, 3, 4, 2, 4, 5]]
        tensors = cond.value.reshape(-1, 3, 3)
        if correct_FSL:
            M = affine[:3, :3] / np.linalg.norm(affine[:3, :3], axis=0)[:, None]
            R = np.eye(3)
            if np.linalg.det(M) > 0:
                R[0, 0] = -1
            M = M.dot(R)
            tensors = tensors.dot(M.T)
            tensors = M.dot(tensors.transpose(2, 1, 0)).transpose(2, 0, 1)
        cond.value = tensors.reshape(-1, 9)

        # VN type conductivities
        if normalize:
            for i, t in enumerate(vol_tags):
                # If tissue is isotropic
                if t not in aniso_tissues:
                    c = test_numerical(cond_list[t-1], t-1)
                    cond.value[(mesh.elm.tag1 == t) * (mesh.elm.elm_type == 4)] = \
                        [c, 0, 0, 0, c, 0, 0, 0, c]

                # If tissue is to be normalized
                else:
                    c = test_numerical(cond_list[t-1], t-1)
                    tensors = cond.value[(mesh.elm.tag1 == t) * (mesh.elm.elm_type == 4)]
                    tensors = tensors.reshape(-1, 3, 3)
                    tensors = _fix_zeros(tensors, c)
                    eigval, eigvec = _get_sorted_eigenv(tensors)
                    # Normalize
                    eigval /= (np.abs(eigval).prod(axis=1) ** (1./3.))[:, None]
                    # Fix
                    eigval = _fix_eigv(eigval, max_cond, max_ratio, c)
                    # Normalize again
                    eigval /= (eigval.prod(axis=1) ** (1./3.))[:, None]
                    # Fix again
                    eigval = _fix_eigv(eigval, max_cond, max_ratio, c)
                    if excentricity_scaling is not None:
                        eigval = _adjust_excentricity(eigval, excentricity_scaling)
                    tensors = _form_tensors(eigval, eigvec)
                    cond.value[(mesh.elm.tag1 == t) * (mesh.elm.elm_type == 4)] = \
                        tensors.reshape(-1, 9) * c

        # MC / DIR type conductivities
        else:
            mean_vol = np.zeros(len(vol_tags))
            eigval_list = [[] for i in vol_tags]
            eigvec_list = [[] for i in vol_tags]
            for i, t in enumerate(vol_tags):
                # If tissue is isotropic
                if t not in aniso_tissues:
                    c = test_numerical(cond_list[t-1], t-1)
                    cond.value[(mesh.elm.tag1 == t) * (mesh.elm.elm_type == 4)] = \
                        [c, 0, 0, 0, c, 0, 0, 0, c]

                else:
                    c = test_numerical(cond_list[t-1], t-1)
                    tensors = cond.value[(mesh.elm.tag1 == t) * (mesh.elm.elm_type == 4)]
                    tensors = tensors.reshape(-1, 3, 3)
                    eigval, eigvec = _get_sorted_eigenv(tensors)
                    if correct_intensity:
                        # First fix, does not apply upper bound
                        eigval = _fix_eigv(eigval, 1e10, max_ratio, -1e-6)
                        eigval_list[i] = eigval
                        eigvec_list[i] = eigvec
                    else:
                        eigval = _fix_eigv(eigval, max_cond, max_ratio, c)
                        if excentricity_scaling is not None:
                            eigval = _adjust_excentricity(eigval, excentricity_scaling)
                        tensors = _form_tensors(eigval, eigvec)
                        tensors = _fix_zeros(tensors, c)
                        cond.value[(mesh.elm.tag1 == t) * (mesh.elm.elm_type == 4)] = \
                            tensors.reshape(-1, 9)

                    mean_vol[i] = np.average(eigval.prod(axis=1))
        # Correct the intensity
        if correct_intensity and not normalize:
            num = 0
            denom = 0
            for i, t in enumerate(vol_tags):
                if t in aniso_tissues:
                    c = test_numerical(cond_list[t-1], t-1)
                    num += c * mean_vol[i] ** (1./3.)
                    denom += mean_vol[i] ** (2./3.)
            s = num / denom
            logger.info('Scaling the tensors with the factor: {0}'.format(s))
            # We need to fix the eigenvalues a second time
            for i, t in enumerate(vol_tags):
                if t in aniso_tissues:
                    c = test_numerical(cond_list[t-1], t-1)
                    eigval = s * eigval_list[i]
                    eigvec = eigvec_list[i]
                    eigval = _fix_eigv(eigval, max_cond, max_ratio, c)
                    logger.info(
                        'Reference conductivity: {0} Mean conductivity: {1}'.format(
                            c, np.average(eigval.prod(axis=1)) ** (1./3.)))
                    if excentricity_scaling is not None:
                        eigval = _adjust_excentricity(eigval, excentricity_scaling)
                    tensors = _form_tensors(eigval, eigvec)
                    tensors = _fix_zeros(tensors, c)
                    cond.value[(mesh.elm.tag1 == t) * (mesh.elm.elm_type == 4)] = \
                        tensors.reshape(-1, 9)

    cond.assign_triangle_values()

    return cond


def _get_sorted_eigenv(tensors):
    eig_val, eig_vec = np.linalg.eig(tensors)
    eig_val = np.real(eig_val)
    eig_vec = np.real(eig_vec)
    sort = eig_val.argsort(axis=1)[:, ::-1]
    eig_val = np.sort(eig_val, axis=1)[:, ::-1]
    s_eig_vec = np.zeros_like(tensors)
    for i in range(3):
        s_eig_vec[:, :, i] = eig_vec[np.arange(len(eig_vec)), :, sort[:, i]]
    eig_vec = s_eig_vec
    return eig_val, eig_vec


def _form_tensors(eigval, eigvec):
    # What we are doing: eig_vec[0].dot((ls[0, None] * eig_vec[0]).T)
    tensors = np.einsum('aij,akj -> aik',
                        eigvec, eigval[:, None] * eigvec)
    return tensors


def _adjust_excentricity(eig_val, scaling):
    '''' Tensor excentricity '''
    if scaling >= 1. or scaling < 0.:
        raise ValueError('Invalid scaling factor for '
                         'tensor excentricity: {0}'.format(scaling))

    if np.any(eig_val) < 0:
        raise ValueError('Found a negative eigenvalue!')
    if np.any(np.isclose(eig_val, 0)):
        raise ValueError('Found a zero eigenvalue!')
    # Excentricity
    ex = np.sqrt(1 - (eig_val[:, [1, 2, 2]]/eig_val[:, [0, 0, 1]]) ** 2)
    # Scale excentricity
    if scaling < .5:
        es = 2.0 * ex * scaling
    elif scaling > .5:
        es = 2.0 * (1 - ex) * scaling + 2 * ex - 1
    else:
        es = ex
    # New eigenvalues
    ls = np.ones_like(eig_val)
    ls[:, 1] = np.sqrt(1 - es[:, 0]**2)
    ls[:, 2] = np.sqrt(1 - es[:, 1]**2)
    # Scale the new eigenvalues so that they keep the volume of the tensor
    ls *= (np.prod(eig_val, axis=1)/np.prod(ls, axis=1))[:, None]**(1.0/3.0)
    # If the tensor is isotropic, don't scale it
    iso = np.isclose(eig_val[:, 0], eig_val[:, 2], rtol=1e-2)
    ls[iso] = eig_val[iso]
    return ls


def _fix_zeros(tensors, c):
    tensors = tensors.reshape(-1, 9)
    negative = np.all(np.isclose(tensors, 0), axis=1)
    frac = np.sum(negative)/len(tensors)
    if frac > .1:
        log = logger.critical
    else:
        log = logger.info
    log('Found {0} ({1:.1%}) Zero tensors in the volume. '
        ' Fixing it.'.format(np.sum(negative), frac))
    tensors[negative] = c * np.eye(3).reshape(-1)
    tensors = tensors.reshape(-1, 3, 3)
    return tensors


def _fix_eigv(eig_val, max_value, max_ratio, c):
    negative = np.all(eig_val <= 0.0, axis=1)
    eig_val[negative] = c
    frac = np.sum(negative)/len(eig_val)
    if frac > .1:
        log = logger.critical
    else:
        log = logger.info
    log('Found {0} ({1:.1%}) Negative Semi-definite tensors in the volume. '
        ' Fixing it.'.format(np.sum(negative), frac))

    large = eig_val > max_value
    eig_val[large] = max_value

    frac = np.sum(large)/(3. * len(eig_val))
    if frac > .1:
        log = logger.critical
    else:
        log = logger.info
    log('Found {0} ({1:.1%}) too large eigenvalues in the volume. '
        ' Fixing it.'.format(np.sum(large), frac))

    small = eig_val < (eig_val[:, 0] / max_ratio)[:, None]
    eig_val[small[:, 1], 1] = eig_val[small[:, 1], 0] / max_ratio
    eig_val[small[:, 2], 2] = eig_val[small[:, 2], 0] / max_ratio
    frac = np.sum(small)/(3 * len(eig_val))
    if frac > .1:
        log = logger.critical
    else:
        log = logger.info
    log('Found {0} ({1:.1%}) too small or negative eigenvalues in the volume. '
        ' Fixing it.'.format(np.sum(small), frac))

    return eig_val


def TensorVisualization(cond, mesh, all_compoents=False):
    ''' Creates a visualization of the tensors for plotting

    Parameters
    -----------
    cond: Nx9 numpy array
        Array with conductivity values at each element
    mesh: simnibs.msh.Msh
        Mesh file
    all_componets: bool (optional)
        Whether or nor to plot the middle and minimum eigenvalue, as well as the mean
        conductivity

    Returns
    ---------
    data: list
        List with the fields for visualization
    '''
    assert cond.nr_comp == 9, 'This function can only be used with tensor conductivities'
    tensors = cond.value.reshape(-1, 3, 3)
    eig_val, eig_vec = _get_sorted_eigenv(tensors)
    mc = np.prod(eig_val, axis=1) ** (1.0/3.0)
    data = []
    data.append(
        mesh_io.ElementData(
            eig_vec[:, :, 0] * eig_val[:, 0, None],
            name='max_eig', mesh=mesh))
    if all_compoents:
        data.append(
            mesh_io.ElementData(
                eig_vec[:, :, 1] * eig_val[:, 1, None],
                name='mid_eig', mesh=mesh))
        data.append(
            mesh_io.ElementData(
                eig_vec[:, :, 2] * eig_val[:, 2, None],
                name='min_eig', mesh=mesh))

        data.append(
            mesh_io.ElementData(mc,
                 name='mean_conductivity', mesh=mesh))
    return data

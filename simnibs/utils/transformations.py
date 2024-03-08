# -*- coding: utf-8 -*-\
'''
    linear and non-linear transfomations for SimNIBS
'''
'''
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.
    Copyright (C) 2020 Kristoffer H Madsen, Guilherme B Saturnino, Axel Thielscher

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

from pathlib import Path
import warnings
import os
from functools import partial
import gc
import copy
import nibabel as nib
import numpy as np
import scipy.ndimage
import scipy.spatial
from scipy.sparse import coo_matrix, csr_matrix
from typing import Union

from simnibs.utils.mesh_element_properties import ElementTags
from ..utils.simnibs_logger import logger
from ..utils.file_finder import templates, SubjectFiles, get_reference_surf
from ..utils.csv_reader import write_csv_positions, read_csv_positions

__all__ = [
    'warp_volume',
    'warp_coordinates',
    'subject2mni_coords',
    'mni2subject_coords',
    'subject_atlas',
    'middle_gm_interpolation'

]


def volumetric_nonlinear(image, deformation, target_space_affine=None,
                         target_dimensions=None, intorder=1, cpus=1,
                         inverse_deformation=None,
                         keep_vector_length=True):
    ''' Applies a volumetric non-linear transformation

    Parameters
    --------
    image: tuple
        tuple (ndarray, affine) with the data to be transformed
    deformation: tuple
        tuple (ndarray, affine) with defomation field. The deformation field is in the
        target space and specifies in each voxel the equivalent coordinates (x, y, z)
        of the original space
    target_space_affine: ndarray (Optional)
        4x4 matrix with the voxel indices to (x, y, z) coordinates in the target space.
        Default: same as the deformation field
    target_dimensions: array (Optional)
        Dimensions of volume in the final image. Default: s
    intorder: int (Optiona)
        Order of interpolation. See the documentation for the
        scipy.ndimage.map_coordinates function. Default=1
    inverse_deformation: tuple
        tuple (ndarray, affine) with inverse field. The inverse deformation field is in the
        original space and specifies in each voxel the equivalent coordinates (x, y, z)
        of the target space. This is used ony to rotate vectors accordingly. If not set,
        vectors wil not be rotated
    keep_vector_length: bool
        Whether to keep the length of the vectors unchanged. Only comes into effect if the
        image is 3 dimensional and inverse_deformation is set.
    Returns
    ------
    outdata: ndarray
        ndarray corresponding to the transormed image
    '''
    # read volumetric data
    im_data, im_affine = image
    df_data, df_affine = deformation
    if len(df_data.shape) > 4:
        df_data = df_data.squeeze()

    # If the resolution is to be changed
    if target_dimensions is None:
        target_dimensions = df_data.shape[:3]
    if target_space_affine is not None:
        # Create grid in target space
        xyzvox = np.array(np.meshgrid(
            np.arange(target_dimensions[0], dtype=float),
            np.arange(target_dimensions[1], dtype=float),
            np.arange(target_dimensions[2], dtype=float),
            indexing='ij'))

        # Bring points from the grid in target space to x, y,z in target space
        # Applyies the inverse transformation of the warp to brin the points from x, y, z
        # to volxels to later interpolation
        iM = np.linalg.inv(df_affine).dot(target_space_affine)
        t = iM[:3, :3].dot(xyzvox.reshape(3, -1)) + iM[:3, 3, None]
        # Figure out the x, y, z coordinates in the original space
        f = partial(scipy.ndimage.map_coordinates,
                    coordinates=t,
                    output=np.float32, order=1,
                    mode='constant', cval=np.nan)
        voxvals = np.array([f(df_data[..., i]) for i in range(3)])
        # Figure out the voxel coordinates in original space
        iM = np.linalg.inv(im_affine)
        coords = iM[:3, :3].dot(voxvals) + iM[:3, 3, None]
        del t
        del xyzvox
        del voxvals
        gc.collect()

    # If we are using the warp apace
    else:
        iM = np.linalg.inv(im_affine)
        coords = iM[:3, :3].dot(df_data.reshape(-1, 3).T) + iM[:3, 3, None]
        target_space_affine = df_affine

    # If it's a vector and the inverse is defined
    # Rotate the vectors while still in the original space
    if len(im_data.shape) == 4 and im_data.shape[3] == 3:
        if inverse_deformation is not None:
            inv_df_data, inv_df_affine = inverse_deformation
            if len(inv_df_data.shape) > 4:
                inv_df_data = inv_df_data.squeeze()
            im_data_rotated = np.zeros_like(im_data)
            # Spacing in voxel space, for gradient transformations
            spacing = np.linalg.norm(inv_df_affine[:3, :3], axis=1)
            # Transformation from the target voxel space to the transformation voxel
            # space
            t = np.linalg.inv(inv_df_affine).dot(im_affine)
            # Calculate gradient and multiply with vectors
            for i in range(3):
                grad = np.gradient(inv_df_data[..., i], *spacing)
                grad = np.nan_to_num(grad, copy=False)
                for j in range(3):
                    grad_t = scipy.ndimage.affine_transform(
                        grad[j], t[:3, :3], t[:3, 3],
                        output_shape=im_data.shape[:3], order=1,
                        cval=0.0, mode='constant')
                    im_data_rotated[..., i] += grad_t * im_data[..., j]
                    del grad_t
                del grad
                gc.collect()

            # The gradiends were calculated in the voxel space, so we still need to
            # rotate them to the xyz space
            im_data_rotated = np.einsum(
                'ij,klmj->klmi', inv_df_affine[:3, :3] / spacing[:, None], im_data_rotated)

            if keep_vector_length:
                im_data_rotated = im_data_rotated * \
                    (np.linalg.norm(im_data, axis=3))[:, :, :, None] /\
                    (1e-32 +
                     np.linalg.norm(im_data_rotated, axis=3))[:, :, :, None]
            im_data = im_data_rotated

    # Interpolate to the target space
    outdata = _interpolate(im_data, coords, intorder, target_dimensions, cpus)
    del im_data
    gc.collect()
    return outdata


def volumetric_affine(image, affine, target_space_affine,
                      target_dimensions, intorder=1, cpus=1,
                      keep_vector_length=True):
    ''' Applies an affine transformation

    Parameters
    --------
    image: tuple
        tuple (ndarray, affine) with the data to be transformed
    affine: ndarray
        4x4 matrix specfying the transformation from (x, y, z) in the target space to
        (x, y, z) in the original space
    target_space_affine: ndarray
        4x4 matrix with the voxel indices to (x, y, z) coordinates in the target space
    target_dimensions: array
        Dimensions of volume in the final image
    intorder: int (Optional)
        Order of interpolation. See the documentation for the
        scipy.ndimage.map_coordinates function. Default: 1
    keep_vector_length: bool
        Whether to keep the length of the vectors unchanged. Only comes into effect if the
        image is 3 dimensional. Default: True
    Returns
    ------
    outdata: ndarray
        ndarray corresponding to the transormed image
    '''
    # read volumetric data
    im_data, im_affine = image

    # Create the grid in the target space
    iM = np.linalg.inv(im_affine).dot(affine.dot(target_space_affine))
    if im_data.ndim == 4:
        outdata = np.stack(
            [scipy.ndimage.affine_transform(
                im_data[..., i], iM[:3, :3], iM[:3, 3], target_dimensions,
                order=intorder) for i in range(im_data.shape[3])], axis=-1)
    else:
        outdata = scipy.ndimage.affine_transform(
            im_data, iM[:3, :3], iM[:3, 3], target_dimensions,
            order=intorder)

    # Rotate vectors
    if len(outdata.shape) == 4 and outdata.shape[3] == 3:
        outdata_rotated = np.einsum(
            'ij,klmj->klmi', np.linalg.inv(affine[:3, :3]), outdata)
        if keep_vector_length:
            outdata_rotated = outdata_rotated * \
                (np.linalg.norm(outdata, axis=3))[:, :, :, None] /\
                (1e-32 +
                 np.linalg.norm(outdata_rotated, axis=3))[:, :, :, None]

        outdata = outdata_rotated

    return outdata


def _interpolate(im_data, coords, intorder, outdim, cpus=1):
    if len(im_data.shape) > 3:
        indim = im_data.shape[3]
        squeeze = False
    else:
        indim = 1
        im_data = im_data[..., None]
        squeeze = True
    f = partial(
        scipy.ndimage.map_coordinates, coordinates=coords,
        output=im_data.dtype, order=intorder, mode='constant',
        cval=0.0)
    outdata = np.array([f(im_data[..., i]) for i in range(indim)]).T
    outdata = outdata.reshape(tuple(outdim)+(indim,))
    if squeeze:
        outdata = outdata.squeeze()
    return outdata


def nifti_transform(image, warp, ref, out=None, mask=None, order=1, inverse_warp=None,
                    binary=False):
    ''' Transforms a nifti file to the reference space usign the defined warp

    Parameters
    --------
    image: str or tuple
        Name of original file or tuple with (image, affine) pair
    warp: str
        Name of file with the transformation. Can either be a nifti file
        where each voxel corresponds to coordinates (x, y, z) in the original space
        or an affine transformation defined from the target space to the original
        space. In the later case, the name must finish in ".mat", and it will be read
        with the numpy loadtxt file
    ref: str
        Name of reference file. The output will be in the same space
        as the reference (same affine transfomation and same dimensions)
    out: str (optional)
        If not None, the result will be written to this file as a nifti
    mask: str (optional)
        Mak to be applied before transformation
    order: int (optional)
        Interpolation order to be used. Default: 1
    inverse_warp: str
        Name of nifti file with inverse the transformation. Used to rotate vectors to the
        target space in the case of non-liner transformations. If the transformation is
        linear, the inverse matrix is used.
    Returns
    ------
    img: nibabel.Nifti1Pair
        Nibabel image object with tranformed field
    '''
    if isinstance(image, str):
        image = nib.load(image)
        im_data = np.asanyarray(image.dataobj)
        im_affine = image.affine.copy()
        im_data = np.nan_to_num(im_data, copy=False)
    else:
        im_data, im_affine = image
        im_data = np.nan_to_num(im_data, copy=False)
    reference_nifti = nib.load(ref)
    if mask is not None:
        mask = nib.load(mask)
        mask_data = np.asanyarray(mask.dataobj)
        mask_data[np.isnan(mask_data)] = 0
        if not np.allclose(mask.affine, im_affine):
            raise ValueError('Could not apply mask: the affine transformation must be '
                             'the same as the one in the input image')
        if not mask_data.shape[:2] == im_data.shape[:2]:
            raise ValueError('Could not apply mask: The shape of the mask '
                             'must be the same as the image shape')
        im_data *= mask_data.astype(im_data.dtype)

    if isinstance(warp, np.ndarray):
        affdeform = warp
        image = volumetric_affine(
            (im_data, im_affine),
            affdeform,
            target_space_affine=reference_nifti.affine,
            target_dimensions=reference_nifti.header['dim'][1:4],
            intorder=order)

    elif os.path.splitext(warp)[1] in ['.mat', '.txt']:
        affdeform = np.loadtxt(warp)
        image = volumetric_affine(
            (im_data, im_affine),
            affdeform,
            target_space_affine=reference_nifti.affine,
            target_dimensions=reference_nifti.header['dim'][1:4],
            intorder=order)

    else:
        warp_nifti = nib.load(warp)
        if inverse_warp is not None:
            inverse_warp = nib.load(inverse_warp)
            inverse_warp = (np.asanyarray(inverse_warp.dataobj),
                            inverse_warp.affine)
        image = volumetric_nonlinear(
            (im_data, im_affine),
            (np.asanyarray(warp_nifti.dataobj), warp_nifti.affine),
            target_space_affine=reference_nifti.affine,
            target_dimensions=reference_nifti.header['dim'][1:4],
            intorder=order,
            inverse_deformation=inverse_warp)

    if binary:
        image = (image >= .5).astype(np.uint8)
    img = nib.Nifti1Image(image, reference_nifti.affine)
    img.header.set_xyzt_units(reference_nifti.header.get_xyzt_units()[0])
    img.header.set_qform(reference_nifti.header.get_qform(), code=2)
    img.update_header()

    if image.dtype == np.float64:
        img.set_data_dtype(np.float32)

    if out is not None:
        nib.save(img, out)

    return img


def get_names_from_folder_structure(m2m_folder):
    ''' Gets the names of the transformations and reference file from the folder
    structure

    Parameters
    ---------
    m2m_folder: str
        Path to the m2m_{subject_id} folder, generated during the segmantation

    Returns
    -------
    names: dict
        Dictionary with path to transformations and the reference file
    '''
    if not os.path.isfile(templates.mni_volume):
        warnings.warn(
            'Could not find reference volume in mni space at: {0}'.format(
                templates.mni_volume))

    sub_files = SubjectFiles(subpath=m2m_folder)
    if not os.path.isdir(sub_files.subpath):
        raise IOError(
            'The given m2m folder name does not correspond to a directory')

    ref_conf = sub_files.reference_volume
    if not os.path.isfile(ref_conf):
        warnings.warn(
            'Could not find reference volume in subject space at: {0}'.format(
                ref_conf))

    mesh = sub_files.fnamehead
    if not os.path.isfile(mesh):
        warnings.warn('Could not find head mesh at: {0}'.format(mesh))

    mni2conf_nonl = sub_files.mni2conf_nonl
    if not os.path.isfile(mni2conf_nonl):
        warnings.warn('Could not find MNI2Conform non-linear transform at: {0}'.format(
            mni2conf_nonl))

    conf2mni_nonl = sub_files.conf2mni_nonl
    if not os.path.isfile(conf2mni_nonl):
        warnings.warn('Could not find Conform2MNI non-linear transform at: {0}'.format(
            conf2mni_nonl))

    mni2conf_6dof = sub_files.mni2conf_6dof
    mni2conf_12dof = sub_files.mni2conf_12dof

    names = {
        'mesh': mesh,
        'reference_mni': templates.mni_volume,
        'reference_conf': ref_conf,
        'mni2conf_nonl': mni2conf_nonl,
        'conf2mni_nonl': conf2mni_nonl,
        'mni2conf_6dof': mni2conf_6dof,
        'mni2conf_12dof': mni2conf_12dof}
    return names


def warp_volume(image_fn, m2m_folder, out_name,
                transformation_direction='subject2mni',
                transformation_type='nonl',
                reference=None, mask=None,
                labels=None,
                out_original=None,
                order=1,
                method='linear',
                continuous=False,
                binary=False,
                keep_tissues=None):
    ''' Warps a nifti image or a mesh using a linear or non-linar transform, writes out
    the output as a nifti file

    Parameters
    --------
    image_fn: str
        path to mesh or volume. If it is a mesh, all the ElementData fields will be
        warped
    m2m_folder: str
        Path to the m2m_{subject_id} folder, generated during the segmantation
    out_name: str
        Name of output nifti file. If the input is a mesh, the names will be appended with the field
        names.
    transfomation_direction: {'subject2mni' or 'mni2subject'}
        Direction of the tranformation. If 'subject2mni', assumes that the input is in
        subject space, and if 'mni2subject', assumes the input is in MNI space
    transformation_type: {'nonlinear', '6dof', '12dof'}
        Type of tranformation
    reference: str (optional)
        Path to reference nifti. Default: look it up
    mask: str (optional)
        Path to a mask, only applied if input is a nifti volume
    labels: list (optional)
        List of volume labels that should be warped. Only valid if the input is a mesh
    out_original: str (optional)
        If not None, the volume in the original grid be written to a nifti. Only valid
        for mesh inputs
    order: int (optional)
        Interpolation order
    method: {'assign' or 'linear'} (Optional)
        Method for gridding the data. If 'assign', gives to each voxel the value of the element that contains
        it. If linear, first assign fields to nodes, and then perform
        baricentric interpolatiom. Default: linear
    continuous: bool (option)
        Wether fields is continuous across tissue boundaries. Default: False
    keep_tissues: list of tissue tags (Optional)
        Only the fields for the listed tissues are interpolated, rest is set to
        zero. Only applied to mesh inputs. Default: None (all tissues are kept)
    '''
    from ..mesh_tools.mesh_io import read_msh
    names = get_names_from_folder_structure(m2m_folder)

    if transformation_direction not in ['subject2mni', 'mni2subject']:
        raise ValueError('Invalid transformation direction')

    if transformation_type not in ['nonl', '6dof', '12dof']:
        raise ValueError('Invalid transformation type')

    # This is right. Double checked
    if transformation_direction == 'subject2mni':
        if reference is None:
            reference = names['reference_mni']

        if transformation_type == 'nonl':
            warp = names['mni2conf_nonl']
            inverse_warp = names['conf2mni_nonl']

        elif transformation_type == '6dof':
            warp = names['mni2conf_6dof']
            inverse_warp = None

        elif transformation_type == '12dof':
            warp = names['mni2conf_12dof']
            inverse_warp = None

    if transformation_direction == 'mni2subject':
        if reference is None:
            reference = names['reference_conf']

        if transformation_type == 'nonl':
            warp = names['conf2mni_nonl']
            inverse_warp = names['mni2conf_nonl']

        elif transformation_type == '6dof':
            warp = np.linalg.inv(np.loadtxt(names['mni2conf_6dof']))
            inverse_warp = None

        elif transformation_type == '12dof':
            warp = np.linalg.inv(np.loadtxt(names['mni2conf_12dof']))
            inverse_warp = None

    def append_name(fn, append):
        if fn.endswith('.nii.gz'):
            name, end = (fn[:-len('.nii.gz')], '.nii.gz')
        else:
            name, end = os.path.splitext(fn)
        return name + '_' + append + end

    # Ugly hack to get MNI at the end of file names
    if transformation_direction == 'subject2mni':
        if out_name.endswith('.nii.gz'):
            name, end = (out_name[:-len('.nii.gz')], '.nii.gz')
        else:
            name, end = os.path.splitext(out_name)
        if end == '':
            end = '.nii.gz'
        out_name = name + '_MNI' + end

    if os.path.splitext(image_fn)[1] == '.msh':
        m = read_msh(image_fn)
        if keep_tissues is not None:
            m = m.crop_mesh(tags=keep_tissues)

        logger.info('Warping mesh: {0}'.format(image_fn))
        for ed in m.elmdata + m.nodedata:
            name = append_name(out_name, ed.field_name)
            if out_original is not None:
                name_original = append_name(out_original, ed.field_name)
            else:
                name_original = None
            logger.info('Warping field: {0}'.format(ed.field_name))
            logger.info('To file: {0}'.format(name))
            logger.debug('Transformation type: {0}'.format(
                transformation_type))
            logger.debug('Method: {0}'.format(method))
            logger.debug('Labels: {0}'.format(labels))
            logger.debug('Continuous: {0}'.format(continuous))
            ed.to_deformed_grid(warp, reference,
                                out=name,
                                out_original=name_original,
                                tags=labels,
                                order=order,
                                method=method,
                                continuous=continuous,
                                inverse_warp=inverse_warp,
                                reference_original=names['reference_conf'],
                                binary=binary)

    else:
        logger.info('Warping nifti file: {0}'.format(image_fn))
        logger.info('To file: {0}'.format(out_name))
        logger.debug('Transformation type: {0}'.format(transformation_type))
        logger.debug('Mask: {0}'.format(mask))
        nifti_transform(image_fn, warp, reference, out=out_name,
                        mask=mask, order=order, inverse_warp=inverse_warp,
                        binary=binary)


def interpolate_to_volume(fn_mesh, reference, fn_out, create_masks=False,
                          method='linear', continuous=False, create_label=False,
                          keep_tissues=None):
    ''' Interpolates the fields in a mesh and writem them to nifti files

    Parameters:
    -----------
    fn_mesh: str, or simnibs.msh.Msh
        Name of mesh file, or a preloaded mesh
    reference: str
        Path to the m2m_{subject_id} folder, generated during the segmantation, or to a
        reference nifti file
    out_name: str
        Name of output nifti file. If the input is a mesh, the names will be appended with the field
        names.
    create_masks: bool
        Mask mode: write tissue masks instead of fields
    method: {'assign' or 'linear'} (Optional)
        Method for gridding the data. If 'assign', gives to each voxel the value of the element that contains
        it. If linear, first assign fields to nodes, and then perform
        baricentric interpolatiom. Default: linear
    continuous: bool
        Wether fields is continuous across tissue boundaries. Default: False
    create_label: bool
        Label mode: write label image from mesh instead of fields
    keep_tissues: list of tissue tags (Optional)
        Only the fields for the listed tissues are interpolated, rest is set to
        zero. Default: None (all tissues are kept)
    '''
    from ..mesh_tools.mesh_io import read_msh, ElementData
    if os.path.isdir(reference):
        names = get_names_from_folder_structure(reference)
        reference = names['reference_conf']
    if not os.path.isfile(reference):
        raise IOError('Could not find reference file: {0}'.format(reference))
    if isinstance(fn_mesh, str):
        if not os.path.isfile(fn_mesh):
            raise IOError('Could not find mesh file: {0}'.format(fn_mesh))
        mesh = read_msh(fn_mesh)
    else:
        mesh = copy.deepcopy(fn_mesh)

    if keep_tissues is not None:
        mesh = mesh.crop_mesh(tags=keep_tissues)

    image = nib.load(reference)
    affine = image.affine
    n_voxels = image.header['dim'][1:4]
    fn, ext = os.path.splitext(fn_out)
    if ext == '':
        ext = '.nii.gz'
        fn_out += '.nii.gz'

    def append_name(fn, append):
        if fn.endswith('.nii.gz'):
            name, end = (fn[:-len('.nii.gz')], '.nii.gz')
        else:
            name, end = os.path.splitext(fn)
        return name + '_' + append + end

    if create_masks:
        mesh = mesh.crop_mesh(elm_type=4)
        vol_tags = np.unique(mesh.elm.tag1)
        for v in vol_tags:
            name = fn + '_mask_{0}'.format(v) + ext
            field = np.zeros(mesh.elm.nr, dtype=np.uint8)
            field[mesh.elm.tag1 == v] = 1
            ed = ElementData(field)
            ed.mesh = mesh
            ed.to_nifti(n_voxels, affine, fn=name, qform=image.header.get_qform(),
                        method='assign')
    elif create_label:
        mesh = mesh.crop_mesh(elm_type=4)
        field = np.zeros(mesh.elm.nr, dtype=np.int32)
        field = mesh.elm.tag1
        ed = ElementData(field)
        ed.mesh = mesh
        ed.to_nifti(n_voxels, affine, fn=fn_out, qform=image.header.get_qform(),
                    method='assign')
    else:
        if len(mesh.elmdata) + len(mesh.nodedata) == 0:
            warnings.warn('No fields found in mesh!')
        for ed in mesh.elmdata + mesh.nodedata:
            name = append_name(fn_out, ed.field_name)
            logger.info('Field: {0}'.format(ed.field_name))
            logger.info('To file: {0}'.format(name))
            logger.debug('Method: {0}'.format(method))
            logger.debug('Continuous: {0}'.format(continuous))
            ed.to_nifti(n_voxels, affine, fn=name, method=method,
                        continuous=continuous)
            gc.collect()


def coordinates_nonlinear(coordinates, deformation, intorder=1, vectors=None,
                          keep_length=True):
    ''' Applies a volumetric non-linear transformation

    Parameters
    --------
    coordinates: nx3 ndarray
        ndarray with coordinates to be transfored
    deformation: tuple
        tuple (ndarray, affine) with defomation field. The deformation field specifies
        in each voxel the coordinates (x, y, z) of the target space. This is the opposite
        transformation from the one used in the volumetric transform
    intorder: int
        Order of interpolation. See the documentation for the
        scipy.ndimage.map_coordinates function
    vectors: nx3 ndattay (optional)
        ndarray with vectors to be transformed. Needs to have the same lenght as
        "coordinates".
    keep_length: bool (optional)
        Wether to preseve the lengths of the vectors. Default: True

    Returns
    ------
    coordinates_transformed: nx3 ndarray
        ndarray corresponding to the transormed coordinates. If outside the volume, will
        return the transformation of the closest coordinate inside the volume
    '''
    if len(coordinates.shape) == 1:
        coordinates = coordinates[None, :]
    if coordinates.shape[1] != 3:
        raise ValueError('Input coordinates should have a nx3 format')
    if vectors is not None and len(vectors.shape) == 1:
        vectors = vectors[None, :]

    df_data, df_affine = deformation
    df_data = df_data.squeeze()

    iM = np.linalg.inv(df_affine)
    t = iM[:3, :3].dot(coordinates.T) + iM[:3, 3, None]
    # Figure out the x, y, z coordinates in the original space
    f = partial(scipy.ndimage.map_coordinates,
                coordinates=t,
                output=np.float32, order=intorder,
                mode='nearest')
    coords = np.array([f(df_data[..., i]) for i in range(3)])

    if vectors is not None:
        original_l = np.linalg.norm(vectors, axis=1)
        v = coordinates + (vectors / original_l[:, None]) * .1
        t = iM[:3, :3].dot(v.T) + iM[:3, 3, None]
        f = partial(scipy.ndimage.map_coordinates,
                    coordinates=t,
                    output=np.float32, order=intorder,
                    mode='nearest')
        v_trafo = np.array([f(df_data[..., i]) for i in range(3)])
        vec_transf = (v_trafo - coords) / .1
        if keep_length:
            vec_transf /= np.linalg.norm(vec_transf, axis=0)[None, :]
        vec_transf *= original_l[None, :]

        return coords.T, vec_transf.T
    return coords.T


def coordinates_affine(coordinates, affine):
    ''' Applies an affine transformation

    Parameters
    --------
    coordinates: nx3 ndarray
        coordinates to be transformed
    affine: ndarray
        4x4 matrix specfying the transformation from (x, y, z) in the original space to
        (x, y, z) in the target space. This is the inverse of the matrix used to
        transform volumes
    Returns
    ------
    coordinates_transformed: nx3 ndarray
        ndarray corresponding to the transormed coordinates
    '''
    # read volumetric data
    if len(coordinates.shape) == 1:
        coordinates = coordinates[None, :]
    if coordinates.shape[1] != 3:
        raise ValueError('Input coordinates should have a nx3 format')

    M = affine
    coords = M[:3, :3].dot(coordinates.T) + M[:3, 3, None]
    return coords.T


def vectors_affine(vectors, affine, keep_length=True):
    ''' Applies an affine transformation to vector

    Parameters
    --------
    vectors: nx3 ndarray
        vectors to be transformed
    affine: ndarray
        4x4 matrix specfying the transformation from (x, y, z) in the original space to
        (x, y, z) in the target space. This is the inverse of the matrix used to
        transform volumes
    keep_length: Bool
        Wether to keep the length of the vectors
    Returns
    ------
    coordinates_transformed: nx3 ndarray
        ndarray corresponding to the transormed vectors
    '''
    # read volumetric data
    if len(vectors.shape) == 1:
        vectors = vectors[None, :]
    if vectors.shape[1] != 3:
        raise ValueError('Input coordinates should have a nx3 format')

    M = affine
    vectors_tr = M[:3, :3].dot(vectors.T)
    vectors_tr = vectors_tr.T
    if keep_length:
        original_l = np.linalg.norm(vectors, axis=1)
        vectors_tr /= np.linalg.norm(vectors_tr, axis=1)[:, None]
        vectors_tr *= original_l[:, None]
    return vectors_tr


def project_points_on_surface(mesh, pts, surface_tags = None, distance = 0.):
    ''' Find the position on the surface closest to each coordinate

     Parameters
     -------
     coords: nx3 ndarray
         Vectors to be transformed
     mesh: simnibs.msh.mesh_io.Msh
         Mesh structure
     surface_tags: int (optional)
             Tag of the target surface.
             Default: None (positions will be projected on closest surface)
     distance: float or nx1 ndarray (optional)
         Distance (normal) to the surface to be enforced. Default: 0
         Note: negative values will move the point inside, positive values
               outside the volume defined by the surface

     Returns
     ------
     coords_projected: nx3 ndarray
         coordinates projected on the surface
    '''
    pts = np.array(pts)
    if pts.ndim == 1:
        pts = pts.reshape((1,len(pts)))

    # get surface
    if surface_tags is not None:
        tr_of_interest = (mesh.elm.elm_type == 2) * (np.in1d(mesh.elm.tag1, surface_tags))
    else:
        tr_of_interest = mesh.elm.elm_type == 2
    tri_node_list = mesh.elm.node_number_list[tr_of_interest, :3] - 1
    tri_nodes = np.unique(tri_node_list)

    old2new = np.zeros(tri_nodes[-1] + 1, dtype = 'int32')
    old2new[tri_nodes] = np.arange(len(tri_nodes))
    surf = {
        'tris': old2new[tri_node_list],
        'points': mesh.nodes.node_coord[tri_nodes]
    }

    # get indices of close-by surface nodes and their connected triangles
    pttris = _get_nearest_triangles_on_surface(pts, surf, n = 3)

    # project points on triangles
    tris, _, projs, dists = _project_points_to_surface(pts, surf, pttris)

    # ensure distance (optional)
    if not np.all(np.isclose(distance, 0)):
        distance = np.array(distance)
        if distance.ndim == 0:
            distance = distance.reshape(1)

        tri_pts = surf['points'][surf['tris'][tris]]
        sideA = tri_pts[:, 1] - tri_pts[:, 0]
        sideB = tri_pts[:, 2] - tri_pts[:, 0]
        n = np.cross(sideA, sideB)
        n /= np.linalg.norm(n, axis=1)[:, None]

        projs += distance[:, None]*n

    return projs


def transform_tdcs_positions(coords, transf_type, transf_def, pos_y=None, mesh=None):
    ''' Transform positions of tDCS electrodes

    Parameters
    ------
    coords: nx3 ndarray
        Array with electrode centers
    transf_type: {'affine', 'nonl'}
        Type of transformation
    transf_def:
        if transf_def == 'affine':
            4x4 matrix defining the transform
        if transf_type == 'nonl':
            tuple (ndarray, affine) with defomation field. The deformation field specifies
            in each voxel the coordinates (x, y, z) of the target space. This is the opposite
            transformation from the one used in the volumetric transform
    pos_y: nx3 ndarray (optional)
        Array with reference points
    mesh: simnibs.msh.mesh_io.Msh
        Mesh structure
    Returns
    ------
    coords_transf: nx3 ndarray
        Array with electrodes transformed and projected to the scalp
    pos_y_transf: nx3 ndarray
        Array with reference points transformed and projected to the scalp
    '''
    if len(coords.shape) == 1:
        coords = coords[None, :]
    if coords.shape[1] != 3:
        raise ValueError('Input coordinates should have a nx3 format')

    if transf_type == 'affine':
        coords_transf = coordinates_affine(coords, transf_def)
        if pos_y is not None:
            pos_y_transf = coordinates_affine(pos_y, transf_def)

    elif transf_type == 'nonl':
        coords_transf = coordinates_nonlinear(coords, transf_def)
        if pos_y is not None:
            pos_y_transf = coordinates_nonlinear(pos_y, transf_def)

    else:
        raise ValueError('Invalid transformation type')

    if mesh:
        coords_transf = project_points_on_surface(mesh, coords_transf, surface_tags = 1005)
        if pos_y is not None:
            pos_y_transf = project_points_on_surface(mesh, pos_y_transf, surface_tags = 1005)

    if pos_y is not None:
        return coords_transf, pos_y_transf
    else:
        return coords_transf


def transform_tms_positions(coords, v_y, v_z, transf_type, transf_def,
                            mesh=None, distances=None):
    ''' Transform positions and vectors for TMS positions
    Enforces the orthonomality of v_z and v_y by orthogonalizing the transformed v_y
    w.r.t v_z

    Parameters
    ------
    coords: nx3 ndarray
        Array with coil centers
    v_y: nx3 ndarray
        Array with the "y" axis of the coil coordinate system (allong handle)
    v_z: nx3 ndarray
        Array with the "z" axis of the coil coordinate system (top - bottom)
    transf_type: {'affine', 'nonl'}
        Type of transformation
    transf_def:
        if transf_def == 'affine':
            4x4 matrix defining the transform
        if transf_type == 'nonl':
            tuple (ndarray, affine) with defomation field. The deformation field specifies
            in each voxel the coordinates (x, y, z) of the target space. This is the opposite
            transformation from the one used in the volumetric transform
    mesh: simnibs.msh.mesh_io.Msh
        Mesh structure. Used to project and get the right distance to skin. Default: do not project
    distances: nx1 ndarray
        Array with coil distances
    Returns
    ------
    coords_t: nx3 ndarray
        Array with transformed coil centers
    v_y_t: nx3 ndarray
        Array with the transformed "y" axis of the coil coordinate system (allong handle)
    v_z_z: nx3 ndarray
        Array with the transformd "z" axis of the coil coordinate system (top - bottom)
    '''
    if len(coords.shape) == 1:
        coords = coords[None, :]
    if coords.shape[1] != 3:
        raise ValueError('Input coordinates should have a nx3 format')

    if transf_type == 'affine':
        coords_transf = coordinates_affine(coords, transf_def)
        if mesh:
            coords_transf = project_points_on_surface(mesh, coords_transf,
                                                      surface_tags = 1005,
                                                      distance=distances)
        vy_transf = vectors_affine(v_y, transf_def)
        vz_transf = vectors_affine(v_z, transf_def)

    elif transf_type == 'nonl':
        coords_transf, vy_transf = coordinates_nonlinear(coords, transf_def,
                                                         vectors=v_y)
        coords_transf, vz_transf = coordinates_nonlinear(coords, transf_def,
                                                         vectors=v_z)
        if mesh:
            coords_transf = project_points_on_surface(mesh, coords_transf,
                                                      surface_tags = 1005,
                                                      distance=distances)
    else:
        raise ValueError('Invalid transformation type')

    vz_transf /= np.linalg.norm(vz_transf, axis=1)[:, None]
    vy_transf = vy_transf - vz_transf * \
        np.sum(vz_transf * vy_transf, axis=1)[:, None]
    vy_transf /= np.linalg.norm(vy_transf, axis=1)[:, None]

    return coords_transf, vy_transf, vz_transf


def subject2mni_coords(coordinates, m2m_folder, transformation_type='nonl'):
    ''' Warps a set of coordinates in subject space to MNI space

    Parameters
    ------------
    coordinates: list or numpy array in Nx3 format
        Coordinates to be transformd

    m2m_folder: str
        Path to the m2m_{subject_id} folder, generated during the segmantation

    transformation_type: {'nonl', '6dof', '12dof'}
        Type of tranformation, non-linear, 6 or 12 degrees of freedom

    Returns
    ----------
    transformed_coords: Nx3 numpy array
        Array with transformed coordinates

    '''
    transformed = warp_coordinates(
        coordinates, m2m_folder,
        transformation_direction='subject2mni',
        transformation_type=transformation_type,
        out_name=None,
        out_geo=None)[1]
    if np.array(coordinates).ndim == 1:
        transformed = np.squeeze(transformed)
    return transformed


def mni2subject_coords(coordinates, m2m_folder, transformation_type='nonl'):
    ''' Warps a set of coordinates in MNI space to subject space

    Parameters
    ------------
    coordinates: list or numpy array in Nx3 format
        Coordinates to be transformd

    m2m_folder: str
        Path to the m2m_{subject_id} folder, generated during the segmantation

    transformation_type: {'nonl', '6dof', '12dof'}
        Type of tranformation, non-linear, 6 or 12 degrees of freedom

    Returns
    ----------
    transformed_coords: Nx3 numpy array
        Array with transformed coordinates

    '''
    transformed = warp_coordinates(
        coordinates, m2m_folder,
        transformation_direction='mni2subject',
        transformation_type=transformation_type,
        out_name=None,
        out_geo=None)[1]
    if np.array(coordinates).ndim == 1:
        transformed = np.squeeze(transformed)
    return transformed


def warp_coordinates(coordinates, m2m_folder,
                     transformation_direction='subject2mni',
                     transformation_type='nonl',
                     out_name=None,
                     out_geo=None,
                     mesh_in=None):
    ''' Warps a set of coordinates
    For simpler calls, please see subject2mni_coords and mni2subject_coords

    Parameters
    --------
    coordinates: str, list or ndarray
        if list or ndarray (Nx3 format):
            Will do a simple transformation of the
        If path to csv file:
            The CSV file must have at least 4 columns, dependind on the type of data to
            be transformed
            Generic:
                Generic, pos_x, pos_y, pos_z, name, ...
            Positions will not be changed after transformation.

            Fiducial, Electrode, ReferenceElectrode:
                 Type, pos_x, pos_y, pos_z, name, whatever
                 Type, pos_x, pos_y, pos_z, pos2_x, pos2_y, pos2_z, name, ...
            if the direction is mni2subject: Positions will be projected on skin after
            transformation. Type must be Fiducial, Electrode, or ReferenceElectrode.

            CoilPos:
                Type, pos_x, pos_y, pos_z, ez_x, ez_y, ez_z, ey_x, ey_y, ey_z, dist, name, ...
            if the direction is mni2subject: position will be adjusted after transformation
            to have specified distance to skin (mni2subject_coords ignores distances)

            You can also input a list with the same structure as the csv files

    m2m_folder: str
        Path to the m2m_{subject_id} folder, generated during the segmantation

    transfomation_direction: {'subject2mni' or 'mni2subject'}
        Direction of the tranformation. If 'subject2mni', assumes that the input is in
        subject space, and if 'mni2subject', assumes the input is in MNI space

    out_name: str
        Name of output file. If defined, the function will write to this file as a csv

    transformation_type: {'nonlinear', '6dof', '12dof'}
        Type of tranformation

    out_geo: str
        Writes out a geo file for visualization. Only works when out_name is also set

    mesh_in : mesh file, optional
        scalp or head mesh, used to project electrodes onto scalp or ensure
        TMS coil distances from scalp; is used only for direction 'mni2subject'
        (standard: None; the head mesh file will then be automatically loaded,
         which is slower for repeated applications)

    Returns
    ----------
    type: list
        List with point types. Can be 'Generic', 'Fiducial', 'Electrode',
        'ReferenceElectrode', 'CoilPos'.
    transformed_coords: nx3 ndarray
        Array with transformed coordinates
    transformed_extra: list
        Transformed extra information for each position (pos_y_dir) for Electrode and
        ReferenceElectrode, coil axis and distance for CoilPos
    name: list
        Names of positions
    '''
    from ..mesh_tools.mesh_io import read_msh
    names = get_names_from_folder_structure(m2m_folder)
    # Read CSV
    if isinstance(coordinates, str):
        type_, coordinates, extra, name, extra_cols, header = read_csv_positions(
            coordinates)
    else:
        try:
            type_, coordinates, extra, name, extra_cols, header = coordinates
        except ValueError:
            coordinates = np.asarray(coordinates, dtype=float)
            if len(coordinates.shape) == 1:
                coordinates = coordinates[None, :]
            if coordinates.shape[1] != 3:
                raise ValueError('Input coordinates should have a nx3 format')
            type_ = ['Generic'] * len(coordinates)
            header = []
            extra_cols = [None] * len(coordinates)
            extra = [None] * len(coordinates)
            name = [None] * len(coordinates)

    transformed_coords = np.zeros_like(coordinates)
    transformed_extra = extra

    if transformation_direction not in ['subject2mni', 'mni2subject']:
        raise ValueError('Invalid transformation direction')

    if transformation_type not in ['nonl', '6dof', '12dof']:
        raise ValueError('Invalid transformation type')

    # Here we use the inverse transforms from the ones in the volumetric transforms
    if transformation_direction == 'subject2mni':
        if transformation_type == 'nonl':
            warp = names['conf2mni_nonl']
        elif transformation_type == '6dof':
            warp = np.linalg.inv(np.loadtxt(names['mni2conf_6dof']))
        elif transformation_type == '12dof':
            warp = np.linalg.inv(np.loadtxt(names['mni2conf_12dof']))
        mesh = None

    if transformation_direction == 'mni2subject':
        if transformation_type == 'nonl':
            warp = names['mni2conf_nonl']
        elif transformation_type == '6dof':
            warp = np.loadtxt(names['mni2conf_6dof'])
        elif transformation_type == '12dof':
            warp = np.loadtxt(names['mni2conf_12dof'])

    # Apply transformation
    if transformation_type == 'nonl':
        ttype = 'nonl'
        image = nib.load(warp)
        warp = (np.asanyarray(image.dataobj), image.affine)
        simple = coordinates_nonlinear
    else:
        ttype = 'affine'
        simple = coordinates_affine

    # Transforms all generic type rows
    generic = [i for i, t in enumerate(type_) if t == 'Generic']
    if len(generic) > 0:
        transformed_coords[generic, :] = simple(coordinates[generic, :], warp)

    # delay reading the mesh
    if len(generic) != len(type_) and transformation_direction == 'mni2subject':
        if mesh_in is None:
            mesh = read_msh(names['mesh'])
        else:
            mesh = mesh_in

    # Transform all electrode types
    electrode = [i for i, t in enumerate(
        type_) if t in ['Fiducial', 'Electrode', 'ReferenceElectrode']]
    if len(electrode) > 0:
        transformed_coords[electrode, :] = transform_tdcs_positions(
            coordinates[electrode, :], ttype, warp, mesh=mesh)

    # Transform the pos_y coordinates, if any
    electrode_posy = [i for i in electrode if extra[i] is not None]
    posy = np.array([extra[i] for i in electrode_posy])
    if len(posy.shape) == 1:
        posy = posy[None, :]

    if len(electrode_posy) > 0:
        posy = transform_tdcs_positions(posy, ttype, warp, mesh=mesh)
    for el, y in zip(electrode_posy, posy):
        transformed_extra[el] = y

    # Transforms the cois
    coil = [i for i, t in enumerate(type_) if t == 'CoilPos']
    coil_extra = []
    if len(coil) > 0:
        e = np.vstack([extra[i] for i in coil])
        if len(e.shape) == 1:
            e = e[None, :]
        vz = e[:, :3]
        vy = e[:, 3:6]
        d = e[:, 6]
        transformed_coords[coil, :], vy, vz = transform_tms_positions(
            coordinates[coil, :], vy, vz, ttype, warp, mesh=mesh, distances=d)

        coil_extra = np.hstack((vz, vy, d[None, :].T))

    for c, e in zip(coil, coil_extra):
        transformed_extra[c] = e

    for t in type_:
        if t not in ['Fiducial', 'Electrode', 'ReferenceElectrode', 'Generic',
                     'CoilPos']:
            warnings.warn('Unrecogized column type: {0}, '
                          'valid types are:\n'
                          'Fiducial, Electrode, ReferenceElectrode, Generic, '
                          'CoilPos'.format(t))

    # Write CSV
    if out_name is not None:
        write_csv_positions(
            out_name, type_, transformed_coords, name, transformed_extra, extra_cols, header)
        if out_geo is not None:
            csv_to_geo(out_name, out_geo)

    return type_, transformed_coords, transformed_extra, name


def csv_to_geo(fn_csv, fn_out):
    ''' Writes a .geo file based on a .csv file

    Parameters
    -------
    fn_csv: str
       name of csv file
    fn_out: str
        name of output geo file
    '''
    type_, coordinates, extra, name, _, _ = read_csv_positions(fn_csv)
    file_text = 'View"' + '' + '"{\n'

    def write_electrode(file_text, coordinates, extra, name):
        file_text += "SP(" + ", ".join([str(c)
                                        for c in coordinates]) + "){0};\n"
        if extra is not None:
            file_text += "SL(" + ", ".join([str(c) for c in coordinates]) +\
                ", " + ", ".join([str(c) for c in extra]) + "){0, 0};\n"
        if name is not None:
            text_coords = coordinates + [0., 0., 5.]
            file_text += 'T3(' + ', '.join([str(c) for c in text_coords]) +\
                ', 0){"' + name + '"};\n'
        return file_text

    def write_coil(file_text, coordinates, extra, name):
        file_text += "SP(" + ", ".join([str(c)
                                        for c in coordinates]) + "){0};\n"
        ez = coordinates + extra[:3] * 30
        ey = coordinates + extra[3:6] * 10
        file_text += "SL(" + ", ".join([str(c) for c in coordinates]) +\
            ", " + ", ".join([str(c) for c in ez]) + "){0., 0.};\n"
        file_text += "SL(" + ", ".join([str(c) for c in coordinates]) +\
            ", " + ", ".join([str(c) for c in ey]) + "){0., 0.};\n"
        if name is not None:
            text_coords = coordinates + [0., 0., 5.]
            file_text += 'T3(' + ', '.join([str(c) for c in text_coords]) +\
                ', 0){"' + name + '"};\n'
        return file_text

    for t, c, e, n in zip(type_, coordinates, extra, name):
        if t in ['Fiducial', 'Electrode', 'ReferenceElectrode', 'Generic']:
            file_text = write_electrode(file_text, c, e, n)
        if t in ['CoilPos']:
            file_text = write_coil(file_text, c, e, n)

    file_text += "};\n"
    file_text += "myView = PostProcessing.NbViews-1;\n"
    file_text += "View[myView].PointType=1;\n"
    file_text += "View[myView].PointSize=6;\n"
    file_text += "View[myView].LineType=1;\n"
    file_text += "View[myView].LineWidth=2;\n"

    with open(fn_out, 'w') as f:
        f.write(file_text)


def get_vox_size(affine):
    ''' Determines the voxel size based on the affine transformation

    Parameters
    ------------
    affine: 4x4 np.ndarray
        Affine transformation

    Returs
    ----------
    voxsize: 3x1 np.ndarray
        Voxel sizes
    '''
    return np.linalg.norm(affine[:3, :3], axis=0)


def crop_vol(vol, affine, mask, thickness_boundary=0):
    ''' Crops a volume from vol based on the bounding box defined by the mask
    Also fixes the affine transformation based on the bounding box

    Parameters
    ------------
    vol: 3D np.ndarray
        Volume to be cropped

    Affine: 4x4 np.ndarray
        Affine transfromation from the original volume to world space

    mask: 3D np.ndarray of floats
        Volume defining the bouding box

    thickness_boundary: float (optional)
        Thickness (in mm) to add to each side of the bounding box defined by the mask.
        Default: 0

    Returns
    ---------
    cropped_vol: 3D np.ndarray
        Cropped volume

    new_affine: 4x4 np.ndarray
        Affine transformation from the cropped volume to world space

    bound_box: 3x2 np.ndarray
        Bound box used to crop the mesh
    '''
    # Figure out bound box from mask
    non_zero = np.where(mask)
    bound_box = np.array([
        (np.min(nz), np.max(nz))
        for nz in non_zero
    ], dtype=int)

    vx_size = get_vox_size(affine)
    thickness_boundary = (thickness_boundary/vx_size).astype(int)

    # Add thickness boundary in the left/bottom/front and fix
    bound_box[:, 0] -= thickness_boundary
    bound_box[bound_box[:, 0] < 0, 0] = 0

    # Add thickness boundary in the right/top/back and fix
    bound_box[:, 1] += thickness_boundary
    out_of_bound = bound_box[:, 1] >= vol.shape
    bound_box[out_of_bound, 1] = np.array(vol.shape)[out_of_bound] - 1

    cropped_vol = vol[
        bound_box[0, 0]:bound_box[0, 1] + 1,
        bound_box[1, 0]:bound_box[1, 1] + 1,
        bound_box[2, 0]:bound_box[2, 1] + 1
    ]

    new_affine = affine.copy()
    new_affine[:3, 3] += affine[:3, :3].dot(bound_box[:, 0])
    return cropped_vol, new_affine, bound_box

def pad_vol(vol, affine, voxel_pad):
    ''' Pads a volume from vol based on the voxel_pad
        Also fixes the affine transformation

        Parameters
        ------------
        vol: 3D np.ndarray
            Volume to be padded

        Affine: 4x4 np.ndarray
            Affine transfromation from the original volume to world space

        voxel_pad: int
            The number of voxels the volume will be padded with all around

        Returns
        ---------
        padded_vol: 3D np.ndarray
            Padded volume

        new_affine: 4x4 np.ndarray
            Affine transformation from the padded volume to world space
        '''
    padded_shape = (
        vol.shape[0] + 2 * voxel_pad,
        vol.shape[1] + 2 * voxel_pad,
        vol.shape[2] + 2 * voxel_pad
    )
    padded_vol = np.zeros(padded_shape, dtype=vol.dtype)
    padded_vol[
        voxel_pad:voxel_pad + vol.shape[0],
        voxel_pad:voxel_pad + vol.shape[1],
        voxel_pad:voxel_pad + vol.shape[2]
    ] = vol

    new_affine = np.copy(affine)
    new_affine[:3, 3] = affine[:3, :3] @ [-voxel_pad, -voxel_pad, -voxel_pad] + affine[:3, 3]

    return padded_vol, new_affine

def resample_vol(vol, affine, target_res, order=1, mode='nearest'):
    ''' Change the resolution of the volume

    Parameters
    -------------
    vol: 3D np.ndarray
        Volume to have the resolution changed

    affine: 4x4 np.ndarray
        Affine tranformation from the original volume to world space

    target_res: 3x1 ndarray
        Target resolution

    order: int (optional)
        Interpolation order. Default: 1
        see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.affine_transform.html#scipy.ndimage.affine_transform

    mode: str (optional)
        How to handle boundaries. Default: 'nearest'
        see https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.affine_transform.html#scipy.ndimage.affine_transform


    Returns
    ------------
    resampled: 3D np.ndarray
        Resampled volume

    new_affine: 4x4 np.ndarray
        Affine tranformation from the resampled volume to world space

    original_res: 3x1 ndarray
        Resolution of the ORIGINAL image (before resampling)
    '''
    target_res = np.squeeze(target_res)
    if target_res.ndim == 0:
        target_res = target_res*np.ones(3)
    original_res = get_vox_size(affine)

    transform = np.squeeze(target_res/original_res)
    new_affine = affine.astype(float)
    new_affine[:3, :3] *= transform
    # We need to change the translation component of the affine
    # to make sure the voxel centers match
    # We just solve
    # R*(-0.5, -0.5, -0.5) + t  = R'*(-0.5, -0.5, -0.5) + t'
    corner = -0.5*np.ones(3)
    offset = (
        affine[:3, :3].dot(corner) -
        new_affine[:3, :3].dot(corner)
    )
    new_affine[:3, 3] += offset
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        resampled = scipy.ndimage.affine_transform(
            vol,
            transform,
            # the voxel coordinates define the center of the voxels
            offset=(transform - 1)/2.,
            output_shape=(vol.shape/transform).astype(int),
            order=order,
            mode=mode
        )

    return resampled, new_affine, original_res


def normalize(arr: np.ndarray, axis=None, inplace: bool = False):
    """Normalize along a particular axis (or axes) of `v` avoiding
    RuntimeWarning due to division by zero.

    PARAMETERS
    ----------
    v : ndarray
        Array to nomalize.
    axis:
        The axis to normalize along, e.g., `axis=1` will normalize rows).
        Like axis argument of numpy.linalg.norm.

    RETURNS
    -------
    v (normalized)
    """
    size = np.linalg.norm(arr, axis=axis, keepdims=True)
    if inplace:
        np.divide(arr, size, where=size != 0, out=arr)
    else:
        return np.divide(arr, size, where=size != 0)


class SurfaceMorph:
    """
    Computes a mapping between surfaces defined on a sphere.
    """

    def __init__(self, surf_from, surf_to, method="linear", n: int = 2):
        assert method in {"nearest", "linear"}
        self.method = method
        self.n = n

        self._fit(surf_from, surf_to)

    def _fit(self, surf_from, surf_to):
        """Create a morph map which allows morphing of values from the nodes in
        `surf_from` to `surf_to` by nearest neighbor or linear interpolation.

        A morph map is a sparse matrix with
        dimensions (n_points_surf_to, n_points_surf_from) where each row has
        exactly three entries that sum to one. It is created by projecting each
        point in `surf_to` onto closest triangle in `surf_from` and determining
        the barycentric coordinates.

        Testing all points against all triangles is expensive and inefficient,
        thus we compute an approximation by finding, for each point in
        `surf_to`, the `self.n` nearest nodes on `surf_from` and the triangles
        to which these points belong. We then test only against these triangles.

        PARAMETERS
        ----------
        surf_from :
            The source mesh (i.e., the mesh to interpolate *from*).
        surf_to :
            The target mesh (i.e., the mesh to interpolate *to*).
        """
        # Ensure on unit sphere
        points_from = normalize(surf_from.nodes.node_coord, axis=1)
        points_to = normalize(surf_to.nodes.node_coord, axis=1)
        n_from = len(points_from)
        n_to = len(points_to)

        if self.method == "nearest":
            # in_v = np.copy(in_surf.nodes.node_coord)
            # in_v /= np.average(np.linalg.norm(in_v, axis=1))
            # kdtree = scipy.spatial.cKDTree(in_v)

            # out_v = np.copy(out_surf.nodes.node_coord)
            # # Normalize the radius of the output sphere
            # out_v /= np.average(np.linalg.norm(out_v, axis=1))
            # _, closest = kdtree.query(out_v)
            kdtree = scipy.spatial.cKDTree(points_from)
            # self.morph_mat = kdtree.query(points_to)[1]
            rows = np.arange(n_to)
            cols = kdtree.query(points_to)[1]
            weights = np.ones(n_to)

        elif self.method == "linear":
            # Find the triangle (in surf) to which each point in points
            # projects and get the associated weights
            # n_points, d_points = points_to.shape
            surf_ = dict(
                points=points_from, tris=surf_from.elm.node_number_list[:, :3] - 1
            )
            pttris = _get_nearest_triangles_on_surface(points_to, surf_, self.n)
            tris, weights, _, _ = _project_points_to_surface(points_to, surf_, pttris)
            rows = np.repeat(np.arange(n_to), points_to.shape[1])
            cols = surf_["tris"][tris].ravel()
            weights = weights.ravel()
        else:
            raise ValueError

        self.morph_mat = csr_matrix((weights, (rows, cols)), shape=(n_to, n_from))

    def transform(self, values):
        return self.morph_mat @ values


def make_cross_subject_morph(
    subject_from: Union[Path, str, SubjectFiles],
    subject_to: Union[Path, str, SubjectFiles],
    subsampling_from: Union[None, int] = None,
    subsampling_to: Union[None, int] = None,
    surface_morph_kwargs: Union[dict, None] = None,
):
    """
    subject_from, subject_to :
        Special subject name is 'fsaverage'.
    subsampling_from, subsampling_to :
        Special subsampling for 'fsaverage' is 10, 40 and None.

    """
    from simnibs.mesh_tools.mesh_io import (
        load_reference_surfaces,
        load_subject_surfaces,
    )

    surface_morph_kwargs = surface_morph_kwargs or {}
    surfs = []
    for subject, subsampling in zip(
        (subject_from, subject_to), (subsampling_from, subsampling_to)
    ):
        if subject == "fsaverage":
            surfaces = load_reference_surfaces("sphere", subsampling)
        else:
            subject_files = (
                subject
                if isinstance(subject, SubjectFiles)
                else SubjectFiles(subpath=str(subject))
            )
            surfaces = load_subject_surfaces(subject_files, "sphere.reg", subsampling)
        surfs.append(surfaces)
    return {h: SurfaceMorph(surfs[0][h], surfs[1][h], **surface_morph_kwargs) for h in surfaces}


def middle_gm_interpolation(
    mesh_fn,
    m2m_folder,
    out_folder,
    out_fsaverage=None,
    depth=0.5,
    quantities=["magn", "normal", "tangent", "angle"],
    fields=None,
    open_in_gmsh=False,
    f_geo=None,
):
    """Interpolates the vector fieds in the middle gray matter surface

    Parameters
    -----------
    mesh_fn: str
        String with file name to mesh
    m2m_folder: str
        Path to the m2m_{subject_id} folder, generated during the segmantation
    out_folder: str
        Name of output folder. Output files will be written to this folder
    out_fsaverage: str (optional)
        Name of output folder for transformed files. In not set, fsaverage transormation
        will not be performed
    depth: float (optional)
        The distance bewteen grey and white matter where the
        interpolation should be done. p = depth * wm + (1 - depth) * gm (default: .5)
    quantities: list with the elements {magn, normal, tangent, angle}
        Quantites to be calculated from vector field
    fields: list of strings (optional)
        Fields to be transformed. Default: all fields
    open_in_gmsh: bool
        If true, opens a Gmsh window with the interpolated fields
    f_geo: str
        String with file name to geo file that accompanies the mesh
    """
    from ..mesh_tools import mesh_io

    m2m_folder = os.path.abspath(os.path.normpath(m2m_folder))
    subject_files = SubjectFiles(subpath=m2m_folder)
    if depth < 0.0 or depth > 1.0:
        raise ValueError("Invalid depth value. Should be between 0 and 1")

    if any([q not in ["magn", "normal", "tangent", "angle"] for q in quantities]):
        raise ValueError("Invalid quanty in {0}".format(quantities))

    def calc_quantities(nd, quantities):
        d = dict.fromkeys(quantities)
        for q in quantities:
            if q == "magn":
                d[q] = nd.norm()
            elif q == "normal":
                d[q] = nd.normal()
                d[q].value *= -1
            elif q == "tangent":
                d[q] = nd.tangent()
            elif q == "angle":
                d[q] = nd.angle()
            else:
                raise ValueError("Invalid quantity: {0}".format(q))
        return d

    m = mesh_io.read_msh(mesh_fn)
    _, sim_name = os.path.split(mesh_fn)
    sim_name = "." + os.path.splitext(sim_name)[0]

    # Crop out WM, GM, and CSF. We add WM and CSF to make the mesh convex.
    m = m.crop_mesh(tags=[ElementTags.WM, ElementTags.GM, ElementTags.CSF])

    # Set the volume to be GM. The interpolation will use only the tetrahedra in the volume.
    th_indices = m.elm.elm_number[m.elm.tag1 == ElementTags.GM]

    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)
    out_folder = os.path.abspath(os.path.normpath(out_folder))

    if out_fsaverage is not None and not os.path.isdir(out_fsaverage):
        os.mkdir(out_fsaverage)
    if out_fsaverage is not None:
        out_fsaverage = os.path.abspath(os.path.normpath(out_fsaverage))

    names_subj = []
    names_fsavg = []

    middle_surf = mesh_io.load_subject_surfaces(subject_files, "central")
    reg_surf = mesh_io.load_subject_surfaces(subject_files, "sphere.reg")

    ref_surf = mesh_io.load_reference_surfaces("sphere")
    avg_surf = mesh_io.load_reference_surfaces("central")

    for hemi in subject_files.hemispheres:
        mesh_io.write_freesurfer_surface(
            middle_surf[hemi],
            os.path.join(out_folder, hemi + ".central")
        )

    if out_fsaverage is not None:
        # Intersubject interpolators
        morph = {
            h: SurfaceMorph(reg_surf[h], ref_surf[h]) for h in subject_files.hemispheres
        }

    h = []
    for name, data in m.field.items():
        for hemi in subject_files.hemispheres:
            if fields is None or name in fields:
                # Interpolate to middle gm
                interpolated = data.interpolate_to_surface(
                    middle_surf[hemi], th_indices=th_indices
                )

                # For vector quantities, calculate quantities (normal, magn, ...)
                if data.nr_comp == 3:
                    q = calc_quantities(interpolated, quantities)
                    for q_name, q_data in q.items():
                        out_subj = os.path.join(
                            out_folder,
                            hemi + sim_name + ".central." + name + "." + q_name,
                        )
                        mesh_io.write_curv(
                            out_subj, q_data.value, middle_surf[hemi].elm.nr
                        )
                        names_subj.append(out_subj)
                        middle_surf[hemi].add_node_field(q_data, f"{name}_{q_name}")
                        h.append(hemi)
                        # Interpolate to fsavg
                        if out_fsaverage is not None:
                            q_transformed = morph[hemi].transform(q_data.value)
                            out_avg = os.path.join(
                                out_fsaverage,
                                hemi + sim_name + ".fsavg." + name + "." + q_name,
                            )
                            mesh_io.write_curv(
                                out_avg, q_transformed, ref_surf[hemi].elm.nr
                            )
                            avg_surf[hemi].add_node_field(
                                q_transformed, f"{name}_{q_name}"
                            )
                            names_fsavg.append(out_avg)

                elif data.nr_comp == 1:
                    field_name = name[-1]
                    q_name = name[:-1]
                    if field_name in m.field.keys() and q_name in quantities:
                        # If we have an equivalent quantity being calculated, skip
                        pass
                    else:
                        out_subj = os.path.join(
                            out_folder, hemi + sim_name + ".central." + name
                        )
                        mesh_io.write_curv(
                            out_subj,
                            interpolated.value.squeeze(),
                            middle_surf[hemi].elm.nr,
                        )
                        names_subj.append(out_subj)
                        h.append(hemi)
                        middle_surf[hemi].add_node_field(interpolated, name)
                        if out_fsaverage is not None:
                            f_transformed = morph[hemi].transform(
                                interpolated.value.squeeze()
                            )
                            out_avg = os.path.join(
                                out_fsaverage, hemi + sim_name + ".fsavg." + name
                            )
                            mesh_io.write_curv(
                                out_avg, f_transformed, ref_surf[hemi].elm.nr
                            )
                            names_fsavg.append(out_avg)
                            avg_surf[hemi].add_node_field(f_transformed, name)

    # Join surfaces, fields and open in gmsh
    # I only work with lh and rh at least for now
    # It also needs to be nicely ordered, otherwise will
    # screw up the atlases

    def join_and_write(surfs, fn_out, open_in_gmsh, f_geo=None):
        mesh = surfs["lh"].join_mesh(surfs["rh"])
        mesh.elm.tag1 = ElementTags.GM_TH_SURFACE * np.ones(mesh.elm.nr, dtype=int)
        mesh.elm.tag2 = ElementTags.GM_TH_SURFACE * np.ones(mesh.elm.nr, dtype=int)
        mesh.nodedata = []
        mesh.elmdata = []
        for k in surfs["lh"].field.keys():
            mesh.add_node_field(
                np.append(surfs["lh"].field[k].value, surfs["rh"].field[k].value), k
            )
        mesh_io.write_msh(mesh, fn_out)

        # write .opt-file
        v = mesh.view(visible_fields=list(surfs["lh"].field.keys())[0])
        if f_geo is not None:
            if not os.path.exists(f_geo):
                raise FileNotFoundError(f"Could not find file: {f_geo}")
            v.add_merge(f_geo, append_views_from_geo=True)
        v.write_opt(fn_out)

        if open_in_gmsh:
            mesh_io.open_in_gmsh(fn_out, True)

    join_and_write(
        middle_surf,
        os.path.join(out_folder, sim_name[1:] + "_central.msh"),
        open_in_gmsh,
        f_geo,
    )
    if out_fsaverage:
        join_and_write(
            avg_surf,
            os.path.join(out_fsaverage, sim_name[1:] + "_fsavg.msh"),
            open_in_gmsh,
        )


def subject_atlas(atlas_name, m2m_dir, hemi="both"):
    """Loads a brain atlas based of the FreeSurfer fsaverage template

    Parameters
    -----------
    atlas_name: 'a2009s', 'DK40' or 'HCP_MMP1'
            Name of atlas to load

            'a2009s': Destrieux atlas (FreeSurfer v4.5, aparc.a2009s)
            Cite: Destrieux, C. Fischl, B. Dale, A., Halgren, E. A sulcal
            depth-based anatomical parcellation of the cerebral cortex.
            Human Brain Mapping (HBM) Congress 2009, Poster #541

            'DK40': Desikan-Killiany atlas (FreeSurfer, aparc.a2005s)
            Cite: Desikan RS, S�gonne F, Fischl B, Quinn BT, Dickerson BC,
            Blacker D, Buckner RL, Dale AM, Maguire RP, Hyman BT, Albert MS,
            Killiany RJ. An automated labeling system for subdividing the
            human cerebral cortex on MRI scans into gyral based regions of
            interest. Neuroimage. 2006 Jul 1;31(3):968-80.

            'HCP_MMP1': Human Connectome Project (HCP) Multi-Modal Parcellation
            Cite: Glasser MF, Coalson TS, Robinson EC, et al. A multi-modal
            parcellation of human cerebral cortex. Nature. 2016;536(7615):171-178.

    m2m_folder: str
        Path to the m2m_{subject_id} folder, generated during the segmantation

    hemi (optional): 'lh', 'rh' or 'both'
        Hemisphere to use. In the case of 'both', will assume that left hemisphere
        nodes comes before right hemisphere nodes

    Returns
    ---------
    atlas: dict
        Dictionary where atlas['region'] = roi
    """
    from ..mesh_tools.mesh_io import read_gifti_surface

    if atlas_name not in ["a2009s", "DK40", "HCP_MMP1"]:
        raise ValueError("Invalid atlas name")

    subject_files = SubjectFiles(subpath=m2m_dir)

    if hemi in ["lh", "rh"]:
        fn_atlas = os.path.join(
            templates.atlases_surfaces, f"{hemi}.aparc_{atlas_name}.annot"
        )
        labels, _, names = nib.freesurfer.io.read_annot(fn_atlas)
        morph = SurfaceMorph(
            read_gifti_surface(get_reference_surf(hemi, "sphere")),
            read_gifti_surface(subject_files.surfaces["sphere.reg"][hemi]),
            method="nearest",  # we are interpolating labels
        )
        labels_sub = morph.transform(labels)
        atlas = {}
        for l, name in enumerate(names):
            atlas[name.decode()] = labels_sub == l

        return atlas

    # If both hemispheres
    elif hemi == "both":
        atlas_lh = subject_atlas(atlas_name, m2m_dir, "lh")
        atlas_rh = subject_atlas(atlas_name, m2m_dir, "rh")
        atlas = {}
        pad_rh = np.zeros_like(list(atlas_rh.values())[0])
        pad_lh = np.zeros_like(list(atlas_lh.values())[0])
        for name, mask in atlas_lh.items():
            atlas[f"lh.{name}"] = np.append(mask, pad_rh)  # pad after
        for name, mask in atlas_rh.items():
            atlas[f"rh.{name}"] = np.append(pad_lh, mask)  # pad after

        return atlas
    else:
        raise ValueError("Invalid hemisphere name")


def _project_points_to_surface(
    points: np.ndarray,
    surf: dict,
    pttris: Union[list, np.ndarray],
    return_all: bool = False,
):
    """Project each point in `points` to the closest point on the surface
    described by `surf` restricted to the triangles in `pttris`.

    PARAMETERS
    ----------
    points : ndarray
        Array with shape (n, d) where n is the number of points and d is the
        dimension.
    surf : dict
        Dictionary with keys 'points' and 'tris' describing the surface mesh on
        which the points are to be projected.
    pttris : ndarray | list
        If a ragged/nested array, the ith entry contains the triangles against
        which the ith point will be tested.
    return_all : bool
        Whether to return all projection results (i.e., the projection of a
        point on each of the triangles which it was tested against) or only the
        projection on the closest triangle.

    RETURNS
    -------
    tris : ndarray
        The index of the triangle onto which a point was projected.
    weights : ndarray
        The linear interpolation weights resulting in the projection of a point
        onto a particular triangle.
    projs :
        The coordinates of the projection of a point on a triangle.
    dists :
        The distance of a point to its projection on a triangle.

    NOTES
    -----
    The cost function to be minimized is the squared distance between a point
    P and a triangle T

        Q(s,t) = |P - T(s,t)|**2 =
               = a*s**2 + 2*b*s*t + c*t**2 + 2*d*s + 2*e*t + f

    The gradient

        Q'(s,t) = 2(a*s + b*t + d, b*s + c*t + e)

    is set equal to (0,0) to find (s,t).

    REFERENCES
    ----------
    https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf

    """

    npttris = list(map(len, pttris))
    pttris = np.concatenate(pttris)

    m = surf["points"][surf["tris"]]
    v0 = m[:, 0]  # Origin of the triangle
    e0 = m[:, 1] - v0  # s coordinate axis
    e1 = m[:, 2] - v0  # t coordinate axis

    # Vector from point to triangle origin (if reverse, the negative
    # determinant must be used)
    rep_points = np.repeat(points, npttris, axis=0)
    w = v0[pttris] - rep_points

    a = np.sum(e0**2, 1)[pttris]
    b = np.sum(e0 * e1, 1)[pttris]
    c = np.sum(e1**2, 1)[pttris]
    d = np.sum(e0[pttris] * w, 1)
    e = np.sum(e1[pttris] * w, 1)
    # f = np.sum(w**2, 1)

    # s,t are so far unnormalized!
    s = b * e - c * d
    t = b * d - a * e
    det = a * c - b**2

    # Project points (s,t) to the closest points on the triangle (s',t')
    sp, tp = np.zeros_like(s), np.zeros_like(t)

    # We do not need to check a point against all edges/interior of a triangle.
    #
    #          t
    #     \ R2|
    #      \  |
    #       \ |
    #        \|
    #         \
    #         |\
    #         | \
    #     R3  |  \  R1
    #         |R0 \
    #    _____|____\______ s
    #         |     \
    #     R4  | R5   \  R6
    #
    # The code below is equivalent to the following if/else structure
    #
    # if s + t <= 1:
    #     if s < 0:
    #         if t < 0:
    #             region 4
    #         else:
    #             region 3
    #     elif t < 0:
    #         region 5
    #     else:
    #         region 0
    # else:
    #     if s < 0:
    #         region 2
    #     elif t < 0
    #         region 6
    #     else:
    #         region 1

    # Conditions
    st_l1 = s + t <= det
    s_l0 = s < 0
    t_l0 = t < 0

    # Region 0 (inside triangle)
    i = np.flatnonzero(st_l1 & ~s_l0 & ~t_l0)
    deti = det[i]
    sp[i] = s[i] / deti
    tp[i] = t[i] / deti

    # Region 1
    # The idea is to substitute the constraints on s and t into F(s,t) and
    # solve, e.g., here we are in region 1 and have Q(s,t) = Q(s,1-s) = F(s)
    # since in this case, for a point to be on the triangle, s+t must be 1
    # meaning that t = 1-s.
    i = np.flatnonzero(~st_l1 & ~s_l0 & ~t_l0)
    aa, bb, cc, dd, ee = a[i], b[i], c[i], d[i], e[i]
    numer = cc + ee - (bb + dd)
    denom = aa - 2 * bb + cc
    sp[i] = np.clip(numer / denom, 0, 1)
    tp[i] = 1 - sp[i]

    # Region 2
    i = np.flatnonzero(~st_l1 & s_l0)  # ~t_l0
    aa, bb, cc, dd, ee = a[i], b[i], c[i], d[i], e[i]
    tmp0 = bb + dd
    tmp1 = cc + ee
    j = tmp1 > tmp0
    j_ = ~j
    k, k_ = i[j], i[j_]
    numer = tmp1[j] - tmp0[j]
    denom = aa[j] - 2 * bb[j] + cc[j]
    sp[k] = np.clip(numer / denom, 0, 1)
    tp[k] = 1 - sp[k]
    sp[k_] = 0
    tp[k_] = np.clip(-ee[j_] / cc[j_], 0, 1)

    # Region 3
    i = np.flatnonzero(st_l1 & s_l0 & ~t_l0)
    cc, ee = c[i], e[i]
    sp[i] = 0
    tp[i] = np.clip(-ee / cc, 0, 1)

    # Region 4
    i = np.flatnonzero(st_l1 & s_l0 & t_l0)
    aa, cc, dd, ee = a[i], c[i], d[i], e[i]
    j = dd < 0
    j_ = ~j
    k, k_ = i[j], i[j_]
    sp[k] = np.clip(-dd[j] / aa[j], 0, 1)
    tp[k] = 0
    sp[k_] = 0
    tp[k_] = np.clip(-ee[j_] / cc[j_], 0, 1)

    # Region 5
    i = np.flatnonzero(st_l1 & ~s_l0 & t_l0)
    aa, dd = a[i], d[i]
    tp[i] = 0
    sp[i] = np.clip(-dd / aa, 0, 1)

    # Region 6
    i = np.flatnonzero(~st_l1 & t_l0)  # ~s_l0
    aa, bb, cc, dd, ee = a[i], b[i], c[i], d[i], e[i]
    tmp0 = bb + ee
    tmp1 = aa + dd
    j = tmp1 > tmp0
    j_ = ~j
    k, k_ = i[j], i[j_]
    numer = tmp1[j] - tmp0[j]
    denom = aa[j] - 2 * bb[j] + cc[j]
    tp[k] = np.clip(numer / denom, 0, 1)
    sp[k] = 1 - tp[k]
    tp[k_] = 0
    sp[k_] = np.clip(-dd[j_] / aa[j_], 0, 1)

    # Distance from original point to its projection on the triangle
    projs = v0[pttris] + sp[:, None] * e0[pttris] + tp[:, None] * e1[pttris]
    dists = np.linalg.norm(rep_points - projs, axis=1)
    weights = np.column_stack((1 - sp - tp, sp, tp))

    if return_all:
        tris = pttris
    else:
        # Find the closest projection
        indptr = [0] + np.cumsum(npttris).tolist()
        i = _sliced_argmin(dists, indptr)
        tris = pttris[i]
        weights = weights[i]
        projs = projs[i]
        dists = dists[i]

    return tris, weights, projs, dists


def _get_nearest_triangles_on_surface(
    points: np.ndarray, surf: dict, n: int = 1, subset=None, return_index: bool = False
):
    """For each point in `points` get the `n` nearest nodes on `surf` and
    return the triangles to which these nodes belong.

    points : ndarray
        Points for which we want to find the candidate triangles. Shape (n, d)
        where n is the number of points and d is the dimension.
    surf : dict
        Dictionary with keys points and tris corresponding to the nodes and
        triangulation of the surface, respectively.
    n : int
        Number of nearest vertices in `surf` to consider for each point in
        `points`.
    subset : array-like
        Use only a subset of the vertices in `surf`. Should be indices *not* a
        boolean mask!
    return_index : bool
        Return the index (or indices if n > 1) of the nearest vertex in `surf`
        for each point in `points`.

    RETURNS
    -------
    pttris : list
        Point to triangle mapping.
    """
    assert isinstance(n, int) and n >= 1

    surf_points = surf["points"] if subset is None else surf["points"][subset]
    tree = scipy.spatial.cKDTree(surf_points)
    _, ix = tree.query(points, n)
    if subset is not None:
        ix = subset[ix]  # ensure ix indexes into surf['points']
    pttris = _get_triangle_neighbors(surf["tris"], len(surf["points"]))[ix]
    if n > 1:
        pttris = list(map(lambda x: np.unique(np.concatenate(x)), pttris))

    return (pttris, ix) if return_index else pttris


def _get_triangle_neighbors(tris: np.ndarray, nr: Union[int, None] = None):
    """For each point get its neighboring triangles (i.e., the triangles to
    which it belongs).

    PARAMETERS
    ----------
    tris : ndarray
        Array describing a triangulation with size (n, 3) where n is the number
        of triangles.
    nr : int
        Number of points. If None, it is inferred from `tris` as tris.max()+1
        (default = None).

    RETURNS
    -------
    pttris : ndarray
        Array of arrays where pttris[i] are the neighboring triangles of the
        ith point.
    """
    n_tris, n_dims = tris.shape
    nr = tris.max() + 1 if nr is None else nr

    rows = tris.ravel()
    cols = np.repeat(np.arange(n_tris), n_dims)
    data = np.ones_like(rows)
    csr = coo_matrix((data, (rows, cols)), shape=(nr, n_tris)).tocsr()
    return np.array(np.split(csr.indices, csr.indptr[1:-1]), dtype=object)


def _sliced_argmin(x: np.ndarray, indptr: np.ndarray):
    """Perform argmin on slices of x.

    PARAMETERS
    ----------
    x : 1-d array
        The array to perform argmin on.
    indptr : 1-d array-like
        The indices of the slices. The ith slice is indptr[i]:indptr[i+1].

    RETURNS
    -------
    res : 1-d array
        The indices (into x) corresponding to the minimum values in each chunk.
    """
    assert x.ndim == 1
    return np.array([x[i:j].argmin() + i for i, j in zip(indptr[:-1], indptr[1:])])


def create_new_connectivity_list_point_mask(points, con, point_mask):
    """
    Creates a new point and connectivity list when applying a point mask (changes indices of points)

    Parameters
    ----------
    points : np.ndarray of float [n_points x 3]
        Point coordinates
    con : np.ndarray of float [n_tri x 3]
        Connectivity of triangles
    point_mask : nparray of bool [n_points]
        Mask of (True/False) which points are kept in the mesh

    Returns
    -------
    points_new : np.ndarray of float [n_points_new x 3]
        New point array containing the remaining points after applying the mask
    con_new : np.ndarray of float [n_tri_new x 3]
        New connectivity list containing the remaining points (includes reindexing)
    """
    con_global = con[point_mask[con].all(axis=1), :]
    unique_points = np.unique(con_global)
    points_new = points[unique_points, :]

    con_new = np.zeros(con_global.shape).astype(int)

    for i, idx in enumerate(unique_points):
        idx_where = np.where(con_global == idx)
        con_new[idx_where[0], idx_where[1]] = i

    return points_new, con_new

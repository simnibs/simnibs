"""
Neuronavigation functions. So far, only for Localite TMS Navigator software.
This code is adapted from A. Thielscher's nnav_read_im.m and nnav_read_localite.m .

Written by Ole Numssen & Konstantin Weise, 2019.
"""
import copy
import io
import os
import nibabel
import numpy as np
from simnibs.utils.simnibs_logger import logger


def get_m_localite(fn_nii,
                   orientation='RAS', flirt_nav2conf=np.eye(4)):
    """Returns global transformation matrix for Localite TMS Navigator -> simnibs.
    
    If the image used for neuronavigation is not the SimNIBS _conform.nii image used for the mesh generation, 
    provide flirt_mat2conf transformation matrix. 

    Parameters
    ----------
    fn_nii: str
        Filename of .nii file the experiments were conducted with.
    orientation: str
        Orientation convention ('RAS' or 'LPS').
        Can be read from localite neuronavigation .xml file under coordinateSpace="RAS".

    flirt_nav2conf: np.ndarray or basestring
        Coregistration matrix for mesh.nii and neuronav.nii (Default: identity matrix).
    """
    if type(flirt_nav2conf) == str:
        flirt_nav2conf = np.loadtxt(flirt_nav2conf)

    fn_exp_nii_ras = splitext_niigz(fn_nii)[0] + '_RAS' + splitext_niigz(fn_nii)[1]

    # transform exp to RAS
    img = to_ras(fn_nii, fn_exp_nii_ras)

    # read q-form matrix from exp
    logger.debug('Constructing transformation matrices:')
    logger.debug(' > q-form matrix of exp')
    m_qform_exp = img.get_qform()

    # invert q-form matrix
    m_qform_exp_inv = np.linalg.inv(m_qform_exp)
    m_vox2mm = get_vox2mm(img.header)
    m_flirt2fs = get_flirt2fs(img.header)
    m_2ras = get_m_2ras(orientation)

    # combine transformation matrices
    m_total = np.dot(np.dot(np.dot(np.dot(m_flirt2fs, flirt_nav2conf), m_vox2mm), m_qform_exp_inv), m_2ras)
    return m_total


def get_vox2mm(exp_hdr):
    """Construct voxel to mm transformation matrix."""

    pixdim_exp = exp_hdr['pixdim'][1:4]
    m_vox2mm = np.eye(4)
    for i in range(3):
        m_vox2mm[i, i] = pixdim_exp[i]

    return m_vox2mm


def get_m_2ras(orientation):
    """Construct transformation matrix to transform into RAS space."""

    if orientation == 'LPS':
        m_2ras = np.array([[-1, 0, 0, 0],
                           [0, -1, 0, 0],
                           [0, 0, 1, 0],
                           [0, 0, 0, 1]])
    elif orientation == 'RAS':
        m_2ras = np.eye(4)

    else:
        raise NotImplementedError(f"orientation {orientation} is not implemented. Use RAS or LPS.")

    return m_2ras


def get_m_2conf(fn_exp_nii, fn_conform_nii, remove_tmp_files=True):
    """
    Calculate headmesh.nii to neuronavigation .nii transformation matrix with FSL FLIRT registration.


    If not remove_tmp_files, the matrix is also stored at fn_exp_nii + _m_2conform.mat

    Parameters
    ----------
    fn_exp_nii: str
        Full filename to .nii used in neuronavigation
    fn_conform_nii: str
        Full filename of .nii produced by SimNIBS for mesh generation.
    remove_tmp_files: bool
        Remove FSL intermediate files. Default: True.

    Returns
    -------
    np.ndarry (4x4)
        Transformation matrix from neuronavigation .nii to mesh conform .nii
    """
    # temp file names to store intermediate results from fsl
    fn_flip = splitext_niigz(fn_exp_nii)[0] + '_flipped_temp.nii'
    fn_out_fslmaths = splitext_niigz(fn_exp_nii)[0] + '_fslmaths_temp'
    fn_mat_m_2conform = splitext_niigz(fn_exp_nii)[0] + '_m_2conform'
    dof = 6

    # flip image of exp along first dimension and save it (to LAS, radiologic)
    exp_nii = nibabel.load(fn_exp_nii)
    data_exp_flipped = np.flip(exp_nii.get_data(), axis=0)
    exp_flipped_nii = nibabel.Nifti1Image(data_exp_flipped, exp_nii.affine, exp_nii.header)
    nibabel.save(exp_flipped_nii, fn_flip)

    # define FSL commands
    cmdstr = [[] for _ in range(4)]
    cmdstr[0] = 'fslorient -setqformcode 1 ' + fn_flip
    cmdstr[1] = 'fslorient -forceradiological ' + fn_flip
    cmdstr[2] = 'fslmaths ' + fn_conform_nii + ' -bin -s 1 ' + fn_out_fslmaths + '.nii.gz'
    cmdstr[3] = 'flirt -in ' + fn_flip + ' -ref ' + fn_conform_nii + ' -refweight ' + fn_out_fslmaths + \
                ' -searchrx -30 30 -searchry -30 30 -searchrz -30 30  -interp sinc -cost mutualinfo ' + \
                '-searchcost mutualinfo -dof ' + str(dof) + ' -omat ' + fn_mat_m_2conform + '.mat -out ' + \
                fn_mat_m_2conform + '.nii.gz'

    # execute FSL commands
    logger.debug('Executing fsl coregistration:')
    for i in range(len(cmdstr)):
        logger.debug('     > {}'.format(cmdstr[i]))
        os.system(cmdstr[i])
    mat2conf = np.loadtxt(fn_mat_m_2conform + '.mat')

    # remove temp files if they had not been there before
    if remove_tmp_files:
        os.remove(fn_flip)
        os.remove(fn_out_fslmaths + '.nii.gz')
        os.remove(fn_mat_m_2conform + '.nii.gz')
        os.remove(fn_mat_m_2conform + '.mat')

    return mat2conf


def get_flirt2fs(img_hdr):
    """Construct flirt to freesurfer transformation matrix"""

    pixdim_conform = img_hdr['pixdim'][1:4]
    dim_conform = img_hdr['dim'][1:4]

    m_flirt2fs = np.eye(4)
    m_flirt2fs[0, 3] = -pixdim_conform[0] * (dim_conform[0] / 2.0 - 1)
    m_flirt2fs[1, 3] = -pixdim_conform[1] * (dim_conform[1] / 2.0)
    m_flirt2fs[2, 3] = -pixdim_conform[2] * (dim_conform[2] / 2.0 - 1)

    return m_flirt2fs


def to_ras(fn_in, fn_out):
    """
    Transforming MRI .nii image to RAS space.

    Parameters
    ----------
    fn_in: str
        Filename of input image .nii file
    fn_out: str
        Filename of output image .nii file in RAS space

    Returns
    -------
    <File>: .nii file
        .nii image in RAS space (fn_out)
    """

    # read image data
    img_in = nibabel.load(fn_in)
    img_in_hdr = img_in.header
    img_out = copy.deepcopy(img_in)

    # read and invert q-form of original image
    m_qform_in = img_in.get_qform()
    m_qform_inv_in = np.linalg.inv(m_qform_in)

    # identify axes to flip
    mathlp = np.sign(m_qform_inv_in)

    ras_dim = np.zeros(3)
    ras_sign = np.zeros(3)

    for i in range(3):
        ras_dim[i] = np.where(np.abs(m_qform_inv_in[:, i]) == np.max(np.abs(m_qform_inv_in[:, i])))[0]
        ras_sign[i] = mathlp[int(ras_dim[i]), i]

    ras_dim = ras_dim.astype(int)

    # apply sorting to qform: first permute, then flip
    m_perm = np.zeros((4, 4))
    m_perm[3, 3] = 1

    for i in range(3):
        m_perm[ras_dim[i], i] = 1

    imgsize = img_in_hdr['dim'][1:4]
    imgsize = imgsize[ras_dim]

    m_flip = np.eye(4)

    for i in range(3):
        if ras_sign[i] < 0:
            m_flip[i, i] = -1
            m_flip[i, 3] = imgsize[i] - 1

    m_qform_out = np.dot(np.dot(m_qform_in, m_perm), m_flip)
    img_out.set_qform(m_qform_out)
    img_out.set_sform(m_qform_out)
    # m_toORG = np.dot(m_perm, m_flip)

    # apply sorting to image: first permute, then flip
    img_out_data = np.transpose(img_in.get_data(), ras_dim)

    for i in range(3):
        if ras_sign[i] < 0:
            img_out_data = np.flip(img_out_data, i)

    # save transformed image in .nii file
    img_out = nibabel.Nifti1Image(img_out_data, img_out.affine, img_out.header)
    nibabel.save(img_out, fn_out)

    return nibabel.load(fn_out)


def splitext_niigz(fn):
    """
    Splitting extension(s) from .nii or .nii.gz file

    Parameters
    ----------
    fn: str
        Filename of input image .nii or .nii.gz file

    Returns
    -------
    path: str
        Path and filename without extension(s)
    ext: str
        Extension, either .nii or .nii.gz
    """

    path, filename = os.path.split(fn)

    file0, ext0 = os.path.splitext(filename)

    if ext0 == '.gz':
        file1, ext1 = os.path.splitext(file0)
        return os.path.join(path, file1), ext1 + ext0
    elif ext0 == '.nii':
        return os.path.join(path, file0), ext0
    else:
        raise ValueError('File extension is neither .nii or .nii.gz!')


def write_tms_navigator_im(ims, xml_fn, overwrite=False):
    """
    Writes a instrument marker .xml file in the fashion of Localite TMS Navigator.

    Parameters
    ----------
    ims: np.ndarray
        Instrument markers with shape (4, 4, n_im).
    xml_fn: str
        Filename.
    overwrite: bool (Default: True)
        Overwrite existing file.

    Returns
    -------
    file: fn_out.xml
    """

    ims = np.atleast_3d(ims)

    assert ims.shape[:2] == (4, 4), 'Expecting array with shape (4, 4, N instrument marker).'

    if not xml_fn.lower().endswith('.xml'):
        xml_fn += '.xml'

    assert not os.path.exists(xml_fn) or overwrite, '.xml file already exists. Remove or set overwrite=True.'

    with io.open(xml_fn, 'w', newline='\n') as f:  # correct windows style would be \r\n, but Localite uses \n
        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<InstrumentMarkerList coordinateSpace="RAS">\n')
        f.write('    <!--This InstrumentMarkerList was written by SimNIBS-->\n')

        for idx in range(ims.shape[-1]):
            im = ims[:, :, idx]
            f.write('    ' + f'<InstrumentMarker alwaysVisible="false" index="{idx}" selected="false">\n')
            f.write('    ' * 2 + f'<Marker additionalInformation="" '
            f'color="#ff0000" description="opt_{idx}" set="true">\n')
            f.write('    ' * 3 + '<Matrix4D \n')
            f.write('    ' * 4 + 'data00="{:+1.17f}" data01="{:+1.17f}" '
                                 'data02="{:+1.17f}" data03="{:+1.17f}"\n'.format(im[0, 0], im[0, 1], im[0, 2],
                                                                                  im[0, 3]))
            f.write('    ' * 4 + 'data10="{:+1.17f}" data11="{:+1.17f}" '
                                 'data12="{:+1.17f}" data13="{:+1.17f}"\n'.format(im[1, 0], im[1, 1], im[1, 2],
                                                                                  im[1, 3]))
            f.write('    ' * 4 + 'data20="{:+1.17f}" data21="{:+1.17f}" '
                                 'data22="{:+1.17f}" data23="{:+1.17f}"\n'.format(im[2, 0], im[2, 1], im[2, 2],
                                                                                  im[2, 3]))
            f.write('    ' * 4 + 'data30="{:+1.17f}" data31="{:+1.17f}" '
                                 'data32="{:+1.17f}" data33="{:+1.17f}"/>\n'.format(0, 0, 0, 1))
            f.write('    ' * 2 + '</Marker>\n')
            f.write('    ' + '</InstrumentMarker>\n')

        f.write('</InstrumentMarkerList>\n')

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

from ..simulation import optim_tms
from simnibs.simulation import sim_struct
from simnibs.utils.simnibs_logger import logger


def get_matsimnibs_from_pos(pos, msh=None):
    """Returns matsimnibs from simnibs.simulation.sim_struct.POSITION as np.ndarray with shape = (4,4)"""
    if pos.matsimnibs is None:
        pos.calc_matsimnibs(msh)

    return np.squeeze(pos.matsimnibs)


def get_matsimnibs_from_tmslist(tmslist, msh=None):
    """Returns matsimnibs from all positions as np.ndarray with shape (4,4,n_pos)"""
    return np.stack([get_matsimnibs_from_pos(pos, msh=msh) for pos in tmslist.pos], axis=2)


def get_matsimnis_from_session(session, msh=None):
    """Returns matsimnibs from all positions from all poslists as np.ndarray with shape (4,4,N_pos)"""
    if msh is None:
        msh = session.fnamehead

    return np.squeeze(np.array([get_matsimnibs_from_tmslist(poslist, msh=msh) for poslist in session.poslists]))


def simnibs2nnav(fn_exp_nii, fn_conform_nii, simnibs_obj,
                 orientation='RAS', fsl_pref='',
                 manufacturer='localite', msh=None, flirt_mat2conf=np.eye(4), skip_flirt='auto'):
    """
    Transforms simnibs positions/orientations to neuronavigation instrument marker space.


    Parameters
    ----------
    fn_exp_nii: str
        Filename of .nii file the experiments were conducted with.
    fn_conform_nii: str
        Filename of .nii file from SimNIBS mri2msh function
        (e.g.: .../fs_subjectID/subjectID_T1fs_conform.nii.gz).
    simnibs_obj: simnibs.simulation.sim_struct object with position/s (POSITION | SESSION | TMSLIST | ...) or np.ndarray
        N positions/orientations of coil.
    orientation: str
        Orientation convention ('RAS' or 'LPS').
        Can be read from localite neuronavigation .xml file under coordinateSpace="RAS".
    fsl_pref: str
        Bash prefix needed to start FSL environment (Default: '').
    manufacturer: str (Default: 'localite')
        Do transforms for which software.
    msh: simnibs.simulation.sim_struct.Mesh or None (Default: None)
        Headmesh. Only needed if matsimnibs is not present in pos.
    flirt_mat2conf: np.ndarray or basestring
        Coregistration matrix for mesh.nii and neuronav.nii (Default: identity matrix).
    skip_flirt: str or bool
        Skip flirt coregistration between nnav.nii und mesh.nii. True, False, 'auto'. (Default: 'auto')

    Returns
    -------
    m_simnibs : nparray of float [4 x 4 x N]
    """

    if type(simnibs_obj) == np.ndarray:
        simnibs_obj_mat = simnibs_obj

    elif type(simnibs_obj) == sim_struct.TMSLIST:
        simnibs_obj_mat = get_matsimnibs_from_tmslist(simnibs_obj, msh=msh)

    elif type(simnibs_obj) == sim_struct.POSITION:
        simnibs_obj_mat = get_matsimnibs_from_pos(simnibs_obj, msh=msh)

    elif type(simnibs_obj) == optim_tms.TMSOPTIMIZATION:
        if msh is None:
            msh = simnibs_obj.fnamehead
        simnibs_obj_mat = get_matsimnibs_from_tmslist(simnibs_obj.optimlist, msh=msh)

    elif type(simnibs_obj) == sim_struct.SESSION:
        simnibs_obj_mat = get_matsimnis_from_session(simnibs_obj, msh=msh)

    else:
        raise NotImplementedError("{} not implemented in simnibs2nnav().".format(type(simnibs_obj)))

    if len(simnibs_obj_mat.shape) == 2:
        simnibs_obj_mat = simnibs_obj_mat[:, :, np.newaxis]

    if manufacturer == 'localite':

        m_total = get_m_localite(fn_exp_nii, fn_conform_nii, fsl_pref, orientation,
                                 skip_flirt=skip_flirt, flirt_mat2conf=flirt_mat2conf)
        m_total = np.linalg.inv(m_total)

        # construct flip matrix
        m_flip = np.array([[0, 0, 1, 0],
                           [0, -1, 0, 0],
                           [1, 0, 0, 0],
                           [0, 0, 0, 1]])

        # transform coil position from neuronavigation to simnibs space
        m_nnav = np.dot(np.dot(m_total,
                               simnibs_obj_mat.transpose([2, 0, 1])
                               ).transpose([1, 0, 2]),
                        m_flip).transpose([1, 2, 0])

    else:
        raise NotImplementedError

    return m_nnav


def nnav2simnibs(fn_exp_nii, fn_conform_nii, m_nnav, orientation='RAS',
                 fsl_pref='', manufacturer='localite'):
    """
    Transforms TMS coil instrument markers from neuronavigation to simnibs space.

    Parameters
    ----------
    fn_exp_nii: str
        Filename of .nii file the experiments were conducted with.
    fn_conform_nii: str
        Filename of .nii file from SimNIBS mri2msh function
        (e.g.: .../fs_subjectID/subjectID_T1fs_conform.nii.gz).
    m_nnav: nparray [4 x 4 x N]
        N position matrices from neuronavigation.
    orientation: str
        Orientation convention ('RAS' or 'LPS').
        Can be read from localite neuronavigation .xml file under coordinateSpace="RAS".
    fsl_pref: str
        Bash prefix needed to start FSL environment (Default: '').
    manufacturer: str (Default: 'localite')
        Do transforms for which software.

    Returns
    -------
    m_simnibs: nparray of float [4 x 4 x N]

    """

    if len(m_nnav.shape) == 2:
        m_nnav = m_nnav[:, :, np.newaxis]

    if manufacturer == 'localite':

        m_total = get_m_localite(fn_exp_nii, fn_conform_nii, fsl_pref, orientation)

        # construct flip matrix
        m_flip = np.array([[0, 0, 1, 0],
                           [0, -1, 0, 0],
                           [1, 0, 0, 0],
                           [0, 0, 0, 1]])

        # transform coil position from neuronavigation to simnibs space
        m_simnibs = np.dot(np.dot(m_total,
                                  m_nnav.transpose([2, 0, 1])
                                  ).transpose([1, 0, 2]),
                           m_flip).transpose([1, 2, 0])

    else:
        raise NotImplementedError

    return m_simnibs


def get_m_localite(fn_exp_nii, fn_conform_nii,
                   fsl_pref='', orientation='RAS',
                   skip_flirt='auto', flirt_mat2conf=np.eye(4)):
    """Returns global transformation matrix for Localite TMS Navigator -> simnibs.

    Parameters
    ----------
    fn_exp_nii: str
        Filename of .nii file the experiments were conducted with.
    fn_conform_nii: str
        Filename of .nii file from SimNIBS mri2msh function
        (e.g.: .../fs_subjectID/subjectID_T1fs_conform.nii.gz).
    fsl_pref: str
        Bash prefix needed to start FSL environment (Default: '').
    orientation: str
        Orientation convention ('RAS' or 'LPS').
        Can be read from localite neuronavigation .xml file under coordinateSpace="RAS".
    skip_flirt: 'auto' or bool
        Skip coregistration between headmeash.nii and neuronavigation.nii
        'auto': same base-filenames and same qforms -> Skip
        True: skip. False: do flirt coregistration.
    flirt_mat2conf: np.ndarray or basestring
        Coregistration matrix for mesh.nii and neuronav.nii (Default: identity matrix).
    """
    if type(flirt_mat2conf) == str:
        flirt_mat2conf = np.loadtxt(flirt_mat2conf)

    # get original qform without RAS
    exp_nii_original = nibabel.load(fn_exp_nii)
    conform_nii_original = nibabel.load(fn_conform_nii)
    m_qform_exp_original = exp_nii_original.get_qform()
    m_qform_conform_original = conform_nii_original.get_qform()

    # check if conform_nii and exp_nii are the same and have the same q-form
    if skip_flirt == 'auto':
        skip_flirt = (os.path.split(splitext_niigz(fn_exp_nii)[0])[1] ==
                      os.path.split(splitext_niigz(fn_conform_nii)[0])[1]) \
                     and np.all(np.isclose(m_qform_conform_original, m_qform_exp_original)) \
                     and np.all(flirt_mat2conf == np.eye(4))

    fn_exp_nii_ras = splitext_niigz(fn_exp_nii)[0] + '_RAS' + splitext_niigz(fn_exp_nii)[1]

    # transform exp to RAS
    exp_nii = to_ras(fn_exp_nii, fn_exp_nii_ras)
    conform_nii = nibabel.load(fn_conform_nii)
    logger.debug('Gathering header information...')

    # extract header information
    conform_hdr = conform_nii.header
    exp_hdr = exp_nii.header

    # read q-form matrix from exp
    logger.debug('Constructing transformation matrices:')
    logger.debug(' > q-form matrix of exp')
    m_qform_exp = exp_nii.get_qform()

    # read q-form matrix from conform
    # m_qform_conform = conform_nii.get_qform()

    # invert q-form matrix
    m_qform_exp_inv = np.linalg.inv(m_qform_exp)
    m_vox2mm = get_vox2mm(exp_hdr)
    m_flirt2fs = get_flirt2fs(conform_hdr)
    m_2ras = get_m_2ras(orientation)

    # construct flirt transformation matrix if necessary
    if skip_flirt:
        logger.debug('Skipping flirt coregistration of headmesh.nii and neuronavigation.nii.')
        m_2conf = flirt_mat2conf
    else:
        m_2conf = get_m_2conf(exp_nii, fn_conform_nii, fn_exp_nii_ras, fsl_pref)

    # combine transformation matrices
    m_total = np.dot(np.dot(np.dot(np.dot(m_flirt2fs, m_2conf), m_vox2mm), m_qform_exp_inv), m_2ras)
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
    else:
        m_2ras = np.eye(4)

    return m_2ras


def get_m_2conf(exp_nii, fn_conform_nii, fn_exp_nii_ras, fsl_cmd):
    """
    Construct headmesh.nii to neuronavigation .nii transformation matrix


    """
    # temo file names to store intermediate results from fsl
    fn_flip = splitext_niigz(fn_exp_nii_ras)[0] + '_flipped_temp.nii'
    fn_out_fslmaths = splitext_niigz(fn_exp_nii_ras)[0] + '_fslmaths_temp'
    fn_mat_m_2conform = splitext_niigz(fn_exp_nii_ras)[0] + '_m_2conform_temp'
    if any([os.path.exists(fn_flip), os.path.exists(fn_out_fslmaths), os.path.exists(fn_mat_m_2conform)]):
        remove_tmp_files = False
    else:
        remove_tmp_files = True

    dof = 6

    # flip image of exp along first dimension and save it (to LAS, radiologic)
    data_exp_flipped = np.flip(exp_nii.get_data(), axis=0)
    exp_flipped_nii = nibabel.Nifti1Image(data_exp_flipped, exp_nii.affine, exp_nii.header)
    nibabel.save(exp_flipped_nii, fn_flip)

    # define FSL commands
    cmdstr = [[] for _ in range(4)]
    cmdstr[0] = fsl_cmd + ' fslorient -setqformcode 1 ' + fn_flip
    cmdstr[1] = fsl_cmd + ' fslorient -forceradiological ' + fn_flip
    cmdstr[2] = fsl_cmd + ' fslmaths ' + fn_conform_nii + ' -bin -s 1 ' + fn_out_fslmaths + '.nii.gz'
    cmdstr[3] = fsl_cmd + ' flirt -in ' + fn_flip + ' -ref ' + fn_conform_nii + ' -refweight ' + fn_out_fslmaths + \
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


def get_flirt2fs(conform_hdr):
    """Construct flirt to freesurfer transformation matrix"""

    pixdim_conform = conform_hdr['pixdim'][1:4]
    dim_conform = conform_hdr['dim'][1:4]

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

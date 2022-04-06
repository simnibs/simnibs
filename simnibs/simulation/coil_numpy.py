# -*- coding: utf-8 -*-\
'''
    Functions for calculating dA/dt fields in mesh nodes from .ccd or nifti coil files
'''

import numpy as np
import nibabel as nib
import os
import re
from .. import __version__
from ..mesh_tools import mesh_io
from ..utils.file_finder import Templates


try:
    import fmm3dpy
except ImportError:
    FMM3D = False
else:
    FMM3D = True


'''
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2013-2019  Andre Antunes, Guilherme B Saturnino, Kristoffer H Madsen


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

def parseccd(ccd_file):
    '''
    Parse ccd file, and return intended bounding box, resolution and dIdtmax
        and other fields if available

    Parameters
    ----------
    ccd_file : string
        ccd file to parse

    On ccd file format version 1.1:
    1) First line is a header line escaped with # which can contain any number
       of variables in the form variable=value, a few of these variables are
       reserved as can be seen below, they are separated by semicolons (;).
    2) The second line is the contains the number of dipoles expected
       (this is actually not used in practice).
    3) Third line contains a header text excaped by #, typically:
       # centers and weighted directions of the elements (magnetic dipoles)
    4-end) Remaining lines are space separated dipole positions and dipole
       moments in a string number format readable by numpy. E.g. each line must contain
       six values: x y z mx my mz, where the first three are x,y and positions
       in meters, and the remaining three are dipole moments in x,y and z direction
       in Coulumb * meter per 1 A/s input current.
       an example could be:
       0 0 0 0 0 1.0e-03
       indicating a dipole at position 0,0,0 in z direction with strength
       0.001 C*m*s/A

    The variables are used to encode additonal optional information in text, specifically:
    dIdtmax=147,100
        which would indicate a max dI/dt (at 100% MSO) of 146.9 A/microsecond
        for first stimulator and 100 A/microsecond for the second stimulator.
    dIdtstim=162
        Indicating the max dI/dt reported on the stimulation display of 162,
        typically this is used to create a rescaled version of the ccd file
        such that the stimulator reported dI/dt max can be used directly.
        This is currently only supported for one stimulator.
    stimulator=Model name 1,Model name 2
    brand=Brand name 1,Brand name 2
    coilname=name of coil
        Indicates the name of the coil for display purposes
    Some variables are used for expansion in to nifti1 format:
    x=-300,300
        Indicates that the ccd file should be expanded into a FOV from
        x=-300mm to x=300mm, this could also be indicated as x=300
    y=-300,300
        The same for y
    z=-200,200
        The same for z
    resolution=3,3,3
        Indicates that the resolution should be 3mm in x,y and z directions,
        this could also be given as resolution=3

    The below is an example header line:
    #Test CCD file;dIdtmax=162;x=300;y=300;z=200;resolution=3;stimulator=MagProX100;brand=MagVenture;


    '''
    def parseField(info, field):
        try:
            a = np.fromstring(info[field],sep=',')
        except:
            a = None
        return a
    if os.path.splitext(ccd_file)[1]=='.ccl':
        d_position, d_moment = read_ccl(ccd_file)
    else:
        d_position, d_moment = read_ccd(ccd_file)

    #reopen to read header
    f = open(ccd_file,'r')
    data = f.readline()
    fields = re.findall(r'(\w*=[^;]*)', data + ';')
    labels = [f.split('=')[0] for f in fields]
    values = [f.split('=')[1].rstrip('\n') for f in fields]
    info = {}
    for i,label in enumerate(labels):
        info[label]=values[i]

    #parse bounding box for nii
    bb = []
    for dim in ('x','y','z'):
        a = parseField(info, dim)
        if a is None:
            bb.append(None)
        else:
            if len(a)<2:
                bb.append((-np.abs(a),np.abs(a)))
            else:
                bb.append(a)

    #parse resolution
    res = []
    a = parseField(info, 'resolution')
    if a is None:
        res.append(None)
    else:
        if len(a)<3:
            for i in range(len(a),3):
                a = np.concatenate((a, (a[i-1],)))
        res = a
    return d_position, d_moment, bb, res, info

def read_ccd(fn):
    """ reads a ccl file, this format is similar to the ccd format. However,
    only line segments positions are included (first 3 columns) and an optional
    weighting in the forth column

    Parameters
    -----------
    fn: str
        name of ccl file

    Returns
    ----------
    [pos, m]: list
        positions of line segments
    """
    ccl_file = np.loadtxt(fn, skiprows=2)

    # if there is only 1 dipole, loadtxt return as array of the wrong shape
    if (len(np.shape(ccl_file)) == 1):
        a = np.zeros([1, 4])
        a[0, 0:3] = ccd_file[0:3]
        a[0, 3:] = ccd_file[3:]
        ccd_file = a

    return ccd_file[:, 0:3], ccd_file[:, 3:]

def read_ccd(fn):
    """ reads a ccd file

    Parameters
    -----------
    fn: str
        name of ccd file

    Returns
    ----------
    [pos, m]: list
        position and moment of dipoles
    """
    ccd_file = np.loadtxt(fn, skiprows=2)

    # if there is only 1 dipole, loadtxt return as array of the wrong shape
    if (len(np.shape(ccd_file)) == 1):
        a = np.zeros([1, 6])
        a[0, 0:3] = ccd_file[0:3]
        a[0, 3:] = ccd_file[3:]
        ccd_file = a

    return ccd_file[:, 0:3], ccd_file[:, 3:]


def A_from_dipoles(d_moment, d_position, target_positions, eps=1e-3, direct='auto'):
    '''
    Get A field from dipoles using FMM3D

    Parameters
    ----------
    d_moment : ndarray
        dipole moments (Nx3).
    d_position : ndarray
        dipole positions (Nx3).
    target_positions : ndarray
        position for which to calculate the A field.
    eps : float
        Precision. The default is 1e-3
    direct : bool
        Set to true to force using direct (naive) approach or False to force use of FMM.
        If set to auto direct method is used for less than 300 dipoles which appears to be faster i these cases.
        The default is 'auto'

    Returns
    -------
    A : ndarray
        A field at points (M x 3) in Tesla*meter.

    '''
    #if set to auto use direct methods if # dipoles less than 300
    if direct=='auto':
        if d_moment.shape[0]<300:
            direct = True
        else:
            direct = False
    if direct is True:
        out = fmm3dpy.l3ddir(charges=d_moment.T, sources=d_position.T,
                  targets=target_positions.T, nd=3, pgt=2)
    elif direct is False:
        #use fmm3dpy to calculate expansion fast
        out = fmm3dpy.lfmm3d(charges=d_moment.T, eps=eps, sources=d_position.T,
                  targets=target_positions.T, nd=3, pgt=2)
    else:
        print('Error: direct flag needs to be either "auto", True or False')
    A = np.empty((target_positions.shape[0], 3), dtype=float)
    #calculate curl
    A[:, 0] = (out.gradtarg[1][2] - out.gradtarg[2][1])
    A[:, 1] = (out.gradtarg[2][0] - out.gradtarg[0][2])
    A[:, 2] = (out.gradtarg[0][1] - out.gradtarg[1][0])
    #scale
    A *= -1e-7
    return A

def B_from_dipoles(d_moment, d_position, target_positions, eps=1e-3, direct='auto'):
    '''
    Get B field from dipoles using FMM3D

    Parameters
    ----------
    d_moment : ndarray
        dipole moments (Nx3).
    d_position : ndarray
        dipole positions (Nx3).
    target_positions : ndarray
        position for which to calculate the B field.
    eps : float
        Precision. The default is 1e-3
    direct : bool
        Set to true to force using direct (naive) approach or False to force use of FMM.
        If set to auto direct method is used for less than 300 dipoles which appears to be faster i these cases.
        The default is 'auto'

    Returns
    -------
    B : ndarray
        B field at points (M x 3) in Tesla.

    '''
    #if set to auto use direct methods if # dipoles less than 300
    if direct=='auto':
        if d_moment.shape[0]<300:
            direct = True
        else:
            direct = False
    if direct is True:
        out = fmm3dpy.l3ddir(dipvec=d_moment.T, sources=d_position.T,
                  targets=target_positions.T, nd=1, pgt=2)
    elif direct is False:
        out = fmm3dpy.lfmm3d(dipvec=d_moment.T, eps=eps, sources=d_position.T,
                  targets=target_positions.T, nd=1, pgt=2)
    else:
        print('Error: direct flag needs to be either "auto", True or False')
    B = out.gradtarg.T
    B *= -1e-7
    return B

def A_biot_savart_path(lsegments, points, I=1.):
    '''
    Calculate A field naively (without FMM) from current path

    Parameters
    ----------
    lsegments : ndarray
        line segments, definition of M wire path as 3D points (3 x M).
    points : ndarray
        3D points where field is calculated (3 x N).
    I : float, optional
        Current strenth in A. The default is 1..

    Returns
    -------
    A : ndarray
        A field at points (3 x N) in Tesla*meter..

    '''
    r1 = np.sqrt(np.sum((points[:,:,None]-lsegments[:,None,:])**2,axis=0))
    dl = np.zeros(lsegments.shape)
    dl[:,:-1] = np.diff(lsegments,axis=1)
    dl[:,-1] = lsegments[:,0]-lsegments[:,-1]
    A = np.sum(dl[:,None,:]/r1[None],axis=2)
    A *= I*1e-7
    return A

def A_biot_savart_path_fmm(lsegments, points, I=1., eps=1e-3, direct='auto'):
    '''
    Calculate A field using FMM3D from current path

    Parameters
    ----------
    lsegments : ndarray
        line segments, definition of M wire path as 3D points (3 x M).
    points : ndarray
        3D points where field is calculated (3 x N).
    I : float, optional
        Current strenth in A. The default is 1..
    eps : float
        Precision. The default is 1e-3
    direct : bool
        Set to true to force using direct (naive) approach or False to force use of FMM.
        If set to auto direct method is used for less than 300 line segments which appears to be faster i these cases.
        The default is 'auto'

    Returns
    -------
    A : ndarray
        A field at points (3 x N) in Tesla*meter.

    '''
    dl = np.zeros(lsegments.shape)
    dl[:,:-1] = np.diff(lsegments,axis=1)
    dl[:,-1] = lsegments[:,0]-lsegments[:,-1]
    if direct == 'auto':
        if dl.shape[1] >= 300:
            direct = False
        else:
            direct = True
    if direct:
        A = fmm3dpy.l3ddir(sources=lsegments, charges=dl, targets=points, nd=3, pgt=1)
    else:
        A = fmm3dpy.lfmm3d(sources=lsegments, charges=dl, targets=points, nd=3, eps=eps, pgt=1)
    A = I*1e-7*A.pottarg
    return A

def B_biot_savart_path(lsegments, points, I=1.):
    '''
    Calculate B field naively (without FMM) from current path

    Parameters
    ----------
    lsegments : ndarray
        line segments, definition of M wire path as 3D points (3 x M).
    points : ndarray
        3D points where field is calculated (3 x N).
    I : float, optional
        Current strenth in Ampere. The default is 1..

    Returns
    -------
    B : ndarray
        B field at points (3 x N) in Tesla.

    '''
    r = points[:,:,None]-lsegments[:,None,:]
    r3 = np.sum(r**2,axis=0)**1.5
    dl = np.zeros(lsegments.shape)
    dl[:,:-1] = np.diff(lsegments,axis=1)
    dl[:,-1] = lsegments[:,0]-lsegments[:,-1]
    dl = dl[:,None,:]
    B = np.sum(np.cross(dl, r, axis=0)/r3[None],axis=2)
    B *= I*1e-7
    return B

def _rotate_coil(ccd_file, coil_matrix):
    # read ccd file
    d_position, d_moment = read_ccd(ccd_file)
    # transfrom positions to mm
    d_position *= 1e3
    # add a column to the position in order to apply the transformation matrix
    d_position = np.hstack([d_position, np.ones((d_position.shape[0], 1))])
    d_position = coil_matrix.dot(d_position.T).T[:, :3]
    # rotate the moment
    d_moment = coil_matrix[:3, :3].dot(d_moment.T).T
    return d_position, d_moment

def _calculate_dadt_ccd(msh, ccd_file, coil_matrix, didt, geo_fn):
    """ auxiliary function to calculate the dA/dt field from a ccd file """
    # read ccd file
    d_position, d_moment = _rotate_coil(ccd_file, coil_matrix)
    A = np.zeros((msh.nodes.nr, 3), dtype=float)
    for p, m in zip(d_position, d_moment):
        # get distance of point to dipole, transform back to meters
        r = (msh.nodes.node_coord - p) * 1e-3
        A += 1e-7 * didt * np.cross(m, r) / (np.linalg.norm(r, axis=1)[:, None] ** 3)
    node_data = mesh_io.NodeData(A)

    if geo_fn is not None:
        mesh_io.write_geo_spheres(
            d_position, geo_fn,
            np.linalg.norm(d_moment, axis=1),
            'coil_dipoles'
        )

    return node_data

def _calculate_dadt_ccd_FMM(msh, ccd_file, coil_matrix, didt, geo_fn, eps=1e-3):
    """ auxiliary function to calculate the dA/dt field from a ccd file using FMM """
    d_position, d_moment = _rotate_coil(ccd_file, coil_matrix)
    # bring everything to SI
    d_position *= 1e-3
    pos = msh.nodes[:] * 1e-3
    A = np.zeros((len(pos), 3), dtype=float)
    out = fmm3dpy.lfmm3d(
            eps=eps,
            sources=d_position.T,
            charges=d_moment.T,
            targets=pos.T,
            pgt=2,
            nd=3
        )

    A[:, 0] = (out.gradtarg[1][2] - out.gradtarg[2][1])
    A[:, 1] = (out.gradtarg[2][0] - out.gradtarg[0][2])
    A[:, 2] = (out.gradtarg[0][1] - out.gradtarg[1][0])

    A *= -1e-7 * didt
    if geo_fn is not None:
        mesh_io.write_geo_spheres(
            d_position*1e3, geo_fn,
            np.linalg.norm(d_moment, axis=1),
            'coil_dipoles'
        )

    return mesh_io.NodeData(A)


def _calculate_dadt_nifti(msh, nifti_image, coil_matrix, didt, geo_fn):

    """ auxiliary function that interpolates the dA/dt field from a nifti file """
    if isinstance(nifti_image, str):
        nifti_image = nib.load(nifti_image)
    elif isinstance(nifti_image, nib.nifti1.Nifti1Image):
        pass
    else:
        raise NameError('Failed to parse input volume (not string or nibabel nifti1 volume)')
    coords = msh.nodes.node_coord

    out = _get_field(nifti_image, coords, coil_matrix)
    out = out * didt

    node_data = mesh_io.NodeData(out.T)

    if geo_fn is not None:
        y_axis = np.arange(1, 10, dtype=float)[:, None] * (0, 1, 0)
        z_axis = np.arange(1, 30, dtype=float)[:, None] * (0, 0, 1)
        pos = np.vstack((((0, 0, 0)), y_axis, z_axis))
        pos = (coil_matrix[:3, :3].dot(pos.T) + coil_matrix[:3, 3][:, None]).T
        mesh_io.write_geo_spheres(pos, geo_fn,
                               name='coil_directions')

    return node_data


def _get_field(nifti_image, coords, coil_matrix, get_norm=False):
    ''' This function is also used in the GUI '''
    from scipy.ndimage import interpolation
    if isinstance(nifti_image, str):
        nifti_image = nib.load(nifti_image)
    elif isinstance(nifti_image, nib.nifti1.Nifti1Image):
        pass
    else:
        raise NameError('Failed to parse input volume (not string or nibabel nifti1 volume)')
    iM = np.dot(np.linalg.pinv(nifti_image.affine),
                np.linalg.pinv(coil_matrix))

    # gets the coordinates in voxel space
    voxcoords = np.dot(iM[:3, :3], coords.T) + iM[:3, 3][:, np.newaxis]

    # Interpolates the values of the field in the given coordinates
    out = np.zeros((3, voxcoords.shape[1]))
    for dim in range(3):
        out[dim] = interpolation.map_coordinates(
            np.asanyarray(nifti_image.dataobj)[..., dim],
            voxcoords, order=1
        )

    # Rotates the field
    out = np.dot(coil_matrix[:3, :3], out)
    if get_norm:
        out = np.linalg.norm(out, axis=0)
    return out


def _add_logo(msh_stl):
    ''' adds the simnibs logo to the coil surface '''

    msh_logo=mesh_io.read_msh(Templates().simnibs_logo)

    # 'simnibs' has tag 1, '3' has tag 2, '4' has tag 3
    # renumber tags, because they will be converted to color:
    # 0 gray, 1 red, 2 lightblue, 3 blue
    major_version = __version__.split('.')[0]
    if  major_version == '3':
        msh_logo=msh_logo.crop_mesh(tags = [1,2])
        msh_logo.elm.tag1[msh_logo.elm.tag1==2] = 3 # version in blue
    elif major_version == '4':
        msh_logo=msh_logo.crop_mesh(tags=[1,3])
    else:
        msh_logo=msh_logo.crop_mesh(tags=1)
    msh_logo.elm.tag1[msh_logo.elm.tag1==1] = 2 # 'simnibs' in light blue

    # center logo in xy-plane, mirror at yz-plane and scale
    bbox_coil=np.vstack([np.min(msh_stl.nodes[:],0),
                         np.max(msh_stl.nodes[:],0)])
    bbox_logo=np.vstack([np.min(msh_logo.nodes[:],0),
                         np.max(msh_logo.nodes[:],0)])
    bbox_ratio=np.squeeze(np.diff(bbox_logo,axis=0)/
                          np.diff(bbox_coil,axis=0))
    bbox_ratio=max(bbox_ratio[0:2]) # maximal size ratio in xy plane

    msh_logo.nodes.node_coord[:,0:2] -= np.mean(bbox_logo[:,0:2],axis=0)
    msh_logo.nodes.node_coord[:,0] = -msh_logo.nodes.node_coord[:,0]
    msh_logo.nodes.node_coord[:,0:2] *= 1/(4*bbox_ratio)

    # shift logo along negative z to the top side of coil
    msh_logo.nodes.node_coord[:,2] += bbox_coil[0,2] - bbox_logo[0,2] - 5

    msh_stl = msh_stl.join_mesh(msh_logo)
    return msh_stl


def _transform_coil_surface(msh_stl, coil_matrix, msh_skin=None, add_logo=False):
    """ Applies a 4x4 transformation matrix to mesh nodes.

        Parameters
        -----------
        msh_stl: simnibs.mesh_io.Msh
            mesh structure with the coil surface
        coil_matrix: np.ndarray (4x4)
            transformation matrix
        msh_skin: simnibs.mesh_io.Msh (optional)
            Mesh structure with the skin. When given, it will be tested
            whether any nodes of msh_stl are inside msh_skin. The triangles
            with these nodes will be marked by setting tag1=1 (other: tag1=0).
        add_logo:
            if set to True, a SimNIBS logo will be added to the coil surface.
            The tags will be marked as 2 and 3.

        Returns
        -------
        msh_stl: simnibs.mesh_io.Msh
            transfomred coil surface
    """
    msh_stl.elm.tag1[:]=0 # tag1=0 will be shown as gray
    if add_logo:
        msh_stl = _add_logo(msh_stl)
    d_position = np.hstack([msh_stl.nodes[:], np.ones((msh_stl.nodes.nr, 1))])
    msh_stl.nodes.node_coord = coil_matrix.dot(d_position.T).T[:, :3]

    if msh_skin is not None:
        idx_inside = msh_skin.pts_inside_surface(msh_stl.nodes[:])
        if len(idx_inside):
            idx_hlp = np.zeros((msh_stl.nodes.nr, 1),dtype=bool)
            idx_hlp[idx_inside] = True
            idx_hlp = np.any(np.squeeze(idx_hlp[msh_stl.elm[:,:3]-1]), axis=1)
            msh_stl.elm.tag1[idx_hlp & (msh_stl.elm.tag1 == 0)] = 1
    return msh_stl


def set_up_tms_dAdt(msh, coil_file, coil_matrix, didt=1e6, fn_geo=None, fn_stl=None):

    add_logo = True # add simnibs logo to coil surface visulization

    coil_matrix = np.array(coil_matrix)
    if isinstance(coil_file, nib.nifti1.Nifti1Image) or\
        coil_file.endswith('.nii.gz') or\
        coil_file.endswith('.nii'):
        dadt = _calculate_dadt_nifti(msh, coil_file,
                                     coil_matrix, didt, fn_geo)
    elif coil_file.endswith('.ccd'):
        if FMM3D:
            dadt = _calculate_dadt_ccd_FMM(
                msh, coil_file,
                coil_matrix, didt, fn_geo
            )
        else:
            dadt = _calculate_dadt_ccd(
                msh, coil_file,
                coil_matrix, didt, fn_geo
            )
    else:
        raise ValueError('coil file must be either a .ccd file or a nifti file')

    if (fn_geo is not None) and (fn_stl is not None):
        idx = (msh.elm.elm_type == 2)&( (msh.elm.tag1 == 1005) |
                                        (msh.elm.tag1 == 5) )
        msh_skin = msh.crop_mesh(elements = msh.elm.elm_number[idx])
        if not os.path.isfile(fn_stl):
            raise IOError('Could not find stl file: {0}'.format(fn_stl))
        msh_stl = mesh_io.read_stl(fn_stl)
        msh_stl = _transform_coil_surface(msh_stl,coil_matrix,msh_skin,add_logo)
        mesh_io.write_geo_triangles(msh_stl.elm[:,:3]-1, msh_stl.nodes[:],
                                    fn_geo, values=msh_stl.elm.tag1,
                                    name='coil_casing', mode='ba')
    dadt.mesh = msh
    dadt = dadt.node_data2elm_data()
    return dadt


def set_up_tms(msh, coil_file, coil_matrix, didt=1e6, fn_geo=None, fn_stl=None):
    """ sets up a tms simulation

    Parameters
    ------------
    msh: simnibs.mesh_io.Msh
        mesh structure
    coil_file: str
        name of file with coil
    coil_matrix: np.ndarray
        4x4 array with coil rotation and center
    didt(optional): float
        value of dI/dt
   fn_geo(optional): str
        name of .geo file with coil representation
        default: don't output a geo file
   fn_stl(optional): str
        name of .stl file of coil casing
        that will be added to fn_geo
        default: None

   Returns
   ---------
    dAdt: simnibs.mesh_io.ElementData
        dAdt at the elements
    """
    coil_matrix = np.array(coil_matrix, dtype=float)
    return set_up_tms_dAdt(msh, coil_file, coil_matrix, didt, fn_geo, fn_stl)

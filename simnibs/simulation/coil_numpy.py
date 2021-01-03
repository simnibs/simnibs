# -*- coding: utf-8 -*-\
'''
    Functions for calculating dA/dt fields in mesh nodes from .ccd or nifti coil files
'''

import numpy as np
import nibabel as nib
import os

from .. import __version__
from ..msh import mesh_io
from ..utils.simnibs_logger import logger
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
    import fmm3dpy
    d_position, d_moment = _rotate_coil(ccd_file, coil_matrix)
    # bring everything to SI
    d_position *= 1e-3
    pos = msh.nodes[:] * 1e-3
    A = np.zeros((len(pos), 3), dtype=float)
    out = [
        fmm3dpy.lfmm3d(
            eps=eps,
            sources=d_position.T,
            charges=d_m,
            targets=pos.T,
            pgt=2
        )
        for d_m in d_moment.T
    ]
    A[:, 0] = (out[1].gradtarg[2] - out[2].gradtarg[1])
    A[:, 1] = (out[2].gradtarg[0] - out[0].gradtarg[2])
    A[:, 2] = (out[0].gradtarg[1] - out[1].gradtarg[0])

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
        out[dim] = interpolation.map_coordinates(nifti_image.get_data()[..., dim],
                                                 voxcoords,
                                                 order=1)

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



def set_up_tms_dAdt(msh, coil_file, coil_matrix, didt=1e6, fn_geo=None, fn_stl=None):
    
    add_logo = True # add simnibs logo to coil surface visulization
    
    coil_matrix = np.array(coil_matrix)
    if coil_file.endswith('.ccd'):
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
    elif coil_file.endswith('.nii.gz') or coil_file.endswith('.nii'):
        dadt = _calculate_dadt_nifti(msh, coil_file,
                                     coil_matrix, didt, fn_geo)
    else:
        raise ValueError('coil file must be either a .ccd file or a nifti file')
    
    if fn_stl is not None:
        if not os.path.isfile(fn_stl):
            raise IOError('Could not find stl file: {0}'.format(fn_stl))
        msh_stl = mesh_io.read_stl(fn_stl)
        msh_stl.elm.tag1[:]=0 # tag1 will be shown as gray
        if add_logo:
            msh_stl = _add_logo(msh_stl)
        d_position = np.hstack([msh_stl.nodes[:], np.ones((msh_stl.nodes.nr, 1))])
        msh_stl.nodes = coil_matrix.dot(d_position.T).T[:, :3]
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

   Returns
   ---------
    dAdt: simnibs.mesh_io.ElementData
        dAdt at the elements
    """
    coil_matrix = np.array(coil_matrix, dtype=float)
    return set_up_tms_dAdt(msh, coil_file, coil_matrix, didt, fn_geo, fn_stl)

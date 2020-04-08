# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 19:14:50 2020

@author: axthi
"""

import sys
import numpy as np
import base64
if sys.version_info >= (3, ):
    base64decode = base64.decodebytes
else:
    base64decode = base64.decodestring

from simnibs import mesh_io
from simnibs.segmentation import _marching_cubes_lewiner_luts as mcluts
import simnibs.segmentation._marching_cubes_lewiner_cy as _marching_cubes_lewiner_cy


def marching_cube(volume, affine=None, level=None, step_size=1, only_largest_component=False):
    """ wrapper around marching_cube_lewiner, adapted from from scikit-image
    
    PARAMETERS
    ----------
    volume : float32
        3d image to find isosurfaces. Will internally be converted to float32
        if necessary.
    affine : 4x4 array of float
        affine transformation from voxel indices to mm space
        (default = None)
    level : float
        Contour value to search for isosurfaces in volume. If not given or None, 
        the average of the min and max of vol is used.
        (default = None)
    step_size : int
        Step size in voxels. 
        (default = 1)
    only_largest_component : bool
        Extract and return only largest component
        (default = False)
    
    RETURNS
    -------
    surface : simnibs mesh structure
        the reconstructed isosurface
    EC : float
        euler characteristic of reconstruced isosurface

    NOTES
    ----------        
        * Zeros are padded to the image before passing it to marching_cube_lewiner 
        to ensure that all surfaces are closed.
        * marching_cube_lewiner is called with allow_degenerate=False

    """

    step_size = max(1,step_size)
    
    # To make the marching cubes algorithm close all surfaces, pad all
    # dimensions with one zero
    Y = np.pad(volume, 1, mode="constant")
    
    # Extract the surface 
    vertices, faces, _, _ = marching_cubes_lewiner(Y, level=level, step_size=step_size, 
                                                   allow_degenerate=False)

    # Undo the effect of zero padding on coordinates and transform  
    vertices -= 1
    if affine is not None:
        vertices = apply_affine(vertices, affine)

    surface = mesh_io.Msh(mesh_io.Nodes(vertices), mesh_io.Elements(faces+1))
    
    # extract largest component
    if only_largest_component:
        components = surface.elm.connected_components()
        components.sort(key=len,reverse=True)
        surface = surface.crop_mesh(elements=components[0])
    
    # ensure outwards pointing triangles
    surface.fix_surface_orientation()    
    
    # estimate Euler characteristics: EC = #vertices + #faces - #edges
    M = np.sort(surface.elm.node_number_list[:,0:3], axis=1)
    mesh_edges=np.unique(np.vstack((M[:,[0,1]],M[:,[1,2]],M[:,[0,2]])), axis=0)
    EC=surface.nodes.nr + surface.elm.nr - mesh_edges.shape[0]
    
    return surface, EC


def apply_affine(points, affine):
    """Apply 4-by-4 affine transformation to points in 3D. 
    
    PARAMETERS
    ----------
    points : array_like
        Points to transform. The last axis defines the coordinates in 3D, e.g.,
        points is an n-by-3, n-by-m-by-3, etc. array where all but the last
        dimension is assumed to define points.
    affine : array_like
        The affine transformation to apply. It is right-multiplied to points.
    
    RETURNS
    ----------
    p : array_like
        Transformed points.    
    """
    
    dim = points.shape
    assert dim[-1] == 3
    assert affine.shape == (4,4)
    
    if len(dim) == 1:
        points = points[None,:]
    elif len(dim) > 2:
        points = points.reshape((-1,3))
    
    # append column of ones and apply transformation retaining original array
    # size
    p = np.concatenate((points, np.ones(len(points))[:,np.newaxis]), axis=1)
    p = p.dot(affine.T)[:,:3]
    
    if len(dim) == 1 or len(dim) > 2:
        p = p.reshape(dim)
        
    return p


def marching_cubes_lewiner(volume, level=None, spacing=(1., 1., 1.),
                           gradient_direction='descent', step_size=1,
                           allow_degenerate=True, use_classic=False):
    ''' See the documentation for marching_cubes_lewiner from scikit-image'''
    # Check volume and ensure its in the format that the alg needs
    if not isinstance(volume, np.ndarray) or (volume.ndim != 3):
        raise ValueError('Input volume should be a 3D numpy array.')
    if volume.shape[0] < 2 or volume.shape[1] < 2 or volume.shape[2] < 2:
        raise ValueError("Input array must be at least 2x2x2.")
    volume = np.ascontiguousarray(volume, np.float32)  # no copy if not necessary

    # Check/convert other inputs:
    # level
    if level is None:
        level = 0.5 * (volume.min() + volume.max())
    else:
        level = float(level)
        if level < volume.min() or level > volume.max():
            raise ValueError("Surface level must be within volume data range.")
    # spacing
    if len(spacing) != 3:
        raise ValueError("`spacing` must consist of three floats.")
    # step_size
    step_size = int(step_size)
    if step_size < 1:
        raise ValueError('step_size must be at least one.')
    # use_classic
    use_classic = bool(use_classic)

    # Get LutProvider class (reuse if possible)
    L = _get_mc_luts()

    # Apply algorithm
    func = _marching_cubes_lewiner_cy.marching_cubes
    vertices, faces, normals, values = func(volume, level, L, step_size, use_classic)

    if not len(vertices):
        raise RuntimeError('No surface found at the given iso value.')

    # Output in z-y-x order, as is common in skimage
    vertices = np.fliplr(vertices)
    normals = np.fliplr(normals)

    # Finishing touches to output
    faces.shape = -1, 3
    if gradient_direction == 'descent':
        # MC implementation is right-handed, but gradient_direction is left-handed
        faces = np.fliplr(faces)
    elif not gradient_direction == 'ascent':
        raise ValueError("Incorrect input %s in `gradient_direction`, see "
                         "docstring." % (gradient_direction))
    if spacing != (1, 1, 1):
        vertices = vertices * np.r_[spacing]

    if allow_degenerate:
        return vertices, faces, normals, values
    else:
        fun = _marching_cubes_lewiner_cy.remove_degenerate_faces
        return fun(vertices, faces, normals, values)


def _to_array(args):
    shape, text = args
    byts = base64decode(text.encode('utf-8'))
    ar = np.frombuffer(byts, dtype='int8')
    ar.shape = shape
    return ar


EDGETORELATIVEPOSX = np.array([ [0,1],[1,1],[1,0],[0,0], [0,1],[1,1],[1,0],[0,0], [0,0],[1,1],[1,1],[0,0] ], 'int8')
EDGETORELATIVEPOSY = np.array([ [0,0],[0,1],[1,1],[1,0], [0,0],[0,1],[1,1],[1,0], [0,0],[0,0],[1,1],[1,1] ], 'int8')
EDGETORELATIVEPOSZ = np.array([ [0,0],[0,0],[0,0],[0,0], [1,1],[1,1],[1,1],[1,1], [0,1],[0,1],[0,1],[0,1] ], 'int8')


def _get_mc_luts():
    """ Kind of lazy obtaining of the luts.
    """
    if not hasattr(mcluts, 'THE_LUTS'):
        mcluts.THE_LUTS = _marching_cubes_lewiner_cy.LutProvider(
                EDGETORELATIVEPOSX, EDGETORELATIVEPOSY, EDGETORELATIVEPOSZ,
                _to_array(mcluts.CASESCLASSIC), _to_array(mcluts.CASES),
                _to_array(mcluts.TILING1), _to_array(mcluts.TILING2), _to_array(mcluts.TILING3_1), _to_array(mcluts.TILING3_2),
                _to_array(mcluts.TILING4_1), _to_array(mcluts.TILING4_2), _to_array(mcluts.TILING5), _to_array(mcluts.TILING6_1_1),
                _to_array(mcluts.TILING6_1_2), _to_array(mcluts.TILING6_2), _to_array(mcluts.TILING7_1),
                _to_array(mcluts.TILING7_2), _to_array(mcluts.TILING7_3), _to_array(mcluts.TILING7_4_1),
                _to_array(mcluts.TILING7_4_2), _to_array(mcluts.TILING8), _to_array(mcluts.TILING9),
                _to_array(mcluts.TILING10_1_1), _to_array(mcluts.TILING10_1_1_), _to_array(mcluts.TILING10_1_2),
                _to_array(mcluts.TILING10_2), _to_array(mcluts.TILING10_2_), _to_array(mcluts.TILING11),
                _to_array(mcluts.TILING12_1_1), _to_array(mcluts.TILING12_1_1_), _to_array(mcluts.TILING12_1_2),
                _to_array(mcluts.TILING12_2), _to_array(mcluts.TILING12_2_), _to_array(mcluts.TILING13_1),
                _to_array(mcluts.TILING13_1_), _to_array(mcluts.TILING13_2), _to_array(mcluts.TILING13_2_),
                _to_array(mcluts.TILING13_3), _to_array(mcluts.TILING13_3_), _to_array(mcluts.TILING13_4),
                _to_array(mcluts.TILING13_5_1), _to_array(mcluts.TILING13_5_2), _to_array(mcluts.TILING14),
                _to_array(mcluts.TEST3), _to_array(mcluts.TEST4), _to_array(mcluts.TEST6),
                _to_array(mcluts.TEST7), _to_array(mcluts.TEST10), _to_array(mcluts.TEST12),
                _to_array(mcluts.TEST13), _to_array(mcluts.SUBCONFIG13),
                )

    return mcluts.THE_LUTS


if __name__ == '__main__':
    pass

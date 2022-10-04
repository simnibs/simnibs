# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 18:35:47 2019

@author: axthi
"""

import itertools
import gc
import logging
#from math import log, sqrt
#from multiprocessing import Process
import nibabel as nib
import numba
import numpy as np
import os
#from queue import Queue, Empty
import scipy.sparse
from scipy.spatial import cKDTree, ConvexHull
#from scipy.special import erf
import scipy.ndimage as ndimage
import scipy.ndimage.morphology as mrph
from scipy.ndimage.measurements import label
#from subprocess import Popen, PIPE
import sys
#from threading import Thread
import time

from . import _cat_c_utils
from .marching_cube import marching_cube
from ..mesh_tools import mesh_io
from ..utils import file_finder
from ..utils.simnibs_logger import logger
from ..utils.spawn_process import spawn_process
from ..utils.transformations import resample_vol, crop_vol



# --------------- expansion from central to pial surface ------------------

def expandCS(vertices_org, faces, mm2move_total, ensure_distance=0.2, nsteps=5,
             deform="expand", smooth_mesh=True, skip_lastsmooth=True,
             smooth_mm2move=True, despike_nonmove=True, fix_faceflips=True,
             actualsurf='', debug=False):
    """Deform a mesh by either expanding it in the direction of node normals
    or shinking it in the opposite direction of the normals.

    PARAMETERS
    ----------
    vertices : ndarray
        Vertices of the surface mesh.
    faces : ndarray
        Faces of the surface mesh.
    mm2move : ndarray
        The total distance by which to move the vertices.
    ensure_distance : float
        Minimum distance to ensure in the mesh (default = 0.2).
    nsteps : int
        The number of steps used to reach the desired distance (default = 5).
    deform : {"expand", "shrink"}
        Whether to expand or shink the mesh. Expand corresponds to pushing
        the nodes outwards (in the direction of the normals) whereas shrink
        will pull the nodes inwards (in the direction of minus the normals)
        (default = "expand").
    smooth_mesh : bool
        Smoothing of the mesh by Gaussian smoothing, for each vertex which has
        been moved (default = True).
    skip_lastsmooth : bool
        Taubin smoothing instead of Gaussian smoothing is applied during last step
        to ensure that mm2move is exactly reached (default = True).
    smooth_mm2move : bool
        Smoothing of mm2move. Prevents jigjag-like structures around sulci where
        some of the vertices are not moved anymore to keep ensure_distance (default = True).
    despike_nonmove : bool
        The ray-intersection test that is used to test whether other vertices are
        less than ensure_distance away in the expansion direction gives some false
        positives. Thus, single positives are removed (default = True).
    fix_faceflips = bool
        Changes in the face orientations > 90° indicate that an expansion step
        created surface self intersections. These are resolved by local smoothing
        (default = True).
    actualsurf : string
        string added to the beginning of the logged messages (default = '')
    debug : bool
        Write results from intermediate steps to disk (default = False).

    RETURNS
    ----------
    vertices : ndarray
        Vertices describing the expanded mesh.
    """

    # check inputs
    assert deform in ["expand", "shrink"]
    assert isinstance(nsteps, int)
    assert len(mm2move_total) == len(vertices_org), "The length of mm2move must match that of vertices"

    vertices = vertices_org.copy() # prevent modification of input "vertices"
    move = np.ones(len(vertices), dtype=bool)
    v2f = verts2faces(vertices, faces)

    edges1 = vertices[faces[:, 0]] - vertices[faces[:, 1]]
    edges2 = vertices[faces[:, 0]] - vertices[faces[:, 2]]
    edges3 = vertices[faces[:, 1]] - vertices[faces[:, 2]]
    edges = np.vstack([edges1, edges2, edges3])
    avg_edge_len = np.average(np.linalg.norm(edges, axis=1))
    for i in range(nsteps):
        node_normals = mesh_io.Msh(
            nodes=mesh_io.Nodes(vertices),
            elements=mesh_io.Elements(faces+1)).nodes_normals()[:]
        if deform == "shrink":
            node_normals *= -1

        sys.stdout.flush()
        logger.info(actualsurf+': Iteration '+str(i+1)+' of '+str(nsteps))

        # ---------------------------------
        # update mm2move and vertex normals
        # ---------------------------------
        if i == 0:
            mm2move = mm2move_total/float(nsteps)
        else:
            # distance adjustment is needed to account for
            # small vertex shifts caused by smoothing
            dist = np.sum((vertices-vertices_org) * node_normals, axis=1)
            mm2move = (mm2move_total-dist)/float(nsteps-i)
        mm2move[~move] = 0

        # ------------------------------------------------
        # testing for intersections in the direction of movement
        #
        # testing all nodes against all triangles is slow.
        # thus, find triangles within some distance of each node and test only
        # these for intersections. This reduces # of computations dramatically.
        #
        # This is done twice, one time  against the non-shifted triangles
        # and a second time against a temporarily shifted mesh that approximates
        # the new triangle positions after the movement
        # ------------------------------------------------

        # NOTES:
        # * mm2move is not smoothed here. A temporary copy of mm2move could be
        #   smoothed to make testing more similar to final movement
        # * for the temporarily shifted mesh, one intersection is expected,
        #   as each non-shifted vertex will intersect with its shifted version

        mesh = vertices[faces]
        facenormals_pre = get_triangle_normals(mesh)

        intersect_pairs, _ = segment_triangle_intersect(
            vertices, faces,
            vertices[move] + 1e-4 * avg_edge_len * node_normals[move],
            vertices[move] + (mm2move[move, None] + ensure_distance) * node_normals[move]
        )
        n_intersections = np.bincount(intersect_pairs[:, 0], minlength=len(vertices[move]))

        # create temporary shifted mesh and test again for intersections
        vc_tst = vertices.copy()
        vc_tst[move] += node_normals[move]*mm2move[move, None]
        # We need to shift the nodes to that the ray tracing becomes stable
        vc_tst = smooth_vertices(vc_tst, faces, v2f_map=v2f, mask_move=move, taubin=True)

        intersect_pairs, _ = segment_triangle_intersect(
            vc_tst, faces,
            vertices[move] + 1e-4 * avg_edge_len * node_normals[move],
            vertices[move] + (mm2move[move, None] + ensure_distance) * node_normals[move]
        )
        n_intersections2 = np.bincount(intersect_pairs[:, 0], minlength=len(vertices[move]))
        if debug:
            move_backup = move.copy()

        move[move] = (n_intersections == 0) & (n_intersections2 < 2)
        # -------------------------------------
        # remove isolated "non-move" vertices
        # this is needed as the ray-triangle intersection testing
        # returns a few spurious false positives
        # --------------------------------------
        if despike_nonmove:
            Nnomove = np.zeros(len(move), dtype='uint16')
            Nfaces = np.zeros(len(move), dtype='uint16')
            for j in range(len(move)):
                Nnomove[j] = np.sum(~move[faces[v2f[j]]])
                Nfaces[j] = len(v2f[j])
            # a single vertex reoccurs #faces --> Nnomove>Nfaces will be true
            # when more than one vertex is marked "non-move"
            move = ~(~move & (Nnomove > Nfaces))

        # ----------------------------
        # update the vertex positions, fix flipped faces by local smoothing
        # ----------------------------
        mm2move[~move] = 0
        if smooth_mm2move:
            mm2move = smooth_vertices(mm2move, faces, Niterations=1)
            mm2move[~move] = 0
            pass

        if debug:
            vertices_beforemove = vertices.copy()

        vertices += node_normals*mm2move[:, None]

        # test for flipped surfaces
        mesh = vertices[faces]
        facenormals_post = get_triangle_normals(mesh)
        flipped_faces = np.sum(facenormals_post*facenormals_pre, axis=1) < 0
        if fix_faceflips & np.any(flipped_faces):
            logger.debug(f'{actualsurf}: Fixing {np.sum(flipped_faces)} flipped faces')
            vertices = smooth_vertices(
                vertices, faces,
                verts2consider=np.unique(faces[flipped_faces]),
                v2f_map=v2f, Niterations=5, Ndilate=2)
            mesh = vertices[faces]
            facenormals_post = get_triangle_normals(mesh)
            flipped_faces = np.sum(facenormals_post*facenormals_pre,axis=1) < 0


        if smooth_mesh:
            if skip_lastsmooth & (i == nsteps-1):
                logger.debug(f'{actualsurf}: Last iteration: skipping vertex smoothing')
                vertices = smooth_vertices(
                    vertices, faces, v2f_map=v2f, Niterations=10, mask_move=move, taubin=True)
            else:
                vertices = smooth_vertices(
                    vertices, faces, v2f_map=v2f, mask_move=move)

        logger.info(f'{actualsurf}: Moved {np.sum(move)} of {len(vertices)} vertices.')

        if debug:
            tmpmsh = mesh_io.Msh(nodes=mesh_io.Nodes(vertices),
                         elements=mesh_io.Elements(faces+1))
            filename = "mesh_expand_{:d}_of_{:d}"
            filename = filename.format(i+1, nsteps)
            mesh_io.write_freesurfer_surface(tmpmsh, filename+".fsmesh", ref_fs=True)

            tmpmsh.add_node_field(move, 'move')

            hlpvar = np.zeros(move.shape)
            hlpvar[move_backup] = n_intersections
            tmpmsh.add_node_field(hlpvar, 'n_intersections')
            hlpvar[move_backup] = n_intersections2
            tmpmsh.add_node_field(hlpvar, 'n_intersections2')
            tmpmsh.add_node_field(mm2move_total, 'mm2move_total')

            tmpmsh.elm.add_triangles(faces+tmpmsh.nodes.nr+1,3)
            tmpmsh.nodes.node_coord = np.concatenate((tmpmsh.nodes.node_coord, vertices_beforemove))
            tmpmsh.add_element_field(np.concatenate((flipped_faces,flipped_faces)),'flipped_faces')

            tmpmsh.elm.tag2 = tmpmsh.elm.tag1
            tmpmsh.write(filename+".msh")

    return vertices


def smooth_vertices(vertices, faces, verts2consider=None,
                    v2f_map=None, Niterations=1,
                    Ndilate=0, mask_move=None,
                    taubin=False):
    """Simple mesh smoothing by averaging vertex coordinates or other data
    across neighboring vertices.

    PARAMETERS
    ----------
    vertices : ndarray
        Vertices describing the mesh.
    faces : ndarray
        Faces describing the mesh.
    verts2consider: ndarray
        Array of indices of the vertex that will be smoothed (default: all vertices)
    v2f_map: {list, ndarray}
        Mapping from vertices to faces. Optional (to save a bit of time for repeated use),
        will be created if not given as input.
    Niterations: int
        Number of smoothing iterations (default: 1)
    Ndilate: int
        Number of times the surface region(s) defined by verts2consider are dilated
        before smoothing
    taubin: bool
        Wether to use Taubin smoothing. Defaut:False
    RETURNS
    ----------
    vertices : ndarray
        Vertices describing the expanded mesh.
    """

    if verts2consider is None:
        verts2consider = np.arange(len(vertices))
    if v2f_map is None:
        v2f_map = verts2faces(vertices,faces)

    for i in range(Ndilate):
        f2c = [v2f_map[n] for n in verts2consider]
        f2c, f2cok = list2numpy(f2c,dtype=int)
        f2c = f2c[f2cok]  # faces of verts2consider
        verts2consider = np.unique(faces[f2c])

    if mask_move is not None:
        verts2consider = verts2consider[mask_move[verts2consider]]

    smoo = vertices.copy()
    if taubin:
        m = mesh_io.Msh(nodes=mesh_io.Nodes(smoo),
                        elements=mesh_io.Elements(faces + 1))
        vert_mask = np.zeros(len(vertices), dtype=bool)
        vert_mask[verts2consider] = True
        m.smooth_surfaces_simple(Niterations, nodes_mask=vert_mask)
        smoo = m.nodes[:]
    else:
        for n in verts2consider:
            smoo[n] = np.average(vertices[faces[v2f_map[n]]], axis=(0,1))
        for i in range(Niterations-1):
            smoo2 = smoo.copy()
            for n in verts2consider:
                smoo[n] = np.average(smoo2[faces[v2f_map[n]]], axis=(0,1))
    return smoo



def get_element_neighbors(elements, ntol=1e-6):
    """Get the neighbors of each element in elements by comparing barycenters
    of element faces (e.g., if elements are tetrahedra, the faces are
    triangles).

    PARAMETERS
    ----------
    elements : ndarray
        Array of elements (e.g., triangles or tetrahedra) described as an
        N-by-M, with N being number of elements and M being the number of
        vertices of each element (e.g., 3 and 4 for triangles and tetrahedra,
        respectively).
    ntol : float, optional
        Neighbor tolerance. This parameters controls the upper bound for when
        elements are considered neighbors, i.e. the distance between elements
        has to be smaller than this value (default = 1e-6).

    RETURNS
    ----------
    nearest_neighbors : ndarray
        N-by-M array of indices into elements, i.e. for each face, which is its
        neighboring element.
    ok : ndarray (bool)
        This array tells, for each entry in nearest_neighbors, if this is an
        actual neighbor or not. The nearest neighbors are returned as a numpy
        ndarray of shape elements.shape for ease of interpretation and
        efficiency (and not for example as a list of lists of [possibly]
        unequal lengths), hence this is needed.
    """

    elements_idx = np.arange(len(elements))

    # barycenters of the faces making up each element
    barycenters = np.zeros_like(elements)
    num_nodes_per_el = elements.shape[1]
    for i in range(num_nodes_per_el):
        nodes = np.roll(np.arange(num_nodes_per_el),-i)[:-1] # nodes that make up the ith face
        barycenters[:,i,:] = np.average(elements[:,nodes,:], 1)

    bar_tree = cKDTree(barycenters.reshape(np.multiply(*elements.shape[:-1]),
                                           elements.shape[-1]))
    face_dist, face_idx = bar_tree.query(bar_tree.data, 2)

    nonself = (face_idx != np.arange(len(face_idx))[:,np.newaxis]) # get non-self-references

    # Distance to nearest neighbor. Neighbors having a distance shorter than
    # ntol are considered actual neighbors (i.e. sharing a face)
    face_dist = face_dist[nonself]
    ok = face_dist < ntol
    ok = ok.reshape(elements.shape[:2])

    # Index of nearest neigbor. From the tree search, indices are to nearest
    # element face, however, we wish to find neighboring elements. Hence,
    # reindex.
    face_idx = face_idx[nonself]
    nearest_neighbors = elements_idx.repeat(num_nodes_per_el)[face_idx]
    nearest_neighbors = nearest_neighbors.reshape(elements.shape[:2])

    return nearest_neighbors, ok



def verts2faces(vertices, faces, pad_val=0, array_out_type="list"):
    """Generate a mapping from vertices to faces in a mesh, i.e. for each
    vertices, which elements are it a part of.

    PARAMETERS
    ----------
    vertices : ndarray
        Vertices describing the mesh.
    faces : ndarray
        Faces describing the mesh.
    array_out_type : {"list", "numpy_array"}, optional
        Output type. Numpy arrays will enable vectorized operations to be
        performed on the output, however, in the case of variable number of
        elements per vertice, this will have to be padded

    RETURNS
    ----------
    v2f : {list, ndarray}
        The mapping from vertices to faces.
    ok : ndarray
        Array describing which entries in v2f are actual faces and which are
        "artificial". Since in a mesh, different vertices will often be part of
        different numbers of elements, some rows will have to be padded. This
        array is only returned if array_out_type is set to "numpy_array" since
        sublists of a list can be of variable length.
    """
    # Mapping from node to triangles, i.e. which nodes belongs to which
    # triangles
    v2f = [[] for i in range(len(vertices))]
    for t in range(len(faces)):
        for n in faces[t]:
            v2f[n].append(t)

    if array_out_type == "list":
        return v2f
    elif array_out_type == "numpy_array":
        v2f, ok = list2numpy(v2f, pad_val, int)
        return v2f, ok
    else:
        raise ValueError("Array output type must be list or numpy array.")



def list2numpy(L, pad_val=0, dtype=float):
    """Convert a python list of lists (the sublists being of varying length)
    to a numpy array.

    PARAMETERS
    ----------
    L : list
        The list of lists.
    pad_val : float, int
        The value with which to pad numpy array.
    dtype : datatype, optional
        Datatype of the output array.

    RETURNS
    ----------
    narr : ndarray
        L expressed as a numpy array.
    """

    max_neighbors = len(sorted(L, key=len, reverse=True)[0])
    narr = np.array([r+[np.nan]*(max_neighbors-len(r)) for r in L])
    ok = ~np.isnan(narr)
    narr[~ok] = pad_val
    narr = narr.astype(dtype)

    return narr, ok


def get_triangle_normals(mesh):
    """Get normal vectors for each triangle in the mesh.

    PARAMETERS
    ----------
    mesh : ndarray
        Array describing the surface mesh. The dimension are:
        [# of triangles] x [vertices (of triangle)] x [coordinates (of vertices)].

    RETURNS
    ----------
    tnormals : ndarray
        Normal vectors of each triangle in "mesh".
    """

    tnormals = np.cross(mesh[:,1,:]-mesh[:,0,:],mesh[:,2,:]-mesh[:,0,:]).astype(float)
    tnormals /= np.sqrt(np.sum(tnormals**2,1))[:,np.newaxis]
    return tnormals


def segment_triangle_intersect(vertices, faces, segment_start, segment_end):
    ''' Computes the intersection between a line segment and a triangulated surface

    Parameters
    -----------
    vertices: ndarray
        Array with mesh vertices positions
    faces: ndarray
        Array describing the surface triangles
    segment_start: ndarray
        N_lines x 2 array with the start of the line segments
    segment_end: ndarray
        N_lines x 2 array with the end of the line segments

    Returns
    --------
    indices_pairs: ndarray
        Nx2 array of ints with the pair (segment index, face index) for each intersection
    positions: ndarray
        Nx3 array of floats with the position of the intersections
    '''
    m = mesh_io.Msh(
            nodes=mesh_io.Nodes(vertices),
            elements=mesh_io.Elements(faces+1)
    )
    indices_pairs, positions = m.intersect_segment(segment_start, segment_end)
    # Go from 1-indexed to 0-indexed
    indices_pairs[:, 1] -= 1
    return indices_pairs, positions


def _rasterize_surface(vertices, faces, affine, shape, axis='z'):
    ''' Function to rastherize a given surface given by (vertices, faces) to a volume
    '''
    inv_affine = np.linalg.inv(affine)
    vertices_trafo = inv_affine[:3, :3].dot(vertices.T).T + inv_affine[:3, 3].T

    # switch vertices, dimensions to align with rastherization axis
    if axis == 'z':
        out_shape = shape
    elif axis == 'y':
        vertices_trafo = vertices_trafo[:, [0, 2, 1]]
        out_shape = np.array(shape, dtype=int)[[0, 2, 1]]
    elif axis == 'x':
        vertices_trafo = vertices_trafo[:, [2, 1, 0]]
        out_shape = np.array(shape, dtype=int)[[2, 1, 0]]
    else:
        raise ValueError('"axis" should be x, y, or z')

    grid_points = np.array(
        np.meshgrid(
            *tuple(map(np.arange, out_shape[:2])), indexing="ij"
        )
    ).reshape((2, -1)).T
    grid_points_near = np.hstack([grid_points, np.zeros((len(grid_points), 1))])
    grid_points_far = np.hstack([grid_points, out_shape[2] * np.ones((len(grid_points), 1))])

    # This fixes the search are such that if the volume area to rastherize is smaller
    # than the mesh, we will still trace rays that cross the whole extension of the mesh
    if np.min(vertices_trafo[:, 2]) < 0:
        grid_points_near[:, 2] = 1.1 * np.min(vertices_trafo[:, 2])
    if np.max(vertices_trafo[:, 2]) > out_shape[2]:
        grid_points_far[:, 2] = 1.1 * np.max(vertices_trafo[:, 2])

    # Calculate intersections
    pairs, positions = segment_triangle_intersect(
        vertices_trafo, faces, grid_points_near, grid_points_far
    )

    # Select the intersecting lines
    lines_intersecting, uq_indices, _, counts = np.unique(
        pairs[:, 0], return_index=True, return_inverse=True, return_counts=True
    )

    # The count should never be odd
    if np.any(counts % 2 == 1):
        logger.warning(
            'Found an odd number of crossings! This could be an open surface '
            'or a self-intersection'
        )

    # "z" voxels where intersections occurs
    #inter_z = np.around(positions[:, 2]).astype(int)
    inter_z = (positions[:, 2] + 1).astype(int)
    inter_z[inter_z < 0] = 0
    inter_z[inter_z > out_shape[2]] = out_shape[2]

    # needed to take care of last line
    uq_indices = np.append(uq_indices, [len(pairs)])

    # Go through each point in the grid and assign the z coordinates that are in the mesh
    # (between crossings)
    mask = np.zeros(out_shape, dtype=bool)
    for i, l in enumerate(lines_intersecting):
        # We can do this because we know that the "pairs" variables is ordered with
        # respect to the first variable
        crossings = np.sort(inter_z[uq_indices[i]: uq_indices[i+1]])
        for j in range(0, len(crossings) // 2):
            enter, leave = crossings[2*j], crossings[2*j + 1]
            mask[grid_points[l, 0], grid_points[l, 1], enter:leave] = True

    # Go back to the original frame
    if axis == 'z':
        pass
    elif axis == 'y':
        mask = np.swapaxes(mask, 2, 1)
    elif axis == 'x':
        mask = np.swapaxes(mask, 2, 0)

    return mask


def mask_from_surface(vertices, faces, affine, shape):
    """ Creates a binary mask based on a surface

    Parameters
    ----------
    vertices: ndarray
        Array with mesh vertices positions
    faces: ndarray
        Array describing the surface triangles
    affine: 4x4 ndarray
        Matrix describing the affine transformation between voxel and world coordinates
    shape: 3x1 list
        shape of output mask

    Returns
    ----------
    mask : ndarray of shape 'shape'
       Volume mask
    """

    masks = []

    if len(vertices) == 0 or len(faces) == 0:
        logger.warning("Surface if empty! Return empty volume")
        return np.zeros(shape, dtype=bool)

    # Do the rastherization in 3 directions
    for axis in ['x', 'y', 'z']:
        masks.append(_rasterize_surface(vertices, faces, affine, shape, axis=axis))

    # Return all voxels which are in at least 2 of the masks
    # This is done to reduce spurious results caused by bad tolopogy
    return np.sum(masks, axis=0) >= 2
    #return masks[2]


# --------------- central surface creation ------------------

def createCS(Ymf, Yleft, Ymaskhemis, vox2mm, actualsurf,
             surffolder, fsavgDir, vdist=[1.0, 0.75], voxsize_pbt=[0.5, 0.25],
             voxsize_refineCS=[0.75, 0.5], th_initial=0.714,
             no_selfintersections=True, debug=False):
    """ reconstruction of cortical surfaces based on probalistic label image

    PARAMETERS
    ----------
    Ymf : float32
        image prepared for surface reconstruction
    Yleft : uint8
        binary mask to distinguish between left and right parts of brain
    Ymaskhemis : uint8
        volume with the following labels: lh - 1, rh - 2, lc - 3, rc - 4
    vox2mm : 4x4 array of float
        affine transformation from voxel indices to mm space
    actualsurf : str
        surface to be reconstructed ('lh', 'rh', 'lc' or 'rc')
    surffolder : str
        path for storing results
    fsavgDir : str
        directory containing fsaverage templates

    --> optional parameters:
    vdist : list of float
        final surface resolution (distance between vertices)
        for cerebrum (1. value) and cerebellum (2. value) surfaces
        (default = [1.0, 0.75])
    voxsize_pbt : list of float
        internal voxel size used for pbt calculation and creation of
        initial central surface by marching cubes for cerebrum and
        cerebellum surfaces
        (default=[0.5, 0.25])
    voxsize_refineCS : list of float
            internal voxel size of GM percentage position image used during
            expansion of central surface for cerebrum and cerebellum surfaces
            (default=[0.75, 0.5])
    th_initial : float
            Intensity threshold for the initial central surface. Values closer
            to 1 moves the intial central surface closer to the WM boundary
            (WM-GM boundary corresponds to 1, final CS to 0.5, pial surface to 0)
            (default = 0.714)
    no_selfintersections : bool
        if True:
            1. CAT_DeformSurf is run with the option to avoid selfintersections
            during surface expansion. Helps in particular for the cerebellum
            central surfaces, but is not perfect.
            2. meshfix is used to remove selfintersections after
            surface expansion. CAT_FixTopology can cut away large parts
            of the initial surface. CAT_DeformSurf tries to expand the
            surface again, which is instable and can result in selfintersections
            even when CAT_DeformSurf is run with the option listed in 1.
        (default = True)
    debug : bool
        Write results from intermediate steps to disk (default = False).

    RETURNS
    ----------
        Pcentral : string
            Filename (including path) of the final central surface (gifti)
        Pspherereg : string
            Filename of the spherical registration to the fsaverage template (gifti)
        Pthick : string
            Filename of the cortical thickness (stored as freesurfer curvature file)
        EC : int
            surface Euler characteristics of initial, uncorrected surface (should be 2)
        defect_size : float
            overall size of topology defects of initial, uncorrected surface

    NOTES
    ----------
        This function is adapted from cat_surf_createCS.m of CAT12
        (version 2019-03-22, http://www.neuro.uni-jena.de/cat/).
    """

    if sys.platform == 'win32':
        # Make logging to stderrr more talkative to capture
        # all logging output in case of multiprocessing
        logger.handlers[0].setLevel(logging.DEBUG)

    # add surface name to logger
    formatter_list=[]
    for i in range(len(logger.handlers)):
        formatter_list.append(logger.handlers[i].formatter._fmt)
        formatter = logging.Formatter(f'{actualsurf} '+logger.handlers[i].formatter._fmt)
        logger.handlers[i].setFormatter(formatter)

    logger.info(f'Processing {actualsurf}')

    # ------- crop and upsample subvolume -------
    if 'lh' == actualsurf.lower():
        Ymfs=np.multiply(Ymf,(Ymaskhemis == 1))
        Yside=np.array(Yleft, copy=True)

        # logger.info('INFO INFO lh')
        # logger.debug('DEBUG DEBUG lh')
        # logger.info('2222 INFO INFO lh')
        # time.sleep(1)
        # raise ValueError('error lh: debug')

    elif 'rh' == actualsurf.lower():
        Ymfs=np.multiply(Ymf,(Ymaskhemis == 2))
        Yside=np.logical_not(Yleft)

        # logger.debug('debug RH')
        # time.sleep(2)
        # raise ValueError('error rh: debug')

    elif 'lc' == actualsurf.lower():
        Ymfs=np.multiply(Ymf,(Ymaskhemis == 3))
        Yside=np.array(Yleft, copy=True)

    elif 'rc' == actualsurf.lower():
        Ymfs=np.multiply(Ymf,(Ymaskhemis == 4))
        Yside=np.logical_not(Yleft)

    else:
        logger.error(f'Unknown value of actualsurf: {actualsurf}\\n')
        raise ValueError('error: unknown surface name')

    iscerebellum = (actualsurf.lower() == 'lc') or (actualsurf.lower() == 'rc')
    vdist=vdist[int(iscerebellum)]
    voxsize_pbt=voxsize_pbt[int(iscerebellum)]
    voxsize_refineCS=voxsize_refineCS[int(iscerebellum)]
    logger.debug(f'iscerebellum: {iscerebellum}, vdist: {vdist}, voxsize_pbt: {voxsize_pbt}, voxsize_refineCS: {voxsize_refineCS}')

    if debug:
        tmp = ndimage.uniform_filter(Ymfs, 3)
        Ymfs_filtered = nib.Nifti1Image(tmp, vox2mm)
        fname_ymfsfiltered=os.path.join(surffolder, 'Ymfs_filtered_' + actualsurf + '.nii.gz')
        nib.save(Ymfs_filtered, fname_ymfsfiltered)

    # crop volume
    # NOTE: CAT12 uses the threshold of 1.5 in multiple places as estimate
    # of the position of the middle of gray matter. For our segmentation
    # results, tweaking it to 1.2 works better (based on a single example, Ernie)
    mask = ndimage.uniform_filter(Ymfs, 3) > 1.2 # original threshold for CAT12: > 1.5
    Ymfs, vox2mm_cropped, _ = crop_vol(Ymfs, vox2mm, mask, 4)
    Yside = crop_vol(Yside, vox2mm, mask, 4)[0]

    if debug:
        Ymfs_tmp = nib.Nifti1Image(Ymfs, vox2mm_cropped)
        fname_ymfstest=os.path.join(surffolder,'Ymfs_masked_' + actualsurf + '.nii.gz')
        nib.save(Ymfs_tmp, fname_ymfstest)

    # upsample using linear interpolation (linear is better than cubic for small thicknesses)
    Ymfs, vox2mm_upsampled, _ =resample_vol(np.maximum(1,Ymfs), vox2mm_cropped, voxsize_pbt, order=1, mode='nearest')
    Ymfs=np.minimum(3,np.maximum(1,Ymfs))
    Yside=resample_vol(Yside, vox2mm_cropped, voxsize_pbt, order=1, mode='nearest')[0] > 0.5

    if debug:
        Ymfs_upsampled = nib.Nifti1Image(Ymfs, vox2mm_upsampled)
        fname_ymfstest = os.path.join(surffolder,'Ymfs_test_upsampled_'+actualsurf+'.nii.gz')
        nib.save(Ymfs_upsampled, fname_ymfstest)


    #  -------- pbt calculation and some postprocessing of thickness --------
    #  ----------------  and GM percentage position map ----------------
    stimet = time.time()

    # NOTE: Yth1i is the cortical thickness map
    #       Yppi is the percentage position map: 1 is WM, 0 is GM surface
    Yth1i, Yppi = cat_vol_pbt_AT(Ymfs, voxsize_pbt, actualsurf, debug,
                                 vox2mm_upsampled, surffolder)
    del Ymfs
    gc.collect()

    # post-process THICKNESS map and save to disk
    # the thickness will later be interpolated onto the central surface (in refineCS)
    Yth1i[Yth1i > 10]=0
    Yppi[np.isnan(Yppi)]=0
    I=_cat_c_utils.cat_vbdist(Yth1i,Yside)[1]
    Yth1i=Yth1i.T.flatten()[I-1]
    Yth1t, vox2mm_Yth1t, _ = resample_vol(Yth1i, vox2mm_upsampled, voxsize_refineCS, order=1, mode='nearest')
    Vthk = nib.Nifti1Image(Yth1t, vox2mm_Yth1t)
    fname_thkimg=os.path.join(surffolder,actualsurf+'_thk.nii')
    nib.save(Vthk, fname_thkimg)
    del I, Yside, Yth1i, Yth1t, Vthk
    gc.collect()

    # post-process PERCENTAGE POSITION map
    # Replace isolated voxels and holes in Ypp by its median value
    # indicate isolated holes and replace by median of the neighbors
    Yppi[((Yppi < 0.35)  & np.logical_not(lab(Yppi < 1)))]=1
    Ymsk=(Yppi == 0) & (dilate(Yppi > 0.9,1))
    Yppi=_cat_c_utils.cat_vol_median3(np.float32(Yppi),Ymsk,np.logical_not(Ymsk))

    # indicate isolated objects and replace by median of the neighbors
    Yppi[((Yppi > 0.65) & lab(Yppi == 0))]=0
    Ymsk=((Yppi > 0.95) & dilate((Yppi < 0.1),1))
    Yppi=_cat_c_utils.cat_vol_median3(np.float32(Yppi),Ymsk,np.logical_not(Ymsk))
    del Ymsk
    gc.collect()

    if debug:
        Vppi = nib.Nifti1Image(Yppi, vox2mm_upsampled)
        fname=os.path.join(surffolder,actualsurf+'_Yppi.nii')
        nib.save(Vppi, fname)
        del Vppi
        gc.collect()

    # save to disk
    # this image will later be used by the CAT12 binaries (in refineCS)
    Yppt, vox2mm_Yppt, _  = resample_vol(Yppi, vox2mm_upsampled, voxsize_refineCS, order=1, mode='nearest')
    Yppt[np.isnan(Yppt)]=1
    Yppt[Yppt > 1] = 1
    Yppt[Yppt<0]=0
    Vpp = nib.Nifti1Image(np.uint8(np.rint(Yppt*255)), vox2mm_Yppt)
    Vpp.header['scl_slope'] = 1/255
    Vpp.header['scl_inter'] = 0
    fname_ppimg=os.path.join(surffolder,actualsurf+'_pp.nii')
    nib.save(Vpp, fname_ppimg)
    del Yppt, Vpp
    gc.collect()
    logger.info(f'Thickness estimation ({"{0:.2f}".format(voxsize_pbt)} mm{chr(179)})...'+time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))


    # ---------- generation of initial central surface (using Yppi) -------
    stimet = time.time()
    logger.info('Calling marching cubes')
    CS, EC = marching_cube(Yppi, affine=vox2mm_upsampled, level=th_initial,
                           step_size=round(vdist/voxsize_pbt), only_largest_component=True, n_uniform=2)

    Praw=os.path.join(surffolder,actualsurf+'.central.nofix.gii')
    mesh_io.write_gifti_surface(CS, Praw)
    del Yppi, CS
    gc.collect()
    logger.info(f'Create initial surface: '+time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))


    # -------- refine intial surface and register to fsaverage template --------
    Pcentral, Pspherereg, Pthick, defect_size = refineCS(Praw, fname_thkimg, fname_ppimg,
                                                         fsavgDir, vdist, no_selfintersections, debug)
    logger.info(f'Surface Euler number: {EC}')
    logger.info(f'Overall size of topology defects: {defect_size}')


    # -------- remove temporary files --------
    if not debug:
        if os.path.isfile(Praw):
            os.remove(Praw)
        if os.path.isfile(fname_ppimg):
            os.remove(fname_ppimg)
        if os.path.isfile(fname_thkimg):
            os.remove(fname_thkimg)

    # restore logger
    for i in range(len(logger.handlers)):
        formatter = logging.Formatter(formatter_list[i])
        logger.handlers[i].setFormatter(formatter)

    if sys.platform == 'win32':
        # Make logging to stderrr less talkative again
        logger.handlers[0].setLevel(logging.INFO)

    return Pcentral, Pspherereg, Pthick, EC, defect_size


def cat_vol_pbt_AT(Ymf, resV, actualsurf, debug=False, vox2mm=None, surffolder=None):
    """ Estimate cortical thickness and surface position using pbt2x

    PARAMETERS
    ----------
       Ymf    : tissue segment image or better the noise, bias, and
                intensity corrected
       resV   : voxel resolution (only isotropic)
       actualsurf: the side being processed (lh,rh,rc,lc)

       --> optional parameters:
       debug  : bool
                (default = False)

       vox2mm : image-to-world transformation for debugging purposes

    RETURNS
    ----------
       Ygmt   : GM thickness map
       Ypp    : percentage position map

    NOTES
    ----------
       This function is adapted from cat_vol_pbt.m of CAT12
       (version 2019-03-22, http://www.neuro.uni-jena.de/cat/).

       This python version fixed a side effect caused by unintended
       modification of variables. The problem is due to the C language
       does not use the passby-value scheme as Matlab when passing arrays,
       so an unintended modification to a dummy array in the C function
       can cause side effects. See the parameters Ywmd, Ycsfd, Ygmt, Ygmt1,
       Ygmt2 in lines 86, 130, 134, 156, 157, 172, 179, 182, 189 in the matlab
       function cat_vol_pbt_AT.m when calling the C function
       "cat_vol_localstat.c".

       The python version also fixed a bug called array index out of bound.
       The index used to address array items in the C function
       "cat_vol_pbtp.cpp" exceeds the allowed value by 1. It causes
       undefined behavior in the C function "cat_vol_pbtp.cpp". See lines 52
       and 63 in the updated file "cat_vol_pbtp.cpp".

     Reference
     ----------
       Dahnke, R; Yotter R; Gaser C.
       Cortical thickness and central surface estimation.
       NeuroImage 65 (2013) 226-248.

    """
    debug=int(debug)

    if (np.sum(np.round(np.asanyarray(Ymf).reshape(-1, 1)) == np.asanyarray(Ymf).reshape(-1, 1)) / np.asarray(Ymf).size) > 0.9:
        binary = True
    else:
        binary = False

    minfdist = 2
    # NOTE: CAT12 uses the threshold of 1.5 in multiple places as estimate
    # of the position of the middle of gray matter. For our segmentation
    # results, tweaking it to 1.2 works better (based on a single example, Ernie)
    thres_magic = 1.2 # CAT12's original magical threshold is 1.5

    #  WM distance
    #  Estimate WM distance Ywmd and the outer CSF distance Ycsfdc to correct
    #  the values in CSF area are to limit the Ywmd to the maximum value that
    #  is possible within the cortex.
    #  The increasement of this area allow a more accurate and robust projection.
    #  cat_vol_eidist used speed map to align voxel to the closer gyrus
    #  that is not required for the correction map.

    #  RD 201803:
    #  The speed map weighting "max(0.5,min(1,Ymf/2))" is not strong enough to
    #  support asymmetric structures. The map "max(eps,min(1,((Ymf-1)/1.1).^4))"
    #  works much better but it leads to much higher thickness results (eg. in
    #  the Insula).

    stimet = time.time()
    stimet2 = stimet

    YMM = np.logical_or(erosion(Ymf < thres_magic, 1), np.isnan(Ymf))
    F = np.fmax(0.5, np.fmin(1, Ymf / 2))
    YM = np.fmax(0, np.fmin(1, (Ymf - 2)))
    # In CAT12 this extends quite far into gray matter, as there are
    # values that are > 2 quite close to CSF. In our case the values drop
    # off faster, so could change the 2 into 1.9
    # YM = np.fmax(0, np.fmin(1, (Ymf - 1.9))) --> didn't do much
    YM[YMM] = np.nan
    Ywmd = _cat_c_utils.cat_vol_eidist(
            YM, F, np.array([1, 1, 1]), 1, 1, 0, debug)[0]

    ## DEBUG
    if debug:
        F_image = nib.Nifti1Image(F, vox2mm)
        fname_F=os.path.join(surffolder,'F_' + actualsurf + '.nii.gz')
        nib.save(F_image, fname_F)

        YM_image = nib.Nifti1Image(YM, vox2mm)
        fname_YM=os.path.join(surffolder,'YM_' + actualsurf + '.nii.gz')
        nib.save(YM_image, fname_YM)

        Ywmd_image = nib.Nifti1Image(Ywmd, vox2mm)
        fname_Ywmd=os.path.join(surffolder,'Ywmd_' + actualsurf + '.nii.gz')
        nib.save(Ywmd_image, fname_Ywmd)

    F = np.fmax(1.0, np.fmin(1, Ymf / 2))
    YM = np.fmax(0, np.fmin(1, (Ymf - 1)))
    YM[YMM] = np.nan
    Ycsfdc = _cat_c_utils.cat_vol_eidist(
        YM, F, np.array([1, 1, 1]), 1, 1, 0, debug)[0]
    del F, YMM
    gc.collect()

    ## DEBUG
    if debug:
        Ycsfdc_image = nib.Nifti1Image(Ycsfdc, vox2mm)
        fname_Ycsfdc=os.path.join(surffolder,'Ycsfdc_' + actualsurf + '.nii.gz')
        nib.save(Ycsfdc_image, fname_Ycsfdc)

    if not binary:
        # limit the distance values outside the GM/CSF boudary to the distance possible in the GM
        notnan = ~np.isnan(Ywmd)
        YM = np.full(Ywmd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ywmd[notnan] > minfdist), (Ymf[notnan] <= thres_magic))

        # It might happen that there are inf-infs here which triggers a
        # runtime warning. inf-inf produces a nan and those seem to be
        # masked out later on, so I'll ignore the warning here
        with np.errstate(invalid='ignore'):
            Ywmd[YM] = Ywmd[YM] - Ycsfdc[YM]
        Ywmd[np.isinf(Ywmd)] = 0
        del Ycsfdc
        gc.collect()

        # smoothing of distance values inside the GM
        notnan = ~np.isnan(Ywmd)
        YM = np.full(Ywmd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ywmd[notnan] > minfdist), (Ymf[notnan] > thres_magic))
        YwmdM = np.array(Ywmd, copy=True)
        YwmdM = _cat_c_utils.cat_vol_localstat(YwmdM, YM, 1, 1)[0]
        Ywmd[YM] = YwmdM[YM]

        # smoothing of distance values outside the GM
        notnan = ~np.isnan(Ywmd)
        YM = np.full(Ywmd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ywmd[notnan] > minfdist), (Ymf[notnan] <= thres_magic))
        YwmdM = np.array(Ywmd, copy=True)
        for i in np.arange(1, 3):
            YwmdM = _cat_c_utils.cat_vol_localstat(YwmdM, YM, 1, 1)[0]
        Ywmd[YM] = YwmdM[YM]

        # reducing outliers in the GM/CSF area
        notnan = ~np.isnan(Ywmd)
        YM = np.full(Ywmd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ywmd[notnan] > minfdist), (Ymf[notnan] < 2.0))
        YwmdM = np.array(Ywmd, copy=True)
        YwmdM = _cat_c_utils.cat_vol_median3(YwmdM, YM, YM)
        Ywmd[YM] = YwmdM[YM]
        del YwmdM, YM
        gc.collect()

    ## DEBUG
    if debug:
        Ywmd_image = nib.Nifti1Image(Ywmd, vox2mm)
        fname_Ywmd=os.path.join(surffolder,'Ywm_after_a_million_steps_' + actualsurf + '.nii.gz')
        nib.save(Ywmd_image, fname_Ywmd)

    logger.info(f'WM distance: ' +
                time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))

    #  CSF distance
    #  Similar to the WM distance, but keep in mind that this map is
    #  incorrect in blurred sulci that is handled by PBT
    stimet = time.time()
    YMM = np.any((erosion(Ymf < thres_magic, 1), erosion(
        Ymf > 2.5, 1), np.isnan(Ymf)), axis=0)
    F = np.fmax(0.5, np.fmin(1, (4 - Ymf) / 2))

    tmp_param10 = 1.6 # in CAT12 originally set to 2
    YM = np.fmax(0, np.fmin(1, (tmp_param10 - Ymf)))
    YM[YMM] = np.nan
    Ycsfd = _cat_c_utils.cat_vol_eidist(
        YM, F, np.array([1, 1, 1]), 1, 1, 0, debug)[0]
    F = np.fmax(1, np.fmin(1, (4 - Ymf) / 2))

    YM = np.fmax(0, np.fmin(1, (3 - Ymf)))
    YM[YMM] = np.nan
    Ywmdc = _cat_c_utils.cat_vol_eidist(
        YM, F, np.array([1, 1, 1]), 1, 1, 0, debug)[0]

    YM = np.fmax(0, np.fmin(1, (2.7 - Ymf)))
    YM[YMM] = np.nan
    Ywmdx = _cat_c_utils.cat_vol_eidist(
            YM, F, np.array([1, 1, 1]), 1, 1, 0, debug)[0] + 0.3
    del F, YMM
    gc.collect()
    Ywmdc = np.fmin(Ywmdc, Ywmdx)

    ## DEBUG
    if debug:
        Ywmdx_image = nib.Nifti1Image(Ywmdx, vox2mm)
        fname_Ywmdx=os.path.join(surffolder,'Ywmdx_'+actualsurf+'.nii.gz')
        nib.save(Ywmdx_image, fname_Ywmdx)

        Ywmdc_image = nib.Nifti1Image(Ywmdc, vox2mm)
        fname_Ywmdc=os.path.join(surffolder,'Ywmdc_'+ actualsurf+ '.nii.gz')
        nib.save(Ywmdc_image, fname_Ywmdc)

        Ycsfdc_image = nib.Nifti1Image(Ycsfd, vox2mm)
        fname_Ycsfdc=os.path.join(surffolder,'Ycsfdc_before_a_million_steps_' + actualsurf + '.nii.gz')
        nib.save(Ycsfdc_image, fname_Ycsfdc)

    if not binary:
        notnan = ~np.isnan(Ycsfd)
        YM = np.full(Ycsfd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ycsfd[notnan] > minfdist), (Ymf[notnan] >= 2.5))
        Ycsfd[YM] = Ycsfd[YM] - Ywmdc[YM]
        Ycsfd[np.isinf(- Ycsfd)] = 0
        del Ywmdc
        gc.collect()
        notnan = ~np.isnan(Ycsfd)
        YM = np.full(Ycsfd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ycsfd[notnan] > minfdist), (Ymf[notnan] < 2.5))
        YcsfdM = np.array(Ycsfd, copy=True)
        YcsfdM = _cat_c_utils.cat_vol_localstat(YcsfdM, YM, 1, 1)[0]
        Ycsfd[YM] = YcsfdM[YM]
        notnan = ~np.isnan(Ycsfd)
        YM = np.full(Ycsfd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ycsfd[notnan] > minfdist), (Ymf[notnan] >= 2.5))
        YcsfdM = np.array(Ycsfd, copy=True)
        for i in np.arange(1, 3):
            YcsfdM = _cat_c_utils.cat_vol_localstat(YcsfdM, YM, 1, 1)[0]
        Ycsfd[YM] = YcsfdM[YM]
        notnan = ~np.isnan(Ycsfd)
        YM = np.full(Ycsfd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ycsfd[notnan] > minfdist), (Ymf[notnan] > 2.0))
        YcsfdM = np.array(Ycsfd, copy=True)
        YcsfdM = _cat_c_utils.cat_vol_median3(YcsfdM, YM, YM)
        Ycsfd[YM] = YcsfdM[YM]
        del YcsfdM, YM
        gc.collect()

    ## DEBUG
    if debug:
        Ycsfdc_image = nib.Nifti1Image(Ycsfd, vox2mm)
        fname_Ycsfdc=os.path.join(surffolder,'Ycsfdc_a_million_steps_' + actualsurf + '.nii.gz')
        nib.save(Ycsfdc_image, fname_Ycsfdc)

    logger.info(f'CSF distance: ' +
                time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))

    # PBT thickness mapping using pbt2x
    # --------------------
    stimet = time.time()

    # add 1 to keep the iteration number the same as matlab
    iterator = 1 / np.mean((np.squeeze(resV)).flatten()) + 1

    # Estimation of the cortical thickness with sulcus (Ygmt1) and gyri
    # correction (Ygmt2) to create the final thickness as the minimum map
    # of both.

    # estimate thickness with PBT approach
    YcsfdM = np.array(Ycsfd, copy=True)
    Ygmt1 = _cat_c_utils.cat_vol_pbtp(Ymf, Ywmd, YcsfdM)[0]
    YwmdM = np.array(Ywmd, copy=True)
    Ygmt2 = _cat_c_utils.cat_vol_pbtp(4 - Ymf, Ycsfd, YwmdM)[0]
    del YcsfdM, YwmdM
    gc.collect()

    # avoid meninges !
    Ygmt1 = np.fmin(Ygmt1, Ycsfd + Ywmd)
    Ygmt2 = np.fmin(Ygmt2, Ycsfd + Ywmd)

    # median filter to remove outliers
    notnan = ~np.isnan(Ygmt1)
    YM = np.full(Ygmt1.shape, False, dtype=bool)
    YM[notnan] = Ygmt1[notnan] > 0
    Ygmt1 = _cat_c_utils.cat_vol_median3(Ygmt1, YM, YM)

    notnan = ~np.isnan(Ygmt2)
    YM = np.full(Ygmt2.shape, False, dtype=bool)
    YM[notnan] = Ygmt2[notnan] > 0
    Ygmt2 = _cat_c_utils.cat_vol_median3(Ygmt2, YM, YM)

    ## DEBUG
    if debug:
        Ygmt1_image = nib.Nifti1Image(Ygmt1, vox2mm)
        fname_Ygmt1=os.path.join(surffolder,'Ygmt1_' + actualsurf + '.nii.gz')
        nib.save(Ygmt1_image, fname_Ygmt1)
        Ygmt2_image = nib.Nifti1Image(Ygmt2, vox2mm)
        fname_Ygmt2=os.path.join(surffolder,'Ygmt2_' + actualsurf + '.nii.gz')
        nib.save(Ygmt2_image, fname_Ygmt2)

    # estimation of Ypp for further GM filtering without sulcul blurring
    Ygmt = np.fmin(Ygmt1, Ygmt2)
    YM = np.logical_and((Ymf >= thres_magic), (Ymf < 2.5))
    Ypp = np.zeros(Ymf.shape, dtype=np.float32)
    Ypp[Ymf >= 2.5] = 1
    eps = np.finfo(float).eps
    Ypp[YM] = np.fmin(Ycsfd[YM], Ygmt[YM] - Ywmd[YM]) / (Ygmt[YM] + eps)
    Ypp[Ypp > 2] = 0
    notnan = ~np.logical_or(np.isnan(Ywmd), np.isnan(Ygmt))
    YM = np.full(Ywmd.shape, False, dtype=bool)
    YM[notnan] = np.squeeze((Ygmt[notnan] <= resV) & (
        Ywmd[notnan] <= resV) & (Ygmt[notnan] > 0))
    Ypp[YM] = (Ymf[YM] - 1) / 2
    Ygmts = np.array(Ygmt, copy=True)
    for i in np.arange(1, iterator):
        Ygmts = _cat_c_utils.cat_vol_localstat(Ygmts, Ygmt1 > 0, 1, 1)[0]

    Ygmt[Ygmts > 0] = Ygmts[Ygmts > 0]
    Ygmts = np.array(Ygmt1, copy=True)
    for i in np.arange(1, iterator):
        Ygmts = _cat_c_utils.cat_vol_localstat(Ygmts, (((Ygmt > 1) | (Ypp > 0.1)) & (
            Ygmt > 0) & ((Ygmt > 1) | (Ymf > 1.8))), 1, 1)[0]

    Ygmt1[Ygmts > 0] = Ygmts[Ygmts > 0]
    Ygmts = np.array(Ygmt2, copy=True)
    for i in np.arange(1, iterator):
        Ygmts = _cat_c_utils.cat_vol_localstat(Ygmts, (((Ygmt > 1) | (Ypp > 0.1)) & (
            Ygmt > 0) & ((Ygmt > 1) | (Ymf > 1.8))), 1, 1)[0]

    Ygmt2[Ygmts > 0] = Ygmts[Ygmts > 0]

    # mix result
    # only minimum possible, because Ygmt2 is incorrect in blurred sulci
    Ygmt = np.fmin(Ygmt1, Ygmt2)

    ## DEBUG
    if debug:
        Ygmt_image = nib.Nifti1Image(Ygmt, vox2mm)
        fname_Ygmt=os.path.join(surffolder,'Ygmt_' + actualsurf + '.nii.gz')
        nib.save(Ygmt_image, fname_Ygmt)

    Ygmts = np.array(Ygmt, copy=True)
    for i in np.arange(1, iterator):
        Ygmts = _cat_c_utils.cat_vol_localstat(Ygmts, (((Ygmt > 1) | (Ypp > 0.1)) & (
            Ygmts > 0) & ((Ygmt > 1) | (Ymf > 1.8))), 1, 1)[0]

    Ygmt[Ygmts > 0] = Ygmts[Ygmts > 0]

    # Estimation of a mixed percentual possion map Ypp.
    YM = ((Ymf >= thres_magic) & (Ymf < 2.5) & (Ygmt > eps))
    Ycsfdc = np.array(Ycsfd, copy=True)
    Ycsfdc[YM] = np.fmin(Ycsfd[YM], Ygmt[YM] - Ywmd[YM])
    Ypp = np.zeros(Ymf.shape, dtype=np.float32)
    Ypp[Ymf >= 2.5] = 1
    Ypp[YM] = Ycsfdc[YM] / (Ygmt[YM] + eps)
    Ypp[Ypp > 2] = 0
    notnan = ~np.logical_or(np.isnan(Ywmd), np.isnan(Ygmt))
    YM = np.full(Ywmd.shape, False, dtype=bool)
    YM[notnan] = np.squeeze((Ygmt[notnan] <= resV) & (
        Ywmd[notnan] <= resV) & (Ygmt[notnan] > 0))

    Ypp[YM] = (Ymf[YM] - 1) / 2 - 0.2
    Ypp[np.isnan(Ypp)] = 0
    Ypp[Ypp < 0] = 0

    ## DEBUG
    if debug:
        Ypp_image = nib.Nifti1Image(Ypp, vox2mm)
        fname_Ypp=os.path.join(surffolder,'Ypp_pbt_' + actualsurf + '.nii.gz')
        nib.save(Ypp_image, fname_Ypp)

    # Final corrections for thickness map with thickness limit of 10 mm.
    # Resolution correction of the thickness map after all other operations,
    # because PBT actually works only with the voxel-distance (isotropic 1 mm)
    Ygmt = Ygmt * resV
    Ygmt[Ygmt > 10] = 10

    logger.info(f'PBT2x thickness: ' +
                time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))

    logger.info(f'Cortical thickness and surface position estimation: {debug}: ' + time.strftime(
                '%H:%M:%S', time.gmtime(time.time() - stimet2)))

    return Ygmt, Ypp


def refineCS(Praw, fname_thkimg, fname_ppimg, fsavgDir, vdist=1.0,
        no_selfintersections=True, debug=False):
    """ wrapper around the CAT12 binaries to refine the initial central surface
        and register it to the fsaverage template.

    Roughly, it goes through the following steps:
        * detection and removal of topological defects
        * surface expansion until 0.5 in the perc. position image is reached
        * interpolation of the cortical thickness image at the surface nodes
        * spherical registration to the fsaverage template

    PARAMETERS
    ----------
    Praw : string
        Filename (including path) of the initial central surface created by
        marching cubes. File has to be in gifti format. The letters until
        the first '.' have to indicate the surface name.
        (e.g. path_to_file/lh.rest_of_name)
    fname_thkimg : string
        Filename of cortical thickness image (nifti)
    fname_ppimg : string
        Filename of percentage position image
    fsavgDir : string
        Directory containing the fsaverage template surfaces
    vdist : float
        Desired vertex distance. Smaller numbers give denser meshes.
        (default = 1.0)
    no_selfintersections : bool
        if True:
            1. CAT_DeformSurf is run with the option to avoid selfintersections
            during surface expansion. Helps in particular for the cerebellum
            central surfaces, but is not perfect.
            2. meshfix is used to remove selfintersections after
            surface expansion. CAT_FixTopology can cut away large parts
            of the initial surface. CAT_DeformSurf tries to expand the
            surface again, which is instable and can result in selfintersections
            even when CAT_DeformSurf is run with the option listed in 1.
        (default = True)
    debug : bool
        Save intermediate results for debugging
        (default = False)

    RETURNS
    -------
    Pcentral : string
        Filename (including path) of the final central surface (gifti)
    Pspherereg : string
        Filename of the spherical registration to the fsaverage template (gifti)
    Pthick : string
        Filename of the cortical thickness (stored as freesurfer curvature file)
    defect_sizeOut : float
        total size of areas with topological defects

    """

    # ---------------- get surface filenames ----------------
    [surffolder,actualsurf]=os.path.split(Praw)
    actualsurf=actualsurf.split('.',1)[0]

    Pcentral=os.path.join(surffolder,actualsurf+'.central.gii')

    Pthick=os.path.join(surffolder,actualsurf+'.thickness')
    Pdefects0=os.path.join(surffolder,actualsurf+'.defects')

    Psphere0=os.path.join(surffolder,actualsurf+'.sphere.nofix.gii')
    Psphere=os.path.join(surffolder,actualsurf+'.sphere.gii')
    Pspherereg=os.path.join(surffolder,actualsurf+'.sphere.reg.gii')

    Pfsavg=os.path.join(fsavgDir,actualsurf+'.central.freesurfer.gii')
    Pfsavgsph=os.path.join(fsavgDir,actualsurf+'.sphere.freesurfer.gii')

    if debug:
        Pdebug=os.path.join(surffolder,actualsurf+'.debug.msh')
        # contains:
        # region 1: initial surface with defects (as node data)
        # region 2: spherical version of initial surface with defects (as node data)
        # region 3: surface after topology correction
        # region 4: final surface with thickness and perc. positions (as node data)


    # ------- mark topological defects --------
    stimet = time.time()

    # spherical surface mapping 1 of the uncorrected surface for topology correction
    cmd = [file_finder.path2bin("CAT_Surf2Sphere"), Praw, Psphere0, '5']
    spawn_process(cmd)

    # estimate size of topology defects (in relation to number of vertices and mean brain with 100000 vertices)
    cmd = [file_finder.path2bin("CAT_MarkDefects"), Praw, Psphere0, Pdefects0]
    spawn_process(cmd)

    defect_sizes=mesh_io.read_curv(Pdefects0)
    defect_sizeOut=np.rint( 100000 * np.sum(defect_sizes.flatten() > 0) / defect_sizes.flatten().size)

    if debug:
        # add initial surface and spherical version with defects to .msh for later viewing in gmsh
        CS_dbg = mesh_io.read_gifti_surface(Praw)
        CS = mesh_io.read_gifti_surface(Psphere0)
        CS_dbg.elm.add_triangles(CS.elm.node_number_list[:,0:3]+CS_dbg.nodes.nr,2)
        CS_dbg.nodes.node_coord = np.concatenate((CS_dbg.nodes.node_coord, CS.nodes.node_coord))
        CS_dbg.add_node_field(np.hstack((defect_sizes,defect_sizes)),'defects on surfs 1 and 2')
        del CS
    del defect_sizes
    gc.collect()

    if os.path.isfile(Pdefects0):
        os.remove(Pdefects0)

    logger.info(f'Preparing surface improvment: '+time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))


    # --------- topology correction ---------
    stimet = time.time()

    cmd=[file_finder.path2bin("CAT_FixTopology"), '-lim', '128', '-bw', '512',
         '-n', '81920', '-refine_length', "{:.2f}".format(2 * vdist),
         Praw, Psphere0, Pcentral]
    spawn_process(cmd)

    if debug:
        # add surface after topology correction to .msh for later viewing in gmsh
        CS = mesh_io.read_gifti_surface(Pcentral)
        CS_dbg.elm.add_triangles(CS.elm.node_number_list[:,0:3]+CS_dbg.nodes.nr,3)
        CS_dbg.nodes.node_coord = np.concatenate((CS_dbg.nodes.node_coord, CS.nodes.node_coord))
        del CS
        gc.collect()

    logger.info(f'Topology correction: '+time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))


    # --------- surface refinement by deformation based on the PP map ---------
    stimet = time.time()

    if no_selfintersections:
        force_no_selfintersections = '1'
    else:
        force_no_selfintersections = '0'

    #cmd=f'\"{file_finder.path2bin("CAT_DeformSurf")}\" \"{fname_ppimg}\" none 0 0 0 \"{Pcentral}\" \"{Pcentral}\" none 0 1 -1 .1 avg -0.1 0.1 .2 .1 5 0 0.5 0.5 n 0 0 0 150 0.01 0.0 {force_no_selfintersections}'
    cmd=[file_finder.path2bin("CAT_DeformSurf"), fname_ppimg, 'none', '0', '0', '0',
         Pcentral, Pcentral, 'none', '0', '1', '-1', '.1',
         'avg', '-0.1', '0.1', '.2', '.1', '5', '0', '0.5', '0.5',
         'n', '0', '0', '0', '150', '0.01', '0.0', force_no_selfintersections]
    spawn_process(cmd)

    if no_selfintersections:
        # remove self-intersections using meshfix
        CS = mesh_io.read_gifti_surface(Pcentral)
        mesh_io.write_off(CS, Pcentral+'.off')

        cmd=[file_finder.path2bin("meshfix"), Pcentral+'.off', '-o', Pcentral+'.off']
        spawn_process(cmd)

        CS = mesh_io.read_off(Pcentral+'.off')
        mesh_io.write_gifti_surface(CS,Pcentral)
        if os.path.isfile(Pcentral+'.off'):
            os.remove(Pcentral+'.off')
        if os.path.isfile('meshfix_log.txt'):
            os.remove('meshfix_log.txt')
        del CS
        gc.collect()

    # need some more refinement because some vertices are distorted after CAT_DeformSurf
    cmd=[file_finder.path2bin("CAT_RefineMesh"), Pcentral, Pcentral, "{:.2f}".format(1.5 * vdist ), '0']
    spawn_process(cmd)

    cmd=[file_finder.path2bin("CAT_DeformSurf"), fname_ppimg, 'none', '0', '0', '0',
          Pcentral, Pcentral, 'none', '0', '1', '-1', '.2',
          'avg', '-0.05', '0.05', '.1', '.1', '5', '0', '0.5', '0.5',
          'n', '0', '0', '0', '50', '0.01', '0.0', force_no_selfintersections]
    spawn_process(cmd)

    # map thickness data on final surface
    CS = mesh_io.read_gifti_surface(Pcentral)
    Vthk=nib.load(fname_thkimg)
    nd = mesh_io.NodeData.from_data_grid(CS, Vthk.get_data(), Vthk.affine, 'thickness')
    mesh_io.write_curv(Pthick, nd.value, nd.nr)

    if debug:
        # add prefinal surface with thickness and pp data to .msh
        thickness=np.hstack((np.zeros_like(CS_dbg.nodes.node_number),nd.value))
        # Yppt sampled on the final surface should be distributed sharply around 0.5
        Vpp = nib.load(fname_ppimg)
        nd = mesh_io.NodeData.from_data_grid(CS, Vpp.get_data(), Vpp.affine, 'pp')
        pponsurf=np.hstack((np.zeros_like(CS_dbg.nodes.node_number),nd.value))

        CS_dbg.elm.add_triangles(CS.elm.node_number_list[:,0:3]+CS_dbg.nodes.nr,4)
        CS_dbg.nodes.node_coord = np.concatenate((CS_dbg.nodes.node_coord, CS.nodes.node_coord))
        CS_dbg.add_node_field(thickness,'thickness on surf 4')
        CS_dbg.add_node_field(pponsurf,'perc. position on surf 4')
        del Vpp, thickness, pponsurf
    del CS, Vthk, nd
    gc.collect()

    logger.info(f'Refine central surface: '+time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))


    # AT: this part can create self intersections when neighboring surfaces are close to each other
    # which lead to artifacts during expansion to pial surfaces.
    # # ---- final correction of central surface in highly folded areas --------
    # # ----------------  with high mean curvature --------------
    # stimet = time.time()

    # cmd = [file_finder.path2bin("CAT_Central2Pial"), '-equivolume',
    #         '-weight', '0.3', Pcentral, Pthick, Pcentral, '0']
    # spawn_process(cmd)

    # if debug:
    #     # add final central surface to .msh and save .msh
    #     CS = mesh_io.read_gifti_surface(Pcentral)
    #     CS_dbg.elm.add_triangles(CS.elm.node_number_list[:,0:3]+CS_dbg.nodes.nr,5)
    #     CS_dbg.nodes.node_coord = np.concatenate((CS_dbg.nodes.node_coord, CS.nodes.node_coord))
    #     mesh_io.write_msh(CS_dbg,Pdebug)
    #     del CS
    #     gc.collect()

    # logger.info(f'Correction of central surface in highly folded areas 2: '+time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))


    # -------- registration to FSAVERAGE template --------
    stimet = time.time()

    # spherical surface mapping 2 of corrected surface
    cmd = [file_finder.path2bin("CAT_Surf2Sphere"), Pcentral, Psphere, '10']
    spawn_process(cmd)

    # spherical registration to fsaverage template
    cmd = [file_finder.path2bin("CAT_WarpSurf"), '-steps', '2', '-avg',
           '-i', Pcentral, '-is', Psphere, '-t', Pfsavg, '-ts', Pfsavgsph, '-ws', Pspherereg]
    spawn_process(cmd)

    logger.info(f'Registration to FSAVERAGE template: '+time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))


    # -------- remove unnecessary files --------
    if not debug:
        if os.path.isfile(Psphere0):
            os.remove(Psphere0)

        if os.path.isfile(Pdefects0):
            os.remove(Pdefects0)

        # if os.path.isfile(Psphere):
        #     os.remove(Psphere)

    return Pcentral, Pspherereg, Pthick, defect_sizeOut


def dilate(image,n):
    nan_inds = np.isnan(image)
    image[nan_inds] = 0
    image = image > 0.5
    se = np.ones((2*n+1,2*n+1,2*n+1),dtype=bool)
    return mrph.binary_dilation(image,se)>0


def erosion(image,n):
    nan_inds = np.isnan(image)
    image[nan_inds] = 0
    image = image > 0.5
    return ~dilate(~image,n)


def lab(image):
    labels, _ = label(image)
    return (labels == np.argmax(np.bincount(labels.flat)[1:])+1)


def close(image,n):
    nan_inds = np.isnan(image)
    image[nan_inds] = 0
    image = image > 0.5
    image_padded = np.pad(image,n,'constant')
    image_padded = dilate(image_padded,n)
    image_padded = erosion(image_padded,n)
    return image_padded[n:-n,n:-n,n:-n]>0


def labclose(image,n):
    nan_inds = np.isnan(image)
    image[nan_inds] = 0
    image = image > 0.5
    tmp = close(image,n)
    return ~lab(~tmp)


# Subsampling of CAT surfaces

def subsample_surfaces(m2m_dir, n_points):
    """Subsample the central gray matter surfaces and spherical registrations
    for each hemisphere generated by CAT. The original surfaces contain
    approximately 100,000 nodes per hemisphere.

    The subsampled surfaces are written to files with an added suffix, e.g.,
    if n_points=10000

        lh.central.gii      -> lh.central.10000.gii
        lh.sphere.reg.gii   -> lh.sphere.reg.10000.gii

    Additionally, the normals from the full resolution surfaces corresponding
    to the subsampled nodes are written to a txt file, e.g.,

        lh.central.10000.normals.txt

    PARAMETERS
    ----------
    m2m_dir : str
        Path to m2m subject folder.
    n_points : int
        Number of nodes in each subsampled hemisphere.

    RETURNS
    -------
    sub : dict
        Dictionary of the subsampled surfaces with entries (hemi, surf).
        Each entry contains rr and tris corresponding to nodes and
        triangulation, respectively. (hemi, 'central') also contains nn, the
        normals of the subsampled points (from the original surfaces).

    NOTES
    -----
    The actual number of points in the subsampled surfaces may be less than the
    requested number since sometimes multiple subsampled points may map to the
    same point on the original surface in which case only one will be kept.
    However, for consistency the filename will reflect the *requested* number
    of points rather than the *actual* number of points.
    """
    surfs = ('central', 'sphere')

    surfs_full, surfs_sub = {}, {s: {} for s in surfs}
    sf = file_finder.SubjectFiles(subpath=m2m_dir)

    for s in surfs:
        surfs_full[s] = sf.v2v2_read_surface(s)
    hemis = surfs_full[surfs[0]].keys()

    for h in hemis:
        surfs_sub["central"][h], surfs_sub["sphere"][h] = subsample_surface(surfs_full["central"][h], surfs_full["sphere"][h], n_points)

    for s in surfs:
        sf.v2v2_write_surface(surfs_sub[s], s, n_points)

    f = sf.v2v2_get_morph_data_file("normals", n_points)
    for h in hemis:
        np.savetxt(f[h].with_suffix(f"{f[h].suffix}.txt"), surfs_sub["central"][h]["normals"])
        # np.save(f[h], surfs_sub["central"][h]["normals"]) # .npy
    return surfs_sub

def subsample_surface(central_surf, sphere_surf, n_points, refine=True):
    """Subsample a hemisphere surface using its spherical registration.

    PARAMETERS
    ----------
    central_surf : dict
        Dictionary representing the surface of a hemisphere containing the keys
        "points" (points) and "tris" (triangulation).
    sphere_surf : dict
        Dictionary representing the surface of a hemisphere containing the keys
        "points" (points) and "tris" (triangulation).
    n_points : int
        Number of points (source positions) in the subsampled surface.

    RETURNS
    -------
    central_surf_sub : dict
        The subsampled surface. Also contains the key 'nn' representing the
        normals from the original surface.
    sphere_surf_sub : dict
        The subsampled surface.
    """
    assert isinstance(n_points, int)
    assert isinstance(central_surf, dict)
    assert isinstance(sphere_surf, dict)
    assert n_points < (n_full := central_surf['points'].shape[0])

    visualize = False # temporary debug flag for testing; requires pyvista

    sphere_points = sphere_surf["points"] / np.linalg.norm(sphere_surf["points"], axis=1, keepdims=True)
    tree = cKDTree(sphere_points)

    points, _ = fibonacci_sphere(n_points)
    _, idx = tree.query(points)

    # If multiple points map to the same points on the original surface then
    # keep only the unique ones.
    uniq = np.unique(idx)
    used = np.zeros(n_full, dtype=bool)
    used[uniq] = True

    n_smooth = get_n_smooth(n_points / n_full)
    basis = compute_gaussian_basis_functions(central_surf, n_smooth)
    # rescale to avoid numerical problems?
    basis.data /= basis.mean(1).mean()
    coverage = basis[:, used].sum(1) # i.e., basis @ x where x is indicator vector

    if visualize:
        surfs = {}
        add_surfs(surfs, central_surf, sphere_surf, coverage.copy(), used, "init")

    # updates `coverage` in-place
    used, unused = maximize_coverage_by_addition(used, coverage, basis, n_points - uniq.size)
    if visualize:
        add_surfs(surfs, central_surf, sphere_surf, coverage.copy(), used, "add")

    if refine:
        # updates `used`, `ununsed`, and `coverage` in-place
        indptr, indices, data = basis.indptr, basis.indices, basis.data
        equalize_coverage_by_swap(used, unused, coverage, indptr, indices, data)
        # covs, pairs, hard_swaps = equalize_coverage_by_swap(used, unused, coverage, indptr, indices, data)

    # Triangulate
    hull = ConvexHull(sphere_surf['points'][used])
    points, tris = hull.points, hull.simplices
    ensure_orientation_consistency(points, tris)

    if visualize:
        add_surfs(surfs, central_surf, sphere_surf, coverage.copy(), used, "swap")

        # visualize with pyvista
        clim = np.percentile(surfs["cent_init"]['coverage'], [5,95])
        names = ["init", "add", "swap"]

        p = pv.Plotter(shape=(2,3), window_size=(1600,1000), notebook=False)
        for i, name in enumerate(names):
            p.subplot(0,i)
            p.add_mesh(surfs[f"cent_{name}_sub"], show_edges=True)
            p.subplot(1,i)
            p.add_mesh(surfs[f"cent_{name}"], clim=clim, show_edges=True)
        # p.add_points(surfs[f'cent_{name}_sub'].points)
        p.link_views()
        p.show()

    # Use the normals from the original (high resolution) surface as this
    # should be more accurate
    normals = mesh_io.Msh(
            mesh_io.Nodes(central_surf["points"]),
            mesh_io.Elements(central_surf["tris"]+1)
            ).nodes_normals().value
    central_surf_sub = dict(
        points = central_surf["points"][used],
        tris = tris,
        normals = normals[used]
    )
    sphere_surf_sub = dict(
        points = sphere_surf["points"][used],
        tris = tris
    )

    return central_surf_sub, sphere_surf_sub


def fibonacci_sphere(n, radius=1):
    """Generate a triangulated sphere with n vertices centered on (0, 0, 0).

    PARMETERS
    ---------
    n : int
        Number of vertices of the sphere.

    RETURNS
    -------
    rr : ndarray (n, 3)
        Point coordinates.
    tris : ndarray (m, 3)
        Array describing the triangulation.

    """
    points = fibonacci_sphere_points(n) * radius
    hull = ConvexHull(points)
    points, tris = hull.points, hull.simplices
    ensure_orientation_consistency(points, tris)
    return points, tris


def fibonacci_sphere_points(n):
    """Evenly distribute n points on a unit sphere using a Fibonacci lattice.

    PARAMETERS
    ----------
    n : int
        The desired number of points.

    RETURNS
    -------
    (n, 3) array with point coordinates in rows.

    NOTES
    -----
    Based on

        http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/

    """

    # Optimize average (instead of minimum) nearest neighbor distance by
    # introducing an offset at the poles
    epsilon = 0.36

    golden_ratio = 0.5 * (1 + np.sqrt(5))
    i = np.arange(0, n, dtype=float)

    # Original fibonacci lattice
    # x2, _ = np.modf(i / golden_ratio)
    # y2 = i / n_points

    # (MODIFIED) FIBONACCI LATTICE

    x2, _ = np.modf(i / golden_ratio)
    y2 = (i + epsilon) / (n - 1 + 2*epsilon)

    # Project to fibonacci spiral via equal area projection
    # theta = 2 * np.pi * x2
    # r = np.sqrt(y2)

    # fig = plt.figure()
    # ax = fig.add_subplot()
    # ax.scatter(x2, y2)
    # ax.set_title('Fibonacci Lattice')

    # fig = plt.figure()
    # ax = fig.add_subplot(projection='polar')
    # ax.scatter(theta, r)
    # ax.set_title('Fibonacci Spiral')

    # FIBONACCI SPHERE

    # Spherical coordinates (r = 1 is implicit because it is the unit sphere)
    # theta : longitude (around sphere, 0 <= theta <= 2*pi)
    # phi   : latitude (from pole to pole, 0 <= phi <= pi)
    theta = 2*np.pi*x2
    phi = np.arccos(1 - 2*y2)

    # Cartesian coordinates
    x3 = np.cos(theta) * np.sin(phi)
    y3 = np.sin(theta) * np.sin(phi)
    z3 = np.cos(phi)

    return np.array([x3, y3, z3]).T


def ensure_orientation_consistency(points, tris):
    """Fix orientation of normals so that they all point outwards. Operates
    in-place on tris.

    PARAMETERS
    ----------
    rr : ndarray (n, 3)
        Point coordinates.
    tris : ndarray (m, 3)
        Array describing the triangulation.

    RETURNS
    -------
    None, operates in-place on tris.
    """
    # centroid_tris: vector from global centroid to centroid of each triangle
    # (in this case the global centroid is [0,0,0] and so can be ignored)
    n = mesh_io.Msh(mesh_io.Nodes(points),
                    mesh_io.Elements(tris+1)
                    ).triangle_normals().value
    centroid_tris = points[tris].mean(1)
    orientation = np.sum(centroid_tris * n, axis=1)
    swap_select_columns(tris, orientation < 0, [1, 2])

def swap_select_columns(arr, rows, cols):
    """Swap the columns (cols) of the selected rows (rows) in the array (arr).
    Operates in-place on arr.
    """
    assert not isinstance(rows, tuple)
    assert len(cols) == 2
    c0, c1 = cols
    arr[rows, c1], arr[rows, c0] = arr[rows, c0], arr[rows, c1]


def recursive_matmul(X, n):
    assert isinstance(n, int) and  n >= 1
    return X if n == 1 else recursive_matmul(X, n-1) @ X


def compute_adjacency_matrix(tris):
    """Make sparse adjacency matrix for vertices with connections `tris`.
    """
    N = tris.max() + 1

    pairs = list(itertools.combinations(np.arange(tris.shape[1]), 2))
    row_ind = np.concatenate([tris[:, i] for p in pairs for i in p])
    col_ind = np.concatenate([tris[:, i] for p in pairs for i in p[::-1]])

    data = np.ones_like(row_ind)
    return scipy.sparse.csr_array((data / 2, (row_ind, col_ind)), shape=(N, N))


def compute_gaussian_basis_functions(surf, degree):
    A = compute_adjacency_matrix(surf['tris'])

    # Estimate standard deviation for Gaussian distance computation
    Acoo = A.tocoo()
    v = surf['points']
    data = np.linalg.norm(v[Acoo.row]-v[Acoo.col], axis=1)
    # set sigma to half average distance to neighbors. This is arbitrary but
    # seems to work. Large sigmas do not seems to work well as they tend to
    # created "striped" pattern in dense regions
    sigma = 0.5*data.mean()

    # smooth to `neighborhood` degree neighbors
    A = recursive_matmul(A, degree)

    # # remove point itself
    # A = A.tolil()
    # A.setdiag(0)
    # A = A.tocsc()

    # compute gaussian distances
    Acoo = A.tocoo()
    v = surf['points']
    data = np.linalg.norm(v[Acoo.row]-v[Acoo.col], axis=1)
    A.data = np.exp(-data**2/(2*sigma**2))

    return A


def get_n_smooth(frac):
    """How much to smooth adjacency matrix.

    These numbers are pretty heuristic at the moment.
    """
    n_smooth = 3
    if frac < 0.2:
        n_smooth += 1
    if frac < 0.05:
        n_smooth += 1
    return n_smooth


@numba.jit(nopython=True, fastmath=True)
def update_vec(vec, indptr, indices, data, addi, subi):
    addi_indptr0, addi_indptr1 = indptr[addi:addi+2]
    subi_indptr0, subi_indptr1 = indptr[subi:subi+2]
    addi_indices = indices[addi_indptr0:addi_indptr1]
    subi_indices = indices[subi_indptr0:subi_indptr1]
    addi_data = data[addi_indptr0:addi_indptr1]
    subi_data = data[subi_indptr0:subi_indptr1]

    vec[addi_indices] += addi_data
    vec[subi_indices] -= subi_data


@numba.jit(nopython=True, fastmath=True)
def masked_indexed_argmin(x, index, mask, range_of_index):
    """Equivalent to (but faster than)

        ixm = index[mask]
        iargmin = x[ixm].argmin()
        ixargmin = ixm[iargmin]

    the more sparse/irregular `mask` is (since the boolean indexing operation
    becomes slow).

    """
    iargmin, ixargmin = 0, 0
    minval = np.inf
    for i in range_of_index:
        if mask[i]:
            ixi = index[i]
            val = x[ixi]
            if val < minval:
                minval = val
                iargmin = i
                ixargmin = ixi
    return iargmin, ixargmin, minval


# numba complains about the np.ones stuff. Also, A would have to be unpacked
# before. However, doesn't seem to make much difference
# @numba.jit(nopython=True, fastmath=True)
def maximize_coverage_by_addition(used, coverage, basis, n_add):
    """Add `n_add` points (move from unused to used) and update `coverage`
    accordingly (this is done in-place).

    PARAMTERS
    ---------
    used : ndarray of bool

    coverage : ndarray
        Array describing the coverage of the original points.

    n_add : int
        Number of points to add

    RETURNS
    -------
    used : ndarray
        Indices of the used points.
    unused : ndarray
        Indices of the unused points.
    """

    indptr, indices, data = basis.indptr, basis.indices, basis.data

    unused = np.where(~used)[0]
    added = -np.ones(n_add, dtype=int)
    mask = np.ones(unused.size, dtype=bool)
    range_of_index = np.arange(unused.size)

    for i in np.arange(n_add):
        # identify vertex to add
        irem, iadd, _ = masked_indexed_argmin(coverage, unused, mask, range_of_index)
        mask[irem] = False
        added[i] = iadd

        # update coverage accordingly
        start, stop = indptr[iadd:iadd+2]
        coverage[indices[start:stop]] += data[start:stop]

    used[added] = True
    unused = np.where(~used)[0]
    used = np.where(used)[0]

    return used, unused


# argpartition is not recognized by numba...
# @numba.jit(nopython=True, fastmath=True)
def equalize_coverage_by_swap(used, unused, coverage, indptr, indices, data, k=10, max_iter=5000):
    """Try to equalize coverage by swapping used and unused points. In
    particular, the objective is to minimize the variance of the coverage.

    First, it tries to swap the points with minimum and maximum coverage (add
    and remove, respectively). If this does not decrease the variance, it tries
    the next k-1 pairs. If neither decrease the variance, the process
    terminates.

    Updates `used`, `unused`, and `coverage` in-place.

    PARAMTERS
    ---------
    used : ndarray
        Indices of the used points.
    unused : ndarray
        Indices of the unused points.
    coverage : ndarray
        Array describing the coverage of the original points.

    k : partition index
        Maximum number of pairs to try.


    """
    coverage_var = coverage.var()

    # n_swaps = 0
    # covs = []
    # pairs = []
    # hard_swaps = []

    for iteration in np.arange(max_iter):

        cov_un = coverage[unused]
        cov_us = coverage[used]
        addii = cov_un.argmin()
        subii = cov_us.argmax()
        addi = unused[addii]
        subi = used[subii]

        update_vec(coverage, indptr, indices, data, addi, subi)

        # if iteration % 1000 == 0:
        #     print(iteration)

        if (coverage_swap_var := coverage.var()) < coverage_var:
            unused[addii] = subi
            used[subii] = addi
            coverage_var = coverage_swap_var

            # covs.append(coverage_swap_var)
            # n_swaps += 1
        else:
            update_vec(coverage, indptr, indices, data, subi, addi) # undo

            # hard_swaps.append(iteration)

            # this is the "slow" bit
            addiis = cov_un.argpartition(k)[:k] # [::10]
            subiis = cov_us.argpartition(-k)[-k:] # [::10]

            # order in pairs (skip the 1st as we already checked that)
            addiis = addiis[cov_un[addiis].argsort()][1:]
            subiis = subiis[cov_us[subiis].argsort()[::-1]][1:]

            found = False
            for i, (addii, subii) in enumerate(zip(addiis, subiis)):

                addi = unused[addii]
                subi = used[subii]

                update_vec(coverage, indptr, indices, data, addi, subi)

                if (coverage_swap_var := coverage.var()) < coverage_var:
                    found = True

                    unused[addii] = subi
                    used[subii] = addi

                    coverage_var = coverage_swap_var

                    # covs.append(coverage_swap_var)
                    # pairs.append(i)

                    break # next iteration
                else:
                    update_vec(coverage, indptr, indices, data, subi, addi) # undo

            if not found:
                break # if nothing to swap, terminate

    # return covs, pairs, hard_swaps

# some temporary stuff for plotting results of subsampling...
def add_surfs(surfs, central_surf, sphere_surf, coverage, used, name):
    surfs[f"cent_{name}"] = pv.make_tri_mesh(central_surf["points"], central_surf['tris'])
    surfs[f"cent_{name}"]['coverage'] = coverage
    surfs[f"sphe_{name}"] = pv.make_tri_mesh(sphere_surf["points"], sphere_surf['tris'])
    surfs[f"sphe_{name}"]['coverage'] = coverage
    hull = ConvexHull(sphere_surf["points"][used])
    rr, tris = hull.points, hull.simplices
    ensure_orientation_consistency(rr, tris)
    surfs[f"cent_{name}_sub"] = pv.make_tri_mesh(central_surf["points"][used],tris)
    surfs[f"sphe_{name}_sub"] = pv.make_tri_mesh(sphere_surf["points"][used],tris)

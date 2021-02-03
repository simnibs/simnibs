# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 18:35:47 2019

@author: axthi
"""


import gc
import logging
from math import log, sqrt
#from multiprocessing import Process
import nibabel as nib
import numpy as np
import os
#from queue import Queue, Empty
from scipy.spatial import cKDTree
from scipy.special import erf
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
             actualsurf='', ref_fs=None):
             #log_level=logging.INFO, actualsurf='', ref_fs=None): # log_level unused
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
        Smoothing of the mesh by local averaging, for each vertex which has
        been moved, the coordinates of itself and all other vertices to which
        it is connected weighting the vertex itself higher than the surrounding
        vertices (default = True).
    skip_lastsmooth : bool
        No smoothing is applied during last step to ensure that mm2move is exactly
        reached (default = True).
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

    RETURNS
    ----------
    vertices : ndarray
        Vertices describing the expanded mesh.
    """

    DEBUG = False # controls writing of additional meshes for debugging
    #  Note: ref_fs needed for debugging to have correct header information
    #        when writing FreeSurfer surfaces

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
        if DEBUG:
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

        if DEBUG:
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
            else:
                vertices = smooth_vertices(
                    vertices, faces, v2f_map=v2f, mask_move=move, taubin=True)

        logger.info(f'{actualsurf}: Moved {np.sum(move)} of {len(vertices)} vertices.')

        if DEBUG:
            tmpmsh = mesh_io.Msh(nodes=mesh_io.Nodes(vertices),
                         elements=mesh_io.Elements(faces+1))
            filename = "mesh_expand_{:d}_of_{:d}"
            filename = filename.format(i+1, nsteps)
            mesh_io.write_freesurfer_surface(tmpmsh, filename+".fsmesh", ref_fs=ref_fs)

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
        f2c, f2cok = list2numpy(f2c,dtype=np.int)
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
        m.smooth_surfaces(Niterations, nodes_mask=vert_mask)
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
        v2f, ok = list2numpy(v2f, pad_val, np.int)        
        return v2f, ok
    else:
        raise ValueError("Array output type must be list or numpy array.")    


        
def list2numpy(L, pad_val=0, dtype=np.float):
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

    tnormals = np.cross(mesh[:,1,:]-mesh[:,0,:],mesh[:,2,:]-mesh[:,0,:]).astype(np.float)
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
    #inter_z = np.around(positions[:, 2]).astype(np.int)
    inter_z = (positions[:, 2] + 1).astype(np.int)
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

def createCS(Ymf, Yleft, Ymaskhemis, Ymaskparahipp, vox2mm, actualsurf, 
             surffolder, fsavgDir, vdist=[1.0, 0.75], voxsize_pbt=[0.5, 0.25],
             voxsize_refineCS=[0.75, 0.5], th_initial=0.714, no_selfintersections=True, 
             add_parahipp=0.1, close_parahipp=False):
    """ reconstruction of cortical surfaces based on probalistic label image
    
    PARAMETERS
    ----------
        Ymf : float32
            image prepared for surface reconstruction
        Yleft : uint8
            binary mask to distinguish between left and right parts of brain
        Ymaskhemis : uint8
            volume with the following labels: lh - 1, rh - 2, lc - 3, rc - 4 
        Ymaskparahipp : uint8
            binary volume mask of hippocampi and gyri parahippocampalis, used to
            "fix" the cutting of the surfaces in these regsions
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
            
        --> these parameters will likely be skipped ...
        add_parahipp : float 
            increases the intensity values in the parahippocampal area to prevent 
            large cuts in the parahippocampal gyrus. Higher values for 
            add_parahipp will shift the initial surface in this area
            closer to GM/CSF border.
            If the parahippocampal gyrus is still cut you can try to increase 
            this value (start with 0.15).
            (default = 0.1)
        close_parahipp : bool
            optionally applies  closing inside mask for parahippocampal gyrus 
            to get rid of the holes that lead to large cuts in gyri after 
            topology correction. However, this may also lead to poorer quality 
            of topology correction for other data and should be only used if 
            large cuts in the parahippocampal areas occur.
            (AT: improved results in some parts of that area, but made it less
             accurate in others for ernie)
            (default = False)
            
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
    debug=False # keep intermediate results if set to True
    
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

    elif 'rh' == actualsurf.lower():
        Ymfs=np.multiply(Ymf,(Ymaskhemis == 2))
        Yside=np.logical_not(Yleft)

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
            
    # crop volume
    if debug:
        tmp = ndimage.uniform_filter(Ymfs, 3)
        Ymfs_filtered = nib.Nifti1Image(tmp, vox2mm)
        fname_ymfsfiltered=os.path.join(surffolder, 'Ymfs_filtered_' + actualsurf + '.nii.gz')
        nib.save(Ymfs_filtered, fname_ymfsfiltered)
    
    # Okay CAT12 uses the magical threshold of 1.5 in multiple places,
    # based on a single example (Ernie) this is not the optimal magical
    # threshold for charm. So from hereon wherever there read thres_magic,
    # it will be 1.3 for charm whereas the original one in CAT12 is 1.5
    iscat = True
    if iscat:
        thres_magic = 1.5
    else:
        thres_magic = 1.2
    thres_magic = 1.2     
    mask = ndimage.uniform_filter(Ymfs, 3) > thres_magic
    
    Ymfs, vox2mm_cropped, _ = crop_vol(Ymfs, vox2mm, mask, 4)
    
    if debug:
        Ymfs_tmp = nib.Nifti1Image(Ymfs, vox2mm_cropped)
        fname_ymfstest=os.path.join(surffolder,'Ymfs_masked_' + actualsurf + '.nii.gz')
        nib.save(Ymfs_tmp, fname_ymfstest)   
    
    Yside = crop_vol(Yside, vox2mm, mask, 4)[0]
    if not iscerebellum:
        mask_parahipp = dilate(Ymaskparahipp,6)
        mask_parahipp = crop_vol(mask_parahipp, vox2mm, mask, 4)[0]

    # upsample using linear interpolation (linear is better than cubic for small thicknesses)
    Ymfs, vox2mm_upsampled, _ =resample_vol(np.maximum(1,Ymfs), vox2mm_cropped, voxsize_pbt, order=1, mode='nearest')
    Ymfs=np.minimum(3,np.maximum(1,Ymfs))
    
    # TESTING
    if debug:
        Ymfs_upsampled = nib.Nifti1Image(Ymfs, vox2mm_upsampled)
        fname_ymfstest = os.path.join(surffolder,'Ymfs_test_upsampled_'+actualsurf+'.nii.gz')
        nib.save(Ymfs_upsampled, fname_ymfstest)
        
    Yside=resample_vol(Yside, vox2mm_cropped, voxsize_pbt, order=1, mode='nearest')[0] > 0.5
    if not iscerebellum:
        mask_parahipp = resample_vol(mask_parahipp, vox2mm_cropped, voxsize_pbt, order=1, mode='nearest')[0] > 0.5


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

    # some more tweaking of parahippocampal area in Yppi
    if not iscerebellum:
        Yppi = Yppi + add_parahipp * spm_smooth(np.float32(mask_parahipp), [8., 8., 8.])
    Yppi[Yppi < 0] = 0

    # optionally apply closing inside mask for parahippocampal gyrus to get rid of the holes that lead to large cuts in gyri
    # after topology correction
    if close_parahipp and not iscerebellum:
        tmp=labclose(Yppi, 1)
        Yppi[mask_parahipp]=tmp[mask_parahipp]
        del tmp
        gc.collect()

    if debug:
        Vppi = nib.Nifti1Image(Yppi, vox2mm_upsampled)
        fname=os.path.join(surffolder,actualsurf+'_Yppi.nii')
        nib.save(Vppi, fname)
        del Vppi
        gc.collect()
     
    # generate intial surface and save as gifti  
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

    if (np.sum(np.round(np.asanyarray(Ymf).reshape(-1, 1)) == np.asanyarray(Ymf).reshape(-1, 1)) / np.asarray(Ymf).size) > 0.9:
        binary = True
    else:
        binary = False

    minfdist = 2

    debug=int(debug)

    # CAT12's original magical threshold is 1.5
    iscat = True
    
    if iscat:
        thres_magic = 1.5
    else:
        thres_magic = 1.2
    
    thres_magic = 1.2
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
    
    ## DEBUG
    if debug:
        F_image = nib.Nifti1Image(F, vox2mm)
        fname_F=os.path.join(surffolder,'F_' + actualsurf + '.nii.gz')
        nib.save(F_image, fname_F)
    
    # In CAT12 this extends quite far into gray matter, as there are
    # values that are > 2 quite close to CSF. In our case the values drop
    # off faster, so could change the 2 into 1.9
    # YM = np.fmax(0, np.fmin(1, (Ymf - 1.9)))
    if iscat:
        tmp_param = 2
    else:
        tmp_param = 1.9
        
    YM = np.fmax(0, np.fmin(1, (Ymf - tmp_param)))
    YM[YMM] = np.nan
    
    ## DEBUG
    if debug:
        YM_image = nib.Nifti1Image(YM, vox2mm)
        fname_YM=os.path.join(surffolder,'YM_' + actualsurf + '.nii.gz')
        nib.save(YM_image, fname_YM)
    
    Ywmd = _cat_c_utils.cat_vol_eidist(
        YM, F, np.array([1, 1, 1]), 1, 1, 0, debug)[0]
    
    ## DEBUG
    if debug:
        Ywmd_image = nib.Nifti1Image(Ywmd, vox2mm)
        fname_Ywmd=os.path.join(surffolder,'Ywmd_' + actualsurf + '.nii.gz')
        nib.save(Ywmd_image, fname_Ywmd)
    
    
    F = np.fmax(1.0, np.fmin(1, Ymf / 2))
    
    YM = np.fmax(0, np.fmin(1, (Ymf - 1)))
    YM[YMM] = np.nan
    Ycsfdc = _cat_c_utils.cat_vol_eidist(
        YM, F, np.array([1, 1, 1]), 1, 1, 0, debug)[0]

    ## FOR DEBUGGING
    if debug:
        Ycsfdc_image = nib.Nifti1Image(Ycsfdc, vox2mm)
        fname_Ycsfdc=os.path.join(surffolder,'Ycsfdc_' + actualsurf + '.nii.gz')
        nib.save(Ycsfdc_image, fname_Ycsfdc)
    
    del F, YMM
    gc.collect()

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

        if iscat:
            tmp_param9 = 2.0
        else:
            tmp_param9 = 1.5

        # reducing outliers in the GM/CSF area
        notnan = ~np.isnan(Ywmd)
        YM = np.full(Ywmd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ywmd[notnan] > minfdist), (Ymf[notnan] < tmp_param9))
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
    if iscat:
        tmp_param6 = 2.5
    else:
        tmp_param6 = 2.2
        

    YMM = np.any((erosion(Ymf < thres_magic, 1), erosion(
        Ymf > tmp_param6, 1), np.isnan(Ymf)), axis=0)

    
    F = np.fmax(0.5, np.fmin(1, (4 - Ymf) / 2))

    if iscat:
        tmp_param10 = 2
    else:
        tmp_param10 = 1.5
        
    tmp_param10 = 1.6 #seems to be better than 1.5
    
    YM = np.fmax(0, np.fmin(1, (tmp_param10 - Ymf)))
    YM[YMM] = np.nan
    Ycsfd = _cat_c_utils.cat_vol_eidist(
        YM, F, np.array([1, 1, 1]), 1, 1, 0, debug)[0]
    F = np.fmax(1, np.fmin(1, (4 - Ymf) / 2))

    YM = np.fmax(0, np.fmin(1, (3 - Ymf)))
    YM[YMM] = np.nan
    Ywmdc = _cat_c_utils.cat_vol_eidist(
        YM, F, np.array([1, 1, 1]), 1, 1, 0, debug)[0]
    
    # The original CAT12 command was: 
    # YM = np.fmax(0, np.fmin(1, (2.7 - Ymf)))
    # This looks slightly fishy as it's not really a half
    # Let's try to set this to 3 and see what happens
    if iscat:
        tmp_param2 = 2.7
    else:
        tmp_param2 = 2.2
        
    YM = np.fmax(0, np.fmin(1, (tmp_param2 - Ymf)))
    
    YM[YMM] = np.nan
    
    # I presume the 0.3 links to the 2.7 above so let's remove it
    # Here the original CAT12 call was:
    if iscat:
        Ywmdx = _cat_c_utils.cat_vol_eidist(
            YM, F, np.array([1, 1, 1]), 1, 1, 0, debug)[0] + 0.3
    else:
        Ywmdx = _cat_c_utils.cat_vol_eidist(
            YM, F, np.array([1, 1, 1]), 1, 1, 0, debug)[0] + 0.8
    
    del F, YMM
    gc.collect()
    
    ## DEBUG
    if debug:
        Ywmdx_image = nib.Nifti1Image(Ywmdx, vox2mm)
        fname_Ywmdx=os.path.join(surffolder,'Ywmdx_'+actualsurf+'.nii.gz')
        nib.save(Ywmdx_image, fname_Ywmdx)

    Ywmdc = np.fmin(Ywmdc, Ywmdx)
    
    ## DEBUG
    if debug:
        Ywmdc_image = nib.Nifti1Image(Ywmdc, vox2mm)
        fname_Ywmdc=os.path.join(surffolder,'Ywmdc_'+ actualsurf+ '.nii.gz')
        nib.save(Ywmdc_image, fname_Ywmdc)
    
    ## DEBUG
    if debug:
        Ycsfdc_image = nib.Nifti1Image(Ycsfd, vox2mm)
        fname_Ycsfdc=os.path.join(surffolder,'Ycsfdc_before_a_million_steps_' + actualsurf + '.nii.gz')
        nib.save(Ycsfdc_image, fname_Ycsfdc)
    
    if not binary:
        notnan = ~np.isnan(Ycsfd)
        YM = np.full(Ycsfd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ycsfd[notnan] > minfdist), (Ymf[notnan] >= tmp_param6))
        Ycsfd[YM] = Ycsfd[YM] - Ywmdc[YM]
        Ycsfd[np.isinf(- Ycsfd)] = 0
        del Ywmdc
        gc.collect()
        notnan = ~np.isnan(Ycsfd)
        YM = np.full(Ycsfd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ycsfd[notnan] > minfdist), (Ymf[notnan] < tmp_param6))
        YcsfdM = np.array(Ycsfd, copy=True)
        YcsfdM = _cat_c_utils.cat_vol_localstat(YcsfdM, YM, 1, 1)[0]
        Ycsfd[YM] = YcsfdM[YM]
        notnan = ~np.isnan(Ycsfd)
        YM = np.full(Ycsfd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ycsfd[notnan] > minfdist), (Ymf[notnan] >= tmp_param6))
        YcsfdM = np.array(Ycsfd, copy=True)
        for i in np.arange(1, 3):
            YcsfdM = _cat_c_utils.cat_vol_localstat(YcsfdM, YM, 1, 1)[0]
        Ycsfd[YM] = YcsfdM[YM]
        notnan = ~np.isnan(Ycsfd)
        YM = np.full(Ycsfd.shape, False, dtype=bool)
        YM[notnan] = np.logical_and(
            (Ycsfd[notnan] > minfdist), (Ymf[notnan] > tmp_param9))
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
    # 2.2 was 2.5 originally
    if iscat:
        tmpparam3 = 2.5
    else:
        tmpparam3 = 2.2
        
    Ygmt = np.fmin(Ygmt1, Ygmt2)
    YM = np.logical_and((Ymf >= thres_magic), (Ymf < tmpparam3))
    Ypp = np.zeros(Ymf.shape, dtype=np.float32)
    Ypp[Ymf >= tmpparam3] = 1
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

    # filter result 1.3 was 1.8 in original CAT12 code
    if iscat:
        tmpparam4 = 1.8
    else:
        tmpparam4 = 1.2
        
    Ygmts = np.array(Ygmt1, copy=True)
    for i in np.arange(1, iterator):
        Ygmts = _cat_c_utils.cat_vol_localstat(Ygmts, (((Ygmt > 1) | (Ypp > 0.1)) & (
            Ygmt > 0) & ((Ygmt > 1) | (Ymf > tmpparam4))), 1, 1)[0]

    Ygmt1[Ygmts > 0] = Ygmts[Ygmts > 0]
    Ygmts = np.array(Ygmt2, copy=True)
    for i in np.arange(1, iterator):
        Ygmts = _cat_c_utils.cat_vol_localstat(Ygmts, (((Ygmt > 1) | (Ypp > 0.1)) & (
            Ygmt > 0) & ((Ygmt > 1) | (Ymf > tmpparam4))), 1, 1)[0]

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
    
    # The original CAT12 version of this loop was:
    # for i in np.arange(1, iterator):
    #     Ygmts = _cat_c_utils.cat_vol_localstat(Ygmts, (((Ygmt > 1) | (Ypp > 0.1)) & (
    #         Ygmts > 0) & ((Ygmt > 1) | (Ymf > 1.8))), 1, 1)[0]
    # Here the 1.8 looks suspicious, so changing that to 2
    
    for i in np.arange(1, iterator):
        Ygmts = _cat_c_utils.cat_vol_localstat(Ygmts, (((Ygmt > 1) | (Ypp > 0.1)) & (
            Ygmts > 0) & ((Ygmt > 1) | (Ymf > tmpparam4))), 1, 1)[0]

    Ygmt[Ygmts > 0] = Ygmts[Ygmts > 0]

    # Estimation of a mixed percentual possion map Ypp. 2.3 should be 2.5
    if iscat:
        tmpparam5 = 2.5
    else:
        tmpparam5 = 2.2
        
    YM = ((Ymf >= thres_magic) & (Ymf < tmpparam5) & (Ygmt > eps))
    Ycsfdc = np.array(Ycsfd, copy=True)
    Ycsfdc[YM] = np.fmin(Ycsfd[YM], Ygmt[YM] - Ywmd[YM])
    Ypp = np.zeros(Ymf.shape, dtype=np.float32)
    Ypp[Ymf >= tmpparam5] = 1
    Ypp[YM] = Ycsfdc[YM] / (Ygmt[YM] + eps)
    Ypp[Ypp > 2] = 0
    notnan = ~np.logical_or(np.isnan(Ywmd), np.isnan(Ygmt))
    YM = np.full(Ywmd.shape, False, dtype=bool)
    YM[notnan] = np.squeeze((Ygmt[notnan] <= resV) & (
        Ywmd[notnan] <= resV) & (Ygmt[notnan] > 0))
    
    # Here the original CAT12 call was 
    # Ypp[YM] = (Ymf[YM] - 1) / 2 - 0.2
    # Maybe the 0.2 is linked to the 1.8 above, 
    # so getting rid of that
    if iscat:
        Ypp[YM] = (Ymf[YM] - 1) / 2 - 0.2
    else:
        Ypp[YM] = (Ymf[YM] - 1) / 2

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


def refineCS(Praw, fname_thkimg, fname_ppimg, fsavgDir, vdist=1.0, no_selfintersections=True, debug=False):
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
        # region 4: pre-final surface with thickness and perc. positions (as node data)
        # region 5: final surface
        

    # ------- mark topological defects -------- 
    stimet = time.time()
    
    # spherical surface mapping 1 of the uncorrected surface for topology correction
    cmd = f'\"{file_finder.path2bin("CAT_Surf2Sphere")}\" \"{Praw}\" \"{Psphere0}\" 5' 
    spawn_process(cmd)

    # estimate size of topology defects (in relation to number of vertices and mean brain with 100000 vertices)
    cmd = f'\"{file_finder.path2bin("CAT_MarkDefects")}\" \"{Praw}\" \"{Psphere0}\" \"{Pdefects0}\"'     
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
    
    cmd=f'\"{file_finder.path2bin("CAT_FixTopology")}\" -lim 128 -bw 512 -n 81920 -refine_length {2 * vdist} \"{Praw}\" \"{Psphere0}\" \"{Pcentral}\"' 
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
        force_no_selfintersections = 1
    else:
        force_no_selfintersections = 0
        
    cmd=f'\"{file_finder.path2bin("CAT_DeformSurf")}\" \"{fname_ppimg}\" none 0 0 0 \"{Pcentral}\" \"{Pcentral}\" none 0 1 -1 .1 avg -0.1 0.1 .2 .1 5 0 0.5 0.5 n 0 0 0 150 0.01 0.0 {force_no_selfintersections}' 
    spawn_process(cmd)
    
    if no_selfintersections:
        # remove self-intersections using meshfix
        CS = mesh_io.read_gifti_surface(Pcentral)
        mesh_io.write_off(CS, Pcentral+'.off')
        
        cmd=f'\"{file_finder.path2bin("meshfix")}\"  \"{Pcentral+".off"}\" -o \"{Pcentral+".off"}\"'
        spawn_process(cmd)
        
        CS = mesh_io.read_off(Pcentral+'.off')
        mesh_io.write_gifti_surface(CS,Pcentral)
        if os.path.isfile(Pcentral+'.off'):
            os.remove(Pcentral+'.off')
        del CS
        gc.collect()
        
    # need some more refinement because some vertices are distorted after CAT_DeformSurf
    cmd=f'\"{file_finder.path2bin("CAT_RefineMesh")}\" \"{Pcentral}\" \"{Pcentral}\" {"%.2f" % (1.5 * vdist )} 0' 
    spawn_process(cmd)

    cmd=f'\"{file_finder.path2bin("CAT_DeformSurf")}\" \"{fname_ppimg}\" none 0 0 0 \"{Pcentral}\" \"{Pcentral}\" none 0 1 -1 .2 avg -0.05 0.05 .1 .1 5 0 0.5 0.5 n 0 0 0 50 0.01 0.0 {force_no_selfintersections}' 
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


    # ---- final correction of central surface in highly folded areas -------- 
    # ----------------  with high mean curvature -------------- 
    stimet = time.time()

    cmd = f'\"{file_finder.path2bin("CAT_Central2Pial")}\" -equivolume -weight 0.3 \"{Pcentral}\" \"{Pthick}\" \"{Pcentral}\" 0' 
    spawn_process(cmd)
    
    if debug:
        # add final central surface to .msh and save .msh
        CS = mesh_io.read_gifti_surface(Pcentral)
        CS_dbg.elm.add_triangles(CS.elm.node_number_list[:,0:3]+CS_dbg.nodes.nr,5)
        CS_dbg.nodes.node_coord = np.concatenate((CS_dbg.nodes.node_coord, CS.nodes.node_coord))
        mesh_io.write_msh(CS_dbg,Pdebug)
        del CS
        gc.collect()

    logger.info(f'Correction of central surface in highly folded areas 2: '+time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))


    # -------- registration to FSAVERAGE template --------
    stimet = time.time()
    
    # spherical surface mapping 2 of corrected surface    
    cmd = f'\"{file_finder.path2bin("CAT_Surf2Sphere")}\" \"{Pcentral}\" \"{Psphere}\" 10' 
    spawn_process(cmd)
    
    # spherical registration to fsaverage template
    cmd = f'\"{file_finder.path2bin("CAT_WarpSurf")}\" -steps 2 -avg -i \"{Pcentral}\" -is \"{Psphere}\" -t \"{Pfsavg}\" -ts \"{Pfsavgsph}\" -ws \"{Pspherereg}\"' 
    spawn_process(cmd)

    logger.info(f'Registration to FSAVERAGE template: '+time.strftime('%H:%M:%S', time.gmtime(time.time() - stimet)))

    
    # -------- remove unnecessary files -------- 
    if not debug:
        if os.path.isfile(Psphere0):
            os.remove(Psphere0)    
            
        if os.path.isfile(Pdefects0):
            os.remove(Pdefects0)
            
        if os.path.isfile(Psphere):
            os.remove(Psphere)
            
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


def spm_smoothkern(fwhm, t=1):
    """  Generate a Gaussian smoothing kernel

    PARAMETERS
    ----------
        fwhm : full width at half maximum

        --> optional parameters:
        t    : either 0 (nearest neighbour) or 1 (linear).
            [Default: 1]

    RETURNS
    ----------
        krn  : value of kernel at position x

    NOTES
    ----------      
        This function returns a Gaussian convolved with a triangular (1st 
        degree B-spline) basis function.

        It is adapted from spm_smoothkern.m of SPM12
        (version 2019-03-22, https://www.fil.ion.ucl.ac.uk/spm/software/spm12/).

     Reference
     ----------
        John Ashburner
        $Id: spm_smoothkern.m 7460 2018-10-29 15:55:12Z john $

        Martin TB, Prunet S, Drissen L. Optimal fitting of Gaussian-apodized
        or under-resolved emission lines in Fourier transform spectra
        providing new insights on the velocity structure of NGC 6720. Monthly
        Notices of the Royal Astronomical Society. 2016 Sep 14;463(4):4223-38.

    """

    length = np.rint(3 * fwhm / sqrt(2*log(2)))
    x = np.arange(-length, length+1, 1)

    # Variance from FWHM
    s = fwhm ** 2 / (8*log(2)) + np.finfo(float).eps

    if t == 0:
        # Gaussian convolved with 0th degree B-spline
        w1 = 1.0 / sqrt(2*s)
        krn = 0.5*(erf(w1 * (x + 0.5)) - erf(w1 * (x - 0.5)))

    elif t == 1:
        # Gaussian convolved with 1st degree B-spline
        w1 = 0.5 * sqrt(2/s)
        w2 = - 0.5 / s
        w3 = sqrt(s * 0.5 / np.pi)
        krn = 0.5*(erf(w1 * (x + 1)) * (x + 1) + erf(w1 * (x - 1)) * (x - 1) - 2 * erf(w1 * x) * x) + \
            w3 * (np.exp(w2 * ((x + 1) ** 2)) + np.exp(w2 *
                                                       ((x - 1) ** 2)) - 2*np.exp(w2 * (x ** 2)))
    else:
        logger.error(
            f'Only defined for nearest neighbour and linear interpolation.')
        raise ValueError('spm_smoothkern only supports zero and first order')

    krn[krn < 0] = 0
    krn = krn / np.sum(krn)
    return krn
        

def spm_smooth(P, s):
    """ Convolve a three dimensional image

    PARAMETERS
    ----------    
        P     : 3D array to be smoothed
        s     : [sx sy sz] Gaussian filter width {FWHM} in edges

    RETURNS
    ----------
        Q     : 3D array of the smoothed image (or 3D array)    

    NOTES
    ----------   
        spm_smooth is used to smooth or convolve image. This function is
        adapted from spm_smooth.m and smooth1.m of SPM12
        (version 2019-03-22, https://www.fil.ion.ucl.ac.uk/spm/software/spm12/).

        The sum of kernel coeficients are set to unity.  Boundary
        conditions assume data does not exist outside the image in z (i.e.
        the kernel is truncated in z at the boundaries of the image space). S
        can be a vector of 3 FWHM values that specifiy an anisotropic
        smoothing.

        The inconsistencies in dealing with the convolution in the x/y 
        direction and z direction in the C function spm_conv_vol.c is removed 
        from the python version spm_smooth(P, s).

     Reference
     ----------        
        John Ashburner & Tom Nichols
        $Id: spm_smooth.m 4419 2011-08-03 18:42:35Z guillaume $

    """

    x = spm_smoothkern(s[0], 1)
    y = spm_smoothkern(s[1], 1)
    z = spm_smoothkern(s[2], 1)

    # Convolve the image in dimension x, y and z iteratively
    Q = P.copy()
    Q[np.isinf(Q)] = 0
    for i, k in enumerate((x.flatten(), y.flatten(), z.flatten())):
        Q = ndimage.convolve1d(
            Q, k, axis=i, mode='constant', cval=0.0, origin=0)

    return Q


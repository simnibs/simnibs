# -*- coding: utf-8 -*-
"""
Created on Sat Nov 23 18:35:47 2019

@author: axthi
"""

from threading import Thread
from multiprocessing import Process
from queue import Queue, Empty

import logging
import numpy as np
from scipy.spatial import cKDTree
import sys
from subprocess import Popen, PIPE
import scipy.ndimage.morphology as mrph
from scipy.ndimage.measurements import label

from ..mesh_tools import mesh_io
from ..utils.simnibs_logger import logger


def expandCS(vertices_org, faces, mm2move_total, ensure_distance=0.2, nsteps=5,
             deform="expand", smooth_mesh=True, skip_lastsmooth=True,
             smooth_mm2move=True, despike_nonmove=True, fix_faceflips=True,
             log_level=logging.INFO, actualsurf='', ref_fs=None):
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
    log_level : integer
        logging level (default = logging.INFO).
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
            *tuple(map(np.arange, out_shape[:2])), indexing="ij")
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
    lines_intersecting, uq_indices, uq_inverse, counts = np.unique(
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

def dilate(image, n):
    se = np.ones((2*n+1, 2*n+1, 2*n+1), dtype=bool)
    return mrph.binary_dilation(image, se) > 0

def erosion(image, n):
    return ~dilate(~image, n)

def lab(image):
    labels, num_features = label(image)
    return (labels == np.argmax(np.bincount(labels.flat)[1:])+1)

def close(image, n):
    image_padded = np.pad(image, n, 'constant')
    image_padded = dilate(image_padded, n)
    image_padded = erosion(image_padded, n)
    return image_padded[n:-n, n:-n, n:-n] > 0


def labclose(image, n):
    tmp = close(image, n)
    return ~lab(~tmp)


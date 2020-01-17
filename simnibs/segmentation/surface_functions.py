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
import os
from scipy.spatial import cKDTree
import struct
import sys
from subprocess import Popen, PIPE
from itertools import islice

from simnibs import mesh_io


logger = logging.getLogger("py.warnings")



# def deform_mesh(vertices, faces, mm2move, ensure_distance, nsteps,
#                 deform="expand", plane_type="one-sided", smooth_mesh=False,
#                 save=None, file_format="off", eps=1e-6, verbose=False):
#     """Deform a mesh by either expanding it in the direction of node normals
#     or shinking it in the opposite direction of the normals.
    
#     PARAMETERS
#     ----------
#     vertices : ndarray
#         Vertices describing the mesh.
#     faces : ndarray
#         Faces describing the mesh.    
#     mm2move : float
#         The amount by which to move a vertice each iteration.
#     ensure_distance : float
#         Minimum distance to ensure in the mesh.
#     nsteps : int
#         The number of steps used to reach the desired distance to move.
#     deform : {"expand", "shrink"}
#         Whether to expand or shink the mesh. Expand corresponds to pushing
#         the nodes outwards (in the direction of the normals) whereas shrink
#         will pull the nodes inwards (in the direction of minus the normals)
#         (default = "expand").
#     plane_type : {"one-sided", "two-sided"}, optional
#         Which type of intersections to consider. If "one-sided", only
#         intersections with the face of a triangle are considered (i.e. the
#         direction of movement and triangle normal are opposite). If
#         "two-sided", all intersections are considered (default = "one-sided").
#     smooth_mesh : bool, optional
#         Local smoothing of the mesh by averaging, for each vertex which has
#         been moved, the coordinates of itself and all other vertices to which
#         it is connected weighting the vertex itself higher than the surrounding
#         vertices (default = True).
#     save : {"all", "final", None}, optional
#         Which steps of the process to save (default = None).
#     file_format : {"off", "stl"}, optional
#         Which file format to save the mesh to (default = "off").    
#     eps : float, optional
#         Epsilon. Error tolerance margin. For numerical stability of results
#         (default = 1e-6).
#     verbose : bool
#         Whether or not to print to console (default = False).
    
#     RETURNS
#     ----------
#     vertices : ndarray
#         Vertices describing the expanded mesh.
#     """
    
#     verbosity = logging.DEBUG
    
#     # check inputs
#     assert deform in ["expand", "shrink"]
#     assert isinstance(nsteps,int)
#     assert plane_type in ["one-sided", "two-sided"]
#     assert save in ["all", "final", None]
#     assert file_format in ["off", "stl"]
    
#     mm2move = mm2move/float(nsteps)
#     n_nodes = len(vertices)
#     if isinstance(mm2move, (list, tuple, np.ndarray)):
#         assert len(mm2move) == len(vertices), "The length of mm2move must match that of vertices"
#         ballsize = max(mm2move)
#         isarray = True
#         verts2consider = np.where(mm2move!=0)[0] # no need to consider others
#     else:
#         assert isinstance(mm2move,(int,float))
#         ballsize = mm2move
#         isarray = False
#         verts2consider = np.arange(n_nodes)
        
#     vertices = vertices.copy() # prevent modification of input "vertices"!
#     v2f = verts2faces(vertices,faces)
      
#     for i in range(nsteps):
#         #if verbose:
#         #    log("Step {0}",i+1)
#         sys.stdout.flush()
        
#         mesh = vertices[faces]
        
#         # get node normals by averaging normals of the triangles to which the
#         # node belongs
#         node_normals = get_node_normals(vertices, faces, verts2consider)
#         if deform == "shrink":
#             node_normals *= -1
        
#         # testing all nodes against all triangles is slow if the mesh is large,
#         # thus, find triangles within some distance of each node and test only
#         # these for intersections. This reduces # of computations dramatically.    
#         barycenters = np.average(mesh, axis=1)
#         bar_tree = cKDTree(barycenters)
#         ver_tree = cKDTree(vertices[verts2consider])
#         res = ver_tree.query_ball_tree(bar_tree, r=(2+eps)*ballsize+2*ensure_distance,
#                                         p=2)
#         # possible intersections
#         pi, piok = list2numpy(res, dtype=np.int)

#         intersect = ray_triangle_intersect(mesh, vertices[verts2consider],
#                         node_normals, pi, piok, plane_type, mindist2int=0,
#                         maxdist2int=ballsize+ensure_distance)
#         move = ~intersect.any(1)
#         verts2consider = verts2consider[move]
        
#         if verbose:
#             pad = str(len(str(n_nodes)))
#             msg = "Moved {:"+pad+"d} of {:d} vertices"
#             log(msg, (np.sum(move), n_nodes), level=verbosity)
#             sys.stdout.flush()

#         # update the vertices
#         if isarray:
#             vertices[verts2consider] += node_normals[move]*mm2move[verts2consider,None]
#         else:
#             vertices[verts2consider] += node_normals[move]*mm2move
        
#         if smooth_mesh:
#             # smooth verts2consider and neighbors
#             if len(verts2consider) < len(vertices):
#                 f2c = [v2f[n] for n in verts2consider]
#                 f2c,f2cok = list2numpy(f2c,dtype=np.int)
#                 f2c = f2c[f2cok] # faces of verts2consider
#                 nn, nnok = get_element_neighbors(vertices[faces]) # get neighboring elements
#                 neighbors = nn[f2c][nnok[f2c]]
#                 # get neighboring vertices and verts2consider themselves
#                 v2smooth = np.unique(np.append(faces[neighbors],verts2consider)).astype(np.int)
#             else:
#                 v2smooth = verts2consider
                
#             smoo = np.zeros_like(vertices[v2smooth])
#             for i,n in zip(range(len(v2smooth)),v2smooth):
#                 smoo[i] = np.average(vertices[faces[v2f[n]]], axis=(0,1))
#             vertices[v2smooth] = smoo

#         if save == "all" or (save == "final" and (i+1)==nsteps):
#             pad = str(len(str(nsteps))-1)
#             filename = "mesh_expand_{:0"+pad+"d}_of_{:d}."+file_format
#             filename = filename.format(i+1, nsteps)
#             mesh_save(vertices, faces, filename)
    
#     return vertices




 
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
    
    DEBUG=False # controls writing of additional meshes for debugging
    #  Note: ref_fs needed for debugging to have correct header information 
    #        when writing FreeSurfer surfaces

    # check inputs
    assert deform in ["expand", "shrink"]
    assert isinstance(nsteps,int)
    assert len(mm2move_total) == len(vertices_org), "The length of mm2move must match that of vertices"

    vertices = vertices_org.copy() # prevent modification of input "vertices"
    move=np.ones(len(vertices),dtype='bool')
    v2f = verts2faces(vertices,faces)
    for i in range(nsteps):
        sys.stdout.flush()
        log(actualsurf+': Iteration '+str(i+1)+' of '+str(nsteps), level=logging.INFO)
        
        # ---------------------------------
        # update mm2move and vertex normals
        # ---------------------------------
        if i==0:
            mm2move=mm2move_total/float(nsteps)
        else:
            # distance adjustment is needed to account for small vertex shifts caused by smoothing
            dist=np.sqrt(np.sum((vertices-vertices_org)**2, axis=1)) 
            mm2move=(mm2move_total-dist)/float(nsteps-i)
        mm2move[~move] = 0
        
        node_normals = get_node_normals(vertices, faces)
        if deform == "shrink":
            node_normals *= -1
        
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
        # * ray_triangle_intersect works only properly for smoothed vc_tst below, 
        #   unclear why!? --> enforcing smoothing here, independent of option settings
        # * mm2move is not smoothed here. A temporary copy of mm2move could be 
        #   smoothed to make testing more similar to final movement
        # * for the temporarily shifted mesh, one intersection is expected,
        #   as each non-shifted vertex will intersect with its shifted version
            
        mesh=vertices[faces]
        barycenters = np.average(mesh, axis=1)
        bar_tree = cKDTree(barycenters)
        ver_tree = cKDTree(vertices[move])
        res = ver_tree.query_ball_tree(bar_tree, r=1.001*max(mm2move)+ensure_distance, p=2)
        pi, piok = list2numpy(res, dtype=np.int) # possible intersections
        
        facenormals_pre = get_triangle_normals(mesh)
        
        intersect = ray_triangle_intersect(mesh, vertices[move],
                        node_normals[move], pi, piok, plane_type="two-sided", mindist2int=0,
                        maxdist2int=mm2move[move,None]+ensure_distance)
        
        # create temporary shifted mesh and test again for intersections
        vc_tst = vertices.copy()
        vc_tst[move] += node_normals[move]*mm2move[move,None]
        vc_tst=smooth_vertices(vc_tst,faces,v2f_map=v2f,mask_move=move)
        
        mesh=vc_tst[faces]
        barycenters = np.average(mesh, axis=1)
        bar_tree = cKDTree(barycenters)
        ver_tree = cKDTree(vertices[move])
        res = ver_tree.query_ball_tree(bar_tree, r=1.001*max(mm2move)+ensure_distance, p=2)
        pi, piok = list2numpy(res, dtype=np.int) # possible intersections
        
        intersect2 = ray_triangle_intersect(mesh, vertices[move],
                        node_normals[move], pi, piok, plane_type="two-sided", mindist2int=0,
                        maxdist2int=mm2move[move,None]+ensure_distance)
                
        if DEBUG:
            move_backup=move.copy()        
        move[move]=(intersect.sum(1)==0)&(intersect2.sum(1)<2)

        # -------------------------------------
        # remove isolated "non-move" vertices
        # this is needed as the ray-triangle intersection testing
        # returns a few spurious false positives
        # --------------------------------------
        if despike_nonmove:        
            Nnomove=np.zeros(len(move),dtype='uint16')
            Nfaces=np.zeros(len(move),dtype='uint16')
            for j in range(len(move)):
                Nnomove[j]=np.sum(~move[faces[v2f[j]]])
                Nfaces[j]=len(v2f[j])
            # a single vertex reoccurs #faces --> Nnomove>Nfaces will be true 
            # when more than one vertex is marked "non-move"
            move=~(~move&(Nnomove>Nfaces))
    
        # ----------------------------
        # update the vertex positions, fix flipped faces by local smoothing
        # ----------------------------
        mm2move[~move]=0
        if smooth_mm2move:
            mm2move=smooth_vertices(mm2move,faces,Niterations=1)
        
        if DEBUG:
            vertices_beforemove=vertices.copy()
            
        vertices += node_normals*mm2move[:,None]
        
        # test for flipped surfaces
        mesh=vertices[faces]
        facenormals_post = get_triangle_normals(mesh)
        flipped_faces=np.sum(facenormals_post*facenormals_pre,axis=1)<0
        if fix_faceflips&np.any(flipped_faces):
            log(actualsurf+': Fixing '+str(np.sum(flipped_faces))+' flipped faces', level=logging.DEBUG)
            vertices=smooth_vertices(vertices,faces,verts2consider=np.unique(faces[flipped_faces]),
                                     v2f_map=v2f,Niterations=5,Ndilate=2)
            mesh=vertices[faces]
            facenormals_post = get_triangle_normals(mesh)
            flipped_faces=np.sum(facenormals_post*facenormals_pre,axis=1)<0
            log(actualsurf+': '+str(np.sum(flipped_faces))+' flipped faces remaining', level=logging.DEBUG)
            
        if smooth_mesh:
            if skip_lastsmooth&(i==nsteps-1):
                log(actualsurf+': Last iteration: skipping vertex smoothing', level=logging.DEBUG)
            else:
                vertices=smooth_vertices(vertices,faces,v2f_map=v2f,mask_move=move)
                
        log(actualsurf+': Moved '+str(np.sum(move))+' of '+str(len(vertices))+' vertices.', level=logging.INFO)
        
        if DEBUG:
            tmpmsh = mesh_io.Msh(nodes=mesh_io.Nodes(vertices),
                             elements=mesh_io.Elements(faces+1))
            filename = "mesh_expand_{:d}_of_{:d}"
            filename = filename.format(i+1, nsteps)
            mesh_io.write_freesurfer_surface(tmpmsh,filename+".fsmesh", ref_fs=ref_fs)
            
            tmpmsh.add_node_field(move,'move')
            
            hlpvar=np.zeros(move.shape)
            hlpvar[move_backup]=intersect.sum(1)
            tmpmsh.add_node_field(hlpvar,'intersect.sum(1)')
            hlpvar[move_backup]=intersect2.sum(1)
            tmpmsh.add_node_field(hlpvar,'intersect2.sum(1)')
            tmpmsh.add_node_field(mm2move_total,'mm2move_total')
                        
            tmpmsh.elm.add_triangles(faces+tmpmsh.nodes.nr+1,3)
            tmpmsh.nodes.node_coord=np.concatenate((tmpmsh.nodes.node_coord,vertices_beforemove))
            tmpmsh.add_element_field(np.concatenate((flipped_faces,flipped_faces)),'flipped_faces')
            
            tmpmsh.elm.tag2=tmpmsh.elm.tag1
            tmpmsh.write(filename+".msh")
                    
    return vertices
    


def smooth_vertices(vertices,faces,verts2consider=None,v2f_map=None,Niterations=1,Ndilate=0,mask_move=None):
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
        f2c,f2cok = list2numpy(f2c,dtype=np.int)
        f2c = f2c[f2cok] # faces of verts2consider
        verts2consider=np.unique(faces[f2c])

    if not mask_move is None:
         verts2consider = verts2consider[mask_move[verts2consider]]
        
    smoo = vertices.copy()
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



def get_node_normals(vertices, faces, mask=None, smooth=False,
                     weigh_by_area=False):
    """Compute the normal of vertices in a mesh by averaging normals of the 
    triangles to which the node belongs
    
    PARAMETERS
    ----------
    vertices : ndarray
        Vertices describing the mesh.
    faces : ndarray
        Faces describing the mesh.
    mask : array_like
        Array of indices or boolean array describing which vertices to return
        the normals of. Returns normals for all vertices by default.
    smooth : bool, optional
        Whether or not to smooth vertice normals. A vertice normal is smoothed
        by taking into consideration also the normals of the neighboring
        vertices when computing the normal (default = False).
    weigh_by_area : bool, optional
        When computing a vertice normal, weight each triangle normal by the
        area of the triangle so that large triangles contribute more and vice
        versa (default = False).
    
    RETURNS
    ----------
    node_normals : ndarray
        The vertex normals.
    """
    verts2consider = np.arange(len(vertices))
    if not mask is None:
        verts2consider = verts2consider[mask]

    if smooth:
        assert len(verts2consider) == len(vertices),"Cannot smooth normals when using mask"
    v2f = verts2faces(vertices, faces)    
    triangle_normals = get_triangle_normals(vertices[faces])
    node_normals = []
    
    # Get the node normals
    if weigh_by_area:
        mesh = vertices[faces]
        v0 = mesh[:,0,:]
        e1 = mesh[:,1,:]-v0
        e2 = mesh[:,2,:]-v0
            
        area = np.linalg.norm(np.cross(e1,e2),axis=1)/2.
        for n in verts2consider:
            w = area[v2f[n]]
            node_normals.append(np.average(triangle_normals[v2f[n]], axis=0,
                                           weights=w/np.sum(w)))
    else:
        for n in verts2consider:
            node_normals.append(np.average(triangle_normals[v2f[n]], axis=0))
    
    # Normalize
    node_normals  = np.array(node_normals)
    node_normals /= np.sqrt(np.sum(node_normals**2,1))[:,np.newaxis]
    
    # Smooth node normals by averaging neighboring node normals
    if smooth:
        s = np.zeros_like(node_normals)
        for n in verts2consider:
            s[n] = np.average(node_normals[faces[v2f[n]]], axis=(0,1))
                
        node_normals  = s
        node_normals /= np.sqrt(np.sum(node_normals**2,1))[:,np.newaxis]
    
    return node_normals



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


    
def ray_triangle_intersect(triangles, ray_origin, ray_direction, posint=None,
                           posint_ok=None, plane_type="two-sided",
                           mindist2int=-np.inf, maxdist2int=np.inf,
                           return_int_points=False, return_int_dists=False, eps=1e-6):
    """Test intersection between rays and triangles.
    
    PARAMETERS
    ----------
    triangles : ndarray
        Nx3x3 where N is # of triangles each of which are described by a 3-by-3
        (vertices-by-coordinates) array.
    ray_origin : ndarray
        Points of ray origin.
    ray_direction : ndarray
        Vector describing the direction of the ray originating from ray_origin.
        This may be specified (1) one vecetor per point in ray_origin, or (2)
        a single vector in which case the ray direction is assumed to be the
        same for all origins.
    posint : ndarray, optional
        An array describing, for each point in ray_origin, which triangles the
        ray could possibly intersect. If testing a large number of points
        and/or triangles it is highly recommended to do an initial (coarse)
        check for possible point-triangle intersection pairs. Otherwise, all
        points will need to be checked against all triangles, which (1) is
        likely to be very redundant, and (2) may even result in a MemoryError
        being raised.
        Specified either as a numpy.ndarray or list of lists (in the latter
        case, posint_ok need not be specified).
    posint_ok ndarray, optional
        An array of shape posint describing which entries in posint are real
        (actual) possible intersection pairs. This array is needed since it is
        highly unlikely that all ray origin points will have the exact same
        number of possible intersections, hence the need to somehow pad the
        array. posint_ok specified which values are 'real' and which are the
        result of padding. Note, however, that this is not the case if posint
        is specified as a list of lists containing, for each ray origin, all
        possible triangles it might intersect. Lists, however, does not support
        vectorized operations, and so will need to be converted to a
        numpy.ndarray in any case.
    plane_type : {"one-sided", "two-sided"}, optional
        The point of intersection between a ray and the plane spanned by the
        triangle sides may be computed as some multiple of the ray_direction
        (from ray_origin). If one-sided, accept only intersections in the
        actual direction of the ray (positive multiple). If two-sided, accept
        all intersections no matter of there are in the actual direction of the
        ray or in the opposite direction (default = "two-sided").
    mindist2int : float
        Minimum distance for which to consider intersections
        (default = -np.inf).
    maxdist2int : float
        Maximum distance for which to consider intersections
        (default = np.inf).
    return_int_points : bool, optional
        Return the exact points where each ray crosses the plane of each
        triangle (default = False).
    return_int_dists : bool, optional
        Return the distance from ray_origin to the point of intersection with
        the plane spanned by each triangle (default = False).
    eps : float, optional
        Error tolerance (default = 1e-6).
    
    RETURNS
    ----------
    intersect : ndarray
        Array which triangles are intersected by the ray from each ray_origin.
    intersect_points : ndarray, optional
        The points where each ray crosses the plane of a particular triangle.
    intersect_distances : ndarray, optional
        The distance from ray_origin to the point of interscetion with the
        plane spanned by each triangle.
        
    NOTES
    ----------
    
    """
    # Check inputs
    ray_direction = np.array(ray_direction)
    if ray_direction.shape == ray_origin.shape:
        ray_direction = ray_direction[:,None,:]
    elif ray_direction.size == 3: # broadcast to all ray origins     
        ray_direction = ray_direction[None,None,:]
    else:
        raise ValueError("ray_direction must be either a single vector or have same dimensions as ray_origin.")
    
    try:
        assert plane_type in ["one-sided", "two-sided"]
    except AssertionError:
        raise ValueError("Choose 'one-sided' or 'two-sided'.")
    if posint is None:
        posint = []
    if posint_ok is None:
        posint_ok = []
        
    if (len(posint) > 0) & isinstance(posint,list):
        posint, posint_ok = list2numpy(posint, dtype=np.int)
        
    # Get vectors (e1, e2) spanning each triangle from point v0
    v0 = triangles[:,0,:]
    e1 = triangles[:,1,:]-v0
    e2 = triangles[:,2,:]-v0
    if len(posint)>0:
        try:
            assert posint.shape == posint_ok.shape
        except AssertionError:
            raise ValueError("posint and posint_ok must have same dimensions!")
    
        v0=v0[posint]
        e1=e1[posint]
        e2=e2[posint]
    else:
        v0=v0[None,...]
        e1=e1[None,...]
        e2=e2[None,...]
    
    # if not specified, consider all intersections possible
    if len(posint_ok) == 0:
        posint_ok = True
        
    # RAY TESTING
    # =======================
    # Implementation of the algorithm presented in
    # 
    # Moeller & Trumbore (1997). Fast, Minimum Storage Ray/Triangle
    # Intersection. Journal of Graphics Tools, 2(1):21--28.
    
    # Vector perpendicular to e2 and ray_direction
    #   P = CROSS(dir, e2)
    P = np.cross(ray_direction, e2)

    # Determinant
    #   D = DOT(P, e1)
    # D describes the relationship between triangle normal (face of the
    # triangle) and ray direction: if det>0 then the ray points towards the
    # outer face of the triangle (i.e. the ray and the triangle normal point
    # towards each other); if det<0, they point in the same direction; if
    # det=0, the ray and triangle are parallel.
    det = np.einsum("ijk,ijk->ij", P, e1)    
    inv_det = 1./(det+1e-10)
    
    # Vector from v0 to ray_origin
    #   T = O-V0
    tvec = ray_origin[:,None,:] - v0

    # 1st barycentric (e.g., x) coordinate of intersection point, i.e. where 
    # the ray intersects the plane spanned by e1 and e2.
    #   u = DOT(T, P)
    u = np.einsum("ijk,ijk->ij", P, tvec)*inv_det
    
    Q = np.cross(tvec, e1)
    
    # 2nd barycentric (e.g., y) coordinate of intersection point
    #   v = DOT(dir, Q)
    v = np.einsum("ijk,ijk->ij", ray_direction, Q)*inv_det
    
    # Distance from ray_origin to point of intersection in units of
    # ray_direction
    t = np.einsum("ijk,ijk->ij", e2, Q)*inv_det

    # CHECK CONDITIONS
    # =======================
    
    # Check that the ray crosses the plane spanned by vectors E1 and E2 in the
    # direction of the ray (one-sided) or if it crosses on either side ()
    if plane_type == "one-sided":
        posint_ok = posint_ok & (det >= eps)
    elif plane_type == "two-sided":
        posint_ok = posint_ok & (np.abs(det) >= eps)
    
    # Check that the ray actually crosses a triangle by testing if the
    # barycentric coordinates fulfills the equation
    #   T(u,v) = (1-u-v)*v0+u*v1+v*v2
    # which they do when u>=0, v>=0, (u+v)<=1.
    # Additionally, check if the intersection is within the allowed distance
    # (maxdist2int) by testing t.
    posint_ok = posint_ok & (u>=-eps) & (u<=1+eps)
    posint_ok = posint_ok & (v>=-eps) & ((u+v)<=1+eps)
    intersect = posint_ok & (t>(mindist2int+eps)) & (t<(maxdist2int+eps))

    t[~posint_ok]=np.inf
    if return_int_points:
        # The points where each ray cross the plane of each triangle
        intersect_points = ray_origin[:,np.newaxis,:]+t[...,np.newaxis]*ray_direction        
        if return_int_dists:      
            return intersect, intersect_points, t
        else:
            return intersect, intersect_points
    elif return_int_dists:
        return intersect, t
    else:
        return intersect



def mesh_load(fname):
    """Load a surface mesh in either .off or .stl file format. Return the 
    vertices and faces of the mesh. If .stl file, assumes only one solid, i.e. 
    only one mesh per file.
    
    PARAMETERS
    ----------
    fname : str 
        Name of the file to be read (.off or .stl file).
    
    RETURNS
    ----------
    vertices : ndarray
        Triangle vertices.
    faces : ndarray
        Triangle faces (indices into "vertices").
    """
    file_format = os.path.splitext(fname)[1].lower()
    
    if file_format == ".off":
        with open(fname, "rb") as f:
            # Read header
            hdr = f.readline().decode().rstrip("\n").lower()
            assert hdr == "off", ".off files should start with OFF"
            while hdr.lower() == "off" or hdr[0] == "#" or hdr == "\n":
                hdr = f.readline().decode()
            hdr = [int(i) for i in hdr.split()]
            
            # Now read the data
            vertices = np.genfromtxt(islice(f,0,hdr[0]))
            faces    = np.genfromtxt(islice(f,0,hdr[1]),
                                     usecols=(1,2,3)).astype(np.uint)
        return vertices,faces
        
    elif file_format == ".stl":
        # test if ascii. If not, assume binary        
        with open(fname, "rb") as f:
            try:
                if f.readline().decode().split()[0] == "solid":
                    is_binary = False
                else:
                    is_binary = True
            except:
                is_binary = True
                
        if is_binary:
            with open(fname, "rb") as f:
                # Skip the header (80 bytes), read number of triangles (1
                # byte). The rest is the data.
                np.fromfile(f, dtype=np.uint8, count=80)         
                np.fromfile(f, dtype=np.uint32, count=1)[0]
                data = np.fromfile(f, dtype=np.uint16, count=-1)
            data = data.reshape((-1,25))[:,:24].copy().view(np.float32)
            vertices = data[:,3:].reshape(-1,3) # discard the triangle normals
            
        else:
            vertices = []
            with open(fname,"rb") as f:
                for line in f:
                    line = line.decode().lstrip().split()
                    if line[0] == "vertex":
                        vertices.append(line[1:])
            vertices = np.array(vertices, dtype=np.float)
        
        # The stl format does not contain information about the faces, hence we
        # will need to figure this out.
        faces = np.arange(len(vertices)).reshape(-1,3)

        # Remove vertice duplicates and sort rows by sum    
        sv = np.sum(vertices+vertices*(100*np.random.random(3))[np.newaxis,:],
                    axis=1)
        sv_arg = np.argsort(sv)
        sv_arg_rev = np.argsort(sv_arg) # reverse indexing for going back
        
        # Get unique rows, indices of these, and counts. Create the new indices
        # and repeat them
        u, u_idx, u_count = np.unique(sv[sv_arg],return_index=True,
                                      return_counts=True)
        repeat_idx = np.repeat(np.arange(len(u)), u_count)
        
        # Retain only unique vertices and modify faces accordingly
        vertices = vertices[sv_arg][u_idx]
        faces = repeat_idx[sv_arg_rev][faces]
        
        return vertices, faces
            
    else:
        raise IOError("Invalid file format. Only files of type .off and .stl are supported.")

        
      
def mesh_save(vertices,faces,fname,file_format="off",binary=True):
    """Save a surface mesh described by points in space (vertices) and indices
    into this array (faces) to an .off or .stl file.
    
    PARAMETERS
    ----------
    vertices : ndarray
        Array of vertices in the mesh.
    faces : ndarray, int
        Array describing the faces of each triangle in the mesh.
    fname : str
        Output filename.
    file_format : str, optional
        Output file format. Choose between "off" and "stl" (default = "off").
    binary : bool
        Only used when file_format="stl". Whether to save file as binary (or
        ascii) (default = True).

    RETURNS
    ----------
    Nothing, saves the surface mesh to disk.
    """
    nFaces = len(faces)
    file_format = file_format.lower()
    
    # if file format is specified in filename, use this
    if fname.split(".")[-1] in ["stl","off"]:
        file_format = fname.split(".")[-1].lower()
    else:
        fname = fname+"."+file_format
    
    if file_format == "off":
        nVertices = len(vertices)
        with open(fname, "wb") as f:
            f.write("OFF\n".encode())
            f.write("# File created by ... \n\n".encode())
            np.savetxt(f,np.array([nVertices,nFaces,0])[np.newaxis,:],fmt="%u")
            np.savetxt(f,vertices,fmt="%0.6f")
            np.savetxt(f,np.concatenate((np.repeat(faces.shape[1],nFaces)[:,np.newaxis],faces),axis=1).astype(np.uint),fmt="%u")
    
    elif file_format == "stl":
        mesh = vertices[faces]
        tnormals  = get_triangle_normals(mesh)
        data = np.concatenate((tnormals, np.reshape(mesh, [nFaces,9])),
                              axis=1).astype(np.float32)
        
        if binary:
            with open(fname, "wb") as f:
                f.write(np.zeros(80, dtype=np.uint8))
                f.write(np.uint32(nFaces))
                f.write(np.concatenate((data.astype(np.float32,order="C",copy=False).view(np.uint16),np.zeros((data.shape[0],1),dtype=np.uint16)),axis=1).reshape(-1).tobytes())
        else:
            with open(fname, "w") as f:
                f.write("solid MESH\n")
                for t in range(len(data)):
                    f.write(" facet normal {0} {1} {2}\n  outer loop\n   vertex {3} {4} {5}\n   vertex {6} {7} {8}\n   vertex {9} {10} {11}\n  endloop\n endfacet\n"\
                    .format(*data[t,:]))
                f.write("endsolid MESH\n")
    else:
        raise IOError("Invalid file format. Please choose off or stl.")



def read_curv(fname):
    """Read a FreeSurfer curvature file. This filetype associates one value
    with each node (e.g., gray matter thickness). Uses big-endian format.
    
    PARAMETERS
    ----------
    fname : str
        Name of the file to be read.
    
    RETURNS
    ----------
    data : ndarray
        The data in the file.
    """
    # if first three bytes equal this number, then it is the "new" curvature
    # format    
    NEW_VERSION_MAGIC_NUMBER = 16777215
    
    with open(fname,"rb") as f:
        n_vertices = read_int3(f)
        if n_vertices == NEW_VERSION_MAGIC_NUMBER:
            n_vertices, n_faces, vals_per_vertex = np.fromfile(f,">i4",3)
            data = np.fromfile(f, ">f4", n_vertices)
        else:
            read_int3(f) # number of faces
            data = np.fromfile(f, ">i4", n_vertices)/100.
    
    return data



def read_int3(f):
    """Read three bytes from the current position in the file as an integer,
    i.e. int24.
    
    PARAMETERS
    ----------
    f : file object
        The file object from which to read.
        
    RETURNS
    ----------
    val : int
        The integer read from three bytes in the file.
    """
    b1 = struct.unpack('B',f.read(1))[0]
    b2 = struct.unpack('B',f.read(1))[0]
    b3 = struct.unpack('B',f.read(1))[0]

    val = (b1 << 16) + (b2 << 8) + b3
    
    return val
  


def make_handler(handler_type="stream", level=logging.INFO, filename=None):
    """Create a stream handler (logging to console) or a file handler (logging
    to a file) for a logger.

    PARAMETERS
    ----------
    handler_type : str
        The handler type, either "stream" or "file" (default = "stream").
    level : logging level
        Logging level (default = logging.INFO)
    filename : str
        The name of the file to which to write logs (only applicable if
        handler_type is "file")
        
    RETURNS
    ---------
    handler : handler object
        Handler object which can be added to a logger.
    """
    
    if handler_type == "stream":
        formatter = logging.Formatter("%(message)s")
        handler = logging.StreamHandler()

    if handler_type == "file":
        formatter = logging.Formatter(fmt="%(asctime)s [%(levelname)-5.5s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        if filename == None:
            filename = os.path.join(os.getcwd(), "log.log")
        handler = logging.FileHandler(filename, mode='a')

    handler.setLevel(level)
    handler.setFormatter(formatter)

    return handler



def spawn_process(cmd, return_exit_status=False, new_thread=False,
                  verbose=False, shell=False, new_process=False):
    """Spawn a new process and communicate its output.
    
    PARAMETERS
    ----------
    cmd : str
        The command to be executed.
    return_exit_status : bool, optional
        Return exit status of process. Only applies if new_thread == False
        (default = False).
    new_thread : bool, optional
        By default, the child process blocks the main process while running,
        however, by setting this option to true the child progress will be
        spawned in a new thread, hence not blocking the main process (default =
        False).
    verbose : bool, optional
        Communicate stdout to console or only to logfile. Only applies if 
        new_thread == False (default = False).
    shell: bool, optional
        Whether to use shell mode. Default: False (forced to be True on Windows)
    new_process: bool, optional
        Starts command in a new process. Default: False
    RETURNS
    ----------
    p.returncode : int (optional)
        Return code of the process.
    """
    log("Running {0}",cmd, level=logging.DEBUG)
    if verbose:
        lvl = logging.INFO
    else:
        lvl = logging.DEBUG

    if new_thread or new_process:
        ON_POSIX = "posix" in sys.builtin_module_names

        def enqueue_output(out, queue):
            for line in iter(out.readline, b''):
                try:
                    queue.put(line.decode())
                except UnicodeDecodeError:
                    queue.put('Could not print line')
            out.close()

        p = Popen(cmd, stdout=PIPE,
                  bufsize=1, close_fds=ON_POSIX,
                  shell=True)

        q = Queue()
        if new_thread:
            t = Thread(target=enqueue_output, args=(p.stdout, q))
        if new_process:
            t = Process(target=enqueue_output, args=(p.stdout, q))
        t.daemon = True  # thread dies with the program
        t.start()
        return t

    else: 
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        for line in iter(p.stdout.readline, b''):
            # to prevent massive output to the logfile retain only the last
            # line starting with \r
            try:
                line = line.decode().split("\r")[-1]
            except UnicodeDecodeError:
                log('Could not print line', level=logging.DEBUG)
                continue
            
            # remove one \n since logger.log adds one itself
            if line[-1]=="\n":
                line = line[:-1]
            
            try:
                log(line, level=lvl)
            except KeyError: # in case line includes ".. {..} .."
                log(line.replace("{","{{").replace("}","}}"), level=lvl)
                
        if return_exit_status:
            return p.wait()



def log(msg, args=None, level=logging.INFO, width=72):
    """Log message with formatting.
    
    PARAMETERS
    ----------
    msg : str
        The message to log.
    args : list or str, optional
        List of arguments. Arguments are applied using .format, hence supplying
        arguments, msg should contain proper placeholders for these arguments
        (e.g., {0}, {1}, {:d} or with options such as {:06.2f}) (default = [],
        i.e. no arguments).        
    width : int, optional
        Wrapping width of msg (default = 72)
        GUILHERME: deprecated this argument, it only makes the log more confusing
    level : logging level, optional
        Specify logging level (default = logging.INFO)
        
    RETURNS
    ----------
    The message logged on logger named "logger".    
    """
    if args is None:
        args = []
    if type(args) not in [list,tuple]:
        args = [args]
        
    try:
        logger.log(level, msg.format(*args))
    # If there is an error in the message formating, logs the raw message
    except ValueError:
        logger.log(level, msg)
        

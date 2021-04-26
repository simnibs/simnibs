import os
import tempfile
import logging
import numpy as np
import scipy.sparse
import scipy.ndimage
import time
import numba

from . import mesh_io
from . import cgal
from ..utils.simnibs_logger import logger, format_time
from ..segmentation.brain_surface import dilate, erosion
from ..segmentation._thickness import _calc_thickness
from ..utils.transformations import get_vox_size


class MeshingError(ValueError):
    pass


def _write_inr(image, voxel_dims, fn_out):
    if image.ndim != 3:
        raise MeshingError('Nifti volume must have 3 dimensions')
    if len(voxel_dims) != 3:
        raise MeshingError('Voxel_dims must have 3 dimensions')

    # Figure out data type and translate it to inria data types
    if image.dtype in [np.uint8, np.uint16]:
        inria_dtype = 'unsigned fixed'
    elif image.dtype in [np.float32, np.float64]:
        inria_dtype = 'float'
    else:
        raise MeshingError(
            f'Image data type: {image.dtype} not compatible. '
            f'compatible types are: uint8, uint16, float32, float64')

    n_bits = 8 * image.itemsize
    with open(fn_out, 'bw') as f:
        header = b'#INRIMAGE-4#{\n'
        header += 'XDIM={0}\nYDIM={1}\nZDIM={2}\nVDIM=1\n'.format(*image.shape).encode()
        header += f'TYPE={inria_dtype}\nPIXSIZE={n_bits} bits\nCPU=decm\n'.encode()
        header += 'VX={0}\nVY={1}\nVZ={2}\n'.format(*voxel_dims).encode()
        # fill out remaining part of the header
        header += (252-len(header)) * b'\n' + b'##}\n'
        f.write(header)
        f.write(image.tobytes(order='F'))


def _mesh_image(image, voxel_dims, facet_angle,
                facet_size, facet_distance,
                cell_radius_edge_ratio, cell_size,
                optimize):

    with tempfile.TemporaryDirectory() as tmpdir:
        fn_image = os.path.join(tmpdir, 'image.inr')
        fn_mesh = os.path.join(tmpdir, 'mesh.mesh')
        _write_inr(image, voxel_dims, fn_image)
        if type(cell_size) is np.ndarray:
            ret = cgal.mesh_image_sizing_field(
                    fn_image.encode(), fn_mesh.encode(),
                    facet_angle, facet_size, facet_distance,
                    cell_radius_edge_ratio, cell_size,
                    optimize
                 )
        else:
            ret = cgal.mesh_image(
                    fn_image.encode(), fn_mesh.encode(),
                    facet_angle, facet_size, facet_distance,
                    cell_radius_edge_ratio, cell_size,
                    optimize
                 )

        if ret != 0:
            raise MeshingError('There was an error while meshing the image')

        mesh = mesh_io.read_medit(fn_mesh)
    # In concurrent meshing there might be some spurious nodes
    used_nodes = np.unique(mesh.elm[:])[1:]
    mesh = mesh.crop_mesh(nodes=used_nodes)

    return mesh

def _decompose_affine(affine):
    ''' Decompose affine transformation into the form A = RZS
    Where R is a rotation matriz, Z a scaling matrix and S a shearing matrix
    '''
    if np.isclose(np.linalg.det(affine), 0):
        raise ValueError('Affine matrix is singular!')
    rot, z = np.linalg.qr(affine[:3, :3])
    scaling = np.diagonal(z)
    rot *= np.sign(scaling)
    z *= np.sign(scaling)[:, None]
    scaling = np.abs(scaling)
    shearing = np.diag(1/scaling).dot(z)
    return rot, scaling, shearing


def _resample2iso(image, affine, sampling_rate=1, order=1):
    ## DEPRECATED ##
    # R is a rotation matrix, K a scaling matrix and S a shearing matrix
    R, Z, S = _decompose_affine(affine[:3, :3])
    # Definew new affine matrix
    new_res = sampling_rate * np.min(Z)
    new_affine = R.dot(new_res * np.eye(3))
    # Transformation from original image to new image
    M = affine[:3, :3].dot(np.linalg.inv(new_affine))
    # Add translation component
    aff = np.eye(4)
    aff[:3, :3] = new_affine
    aff[:3, 3] = affine[:3, 3]
    new_affine = aff
    if np.allclose(M, np.eye(3)):
        return image.copy(), new_affine
    # New image shape
    target_shape = np.ceil(M.dot(image.shape)).astype(int)
    # Resample
    image_resampled = scipy.ndimage.affine_transform(
        image,
        np.linalg.inv(M),
        output_shape=target_shape,
        output=image.dtype,
        order=order
    )
    return image_resampled, new_affine

def image2mesh(image, affine, facet_angle=30,
               facet_size=None, facet_distance=None,
               cell_radius_edge_ratio=3, cell_size=None,
               optimize=True):
    ''' Creates a mesh from a 3D image

    Parameters
    ------------
    image: ndarray
        3D ndarray corresponding to time image. Must be of type uint8 or uint16

    affine: 4x4 ndarray
        Array describing the affine transformation from voxel space to world space. Must
        be decomposable in the format (R*Z) where R is a rotation matrix and Z a diagonal
        scaling matrix (shearing not accepted).

    facet_angle: float (optional)
        This parameter controls the shape of surface facets. Specifically, it is a lower
        bound for the angle (in degrees) of surface facets. When boundary surfaces are
        smooth, the termination of the meshing process is guaranteed if this angular bound is
        at most 30 degrees. Default: 30

    facet_size: float or ndarray wih same shape as image(optional)
        This parameter controls the size of surface facets. Each surface facet has a surface
        Delaunay ball which is a ball circumscribing the surface facet and centered on the
        surface patch. The parameter facet_size is either a constant or a spatially variable
        scalar field, providing an upper bound for the radii of surface Delaunay balls.
        Default: minimum voxel size (very low!)

    facet_distance: float or ndarray with same shape as image(optional)
        This parameter controls the approximation error of boundary and subdivision surfaces.
        Specifically, it is either a constant or a spatially variable scalar field. It
        provides an upper bound for the distance between the circumcenter of a surface facet
        and the center of a surface Delaunay ball of this facet. Default: minimum voxel size

    cell_radius_edge_ratio: float (optional)
        This parameter controls the shape of mesh cells (but can't filter slivers, as we
        discussed earlier). It is an upper bound for the ratio between the circumradius of a
        mesh tetrahedron and its shortest edge. There is a theoretical bound for this
        parameter: the Delaunay refinement process is guaranteed to terminate for values of
        cell_radius_edge_ratio bigger than 2. Default: 3

    cell_size: float or ndarray with same shape as image(optional)
        This parameter controls the size of mesh tetrahedra. It is either a scalar or a
        spatially variable scalar field. It provides an upper bound on the circumradii of the
        mesh tetrahedra. Default: minimum voxel size (very low!)

    optimize: bool (optional)
        Tunrn on Lloyd optimization. Sliver perturbation and exudation is always done. Default: True

    Returns
    ----------
    msh: simnibs.Msh
        Mesh structure

    References
    -----------
    https://doc.cgal.org/latest/Mesh_3/index.html
    '''

    if image.dtype not in [np.uint8, np.uint16]:
        raise MeshingError('Image must be of type uint8 or uint16')

    rot, voxel_dims, shearing = _decompose_affine(affine)

    if not np.allclose(shearing, np.eye(3), atol=1.e-5):
        raise ValueError('Affine matrix has a shearing component')

    if facet_size is None:
        facet_size = min(voxel_dims)

    if facet_distance is None:
        facet_distance = min(voxel_dims)

    if cell_size is None:
        cell_size = min(voxel_dims)

    if type(cell_size) is np.ndarray or\
       type(facet_size) is np.ndarray or \
       type(facet_distance) is np.ndarray:

        if type(cell_size) is np.ndarray:
            assert cell_size.shape == image.shape
        else:
            cell_size = cell_size * np.ones_like(image, dtype=np.float32, order='F') 

        if type(facet_size) is np.ndarray:
            assert facet_size.shape == facet_size.shape
        else:
            facet_size = facet_size * np.ones_like(image, dtype=np.float32, order='F') 

        if type(facet_distance) is np.ndarray:
            assert facet_distance.shape == facet_distance.shape
        else:
            facet_distance = facet_distance * np.ones_like(image, dtype=np.float32, order='F') 

    mesh = _mesh_image(
        image, voxel_dims,
        facet_angle, facet_size, facet_distance,
        cell_radius_edge_ratio, cell_size,
        optimize
    )
    # Rotate nodes
    mesh.nodes.node_coord = rot.dot(mesh.nodes.node_coord.T).T
    # Translate nodes
    mesh.nodes.node_coord += affine[:3, 3]
    # Fix node orderings
    if np.linalg.det(rot) < 0:
        mesh.elm.node_number_list = mesh.elm.node_number_list[:, [1, 0, 2, 3]]
    return mesh


def _mesh_surfaces(surfaces, subdomains, facet_angle,
                   facet_size, facet_distance,
                   cell_radius_edge_ratio, cell_size,
                   optimize):
    if len(surfaces) != len(subdomains):
        raise MeshingError('Please define one surface per subdomain')

    with tempfile.TemporaryDirectory() as tmpdir:
        fn_surfaces = []
        sd_formated = []
        for surf, sd in zip(surfaces, subdomains):
            if len(sd) != 2:
                raise MeshingError('Subdomain speficifation must have exactly 2 entries')
            fn = os.path.join(tmpdir, f'{sd[0]:03d}-{sd[1]:03d}.off')
            mesh_io.write_off(surf, fn)
            fn_surfaces.append(fn.encode())
            sd_formated.append((sd[0], sd[1]))
        fn_mesh = os.path.join(tmpdir, 'mesh.mesh')
        ret = cgal.mesh_surfaces(
                fn_surfaces, sd_formated, fn_mesh.encode(),
                facet_angle, facet_size, facet_distance,
                cell_radius_edge_ratio, cell_size,
                optimize
               )
        if ret != 0:
            raise MeshingError('There was an error while meshing the surfaces')
        mesh = mesh_io.read_medit(fn_mesh)

    # In concurrent meshing there might be some spurious nodes
    used_nodes = np.unique(mesh.elm[:])[1:]
    mesh = mesh.crop_mesh(nodes=used_nodes)

    return mesh


def remesh(mesh, facet_size, cell_size,
           facet_angle=30, facet_distance=0.1,
           cell_radius_edge_ratio=2, optimize=True):
    ''' Extracts the mesh subdomains and recreate the mesh usign CGAL

    Parameters
    ------------
    msh: simnibs.Msh
        Mesh structure, to be reconstructed

    facet_size: float
        See https://doc.cgal.org/latest/Mesh_3/index.html#title10.

    cell_size: float (optional)
        See https://doc.cgal.org/latest/Mesh_3/index.html#title10.

    facet_angle: float (optional)
        See https://doc.cgal.org/latest/Mesh_3/index.html#title10. Default: 30

    facet_distance: float (optional)
        See https://doc.cgal.org/latest/Mesh_3/index.html#title10. Default: 0.1

    cell_radius_edge_ratio: float (optional)
        See https://doc.cgal.org/latest/Mesh_3/index.html#title10. Default: 2

    Returns
    ------------
    msh: simnibs.Msh
        Reconstructed mesh
    '''
    # Reconstruct surfaces from tetrahedra
    mesh = mesh.crop_mesh(elm_type=4)
    mesh.reconstruct_surfaces()
    mesh.elm.tag1[mesh.elm.tag1 > 1000] -= 1000
    # Find the tetrahedra adjacent to each surface
    adj_th = mesh.elm.find_adjacent_tetrahedra()
    adj_labels = mesh.elm.tag1[adj_th - 1]
    adj_labels[adj_th == -1] = 0
    triangles = mesh.elm.elm_type == 2
    adj_labels = adj_labels[triangles, :2]
    tr_labels = mesh.elm.tag1[triangles]
    # fix the label order so that the first is the triangle label
    # and the second the label of the adjacent tissue
    adj_labels_fixed = np.zeros_like(adj_labels)
    adj_labels_fixed[:, 0] = tr_labels
    for i in range(2):
        adj_labels_fixed[
            adj_labels[:, i] != tr_labels, 1
        ] = adj_labels[adj_labels[:, i] != tr_labels, i]
    adj_labels = adj_labels_fixed
    # Go throu all the
    tags = np.unique(tr_labels)
    tags = np.append(tags, 0)
    surfaces = []
    subdomains = []
    for i, t1 in enumerate(tags):
        for j, t2 in enumerate(tags[i+1:]):
            tr_interface = \
                (adj_labels[:, 0] == t1) *\
                (adj_labels[:, 1] == t2)
            if np.any(tr_interface):
                surfaces.append(
                    mesh.crop_mesh(elements=mesh.elm.triangles[tr_interface])
                )
                subdomains.append((t1, t2))

    mesh = _mesh_surfaces(
        surfaces, subdomains, facet_angle,
        facet_size, facet_distance,
        cell_radius_edge_ratio, cell_size,
        optimize
    )
    return mesh

def relabel_spikes(elm, tag, with_labels, adj_labels, label_a, label_b, 
                   target_label, labels, nodes_label, adj_th, adj_threshold=2,
                   log_level=logging.DEBUG, relabel_tol=1e-6, max_iter=20,
                   nlist=None,maxn=None):
    ''' Relabels the spikes in a mesh volume, in-place

    A spike is defined as a tetrahedron in "label_a" or "label_b"
    which has at least one node in the other volume and
    at least "adj_threshold" faces adjacent to tetrahedra in
    "target_label".

    Parameters
    -----------
    elm: ndarray
       simnibs.Msh.elm[:] mesh structure
    tag: ndarray
        labels for elements (simnibs.Msh.elm.tag1)
    with_labels: ndarray (ntag x nelements) bool
        indicates if element contains each of the labels
    adj_labels:  ndarray int
        labels for adjacent elements
    label_a: int
        index for volume label with spikes
    label_b: int
        index for second volume label with spikes
    target_label: int
        index for volume label to relabel spikes to
    labels: ndarray int
        list of labels
    nodes_label: ndarray (ntag x nelements) int
        count of how many nodes in each element have each label, used when updating
    adj_th: list
       value of m.elm.find_adjacent_tetrahedra()
    adj_threshold: int (optional)
        Threshhold of number of adjacent faces for being considered a spike
    relabel_tol: float (optional)
        Fraction of the elements that indicates convergence
    max_iter: int
        Maximum number of relabeling iterations
    nlist : ndarray int
        list of which elements each nodes is connected to
    maxn : ndarray int
        number of elements that each nodes is connected to 
        (needed for lookup in nlist)
    '''
    logger.log(
        log_level,
        f'Relabeling spikes in {labels[label_a]} and {labels[label_b]} to {labels[target_label]}'
    )
    if not np.any(with_labels[label_a] * with_labels[label_b]):
        return

    for i in range(max_iter):
        # Relabel tissue A
        # Find spikes
        A_to_relabel, frac_A_relabeled = _find_spikes(tag, label_a, label_b,
                  with_labels, adj_labels, target_label, labels, adj_threshold)
        # Update tags and adjlabels, with_labels and nodes_label in place
        _update_tags(tag, elm, adj_th, with_labels, adj_labels, A_to_relabel,
                     label_a, target_label, labels, nodes_label,nlist,maxn)

        # Relabel tissue B
        # Find spikes
        B_to_relabel, frac_B_relabeled = _find_spikes(tag, label_b, label_a,
                  with_labels, adj_labels, target_label, labels, adj_threshold)
        # Update tags and adjlabels, with_labels and nodes_label in place
        _update_tags(tag, elm, adj_th, with_labels, adj_labels, B_to_relabel,
                     label_b, target_label, labels, nodes_label,nlist,maxn)
        
        logger.log(log_level,
                   f'Relabeled {np.sum(A_to_relabel)} from {labels[label_a]} '
                   f'and {np.sum(B_to_relabel)} from {labels[label_b]}'
                   )

        # Stop if converge has been reached
        if frac_A_relabeled < relabel_tol and frac_B_relabeled < relabel_tol:
            break
        
    # A_to_relabel, frac_A_relabeled = _find_spikes(tag, label_a, label_b,
    #           with_labels, adj_labels, target_label, labels, adj_threshold)
    # B_to_relabel, frac_B_relabeled = _find_spikes(tag, label_b, label_a,
    #           with_labels, adj_labels, target_label, labels, adj_threshold)
    # print(f'converged after {i+1} iterations, {np.sum(A_to_relabel)+np.sum(B_to_relabel)} left to relabel')
        
@numba.njit(parallel=True, fastmath=True)
def _find_spikes(tag, label, label2, with_labels, adj_labels, target_label, labels, adj_threshold=2):
    '''
    Find spikes
    
    Parameters
    ----------
    tag : ndarray int
        labels for elements
    label : int
            index for volume label with spikes
    label2 : int
            index for second volume label with spikes
    with_labels : ndarray (ntag x nelements) bool
            indicates if element contains each of the labels
    adj_labels : ndarray int
        labels for adjacent elements
    target_label : int
        target label index
    labels : ndarray int
            list of labels.
    adj_threshold : list, optional
        Threshold of number of adjacent faces for being considered a spike. The default is 2.
    
    Returns
    -------
    ndarray bool
        Indicates if spikes were found
    
    '''
    # Initialize output
    found_spikes = np.zeros(tag.shape[0], dtype='bool')
    # initialize number of elements with relevant label
    na = 0
    nspikes = 0
    # Parallel loop over elements
    for i in numba.prange(tag.shape[0]):
    # if element has the label
        if tag[i] == labels[label]:
            # increment element count with label
            na += 1
            # if we have the other label too
            if with_labels[label2, i]:
                # local variable (for thread) holding spikes count
                spikes = 0
                # loop over adjacent element
                for j in range(adj_labels.shape[1]):
                    # if they have the target label too
                    if adj_labels[i, j] == labels[target_label]:
                        # increment spike count
                        spikes += 1
                        # if above threshold indicate spikes
                        if spikes >= adj_threshold:
                            found_spikes[i] = True
                            # increment number of spikes found
                            nspikes += 1
                            # it is already a spike no need to go on
                            break
    # return found spikes and ratio of spikes (to check convergence)
    return found_spikes, float(nspikes) / float(na)

@numba.njit
def _update_tags(tag, elm, adj_th, with_labels, adj_labels, to_relabel, label, 
                 target_label, labels, nodes_label,nlist,maxn):
    '''
    Update attributes needed for identifying spikes in place
    Parameters
    ----------
    tag : ndarray int
        labels for elements
    elm : ndarray int
        elements
    adj_th : ndarray int
        value of m.elm.find_adjacent_tetrahedra()
    with_labels : ndarray bool
    
    adj_labels : ndarray  int
        labels for adjacent elements
    to_relabel : ndarray bool
        Elements to relabel.
    label : int
        original label index
    target_label : int
        new label
    labels : ndarray int
        list of labels
    nodes_label : ndarray (ntag x nelements) int
        count of how many nodes in each element have each label
    nlist : ndarray int
            list of which elements each nodes is connected to
    maxn : ndarray int
        number of elements that each nodes is connected to 
        (needed for lookup in nlist)
    Returns
    -------
    None.
    
    '''
    # Loop over elements
    for i in range(tag.shape[0]):
        # relabel to intended label if indicated
        if to_relabel[i]:
            tag[i] = labels[target_label]
            # update count of nodes with labels
            for j in range(elm.shape[1]):
                nodes_label[label, elm[i, j]] -= 1  # decrease original label
                nodes_label[target_label, elm[i, j]] += 1  # increase intended label
                # loop over elements this node is connected to update with_labels
                # avoids looping over the entire mesh
                for m in range(maxn[elm[i,j]-1],maxn[elm[i,j]]):
                    # Update only relevant if the element actually had the node
                    if with_labels[label, nlist[m]]:
                        with_labels[label, nlist[m]] = False
                        for n in range(elm.shape[1]):
                            if nodes_label[label, elm[nlist[m], n]] > 0:
                                with_labels[label, nlist[m]] = True
                                # no need to go on it is already True
                                break
                    # Update can also be relevant if element does not have target label
                    if not with_labels[target_label, nlist[m]]:
                        for n in range(elm.shape[1]):
                            # Update target label if nodes has it
                            if nodes_label[target_label, elm[nlist[m], n]] > 0:
                                with_labels[target_label, nlist[m]] = True
                                # no need to go on it is already True
                                break
            # loop over neighboors (4)
            for k in range(adj_th.shape[1]):
                adj_i = adj_th[i, k] - 1  # indexing is from 1, 0 indicate no neighbor
                if adj_i >= 0:  # only if neighbor
                    for j in range(adj_th.shape[1]):
                        # if this element is the one being updated
                        if adj_th[adj_i, j] - 1 == i:
                            # update adjacent label
                            adj_labels[adj_i, j] = labels[target_label]
                            # no need to test more there is only one
                            break

@numba.njit
def _with_label_numba_all(elm, tag1, labels, N):
    ''' Returns a count of how many labels of each type belongs to each node
        and all elements in the mesh which have a node that is in
        a region with each label
    '''
    # output for counting label types for each node indexing starts from 1
    # so first element is empty.
    nodes_label = np.zeros((len(labels), N+1), dtype='uint16')
    # output for boolean array indicating if element touches each node type
    with_labels = np.zeros((len(labels), elm.shape[0]), dtype='bool')
    maxn = np.zeros((N+1), dtype=elm.dtype)
    # loop over labels
    for k in range(len(labels)):
        # Loop over elements
        for i in range(elm.shape[0]):
            # increment count if elements has current label
            if tag1[i] == labels[k]:
                for j in range(elm.shape[1]):
                    nodes_label[k, elm[i, j]] += 1
                    maxn[elm[i,j]] += 1
    # cummulative count of how many elements is connected nodes
    # equivalent to np.cumsum(maxn) but in-place
    for i in range(2,N+1):
        maxn[i] += maxn[i-1]
    # create an array containing the element number that each node is connected to
    nlist = np.zeros((elm.shape[0]*elm.shape[1]),dtype=elm.dtype)
    # list for counting how many elements a nodes has currently been asigned to
    ncount = np.zeros((N+1,),dtype='uint16')
    for i in range(elm.shape[0]):
        for j in range(elm.shape[1]):
            # element index (starting from 0 here)
            n = elm[i,j]-1
            # set the n'th element that the node is connected to 
            # ncount[n] counts how many has already been set
            nlist[maxn[n] + ncount[n]] = i
            #increment the elements that the node is connected to
            ncount[n] += 1
    # Loop over labels
    for k in range(len(labels)):
        # Loop over elements
        for i in range(elm.shape[0]):
        # Loop over nodes within label (4)
            for j in range(elm.shape[1]):
                # Indicate if element contains label
                if nodes_label[k, elm[i, j]] > 0:
                    with_labels[k, i] = True
                    break
    return nodes_label, with_labels, nlist, maxn

def despike(msh, adj_threshold=2, relabel_tol=1e-6, max_iter=20,
            log_level=logging.DEBUG):

    ''' Goes through the mesh removing spiles
    A spike is defined as a tetrahedron in a volume "a"
    which has at least one node in the other volume "b" and
    at least "adj_threshold" faces adjacent to tetrahedra in a volume "c"

    Parameters
    -----------
    m: simnibs.Msh
       Mesh structure
    adj_threshold: int (optional)
        Threshhold of number of adjacent faces for being considered a spike
    relabel_tol: float (optional)
        Fraction of the elements that indicates convergence
    max_iter: int
        Maximum number of relabeling iterations
    '''
    
    if np.any(msh.elm.elm_type != 4):
        logger.log(log_level,
                   'Error: Attempting to despike mesh containing not only'
                   'tetrahedra. Please consider cropping the mesh first.' 
                   )
        raise ValueError()
        return

    tags = np.unique(msh.elm.tag1)
    adj_th = msh.elm.find_adjacent_tetrahedra()
    elm = msh.elm[:]
    tag = msh.elm.tag1
    adj_labels = tag[adj_th - 1]
    adj_labels[adj_th == -1] = -1

    # Total number of nodes
    N = msh.nodes.nr
    # Count how many labels of each type belongs to each node - first output
    # and determine if elements touch each labels (has a node with that label) - second output
    # and determine which element each node is connected too - third output
    nodes_label, with_labels, nlist, maxn = _with_label_numba_all(elm, tag, tags, N)
    # Loop over labels
    for i, t1 in enumerate(tags):
        for j, t2 in enumerate(tags[i + 1:]):
            # Only if the labels are not the same
            if t1 == t2:
                continue
            #only if at least one elements have this label combination
            if not np.any((nodes_label[i] > 0) * (nodes_label[j + i + 1] > 0)):
                continue
            for k, t3 in enumerate(tags):
                # Only if target label is different
                if t1 == t3 or t2 == t3:
                    continue
                #only if at least one elements have this label combination
                if not np.any((nodes_label[i] > 0) * 
                              (nodes_label[j + i + 1] > 0) * nodes_label[k]):
                    continue
                #call relabel function
                relabel_spikes(elm = elm, tag = tag, with_labels = with_labels,
                               adj_labels = adj_labels, label_a = i, 
                               label_b = j + i + 1, target_label = k,
                               relabel_tol = relabel_tol, labels = tags,
                               adj_threshold = adj_threshold, adj_th = adj_th,
                               max_iter = max_iter, log_level = log_level, 
                               nodes_label = nodes_label,nlist=nlist,maxn=maxn)
    #set tag1/tag2 in msh structure
    msh.elm.tag1 = tag
    msh.elm.tag2 = tag



def _sizing_field_from_thickness(label, thickness, elem_sizes):
    ''' Calculates a sizing field from thickness,
        with the option to define label-wise sizing rules'''
    # standard mapping
    slope = elem_sizes['standard']['slope']
    ranges = elem_sizes['standard']['range']
    if ranges[0] > ranges[1]:
        raise ValueError('Ranges value should be in format (min, max)')
    field = np.array(slope*thickness, dtype=np.float32, order='F')
    field[field < ranges[0]] = ranges[0]
    field[field > ranges[1]] = ranges[1]
    
    # label-specific mappings
    if len(elem_sizes)>1:
        tissues = list(elem_sizes)
        for i in range(len(tissues)):
            if tissues[i] == 'standard':
                continue
            slope = elem_sizes[tissues[i]]['slope']
            ranges = elem_sizes[tissues[i]]['range']
            if ranges[0] > ranges[1]:
                raise ValueError('Ranges value should be in format (min, max)')
            idx = label==int(tissues[i])
            field[idx] = slope*thickness[idx]
            idx2 = idx*(field < ranges[0])
            field[idx2] = ranges[0]
            idx2 = idx*(field > ranges[1])
            field[idx2] = ranges[1]
    return field


def create_mesh(label_img, affine, 
                elem_sizes={"standard": {"range": [1, 5], "slope": 1.0}},
                smooth_size_field = 2,
                skin_facet_size=2.0, 
                facet_distances={"standard": {"range": [0.1, 3], "slope": 0.5}},
                optimize=True, remove_spikes=True, skin_tag=1005,
                remove_twins=True, hierarchy=None, smooth_steps=5):
    """Create a mesh from a labeled image.

    The maximum element sizes (CGAL facet_size and cell_size) are controlled 
    by elem_sizes:
        size = slope * thickness
        size[size < range[0] = range[0]
        size[size > range[1] = range[1]
    where "thickness" is the local tissue thickness, 
    "range" is the size range (label-specific if label is added to elem_sizes, 
                                otherwise the "standard" range is used)
    The distance (CGAL facet_distance) parameter is calcualted in a similar way.
    This allows for the meshing to adjust sizes according to local needs.

    Parameters
    ----------
    label_img: 3D np.ndarray in uint8 format
        Labeled image from segmentation
    affine: 4x4 np.ndarray
        Affine transformation from voxel coordinates to world coordinates
    elem_sizes: dictionary (optional)
        Lists the relationship between thickness and elem_sizes.
        Label-specific relationships can be added if needed, e.g. for label 2:
            {"standard": {"range": [1, 5], "slope": 1.0},
                    "2": {"range": [1, 2], "slope": 0.7}}
        "range" determines the minimum and maximum values for element sizes.
        "slope" determines relationship between thickness and element sizes.      
        Note: Label indices are used as keys, and need to be provided as string
        Default: {"standard": {"range": [1, 5], "slope": 1.0}}
    smooth_size_field: int (optional)
        Defines the size of a triangular kernel to smooth the size field. A bit
        of smoothing helps to remove the effect of a few outliers in the 
        thickness estimates on the size field. Set to 0 to disable.
        The kernel size is 2*smooth_size_field+1 in voxels. Default:  2
    skin_facet_size: float (optional)
        Maximum size for the triangles of the outermost surface. If set to None,
        the elements of the outermost surface will be scaled the same way as
        the other surfaces. Default:  2.0
    facet_distances: dictionary (optional)
        Relationship between thickness and facet_distance. For small
        facet_distance values, the meshing will follow the label boundaries 
        in the original image more strictly. This also means more elements. 
        Label-specific relationships can be added if needed.
        "range": Minimum and maximum values for facet_distance
        "slope": Steepness of relationship between thickness and facet_distance.
        Default: {"standard": {"range": [0.1, 3], "slope": 0.5}}
    optimize: bool (optional)
        Whether to run lloyd optimization on the mesh. Default: True
    remove_spikes: bool (optional)
        Whether to remove spikes to create smoother meshes. Default: True
    skin_tag: float (optional)
        1) Restrict effects of skin_facet_size to the outer boundary of the 
           region with label skin_tag-1000 (i.e. 5 for 1005)
        2) Add outer surface to mesh using given tag. Set to None to disable.
        NOTE: This surface will replace any other surface with the same tag.
        Default: 1005
    remove_twins: bool (optional)
        Remove triangle twins created during surface reconstruction.
        Default: True
    hierarchy: list of ints or None (optional)
        List of surface tags that determines the order in which triangles 
        are kept for twin pairs (for remove_twins=True). Default for hierarchy=None:
        (1005, 1001, 1002, 1009, 1003, 1004, 1008, 1007, 1006, 1010)
        i.e. Skin (1005) has highest priority, WM (1001) comes next, etc; 
    smooth_steps: int (optional)
        Number of smoothing steps to apply to the final mesh surfaces. Default: 5

    Returns
    -------
    msh: simnibs.Msh
        Mesh structure
    """
    if hierarchy is None:
        hierarchy = (1005, 1001, 1002, 1009, 1003, 1004, 1008, 1007, 1006, 1010)
    if not 'standard' in elem_sizes:
        raise ValueError('elem_sizes needs a \"standard\" entry')
    if not 'standard' in facet_distances:
        raise ValueError('facet_distances needs a \"standard\" entry')

    # Calculate thickness
    logger.info('Calculating tissue thickness')
    start = time.time()
    thickness = _calc_thickness(label_img)
    thickness[thickness < .5] = 100 # set background thickness to some large value
    voxel_size = get_vox_size(affine) 
    if not np.allclose(np.diff(voxel_size), 0):
        logger.warn('Anisotropic image, meshing may contain extra artifacts')
    thickness *= np.average(voxel_size) # Scale thickness with voxel size
    
    # Define size fields and distance field
    logger.info('Calculating sizing fields')
    size_field = _sizing_field_from_thickness(
        label_img, thickness, elem_sizes
    )
    distance_field = _sizing_field_from_thickness(
        label_img, thickness, facet_distances
    )
    del thickness
    
    # Smooth size field a bit to reduce effect of a few outliers in the thickness
    # map on the mesh; the outliers show up as localized small thickness values at 
    # some of the tissue boundaries
    if smooth_size_field:
        size_field = size_field**(1/3) # compress high values to preserve edges a bit better
        kernel = smooth_size_field+1-np.abs(np.arange(-smooth_size_field, 
                                                      smooth_size_field+1, 1))
        kernel = kernel/np.sum(kernel)
        for i in range(3):
            size_field = scipy.ndimage.convolve1d(size_field, kernel, axis=i, 
                                            mode='constant', cval=0.0, origin=0)
        size_field = size_field**3
    
    # Control triangle size of outer surface to ensure eletrode meshing works OK
    if skin_facet_size is not None:      
        boundary = (label_img > 0).astype('int8')
        boundary = boundary-erosion(boundary,1)
        if skin_tag is not None:
            # keep boundary only at regions with label skin_tag-1000
            # to save some tetrahedra
            skin_tet_tag = skin_tag-1000
            boundary *= (label_img == skin_tet_tag)
            boundary = dilate(boundary,1)
            boundary *= (label_img == skin_tet_tag)
        else:
            boundary = dilate(boundary,1)
        size_field = size_field.flatten()
        size_field[boundary.flatten()] = skin_facet_size
        size_field = size_field.reshape(label_img.shape)
        del boundary
    logger.info(
        'Time to prepare meshing: ' +
        format_time(time.time()-start)
    )
    
    # Run meshing
    logger.info('Meshing')
    start = time.time()
    mesh = image2mesh(
        label_img,
        affine,
        facet_size=size_field,
        facet_distance=distance_field,
        cell_size=size_field,
        optimize=optimize
    )
    del size_field, distance_field
    logger.info(
        'Time to mesh: ' +
        format_time(time.time()-start)
    )
    
    # Separate out tetrahedron (will reconstruct triangles later)
    start = time.time()
    mesh = mesh.crop_mesh(elm_type=4)
    # Assign the right labels to the mesh as CGAL modifies them
    indices_seg = np.unique(label_img)[1:]
    new_tags = np.copy(mesh.elm.tag1)
    for i, t in enumerate(indices_seg):
        new_tags[mesh.elm.tag1 == i+1] = t
    mesh.elm.tag1 = new_tags
    mesh.elm.tag2 = new_tags.copy()
    
    # Remove spikes from mesh
    if remove_spikes:
        logger.info('Removing Spikes')
        despike(
            mesh, relabel_tol=1e-5,
            adj_threshold=2
        )
        
    # Reconctruct mesh surfaces
    logger.info('Reconstructing Surfaces')
    mesh.fix_th_node_ordering()
    idx=mesh.elm.connected_components()
    mesh = mesh.crop_mesh(elements=max(idx,key=np.size))
    mesh.reconstruct_surfaces(add_outer_as=skin_tag)
    if remove_twins:
        mesh = mesh.remove_triangle_twins(hierarchy=hierarchy)
        
    # Smooth mesh
    if smooth_steps > 0:
        logger.info('Smoothing Mesh Surfaces')
        mesh.smooth_surfaces(smooth_steps, step_size=0.3, max_gamma=10)

    logger.info(
        'Time to post-process mesh: ' +
        format_time(time.time()-start)
    )
    return mesh


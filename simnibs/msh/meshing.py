import os
import tempfile
import logging
import numpy as np
from simnibs.msh import mesh_io
import scipy.sparse
import scipy.ndimage

from ..utils.simnibs_logger import logger
from .._compiled import create_mesh


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
        ret = create_mesh.mesh_image(
                fn_image.encode(), fn_mesh.encode(),
                facet_angle, facet_size, facet_distance,
                cell_radius_edge_ratio, cell_size,
                optimize
             )
        if ret != 0:
            raise MeshingError('There was an error while meshing the image')
        mesh = mesh_io.read_medit(fn_mesh)
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
               cell_radius_edge_ratio=2, cell_size=None,
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
        See https://doc.cgal.org/latest/Mesh_3/index.html#title10. Default: 30

    facet_size: float (optional)
        See https://doc.cgal.org/latest/Mesh_3/index.html#title10. Default: minimum voxel
        size

    facet_distance: float (optional)
        See https://doc.cgal.org/latest/Mesh_3/index.html#title10. Default: minimum voxel
        size

    cell_radius_edge_ratio: float (optional)
        See https://doc.cgal.org/latest/Mesh_3/index.html#title10. Default: 2

    cell_size: float (optional)
        See https://doc.cgal.org/latest/Mesh_3/index.html#title10. Default: minimum voxel
        size

    Returns
    ----------
    msh: simnibs.Msh
        Mesh structure
    '''

    if image.dtype not in [np.uint8, np.uint16]:
        raise MeshingError('Image must be of type uint8 or uint16')

    rot, voxel_dims, shearing = _decompose_affine(affine)

    if not np.allclose(shearing, np.eye(3)):
        raise ValueError('Affine matrix has a shearing component')

    if facet_size is None:
        facet_size = min(voxel_dims)

    if facet_distance is None:
        facet_distance = min(voxel_dims)

    if cell_size is None:
        cell_size = min(voxel_dims)

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
        ret = create_mesh.mesh_surfaces(
                fn_surfaces, sd_formated, fn_mesh.encode(),
                facet_angle, facet_size, facet_distance,
                cell_radius_edge_ratio, cell_size,
                optimize
               )
        if ret != 0:
            raise MeshingError('There was an error while meshing the surfaces')
        mesh = mesh_io.read_medit(fn_mesh)
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


def smooth_surfaces(mesh, n_steps, step_size=.3):
    ''' Smoothes the mesh surfaces using Taubin smoothing
    IN-PLACE
    
    Can severely decrease tetrahedra quality

    Parameters
    ------------
    msh: simnibs.msh.Msh
        Mesh structure, must have inner triangles
    n_steps: int
        number of smoothing steps to perform
    step_size (optional): float 0<step_size<1
        Size of each smoothing step. Default: 0.3
    '''
    assert step_size > 0 and step_size < 1
    # Triangle neighbourhood information
    adj_tr = scipy.sparse.csr_matrix((mesh.nodes.nr, mesh.nodes.nr), dtype=bool)
    tr = mesh.elm[mesh.elm.triangles, :3] - 1
    ones = np.ones(len(tr), dtype=bool)
    for i in range(3):
        for j in range(3):
            adj_tr += scipy.sparse.csr_matrix(
                (ones, (tr[:, i], tr[:, j])),
                shape=adj_tr.shape
            )
    adj_tr -= scipy.sparse.dia_matrix((np.ones(adj_tr.shape[0]), 0), shape=adj_tr.shape)

    # Tetrahedron neighbourhood information
    th = mesh.elm[mesh.elm.tetrahedra] - 1
    adj_th = scipy.sparse.csr_matrix((mesh.nodes.nr, len(th)), dtype=bool)
    th_indices = np.arange(len(th))
    ones = np.ones(len(th), dtype=bool)
    for i in range(4):
        adj_th += scipy.sparse.csr_matrix(
            (ones, (th[:, i], th_indices)),
            shape=adj_th.shape
        )

    surf_nodes = np.unique(tr)
    nodes_coords = np.ascontiguousarray(mesh.nodes.node_coord, np.float)
    for i in range(n_steps):
        create_mesh.gauss_smooth(
            surf_nodes.astype(np.uint),
            nodes_coords,
            np.ascontiguousarray(th, np.uint),
            np.ascontiguousarray(adj_tr.indices, np.uint),
            np.ascontiguousarray(adj_tr.indptr, np.uint),
            np.ascontiguousarray(adj_th.indices, np.uint),
            np.ascontiguousarray(adj_th.indptr, np.uint),
            float(step_size)
        )
        # Taubin step
        create_mesh.gauss_smooth(
            surf_nodes.astype(np.uint),
            nodes_coords,
            np.ascontiguousarray(th, np.uint),
            np.ascontiguousarray(adj_tr.indices, np.uint),
            np.ascontiguousarray(adj_tr.indptr, np.uint),
            np.ascontiguousarray(adj_th.indices, np.uint),
            np.ascontiguousarray(adj_th.indptr, np.uint),
            -1.05 * float(step_size)
        )
    mesh.nodes.node_coord = nodes_coords


def relabel_spikes(m, label_a, label_b, target_label, adj_threshold=2, log_level=logging.DEBUG):
    ''' Relabels the spikes in a mesh volume, in-place

    A spike is defined as a tetrahedron in "label_a" or "label_b" which has at least one
    node in the other volume and at least "adj_threshold" faces adjacent to tetrahedra in
    "tarrget_label". For example, for relabeling GM and Skull spikes going through CSF,
    one can use
    relabel_spikes(m, 2, 4, 3)

    Parameters
    -----------
    m: simnibs.Msh
       Mesh structure
    label_a: int
        Volume label with spikes
    label_b: int
        Second volume label with spikes are located
    target_label: int
        Volume label where the spikes are locate
    adj_threshold: int (optional)
        Threshhold of number of adjacent faces for being considered a spike 
    '''
    logger.log(log_level, f'relabeling sipkes in {label_a} and {label_b} to {target_label}')
    def find_spikes():
        # First, get all tetrahedra which have a node that is shared between the tissues
        shared_nodes = m.find_shared_nodes([label_a, label_b])
        with_shared_nodes = (
            np.in1d(m.elm[:], shared_nodes)
            .reshape(-1, 4)
            .sum(axis=1, dtype=bool)
        )
        # Now, select the tetrahedra with at least adj_th faces adjacent to the tissue
        # with holes
        adj_labels = m.elm.tag1[adj_th - 1]
        adj_labels[adj_th == -1] = -1
        spikes = with_shared_nodes *\
            (np.sum(adj_labels == target_label, axis=1) >= adj_threshold)
        return spikes

    adj_th = m.elm.find_adjacent_tetrahedra()
    relabeled = True
    while np.any(relabeled):
        # Relabel tissue A
        A_to_relabel = (m.elm.tag1 == label_a) * find_spikes()
        frac_A_relabeled = np.sum(A_to_relabel)/np.sum(m.elm.tag1 == label_a)
        m.elm.tag1[A_to_relabel] = target_label
        m.elm.tag2[A_to_relabel] = target_label
        # Relabel tissue B
        B_to_relabel = (m.elm.tag1 == label_b) * find_spikes()
        frac_B_relabeled = np.sum(B_to_relabel)/np.sum(m.elm.tag1 == label_b)
        m.elm.tag1[B_to_relabel] = target_label
        m.elm.tag2[B_to_relabel] = target_label
        relabeled = A_to_relabel + B_to_relabel
        logger.log(log_level,
            f'Relabeled {frac_A_relabeled:.2%} of elements from {label_a} '
            f'and {frac_B_relabeled:.2%} of elements from {label_b}'
        )

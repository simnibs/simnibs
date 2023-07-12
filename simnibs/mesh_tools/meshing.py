import os
import tempfile
import logging
import numpy as np
import nibabel as nib
import scipy.sparse
import scipy.ndimage
import time

from simnibs.utils.mesh_element_properties import ElementTags

from . import mesh_io
from . import cgal
from ..utils import file_finder
from ..utils.spawn_process import spawn_process
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
               optimize=False):
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
    mesh.elm.tag1[mesh.elm.tag1 > ElementTags.TH_END] -= ElementTags.TH_SURFACE_START
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


def _get_connectivity_matrix(faces, nr_nodes = None):
    ''' get connectivity matrix of surface nodes from triangles

        Parameters
        ----------
        faces:
            n_facesx3 ndarray of triangle nodes
        nr_nodes:
            total number of nodes (set to np.max(faces)+1 if left empty)

        Returns
        -------
        boolean csc sparse matrix (nr_nodes x nr_nodes)
    '''
    if nr_nodes == None:
        nr_nodes = np.max(faces)+1
    nr_faces = faces.shape[0]

    i = np.tile( np.arange(nr_faces).reshape(-1,1), (1,3) )
    i = i.reshape(-1)
    j = faces.reshape(-1)
    ones = np.ones( i.shape[0], dtype=bool)

    C = scipy.sparse.csr_matrix( (ones, (i,j)), shape=(nr_faces, nr_nodes), dtype=bool )
    C=C.T*C

    nonzero, = C.diagonal().nonzero()
    C[nonzero, nonzero] = 0
    C.eliminate_zeros()
    return C


def _get_surfaces(faces, tet_faces, adj_tets, tag, nr_nodes):
    ''' reconstructs surfaces between tets of different labels

        Parameters
        ----------
        faces:
            n_facesx3 ndarray of triangle nodes
        tet_faces:
            n_tetsx4 ndarray of tet faces (indices into the faces array)
        adj_tets:
            n_tetx4 ndarray of tet neighbors (-1 in case of "air")
        tag:
            n_tetx1 ndarray of tet labels
        nr_nodes:
            total number of nodes

        Returns
        -------
        idx_surface_tri:
            ndarray of surface triangles (indices into the faces array)
        face_node_diff:
            nr_nodesx1 ndarray indicating difference between number of surface
            faces connected to a node and number of neigbhor surface nodes
        nneighb:
            nr_nodesx1 ndarray indicating number of neigbhor surface nodes
        conn_nodes:
            connectivity matrix of surface nodes (boolean csc sparse matrix,
                                                  nr_nodes x nr_nodes)
    '''
    adj_labels = tag[adj_tets]
    adj_labels[adj_tets == -1] = -1
    adj_diff = adj_labels - tag.reshape((len(tag),1)) != 0
    adj_diff[tag == -1] = False # no triangles for "air" tetrahedra

    # index of surface faces
    idx_surface_tri = np.unique(tet_faces[adj_diff])
    # node connectivity matrix
    conn_nodes = _get_connectivity_matrix(faces[idx_surface_tri], nr_nodes)
    # number of neighbor nodes
    nneighb=np.asarray( conn_nodes.sum(axis=1) ).reshape(-1)
    # number of surface faces connected to each node
    nfaces=np.bincount(faces[idx_surface_tri].reshape(-1),minlength=nr_nodes)
    # face_node_diff > 0 indicates a more complex surface topology (T-junction, ...)
    face_node_diff = nfaces-nneighb

    return idx_surface_tri, face_node_diff, nneighb, conn_nodes


def _get_elm_and_new_tag(tag, adj_tets, nr_diff, return_diffmat = False):
    ''' get index of all tets with nr_diff (2, 3 or 4) neigboring tets
        with different labels and returns also the labels of these neighbor tets

        Parameters
        ----------
        tag:
            n_tetx1 ndarray of tet labels
        adj_tets:
            n_tetx4 ndarray of tet neighbors (-1 in case of "air")
        nr_diff:
            required number of neighbors with different labels (2, 3 or 4)
        return_diffmat (optional; Default: False):
            whether to return a boolean matrix indicating neighbors with different labels

        Returns
        -------
        idx_elm:
            ndarray of selected tets
        new_tag:
            new tag based on the neighbor labels
        adj_diff (optional):
            n_tetx4 boolean ndarray indicating neighbor tets with different labels

        Notes
        -------
        * nr_diff == 4: the most frequent neighbor label is selected
        * nr_diff == 3: only tets are returned where at least 2 neighbors
                        have the same label; the most frequent neighbor label is
                        returned
        * nr_diff == 2: only tets are returned where the two neighbors
                        have the same label
    '''
    adj_labels = tag[adj_tets]
    adj_labels[adj_tets == -1] = -1
    adj_diff = adj_labels - tag.reshape((len(tag),1)) != 0

    idx_elm = np.where(np.sum(adj_diff, axis=1) == nr_diff)[0]
    adj_labels = adj_labels[idx_elm]
    adj_labels = adj_labels[adj_diff[idx_elm]].reshape((-1,nr_diff))

    if nr_diff == 4:
        new_tag = scipy.stats.mode(adj_labels, axis=1)[0].flatten()
    elif nr_diff == 3:
        new_tag, n_occ = scipy.stats.mode(adj_labels, axis=1)
        # ensure that at least 2 neighbors have the same label
        idx_relabel = (n_occ > 1).flatten()
        idx_elm = idx_elm[idx_relabel]
        new_tag = new_tag[idx_relabel].flatten()
    elif nr_diff == 2:
        # ensure that the two neighbors have the same label
        idx_relabel = np.diff(adj_labels).flatten() == 0
        idx_elm = idx_elm[idx_relabel]
        new_tag = adj_labels[idx_relabel,0]
    else:
        raise ValueError('nr_diff has to be 2,3 or 4')

    if return_diffmat:
        return idx_elm, new_tag, adj_diff
    return idx_elm, new_tag


def _get_test_nodes(faces, tet_faces, adj_tets, idx_surface_tri, face_node_diff,
                    nneighb, node_coord, tag, node_number_list, fast_track=False):
    ''' returns indices of the to-be-tested surface nodes that might potentially
        be a spike

        Parameters
        ----------
        faces:
            n_facesx3 ndarray of triangle nodes
        tet_faces:
            n_tetsx4 ndarray of tet faces (indices into the faces array)
        adj_tets:
            n_tetx4 ndarray of tet neighbors (-1 in case of "air")
        idx_surface_tri:
            ndarray of surface triangles (indices into the faces array)
        face_node_diff:
            nr_nodesx1 ndarray indicating difference between number of surface
            faces connected to a node and number of neigbhor surface nodes
        nneighb:
            nr_nodesx1 ndarray indicating number of neigbhor surface nodes
        node_coord:
            n_nodesx3 ndarray of node positions
        tag:
            n_tetx1 ndarray of tet labels
        node_number_list:
                n_tetx4 ndarray of node indices
        fast_track: bool (standard: False)
                when set to True, only nodes connected to three different
                regions are detected as putative spike nodes. This leads
                to far less nodes that need to be tested, making the next
                steps faster. However, spikes, e.g. in GM sulci will then
                not be detected.

        Returns
        -------
        idx_test_nodes :
            ndarray of selected candidate nodes

        Notes
        -------
        * uses heuristics to lower the number of candidate nodes
          while not loosing too many real spike nodes
        * required as _get_spikes_from_conn_matrix is rather slow
    '''
    idx_test_nodes = np.zeros(len(nneighb), dtype=bool)

    if not fast_track:
        # get nodes belonging to tets with 3 different neighbors
        idx_elm, _, adj_diff = _get_elm_and_new_tag(tag, adj_tets, 3, return_diffmat = True)
        # get the node that is shared by the three tet faces facing the different neighbor tets
        faces_elm = tet_faces[idx_elm]
        faces_elm = faces_elm[adj_diff[idx_elm]].reshape((-1,3))
        idx_node = scipy.stats.mode(faces[faces_elm].reshape((-1,9)), axis=1)[0]
        if len(idx_node) > 0:
            idx_test_nodes[idx_node] = True

        # get nodes belonging to tets with 2 different neighbors
        idx_elm = _get_elm_and_new_tag(tag, adj_tets, 2)[0]
        # get the 2 nodes that are shared by the two tet faces facing the different neighbor tets
        faces_elm = tet_faces[idx_elm]
        faces_elm = faces_elm[adj_diff[idx_elm]].reshape((-1,2))
        idx_node = np.sort(faces[faces_elm].reshape((-1,6)))
        idx_node = idx_node[:,:5][np.diff(idx_node) == 0]
        if len(idx_node) > 0:
            idx_test_nodes[idx_node] = True

        # the above steps are quite liberal and add too many nodes
        # therefore, exclude nodes with too low surface curvature
        # even though a few true spikes can get lost
        m_surf = mesh_io.Msh()
        m_surf.nodes.node_coord = node_coord
        m_surf.elm.add_triangles(faces[idx_surface_tri,:]+1,
                                  np.ones(len(idx_surface_tri),dtype=int))
        nd = m_surf.gaussian_curvature()
        nd.value = np.abs(nd.value)
        idx_test_nodes *= nd.value > 0.1

    # add nodes that belong to tets of three different regions
    C = scipy.sparse.lil_matrix((len(nneighb), np.max(tag)+2), dtype=bool)
    for i in range(4):
        C[node_number_list[:,i]-1, tag] = True
    C[:,-1] = False # do not count tags of air tets
    n_node_tags = np.asarray( C.sum(axis=1) ).reshape(-1)
    idx_test_nodes += n_node_tags == 3

    # a spike node needs at least 6 neighbor nodes
    idx_test_nodes *= nneighb > 5
    # exclude more complex surface geometries (e.g. T-junctions)
    idx_test_nodes *= face_node_diff == 0
    # exclude nodes connected to outer faces
    idx = tet_faces[adj_tets == -1]
    idx = np.unique(faces[idx])
    idx_test_nodes[idx] = False

    return idx_test_nodes


def _get_spikes_from_conn_matrix(conn_nodes, idx_test_nodes, nneighb):
    ''' test which nodes in idx_test_nodes are spike nodes

        Parameters
        ----------
        conn_nodes: boolean csc sparse matrix (nr_nodes x nr_nodes)
            node connectivity matrix
        idx_test_nodes: ndarray
            indices of candidate nodes taht should be tested
        nneighb: nr_nodesx1 ndarray
            number of neighbor surface nodes

        Returns
        -------
        ndarray of spike node indices

        Notes
        ------
        * nodes belonging to more complex surface configurations
         (i.e. face_node_diff > 0) must not be candidate nodes
    '''
    # restrict connectivty matrix to nodes on surfaces
    # and being test nodes or neighbors of test nodes
    # to gain some speed up

    def _get_spikes(conn_nodes, idx_test_nodes, nneighb):
        idx_surface_nodes = nneighb != 0
        idx_surface_nodes *= conn_nodes.dot(idx_test_nodes)
        idx_surface_nodes += idx_test_nodes
        conn_nodes = conn_nodes[:,idx_surface_nodes][idx_surface_nodes]

        map_nodes_new_old=np.where(idx_surface_nodes)[0]
        map_nodes_old_new=-1*np.ones(len(nneighb), dtype=int)
        map_nodes_old_new[idx_surface_nodes] = np.arange(np.sum(idx_surface_nodes))

        idx_test_nodes = map_nodes_old_new[idx_test_nodes]
        idx_spike_nodes = np.zeros_like(idx_test_nodes,dtype=bool)

        # loop over all test nodes and test whether their neighbors
        # are all connected to each other
        for (node, k) in zip(idx_test_nodes, range(len(idx_test_nodes))):
            idx = conn_nodes[:,node].nonzero()[0]
            c=conn_nodes[:,idx][idx]

            # try to resolve topologies where some nodes have
            # more than 2 neighbors
            nnz=c.getnnz(axis=0)
            while np.any(nnz > 2):
                idx_2 = nnz==2
                # 1) select nodes with more than 2 neighbors and which have
                # at least two neighboring nodes with excatly 2 neighbors
                idx_two2neighb = np.where( ( c.astype(int).dot(idx_2)>1 ) * ~idx_2 )[0]
                idx_not2 = np.where(~idx_2)[0]
                # 2) remove the connection(s) between these nodes and
                # other nodes with more than two neighboars
                if len(idx_two2neighb)>0:
                    ix, iy = np.meshgrid(idx_two2neighb, idx_not2)
                    ix = ix.ravel()
                    iy = iy.ravel()

                    idx_set = np.ravel(c[ix,iy])
                    c[ix[idx_set],iy[idx_set]]=False

                    idx_set = np.ravel(c[iy,ix])
                    c[iy[idx_set],ix[idx_set]]=False

                    c.eliminate_zeros()
                else:
                    break
                nnz=c.getnnz(axis=0)

            # test whether all nodes are connected to the first node
            a=c.getcol(0)
            nnodes = a.shape[0]
            for i in range(int(nnodes/2)-1):
                a += c.dot(a)
            idx_spike_nodes[k] = nnodes != a.nnz

        return map_nodes_new_old[idx_test_nodes[idx_spike_nodes]]

    # iterate over smaller smaller chunks
    # to speeds up slicing of sparse matrix
    idx_test_nodes = np.where(idx_test_nodes)[0]
    start_idx = np.arange(0, len(idx_test_nodes), 2500)
    stop_idx  = np.append(start_idx[1:]-1, len(idx_test_nodes)-1)
    idx_spike_nodes = np.array([],dtype=int)
    idx_hlp = np.zeros(len(nneighb),dtype = bool)

    for i, j in zip(start_idx, stop_idx):
        idx_hlp[:] = False
        idx_hlp[idx_test_nodes[i:j+1]] = True
        idx_spike_nodes = np.append(idx_spike_nodes,
                                    _get_spikes(conn_nodes, idx_hlp, nneighb))

    return idx_spike_nodes


def _get_new_tag_for_spikes(idx_spike_nodes, adj_tets, elm, tag, node_nr):
    ''' resolve spike by relabeling some of the tetrahedra connected
        to a spike node.

        The tetrahedra connected to a spike node are sorted into
        groups of connected tetrahedra (connected by sharing faces
        and having the same label). The smallest group is then relabelled
        to the label of the immediate neighbor group (i.e. sharing faces with
        the spike group)

        Parameters
        ----------
        idx_spike_nodes : ndarray
            indices of the spike nodes
        adj_tets:
            n_tetx4 ndarray of tet neighbors (-1 in case of "air")
        elm:
            n_tetx4 ndarray of node indices
        tag:
            n_tetx1 ndarray of tet labels
        nr_nodes:
            total number of nodes

        Returns
        -------
        tag:
            updated n_tetx1 ndarray of tet labels
        spike_data: list of tuples of the form
            (idx_spikenode (1-based), idx_spike_tets, tag_old, tag_new)

        Notes
        ------
        * elm must contain only tetrahedra
    '''
    # restrict to tets neighboring spike nodes
    # for speed up
    idx = np.zeros(node_nr,dtype=bool)
    idx[idx_spike_nodes] = True
    elm = np.copy(elm)-1

    map_new_old = np.where(np.any(idx[elm],axis=1))[0]
    map_old_new=-1*np.ones(len(elm), dtype=int)
    map_old_new[map_new_old] = np.arange(len(map_new_old))

    elm = elm[map_new_old]
    tag_spk = tag[map_new_old]
    adj_tets_spk = map_old_new[adj_tets[map_new_old]]

    tag_buff = -1*np.ones_like(tag_spk)
    idx_spike_tets = np.empty(0,dtype=int)
    spike_data = []
    for idx_node in idx_spike_nodes:
        idx_all_tets = np.where(np.any(elm == idx_node,axis=1))[0]
        idx_tets = np.copy(idx_all_tets)

        len_idx_spike_tets = 1000000 # initialze with some large number
        while len(idx_tets) > 0:
            idx_group = [idx_tets[0]]
            tag_group = tag_spk[idx_group[0]]
            len_idx_group = 0
            while len(idx_group) > len_idx_group:
                len_idx_group = len(idx_group)
                idx_group = np.append(idx_group, adj_tets_spk[idx_group].flatten())
                idx_group = idx_group[ (idx_group != -1)*(tag_spk[idx_group] == tag_group) ]
                idx_group = np.intersect1d(idx_group,idx_tets)
            if len_idx_group < len_idx_spike_tets:
                idx_spike_tets = idx_group
                len_idx_spike_tets = len_idx_group
            idx_tets = np.setdiff1d(idx_tets,idx_group)

        idx_neigh = np.intersect1d(adj_tets_spk[idx_spike_tets], idx_all_tets)
        idx_neigh = np.setdiff1d(idx_neigh, idx_spike_tets)
        tag_neigh = tag_spk[idx_neigh]
        if np.any(np.diff(tag_neigh)):
            logger.warning('ambiguous new tag for node ' + str(idx_node))
        tag_buff[idx_spike_tets] = tag_neigh[0]

        spike_data.append((idx_node+1, map_new_old[idx_spike_tets],
                           tag_spk[idx_spike_tets[0]], tag_neigh[0]))

    new_tag = -1*np.ones_like(tag)
    new_tag[map_new_old] = tag_buff
    tag = np.copy(tag)
    tag[new_tag != -1] = new_tag[new_tag != -1]
    return tag, spike_data


def _get_candidates_for_splitting(sp_dat, node_number_list, idx_surf_nodes):
    """
    Determines all spikes that might be suited for splitting. The basic idea is
    that a spike can be split if there is a second node to which all tets of the
    spike are also connected. Then, the spike can be split by adding a new node
    in the middle of the line between the orignal spike node and the second node.

    Parameters
    ----------
    sp_dat : list of tuples of the form
        (idx_spikenode (1-based), idx_spike_tets, tag_old, tag_new)
    node_number_list :
        n_tetx4 ndarray of node indices
    idx_surf_nodes : boolean ndarray
        True for nodes at region boundaries

    Returns
    -------
    splittest : list of tuples of the form
        (idx_spikenode (1-based), idx_2nd_node, tag_old, tag_new)
    sp2_spike_tets : nx2 ndarray
        tet indices of spikes having 2 tets
        (see _combine_small_spikes for further information)
    sp2_unique_nodes : nx1 ndarray
        indices of the spike nodes for the spikes with 2 tets

    """
    splittest = []
    sp2_spike_tets = np.empty((0,2),dtype='int32')
    sp2_unique_nodes = np.empty((0,2),dtype='int32')

    for sp in sp_dat:
        idx_sp_node = sp[0]
        idx_sp_tets = sp[1]
        old_tag = sp[2]
        new_tag = sp[3]

        if len(idx_sp_tets)<2:
            logger.info(' _get_candidates_for_splitting: spike with one tetrahedron should not occur at this stage')

        elif len(idx_sp_tets)==2:
            un, cn = np.unique(node_number_list[idx_sp_tets], return_counts=True)
            idx = np.where(cn == 1)[0]
            if len(idx) == 2:
                sp2_unique_nodes = np.append(sp2_unique_nodes, un[idx].reshape((1,2)), axis=0)
                sp2_spike_tets = np.append(sp2_spike_tets, idx_sp_tets.reshape((1,2)), axis=0)
            else:
                logger.info(' _get_candidates_for_splitting: weird spike with two tets')
        else:
            # try to find a 2nd node to which all tets of the spike are also connected
            un, cn = np.unique(node_number_list[idx_sp_tets], return_counts=True)
            idx = np.where((un != idx_sp_node) * (cn == len(idx_sp_tets)))[0]
            if len(idx) == 0:
                # if not all tets of the spike are connected to a 2nd node,
                # try to pick one non-surface node as 2nd node
                idx2 = np.where(~idx_surf_nodes[un-1])[0]
                if len(idx2) == 1:
                    splittest.append((idx_sp_node, un[idx2[0]], old_tag, new_tag))
                    # Note: all tets of the spike that are not connected to the
                    # 2nd node will not be tested and automatically get the standard new tag
            elif len(idx) == 1:
                splittest.append((idx_sp_node, un[idx[0]], old_tag, new_tag))
            else:
                logger.info(' _get_candidates_for_splitting: weird spike with two fully connected 2nd nodes')

    return splittest, sp2_spike_tets, sp2_unique_nodes


def _select_splits_from_candidates(splittest, node_number_list, node_coord, tag_org):
    """
    Determines the spikes that will be split from the candidates.
    The spikes have to fulfill two criteria:
        * All tets connected to both the spike node and the 2nd node have to
          have the same tag.
        * When projected on the line between the spike node and the 2nd node,
          all other nodes of these tets have to fall approx. in the middle of
          the line.

    Parameters
    ----------
    splittest : list of tuples of the form
        (idx_spikenode (1-based), idx_2nd_node, tag_old, tag_new)
    node_number_list :
        n_tetx4 ndarray of node indices
    node_coord :
        n_nodesx3 ndarray of node positions.
    tag_org :
        n_tetx1 ndarray of tet labels BEFORE any spike removal

    Returns
    -------
    splitlist : list of tuples of the form
        (idx_spikenode (1-based), idx_2nd_node, tag_old, tag_new)

    """
    splitlist = []
    for sp in splittest:
        idx_n1 = sp[0]
        idx_n2 = sp[1]
        idx_orgtets = np.where( np.any(node_number_list == idx_n1,axis=1) *
                                np.any(node_number_list == idx_n2,axis=1) )[0]
        if len(idx_orgtets) == 0:
            raise ValueError("The two nodes are not connected!")

        if np.max(tag_org[idx_orgtets]) != np.min(tag_org[idx_orgtets]):
            # this can happen when the 2nd node is part of the "ring" of surface nodes connected to the first node
            continue

        un = np.unique(node_number_list[idx_orgtets])
        idx_others = un[(un != idx_n1) * (un != idx_n2)]

        v1 = node_coord[idx_n1-1,:] - node_coord[idx_n2-1,:]
        v2 = node_coord[idx_n1-1,:] - node_coord[idx_others-1,:]
        norm_v1 = np.linalg.norm(v1)
        norm_v2 = np.linalg.norm(v2,axis=1)
        cosalpha = np.sum(v1*v2, axis=1)/(norm_v1*norm_v2)

        if np.all(np.abs(cosalpha-0.5) < 0.3): # 0.3 is a magic number determined during initial testing
            splitlist.append(sp)

    return splitlist


def _combine_small_spikes(sp2_unique_nodes, sp2_spike_tets, adj_tets,
                          node_number_list, tag_org, nr_nodes):
    """
    At thin interfaces, two spikes with each 2 tets can be directly next to each other.
    The function combines them to a common spike with 4 tets and returns them as list
    of spikes for splitting.

    Parameters
    ----------
    sp2_unique_nodes : nx2 ndarray
        tet indices of spikes having 2 tets
    sp2_spike_tets : nx1 ndarray
        indices of the spike nodes for the spikes with 2 tets
    adj_tets :
        n_tetx4 ndarray of tet neighbors (-1 in case of "air")
    node_number_list :
        n_tetx4 ndarray of node indices
    tag_org : TYPE
        DESCRIPTION.
    nr_nodes :
        total number of nodes

    Returns
    -------
    splitlist : list of tuples of the form
        (idx_spikenode (1-based), idx_2nd_node, tag_old, tag_new)

    """
    sp2_unique_nodes = np.sort(np.copy(sp2_unique_nodes), axis=1)
    un, cn = np.unique(sp2_unique_nodes, axis=0, return_counts = True)
    un = un[cn>1]

    splitlist = []
    for i in un:
        idx_spikes = np.where(np.all(i == sp2_unique_nodes,axis=1))[0]
        un2, cn2 = np.unique(node_number_list[sp2_spike_tets[idx_spikes]], return_counts=True)
        idx_sp_nodes = un2[cn2 == 4]

        if len(idx_sp_nodes) == 2:
            # get the correct tags and append to splitlist
            _, sp_dat_hlp = _get_new_tag_for_spikes(idx_sp_nodes-1, adj_tets, node_number_list,
                                                    tag_org, nr_nodes)
            splitlist.append((idx_sp_nodes[0], idx_sp_nodes[1], sp_dat_hlp[1][3], sp_dat_hlp[0][3]))

    return splitlist


def _split_spikes(m, splitlist):
    """
    Splits the spikes in the splitlist and updates the tags
    (works in place on the supplied mesh)

    Parameters
    ----------
    m :
        mesh of meshio.Msh() type
    splitlist : list of tuples of the form
        (idx_spikenode (1-based), idx_2nd_node, tag_old, tag_new)

    Returns
    -------
    idx_splittets :
        indices of the split tets (for visualization and debugging)

    """
    idx_splittets = np.empty((0),dtype='int32')
    for sp in splitlist:
        idx_n1 = sp[0]
        idx_n2 = sp[1]
        old_tag = sp[2]
        new_tag = sp[3]

        idx_tets1, idx_tets2 = m.split_tets_along_line(idx_n1,idx_n2,return_tetindices = True)
        m.elm.tag1[idx_tets1-1] = new_tag
        m.elm.tag1[idx_tets2-1] = old_tag

        idx_splittets=np.append(idx_splittets,idx_tets1)
        idx_splittets=np.append(idx_splittets,idx_tets2)

    m.elm.tag2[:] =  m.elm.tag1
    return idx_splittets


def update_tag_from_label_img(m, adj_tets, vol, affine, label_GM=None, label_CSF=None):
    ''' relables tags when:
            * tetrahedron more likely belongs to another label, based on the label image
            * tetrahdron has more than one face to a neighbor with different label

        Parameters
        ----------
        m:
            mesh of meshio.Msh() type
        adj_tets:
            n_tetx4 ndarray of tet neighbors (-1 in case of "air")
        vol : 3d ndarray
            label image
        affine : 4x4 ndarray
            affine transformation from voxel indices to world mm coordinates
        label_GM : int, optional
            label used for GM tets. The default is None.
        label_CSF : int, optional
            label used for CSF tets. The default is None.

        Returns
        -------
        m:
            mesh with updated m.elm.tag1 and m.elm.tag2

        Notes
        ------
        * The mesh must contain only tetrahedra
        * done repeatedly until no tets are relabeled anymore (max 20 times)
        * when labels for GM and CSF are given, CSF is not relabled to GM (avoids
          a few GM spikes)
        * updates tag1 and tag2 of m.elm
    '''

    # get most likely tag for each tet from label image
    best_tag = np.zeros_like(m.elm.tag1)
    best_tag_p = np.zeros_like(m.elm.tag1, dtype = np.float32)
    for i in np.unique(m.elm.tag1):
         ed = mesh_io.ElementData.from_data_grid(m, (vol[:]==i).astype(np.float32),
                                                 affine, '', order=1)
         idx = ed.value > best_tag_p
         best_tag[idx] = i
         best_tag_p[idx] = ed.value[idx]

    # relabel tets having more than one neighbor with different label
    def get_nr_diff_tag(tag, adj_tets):
        adj_labels = tag[adj_tets]
        adj_labels[adj_tets == -1] = -1
        adj_diff = adj_labels - m.elm.tag1.reshape((len(tag),1)) != 0
        return np.sum(adj_diff, axis=1)

    m.elm.tag2 = np.copy(m.elm.tag1)
    for i in range(20):
        nr_diff_pre = get_nr_diff_tag(m.elm.tag1, adj_tets)
        idx_relabel = (nr_diff_pre > 1)*(best_tag != m.elm.tag1)
        if label_GM is not None:
            idx_relabel[(m.elm.tag1 == label_CSF)*(best_tag == label_GM)] = False
        m.elm.tag1[idx_relabel] = best_tag[idx_relabel]

        # undo relabeling for tets that got more "spiky"
        nr_diff_post = get_nr_diff_tag(m.elm.tag1, adj_tets)
        diff_pre_post = nr_diff_pre - nr_diff_post
        idx_undo = idx_relabel*(diff_pre_post<0)
        m.elm.tag1[idx_undo] = m.elm.tag2[idx_undo]
        if (np.sum(idx_relabel) == np.sum(idx_undo)):
            break

    idx_relabel = m.elm.tag1 == 0 # set back tets relabled to 0 to their original label
    m.elm.tag1[idx_relabel] = m.elm.tag2[idx_relabel]
    logger.info('   Relabled ' + str(np.sum(m.elm.tag1 != m.elm.tag2)) + ' tets')
    m.elm.tag2[:] = m.elm.tag1
    return m


def update_tag_from_tet_neighbors(m, faces, tet_faces, adj_tets, nr_iter = 12):
    ''' relables tetrahedra when they are surrounded by
            * 4 neighbors all having a different label
            * at least 3 neighbors with a different label, whereby 2 of these
              3 neighbors need to share the same label (e.g. a tet with label 2
              surrounded by 2, 3, 4, 4 is relabeled to 4)
            * at least 2 neighbors with a different label, whereby these 2
              neighbors need to have the same label and a simple surface
              analysis indicates a surface defect

        this is done iteratively (standard: 12 times), whereby (most) tetrahedra
        that get repeatedly relabled are blocked from further relabeling;
        full convergence is not guaranteed

        Parameters
        ----------
        m:
            mesh of meshio.Msh() type
        faces:
            n_facesx3 ndarray of triangle nodes
        tet_faces:
            n_tetsx4 ndarray of tet faces (indices into the faces array)
        adj_tets:
            n_tetx4 ndarray of tet neighbors (-1 in case of "air")
        nr_iter: int, optional
            number of iteratins. The default is 10.

        Returns
        -------
        m:
            mesh with updated m.elm.tag1 and m.elm.tag2

        Notes
        ------
        * The mesh must contain only tetrahedra
        * updates tag1 and tag2 of m.elm
    '''
    tag = np.copy(m.elm.tag1)
    tag_buffer = np.copy(m.elm.tag1)
    relabeling_allowed = np.ones_like(m.elm.tag1, dtype = bool)
    just_relabelled = np.zeros_like(m.elm.tag1, dtype = bool)
    for i in range(nr_iter):
        just_relabelled[:] = False

        # relabel tets with 4 and 3 different neighbors
        for k in (4,3):
            idx_elm, new_tag = _get_elm_and_new_tag(tag, adj_tets, k)
            tag[idx_elm] = new_tag
            just_relabelled[idx_elm] = True

        # relabel tets with 2 different neighbors
        idx_elm, new_tag, adj_diff = _get_elm_and_new_tag(tag, adj_tets, 2, return_diffmat = True)
        # exclude tets at outer surface and ensure that tets can still be relabeled
        idx = new_tag > -1
        idx *= np.in1d( idx_elm, np.where(relabeling_allowed)[0] )
        idx_elm = idx_elm[idx]
        new_tag = new_tag[idx]
        # get the 2 nodes that are shared by the two tet faces facing the
        # different neighbor tets
        faces_elm = tet_faces[idx_elm]
        faces_elm = faces_elm[adj_diff[idx_elm]].reshape((-1,2))
        facenodes_elm = np.sort(faces[faces_elm].reshape((-1,6)))
        facenodes_elm = facenodes_elm[:,:5][np.diff(facenodes_elm) == 0].reshape((-1,2))
        # test whether these nodes are part of a surface defect
        # (note: get_elm_and_new_tag ensure only tets with diff neighbors that
        #  have the same label --> face_node_diff reveals defect, not T-junction)
        idx_surface_tri, face_node_diff = _get_surfaces(faces, tet_faces, adj_tets, tag, m.nodes.nr)[:2]
        idx = np.max(face_node_diff[facenodes_elm],axis=1)>0 # face_node_diff > 0 indicates a surface defect
        idx_elm = idx_elm[idx]
        new_tag = new_tag[idx]
        tag[idx_elm] = new_tag

        if i > 2:
            # stop flip-flopping between the orginal and a second label
            idx = (tag == m.elm.tag1) * (tag != tag_buffer)
            relabeling_allowed[idx] = False
            # prevent that relabeling tets with two faces converts them to back to tets with 3 or 4 faces
            idx = (tag == tag_buffer) * just_relabelled
            relabeling_allowed[idx] = False
        logger.info('     It. ' + str(i) + ': relabled ' + str(np.sum(tag_buffer != tag)) + ' tets')
        tag_buffer[:] = tag

    logger.info('   Relabeled ' + str(np.sum(m.elm.tag1 != tag)) + ' tets')
    m.elm.tag1 = tag
    m.elm.tag2[:] = m.elm.tag1
    return m


def update_tag_from_surface(m, faces, tet_faces, adj_tets, do_splits = False,
                            fast_track = False):
    ''' relables tetrahedra when they are part of a localized spike
        as detected by analysis of the surface topology

        Parameters
        ----------
        m:
            mesh of meshio.Msh() type
        faces:
            n_facesx3 ndarray of triangle nodes
        tet_faces:
            n_tetsx4 ndarray of tet faces (indices into the faces array)
        adj_tets:
            n_tetx4 ndarray of tet neighbors (-1 in case of "air")
        do_splits: bool
            whether to split spikes that reach deep into two regions
            (standard: False)
        fast_track: bool
            whether to test only nodes connected to three different regions.
            This leads to far less nodes that need to be tested. However,
            spikes, e.g. in GM sulci will then not be removed.
            (standard: False)

        Returns
        -------
        m:
            mesh with updated m.elm.tag1 and m.elm.tag2

        Notes
        ------
        * The mesh must contain only tetrahedra
        * updates tag1 and tag2 of m.elm
        * for do_splits=True: variables faces, tet_faces, adj_tets will
          be outdated after the call, as the function adds new tetrahedra
    '''
    DEBUG=False
    # reconstruct surface, get node-connectivity matrix
    idx_surface_tri, face_node_diff, nneighb, conn_nodes = _get_surfaces(faces, tet_faces, adj_tets,
                                                                         m.elm.tag1, m.nodes.nr)

    # get to-be-tested surface nodes (uses heuristics to lower the number of candidate nodes)
    idx_test_nodes = _get_test_nodes(faces, tet_faces, adj_tets, idx_surface_tri, face_node_diff,
                                     nneighb, m.nodes[:], m.elm.tag1, m.elm.node_number_list,
                                     fast_track = fast_track)

    # detect spikes by analysis of surface topology around each test node
    logger.info('     Testing '+ str(np.sum(idx_test_nodes)) + ' nodes (matrix: '
            + str(conn_nodes.shape) + ', ' + str(conn_nodes.nnz) + ' entries)')
    idx_spike_nodes = _get_spikes_from_conn_matrix(conn_nodes, idx_test_nodes, nneighb)

    # set new tag for spike nodes
    tag_buff = m.elm.tag1
    m.elm.tag1, sp_dat = _get_new_tag_for_spikes(idx_spike_nodes, adj_tets, m.elm[:],
                                                 m.elm.tag1, m.nodes.nr)
    logger.info('   Relabled ' + str(np.sum(m.elm.tag1 != m.elm.tag2)) + ' tets')
    m.elm.tag2[:] = m.elm.tag1

    if do_splits:
        # split spikes that reach deep into two regions
        #   step 1: determine candidate spikes that might be suited for splitting
        idx_surf_nodes = conn_nodes.getnnz(axis=0)>0
        splittest, sp2_tets, sp2_uniquenodes = _get_candidates_for_splitting(sp_dat, m.elm.node_number_list,
                                                                             idx_surf_nodes)
        #   step 2: determine the spikes that will be split from the candidates
        splitlist = _select_splits_from_candidates(splittest, m.elm.node_number_list,
                                                   m.nodes.node_coord, tag_buff)

        #   step 3: at thin interfaces, two spikes with each 2 tets can be directly next
        #   to each other --> combine to a common spike with 4 tets that will be split
        splitlist += _combine_small_spikes(sp2_uniquenodes, sp2_tets, adj_tets,
                                           m.elm.node_number_list, tag_buff, m.nodes.nr)

        #   step 4: split (updates the mesh in place)
        n_tet_pre = m.elm.nr
        idx_splittets = _split_spikes(m, splitlist)
        logger.info('   Split ' + str(m.elm.nr - n_tet_pre) + ' tets')
        m.elm.tag2 = m.elm.tag1.copy()

        if DEBUG:
            ed=np.zeros_like(m.elm.tag1)
            ed[idx_splittets-1] = 1
            ed=mesh_io.ElementData(ed)
            m.add_element_field(ed,'splittets')

    return m


def _remove_spikes(m, label_img, affine, label_GM = 2, label_CSF = 3):
    """
    wrapper function around the three spike removal steps

    Parameters
    ----------
    m : simnibs.Msh
        mesh from cgal, has to be without surfaces!
    label_img: 3D np.ndarray in uint8 format
        Labeled image from segmentation
    affine: 4x4 np.ndarray
        Affine transformation from voxel coordinates to world coordinates
    label_GM : int, optional
        label for GM volume. The default is 2.
    label_CSF : int, optional
        label for CSF volume. The default is 3.

    Returns
    -------
    msh: simnibs.Msh
        Mesh structure

    """
    logger.info('Removing Spikes')
    logger.info(' Step 1: Update tags from label image')
    faces, tet_faces, adj_tets = m.elm._get_tet_faces_and_adjacent_tets()
    tag_buff = m.elm.tag1.copy()
    m = update_tag_from_label_img(m, adj_tets, label_img, affine,
                                  label_GM = label_GM, label_CSF = label_CSF)

    logger.info(' Step 2: Update tags from tet neighbors')
    m = update_tag_from_tet_neighbors(m, faces, tet_faces, adj_tets)

    logger.info(' Step 3: Resolve remaining localized spikes ')
    m = update_tag_from_surface(m, faces, tet_faces, adj_tets, do_splits = True)

    logger.info('Done Removing Spikes: Total number of relabled tets: ' +
                str(np.sum(m.elm.tag1[:len(tag_buff)] != tag_buff)) +
                '; Number of split tets: ' + str(len(m.elm.tag1) - len(tag_buff)))

    # remove "air" tetrahedra with label -1 and corresponding nodes
    idx_keep = np.where(m.elm.tag1 != -1)[0] + 1
    m = m.crop_mesh(elements = idx_keep)

    return m


def _fix_labels(m, label_img):
    ''' Assign the right labels to the mesh as CGAL modifies them '''
    indices_seg, label_counts = np.unique(label_img, return_counts=True)
    indices_seg = indices_seg[1:]
    label_counts = label_counts[1:]
    indices_cgal = np.unique(m.elm.tag1)
    n_dropped = len(indices_seg)-len(indices_cgal)
    if n_dropped:
        idx_keep = np.argsort(label_counts)[::-1]
        idx_keep = idx_keep[:-n_dropped]
        indices_seg = np.sort(indices_seg[idx_keep])
        logger.warn('{} small region(s) dropped during meshing. Check label numbers in mesh!'.format(n_dropped))
    new_tags = np.copy(m.elm.tag1)
    for i, t in enumerate(indices_seg):
        new_tags[m.elm.tag1 == i+1] = t
    m.elm.tag1 = new_tags
    m.elm.tag2 = new_tags.copy()
    return m


def _relabel_microtets(m, el_max = 0.0001):
    """ CGAL can create spurious groups of microscopic tetrahedra
    at region boundaries. In order to get mmg to fix them, they are
    relabled to the most common tag in each group. By that, the
    boundary is moved away and mmg will resolve them even
    when -nosurf is set.

    Parameters
    ----------
    m : simnibs.Msh
        Mesh structure.
    el_max : float, optional
        maximal edge length. Tetrahedra will be relabeled
        when all edges are shorter than el_max.
        The default is 0.0001.

    Returns
    -------
    m : simnibs.Msh
        Mesh structure.

    """
    M = m.nodes[m.elm[:]]
    E = np.array([
        M[:, 0] - M[:, 1],
        M[:, 0] - M[:, 2],
        M[:, 0] - M[:, 3],
        M[:, 1] - M[:, 2],
        M[:, 1] - M[:, 3],
        M[:, 2] - M[:, 3]])
    E = np.swapaxes(E, 0, 1)
    # max Edge length
    Smax = np.max(np.linalg.norm(E, axis=2), axis=1)

    idx = np.argwhere(Smax<el_max).flatten()

    if len(idx) == 0:
        return m

    idx_c = m.elm.connected_components(idx+1)

    for i in idx_c:
        logger.debug(f"Relabeling group of {len(i)} micro-tets")
        m.elm.tag1[i-1] = np.argmax(np.bincount(m.elm.tag1[i-1]))
    m.elm.tag2[:] = m.elm.tag1
    return m


def _run_mmg(m, mmg_noinsert=True, fn_sol=None):
    """
    Wrapper around mmg command line call to improve mesh quality.

    Parameters
    ----------
    m : simnibs.Msh
        Mesh structure.
    mmg_noinsert : bool, optional
        set -noinsert flag to prevent mmg from adding nodes. The default is True.
    fn_sol : str, optional, default: None
        Filename of .sol file containing the sizing field (created with create_sizing_field_sol_file())

    Returns
    -------
    m : simnibs.Msh
        Mesh structure.
    """
    logger.info('Improving Mesh Quality')
    fn_tmp_in_gmsh = tempfile.NamedTemporaryFile().name + ".msh"
    fn_tmp_in_medit = tempfile.NamedTemporaryFile().name + ".mesh"
    fn_tmp_out = tempfile.NamedTemporaryFile().name + ".msh"
    mesh_io.write_msh(m, fn_tmp_in_gmsh)
    del m

    # set meshio convert command
    cmd = ["meshio", "convert", fn_tmp_in_gmsh, fn_tmp_in_medit,  # file_finder.path2bin("meshio")
           "--input-format", "gmsh", "--output-format", "medit"]

    # convert mesh from gmsh (.msh) to medit format (.mesh) for mmg
    spawn_process(cmd, lvl=logging.DEBUG)

    # set MMG command
    if mmg_noinsert:
        cmd = [file_finder.path2bin("mmg3d_O3"), "-v", "6", "-nosurf", "-hgrad", "-1", "-rmc", "-noinsert",
               "-in", fn_tmp_in_medit, "-out", fn_tmp_out]
    else:
        if fn_sol is not None:
            cmd = [file_finder.path2bin("mmg3d_O3"), "-v", "6", "-nosurf", "-hgrad", "-1", "-rmc",
                   "-hsiz", "100.0", "-hmin", "1.3", "-in", fn_tmp_in_medit, "-out", fn_tmp_out, "-sol", fn_sol]
        else:
            cmd = [file_finder.path2bin("mmg3d_O3"), "-v", "6", "-nosurf", "-hgrad", "-1", "-rmc",
                   "-hsiz", "100.0", "-hmin", "1.3", "-in", fn_tmp_in_medit, "-out", fn_tmp_out]

    # run MMG to improve mesh
    spawn_process(cmd, lvl=logging.DEBUG)

    # read mesh written by MMG (msh in ascii format)
    m = mesh_io.read_msh(fn_tmp_out, skip_data=True)

    # remove tmp-files
    try:
        os.remove(fn_tmp_in_gmsh)
    except:
        logger.warning(f'Could not delete {fn_tmp_in_gmsh}')

    try:
        os.remove(fn_tmp_in_medit)
    except:
        logger.warning(f'Could not delete {fn_tmp_in_medit}')

    try:
        os.remove(fn_tmp_out)
    except:
        logger.warning(f'Could not delete {fn_tmp_out}')

    return m


def create_mesh(label_img, affine,
                elem_sizes={"standard": {"range": [1, 5], "slope": 1.0}},
                smooth_size_field = 2,
                skin_facet_size=2.0,
                facet_distances={"standard": {"range": [0.1, 3], "slope": 0.5}},
                optimize=True, remove_spikes=True, skin_tag=1005,
                hierarchy=None, smooth_steps=5, skin_care=20,
                sizing_field=None, mmg_noinsert=False, debug=False):
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
    hierarchy: list of ints or None (optional)
        List of surface tags that determines the order in which triangles
        are kept for twin pairs (for remove_twins=True). Default for hierarchy=None:
        (1, 2, 9, 3, 4, 8, 7, 6, 10, 5)
        i.e. WM (1) has highest priority, GM (2) comes next, etc;
    smooth_steps: int (optional)
        Number of smoothing steps applied to the final mesh surfaces. Default: 5
    skin_care: int (optional)
        Number of addtional smoothing steps applied to the skin. Default: 20
    sizing_field: 3D np.ndarray in float format (optional)
        Sizing field to control the element sizes. Its shape has to be the same
        as label_img.shape. Zeros will be replaced by values from the
        standard sizing field. Default: None
    mmg_noinsert : bool, optional, default: False
        Set this flag to constrain the mesh improvement algorithm of MMG to not insert additional points.
        In this way, the number of elements of the mesh is not increased. (not recommended)

    Returns
    -------
    msh: simnibs.Msh
        Mesh structure
    """
    if hierarchy is None:
        hierarchy = (1, 2, 9, 3, 4, 8, 7, 6, 10, 5)
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

    # Replace values in size_field at positions where sizing_field > 0
    if sizing_field is not None:
        assert sizing_field.shape == label_img.shape
        size_field[sizing_field>0] = sizing_field[sizing_field>0]

    if debug:
        tmp_nii = nib.Nifti1Image(size_field, affine)
        nib.save(tmp_nii, 'size_field.nii.gz')
        tmp_nii = nib.Nifti1Image(distance_field, affine)
        nib.save(tmp_nii, 'distance_field.nii.gz')
        del tmp_nii

    # Run meshing
    logger.info('Meshing')
    logger.info('================================')
    logger.info('USING cell_radius_edge_ratio=2.1')
    logger.info('================================')
    start = time.time()
    m = image2mesh(
        label_img,
        affine,
        facet_size=size_field,
        facet_distance=distance_field,
        cell_size=size_field,
        optimize=optimize,
        cell_radius_edge_ratio=2.1
    )

    del size_field, distance_field
    logger.info(
        'Time to mesh: ' +
        format_time(time.time()-start)
    )

    # separate out tetrahedron (will reconstruct surfaces later)
    start = time.time()
    m = m.crop_mesh(elm_type=4)

    # assign the right labels to the mesh as CGAL modifies them
    m = _fix_labels(m, label_img)

    # relabel groups of microscopic tets to a common tag, so that mmg fixes them
    m = _relabel_microtets(m)

    if debug:
        mesh_io.write_msh(m, 'before_despike.msh')

    # remove spikes from mesh
    if remove_spikes:
        m = _remove_spikes(m, label_img, affine, label_GM=2, label_CSF=3)

    # keep only largest component
    idx = m.elm.connected_components()
    m = m.crop_mesh(elements=max(idx, key=np.size))

    # reconstruct surfaces
    logger.info('Reconstructing Surfaces')
    m.fix_th_node_ordering()
    m.reconstruct_unique_surface(hierarchy=hierarchy, add_outer_as=skin_tag)

    # smooth surfaces
    if smooth_steps > 0:
        logger.info('Smoothing Mesh Surfaces')
        m.smooth_surfaces(smooth_steps, step_size=0.3, max_gamma=10)

    if skin_care > 0:
        logger.info('Extra Skin Care')
        m.smooth_surfaces(skin_care, step_size=0.3, tags=skin_tag, max_gamma=10)

    if debug:
        mesh_io.write_msh(m, 'before_mmg.msh')

    # sizing fields can theoretically also applied by creating a NodeData field but this is currently not working.
    # add sizing field from sizing image to mesh in a NodeData field called "sizing_field:metric" for mmg
    # m.add_sizing_field(sizing_field=sizing_field, affine=affine)

    if sizing_field is not None:
        fn_sol = os.path.splitext(m.fn)[0] + ".sol"
        create_sizing_field_sol_file(mesh=m,
                                     sizing_field=sizing_field,
                                     affine=affine,
                                     fn_sol=os.path.splitext(m.fn)[0] + ".sol")
    else:
        fn_sol = None

    # improve mesh quality using mmg
    m = _run_mmg(m, mmg_noinsert, fn_sol)

    logger.info(
        'Time to post-process mesh: ' +
        format_time(time.time()-start)
    )
    return m


def create_sizing_field_sol_file(mesh, sizing_field, affine, fn_sol=None):
    """

    Parameters
    ----------
    mesh : Msh
        Head mesh the sizing field is computed for.
    sizing_field : str or nifti image or np.ndarray
        Filename of nifti image or nifti image of sizing field, or 3D numpy array with the same shape as
        label_img.shape, which will be applied to the nodes.
    affine : 4x4 ndarray
        Array describing the affine transformation from the data grid to the mesh space.
    fn_sol : str
        Filename of .sol file containing the sizing field in the nodes.
        If nothing is provided, the file will have the same name as the mesh with a .sol extension.
    """
    if fn_sol is None:
        os.path.splitext(mesh.fn)[0] + ".sol"

    # read sizing image
    if type(sizing_field) is str:
        sizing_image = nib.load(sizing_field)
        sizing_field = sizing_image.get_fdata()
    # if an image is passed read the sizing field data out of it
    elif type(sizing_field) is nib.nifti1.Nifti1Image:
        affine = sizing_field.affine
        sizing_field = sizing_field.get_fdata()

    if affine is None:
        raise ValueError("Please provide affine for sizing field.")

    # create a NodeData field containiong the sizing field with ":metric" tag for mmg
    sizing_field_NodeData = mesh_io.NodeData.from_data_grid(mesh=mesh,
                                                            data_grid=sizing_field,
                                                            affine=affine,
                                                            field_name='sizing_field:metric')

    # ensure positive element sizes
    sizing_field_NodeData.value = np.abs(sizing_field_NodeData.value)

    # write .sol file with sizing field
    with open(fn_sol, 'w') as f:
        f.write('MeshVersionFormatted 2\n')
        f.write('\n')
        f.write('Dimension 3\n')
        f.write('\n')
        f.write('SolAtVertices\n')
        f.write(f'{len(sizing_field_NodeData.value)}\n')
        f.write('1 1\n')
        np.savetxt(f, sizing_field_NodeData.value, newline="\n", fmt='%.6f')
        f.write('\n')
        f.write('End')

# def relabel_spikes(elm, tag, with_labels, adj_labels, label_a, label_b,
#                    target_label, labels, nodes_label, adj_th, adj_threshold=2,
#                    log_level=logging.DEBUG, relabel_tol=1e-6, max_iter=20,
#                    nlist=None,maxn=None):
#     ''' Relabels the spikes in a mesh volume, in-place

#     A spike is defined as a tetrahedron in "label_a" or "label_b"
#     which has at least one node in the other volume and
#     at least "adj_threshold" faces adjacent to tetrahedra in
#     "target_label".

#     Parameters
#     -----------
#     elm: ndarray
#        simnibs.Msh.elm[:] mesh structure
#     tag: ndarray
#         labels for elements (simnibs.Msh.elm.tag1)
#     with_labels: ndarray (ntag x nelements) bool
#         indicates if element contains each of the labels
#     adj_labels:  ndarray int
#         labels for adjacent elements
#     label_a: int
#         index for volume label with spikes
#     label_b: int
#         index for second volume label with spikes
#     target_label: int
#         index for volume label to relabel spikes to
#     labels: ndarray int
#         list of labels
#     nodes_label: ndarray (ntag x nelements) int
#         count of how many nodes in each element have each label, used when updating
#     adj_th: list
#        value of m.elm.find_adjacent_tetrahedra()
#     adj_threshold: int (optional)
#         Threshhold of number of adjacent faces for being considered a spike
#     relabel_tol: float (optional)
#         Fraction of the elements that indicates convergence
#     max_iter: int
#         Maximum number of relabeling iterations
#     nlist : ndarray int
#         list of which elements each nodes is connected to
#     maxn : ndarray int
#         number of elements that each nodes is connected to
#         (needed for lookup in nlist)
#     '''
#     logger.log(
#         log_level,
#         f'Relabeling spikes in {labels[label_a]} and {labels[label_b]} to {labels[target_label]}'
#     )
#     if not np.any(with_labels[label_a] * with_labels[label_b]):
#         return

#     for i in range(max_iter):
#         # Relabel tissue A
#         # Find spikes
#         A_to_relabel, frac_A_relabeled = _find_spikes(tag, label_a, label_b,
#                   with_labels, adj_labels, target_label, labels, adj_threshold)
#         # Update tags and adjlabels, with_labels and nodes_label in place
#         _update_tags(tag, elm, adj_th, with_labels, adj_labels, A_to_relabel,
#                      label_a, target_label, labels, nodes_label,nlist,maxn)

#         # Relabel tissue B
#         # Find spikes
#         B_to_relabel, frac_B_relabeled = _find_spikes(tag, label_b, label_a,
#                   with_labels, adj_labels, target_label, labels, adj_threshold)
#         # Update tags and adjlabels, with_labels and nodes_label in place
#         _update_tags(tag, elm, adj_th, with_labels, adj_labels, B_to_relabel,
#                      label_b, target_label, labels, nodes_label,nlist,maxn)

#         logger.log(log_level,
#                    f'Relabeled {np.sum(A_to_relabel)} from {labels[label_a]} '
#                    f'and {np.sum(B_to_relabel)} from {labels[label_b]}'
#                    )

#         # Stop if converge has been reached
#         if frac_A_relabeled < relabel_tol and frac_B_relabeled < relabel_tol:
#             break

#     # A_to_relabel, frac_A_relabeled = _find_spikes(tag, label_a, label_b,
#     #           with_labels, adj_labels, target_label, labels, adj_threshold)
#     # B_to_relabel, frac_B_relabeled = _find_spikes(tag, label_b, label_a,
#     #           with_labels, adj_labels, target_label, labels, adj_threshold)
#     # print(f'converged after {i+1} iterations, {np.sum(A_to_relabel)+np.sum(B_to_relabel)} left to relabel')


# @numba.njit(parallel=True, fastmath=True)
# def _find_spikes(tag, label, label2, with_labels, adj_labels, target_label, labels, adj_threshold=2):
#     '''
#     Find spikes

#     Parameters
#     ----------
#     tag : ndarray int
#         labels for elements
#     label : int
#             index for volume label with spikes
#     label2 : int
#             index for second volume label with spikes
#     with_labels : ndarray (ntag x nelements) bool
#             indicates if element contains each of the labels
#     adj_labels : ndarray int
#         labels for adjacent elements
#     target_label : int
#         target label index
#     labels : ndarray int
#             list of labels.
#     adj_threshold : list, optional
#         Threshold of number of adjacent faces for being considered a spike. The default is 2.

#     Returns
#     -------
#     ndarray bool
#         Indicates if spikes were found

#     '''
#     # Initialize output
#     found_spikes = np.zeros(tag.shape[0], dtype='bool')
#     # initialize number of elements with relevant label
#     na = 0
#     nspikes = 0
#     # Parallel loop over elements
#     for i in numba.prange(tag.shape[0]):
#     # if element has the label
#         if tag[i] == labels[label]:
#             # increment element count with label
#             na += 1
#             # if we have the other label too
#             if with_labels[label2, i]:
#                 # local variable (for thread) holding spikes count
#                 spikes = 0
#                 # loop over adjacent element
#                 for j in range(adj_labels.shape[1]):
#                     # if they have the target label too
#                     if adj_labels[i, j] == labels[target_label]:
#                         # increment spike count
#                         spikes += 1
#                         # if above threshold indicate spikes
#                         if spikes >= adj_threshold:
#                             found_spikes[i] = True
#                             # increment number of spikes found
#                             nspikes += 1
#                             # it is already a spike no need to go on
#                             break
#     # return found spikes and ratio of spikes (to check convergence)
#     return found_spikes, float(nspikes) / float(na)


# @numba.njit
# def _update_tags(tag, elm, adj_th, with_labels, adj_labels, to_relabel, label,
#                  target_label, labels, nodes_label,nlist,maxn):
#     '''
#     Update attributes needed for identifying spikes in place
#     Parameters
#     ----------
#     tag : ndarray int
#         labels for elements
#     elm : ndarray int
#         elements
#     adj_th : ndarray int
#         value of m.elm.find_adjacent_tetrahedra()
#     with_labels : ndarray bool

#     adj_labels : ndarray  int
#         labels for adjacent elements
#     to_relabel : ndarray bool
#         Elements to relabel.
#     label : int
#         original label index
#     target_label : int
#         new label
#     labels : ndarray int
#         list of labels
#     nodes_label : ndarray (ntag x nelements) int
#         count of how many nodes in each element have each label
#     nlist : ndarray int
#             list of which elements each nodes is connected to
#     maxn : ndarray int
#         number of elements that each nodes is connected to
#         (needed for lookup in nlist)
#     Returns
#     -------
#     None.

#     '''
#     # Loop over elements
#     for i in range(tag.shape[0]):
#         # relabel to intended label if indicated
#         if to_relabel[i]:
#             tag[i] = labels[target_label]
#             # update count of nodes with labels
#             for j in range(elm.shape[1]):
#                 nodes_label[label, elm[i, j]] -= 1  # decrease original label
#                 nodes_label[target_label, elm[i, j]] += 1  # increase intended label
#                 # loop over elements this node is connected to update with_labels
#                 # avoids looping over the entire mesh
#                 for m in range(maxn[elm[i,j]-1],maxn[elm[i,j]]):
#                     # Update only relevant if the element actually had the node
#                     if with_labels[label, nlist[m]]:
#                         with_labels[label, nlist[m]] = False
#                         for n in range(elm.shape[1]):
#                             if nodes_label[label, elm[nlist[m], n]] > 0:
#                                 with_labels[label, nlist[m]] = True
#                                 # no need to go on it is already True
#                                 break
#                     # Update can also be relevant if element does not have target label
#                     if not with_labels[target_label, nlist[m]]:
#                         for n in range(elm.shape[1]):
#                             # Update target label if nodes has it
#                             if nodes_label[target_label, elm[nlist[m], n]] > 0:
#                                 with_labels[target_label, nlist[m]] = True
#                                 # no need to go on it is already True
#                                 break
#             # loop over neighboors (4)
#             for k in range(adj_th.shape[1]):
#                 adj_i = adj_th[i, k] - 1  # indexing is from 1, 0 indicate no neighbor
#                 if adj_i >= 0:  # only if neighbor
#                     for j in range(adj_th.shape[1]):
#                         # if this element is the one being updated
#                         if adj_th[adj_i, j] - 1 == i:
#                             # update adjacent label
#                             adj_labels[adj_i, j] = labels[target_label]
#                             # no need to test more there is only one
#                             break


# @numba.njit
# def _with_label_numba_all(elm, tag1, labels, N):
#     ''' Returns a count of how many labels of each type belongs to each node
#         and all elements in the mesh which have a node that is in
#         a region with each label
#     '''
#     # output for counting label types for each node indexing starts from 1
#     # so first element is empty.
#     nodes_label = np.zeros((len(labels), N+1), dtype='uint16')
#     # output for boolean array indicating if element touches each node type
#     with_labels = np.zeros((len(labels), elm.shape[0]), dtype='bool')
#     maxn = np.zeros((N+1), dtype=elm.dtype)
#     # loop over labels
#     for k in range(len(labels)):
#         # Loop over elements
#         for i in range(elm.shape[0]):
#             # increment count if elements has current label
#             if tag1[i] == labels[k]:
#                 for j in range(elm.shape[1]):
#                     nodes_label[k, elm[i, j]] += 1
#                     maxn[elm[i,j]] += 1
#     # cummulative count of how many elements is connected nodes
#     # equivalent to np.cumsum(maxn) but in-place
#     for i in range(2,N+1):
#         maxn[i] += maxn[i-1]
#     # create an array containing the element number that each node is connected to
#     nlist = np.zeros((elm.shape[0]*elm.shape[1]),dtype=elm.dtype)
#     # list for counting how many elements a nodes has currently been asigned to
#     ncount = np.zeros((N+1,),dtype='uint16')
#     for i in range(elm.shape[0]):
#         for j in range(elm.shape[1]):
#             # element index (starting from 0 here)
#             n = elm[i,j]-1
#             # set the n'th element that the node is connected to
#             # ncount[n] counts how many has already been set
#             nlist[maxn[n] + ncount[n]] = i
#             #increment the elements that the node is connected to
#             ncount[n] += 1
#     # Loop over labels
#     for k in range(len(labels)):
#         # Loop over elements
#         for i in range(elm.shape[0]):
#         # Loop over nodes within label (4)
#             for j in range(elm.shape[1]):
#                 # Indicate if element contains label
#                 if nodes_label[k, elm[i, j]] > 0:
#                     with_labels[k, i] = True
#                     break
#     return nodes_label, with_labels, nlist, maxn


# def despike(msh, adj_threshold=2, relabel_tol=1e-6, max_iter=20,
#             log_level=logging.DEBUG):

#     ''' Goes through the mesh removing spiles
#     A spike is defined as a tetrahedron in a volume "a"
#     which has at least one node in the other volume "b" and
#     at least "adj_threshold" faces adjacent to tetrahedra in a volume "c"

#     Parameters
#     -----------
#     m: simnibs.Msh
#        Mesh structure
#     adj_threshold: int (optional)
#         Threshhold of number of adjacent faces for being considered a spike
#     relabel_tol: float (optional)
#         Fraction of the elements that indicates convergence
#     max_iter: int
#         Maximum number of relabeling iterations
#     '''

#     if np.any(msh.elm.elm_type != 4):
#         logger.log(log_level,
#                    'Error: Attempting to despike mesh containing not only'
#                    'tetrahedra. Please consider cropping the mesh first.'
#                    )
#         raise ValueError()
#         return

#     tags = np.unique(msh.elm.tag1)
#     adj_th = msh.elm.find_adjacent_tetrahedra()
#     elm = msh.elm[:]
#     tag = msh.elm.tag1
#     adj_labels = tag[adj_th - 1]
#     adj_labels[adj_th == -1] = -1

#     # Total number of nodes
#     N = msh.nodes.nr
#     # Count how many labels of each type belongs to each node - first output
#     # and determine if elements touch each labels (has a node with that label) - second output
#     # and determine which element each node is connected too - third output
#     nodes_label, with_labels, nlist, maxn = _with_label_numba_all(elm, tag, tags, N)
#     # Loop over labels
#     for i, t1 in enumerate(tags):
#         for j, t2 in enumerate(tags[i + 1:]):
#             # Only if the labels are not the same
#             if t1 == t2:
#                 continue
#             #only if at least one elements have this label combination
#             if not np.any((nodes_label[i] > 0) * (nodes_label[j + i + 1] > 0)):
#                 continue
#             for k, t3 in enumerate(tags):
#                 # Only if target label is different
#                 if t1 == t3 or t2 == t3:
#                     continue
#                 #only if at least one elements have this label combination
#                 if not np.any((nodes_label[i] > 0) *
#                               (nodes_label[j + i + 1] > 0) * nodes_label[k]):
#                     continue
#                 #call relabel function
#                 relabel_spikes(elm = elm, tag = tag, with_labels = with_labels,
#                                adj_labels = adj_labels, label_a = i,
#                                label_b = j + i + 1, target_label = k,
#                                relabel_tol = relabel_tol, labels = tags,
#                                adj_threshold = adj_threshold, adj_th = adj_th,
#                                max_iter = max_iter, log_level = log_level,
#                                nodes_label = nodes_label,nlist=nlist,maxn=maxn)
#     #set tag1/tag2 in msh structure
#     msh.elm.tag1 = tag
#     msh.elm.tag2 = tag


# def create_mesh(label_img, affine,
#                 elem_sizes={"standard": {"range": [1, 5], "slope": 1.0}},
#                 smooth_size_field = 2,
#                 skin_facet_size=2.0,
#                 facet_distances={"standard": {"range": [0.1, 3], "slope": 0.5}},
#                 optimize=True, remove_spikes=True, skin_tag=1005,
#                 remove_twins=True, hierarchy=None, smooth_steps=5, sizing_field=None):
#     """Create a mesh from a labeled image.

#     The maximum element sizes (CGAL facet_size and cell_size) are controlled
#     by elem_sizes:
#         size = slope * thickness
#         size[size < range[0] = range[0]
#         size[size > range[1] = range[1]
#     where "thickness" is the local tissue thickness,
#     "range" is the size range (label-specific if label is added to elem_sizes,
#                                 otherwise the "standard" range is used)
#     The distance (CGAL facet_distance) parameter is calcualted in a similar way.
#     This allows for the meshing to adjust sizes according to local needs.

#     Parameters
#     ----------
#     label_img: 3D np.ndarray in uint8 format
#         Labeled image from segmentation
#     affine: 4x4 np.ndarray
#         Affine transformation from voxel coordinates to world coordinates
#     elem_sizes: dictionary (optional)
#         Lists the relationship between thickness and elem_sizes.
#         Label-specific relationships can be added if needed, e.g. for label 2:
#             {"standard": {"range": [1, 5], "slope": 1.0},
#                     "2": {"range": [1, 2], "slope": 0.7}}
#         "range" determines the minimum and maximum values for element sizes.
#         "slope" determines relationship between thickness and element sizes.
#         Note: Label indices are used as keys, and need to be provided as string
#         Default: {"standard": {"range": [1, 5], "slope": 1.0}}
#     smooth_size_field: int (optional)
#         Defines the size of a triangular kernel to smooth the size field. A bit
#         of smoothing helps to remove the effect of a few outliers in the
#         thickness estimates on the size field. Set to 0 to disable.
#         The kernel size is 2*smooth_size_field+1 in voxels. Default:  2
#     skin_facet_size: float (optional)
#         Maximum size for the triangles of the outermost surface. If set to None,
#         the elements of the outermost surface will be scaled the same way as
#         the other surfaces. Default:  2.0
#     facet_distances: dictionary (optional)
#         Relationship between thickness and facet_distance. For small
#         facet_distance values, the meshing will follow the label boundaries
#         in the original image more strictly. This also means more elements.
#         Label-specific relationships can be added if needed.
#         "range": Minimum and maximum values for facet_distance
#         "slope": Steepness of relationship between thickness and facet_distance.
#         Default: {"standard": {"range": [0.1, 3], "slope": 0.5}}
#     optimize: bool (optional)
#         Whether to run lloyd optimization on the mesh. Default: True
#     remove_spikes: bool (optional)
#         Whether to remove spikes to create smoother meshes. Default: True
#     skin_tag: float (optional)
#         1) Restrict effects of skin_facet_size to the outer boundary of the
#            region with label skin_tag-1000 (i.e. 5 for 1005)
#         2) Add outer surface to mesh using given tag. Set to None to disable.
#         NOTE: This surface will replace any other surface with the same tag.
#         Default: 1005
#     remove_twins: bool (optional)
#         Remove triangle twins created during surface reconstruction.
#         Default: True
#     hierarchy: list of ints or None (optional)
#         List of surface tags that determines the order in which triangles
#         are kept for twin pairs (for remove_twins=True). Default for hierarchy=None:
#         (1005, 1001, 1002, 1009, 1003, 1004, 1008, 1007, 1006, 1010)
#         i.e. Skin (1005) has highest priority, WM (1001) comes next, etc;
#     smooth_steps: int (optional)
#         Number of smoothing steps to apply to the final mesh surfaces. Default: 5
#     sizing_field: 3D np.ndarray in float format (optional)
#         Sizing field to control the element sizes. Its shape has to be the same
#         as label_img.shape. Zeros will be replaced by values from the
#         standard sizing field. Default: None

#     Returns
#     -------
#     msh: simnibs.Msh
#         Mesh structure
#     """
#     if hierarchy is None:
#         hierarchy = (1005, 1001, 1002, 1009, 1003, 1004, 1008, 1007, 1006, 1010)
#     if not 'standard' in elem_sizes:
#         raise ValueError('elem_sizes needs a \"standard\" entry')
#     if not 'standard' in facet_distances:
#         raise ValueError('facet_distances needs a \"standard\" entry')

#     # Calculate thickness
#     logger.info('Calculating tissue thickness')
#     start = time.time()
#     thickness = _calc_thickness(label_img)
#     thickness[thickness < .5] = 100 # set background thickness to some large value
#     voxel_size = get_vox_size(affine)
#     if not np.allclose(np.diff(voxel_size), 0):
#         logger.warn('Anisotropic image, meshing may contain extra artifacts')
#     thickness *= np.average(voxel_size) # Scale thickness with voxel size

#     # Define size fields and distance field
#     logger.info('Calculating sizing fields')
#     size_field = _sizing_field_from_thickness(
#         label_img, thickness, elem_sizes
#     )
#     distance_field = _sizing_field_from_thickness(
#         label_img, thickness, facet_distances
#     )
#     del thickness

#     # Smooth size field a bit to reduce effect of a few outliers in the thickness
#     # map on the mesh; the outliers show up as localized small thickness values at
#     # some of the tissue boundaries
#     if smooth_size_field:
#         size_field = size_field**(1/3) # compress high values to preserve edges a bit better
#         kernel = smooth_size_field+1-np.abs(np.arange(-smooth_size_field,
#                                                       smooth_size_field+1, 1))
#         kernel = kernel/np.sum(kernel)
#         for i in range(3):
#             size_field = scipy.ndimage.convolve1d(size_field, kernel, axis=i,
#                                             mode='constant', cval=0.0, origin=0)
#         size_field = size_field**3

#     # Control triangle size of outer surface to ensure eletrode meshing works OK
#     if skin_facet_size is not None:
#         boundary = (label_img > 0).astype('int8')
#         boundary = boundary-erosion(boundary,1)
#         if skin_tag is not None:
#             # keep boundary only at regions with label skin_tag-1000
#             # to save some tetrahedra
#             skin_tet_tag = skin_tag-1000
#             boundary *= (label_img == skin_tet_tag)
#             boundary = dilate(boundary,1)
#             boundary *= (label_img == skin_tet_tag)
#         else:
#             boundary = dilate(boundary,1)
#         size_field = size_field.flatten()
#         size_field[boundary.flatten()] = skin_facet_size
#         size_field = size_field.reshape(label_img.shape)
#         del boundary
#     logger.info(
#         'Time to prepare meshing: ' +
#         format_time(time.time()-start)
#     )

#     # Replace values in size_field at positions where sizing_field > 0
#     if sizing_field is not None:
#         assert sizing_field.shape == label_img.shape
#         size_field[sizing_field>0] = sizing_field[sizing_field>0]

#     # Run meshing
#     logger.info('Meshing')
#     start = time.time()
#     mesh = image2mesh(
#         label_img,
#         affine,
#         facet_size=size_field,
#         facet_distance=distance_field,
#         cell_size=size_field,
#         optimize=optimize
#     )
#     del size_field, distance_field
#     logger.info(
#         'Time to mesh: ' +
#         format_time(time.time()-start)
#     )

#     # Separate out tetrahedron (will reconstruct triangles later)
#     start = time.time()
#     mesh = mesh.crop_mesh(elm_type=4)
#     # Assign the right labels to the mesh as CGAL modifies them
#     indices_seg = np.unique(label_img)[1:]
#     new_tags = np.copy(mesh.elm.tag1)
#     for i, t in enumerate(indices_seg):
#         new_tags[mesh.elm.tag1 == i+1] = t
#     mesh.elm.tag1 = new_tags
#     mesh.elm.tag2 = new_tags.copy()

#     # Remove spikes from mesh
#     if remove_spikes:
#         logger.info('Removing Spikes')
#         despike(
#             mesh, relabel_tol=1e-5,
#             adj_threshold=2
#         )

#     # Reconctruct mesh surfaces
#     logger.info('Reconstructing Surfaces')
#     mesh.fix_th_node_ordering()
#     idx=mesh.elm.connected_components()
#     mesh = mesh.crop_mesh(elements=max(idx,key=np.size))
#     mesh.reconstruct_surfaces(add_outer_as=skin_tag)
#     if remove_twins:
#         mesh = mesh.remove_triangle_twins(hierarchy=hierarchy)

#     # Smooth mesh
#     if smooth_steps > 0:
#         logger.info('Smoothing Mesh Surfaces')
#         mesh.smooth_surfaces(smooth_steps, step_size=0.3, max_gamma=10)

#     logger.info(
#         'Time to post-process mesh: ' +
#         format_time(time.time()-start)
#     )
#     return mesh

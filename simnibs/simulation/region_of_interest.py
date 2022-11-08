import numpy as np
import scipy.sparse as sparse

from ..utils.simnibs_logger import logger

# TODO: properly implement the following postproc methods to compute e-fields from potential very fast

def _get_gradient(local_dist):
    ''' Calculate the gradient of a function in each tetrahedra
    The way it works: The operator has 2 parts
    local_dist * gradient = project_matrix
    project_matrix is a projection matrix
    project_matrix = [-1, 1, 0, 0]
        [-1, 0, 1, 0]
        [-1, 0, 0, 1]
    And local_dist is the distances between local node [0] and nodes [1], [2] and [3] in each element
    '''

    # define projection matrix
    project_matrix = np.hstack([-np.ones((3, 1)), np.eye(3)])

    # solve local_dist * gradient = project_matrix.
    gradient = np.linalg.solve(local_dist, project_matrix[None, :, :])

    # swapaxes
    return np.swapaxes(gradient, 1, 2)


def _get_edge_lengths(mesh):
    ''' get length of the tetrahedral edges '''

    node_coordinates = mesh.nodes.node_coord.T
    # get the local node coordinates in each element
    # node_coordinates: (3, number_of_nodes), node_numbers: (number_of_elements, 4), local_node_coords: (3, number_of_elements, 4)
    local_node_coords = node_coordinates[:, mesh.elm.node_number_list - 1]

    # get the distances between local node [0] and nodes [1], [2] and [3] in each element
    # local_dist: (3, number_of_elements, 3)
    local_dist = local_node_coords[:, :, 1:] - local_node_coords[:, :, [0]]

    # out: (number_of_elements, 3, 3)
    return np.moveaxis(local_dist, 0, -1), (local_node_coords[:, :, 0].T)[:, :, None]


class RegionOfInterest:
    """ Region of interest class containing methods to compute the electric field from the electric
    potential very fast.

    Parameters
    ----------


    Attributes
    ----------

    """
    def __init__(self, points, msh, con=None, gradient=None, out_fill=0):
        """
        Initializes RegionOfInterest class instance
        """
        self.points = points
        assert self.points.shape[1] == 3
        self.n_points = self.points.shape[0]
        self.con = con

        # crop mesh that only tetrahedra are included
        msh_cropped = msh.crop_mesh(elm_type=4)

        # ensure that the nodes did not change
        assert msh_cropped.nodes.nr == msh.nodes.nr
        assert (msh_cropped.nodes.node_coord == msh.nodes.node_coord).all()

        # compute sF matrix
        self._get_sF_matrix(msh_cropped, self.points, out_fill)

        if not gradient:
            # get the lengths of the tetrahedral edges
            local_dist, _ = _get_edge_lengths(msh_cropped)
            # get gradient operator
            gradient = _get_gradient(local_dist)
        self.gradient = gradient[self.idx]
        self.node_index_list = msh_cropped.elm.node_number_list[self.idx] - 1

    def _get_sF_matrix(self, msh, points, out_fill, tags=None):
        """Create a sparse matrix for SPR interpolation from element data to arbitrary positions (here: the surface nodes)

        Args:
            msh (Msh object): loaded mesh
            points (numpy array): the coordinates of the points that we want to interpolate (number_of_points, 3)
            out_fill (number or None): Value to be given to points outside the volume. If None then use nearest neighbor assigns the nearest value; otherwise assign to out_fill, for example 0)
            tags (list or None, optional): The tissue type defines the selected volume in the loaded mesh. Defaults to None.

        Sets the following object variables:
           sF       sparse.csr_matrix: (number_of_points, number_of_kept_tetrahedra)
           idx      index of the tetrahedra included in sF as columns
           inside   indices of the positions inside the mesh
        """

        assert out_fill in [0, 1]

        # Set the volume to be GM. The interpolation will use only the tetrahedra in the volume.
        if tags is None:
            th_indices = msh.elm.elm_number
        else:
            th_indices = msh.elm.elm_number[np.in1d(msh.elm.tag1, tags)]

        th_with_points, bar = msh.find_tetrahedron_with_points(points, compute_baricentric=True)
        inside = np.isin(th_with_points, th_indices)
        self.inside = inside

        # if any points are inside
        if np.any(inside):
            # get the 'elm_number' of the tetrahedra in 'msh' which contain 'points' in 'points' order
            # assert where_inside.shape == th.shape
            th = th_with_points[inside]

            # interpolate the E field to the points using SPR
            sF = self._get_sF_inside_tissues(msh, th, np.where(inside)[0], bar, points.shape[0])

        # Finally, fill in the unassigned values
        if np.any(~inside):
            if out_fill == 1:  # fill == 'nearest'
                if tags is not None:
                    is_in = np.in1d(msh.elm.elm_number, th_indices)
                    elm_in_volume = msh.elm.elm_number[is_in]
                    m_in_volume = msh.crop_mesh(elements=elm_in_volume)

                    _, nearest = m_in_volume.find_closest_element(points[~inside], return_index=True)

                    sF[np.where(~inside)[0], elm_in_volume[nearest - 1] - 1] = 1
                else:
                    _, nearest = msh.find_closest_element(points[~inside], return_index=True)

                    sF[np.where(~inside)[0], nearest - 1] = 1

        # convert to csc matrix for fast column indexing
        sF = sF.tocsc()
        self.idx = np.nonzero(sF.indptr[:-1] != sF.indptr[1:])[0]

        # delete the all-zero columns in sF
        sF = sF[:, self.idx]

        # convert to csr matrix for fast row indexing in the matrix multiplication
        self.sF = sparse.csr_matrix(sF)

    def _get_sF_inside_tissues(self, msh, th, w, bar, n_points):
        """  Create a sparse matrix to interpolate from element data to
        arbitrary positions (here: the surface nodes) using superconvergent
        patch recovery

        Args:
            msh (Msh object): loaded mesh
            th (numpy array): indices of the elements (start from 0) that contains the points
            w (numpy array): indices of the points that are inside the mesh
            bar (numpy array): barycentric coordinates
            n_points (int): number of points

        Returns:
            sF (sparse matrix): the output sparse matrix

        """

        # initialize the sparse matrix
        sF = sparse.dok_matrix((n_points, msh.elm.nr))

        # get the 'tag1' from 'msh' for every element in 'th' in 'points' order
        tag1_inside = msh.elm.tag1[th - 1]

        for t in np.unique(tag1_inside):
            # find the elements in 'tag1_inside' which equals to 't'
            is_t = tag1_inside == t

            if np.any(is_t):
                logger.debug('points inside volumn with tag1 = {}: {}'.format(t, np.sum(is_t)))

                # 'm_tag' contains only the tetrahedra with 'elm_number == th_with_t'
                m_tag = msh.crop_mesh(tags=t)

                # 'm_with_t' is sorted because 'elm_number' is always sorted
                m_with_t = msh.elm.elm_number[msh.elm.tag1 == t]

                # the 'elm_number' of elements in 'msh'. These elements contain points and 'tag1 == t'
                th_with_t = th[is_t]

                # get the indices of elements in 'm_tag'. The 'elm_number' of the same elements are 'th_with_t' in 'msh'. 'idx' starts from 0, not 1.
                idx = np.searchsorted(m_with_t, th_with_t)

                # convert the 'elmdata' to 'nodedata'
                sM = self._elm2point_SPR(m_tag, bar[w[is_t]], idx)

                i, j, v = sparse.find(sM)
                sF[(w[is_t])[i], m_with_t[j] - 1] = v

        return sF

    def _elm2point_SPR(self, msh, bar, idx):
        """Create the sparse matrix to interpolate from element data to arbitrary positions
        (here: the surface nodes) using superconvergent patch recovery

        Args:
            msh (Msh object): loaded mesh
            bar (numpy array): barycentric coordinates
            idx (numpy array): indices of the elements (start from 0) that need to calculate SPR

        Raises:
            ValueError: msh has to contain volume data

        Returns:
            sparse.csr_matrix: (bar.shape[0], msh.elm.nr)

        References
            Zienkiewicz, Olgierd Cecil, and Jian
            Zhong Zhu. "The superconvergent patch recovery and a posteriori error
            estimates. Part 1: The recovery technique." International Journal for
            Numerical Methods in Engineering 33.7 (1992): 1331-1364.
        """

        if len(msh.elm.tetrahedra) == 0:
            raise ValueError("Can only transform volume data")

        # Get the point in the outside surface
        points_outside = np.unique(msh.elm.get_outside_faces())
        outside_points_mask = np.in1d(msh.elm[msh.elm.tetrahedra],
                                      points_outside).reshape(-1, 4)

        th_indices = msh.elm.tetrahedra
        th_nodes = msh.elm[th_indices] - 1
        masked_th_nodes = np.copy(th_nodes)
        masked_th_nodes[outside_points_mask] = -1
        masked_th_nodes += 1

        # Calculates the quantities needed for the superconvergent patch recovery
        uq_in, th_nodes = np.unique(masked_th_nodes, return_inverse=True)
        baricenters = msh.elements_baricenters()[th_indices]

        volumes = msh.elements_volumes_and_areas()[th_indices]
        baricenters = np.hstack([np.ones((baricenters.shape[0], 1)), baricenters])

        A = np.empty((msh.nodes.nr + 1, 4, 4))
        for i in range(4):
            for j in range(i, 4):
                A[:, i, j] = np.bincount(masked_th_nodes.reshape(-1),
                                         np.repeat(baricenters[:, i], 4) *
                                         np.repeat(baricenters[:, j], 4),
                                         minlength=msh.nodes.nr + 1)

        # This here only ensures we can invert
        outside = np.isclose(A[:, 0, 0], 0)
        for i in range(4):
            A[outside, i, i] = 1

        A[:, 1, 0] = A[:, 0, 1]
        A[:, 2, 0] = A[:, 0, 2]
        A[:, 3, 0] = A[:, 0, 3]
        A[:, 2, 1] = A[:, 1, 2]
        A[:, 3, 1] = A[:, 1, 3]
        A[:, 3, 2] = A[:, 2, 3]

        Ainv = np.linalg.inv(A)

        node_pos = np.hstack(
            [np.ones((msh.nodes.nr, 1)), msh.nodes.node_coord])
        # Added a dummy to the first position
        node_pos = np.vstack([np.ones((1, 4)), node_pos])

        M = sparse.csr_matrix((msh.nodes.nr + 1, msh.elm.nr))
        for i in range(4):
            M += sparse.csr_matrix(
                (np.einsum(
                    'bi, bij, bj -> b',
                    node_pos[masked_th_nodes[:, i]],
                    Ainv[masked_th_nodes[:, i]],
                    baricenters),
                 (masked_th_nodes[:, i], th_indices - 1)),
                shape=M.shape)

        # Assigns the average value to the points in the outside surface
        th_nodes = msh.elm[th_indices] - 1
        masked_th_nodes = np.copy(th_nodes)
        masked_th_nodes[~outside_points_mask] = -1
        masked_th_nodes += 1

        for i in range(4):
            M += sparse.csr_matrix(
                (volumes, (masked_th_nodes[:, i], th_indices - 1)),
                shape=M.shape)

        M = M[1:]

        node_vols = np.bincount(
            th_nodes.reshape(-1),
            np.repeat(volumes, 4),
            minlength=msh.nodes.nr + 1)

        normalization = np.ones(msh.nodes.nr)
        normalization[points_outside - 1] = 1 / (node_vols[points_outside - 1] + np.finfo(float).eps)

        D = sparse.dia_matrix(
            (normalization, 0), shape=(msh.nodes.nr, msh.nodes.nr))
        M = D.dot(M)

        # interpolate from the nodal values to the point values using linear interpolation
        P = bar.shape[0]
        sS = sparse.dok_matrix((P, msh.nodes.nr))
        node_idx = msh.elm.node_number_list[idx] - 1
        sS[np.arange(P)[:, None], node_idx] = bar
        sS = sparse.csr_matrix(sS)
        res = sS @ M

        return res

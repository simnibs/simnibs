import os
import time
import h5py
import logging
import numpy as np
import nibabel as nib

from numpy.linalg import pinv
from scipy import sparse
from scipy.ndimage import zoom

from .fem import FEMSystem, get_dirichlet_node_index_cog
from .sim_struct import SimuList
from .cond import standard_cond
from ..mesh_tools import Msh, mesh_io, read_msh
from ..utils.simnibs_logger import logger
from ..utils.utils_numba import map_coord_nn, map_coord_lin, sumf, bincount_nb
from simnibs.simulation import pardiso


class OnlineFEM:
    """
    OnlineFEM class for fast FEM calculations using the Pardiso solver of the MKL library

    Parameters
    ----------
    """
    def __init__(self, mesh, method, roi, anisotropy_type="scalar", solver_options=None, fn_results=None,
                 useElements=True, grid_spacing=None, order=1, fn_coil=None, didtmax=1,
                 dataType=0):
        """
        Constructor of the OnlineFEM class
        """
        self.grid_spacing = grid_spacing            # grid spacing/interval for A field interpolation
        self.order = order                          # the order of the dadt interpolation
        self.method = method                        # 'TES' or 'TMS'
        self.dataType = dataType                    # Calc. magn. of e-field for dataType=0 otherwise return Ex, Ey, Ez
        self.anisotropy_type = anisotropy_type      # "scalar", "vn"
        self.b = None                               # rhs
        self.v = None                               # electric potential
        self.solver = None                          # solver
        self.fn_coil = fn_coil                      # filename of TMS coil (.nii)
        self.didtmax = didtmax                      # Maximum rate of change of coil current (A/us)
        self.roi = roi                              # list of ROI instances
        self.fn_results = fn_results                # name of output results file (.hdf5)

        if type(self.roi) is not list:
            self.roi = [self.roi]

        self.n_roi = len(self.roi)
        self.e = [0 for _ in range(self.n_roi)]
        self.dadt = None

        # True: interpolate the dadt field using positions (coordinates) of elements centroids
        # False: interpolate the dadt field using positions (coordinates) of nodes
        self.useElements = useElements

        if solver_options is None:
            self.solver_options = "pardiso"
        else:
            self.solver_options = solver_options

        # creating logger
        self.logger = setup_logger(os.path.join(os.path.split(fn_results)[0], "simnibs_simulation_" + time.strftime("%Y%m%d-%H%M%S")))

        # read mesh or store in self
        if type(mesh) is str:
            self.mesh = mesh_io.read_msh(mesh)
        elif type(mesh) == Msh:
            self.mesh = mesh
        else:
            raise TypeError("'mesh' parameter has to be either path to .msh file or SimNIBS mesh object.")

        # get Dirichlet node index where V=0 (close to center of gravity of mesh but not on a surface)
        self.dirichlet_node = get_dirichlet_node_index_cog(mesh=self.mesh)

        # For TMS we use only tetrahedra (elm_type=4). The triangles (elm_type=3) are removed.
        if self.method == "TMS":
            self.mesh = remove_triangles_from_mesh(self.mesh)

        # prepare solver and set self.solver
        self._set_matrices_and_prepare_solver()

        # prepare coil for TMS
        if method == "TMS":
            if fn_coil is None:
                raise AssertionError("Please provide filename of TMS coil (.nii) file.")
            self.coil = Coil(filename=self.fn_coil,
                             grid_spacing=self.grid_spacing,
                             didtmax=self.didtmax,
                             logger=self.logger)
            self.logger.info(f'Loaded coil from file: {self.fn_coil}')

        else:
            self.coil = None

    def update_field(self, electrodes=None, currents=None, matsimnibs=None, sim_idx_hdf5=0):
        """
        Calculating and updating electric field for given coil position (matsimnibs) for TMS or electrode position (TES)

        Parameters
        ----------
        electrodes : list of list of np.ndarray [n_sims][n_electrodes]
            List of list of the surface tags or nodes where the currents to be applied for multiple simulations.
        currents : list of list of np.ndarray [n_sims][n_electrodes]
            List of list of the currents in each surface for multiple simulations.
        matsimnibs : np.array of float [4 x 4 x n_sim]
            Tensor containing the coil positions and orientations in SimNIBS space for multiple simulations.
        sim_idx_hdf5 : int
            Simulation index, to continue to write in .hdf5 file

        Returns
        -------
        e : list of lists [n_sim][n_roi] containing np.arrays of float
            Electric field for queried simulations in ROIs.
        """

        if self.method == "TMS":
            if matsimnibs.ndim < 3:
                matsimnibs = matsimnibs[:, :, np.newaxis]
            n_sim = matsimnibs.shape[2]
        else:
            n_sim = len(electrodes)

        self.e = [[0 for _ in range(self.n_roi)] for _ in range(n_sim)]

        for i_sim in range(n_sim):
            start = time.time()

            # determine RHS
            if self.method == "TMS":
                self.b = self.set_rhs(matsimnibs=matsimnibs[:, :, i_sim])
            elif self.method == "TES":
                self.b = self.set_rhs(electrodes=electrodes[i_sim],
                                      currents=currents[i_sim])

            # solve for potential
            self.v = self.solve(b=self.b)

            # calculate e-field
            for i_roi, r in enumerate(self.roi):
                self.e[i_sim][i_roi] = r.calc_fields(v=self.v, dadt=self.dadt, dataType=self.dataType)

                # store results (overwrite if existing)
                with h5py.File(self.fn_results, "a") as f:
                    try:
                        if f"{sim_idx_hdf5:04}" in f[f"e/roi_{i_roi}"].keys():
                            del f[f"e/roi_{i_roi}/{sim_idx_hdf5:04}"]
                    except KeyError:
                        pass
                    f.create_dataset(name=f"e/roi_{i_roi}/{sim_idx_hdf5:04}", data=self.e[i_sim][i_roi])

            sim_idx_hdf5 += 1

            stop = time.time()

            self.logger.info(f"Finished simulation #{i_sim + 1}/{n_sim} (time: {(stop-start):.3f}s).")

        return self.e

    def set_rhs(self, electrodes=None, currents=None, matsimnibs=None):
        """
        Set up right hand side (force vector) of equation system.

        Parameters
        ----------
        electrodes : list of np.ndarray [n_electrodes]
            list of the surface tags or nodes where the currents to be applied.
            WARNING: should NOT include the ground electrode
        currents : list of np.ndarray [n_electrodes]
            list of the currents in each surface
        matsimnibs : np.array of float [4 x 4]
            Matrix containing the coil position and orientation in SimNIBS space.

        Returns
        -------
        b : np.array of float [n_nodes - 1]
            Right hand side of equation system (without Dirichlet node)
        """
        if self.method == "TES":
            b = self.fem.assemble_tdcs_neumann_rhs(electrodes=electrodes,
                                                   currents=currents,
                                                   input_type='node',
                                                   areas=self.msh_nodes_areas)
        elif self.method == "TMS":
            self.dadt = calculate_dadt(a_affine=self.coil.a_affine,
                                       matsimnibs=matsimnibs,
                                       coordinates=self.coordinates,
                                       a_field=self.coil.a_field,
                                       order=self.order,
                                       node_numbers=self.mesh.elm.node_number_list,
                                       useElements=self.useElements)

            b = assemble_force_vector(force_integrals=self.force_integrals,
                                      reshaped_node_numbers=self.reshaped_node_numbers,
                                      dadt=self.dadt,
                                      dirichlet_node=self.dirichlet_node)
        else:
            raise NotImplementedError("Simulation method not implemented yet. Method is either 'TMS' or 'TES'.")

        return b

    def solve(self, b):
        """
        Solve system of equations Ax=b and add Dirichlet node (V=0) to solution.

        Parameters
        ----------
        b : np.array of float [n_nodes - 1]
            Right hand side of equation system (without Dirichlet node)

        Returns
        -------
        v : np.array of float [n_nodes]
            Solution (including the Dirichlet node at the right position)
        """

        logger.disabled = True

        # solve system
        x = self.solver.solve(b)

        # add Dirichlet node (v=0) to solution (Dirichlet node is defined with node indexing starting with 1)
        v = np.insert(x, self.dirichlet_node - 1, 0)

        logger.disabled = False

        return v

    def _set_matrices_and_prepare_solver(self):
        """
        Set matrices and initialize the pardiso solver for update_position in self.solver.
        """

        if self.method == "TMS":
            node_numbers = self.mesh.elm.node_number_list  # self.node_numbers
            useElements = self.useElements
            node_coordinates = self.mesh.nodes.node_coord.T
            tag1 = self.mesh.elm.tag1

            # get the number of nodes
            number_of_nodes = node_coordinates.shape[1]

            # get the coordinates of nodal points/element centers to prepare for the calculation of dadt.
            # Use self.coordinates to interpolate the field in calculate_dadt().
            self.coordinates = get_coordinates(node_coordinates, node_numbers, useElements)

            # set the reshaped_node_numbers
            self.reshaped_node_numbers = (node_numbers - 1).T.reshape(-1)

            # get the distances between local node [0] and nodes [1], [2] and [3] in each element
            local_dist, _ = get_local_distances(node_numbers, node_coordinates)

            # calculate the volume of a tetrahedron given the coordinates of its four nodal points
            volume = np.abs(np.linalg.det(local_dist)) / 6.

            # get gradient operator
            self.gradient = get_gradient(local_dist)

            # get the conductivity values of each tetrahedra
            # tag1 indexes starting from 1, the array cond_table indexes starting from 0
            # Use isotropic conductivities
            conductivity = np.array([float(c.value) for c in standard_cond()])[tag1 - 1]

            # set the force integrals. We use the force integrals to assemble the right hand side force vector
            self.force_integrals = get_force_integrals(volume, self.gradient, conductivity)

            #  assemble the left hand side stiffness matrix
            stiffmat = assemble_stiffness_matrix(volume, self.gradient, conductivity, node_numbers, number_of_nodes, self.dirichlet_node)

            self.solver = pardiso.Solver(stiffmat)

        elif self.method == "TES":
            # Calculate node areas for whole mesh
            self.msh_nodes_areas = self.mesh.nodes_areas()

            # prepare FEM
            self.simulist = SimuList(mesh=self.mesh)
            self.simulist.anisotropy_type = self.anisotropy_type
            cond = self.simulist.cond2elmdata()
            self.fem = FEMSystem.tdcs_neumann(mesh=self.mesh,
                                              cond=cond,
                                              ground_electrode=self.dirichlet_node,
                                              solver_options=self.solver_options,
                                              input_type='node')
            self.fem.prepare_solver()
            self.solver = self.fem._solver


class Coil:
    """
    Coil class for TMS calculations. Contains precalculated magnetic vector potential and resolution information.

    Parameters
    ----------
    filename : str
        Path to the coil file (.ccd or .nii format)
    grid_spacing : float, optional, default: None
        Grid spacing/interval for A field interpolation
    didtmax : float
        Maximum rate of change of coil current (A/us)
    logger : logger object
        Logger

    Attributes
    ----------
    a_affine : np.array of float [4 x 4]
        Affine matrix describing resolution, location and orientation of untransformed magnetic vector potential
    a_field : np.array of float [3 x N_x x N_y x N_z]
        Magnetic vector potential (A_x, A_y, A_z)

    """
    def __init__(self, filename, grid_spacing=None, didtmax=1, logger=None):
        """
        Constructor of Coil class
        """

        self.a_field, self.a_affine = load_A_from_coil_file(coil_file=filename,
                                                            grid_spacing=grid_spacing,
                                                            logger=logger)
        self.didtmax = didtmax
        self.filename = filename


def load_A_from_coil_file(coil_file, grid_spacing, logger=None):
    """
    Load coil file.

    Parameters
    ----------
    coil_file : str
        Coil file (.ccd or .nii format)
    grid_spacing : float, optional, default: None
        Grid spacing/interval for A field interpolation
    logger : logger object, optional, default: None
        Logger

    Returns
    -------
    a_affine : np.array of float [4 x 4]
        Affine matrix describing resolution, location and orientation of untransformed magnetic vector potential
    a_field : np.array of float [3 x N_x x N_y x N_z]
        Magnetic vector potential (A_x, A_y, A_z)
    """

    if not isinstance(coil_file, str):
        raise NameError(
            'Failed to parse input volume (not string or nibabel nifti1 volume)')

    # load the A field and affine matrix
    if coil_file.endswith('.nii.gz') or coil_file.endswith('.nii'):
        # update the A field and affine matrix
        coil = nib.load(coil_file)
        A = np.moveaxis(np.asanyarray(coil.dataobj), -1, 0)
        affine = coil.affine
    else:
        raise ValueError('Coil file must be a nifti file')

    if grid_spacing:

        # interpolate the A field and recalculate the affine matrix if grid_spacing is not an empty array (None)
        start = time()

        grid_spacing = np.array(grid_spacing, dtype='float64')

        # Linear interpolation on the A field.
        A, affine, zoom_factor, zoom_order = zoom_a_field(A, affine, grid_spacing)

        if logger is not None:
            logger.info('Interpolate A field with zoom: {}s (zoom_factor = {}, zoom_order = {})'.format(time() - start, zoom_factor, zoom_order))

    if logger is not None:
        logger.debug('Load A from coil file: ' + coil_file)

    return np.asfortranarray(A), affine


def zoom_a_field(data, affine, grid_spacing, zoom_order=1):
    """
    Produce a denser regular grid based on interpolating the original data.

    Parameters
    ----------
    data : np.array of float [3 x N_x x N_y x N_z]
        Data to zoom (e.g. magnetic vector potential (A_x, A_y, A_z))
    affine : np.array of float [4 x 4]
        Affine matrix describing resolution, location and orientation of untransformed magnetic vector potential
    grid_spacing : float, optional, default: None
        Grid spacing/interval for A field interpolation
    zoom_order : int
        The order of the spline interpolation. The order has to be in the range 0-5. zoom_order=1 for bilinear.

    Returns
    -------
    data_zoomed : np.array of float [3 x N_x x N_y x N_z]
        Zoomed data.
    affine_zoomed : np.array of float [4 x 4]

    zoom_factor : np.array of float [3]
        zoom factor along the axes.
    zoom_order : int
        The order of the spline interpolation. The order has to be in the range 0-5. zoom_order=1 for bilinear.
    """

    # get the size of the A field
    size = np.array(data.shape[1:4], dtype='float64')

    # the scaling factors in the zooming matrix
    dense_factor = affine.diagonal()[0:3] / grid_spacing

    # calculate the size of the A field on denser grid
    size_zoomed = (size - 1.) * dense_factor + 1

    # calculate the zoom factor along the axes.
    # One value for each axis.
    zoom_factor = size_zoomed/size

    # initialize the array
    data_zoomed = np.zeros(np.rint(np.concatenate(([data.shape[0]], size_zoomed))).astype(np.longlong))

    # interpolate uniformly-spaced 3D array on a finer uniformed-spacing grid
    # Use zoom_order=1 for bilinear
    for i in range(data.shape[0]):
        data_zoomed[i] = zoom(data[i, ...], zoom_factor, order=zoom_order)

    zoom_matrix = np.diag(np.append(np.array([1., 1., 1.])/dense_factor, 1))
    affine_zoomed = np.matmul(affine, zoom_matrix)

    return data_zoomed, affine_zoomed, zoom_factor, zoom_order


def calculate_dadt(a_affine, matsimnibs, coordinates, a_field, order, node_numbers, useElements):
    """
    Interpolates the dadt field from a coil file.

    Parameters
    ----------
    a_affine : np.array of float [4 x 4]
        Affine matrix describing resolution, location and orientation of untransformed magnetic vector potential
    matsimnibs : np.array of float [4 x 4]
        Matrix containing the coil positions and orientations in SimNIBS space.
    coordinates :

    a_field : np.array of float [3 x N_x x N_y x N_z]
        Magnetic vector potential (A_x, A_y, A_z)
    order :

    node_numbers : np.array of size [n_elements x 4]
        Node number list (connectivity list)
    useElements : bool
        Get coordinates in the element centers (True) or of the nodes (False)

    Returns
    -------
    dadt : np.array of float [3 x N_x x N_y x N_z]
        Magnetic vector potential (A_x, A_y, A_z)
    """

    if order == 0:
        fun_interp = map_coord_nn
    elif order == 1:
        fun_interp = map_coord_lin
    else:
        raise NotImplementedError('simnibs.simulation.examples.map_coord_numba currently requires order<=1')

    # Get the affine transformation from the "coordinates" to the coil space (defined in coil file)
    trans = np.dot(pinv(a_affine), pinv(matsimnibs))

    # gets the coordinates in voxel space
    pos = trans[:3, :3] @ coordinates + trans[:3, [3]]

    # Interpolates the values of the field in the given coordinates
    # ensure fortran order - then code is faster - the 32 bit does not really help
    dadt_interp = np.empty((3, pos.shape[1]), dtype='float64', order='F')

    # Interpolates the values of the field in the given coordinates
    # run once with dummy pos (0) to initialize function:
    fun_interp(a_field, pos, dadt_interp)

    # Rotates the field
    dadt_rotated = (matsimnibs[:3, :3] @ dadt_interp)

    # Interpolates the field from a nifti file
    if useElements:
        # if useElements == True, no interpolation because dadt_rotated is defined on elements
        dadt = dadt_rotated

    else:
        # if useElements == False, interpolate dadt_rotated from nodes to element centroids
        dadt = calculate_element_centers(dadt_rotated, node_numbers)

    return dadt.T


def assemble_force_vector(force_integrals, reshaped_node_numbers, dadt, dirichlet_node):
    """
    Assembly of the force vector in a system of linear equations stiffmat * x = forcevec. for TMS.

    Parameters
    ----------
    force_integrals : np.array of float [4, number_of_elements, 3]
        Force_integrals for rhs calculation derived by volume * conductivity * gradient
    reshaped_node_numbers :


    dadt: NodeData or ElementData
        dA/dt field at each node or element
    dirichlet_node : int
        Index of the Dirichlet node (defined with node indexing starting with 1)

    Returns
    -------
    forcevec: np.array
        Right-hand side
    """

    # integrate in each node of each element, the value for repeated nodes will be summed
    # together later
    node_integrals = np.zeros(force_integrals.shape[:2], order='C')

    # force_integrals: (4, number_of_elements, 3); dadt: (number_of_elements, 3); node_integrals: (4, number_of_elements)
    # node_integrals = np.einsum('ijk,jk->ij', force_integrals, dadt)
    sumf(force_integrals, dadt, node_integrals)

    # Assembles the right hand side for TMS.
    # forcevec: (number_of_nodes,), reshaped_node_numbers: (number_of_elements*4,), node_integrals: (4, number_of_elements)
    # forcevec = bincount_nb(reshaped_node_numbers, node_integrals.reshape(-1)).reshape(-1, 1)

    # keep np.bincount to make the testing more stable. May change it back to bincount_nb() for performance reason.
    forcevec = np.bincount(reshaped_node_numbers, node_integrals.reshape(-1)).reshape(-1, 1)

    # Applies the dirichlet BC to the right hand side force_vector
    forcevec = np.delete(forcevec, dirichlet_node-1)

    return forcevec


def calculate_element_centers(node_coordinates, node_numbers):
    """
    Calculate the center of each element. The center is the average of the vertices:
    Centroid = (a + b + c + d) / 4.

    Parameters
    ----------
    node_coordinates : np.array of size [3 x n_nodes]
        Coordinates of the nodes (x, y, z)
    node_numbers : np.array of size [n_elements x 4]
        Node number list (connectivity list)

    Returns
    -------
    element_centers : np.array of float [3 x n_elements]
        Coordinates (x, y, z) of the element centers
    """
    # node_coordinates: (number_of_nodes, 3)
    # node_numbers: (number_of_elements, 4)
    index = node_numbers - 1

    # calculate a/4, b/4, c/4 and d/4 first
    coordinates_normalized = node_coordinates * (1./index.shape[1])

    # locate the coordinates of nodes in each element
    # (3, number_of_elements, 4)
    pos = coordinates_normalized[:, index]

    # then calculate the sum of a/4, b/4, c/4 and d/4
    # (3, number_of_elements)
    element_centers = np.einsum('ijk->ij', pos)

    return element_centers


def get_force_integrals(volume, gradient, conductivity):
    """
    Assembly of the integration part in the force vector without dadt for TMS and conductivity.ndim == 1
    in physical space (unit: cubic metre)

    Parameters
    ----------
    volume : np.array of float [n_elements]
        Volume of the tetrahedra
    gradient : np.array of size [n_elements, 4, 3]
        Gradient in each tetrahedra
    conductivity : np.array of float [n_elements]
        Electrical conductivity value (in S/m) assigned to each tetrahedra

    Returns
    -------
    force_integrals : np.array of float [4, number_of_elements, 3]
        Force_integrals for rhs calculation derived by volume * conductivity * gradient
    """

    if conductivity.ndim != 1:
        raise ValueError('Invalid conductivity array')

    # volume and conductivity are 1D array shape=(M,)
    vol_cond = volume * conductivity * 1e-6 # 1e-6 is to convert 'mm' to 'm'

    # calculate the force_integrals = volume * conductivity * gradient, and rearrange the dimension to 4xMx3
    force_integrals = np.swapaxes(-vol_cond[:, None, None] * gradient, 0, 1)  # (4, number_of_elements, 3)

    return force_integrals


def assemble_stiffness_matrix(volume, gradient, conductivity, node_numbers, number_of_nodes, dirichlet_node):
    """
    Assembly of the l.h.s stiffness matrix. Based in the OptVS algorithm in Cuvelier et. al. 2016.

    Cuvelier, F., Japhet, C., & Scarella, G. (2016). An efficient way to assemble finite element matrices in
    vector languages. BIT Numerical Mathematics, 56(3), 833-864.

    Parameters
    ----------
    volume : np.array of float [n_elements]
        Volume of the tetrahedra
    gradient : np.array of size [n_elements, 4, 3]
        Gradient in each tetrahedra
    conductivity : np.array of float [n_elements]
        Electrical conductivity value (in S/m) assigned to each tetrahedra
    node_numbers : np.array of size [n_elements x 4]
        Node number list (connectivity list)
    number_of_nodes : int
        Number of nodes
    dirichlet_node : int
        Index of the Dirichlet node (defined with node indexing starting with 1)

    Returns
    -------
    stiffmat : scipy.sparse [(n_nodes-1) x (n_nodes-1)]
        Stiffness matrix in sparse (CSR) format (without Dirichlet node)
    """

    # gradient: (number_of_elements, 4, 3),
    # volume: (number_of_elements,),
    # conductivity: (number_of_elements,),
    # node_numbers: (number_of_elements, 4),
    # number_of_nodes: number_of_nodes
    # (number_of_elements, 4)
    index = node_numbers - 1

    dim = np.arange(index.shape[1])
    idx = np.repeat(dim, len(dim))
    idy = np.tile(dim, len(dim))

    # Simplify the integration using commutative law of dot product (elementary-wise product in this case
    # units == 'mm': * 1e6 from the gradient operator, 1e-9 from the volume
    factor = (volume * conductivity * 1e-3)[:, None, None] * gradient

    stiffmat = sparse.coo_matrix((number_of_nodes, number_of_nodes), dtype='float64')

    stiffmat.data = (factor @ np.swapaxes(gradient, 1, 2)).flatten()
    stiffmat.row = (index[:, idx]).flatten()
    stiffmat.col = (index[:, idy]).flatten()

    stiffmat = stiffmat.tocsr()
    stiffmat.eliminate_zeros()

    # Make stiffmat symmetric if necessary
    # stiffmat = (stiffmat + stiffmat.T) * 0.5
    stiffmat = stiffmat.sorted_indices()

    # remove row and column for the dirichlet boundary condition
    stiffmat = delete_row_csr(stiffmat, dirichlet_node-1)
    stiffmat = delete_cols_csr(stiffmat, dirichlet_node-1)

    return stiffmat


def get_coordinates(node_coordinates, node_numbers, useElements):
    """
    Calculate the centers of the elements.
    If useElements is true, get the coordinates of element centers; otherwise get the coordinates of nodes

    Parameters
    ----------
    node_coordinates : np.array of size [3 x n_nodes]
        Coordinates of the nodes (x, y, z)
    node_numbers : np.array of size [n_elements x 4]
        Node number list (connectivity list)
    useElements : bool
        Get coordinates in the element centers (True) or of the nodes (False)

    Returns
    -------
    coordinates : np.array of float of size [n_element] or [n_nodes]
        Centers of the elements or node coordinates.
    """

    # useElements = True: get the coordinates of the element centers
    if useElements:
        # calculate the coordinates of the element centers
        # Notice coordinates: input and output are the same
        coordinates = calculate_element_centers(node_coordinates, node_numbers)

    else:
        # get the coordinates of the nodes
        coordinates = node_coordinates

    return coordinates


def get_local_distances(node_numbers, node_coordinates):
    """
    Get the distances between local node [0] and nodes [1], [2] and [3] in each element

    Parameters
    ----------
    node_numbers : np.array of size [number_of_elements, 4]
        Node number list (connectivity list of tetrahedra)
    node_coordinates : np.array of size [3 x n_nodes]
        Coordinates of the nodes (x, y, z)

    Returns
    -------
    local_dist : np.array of size [number_of_elements, 3, 3]
        Distances between local node [0] and nodes [1], [2] and [3] in each element
    local_node_coords : np.array of size [number_of_elements, 3, 3]
        Node coordinates
    """

    # get the local node coordinates in each element
    # node_coordinates: (3, number_of_nodes)
    # node_numbers: (number_of_elements, 4)
    # local_node_coords: (3, number_of_elements, 4)
    local_node_coords = node_coordinates[:, node_numbers - 1]

    # get the distances between local node [0] and nodes [1], [2] and [3] in each element
    # local_dist: (3, number_of_elements, 3)
    local_dist = local_node_coords[:, :, 1:] - local_node_coords[:, :, [0]]

    # out: (number_of_elements, 3, 3)
    return np.moveaxis(local_dist, 0, -1), (local_node_coords[:, :, 0].T)[:, :, None]


def get_gradient(local_dist):
    """
    Calculate the gradient of a function in each tetrahedra.

    local_dist * gradient = project_matrix

    project_matrix is a projection matrix
    project_matrix = [-1, 1, 0, 0]
                     [-1, 0, 1, 0]
                     [-1, 0, 0, 1]
    and local_dist is the distances between local node [0] and nodes [1], [2] and [3] in each element.

    Parameters
    ----------
    local_dist : np.array of size [n_elements, 3, 3]
        Distances between local node [0] and nodes [1], [2] and [3] in each element (derived from get_local_distances)

    Returns
    -------
    gradient : np.array of size [n_elements, 4, 3]
        Gradient in each tetrahedra
    """

    # define projection matrix
    project_matrix = np.hstack([-np.ones((3, 1)), np.eye(3)])

    # solve local_dist * gradient = project_matrix.
    gradient = np.linalg.solve(local_dist, project_matrix[None, :, :])

    # swapaxes
    return np.swapaxes(gradient, 1, 2)


def remove_triangles_from_mesh(mesh):
    """
    Load the mesh file. For the E field calculation, we use only tetrahedra (elm_type=4).
    The triangles (elm_type=3) are removed.

    Parameters
    ----------
    mesh : Msh object
        Mesh object

    Returns
    -------
    mesh : Msh object
        Mesh object without triangles
    """

    # keep the tetrahedra in the mesh file (triangles are removed)
    mesh = mesh.crop_mesh(elm_type=4)

    assert mesh.elm.node_number_list.max() == mesh.nodes.nr

    # get the nodal index and value with dirichlet boundary condition
    index_naught = mesh.nodes.node_coord[:, 2].argmin()

    # switch the node with dirichlet bc and the last node
    mesh.nodes.node_coord[[index_naught, -1]] = mesh.nodes.node_coord[[-1, index_naught]]

    # update the node_number_list
    mask1 = np.nonzero(mesh.elm.node_number_list == index_naught + 1)
    mask2 = np.nonzero(mesh.elm.node_number_list == mesh.nodes.nr)
    mesh.elm.node_number_list[mask1] = mesh.nodes.nr
    mesh.elm.node_number_list[mask2] = index_naught + 1

    return mesh


def get_mesh_with_tetrahedra(mesh_file, logger=None):
    """
    Load the mesh file. For the E field calculation, we use only tetrahedra (elm_type=4).
    The triangles (elm_type=3) are removed.

    Parameters
    ----------
    mesh_file : str
        Filename (incl. path) to .msh file.
    logger : logger object
        Logger

    Returns
    -------
    mesh : Msh object
        Mesh object
    """
    # read_mesh: load the mesh file; crop_mesh: keep the tetrahedra in the mesh file (triangles are removed)
    mesh = read_msh(mesh_file).crop_mesh(elm_type=4)

    assert mesh.elm.node_number_list.max() == mesh.nodes.nr

    # get the nodal index and value with dirichlet boundary condition
    index_naught = mesh.nodes.node_coord[:, 2].argmin()

    # switch the node with dirichlet bc and the last node
    mesh.nodes.node_coord[[index_naught, -1]] = mesh.nodes.node_coord[[-1, index_naught]]

    # update the node_number_list
    mask1 = np.nonzero(mesh.elm.node_number_list == index_naught + 1)
    mask2 = np.nonzero(mesh.elm.node_number_list == mesh.nodes.nr)
    mesh.elm.node_number_list[mask1] = mesh.nodes.nr
    mesh.elm.node_number_list[mask2] = index_naught + 1

    if logger is not None:
        logger.info('Loaded mesh file: ' + mesh_file)

    return mesh


def setup_logger(logname, filemode='w', format='[ %(name)s ] %(levelname)s: %(message)s', datefmt='%H:%M:%S'):
    """
    Setup logger.

    Parameters
    ----------
    logname : str
        Filename of logfile
    filemode : str, optional, default: 'w'
        'a' append or 'w' overwrite existing logfile.
    format : str, optional, default: '[ %(name)s ] %(levelname)s: %(message)s'
        format of the output message.
    datefmt : str, optional, default: '%H:%M:%S'
        Format of the output time.

    Returns
    -------
    logger : logger instance
        Logger using loging module
    """

    logging.basicConfig(filename=logname,
                        filemode=filemode,
                        format=format,
                        datefmt=datefmt,
                        level=logging.DEBUG)

    logger = logging.getLogger("simnibs")

    return logger


def delete_row_csr(mat, i):
    """
    Delete rows from a sparse matrix in Compressed Sparse Row (CSR) format.

    Parameters
    ----------
    mat : scipy.sparse matrix in csr format
        Sparse matrix (CSR)
    i : int or list of int
        Indices of rows to delete from mat.

    Returns
    -------
    mat : scipy.sparse matrix in csr format
        Sparse matrix with deleted rows.
    """
    if not isinstance(mat, sparse.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    n = mat.indptr[i+1] - mat.indptr[i]
    if n > 0:
        mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
        mat.data = mat.data[:-n]
        mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
        mat.indices = mat.indices[:-n]
    mat.indptr[i:-1] = mat.indptr[i+1:]
    mat.indptr[i:] -= n
    mat.indptr = mat.indptr[:-1]
    mat._shape = (mat._shape[0]-1, mat._shape[1])

    return mat


def delete_cols_csr(mat, i):
    """
    Delete columns from a sparse matrix in Compressed Sparse Row (CSR) format.

    Parameters
    ----------
    mat : scipy.sparse matrix in csr format
        Sparse matrix (CSR)
    i : int or list of int
        Indices of columns to delete from mat.

    Returns
    -------
    mat : scipy.sparse matrix in csr format
        Sparse matrix with deleted columns.
    """
    if not isinstance(mat, sparse.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")

    if type(i) is int:
        i = [i]

    mask = np.ones(mat.shape[1], dtype=bool)
    mask[i] = False

    return mat[:, mask]

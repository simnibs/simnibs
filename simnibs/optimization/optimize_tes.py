import copy
import os
import time
import copy
import numpy as np
import nibabel as nib
import scipy.ndimage.morphology as mrph

from scipy.spatial import ConvexHull, convex_hull_plot_2d
from numpy.linalg import eig, inv

from ..mesh_tools import Msh
from ..mesh_tools import cgal
from ..mesh_tools import mesh_io
from ..mesh_tools import surface
from ..simulation.sim_struct import ELECTRODE
from ..utils.simnibs_logger import logger
from ..utils.file_finder import Templates, SubjectFiles
from ..utils.transformations import subject2mni_coords

SIMNIBSDIR = os.path.split(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))[0]

# TODO: in sim_struct -> filefinder to get filenames of EEG cap for example
# TODO: implement new tes fast optimize class
# TODO: adapt matlab_tools/opt_struct.m and include new class
# 1:1 Mapping der Elektroden -> von Neumann Kanal 1 -> Sink 1, Kanal 2 -> Sink 2
# 1 Kanal -> mehrere sink elektroden -> fake Dirichlet -> Kanal 1 -> Sink 1,2,3 (zusammengeschaltet) gleiche Spannung


class TESoptimize():
    ''' Defines a TES optimization problem using a direct approach

    Parameters
    --------------
    electrode : Electrode Object
        Electrode object containing ElectrodeArray instances
        (see /simulation/array_layout.py for pre-implemented examples)
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float (optional)
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int (optional)
        Maximum number of active electrodes. Default: no maximum
    init_pos : str or list of str or list of str and np.array of float [3]
        Initial positions of movable Electrode arrays (for each movable array)
    fn_eeg_cap : str, optional, default: 'EEG10-10_UI_Jurak_2007.csv'
        Filename of EEG cap to use for initial position (without path)
        - 'EEG10-10_UI_Jurak_2007.csv'
        - 'easycap_BC_TMS64_X21.csv'
        - 'EEG10-10_Cutini_2011.csv'
        - 'EEG10-10_Neuroelectrics.csv'
        - 'EEG10-20_extended_SPM12.csv'
        - 'EEG10-20_Okamoto_2004.csv'
    plot : bool, optional, default: False
        Plot configurations in output folder for visualization and control

    name: str (optional)
        Name of optimization problem. Default: optimization
    target: list of TDCStarget objects (optional)
        Targets for the optimization. Default: no target
    avoid: list of TDCSavoid objects
        list of TDCSavoid objects defining regions to avoid

    Attributes
    --------------
    electrode : Electrode Object
        Electrode object containing ElectrodeArray instances
        (see /simulation/array_layout.py for pre-implemented examples)
    nodes_areas : np.array of float [n_nodes_skin]
        Areas of skin nodes
    nodes_normals : np.array of float [n_nodes_skin x 3]
        Normals of skin nodes
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int
        Maximum number of active electrodes. Default: no maximum



    The two above are used to define:

    mesh: simnibs.msh.mesh_io.Msh
        Mesh with problem geometry

    leadfield: np.ndarray
        Leadfield matrix (N_elec -1 x M x 3) where M is either the number of nodes or the
        number of elements in the mesh. We assume that there is a reference electrode

    Alternatively, you can set the three attributes above and not leadfield_path,
    mesh_path and leadfield_hdf

    lf_type: None, 'node' or 'element'
        Type of leadfield.

    name: str
        Name for the optimization problem. Defaults tp 'optimization'

    target: list of TDCStarget objects
        list of TDCStarget objects defining the targets of the optimization

    avoid: list of TDCSavoid objects (optional)
        list of TDCSavoid objects defining regions to avoid

    open_in_gmsh: bool (optional)
        Whether to open the result in Gmsh after the calculations. Default: False

    Warning
    -----------
    Changing leadfield_hdf, leadfield_path and mesh_path after constructing the class
    can cause unexpected behaviour
    '''

    def __init__(self,
                 msh=None,
                 electrode=None,
                 roi=None,
                 init_pos=None,
                 fn_eeg_cap=None,
                 plot=False,
                 max_total_current=2e-3,
                 max_individual_current=1e-3,
                 max_active_electrodes=None,
                 output_folder=None,
                 target=None,
                 avoid=None):

        if init_pos is None:
            init_pos = ["C3"]

        if type(init_pos) is str:
            init_pos = list(init_pos)

        if type(init_pos) is np.ndarray:
            init_pos = [init_pos]

        if fn_eeg_cap is None:
            self.fn_eeg_cap = 'EEG10-10_UI_Jurak_2007.csv'
        else:
            self.fn_eeg_cap = os.path.split(fn_eeg_cap)[1]

        assert type(output_folder) is str, "Please prove an output folder to save optimization results in."

        self.electrode = electrode
        self.n_ele_free = len(electrode.electrode_arrays)
        self.max_total_current = max_total_current
        self.max_individual_current = max_individual_current
        self.max_active_electrodes = max_active_electrodes
        self.roi = roi
        self.init_pos = init_pos
        self.init_pos_subject_coords = []
        self.output_folder = output_folder
        self.plot_folder = os.path.join(self.output_folder, "plots")
        self.plot = plot
        self.fn_results_hdf5 = os.path.join(self.output_folder, "opt.hdf5")
        init_pos_list = ["C3", "C4"]

        assert len(init_pos) == self.n_ele_free, "Number of initial positions has to match number of freely movable" \
                                                 "electrode arrays"

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        if not os.path.exists(self.plot_folder):
            os.makedirs(self.plot_folder)

        # read mesh or store in self
        if type(msh) is str:
            self.msh = mesh_io.read_msh(msh)
        elif type(msh) == Msh:
            self.msh = msh
        else:
            raise TypeError("msh has to be either path to .msh file or SimNIBS mesh object.")

        self.ff_templates = Templates()
        self.ff_subject = SubjectFiles(fnamehead=self.msh.fn)
        self.fn_electrode_mask = self.ff_templates.mni_volume_upper_head_mask  # os.path.join(SIMNIBSDIR, "resources", "templates", "MNI152_T1_1mm_electrode_mask.nii.gz")

        # relabel internal air
        self.msh_relabel = relabel_internal_air(m=self.msh,
                                                subpath=os.path.split(self.msh.fn)[0],
                                                label_skin=1005,
                                                label_new=1099,
                                                label_internal_air=501)

        # create skin surface
        self.skin_surface = surface.Surface(mesh=self.msh_relabel, labels=1005)

        # determine point indices where the electrodes may be applied during optimization
        self.skin_surface = self.valid_skin_region(skin_surface=self.skin_surface, mesh=self.msh_relabel)

        # fit optimal ellipsoid to valid skin points
        self.eli_center, self.eli_radii, self.eli_rotmat = fit_ellipsoid(self.skin_surface.nodes)

        # set initial positions to electrode C3 (and C4) if nothing is provided
        if init_pos is None:
            if self.n_ele_free > 2:
                raise NotImplementedError("Please specify initial coordinates or EEG electrode positions for each"
                                          "freely movable electrode array (init_pos)!")
            self.init_pos = [init_pos_list[i] for i in range(self.n_ele_free)]

        # get subject coordinates of initial positions
        if type(self.init_pos[0]) is str:
            for eeg_pos in self.init_pos:
                tmp = ELECTRODE()
                tmp.centre = eeg_pos
                tmp.substitute_positions_from_cap(cap=self.ff_subject.get_eeg_cap(cap_name=self.fn_eeg_cap))
                self.init_pos_subject_coords.append(tmp.centre)
        else:
            self.init_pos_subject_coords = self.init_pos

        if target is None:
            self.target = []
        else:
            self.target = target
        if avoid is None:
            self.avoid = []
        else:
            self.avoid = avoid

        if self.plot:
            theta = np.linspace(0, np.pi, 180)
            phi = np.linspace(0, 2*np.pi, 360)
            coords_sphere = np.array(np.meshgrid(theta, phi)).T.reshape(-1, 2)
            # coords_sphere = np.hstack([theta[:, np.newaxis], phi[:, np.newaxis]])
            eli_coords = ellipsoid2cartesian(coords=coords_sphere,
                                             center=self.eli_center,
                                             radii=self.eli_radii,
                                             rotmat=self.eli_rotmat,
                                             return_normal=False)

            np.savetxt(os.path.join(self.output_folder, "plots", "fitted_ellipsoid.txt"), eli_coords)

            import pynibs
            # write hdf5 _geo files for visualization in paraview
            pynibs.write_geo_hdf5_surf(out_fn=os.path.join(self.output_folder, "plots", "upper_head_region_geo.hdf5"),
                                       points=self.skin_surface.nodes,
                                       con=self.skin_surface.tr_nodes,
                                       replace=True,
                                       hdf5_path='/mesh')

            pynibs.write_data_hdf5_surf(data=[np.zeros(self.skin_surface.tr_nodes.shape[0])],
                                        data_names=["domain"],
                                        data_hdf_fn_out=os.path.join(self.output_folder, "plots", "upper_head_region_data.hdf5"),
                                        geo_hdf_fn=os.path.join(self.output_folder, "plots", "upper_head_region_geo.hdf5"),
                                        replace=True)

        # gauge problem (find node in COG of head and move to end), modifies self.mesh
        self.gauge_mesh()

        # get_surround_pos

        # TODO: solver_options: use PARDISO solver as standard
        # assemble FEM matrix
        # self.fem=FEMSystem.tdcs_neumann(
        #     self.mesh, elm_cond, self.mesh.node.nr,  # last element is 0 was shifted from center
        #     solver_options=solver_options,
        #     input_type='node'
        # )

    def subject2ellipsoid(self, coords: np.ndarray, normals: np.ndarray):
        """
        Transform coordinates from subject to ellipsoid space
        Line intersection of skin surface point along normal and ellipsoid.
        https://math.stackexchange.com/questions/3309397/line-ellipsoid-intersection

        Parameters
        ----------
        coords : np.array of float [n_points x 3]
            Cartesian coordinates in subject space (x, y, z)
        normals : np.array of float [n_points x 3]
            Surface normals of points

        Returns
        -------
        coords : np.array of float [n_coords x 2]
            Spherical coordinates [theta, phi] in radiant.
        """

        R = self.eli_rotmat
        A = np.diag(np.array(1/(self.eli_radii**2)))
        v = coords - self.eli_center
        B = R @ A @ R.transpose()

        p = np.zeros(3)
        p[0] = normals @ B @ normals.transpose()
        p[1] = 2 * v @ B @ normals.transpose()
        p[2] = v @ B @ v.transpose()

        lam = np.roots(p)

    def ellipsoid2subject(self, coords: np.ndarray):
        """
        Transform coordinates from ellipsoid to subject space

        Parameters
        ----------
        coords : np.array of float [n_points x 2]
            Spherical coordinates [theta, phi] in radiant.

        Returns
        -------
        index : int
        coords : np.array of float [n_coords x 3]
            Cartesian coordinates in subject space (x, y, z)
        """

        # get cartesian coordinates and normal
        x0, n = ellipsoid2cartesian(coords=coords,
                                    center=self.eli_center,
                                    radii=self.eli_radii,
                                    rotmat=self.eli_rotmat,
                                    return_normal=True)

        # define near and far point
        near = (x0 - 100*n)[np.newaxis, :]
        far = (x0 + 100*n)[np.newaxis, :]

        indices, points = cgal.segment_triangle_intersection(
            self.skin_surface.nodes_valid,
            self.skin_surface.tr_nodes_valid,
            near, far
        )

        idx_closest_point = np.argmin(np.linalg.norm(points-x0, axis=1))
        index = indices[idx_closest_point, 1]
        points = points[idx_closest_point, :]

        return index, points

    def valid_skin_region(self, skin_surface, mesh):
        """
        Determine the nodes of the scalp surface where the electrode can be applied (not ears and face etc.)

        Parameters
        ----------
        skin_surface : Surface object
            Surface of the mesh (mesh_tools/surface.py)
        mesh : Msh object
            Mesh object created by SimNIBS (mesh_tools/mesh_io.py)
        """
        # load mask of valid electrode positions (in MNI space)
        mask_img = nib.load(self.fn_electrode_mask)
        mask_img_data = mask_img.get_fdata()

        # transform skin surface points to MNI space
        skin_nodes_mni_ras = subject2mni_coords(coordinates=skin_surface.nodes,
                                                m2m_folder=os.path.split(mesh.fn)[0],
                                                transformation_type='nonl')

        # transform coordinates to voxel space
        skin_nodes_mni_voxel = np.floor(np.linalg.inv(mask_img.affine) @
                                        np.hstack((skin_nodes_mni_ras,
                                                   np.ones(skin_nodes_mni_ras.shape[0])[:, np.newaxis]))
                                        .transpose())[:3, :].transpose().astype(int)
        skin_nodes_mni_voxel[skin_nodes_mni_voxel[:, 0] >= mask_img.shape[0], 0] = mask_img.shape[0] - 1
        skin_nodes_mni_voxel[skin_nodes_mni_voxel[:, 1] >= mask_img.shape[1], 1] = mask_img.shape[1] - 1
        skin_nodes_mni_voxel[skin_nodes_mni_voxel[:, 2] >= mask_img.shape[2], 2] = mask_img.shape[2] - 1
        # skin_nodes_mni_voxel[skin_nodes_mni_voxel < 0] = 0

        # get boolean mask of valid skin points
        skin_surface.mask_valid_nodes = mask_img_data[skin_nodes_mni_voxel[:, 0],
                                                      skin_nodes_mni_voxel[:, 1],
                                                      skin_nodes_mni_voxel[:, 2]].astype(bool)

        # remove points outside of MNI space (lower neck)
        skin_surface.mask_valid_nodes[(skin_nodes_mni_voxel < 0).any(axis=1)] = False

        skin_surface.mask_valid_tr = np.zeros(skin_surface.tr_centers.shape).astype(bool)

        unique_points = np.unique(skin_surface.tr_nodes[skin_surface.mask_valid_nodes[skin_surface.tr_nodes].all(axis=1), :])
        for point in unique_points:
            idx_where = np.where(skin_surface.tr_nodes == point)
            skin_surface.mask_valid_tr[idx_where[0], idx_where[1]] = True
        skin_surface.mask_valid_tr = skin_surface.mask_valid_tr.all(axis=1)

        # determine connectivity list of valid skin region (creates new node and connectivity list)
        skin_surface.nodes, skin_surface.tr_nodes = \
            self.create_new_connectivity_list_point_mask(points=skin_surface.nodes,
                                                         con=skin_surface.tr_nodes,
                                                         point_mask=skin_surface.mask_valid_nodes)

        # identify spurious skin patches inside head and remove them
        tri_domain = np.ones(skin_surface.tr_nodes.shape[0]).astype(int) * -1
        point_domain = np.ones(skin_surface.nodes.shape[0]).astype(int) * -1

        domain = 0
        while (tri_domain == -1).any():
            nodes_idx_of_domain = np.array([])
            tri_idx_of_domain = np.where(tri_domain == -1)[0][0]

            n_current = -1
            n_last = 0
            while n_last != n_current:
                n_last = copy.deepcopy(n_current)
                nodes_idx_of_domain = np.unique(np.append(nodes_idx_of_domain, skin_surface.tr_nodes[tri_idx_of_domain, :])).astype(int)
                tri_idx_of_domain = np.isin(skin_surface.tr_nodes, nodes_idx_of_domain).any(axis=1)
                n_current = np.sum(tri_idx_of_domain)
                # print(f"domain: {domain}, n_current: {n_current}")

            tri_domain[tri_idx_of_domain] = domain
            point_domain[nodes_idx_of_domain] = domain
            domain += 1

        domain_idx_main = np.argmax([np.sum(point_domain == d) for d in range(domain)])

        skin_surface.nodes, skin_surface.tr_nodes = \
            self.create_new_connectivity_list_point_mask(points=skin_surface.nodes,
                                                         con=skin_surface.tr_nodes,
                                                         point_mask=point_domain == domain_idx_main)

        skin_surface.nodes_areas = skin_surface.nodes_areas[skin_surface.mask_valid_nodes]
        skin_surface.nodes_normals = skin_surface.nodes_normals[skin_surface.mask_valid_nodes, :]
        skin_surface.surf2msh_nodes = skin_surface.surf2msh_nodes[skin_surface.mask_valid_nodes, :]
        skin_surface.surf2msh_triangles = skin_surface.surf2msh_triangles[skin_surface.mask_valid_tr]
        skin_surface.tr_areas = skin_surface.tr_areas[skin_surface.mask_valid_tr]
        skin_surface.tr_centers = skin_surface.tr_centers[skin_surface.mask_valid_tr, :]
        skin_surface.tr_normals = skin_surface.tr_normals[skin_surface.mask_valid_tr, :]

        return skin_surface

    def create_new_connectivity_list_point_mask(self, points, con, point_mask):
        """
        Creates a new point and connectivity list when applying a point mask (changes indices of points)

        Parameters
        ----------
        points : np.array of float [n_points x 3]
            Point coordinates
        con : np.array of float [n_tri x 3]
            Connectivity of triangles
        point_mask : nparray of bool [n_points]
            Mask of (True/False) which points are kept in the mesh

        Returns
        -------
        points_new : np.array of float [n_points_new x 3]
            New point array containing the remaining points after applying the mask
        con_new : np.array of float [n_tri_new x 3]
            New connectivity list containing the remaining points (includes reindexing)
        """
        con_global = con[point_mask[con].all(axis=1), :]
        unique_points = np.unique(con_global)
        points_new = points[unique_points, :]

        con_new = np.zeros(con_global.shape).astype(int)

        for i, idx in enumerate(unique_points):
            idx_where = np.where(con_global == idx)
            con_new[idx_where[0], idx_where[1]] = i

        return points_new, con_new

    def assign_skin_points_electrode(self, coords, electrode_array):
        """
        Assigns the skin points of the electrodes in electrode array and writes the points in
        electrode_array.electrodes[i].points and electrode_array.electrodes[i].points_area

        Parameters
        ----------
        coords : np.array of float [3]
            Spherical coordinates (theta, phi) and orientation angle (alpha) of electrode array (in this order)
        electrode_array : ElectrodeArray instance
            ElectrodeArray instance containing Electrode instances
        """
        # get
        index, point = self.ellipsoid2subject(coords)

    def gauge_mesh(self):
        """

        :return:
        """
        pass

    def update_electrode(self, location_parameters):
        """
        Updates
        :return:
        """
        # this function is currently in sim_struct and rather slow
        # self.array_layout.get_surround_pos(center_pos, fnamehead, radius_surround=50, N=4,
        #                  pos_dir_1stsurround=None, phis_surround=None,
        #                  tissue_idx=1005, DEBUG=False)
        self.array_layout_node_idx = 1

    def update_rhs(self):
        """
        Update RHS with new electrode positions
        :return:
        """
        # self.rhs = self.fem.assemble_tdcs_neumann_rhs([np.hstack((self.array_layout_node_idx))], [np.hstack((I))], input_type='node', areas=self.areas)


    # TODO: From Fang (TMS) adapt accordingly
    # TODO: surface can be a ROI class with a .get_fields(v) method, is initialized with either "TMS" (with dAdt or "TES" raw)
    # TODO: change dataType ROI specific in ROI class init
    def update_field(self, matsimnibs, coil_id, amplitude, surfaceIDs, dataType):
        '''Calculate the E field with the different coil positions. The input matsimnibs is a 4x4 affine matrix which defines the coil position in 3D'''

        start = time()
        a_affine = self.coil[coil_id].a_affine
        a_field = self.coil[coil_id].a_field
        didtmax = self.coil[coil_id].didtmax

        assert dataType in [0, 1]



        start = time()

        # update electrode positions
        self.update_electrode()

        # Assemble RHS
        self.update_rhs()

        logger.debug('Assemble the right hand side force vector: {0:.6f}s'.format(time() - start))

        start = time()

        # Solve the linear equation Mx = F. M is the stiffness matrix and F is the force vector.
        # x = self.solve_and_normalize(rhs)

        # TODO: generalize
        # Dirichlet node in case of TMS (last one set to zero)
        v = np.append(x, 0)

        logger.debug('Estimate the x: {0:.2f}s'.format(
            time() - start))

        start = time()

        m = []

        # TODO: this is done in the new ROI class in the .calc_fields() method, TMS needs dAdt !!!
        # loop over the surfaces
        for s in surfaceIDs:

            # Post-processing: estimate the E field
            # step1: get the E field in all tetrahedra
            # v: (number_of_nodes, 1)
            fields = np.einsum('ijk,ij->ik', self.surface[s].gradient, - (v * 1e3)[self.surface[s].node_index_list])

            if dataType == 0:
                # calculate the magnitutde of E
                fields = np.linalg.norm(fields, axis=1, keepdims=True)

            m.append(self.surface[s].sF @ fields)

        logger.debug('Estimate the E field on the surfaces: {0:.2f}s'.format(
            time() - start))

        return m

    def solve_and_normalize(self, node_idx, I, v_norm, I_norm):
        # set RHS (in fem.py, check for speed)
        b = self.fem.assemble_tdcs_neumann_rhs([np.hstack((self.array_layout_node_idx))], [np.hstack((I))], input_type='node')
        # solve
        v = self.fem.solve(b)
        v_elec = [v[self.array_layout_node_idx[k] - 1] for k in range(len(I))]

        # TODO: add postprocessing step from Fang to get e-fields

        v_mean = [np.mean(v_elec[k]) for k in range(len(I))]
        for k in range(len(I)):
            vn = (v_elec[k] - np.mean(v_elec[k])) / np.std(v_mean)
            In = (I[k] - np.mean(I[k])) / np.mean(I[k])

            v_norm[k] = np.append(v_norm[k], vn.reshape(1, len(vn)), axis=0)
            I_norm[k] = np.append(I_norm[k], In.reshape(1, len(In)), axis=0)
        return v, v_elec, v_norm, I_norm

    def run(self, cpus=1, allow_multiple_runs=False, save_mat=True, return_n_max=1):
        """
        Runs the optimization problem

        Parameters
        -------------
        fn_out_mesh: str
            If set, will write out the electric field and currents to the mesh

        fn_out_mesh: str
            If set, will write out the currents and electrode names to a CSV file


        Returns
        ------------
        currents: N_elec x 1 ndarray
            Optimized currents. The first value is the current in the reference electrode
        """
        pass


def fit_ellipsoid(points):
    """
    Determines best fitting ellipsoid to point cloud.

    Parameters
    ----------
    points : np.array of float [n_points x 3]
        Scattered points an ellipsoid is fitted to
    """
    # get convex hull
    hullV = ConvexHull(points)
    lH = len(hullV.vertices)
    hull = np.zeros((lH, 3))
    for i in range(len(hullV.vertices)):
        hull[i] = points[hullV.vertices[i]]
    hull = np.transpose(hull)

    # fit ellipsoid on convex hull
    eansa = ls_ellipsoid(hull[0], hull[1], hull[2])  # get ellipsoid polynomial coefficients
    center, radii, rotmat = polyToParams3D(eansa)    # get ellipsoid 3D parameters

    return center, radii, rotmat


def ls_ellipsoid(xx, yy, zz):
    """
    Finds best fit ellipsoid. (http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html)
    least squares fit to a 3D-ellipsoid:

    Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz  = 1

    Note that sometimes it is expressed as a solution to
    Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz  = 1
    where the last six terms have a factor of 2 in them
    This is in anticipation of forming a matrix with the polynomial coefficients.
    Those terms with factors of 2 are all off diagonal elements.  These contribute
    two terms when multiplied out (symmetric) so would need to be divided by two
    """
    x = xx[:, np.newaxis]
    y = yy[:, np.newaxis]
    z = zz[:, np.newaxis]

    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz = 1
    J = np.hstack((x * x, y * y, z * z, x * y, x * z, y * z, x, y, z))
    K = np.ones_like(x)
    JT = J.transpose()
    JTJ = np.dot(JT, J)
    InvJTJ = np.linalg.inv(JTJ)
    ABC = np.dot(InvJTJ, np.dot(JT, K))

    # Rearrange, move the 1 to the other side
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz - 1 = 0
    #    or
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz + J = 0
    #  where J = -1
    eansa = np.append(ABC, -1)

    return (eansa)


def polyToParams3D(vec):
    """
    Gets 3D parameters of an ellipsoid. (http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html)
    Convert the polynomial form of a 3D-ellipsoid to parameters center, axes, and transformation matrix.
    Algebraic form: X.T * Amat * X --> polynomial form

    Parameters
    ----------
    vec : np.array of float []
        Vector whose elements are the polynomial coefficients A...J

    Returns
    -------
    center : np.array of float [3]
        Center of ellipsoid
    radii : np.array of float [3]
        Axes length
    rotmat : np.array of float [3 x 3]
        Rotation matrix of ellipsoid (direction of principal axes)
    """

    # polynomial
    Amat = np.array(
        [
            [vec[0], vec[3] / 2.0, vec[4] / 2.0, vec[6] / 2.0],
            [vec[3] / 2.0, vec[1], vec[5] / 2.0, vec[7] / 2.0],
            [vec[4] / 2.0, vec[5] / 2.0, vec[2], vec[8] / 2.0],
            [vec[6] / 2.0, vec[7] / 2.0, vec[8] / 2.0, vec[9]]
        ])

    # Algebraic form of polynomial
    # See B.Bartoni, Preprint SMU-HEP-10-14 Multi-dimensional Ellipsoidal Fitting
    # eq. (20) for the following method for finding the center
    A3 = Amat[0:3, 0:3]
    A3inv = inv(A3)
    ofs = vec[6:9] / 2.0

    # center of ellipsoid
    center = -np.dot(A3inv, ofs)

    # Center the ellipsoid at the origin
    Tofs = np.eye(4)
    Tofs[3, 0:3] = center
    R = np.dot(Tofs, np.dot(Amat, Tofs.T))

    # Algebraic form translated to center
    R3 = R[0:3, 0:3]
    s1 = -R[3, 3]
    R3S = R3 / s1
    (el, ec) = eig(R3S)
    recip = 1.0 / np.abs(el)

    # radii
    radii = np.sqrt(recip)

    # rotation matrix
    rotmat = inv(ec)

    return (center, radii, rotmat)


def ellipsoid2cartesian(coords, center, radii, rotmat, return_normal=False):
    """
    Transforms spherical coordinates on ellipsoid to cartesian coordinates and optionally returns surface normal.

    Parameters
    ----------
    coords : np.array of float [n_points x 2]
        Spherical coordinates [theta, phi] in radiant.
    return_normal : bool, optional, default: False
        Additionally return surface normal (pointing outwards)

    Returns
    -------
    coords_cart : np.array of float [n_points x 3]
        Cartesian coordinates of given points
    normal : np.array of float [n_points x 3]
        Surface normal of given points (pointing outwards)
    """

    u = coords[:, 1]
    v = coords[:, 0]

    # Cartesian coordinates that correspond to the spherical angles
    x = radii[0] * np.cos(u) * np.sin(v)
    y = radii[1] * np.sin(u) * np.sin(v)
    z = radii[2] * np.ones_like(u) * np.cos(v)

    # x = radii[0] * np.outer(np.cos(u), np.sin(v))
    # y = radii[1] * np.outer(np.sin(u), np.sin(v))
    # z = radii[2] * np.outer(np.ones_like(u), np.cos(v))

    x_flatten = x.flatten()
    y_flatten = y.flatten()
    z_flatten = z.flatten()

    xyz_flatten = np.hstack((x_flatten[:, np.newaxis], y_flatten[:, np.newaxis], z_flatten[:, np.newaxis])).T
    xyz_flatten_rot = (rotmat @ xyz_flatten + center[:, np.newaxis]).T

    if return_normal:
        normal = 2 * np.hstack(((x_flatten / radii[0] ** 2)[:, np.newaxis],
                                (y_flatten / radii[1] ** 2)[:, np.newaxis],
                                (z_flatten / radii[2] ** 2)[:, np.newaxis]))
        normal = normal / np.linalg.norm(normal, axis=1)[:, np.newaxis]
        normal_rot = (rotmat @ normal.T).T

        return xyz_flatten_rot, normal_rot
    else:
        return xyz_flatten_rot


def relabel_internal_air(m, subpath, label_skin=1005, label_new=1099, label_internal_air=501):
    ''' relabels skin in internal air cavities to something else;
        relevant for charm meshes
    '''
    subject_files = SubjectFiles(subpath=subpath)

    # relabel internal skin to some other label
    label_nifti = nib.load(subject_files.labeling)
    label_affine = label_nifti.affine
    label_img = label_nifti.get_fdata().astype(int)
    label_img = label_img == label_internal_air
    label_img = mrph.binary_dilation(label_img, iterations=2)

    m = copy.copy(m)
    ed = mesh_io.ElementData.from_data_grid(m, label_img, label_affine, order=0)
    idx = ed.value * (m.elm.tag1 == label_skin)
    m.elm.tag1[idx] = label_new
    m.elm.tag2[:] = m.elm.tag1

    return m

# def get_element_intersect_line_surface(p, w, points, con, triangle_center=None, triangle_normals=None):
#     """
#     Get a element indices of where a given line (point p, direction w) intersects with a surface
#     given by points (points) and connectivity list (con).
#
#     Parameters
#     ----------
#     p : np.array of float [3]
#         Base point of ray
#     w : np.array of float [3]
#         Direction of ray
#     points : np.array of float [n_points x 3]
#         Surface points
#     con : np.array of int [n_tri x 3]
#         Connectivity list of triangles
#     triangle_center : np.array of float [n_tri x 3]
#         Center of triangles
#     triangle_normals : np.array of float [n_tri x 3]
#         Normals of triangles
#
#     Returns
#     -------
#     ele_idx : np.array of int [n_tri_intersect]
#         Index of triangles where ray intersects surface
#     dist : np.array of float [n_tri_intersect]
#         Distances between source point p and intersection
#     """
#
#     p = np.tile(p, (con.shape[0], 1))
#     w = np.tile(w, (con.shape[0], 1))
#
#     p1 = points[con[:, 0], :]
#     p2 = points[con[:, 1], :]
#     p3 = points[con[:, 2], :]
#
#     if triangle_center is None:
#         triangle_center = 1 / 3. * (p1 + p2 + p3)
#
#     if triangle_normals is None:
#         triangle_normals = np.cross(p2 - p1, p3 - p1)
#
#     pi = p + w * (-np.sum((p - triangle_center) * triangle_normals, axis=1) / np.sum(w * triangle_normals, axis=1))[:, np.newaxis]
#
#     mask = np.ones(con.shape[0]).astype(bool)
#
#     p1p2 = p1 - p2
#     p3p2 = p3 - p2
#     pp2 = pi - p2
#     p1p3 = p1 - p3
#     p2p3 = p2 - p3
#     pp3 = pi - p3
#
#     l1 = np.sum(p1p2 * pp2, axis=1)
#     l2 = np.sum(p3p2 * pp2, axis=1)
#     l3 = np.sum(p1p3 * pp3, axis=1)
#     l4 = np.sum(p2p3 * pp3, axis=1)
#
#     mask *= (l1 >= 0) * (l2 >= 0) * (l3 >= 0) * (l4 >= 0)
#     ele_idx = np.where(mask)[0]
#     dist = np.linalg.norm(pi[mask] - triangle_center[mask], axis=1)
#
#     return ele_idx, dist

    # def flatten_surface(self):
    #     """
    #     Flatten surface
    #     :return:
    #     """
    #     import flatsurf.halfedge as halfedge
    #     import flatsurf.mccartney1999 as mc1999
    #
    #     # create surface from valid skin region points
    #
    #     triangulation = halfedge.HalfEdgeTriangulation.from_coords_and_simplices(self.skin_surface.nodes_valid, self.skin_surface.tr_nodes_valid.tolist())
    #     flat_trig = mc1999.McCartney1999Flattening(triangulation, Et=1e-4, delta=1e-5)
    #     flat_trig.flatten()
    #     coords_2d = [np.concatenate([vertex.coord_2d, [0.0]]) for vertex in triangulation.vertices]
    #
    #     import pynibs
    #     # write hdf5 _geo file
    #     pynibs.write_geo_hdf5_surf(out_fn="/home/kporzig/tmp/test_geo_head.hdf5",
    #                                points=self.skin_surface.nodes_valid,
    #                                con=self.skin_surface.tr_nodes_valid,
    #                                replace=True,
    #                                hdf5_path='/mesh')
    #
    #     pynibs.write_data_hdf5_surf(data=[np.zeros(self.skin_surface.tr_nodes_valid.shape[0])],
    #                                 data_names=["test"],
    #                                 data_hdf_fn_out="/home/kporzig/tmp/test_data_head.hdf5",
    #                                 geo_hdf_fn="/home/kporzig/tmp/test_geo_head.hdf5",
    #                                 replace=True)




# theta = np.pi / 2# np.linspace(0, np.pi, 20)
#         phi = np.pi / 2 # np.linspace(0, 2*np.pi, 20)
#         # coords_sphere = np.hstack((theta[:, np.newaxis], phi[:, np.newaxis]))
#         coords_sphere = np.array([theta, phi])[np.newaxis, :]
#         index, point = self.ellipsoid2subject(coords_sphere)
#
#         import pynibs
#         pynibs.write_geo_hdf5_surf(out_fn="/home/kporzig/tmp/test_geo_skin_surface_valid_project.hdf5",
#                                    points=self.skin_surface.nodes_valid,
#                                    con=self.skin_surface.tr_nodes_valid,
#                                    replace=True,
#                                    hdf5_path='/mesh')
#
#         data = np.zeros(self.skin_surface.tr_nodes_valid.shape[0])
#         data[index] = 1
#
#         pynibs.write_data_hdf5_surf(data=[data],
#                                     data_names=["intersect"],
#                                     data_hdf_fn_out="/home/kporzig/tmp/test_data_skin_surface_valid_project.hdf5",
#                                     geo_hdf_fn="/home/kporzig/tmp/test_geo_skin_surface_valid_project.hdf5",
#                                     replace=True)
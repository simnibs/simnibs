"""
Functions to find optimal TMS coil positions ior a given cortical target. See examples/tms_optimization.py for
some examples on how to use these.

Written by Ole Numssen & Konstantin Weise, 2019.
Adapted by Guilherme Saturnino, 2019
"""
import numpy as np

import copy
import csv
import re
import os
import time
import glob
import functools
import logging
import gc

import numpy as np
import scipy.spatial
import h5py
import nibabel

from . import optimize_tms
from . import optimization_methods
from . import ADMlib
from ..simulation import fem
from ..simulation import cond
from ..simulation.sim_struct import SESSION, TMSLIST, SimuList, save_matlab_sim_struct
from ..mesh_tools import mesh_io, gmsh_view
from ..utils import transformations
from ..utils.simnibs_logger import logger
from ..utils.file_finder import SubjectFiles
from ..utils.matlab_read import try_to_read_matlab_field, remove_None

# TODO: adapt imports

class TMSoptimize():
    """
    Attributes:
    ------------------------------
    fnamehead: str
        file name of mesh
    subpath (optional): str
        path to m2m folder
    pathfem (optional): str
        path where the leadfield should be saved
    fnamecoil: str
        name of coil file to be used
    cond (optional): list
        list of COND structures with conductivity information
    anisotropy_type (optional): property, can be 'scalar', 'vn' or 'mc'
        type of anisotropy for simulation
    target: (3x1) array
        Optimization target
    target_size (optional): float
        Size of target region, in mm. Defualt: 5
    tissues (optional): list of ints
        Tissues where the target is defined. Default: [2]
    centre(optional): (3x1) array or None
        Position in scalp to use as a reference for the search space. By dafault, will
        project the target to the scalp
    pos_ydir(optional): (3x1) array or None
        Reference position for the coil Y axis, with respect to the target (or the pos
        variable, if it is defined). If left empty, will search positions in a 360 degrees radius. Default: None
    distance (optional): float
        Distance from coil to scalp. Default: 4
    didt (optional): float
        Coil dI/dt value. Default: 1e6
    search_radius (optional): float
        Radius of area where to search for coil positions, in mm. Default: 20
    spatial_resolution (optional): float
        Spatial resolution of search area, in mm. Default: 5
    search_angle (optional): float
        Range of angles to use in search, in degrees. Default: 360
    angle_resolution (optional): float
        Resolution to use for the angles, in degrees. Default: 20
    open_in_gmsh(optional): bool
        Wether to open the results in gmsh. Default: False
    solver_options (optional): str
        Options for the FEM solver. Default: CG+AMG
    method (optional): 'direct' or 'ADM'
        Method to be used. Either 'direct' for running full TMS optimizations or
        'ADM' for using the Auxiliary Dipole Method. 'ADM' is only compatible with ".ccd"
        coil format
    scalp_normals_smoothing_steps (optional): float
        Number of iterations for smoothing the scalp normals to control tangential scalp placement of TMS coil
    keep_hdf5: bool
        Keep intermediate _direct_optimization() ouput. Default: False.
    """

    def __init__(self, matlab_struct=None):
        # : Date when the session was initiated
        self.date = time.strftime("%Y-%m-%d %H:%M:%S")
        self.time_str = time.strftime("%Y%m%d-%H%M%S")
        # Input stuff
        self.fnamehead = None
        self.subpath = None
        self.pathfem = 'tms_optimization/'
        self.fname_tensor = None
        self.mesh = None
        # Name of coil file
        self.fnamecoil = None
        # Conductivity stuff
        self.cond = cond.standard_cond()
        self.anisotropy_type = 'scalar'
        self.aniso_maxratio = 10
        self.aniso_maxcond = 2
        self.anisotropic_tissues = [1, 2]
        # If set, they have priority over fname_tensor
        self.anisotropy_vol = None  # 4-d data with anisotropy information
        self.anisotropy_affine = None  # 4x4 affine transformation from the regular grid
        # Optimization stuff
        self.target = None
        self.target_direction = None
        self.tissues = [2]
        self.target_size = 5
        self.centre = []
        self.pos_ydir = []
        self.distance = 4.
        self.didt = 1e6
        self.search_radius = 20
        self.spatial_resolution = 5
        self.search_angle = 360
        self.angle_resolution = 30
        self.scalp_normals_smoothing_steps = 20
        self.open_in_gmsh = True
        self.solver_options = ''
        self.method = 'direct'
        self.keep_hdf5 = False

        self.name = ''  # This is here only for leagacy reasons, it doesnt do anything

        self._log_handlers = []
        if matlab_struct:
            self.read_mat_struct(matlab_struct)

    @property
    def fn_tensor_nifti(self):
        return self.fname_tensor

    @fn_tensor_nifti.setter
    def fn_tensor_nifti(self, value):
        self.fname_tensor = value

    @property
    def type(self):
        return self.__class__.__name__

    def _set_logger(self):
        SESSION._set_logger(self, summary=False)

    def _finish_logger(self):
        SESSION._finish_logger(self)

    def _get_vol_info(self):
        return SimuList._get_vol_info(self)

    def _prepare(self):
        """Prepares Leadfield for simulations
        relative paths are made absolute,
        empty fields are set to default values,
        check if required fields exist
        """
        sub_files = SubjectFiles(self.fnamehead, self.subpath)
        self.fnamehead = sub_files.fnamehead
        self.subpath = sub_files.subpath

        self.fnamehead = os.path.abspath(os.path.expanduser(self.fnamehead))
        if not os.path.isfile(self.fnamehead):
            raise IOError('Cannot locate head mesh file: %s' % self.fnamehead)

        if not os.path.isdir(self.subpath):
            logger.warning('Cannot locate subjects m2m folder')
            logger.warning('some postprocessing options might fail')
            self.subpath = None

        if not self.fname_tensor:
            self.fname_tensor = sub_files.tensor_file

        logger.info('Head Mesh: {0}'.format(self.fnamehead))
        logger.info('Subject Path: {0}'.format(self.subpath))
        self.pathfem = os.path.abspath(os.path.expanduser(self.pathfem))
        logger.info('Simulation Folder: {0}'.format(self.pathfem))

        self.mesh = mesh_io.read_msh(self.fnamehead)
        TMSLIST.resolve_fnamecoil(self)
        if self.target is None or len(self.target) != 3:
            raise ValueError('Target for optimization not defined')

        assert self.search_radius > 0
        assert self.spatial_resolution > 0
        assert self.angle_resolution > 0
        # fix a few variables
        if len(self.pos_ydir) == 0:
            self.pos_ydir = None
        if len(self.centre) == 0:
            self.centre = np.copy(self.target)
        if self.target_direction is not None and len(self.target_direction) == 0:
            self.target_direction = None

    def sim_struct2mat(self):
        mat = SimuList.cond_mat_struct(self)
        mat['type'] = self.type
        mat['date'] = remove_None(self.date)
        mat['subpath'] = remove_None(self.subpath)
        mat['fnamehead'] = remove_None(self.fnamehead)
        mat['pathfem'] = remove_None(self.pathfem)
        mat['fname_tensor'] = remove_None(self.fname_tensor)
        mat['tissues'] = remove_None(self.tissues)
        mat['fnamecoil'] = remove_None(self.fnamecoil)

        mat['target'] = remove_None(self.target)
        mat['target_direction'] = remove_None(self.target_direction)
        mat['target_size'] = remove_None(self.target_size)
        mat['centre'] = remove_None(self.centre)
        mat['pos_ydir'] = remove_None(self.pos_ydir)
        mat['distance'] = remove_None(self.distance)
        mat['didt'] = remove_None(self.didt)
        mat['search_radius'] = remove_None(self.search_radius)
        mat['spatial_resolution'] = remove_None(self.spatial_resolution)
        mat['search_angle'] = remove_None(self.search_angle)
        mat['angle_resolution'] = remove_None(self.angle_resolution)
        mat['open_in_gmsh'] = remove_None(self.open_in_gmsh)
        mat['solver_options'] = remove_None(self.solver_options)
        mat['method'] = remove_None(self.method)
        mat['scalp_normals_smoothing_steps'] = remove_None(self.scalp_normals_smoothing_steps)
        return mat

    @classmethod
    def read_mat_struct(self, mat):
        """ Reads form matlab structure
        Parameters
        ------------------
        mat: scipy.io.loadmat
            Loaded matlab structure
        """
        self = self()
        SimuList.read_cond_mat_struct(self, mat)
        self.date = try_to_read_matlab_field(
            mat, 'date', str, self.date
        )
        self.subpath = try_to_read_matlab_field(
            mat, 'subpath', str, self.subpath
        )
        self.fnamehead = try_to_read_matlab_field(
            mat, 'fnamehead', str, self.fnamehead
        )
        self.pathfem = try_to_read_matlab_field(
            mat, 'pathfem', str, self.pathfem
        )
        self.fnamecoil = try_to_read_matlab_field(
            mat, 'fnamecoil', str, self.fnamecoil
        )
        self.fname_tensor = try_to_read_matlab_field(
            mat, 'fname_tensor', str, self.fname_tensor
        )
        self.target = try_to_read_matlab_field(
            mat, 'target', list, self.target
        )
        self.target_direction = try_to_read_matlab_field(
            mat, 'target_direction', list, self.target_direction
        )
        self.target_size = try_to_read_matlab_field(
            mat, 'target_size', float, self.target_size
        )
        self.centre = try_to_read_matlab_field(
            mat, 'centre', list, self.centre
        )
        self.centre = try_to_read_matlab_field(
            mat, 'center', list, self.centre
        )
        self.pos_ydir = try_to_read_matlab_field(
            mat, 'pos_ydir', list, self.pos_ydir
        )
        self.distance = try_to_read_matlab_field(
            mat, 'distance', float, self.distance
        )
        self.didt = try_to_read_matlab_field(
            mat, 'didt', float, self.didt
        )
        self.search_radius = try_to_read_matlab_field(
            mat, 'search_radius', float, self.search_radius
        )
        self.spatial_resolution = try_to_read_matlab_field(
            mat, 'spatial_resolution', float, self.spatial_resolution
        )
        self.search_angle = try_to_read_matlab_field(
            mat, 'search_angle', float, self.search_angle
        )
        self.angle_resolution = try_to_read_matlab_field(
            mat, 'angle_resolution', float, self.angle_resolution
        )
        self.open_in_gmsh = try_to_read_matlab_field(
            mat, 'open_in_gmsh', bool, self.open_in_gmsh
        )
        self.solver_options = try_to_read_matlab_field(
            mat, 'solver_options', str, self.solver_options
        )
        self.method = try_to_read_matlab_field(
            mat, 'method', str, self.method
        )
        self.scalp_normals_smoothing_steps = try_to_read_matlab_field(
            mat, 'scalp_normals_smoothing_steps', int, self.scalp_normals_smoothing_steps
        )
        return self

    def run(self, cpus=1, allow_multiple_runs=False, save_mat=True, return_n_max=1):
        ''' Runs the tms optimization

        Parameters
        -----------
        cpus: int (optional)
            Number of cpus to use. Not nescessaraly will use all cpus. Default: 1
        allow_multiple_runs: bool (optinal)
            Wether to allow multiple runs in one folder. Default: False
        save_mat: bool (optional)
            Whether to save the ".mat" file of this structure
        return_n_max: int (optional)
            Return n-th best solutions. Default: 1

        Returns
        --------
        best_results: array_like
            optimal coil positions/orientations. shape  = (return_n_max, 4, 4)
        '''
        self._set_logger()
        dir_name = os.path.abspath(os.path.expanduser(self.pathfem))
        if os.path.isdir(dir_name):
            g = glob.glob(os.path.join(dir_name, 'simnibs_simulation*.mat'))
            if len(g) > 0 and not allow_multiple_runs:
                raise IOError(
                    '\nFound already existing simulation results in directory.'
                    '\nPlease run the simulation in a new directory or delete'
                    ' the simnibs_simulation*.mat files from the folder : {0}'.format(dir_name))
            logger.info(
                'Running simulations in the directory: {0}'.format(dir_name))
        else:
            logger.info('Running simulations on new directory: {0}'.dir_name)
            os.makedirs(dir_name)
        assert return_n_max >= 1

        self._prepare()
        if save_mat:
            save_matlab_sim_struct(
                self,
                os.path.join(
                    dir_name,
                    'simnibs_simulation_{0}.mat'.format(self.time_str)))
        logger.info(str(self))
        pos_matrices = self._get_coil_positions()
        target_region = self._get_target_region()
        cond_field = SimuList.cond2elmdata(self)

        if len(target_region) == 0:
            raise ValueError('Did not find any elements within the defined target region')

        logger.info(f'Searching {len(pos_matrices)} positions')
        # write mesh with target and get file with grid
        m = self.mesh
        fn_target = os.path.join(self.pathfem, 'target.msh')
        target = np.zeros(m.elm.nr)
        target[target_region - 1] = 1
        m.add_element_field(target, 'Target')
        m.write(fn_target)
        v = m.view(
            visible_tags=self.tissues,
            visible_fields='all'
        )
        v.View[0].CustomMax = 1
        v.View[0].CustomMin = 0
        m.elmdata = []

        # Run simulations
        if self.method.lower() == 'direct':
            # Write out the grid
            optimize_tms.plot_matsimnibs_list(
                pos_matrices,
                np.ones(len(pos_matrices)),
                "Grid",
                os.path.join(self.pathfem, 'coil_positions.geo')
            )
            v.add_merge(os.path.join(self.pathfem, 'coil_positions.geo'))
            v.add_view(
                CustomMax=1, CustomMin=1,
                VectorType=4, CenterGlyphs=0,
                Visible=1, ColormapNumber=0
            )
            v.write_opt(fn_target)
            if self.open_in_gmsh:
                mesh_io.open_in_gmsh(fn_target, True)
            E_roi = self._direct_optimize(cond_field, target_region, pos_matrices, cpus)
        elif self.method.lower() == 'adm':
            E_roi, pos_matrices = self._ADM_optimize(cond_field, target_region)
        else:
            raise ValueError("method should be 'direct' or 'ADM'")
        # Update the .geo file with the E values
        optimize_tms.plot_matsimnibs_list(
            pos_matrices,
            E_roi,
            "E at target",
            os.path.join(self.pathfem, 'coil_positions.geo')
        )
        v.add_merge(os.path.join(self.pathfem, 'coil_positions.geo'))
        v.add_view(VectorType=4, CenterGlyphs=0)
        v.write_opt(fn_target)

        # Run one extra simulation with the best position
        logger.info('Re-running best position')
        fn_hdf5 = os.path.join(self.pathfem, self._name_hdf5())
        fn_out = fn_hdf5[:-5] + '.msh'
        fn_geo = fn_hdf5[:-5] + '_coil_pos.geo'
        fem.tms_coil(
            self.mesh, cond_field, self.fnamecoil, 'eEjJ',
            [pos_matrices[np.argmax(E_roi)]],
            [self.didt],
            [fn_out],
            [fn_geo],
        )
        # Write the target to the final simulation file
        m = mesh_io.read_msh(fn_out)
        m.add_element_field(target, 'Target')
        m.write(fn_out)
        v = m.view(
            visible_tags=self.tissues,
            visible_fields=['magnE']
        )
        v.View[-1].CustomMax = 1
        v.View[-1].CustomMin = 0
        v.add_merge(fn_geo)
        v.add_merge(os.path.join(self.pathfem, 'coil_positions.geo'))
        v.add_view(VectorType=4, CenterGlyphs=0)
        v.write_opt(fn_out)
        if self.open_in_gmsh:
            mesh_io.open_in_gmsh(fn_out, True)

        logger.info('\n' + self.summary(pos_matrices[np.argmax(E_roi)]))

        # return optimum coil position(s)
        return_n_max = np.min((return_n_max, E_roi.shape[0]))
        best_results = np.array(pos_matrices)[E_roi.argsort()[-return_n_max:][::-1]]
        return best_results

    def _name_hdf5(self):
        try:
            subid = os.path.splitext(os.path.basename(self.fnamehead))[0]
            subid += '_'
        except TypeError:
            subid = ''
        try:
            coil_name = os.path.splitext(os.path.basename(self.fnamecoil))[0]
            if coil_name.endswith('.nii'):
                coil_name = coil_name[:-4] + '_nii'
            coil_name = '_' + coil_name
        except TypeError:
            coil_name = ''
        name = '{0}TMS_optimize{1}.hdf5'.format(subid, coil_name)
        return name

    def _get_coil_positions(self):
        return optimize_tms.get_opt_grid(
            self.mesh, self.centre,
            handle_direction_ref=self.pos_ydir,
            distance=self.distance, radius=self.search_radius,
            resolution_pos=self.spatial_resolution,
            resolution_angle=self.angle_resolution,
            angle_limits=[-self.search_angle / 2, self.search_angle / 2],
            scalp_normals_smoothing_steps=self.scalp_normals_smoothing_steps
        )

    def _get_target_region(self):
        return optimize_tms.define_target_region(
            self.mesh,
            self.target,
            self.target_size,
            self.tissues
        )

    def _direct_optimize(self, cond_field, target_region, pos_matrices, cpus):
        didt_list = [self.didt for i in pos_matrices]
        fn_hdf5 = os.path.join(self.pathfem, self._name_hdf5())
        if os.path.isfile(fn_hdf5):
            os.remove(fn_hdf5)
        dataset = 'tms_optimization/E_magn'
        volumes = self.mesh.elements_volumes_and_areas()[target_region]
        # Define postporcessing to calculate average field magnitude
        if self.target_direction is None:
            direction = None
        else:
            if len(self.target_direction) != 3:
                raise ValueError('target direction should have 3 elements!')
            direction = np.array(self.target_direction, dtype=float)
            direction /= np.linalg.norm(direction)

        def postprocessing(E, target_region, volumes, direction):
            if direction is None:
                return np.average(
                    np.linalg.norm(E[target_region - 1], axis=1),
                    weights=volumes,
                )
            else:
                return np.average(
                    np.sum(E[target_region - 1] * direction[None, :], axis=1),
                    weights=volumes
                )

        postpro = functools.partial(
            postprocessing,
            target_region=target_region,
            volumes=volumes,
            direction=direction
        )
        fem.tms_many_simulations(
            self.mesh, cond_field,
            self.fnamecoil,
            pos_matrices, didt_list,
            fn_hdf5, dataset,
            post_pro=postpro,
            solver_options=self.solver_options,
            n_workers=cpus
        )
        # Read the fields
        with h5py.File(fn_hdf5, 'a') as f:
            E_roi = f[dataset][:]

        if not hasattr(self, 'keep_hdf5') or not self.keep_hdf5:
            os.remove(fn_hdf5)

        return E_roi

    def _ADM_optimize(self, cond_field, target_region):
        coil_matrices, rotations = optimize_tms.get_opt_grid_ADM(
            self.mesh, self.centre,
            handle_direction_ref=self.pos_ydir,
            distance=self.distance, radius=self.search_radius,
            resolution_pos=self.spatial_resolution,
            resolution_angle=self.angle_resolution,
            angle_limits=[-self.search_angle / 2, self.search_angle / 2],
            scalp_normals_smoothing_steps=self.scalp_normals_smoothing_steps
        )
        # trasnform coil matrix to meters
        coil_matrices[:3, 3, :] *= 1e-3

        baricenters = self.mesh.elements_baricenters()

        th = self.mesh.elm.elm_type == 4
        if not self.fnamecoil.endswith('.ccd'):
            raise ValueError('ADM optimization is only possible with ".ccd" coil files')
        if not np.all(th[target_region - 1]):
            raise ValueError('Target region must contain only tetrahedra')
        ccd_file = np.loadtxt(self.fnamecoil, skiprows=2)
        dipoles, moments = ccd_file[:, 0:3], ccd_file[:, 3:]
        # Run dipole simulations
        S = fem.FEMSystem.electric_dipole(
            self.mesh, cond_field,
            solver_options=self.solver_options
        )

        vols = self.mesh.elements_volumes_and_areas()

        def calc_dipole_J(dipole_dir):
            Jp = mesh_io.ElementData(np.zeros((self.mesh.elm.nr, 3), dtype=float))
            Jp[target_region] = dipole_dir
            b = S.assemble_electric_dipole_rhs(Jp)
            v = mesh_io.NodeData(S.solve(b), mesh=self.mesh)
            m = fem.calc_fields(v, 'J', cond=cond_field)
            J = m.field['J'][:] + Jp[:]
            J /= np.sum(vols[target_region])
            return J

        if self.target_direction is None:
            J_x = calc_dipole_J([1, 0, 0]) * vols[:, None]
            J_y = calc_dipole_J([0, 1, 0]) * vols[:, None]
            J_z = calc_dipole_J([0, 0, 1]) * vols[:, None]
            del S
            gc.collect()
            logger.info('Running ADM')
            # Notice that there is an uknown scale factor
            # as we need to know the pulse angular frequency
            # \Omega and amplitude A
            E_roi = ADMlib.ADMmag(
                baricenters[th].T * 1e-3,
                J_x[th].T, J_y[th].T, J_z[th].T,
                dipoles.T, moments.T,  # .ccd file is already in SI units
                coil_matrices, rotations
            ) * self.didt
        else:
            if len(self.target_direction) != 3:
                raise ValueError('target direction should have 3 elements!')
            direction = np.array(self.target_direction, dtype=float)
            direction /= np.linalg.norm(direction)
            J_d = calc_dipole_J(direction) * vols[:, None]
            E_roi = ADMlib.ADM(
                baricenters[th].T * 1e-3,
                J_d[th].T,
                dipoles.T, moments.T,  # .ccd file is already in SI units
                coil_matrices, rotations
            ) * self.didt

        z = np.array([0., 0., 1.])
        pos_matrices = []
        coil_matrices[:3, 3, :] *= 1e3
        for cm in coil_matrices.transpose(2, 0, 1):
            for r in rotations.T:
                R = np.eye(4)
                R[:3, :3] = np.array([np.cross(r, z), r, z]).T
                pos_matrices.append(cm.dot(R))

        return E_roi.T.reshape(-1), pos_matrices

    def __str__(self):
        string = 'Subject Folder: %s\n' % self.subpath
        string += 'Mesh file name: %s\n' % self.fnamehead
        string += 'Coil file: %s\n' % self.fnamecoil
        string += 'Target: %s\n' % self.target
        string += 'Target Direction: %s\n' % self.target_direction
        string += 'Centre position: %s\n' % self.centre
        string += 'Reference y: %s\n' % self.pos_ydir
        string += 'Coil distance: %s\n' % self.distance
        string += 'Search radius: %s\n' % self.search_radius
        string += 'Spatial resolution: %s\n' % self.spatial_resolution
        string += 'Search angle: %s\n' % self.search_angle
        string += 'Angle resolution: %s\n' % self.angle_resolution
        string += 'method: %s' % self.method
        return string

    def summary(self, position):
        ''' Returns a string with the optimal coil position

        Parameters
        ------------
        position: 4x4 position vector

        Returns
        ------------
        summary: str
            Best coil position
        '''
        s = 'Best coil position\n'
        s += '=============================\n'
        s += '%s\n' % position

        return s


def _create_grid(mesh, pos, distance, radius, resolution_pos, scalp_normals_smoothing_steps=20):
    ''' Creates a position grid '''
    # extract ROI
    msh_surf = mesh.crop_mesh(elm_type=2)
    msh_skin = msh_surf.crop_mesh([5, 1005])
    target_skin = msh_skin.find_closest_element(pos)
    elm_center = msh_skin.elements_baricenters()[:]
    elm_mask_roi = np.linalg.norm(elm_center - target_skin, axis=1) < 1.2 * radius
    elm_center_zeromean = (
        elm_center[elm_mask_roi] -
        np.mean(elm_center[elm_mask_roi], axis=0)
    )
    msh_roi = msh_skin.crop_mesh(elements=msh_skin.elm.elm_number[elm_mask_roi])
    # tangential plane of target_skin point
    u, s, vh = np.linalg.svd(elm_center_zeromean)
    vh = vh.transpose()

    # define regular grid and rotate it to head space
    coords_plane = np.array(
        np.meshgrid(
            np.linspace(-radius, radius, int(2 * radius / resolution_pos + 1)),
            np.linspace(-radius, radius, int(2 * radius / resolution_pos + 1)),
        )
    ).T.reshape(-1, 2)
    coords_plane = coords_plane[np.linalg.norm(coords_plane, axis=1) <= radius]
    coords_plane = np.dot(coords_plane, vh[:, :2].transpose()) + target_skin

    # project grid-points to skin surface
    coords_mapped = []
    coords_normals = []
    normals_roi = msh_roi.triangle_normals(smooth=scalp_normals_smoothing_steps)

    q1 = coords_plane + 1e2 * resolution_pos * vh[:, 2]
    q2 = coords_plane - 1e2 * resolution_pos * vh[:, 2]
    if not q1.size and not q2.size:
        raise ValueError(f"Couldn't determine valid coil positions within search radius. Search radius too small?")
    idx, pos = msh_roi.intersect_segment(q1, q2)
    for i, c in enumerate(coords_plane):
        intersections = idx[:, 0] == i
        if np.any(intersections):
            intersect_pos = pos[intersections]
            intersect_triangles = idx[intersections, 1]
            dist = np.linalg.norm(c[None, :] - intersect_pos, axis=1)
            closest = np.argmin(dist)
            coords_normals.append(normals_roi[intersect_triangles[closest]])
            coords_mapped.append(intersect_pos[closest])

    coords_mapped = np.array(coords_mapped)
    coords_normals = np.array(coords_normals)
    inside = np.linalg.norm(coords_mapped - target_skin, axis=1) <= radius
    coords_mapped += distance * coords_normals
    return coords_mapped[inside], coords_normals[inside]

def _rotate_system(R, angle_limits, angle_res):
    ''' Rotates the vector "y" aroud "z" between the given limits and in the given
    resolution and return rotation matrices'''
    # Define rotation matrix around Z
    n_steps = int((angle_limits[1] - angle_limits[0])/angle_res + 1)
    angles = np.linspace(angle_limits[0], angle_limits[1], n_steps)
    angles = np.deg2rad(angles[(angles > -180.1) * (angles < 180.)])
    matrices = []
    for a in angles:
        Rz = np.array((
            (np.cos(a), -np.sin(a), 0),
            (np.sin(a), np.cos(a), 0),
            (0, 0, 1),
        ))
        matrices.append(R.dot(Rz))
    return matrices


def get_opt_grid(mesh, pos, handle_direction_ref=None, distance=1., radius=20,
                 resolution_pos=1, resolution_angle=20, angle_limits=None, scalp_normals_smoothing_steps=20):
    """ Determine the coil positions and orientations for bruteforce TMS optimization

    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh object
        Simnibs mesh object
    pos: ndarray
        Coordinates (x, y, z) of reference position
    handle_direction_ref (optinal): list of float or np.ndarray
        Vector of handle prolongation direction, in relation to "pos". (Default: do
        not select a handle direction and scan rotations from -180 to 180)
    distance: float or None
        Coil distance to skin surface [mm]. (Default: 1.)
    radius: float or None
        Radius of region of interest around the reference position, where the
        bruteforce simulations are conducted
    resolution_pos: float or None
        Resolution in mm of the coil positions in the region of interest.
    resolution_angle: float or None
        Resolution in deg of the coil positions in the region of interest (Default: 20)
    angle_limits: list of float or None
        Range of angles to get coil rotations for (Default: [-180, 180])
    scalp_normals_smoothing_steps: int
        number of smoothing steps for scalp surface (Default: 20)

    Returns
    -------
    matsimnibs_list: list
        list of MATSIMNIBS matrices
    """
    # creates the spatial grid
    coords_mapped, coords_normals = _create_grid(
        mesh, pos, distance, radius, resolution_pos, scalp_normals_smoothing_steps)
    
    # Determines the seed y direction
    if handle_direction_ref is None:
        y_seed = np.array([0., 1., 0.])
    else:
        y_seed = np.array(handle_direction_ref) - np.array(pos)
        if np.isclose(np.linalg.norm(y_seed), 0.):
            raise ValueError('The coil Y axis reference is too close to the coil center! ')
    if angle_limits is None:
        angle_limits = -180, 180

    matrices = []
    for p, n in zip(coords_mapped, coords_normals):
        z = -n
        y = y_seed - (z * y_seed.dot(z))
        y /= np.linalg.norm(y)
        x = np.cross(y, z)
        R = np.array([x, y, z]).T
        rotated = _rotate_system(R, angle_limits, resolution_angle)
        for r in rotated:
            matrices.append(
                np.vstack((
                    np.hstack((r, p[:, None])),
                    [0, 0, 0, 1]
                ))
            )

    return matrices


def plot_matsimnibs_list(matsimnibs_list, values, field_name, fn_geo):
    ''' Plots the center and the y vector of each matsimnibs matrix as a geo file

    Parameters
    -------------
    matsimnibs_list: list
        list of matsimnibs matrices
    values: array
        Value to assign to each matsimnibs
    field_name: str
        Name of field being printed
    fn_geo: str
        Name of output geo file
    '''
    with open(fn_geo, 'w') as f:
        f.write('View"' + field_name + '"{\n')
        for mat, v in zip(matsimnibs_list, values):
            c = mat[:3, 3]
            y = mat[:3, 1] * v
            f.write(
                "VP(" + ", ".join([str(i) for i in c]) + ")"
                "{" + ", ".join([str(i) for i in y]) + "};\n")
            f.write(
                "SP(" + ", ".join([str(i) for i in c]) + ")"
                "{" + str(v) + "};\n")
        f.write("};\n")


def define_target_region(mesh, target_position, target_radius, tags, elm_type=4):
    ''' Defines a target based on a position, a radius and an element tag

    Paramters
    ------------
    mesh: simnibs.mesh_io.msh
        Mesh
    target_position: array of size (3,)
        Position of target
    target_radius: float
        Size of target
    tags: array of ints
        Tag where the target is located
    elm_type: int
        Type of target element (4 for tetrahedra and 2 for triangles)

    Returns
    -------
    elements: array of size (n,)
        Numbers (1-based) of elements in the tag
    '''
    bar = mesh.elements_baricenters()[:]
    dist = np.linalg.norm(bar - target_position, axis=1)
    elm = mesh.elm.elm_number[
        (dist < target_radius) *
        np.isin(mesh.elm.tag1, tags) *
        np.isin(mesh.elm.elm_type, elm_type)
    ]
    return elm


def get_opt_grid_ADM(mesh, pos, handle_direction_ref=None, distance=1., radius=20,
                 resolution_pos=1, resolution_angle=20, angle_limits=None, scalp_normals_smoothing_steps=20):
    """ Determine the coil positions and orientations for ADM TMS optimization

    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh object
        Simnibs mesh object
    pos: ndarray
        Coordinates (x, y, z) of reference position
    handle_direction_ref (optinal): list of float or np.ndarray
        Vector of handle prolongation direction, in relation to "pos". (Default: do
        not select a handle direction and scan rotations from -180 to 180)
    distance: float or None
        Coil distance to skin surface [mm]. (Default: 1.)
    radius: float or None
        Radius of region of interest around the reference position, where the
        bruteforce simulations are conducted
    resolution_pos: float or None
        Resolution in mm of the coil positions in the region of interest.
    resolution_angle: float or None
        Resolution in deg of the coil positions in the region of interest (Default: 20)
    angle_limits: list of float or None
        Range of angles to get coil rotations for (Default: [-180, 180])
    scalp_normals_smoothing_steps: int
        number of smoothing steps for scalp surface (Default: 20)

    Returns
    -------
    matsimnibs_list: ndarray of size 4x4xNpos
        list of MATSIMNIBS matrices
    coil_dir: ndarray of size 3xNrot
        list of coil directions
    """
    # creates the spatial grid
    coords_mapped, coords_normals = _create_grid(
        mesh, pos, distance, radius, resolution_pos, scalp_normals_smoothing_steps)
    
    # Determines the seed y direction
    if handle_direction_ref is None:
        y_seed = np.array([0., 1., 0.])
    else:
        y_seed = np.array(handle_direction_ref) - np.array(pos)
        if np.isclose(np.linalg.norm(y_seed), 0.):
            raise ValueError('The coil Y axis reference is too close to the coil center! ')
    if angle_limits is None:
        angle_limits = -180, 180

    matrices = []
    for p, n in zip(coords_mapped, coords_normals):
        z = -n
        y = y_seed - (z * y_seed.dot(z))
        y /= np.linalg.norm(y)
        x = np.cross(y, z)
        R = np.array([x, y, z]).T
        A = np.eye(4)
        A[:3, :3] = R
        A[:3, 3] = p
        matrices.append(A)

    matrices = np.array(matrices).transpose(1, 2, 0)

    directions = np.array(
        _rotate_system(np.eye(3), angle_limits, resolution_angle)
    )[:, 1].T

    return matrices, directions

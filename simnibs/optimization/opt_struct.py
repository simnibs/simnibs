'''
    Optimization problem set-up and post-processing in SimNIBS
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.
    Copyright (C) 2023 Guilherme B Saturnino, Konstantin Weise

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''


import copy
import csv
import re
import os
import time
import glob
import h5py
import types
import functools
import datetime
import logging
import gc

import numpy as np
import scipy.spatial
import h5py
import nibabel
import scipy.ndimage.morphology as mrph
from numpy.linalg import eig
from scipy.optimize import direct, Bounds, minimize, differential_evolution, shgo, basinhopping

from . import optimization_methods
from . import ADMlib
from ..simulation import fem
from ..utils import cond_utils
from ..simulation.sim_struct import SESSION, TMSLIST, SimuList, save_matlab_sim_struct, ELECTRODE
from ..simulation.fem import get_dirichlet_node_index_cog
from ..simulation.region_of_interest import RegionOfInterest, RegionOfInterestInitializer
from ..simulation.array_layout import create_tdcs_session_from_array, CircularArray, ElectrodeArrayPair, ElectrodeArrayPairOpt, ElectrodeInitializer
from ..simulation.onlinefem import OnlineFEM, postprocess_e
from ..mesh_tools import mesh_io, gmsh_view, Msh, surface
from ..utils import transformations
from ..utils.simnibs_logger import logger
from ..utils.file_finder import SubjectFiles, Templates
from ..utils.matlab_read import try_to_read_matlab_field, remove_None
from ..utils.mesh_element_properties import ElementTags
from ..utils.transformations import subject2mni_coords, create_new_connectivity_list_point_mask
from ..utils.ellipsoid import Ellipsoid, subject2ellipsoid, ellipsoid2subject
from ..utils.TI_utils import get_maxTI, get_dirTI
from ..utils.measures import AUC, integral_focality, ROC
from simnibs import run_simnibs
from simnibs.optimization import optimize_tms

from simnibs.optimization.tms_flex_optimization import TmsFlexOptimization



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
        self.cond = cond_utils.standard_cond()
        self.anisotropy_type = 'scalar'
        self.aniso_maxratio = 10
        self.aniso_maxcond = 2
        self.anisotropic_tissues = [ElementTags.WM, ElementTags.GM]
        # If set, they have priority over fname_tensor
        self.anisotropy_vol = None  # 4-d data with anisotropy information
        self.anisotropy_affine = None  # 4x4 affine transformation from the regular grid
        # Optimization stuff
        self.target = None
        self.target_direction = None
        self.tissues = [ElementTags.GM]
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
        """ Reads parameters from matlab structure

        Parameters
        ----------
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
            self.mesh, cond_field, self.cond, self.fnamecoil, 'eEjJ',
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
        return_n_max = np.min((return_n_max,E_roi.shape[0]))
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
            angle_limits=[-self.search_angle/2, self.search_angle/2],
            scalp_normals_smoothing_steps=self.scalp_normals_smoothing_steps
        )

    def _get_target_region(self):
        return optimize_tms.define_target_region(
            self.mesh,
            self.target,
            self.target_size,
            self.tissues
        )

    def _direct_optimize(self, cond_field, target_region, pos_matrices, cpus, keep_hdf5=False):
        """
        Parameters:
        -----------
        keep_hdf5 : bool (optional)
            Keep .hdf5 file with target E for all coil positions. Default: False.
        """
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

        if not keep_hdf5:
            os.remove(fn_hdf5)

        return E_roi

    def _ADM_optimize(self, cond_field, target_region):
        coil_matrices, rotations = optimize_tms.get_opt_grid_ADM(
            self.mesh, self.centre,
            handle_direction_ref=self.pos_ydir,
            distance=self.distance, radius=self.search_radius,
            resolution_pos=self.spatial_resolution,
            resolution_angle=self.angle_resolution,
            angle_limits=[-self.search_angle/2, self.search_angle/2],
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
        S = fem.DipoleFEM(self.mesh, cond_field, self.solver_options)
        vols = self.mesh.elements_volumes_and_areas()

        def calc_dipole_J(dipole_dir):
            Jp = mesh_io.ElementData(np.zeros((self.mesh.elm.nr, 3), dtype=float))
            Jp[target_region] = dipole_dir

            dip_pos = baricenters[target_region]
            # `dip_dir` is the desired current density in the target region,
            # therefore; weigh the dipole moments such that the current density
            # is equal to this in all elements
            # (factor 1e-9 converts mm3 to m3)
            b = S.assemble_rhs(
                dip_pos,
                np.atleast_2d(dipole_dir) * 1e-9*np.atleast_1d(vols[target_region])[:, None],
                "partial integration",
            )
            if b.ndim == 2:
                b = b.sum(1)

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
            # Notice that there is an unknown scale factor
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



class TDCSoptimize():
    ''' Defines a tdcs optimization problem

    Parameters
    --------------
    leadfield_hdf: str (optional)
        Name of file with leadfield
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float (optional)
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int (optional)
        Maximum number of active electrodes. Default: no maximum
    name: str (optional)
        Name of optimization problem. Default: optimization
    target: list of TDCStarget objects (optional)
        Targets for the optimization. Default: no target
    avoid: list of TDCSavoid objects
        list of TDCSavoid objects defining regions to avoid


    Attributes
    --------------
    leadfield_hdf: str
        Name of file with leadfield
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int
        Maximum number of active electrodes. Default: no maximum

    ledfield_path: str
        Path to the leadfield in the hdf5 file. Default: '/mesh_leadfield/leadfields/tdcs_leadfield'
    mesh_path: str
        Path to the mesh in the hdf5 file. Default: '/mesh_leadfield/'

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
    def __init__(self, leadfield_hdf=None,
                 max_total_current=2e-3,
                 max_individual_current=1e-3,
                 max_active_electrodes=None,
                 name='optimization/tdcs',
                 target=None,
                 avoid=None,
                 open_in_gmsh=True):
        self.leadfield_hdf = leadfield_hdf
        self.max_total_current = max_total_current
        self.max_individual_current = max_individual_current
        self.max_active_electrodes = max_active_electrodes
        self.leadfield_path = '/mesh_leadfield/leadfields/tdcs_leadfield'
        self.mesh_path = '/mesh_leadfield/'
        self.open_in_gmsh = open_in_gmsh
        self._mesh = None
        self._leadfield = None
        self._field_name = None
        self._field_units = None
        self.name = name
        # I can't put [] in the arguments for weird reasons (it gets the previous value)
        if target is None:
            self.target = []
        else:
            self.target = target
        if avoid is None:
            self.avoid = []
        else:
            self.avoid = avoid


    @property
    def lf_type(self):
        if self.mesh is None or self.leadfield is None:
            return None
        if self.leadfield.shape[1] == self.mesh.nodes.nr:
            return 'node'
        elif self.leadfield.shape[1] == self.mesh.elm.nr:
            return 'element'
        else:
            raise ValueError('Could not find if the leadfield is node- or '
                              'element-based')

    @property
    def leadfield(self):
        ''' Reads the leadfield from the HDF5 file'''
        if self._leadfield is None and self.leadfield_hdf is not None:
            with h5py.File(self.leadfield_hdf, 'r') as f:
                self.leadfield = f[self.leadfield_path][:]

        return self._leadfield

    @leadfield.setter
    def leadfield(self, leadfield):
        if leadfield is not None:
            assert leadfield.ndim == 3, 'leadfield should be 3 dimensional'
            assert leadfield.shape[2] == 3, 'Size of last dimension of leadfield should be 3'
        self._leadfield = leadfield

    @property
    def mesh(self):
        if self._mesh is None and self.leadfield_hdf is not None:
            self.mesh = mesh_io.Msh.read_hdf5(self.leadfield_hdf, self.mesh_path)

        return self._mesh

    @mesh.setter
    def mesh(self, mesh):
        elm_type = np.unique(mesh.elm.elm_type)
        if len(elm_type) > 1:
            raise ValueError('Mesh has both tetrahedra and triangles')
        else:
            self._mesh = mesh

    @property
    def field_name(self):
        if self.leadfield_hdf is not None and self._field_name is None:
            try:
                with h5py.File(self.leadfield_hdf, 'r') as f:
                    self.field_name = f[self.leadfield_path].attrs['field']
            except:
                return 'Field'

        if self._field_name is None:
            return 'Field'
        else:
            return self._field_name


    @field_name.setter
    def field_name(self, field_name):
        self._field_name = field_name

    @property
    def field_units(self):
        if self.leadfield_hdf is not None and self._field_units is None:
            try:
                with h5py.File(self.leadfield_hdf, 'r') as f:
                    self.field_units = f[self.leadfield_path].attrs['units']
            except:
                return 'Au'

        if self._field_units is None:
            return 'Au'
        else:
            return self._field_units

    @field_units.setter
    def field_units(self, field_units):
        self._field_units = field_units


    def to_mat(self):
        """ Makes a dictionary for saving a matlab structure with scipy.io.savemat()

        Returns
        --------------------
        dict
            Dictionaty for usage with scipy.io.savemat
        """
        mat = {}
        mat['type'] = 'TDCSoptimize'
        mat['leadfield_hdf'] = remove_None(self.leadfield_hdf)
        mat['max_total_current'] = remove_None(self.max_total_current)
        mat['max_individual_current'] = remove_None(self.max_individual_current)
        mat['max_active_electrodes'] = remove_None(self.max_active_electrodes)
        mat['open_in_gmsh'] = remove_None(self.open_in_gmsh)
        mat['name'] = remove_None(self.name)
        mat['target'] = _save_TDCStarget_mat(self.target)
        mat['avoid'] = _save_TDCStarget_mat(self.avoid)
        return mat

    @classmethod
    def read_mat_struct(cls, mat):
        '''Reads a .mat structure

        Parameters
        -----------
        mat: dict
            Dictionary from scipy.io.loadmat

        Returns
        ----------
        p: TDCSoptimize
            TDCSoptimize structure
        '''
        t = cls()
        leadfield_hdf = try_to_read_matlab_field(
            mat, 'leadfield_hdf', str, t.leadfield_hdf)
        max_total_current = try_to_read_matlab_field(
            mat, 'max_total_current', float, t.max_total_current)
        max_individual_current = try_to_read_matlab_field(
            mat, 'max_individual_current', float, t.max_individual_current)
        max_active_electrodes = try_to_read_matlab_field(
            mat, 'max_active_electrodes', int, t.max_active_electrodes)
        open_in_gmsh = try_to_read_matlab_field(
            mat, 'open_in_gmsh', bool, t.open_in_gmsh)
        name = try_to_read_matlab_field(
            mat, 'name', str, t.name)
        target = []
        if len(mat['target']) > 0:
            for t in mat['target'][0]:
                target_struct = TDCStarget.read_mat_struct(t)
                if target_struct is not None:
                    target.append(target_struct)
        if len(target) == 0:
            target = None

        avoid = []
        if len(mat['avoid']) > 0:
            avoid = []
            for t in mat['avoid'][0]:
                avoid_struct = TDCSavoid.read_mat_struct(t)
                if avoid_struct is not None:
                    avoid.append(avoid_struct)
        if len(avoid) == 0:
            avoid = None

        return cls(leadfield_hdf, max_total_current,
                   max_individual_current, max_active_electrodes,
                   name, target, avoid, open_in_gmsh)

    def get_weights(self):
        ''' Calculates the volumes or areas of the mesh associated with the leadfield
        '''
        assert self.mesh is not None, 'Mesh not defined'
        if self.lf_type == 'node':
            weights = self.mesh.nodes_volumes_or_areas().value
        elif self.lf_type == 'element':
            weights = self.mesh.elements_volumes_and_areas().value
        else:
            raise ValueError('Cant calculate weights: mesh or leadfield not set')

        weights *= self._get_avoid_field()
        return weights

    def _get_avoid_field(self):
        fields = []
        for a in self.avoid:
            a.mesh = self.mesh
            a.lf_type = self.lf_type
            fields.append(a.avoid_field())

        if len(fields) > 0:
            total_field = np.ones_like(fields[0])
            for f in fields:
                total_field *= f
            return total_field

        else:
            return 1.


    def add_target(self, target=None):
        ''' Adds a target to the current tDCS optimization

        Parameters:
        ------------
        target: TDCStarget (optional)
            TDCStarget structure to be added. Default: empty TDCStarget

        Returns:
        -----------
        target: TDCStarget
            TDCStarget added to the structure
        '''
        if target is None:
            target = TDCStarget(mesh=self.mesh, lf_type=self.lf_type)
        self.target.append(target)
        return target


    def add_avoid(self, avoid=None):
        ''' Adds an avoid structure to the current tDCS optimization

        Parameters:
        ------------
        target: TDCStarget (optional)
            TDCStarget structure to be added. Default: empty TDCStarget

        Returns:
        -----------
        target: TDCStarget
            TDCStarget added to the structure
        '''
        if avoid is None:
            avoid = TDCSavoid(mesh=self.mesh, lf_type=self.lf_type)
        self.avoid.append(avoid)
        return avoid


    def _assign_mesh_lf_type_to_target(self):
        for t in self.target:
            if t.mesh is None: t.mesh = self.mesh
            if t.lf_type is None: t.lf_type = self.lf_type
        for a in self.avoid:
            if a.mesh is None: a.mesh = self.mesh
            if a.lf_type is None: a.lf_type = self.lf_type

    def optimize(self, fn_out_mesh=None, fn_out_csv=None):
        ''' Runs the optimization problem

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
        '''
        assert len(self.target) > 0, 'No target defined'
        assert self.leadfield is not None, 'Leadfield not defined'
        assert self.mesh is not None, 'Mesh not defined'
        if self.max_active_electrodes is not None:
            assert self.max_active_electrodes > 1, \
                    'The maximum number of active electrodes should be at least 2'

        if self.max_total_current is None:
            logger.warning('Maximum total current not set!')
            max_total_current = 1e3

        else:
            assert self.max_total_current > 0
            max_total_current = self.max_total_current

        if self.max_individual_current is None:
            max_individual_current = max_total_current

        else:
            assert self.max_individual_current > 0
            max_individual_current = self.max_individual_current

        self._assign_mesh_lf_type_to_target()
        weights = self.get_weights()
        norm_constrained = [t.directions is None for t in self.target]

        # Angle-constrained optimization
        if any([t.max_angle is not None for t in self.target]):
            if len(self.target) > 1:
                raise ValueError("Can't apply angle constraints with multiple target")
            t = self.target[0]
            max_angle = t.max_angle
            indices, directions = t.get_indexes_and_directions()
            assert max_angle > 0, 'max_angle must be >= 0'
            if self.max_active_electrodes is None:
                opt_problem = optimization_methods.TESLinearAngleConstrained(
                    indices, directions,
                    t.target_mean, max_angle, self.leadfield,
                    max_total_current, max_individual_current,
                    weights=weights, weights_target=t.get_weights()
                )

            else:
                opt_problem = optimization_methods.TESLinearAngleElecConstrained(
                    self.max_active_electrodes, indices, directions,
                    t.target_mean, max_angle, self.leadfield,
                    max_total_current, max_individual_current,
                    weights, weights_target=t.get_weights()
                )

        # Norm-constrained optimization
        elif any(norm_constrained):
            if not all(norm_constrained):
                raise ValueError("Can't mix norm and linear constrained optimization")
            if self.max_active_electrodes is None:
                opt_problem = optimization_methods.TESNormConstrained(
                        self.leadfield, max_total_current,
                        max_individual_current, weights
                )
            else:
                opt_problem = optimization_methods.TESNormElecConstrained(
                        self.max_active_electrodes,
                        self.leadfield, max_total_current,
                        max_individual_current, weights
                )
            for t in self.target:
                if t.intensity < 0:
                    raise ValueError('Intensity must be > 0')
                opt_problem.add_norm_constraint(
                    t.get_indexes_and_directions()[0], t.intensity,
                    t.get_weights()
                )

        # Simple QP-style optimization
        else:
            if self.max_active_electrodes is None:
                opt_problem = optimization_methods.TESLinearConstrained(
                    self.leadfield, max_total_current,
                    max_individual_current, weights)

            else:
                opt_problem = optimization_methods.TESLinearElecConstrained(
                    self.max_active_electrodes, self.leadfield,
                    max_total_current, max_individual_current, weights)

            for t in self.target:
                opt_problem.add_linear_constraint(
                    *t.get_indexes_and_directions(), t.intensity,
                    t.get_weights()
                )

        currents = opt_problem.solve()

        logger.log(25, '\n' + self.summary(currents))

        if fn_out_mesh is not None:
            fn_out_mesh = os.path.abspath(fn_out_mesh)
            m = self.field_mesh(currents)
            m.write(fn_out_mesh)
            v = m.view()
            ## Configure view
            v.Mesh.SurfaceFaces = 0
            v.View[0].Visible = 1
            # Change vector type for target field
            offset=2
            if self.lf_type == 'node':
                offset=3
            for i, t in enumerate(self.target):
                v.View[offset + i].VectorType = 4
                v.View[offset + i].ArrowSizeMax = 60
                v.View[offset + i].Visible = 1
            # Electrode geo file
            el_geo_fn = os.path.splitext(fn_out_mesh)[0] + '_el_currents.geo'
            self.electrode_geo(el_geo_fn, currents)
            v.add_merge(el_geo_fn)
            max_c = np.max(np.abs(currents))
            v.add_view(Visible=1, RangeType=2,
                       ColorTable=gmsh_view._coolwarm_cm(),
                       CustomMax=max_c, CustomMin=-max_c,
                       PointSize=10)
            v.write_opt(fn_out_mesh)
            if self.open_in_gmsh:
                mesh_io.open_in_gmsh(fn_out_mesh, True)


        if fn_out_csv is not None:
            self.write_currents_csv(currents, fn_out_csv)

        return currents

    def field(self, currents):
        ''' Outputs the electric fields caused by the current combination

        Parameters
        -----------
        currents: N_elec x 1 ndarray
            Currents going through each electrode, in A. Usually from the optimize
            method. The sum should be approximately zero

        Returns
        ----------
        E: simnibs.mesh.NodeData or simnibs.mesh.ElementData
            NodeData or ElementData with the field caused by the currents
        '''

        assert np.isclose(np.sum(currents), 0, atol=1e-5), 'Currents should sum to zero'
        E = np.einsum('ijk,i->jk', self.leadfield, currents[1:])

        if self.lf_type == 'node':
            E = mesh_io.NodeData(E, self.field_name, mesh=self.mesh)

        if self.lf_type == 'element':
            E = mesh_io.ElementData(E, self.field_name, mesh=self.mesh)

        return E

    def electrode_geo(self, fn_out, currents=None, mesh_elec=None, elec_tags=None,
                      elec_positions=None):
        ''' Creates a mesh with the electrodes and their currents

        Parameters
        ------------
        currents: N_elec x 1 ndarray (optional)
            Electric current values per electrode. Default: do not print currents
        mesh_elec: simnibs.mesh.Msh (optional)
            Mesh with the electrodes. Default: look for a mesh called mesh_electrodes in
            self.leadfield_hdf
        elec_tags: N_elec x 1 ndarray of ints (optional)
            Tags of the electrodes corresponding to each leadfield column. The first is
            the reference electrode. Default: load at the attribute electrode_tags in the
            leadfield dataset
        elec_positions: N_elec x 3 ndarray of floats (optional)
            Positions of the electrodes in the head. If mesh_elec is not defined, will
            create small sphres at those positions instead.
            Default: load at the attribute electrode_pos in the leadfield dataset

        '''
        # First try to set the electrode visualizations using meshed electrodes
        if mesh_elec is None:
            if self.leadfield_hdf is not None:
                try:
                    mesh_elec = mesh_io.Msh.read_hdf5(self.leadfield_hdf, 'mesh_electrodes')
                except KeyError:
                    pass
            else:
                raise ValueError('Please define a mesh with the electrodes')

        if elec_tags is None and mesh_elec is not None:
            if self.leadfield_hdf is not None:
                with h5py.File(self.leadfield_hdf, 'r') as f:
                    elec_tags = f[self.leadfield_path].attrs['electrode_tags']
            else:
                raise ValueError('Please define the electrode tags')

        # If not, use point electrodes
        if mesh_elec is None and elec_positions is None:
            if self.leadfield_hdf is not None:
                with h5py.File(self.leadfield_hdf, 'r') as f:
                    elec_positions = f[self.leadfield_path].attrs['electrode_pos']
            else:
                raise ValueError('Please define the electrode positions')

        if mesh_elec is not None:
            elec_pos = self._electrode_geo_triangles(fn_out, currents, mesh_elec, elec_tags)
            # elec_pos is used for writing electrode names
        elif elec_positions is not None:
            self._electrode_geo_points(fn_out, currents, elec_positions)
            elec_pos = elec_positions
        else:
            raise ValueError('Neither mesh_elec nor elec_positions defined')
        if self.leadfield_hdf is not None:
            with h5py.File(self.leadfield_hdf, 'r') as f:
                try:
                    elec_names = f[self.leadfield_path].attrs['electrode_names']
                    elec_names = [n.decode() if isinstance(n,bytes) else n for n in elec_names]
                except KeyError:
                    elec_names = None

            if elec_names is not None:
                mesh_io.write_geo_text(
                    elec_pos, elec_names,
                    fn_out, name="electrode_names", mode='ba')

    def _electrode_geo_triangles(self, fn_out, currents, mesh_elec, elec_tags):
        if currents is None:
            currents = np.ones(len(elec_tags))

        assert len(elec_tags) == len(currents), 'Define one current per electrode'

        triangles = []
        values = []
        elec_pos = []
        bar = mesh_elec.elements_baricenters()
        norms = mesh_elec.triangle_normals()
        for t, c in zip(elec_tags, currents):
            triangles.append(mesh_elec.elm[mesh_elec.elm.tag1 == t, :3])
            values.append(c * np.ones(len(triangles[-1])))
            avg_norm = np.average(norms[mesh_elec.elm.tag1 == t], axis=0)
            pos = np.average(bar[mesh_elec.elm.tag1 == t], axis=0)
            pos += avg_norm * 4
            elec_pos.append(pos)

        triangles = np.concatenate(triangles, axis=0)
        values = np.concatenate(values, axis=0)
        elec_pos = np.vstack(elec_pos)
        mesh_io.write_geo_triangles(
            triangles - 1, mesh_elec.nodes.node_coord,
            fn_out, values, 'electrode_currents')

        return elec_pos

    def _electrode_geo_points(self, fn_out, currents, elec_positions):
        if currents is None:
            currents = np.ones(len(elec_positions))

        assert len(elec_positions) == len(currents), 'Define one current per electrode'
        mesh_io.write_geo_spheres(elec_positions, fn_out, currents, "electrode_currents")


    def field_mesh(self, currents):
        ''' Creates showing the targets and the field
        Parameters
        -------------
        currents: N_elec x 1 ndarray
            Currents going through each electrode, in A. Usually from the optimize
            method. The sum should be approximately zero

        Returns
        ---------
        results: simnibs.msh.mesh_io.Msh
            Mesh file
        '''
        target_fields = [t.as_field('target_{0}'.format(i+1)) for i, t in
                         enumerate(self.target)]
        weight_fields = [t.as_field('avoid_{0}'.format(i+1)) for i, t in
                         enumerate(self.avoid)]
        e_field = self.field(currents)
        e_magn_field = e_field.norm()
        if self.lf_type == 'node':
            normals = -self.mesh.nodes_normals()[:]
            e_normal_field = np.sum(e_field[:]*normals, axis=1)
            e_normal_field = mesh_io.NodeData(e_normal_field, 'normal' + e_field.field_name, mesh=self.mesh)
        m = copy.deepcopy(self.mesh)
        if self.lf_type == 'node':
            m.nodedata = [e_magn_field, e_field, e_normal_field] + target_fields + weight_fields
        elif self.lf_type == 'element':
            m.elmdata = [e_magn_field, e_field] + target_fields + weight_fields
        return m


    def write_currents_csv(self, currents, fn_csv, electrode_names=None):
        ''' Writes the currents and the corresponding electrode names to a CSV file

        Parameters
        ------------
        currents: N_elec x 1 ndarray
            Array with electrode currents
        fn_csv: str
            Name of CSV file to write
        electrode_names: list of strings (optional)
            Name of electrodes. Default: will read from the electrode_names attribute in
            the leadfield dataset
        '''
        if electrode_names is None:
            if self.leadfield_hdf is not None:
                with h5py.File(self.leadfield_hdf, 'r') as f:
                    electrode_names = f[self.leadfield_path].attrs['electrode_names']
                    electrode_names = [n.decode() if isinstance(n,bytes) else n for n in electrode_names]
            else:
                raise ValueError('Please define the electrode names')

        assert len(electrode_names) == len(currents)
        with open(fn_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            for n, c in zip(electrode_names, currents):
                writer.writerow([n, c])

    def run(self, cpus=1):
        ''' Interface to use with the run_simnibs function

        Parameters
        ---------------
        cpus: int (optional)
            Does not do anything, it is just here for the common interface with the
            simulation's run function
        '''
        if not self.name:
            if self.leadfield_hdf is not None:
                try:
                    name = re.search(r'(.+)_leadfield_', self.leadfield_hdf).group(1)
                except AttributeError:
                    name = 'optimization'
            else:
                name = 'optimization'
        else:
            name = self.name
        out_folder = os.path.dirname(name)
        os.makedirs(out_folder, exist_ok=True)

        # Set-up logger
        fh = logging.FileHandler(name + '.log', mode='w')
        formatter = logging.Formatter(
            '[ %(name)s - %(asctime)s - %(process)d ]%(levelname)s: %(message)s')
        fh.setFormatter(formatter)
        fh.setLevel(logging.DEBUG)
        logger.addHandler(fh)

        fn_summary = name + '_summary.txt'
        fh_s = logging.FileHandler(fn_summary, mode='w')
        fh_s.setFormatter(logging.Formatter('%(message)s'))
        fh_s.setLevel(25)
        logger.addHandler(fh_s)

        fn_out_mesh = name + '.msh'
        fn_out_csv = name + '.csv'
        logger.info('Optimizing')
        logger.log(25, str(self))
        self.optimize(fn_out_mesh, fn_out_csv)
        logger.log(
            25,
            '\n=====================================\n'
            'SimNIBS finished running optimization\n'
            'Mesh file: {0}\n'
            'CSV file: {1}\n'
            'Summary file: {2}\n'
            '====================================='
            .format(fn_out_mesh, fn_out_csv, fn_summary))

        logger.removeHandler(fh)
        logger.removeHandler(fh_s)

        return fn_out_mesh

    def __str__(self):
        s = 'Optimization set-up\n'
        s += '===========================\n'
        s += 'Leadfield file: {0}\n'.format(self.leadfield_hdf)
        s += 'Max. total current: {0} (A)\n'.format(self.max_total_current)
        s += 'Max. individual current: {0} (A)\n'.format(self.max_individual_current)
        s += 'Max. active electrodes: {0}\n'.format(self.max_active_electrodes)
        s += 'Name: {0}\n'.format(self.name)
        s += '----------------------\n'
        s += 'N targets: {0}\n'.format(len(self.target))
        s += '......................\n'.join(
            ['Target {0}:\n{1}'.format(i+1, str(t)) for i, t in
             enumerate(self.target)])
        s += '----------------------\n'
        s += 'N avoid: {0}\n'.format(len(self.avoid))
        s += '......................\n'.join(
            ['Avoid {0}:\n{1}'.format(i+1, str(t)) for i, t in
             enumerate(self.avoid)])
        return s


    def summary(self, currents):
        ''' Returns a string with a summary of the optimization

        Parameters
        ------------
        field: ElementData or NodeData
            Field of interest

        Returns
        ------------
        summary: str
            Summary of field
        '''
        s = 'Optimization Summary\n'
        s += '=============================\n'
        s += 'Total current: {0:.2e} (A)\n'.format(np.linalg.norm(currents, ord=1)/2)
        s += 'Maximum current: {0:.2e} (A)\n'.format(np.max(np.abs(currents)))
        s += 'Active electrodes: {0}\n'.format(int(np.linalg.norm(currents, ord=0)))
        field = self.field(currents)
        s += 'Field Summary\n'
        s += '----------------------------\n'
        s += 'Peak Value (99.9 percentile): {0:.2f} ({1})\n'.format(
            field.get_percentiles(99.9)[0], self.field_units)
        s += 'Mean field magnitude: {0:.2e} ({1})\n'.format(
            field.mean_field_norm(), self.field_units)
        if np.any(self.mesh.elm.elm_type==4):
            v_units = 'mm3'
        else:
            v_units = 'mm2'
        s += 'Focality: 50%: {0:.2e} 70%: {1:.2e} ({2})\n'.format(
            *field.get_focality(cuttofs=[50, 70], peak_percentile=99.9),
            v_units)
        for i, t in enumerate(self.target):
            s += 'Target {0}\n'.format(i + 1)
            s += '    Intensity specified:{0:.2f} achieved: {1:.2f} ({2})\n'.format(
                t.intensity, t.mean_intensity(field), self.field_units)
            if t.max_angle is not None:
                s += ('    Average angle across target: {0:.1f} '
                      '(max set to {1:.1f}) (degrees)\n'.format(
                      t.mean_angle(field), t.max_angle))
            else:
                s += '    Average angle across target: {0:.1f} (degrees)\n'.format(
                    t.mean_angle(field))

        for i, a in enumerate(self.avoid):
            s += 'Avoid {0}\n'.format(i + 1)
            s += '    Mean field magnitude in region: {0:.2e} ({1})\n'.format(
                a.mean_field_norm_in_region(field), self.field_units)

        return s

class TDCStarget:
    ''' Defines a target for TDCS optimization

    Attributes
    -------------
    positions: Nx3 ndarray
        List of target positions, in x, y, z coordinates and in the subject space. Will find the closes mesh points
    indexes: Nx1 ndarray of ints
        Indexes (1-based) of elements/nodes for optimization. Overwrites positions
    directions: Nx3 ndarray
        List of Electric field directions to be optimized for each mesh point, the string
        'normal' or None (for magnitude optimization), Default: 'normal'
    intensity: float (optional)
        Target intensity of the electric field component in V/m. Default: 0.2
    max_angle: float (optional)
        Maximum angle between electric field and target direction, in degrees. Default:
        No maximum
    radius: float (optional)
        Radius of target. All the elements/nodes within the given radies of the indexes
        will be included.
    tissues: list or None (Optional)
        Tissues included in the target. Either a list of integer with tissue tags or None
        for all tissues. Default: None

    THE ONES BELOW SHOULD NOT BE FILLED BY THE USERS IN NORMAL CIRCUNTANCES:
    mesh: simnibs.msh.mesh_io.Msh (optional)
        Mesh where the target is defined. Set by the TDCSoptimize methods
    lf_type: 'node' or 'element'
        Where the electric field values are defined

    '''
    def __init__(self, positions=None, indexes=None, directions='normal',
                 intensity=0.2, max_angle=None, radius=2, tissues=None,
                 mesh=None, lf_type=None):

        self.lf_type = lf_type
        self.mesh = mesh
        self.radius = radius
        self.tissues = tissues
        self.positions = positions
        self.indexes = indexes
        self.intensity = intensity
        self.max_angle = max_angle
        self.directions = directions

    @property
    def directions(self):
        return self._directions

    @directions.setter
    def directions(self, value):
        if value == 'normal':
            pass
        elif value == 'none':
            value = None
        elif isinstance(value, str):
            raise ValueError(
                'Invalid value for directions: f{directions} '
                'valid arguments are "normal", "none" or an array'
            )
        if value is None and self.max_angle is not None:
            raise ValueError(
                "Can't constrain angle in magnitude optimizations"
            )
        self._directions = value


    @classmethod
    def read_mat_struct(cls, mat):
        '''Reads a .mat structure

        Parameters
        -----------
        mat: dict
            Dictionary from scipy.io.loadmat

        Returns
        ----------
        t: TDCStarget
            TDCStarget structure
        '''
        t = cls()
        positions = try_to_read_matlab_field(mat, 'positions', list, t.positions)
        indexes = try_to_read_matlab_field(mat, 'indexes', list, t.indexes)
        directions = try_to_read_matlab_field(mat, 'directions', list, t.directions)
        try:
            directions[0]
        except IndexError:
            directions = 'normal'
        else:
            if isinstance(directions[0], str):
                directions = ''.join(directions)
            if isinstance(directions[0], bytes):
                directions = ''.join([d.decode() for d in directions])
        intensity = try_to_read_matlab_field(mat, 'intensity', float, t.intensity)
        max_angle = try_to_read_matlab_field(mat, 'max_angle', float, t.max_angle)
        radius = try_to_read_matlab_field(mat, 'radius', float, t.radius)
        tissues = try_to_read_matlab_field(mat, 'tissues', list, t.tissues)
        if positions is not None and len(positions) == 0:
            positions = None
        if indexes is not None and len(indexes) == 0:
            indexes = None
        if tissues is not None and len(tissues) == 0:
            tissues = None

        is_empty = True
        is_empty *= t.positions == positions
        is_empty *= t.indexes == indexes
        is_empty *= t.directions == directions
        is_empty *= t.intensity == intensity
        is_empty *= t.max_angle == max_angle
        is_empty *= t.tissues == tissues
        if is_empty:
            return None

        return cls(positions, indexes, directions, intensity, max_angle, radius, tissues)

    def get_weights(self):
        assert self.lf_type is not None, 'Please set a lf_type'

        if self.lf_type == 'node':
            weights = self.mesh.nodes_volumes_or_areas().value
        elif self.lf_type == 'element':
            weights = self.mesh.elements_volumes_and_areas().value
        else:
            raise ValueError('Invalid lf_type: {0}, should be '
                             '"element" or "node"'.format(self.lf_type))

        return weights

    def get_indexes_and_directions(self):
        ''' Calculates the mesh indexes and directions corresponding to this target
        Returns
        ----------
        indexes: (n,) ndarray of ints
            0-based region indexes

        indexes: (n,3) ndarray of floats
            Target directions
        '''
        indexes, mapping = _find_indexes(self.mesh, self.lf_type,
                                         positions=self.positions,
                                         indexes=self.indexes,
                                         tissues=self.tissues,
                                         radius=self.radius)

        directions = _find_directions(self.mesh, self.lf_type,
                                      self.directions, indexes,
                                      mapping)

        return indexes-1, directions


    def as_field(self, name='target_field'):
        ''' Returns the target as an ElementData or NodeData field

        Parameters
        -----------
        name: str
            Name of the field. Default: 'target_field'
        Returns
        ---------
        target: ElementData or NodeData
            A vector field with a vector pointing in the given direction in the target
        '''
        if (self.positions is None) == (self.indexes is None): # negative XOR operation
            raise ValueError('Please set either positions or indexes')

        assert self.mesh is not None, 'Please set a mesh'

        if self.directions is None:
            nr_comp = 1
        else:
            nr_comp = 3

        if self.lf_type == 'node':
            field = np.zeros((self.mesh.nodes.nr, nr_comp))
            field_type = mesh_io.NodeData
        elif self.lf_type == 'element':
            field = np.zeros((self.mesh.elm.nr, nr_comp))
            field_type = mesh_io.ElementData
        else:
            raise ValueError("lf_type must be 'node' or 'element'."
                             " Got: {0} instead".format(self.lf_type))

        indexes, mapping = _find_indexes(self.mesh, self.lf_type,
                                         positions=self.positions,
                                         indexes=self.indexes,
                                         tissues=self.tissues,
                                         radius=self.radius)

        if self.directions is None:
            field[indexes-1] = self.intensity
        else:
            directions = _find_directions(
                self.mesh, self.lf_type,
                self.directions, indexes,
                mapping
            )
            field[indexes-1] = directions * self.intensity

        return field_type(field, name, mesh=self.mesh)

    def mean_intensity(self, field):
        ''' Calculates the mean intensity of the given field in this target

        Parameters
        -----------
        field: Nx3 NodeData or ElementData
            Electric field

        Returns
        ------------
        intensity: float
            Mean intensity in this target and in the target direction
        '''
        if (self.positions is None) == (self.indexes is None): # negative XOR operation
            raise ValueError('Please set either positions or indexes')

        assert self.mesh is not None, 'Please set a mesh'
        assert field.nr_comp == 3, 'Field must have 3 components'

        indexes, mapping = _find_indexes(self.mesh, self.lf_type,
                                         positions=self.positions,
                                         indexes=self.indexes,
                                         tissues=self.tissues,
                                         radius=self.radius)

        f = field[indexes]
        if self.directions is None:
            components = np.linalg.norm(f, axis=1)

        else:
            directions = _find_directions(self.mesh, self.lf_type,
                                          self.directions, indexes,
                                          mapping)

            components = np.sum(f * directions, axis=1)

        if self.lf_type == 'node':
            weights = self.mesh.nodes_volumes_or_areas()[indexes]
        elif self.lf_type == 'element':
            weights = self.mesh.elements_volumes_and_areas()[indexes]
        else:
            raise ValueError("lf_type must be 'node' or 'element'."
                             " Got: {0} instead".format(self.lf_type))

        return np.average(components, weights=weights)

    def mean_angle(self, field):
        ''' Calculates the mean angle between the field and the target

        Parameters
        -----------
        field: Nx3 NodeData or ElementData
            Electric field

        Returns
        ------------
        angle: float
            Mean angle in this target between the field and the target direction, in
            degrees
        '''
        if (self.positions is None) == (self.indexes is None): # negative XOR operation
            raise ValueError('Please set either positions or indexes')

        assert self.mesh is not None, 'Please set a mesh'
        assert field.nr_comp == 3, 'Field must have 3 components'
        if self.directions is None:
            return np.nan

        indexes, mapping = _find_indexes(self.mesh, self.lf_type,
                                         positions=self.positions,
                                         indexes=self.indexes,
                                         tissues=self.tissues,
                                         radius=self.radius)

        directions = _find_directions(self.mesh, self.lf_type,
                                      self.directions, indexes,
                                      mapping)
        if self.intensity < 0:
            directions *= -1
        f = field[indexes]
        components = np.sum(f * directions, axis=1)
        norm = np.linalg.norm(f, axis=1)
        tangent = np.sqrt(norm ** 2 - components ** 2)
        angles = np.rad2deg(np.arctan2(tangent, components))
        if self.lf_type == 'node':
            weights = self.mesh.nodes_volumes_or_areas()[indexes]
        elif self.lf_type == 'element':
            weights = self.mesh.elements_volumes_and_areas()[indexes]
        else:
            raise ValueError("lf_type must be 'node' or 'element'."
                             " Got: {0} instead".format(self.lf_type))
        weights *= norm
        return np.average(angles, weights=weights)

    def __str__(self):
        s = ('positions: {0}\n'
             'indexes: {1}\n'
             'directions: {2}\n'
             'radius: {3}\n'
             'intensity: {4}\n'
             'max_angle: {5}\n'
             'tissues: {6}\n'
             .format(
                 str(self.positions),
                 str(self.indexes),
                 str(self.directions),
                 self.radius,
                 self.intensity,
                 str(self.max_angle),
                str(self.tissues)))
        return s


class TDCSavoid:
    ''' List of positions to be avoided by optimizer

    Attributes
    -------------
    positions: Nx3 ndarray
        List of positions to be avoided, in x, y, z coordinates and in the subject space.
        Will find the closest mesh points
    indexes: Nx1 ndarray of ints
        Indexes (1-based) of elements/nodes to be avoided. Overwrites positions
    weight : float (optional)
        Weight to give to avoid region. The larger, the more we try to avoid it. Default:
        1e3
    radius: float (optional)
        Radius of region. All the elements/nodes within the given radius of the indexes
        will be included.
    tissues: list or None (Optional)
        Tissues to be included in the region. Either a list of integer with tissue tags or None
        for all tissues. Default: None


    Note
    -------
    If both positions and indexes are set to None, and a tissue is set, it will set the
    given weith to all elements/nodes in the given tissues

    THE ONES BELLOW SHOULD NOT BE FILLED BY THE USERS IN NORMAL CIRCUNTANCES:

    mesh: simnibs.msh.mesh_io.Msh (optional)
        Mesh where the target is defined. Set by the TDCSoptimize methods
    lf_type: 'node' or 'element'
        Where the electric field values are defined

    Warning
    -----------
    Changing positions constructing the class
    can cause unexpected behaviour
    '''
    def __init__(self, positions=None, indexes=None,
                 weight=1e3, radius=2, tissues=None,
                 mesh=None, lf_type=None):
        self.lf_type = lf_type
        self.mesh = mesh
        self.radius = radius
        self.tissues = tissues
        self.positions = positions
        self.indexes = indexes
        self.weight = weight


    @classmethod
    def read_mat_struct(cls, mat):
        '''Reads a .mat structure

        Parameters
        -----------
        mat: dict
            Dictionary from scipy.io.loadmat

        Returns
        ----------
        t: TDCSavoid
            TDCSavoid structure
        '''
        t = cls()
        positions = try_to_read_matlab_field(mat, 'positions', list, t.positions)
        indexes = try_to_read_matlab_field(mat, 'indexes', list, t.indexes)
        weight = try_to_read_matlab_field(mat, 'weight', float, t.weight)
        radius = try_to_read_matlab_field(mat, 'radius', float, t.radius)
        tissues = try_to_read_matlab_field(mat, 'tissues', list, t.tissues)
        if positions is not None and len(positions) == 0:
            positions = None
        if indexes is not None and len(indexes) == 0:
            indexes = None
        if tissues is not None and len(tissues) == 0:
            tissues = None

        is_empty = True
        is_empty *= t.positions == positions
        is_empty *= t.indexes == indexes
        is_empty *= t.weight == weight
        is_empty *= t.tissues == tissues
        if is_empty:
            return None

        return cls(positions, indexes, weight, radius, tissues)


    def _get_avoid_region(self):
        if (self.indexes is not None) or (self.positions is not None):
            indexes, _ = _find_indexes(self.mesh, self.lf_type,
                                       positions=self.positions,
                                       indexes=self.indexes,
                                       tissues=self.tissues,
                                       radius=self.radius)
            return indexes
        elif self.tissues is not None:
            if self.lf_type == 'element':
                return self.mesh.elm.elm_number[
                    np.isin(self.mesh.elm.tag1, self.tissues)]
            elif self.lf_type == 'node':
                return self.mesh.elm.nodes_with_tag(self.tissues)
        else:
            raise ValueError('Please define either indexes/positions or tissues')


    def avoid_field(self):
        ''' Returns a field with self.weight in the target area and
        weight=1 outside the target area

        Returns
        ------------
        w: float, >= 1
            Weight field
        '''
        assert self.mesh is not None, 'Please set a mesh'
        assert self.lf_type is not None, 'Please set a lf_type'
        assert self.weight >= 0, 'Weights must be >= 0'
        if self.lf_type == 'node':
            f = np.ones(self.mesh.nodes.nr)
        elif self.lf_type == 'element':
            f = np.ones(self.mesh.elm.nr)
        else:
            raise ValueError("lf_type must be 'node' or 'element'."
                             " Got: {0} instead".format(self.lf_type))

        indexes = self._get_avoid_region()
        f[indexes - 1] = self.weight
        if len(indexes) == 0:
            raise ValueError('Empty avoid region!')

        return f

    def as_field(self, name='weights'):
        ''' Returns a NodeData or ElementData field with the weights

        Paramets
        ---------
        name: str (optional)
            Name for the field

        Returns
        --------
        f: NodeData or ElementData
            Field with weights
        '''
        w = self.avoid_field()
        if self.lf_type == 'node':
            return mesh_io.NodeData(w, name, mesh=self.mesh)
        elif self.lf_type == 'element':
            return mesh_io.ElementData(w, name, mesh=self.mesh)


    def mean_field_norm_in_region(self, field):
        ''' Calculates the mean field magnitude in the region defined by the avoid structure

        Parameters
        -----------
        field: ElementData or NodeData
            Field for which we calculate the mean magnitude
        '''
        assert self.mesh is not None, 'Please set a mesh'
        assert self.lf_type is not None, 'Please set a lf_type'
        indexes = self._get_avoid_region()
        v = np.linalg.norm(field[indexes], axis=1)
        if self.lf_type == 'node':
            weight = self.mesh.nodes_volumes_or_areas()[indexes]
        elif self.lf_type == 'element':
            weight = self.mesh.elements_volumes_and_areas()[indexes]
        else:
            raise ValueError("lf_type must be 'node' or 'element'."
                             " Got: {0} instead".format(self.lf_type))

        return np.average(v, weights=weight)

    def __str__(self):
        s = ('positions: {0}\n'
             'indexes: {1}\n'
             'radius: {2}\n'
             'weight: {3:.1e}\n'
             'tissues: {4}\n'
             .format(
                 str(self.positions),
                 str(self.indexes),
                 self.radius,
                 self.weight,
                 str(self.tissues)))
        return s

class TDCSDistributedOptimize():
    ''' Defines a tdcs optimization problem with distributed sources

    This function uses the problem setup from

    Ruffini et al. "Optimization of multifocal transcranial current
    stimulation for weighted cortical pattern targeting from realistic modeling of
    electric fields", NeuroImage, 2014

    And the algorithm from

    Saturnino et al. "Accessibility of cortical regions to focal TES:
    Dependence on spatial position, safety, and practical constraints."
    NeuroImage, 2019

    Parameters
    --------------
    leadfield_hdf: str (optional)
        Name of file with leadfield
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float (optional)
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int (optional)
        Maximum number of active electrodes. Default: no maximum
    name: str (optional)
        Name of optimization problem. Default: optimization
    target_image: str or pair (array, affine)
        Image to be "reproduced" via the optimization
    mni_space: bool (optional)
        Wether the image is in MNI space. Default True
    subpath: str (optional)
        Path to the subject "m2m" folder. Needed if mni_space=True
    intensity: float
        Target field intensity
    min_img_value: float >= 0 (optional)
        minimum image (for example t value) to be considered. Corresponds to T_min in
        Ruffini et al. 2014. Default: 0
    open_in_gmsh: bool (optional)
        Whether to open the result in Gmsh after the calculations. Default: False

    Attributes
    --------------
    leadfield_hdf: str
        Name of file with leadfield
    max_total_current: float (optional)
        Maximum current across all electrodes (in Amperes). Default: 2e-3
    max_individual_current: float
        Maximum current for any single electrode (in Amperes). Default: 1e-3
    max_active_electrodes: int
        Maximum number of active electrodes. Default: no maximum
    ledfield_path: str
        Path to the leadfield in the hdf5 file. Default: '/mesh_leadfield/leadfields/tdcs_leadfield'
    mesh_path: str
        Path to the mesh in the hdf5 file. Default: '/mesh_leadfield/'

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

    target_image: str or pair (array, affine)
        Image to be "reproduced" via the optimization

    mni_space: bool (optional)
        Wether the image is in MNI space. Default True

    subpath: str (optional)
        Path to the subject "m2m" folder. Needed if mni_space=True

    intensity: float
        Target field intensity

    min_img_value: float >= 0 (optional)
        minimum image (for example t value) to be considered. Corresponds to T_min in
        Ruffini et al. 2014. Default: 0

    open_in_gmsh: bool (optional)
        Whether to open the result in Gmsh after the calculations. Default: False

    Warning
    -----------
    Changing leadfield_hdf, leadfield_path and mesh_path after constructing the class
    can cause unexpected behaviour
    '''
    def __init__(self, leadfield_hdf=None,
                 max_total_current=2e-3,
                 max_individual_current=1e-3,
                 max_active_electrodes=None,
                 name='optimization/tdcs',
                 target_image=None,
                 mni_space=True,
                 subpath=None,
                 intensity=0.2,
                 min_img_value=0,
                 open_in_gmsh=True):

        self._tdcs_opt_obj = TDCSoptimize(
            leadfield_hdf=leadfield_hdf,
            max_total_current=max_total_current,
            max_individual_current=max_individual_current,
            max_active_electrodes=max_active_electrodes,
            name=name,
            target=[],
            avoid=[],
            open_in_gmsh=open_in_gmsh
        )
        self.max_total_current = max_total_current
        self.max_individual_current = max_individual_current
        self.max_active_electrodes = max_active_electrodes
        self.leadfield_path = '/mesh_leadfield/leadfields/tdcs_leadfield'
        self.mesh_path = '/mesh_leadfield/'
        self.target_image = target_image
        self.mni_space = mni_space
        self.open_in_gmsh = open_in_gmsh
        self.subpath = subpath
        self.name = name

        self.intensity = intensity
        self.min_img_value = min_img_value

        if min_img_value < 0:
            raise ValueError('min_img_value must be > 0')

    @property
    def lf_type(self):
        self._tdcs_opt_obj.mesh = self.mesh
        self._tdcs_opt_obj.leadfield = self.leadfield

        return self._tdcs_opt_obj.lf_type

    @property
    def leadfield_hdf(self):
        return self._tdcs_opt_obj.leadfield_hdf

    @leadfield_hdf.setter
    def leadfield_hdf(self, leadfield_hdf):
        self._tdcs_opt_obj.leadfield_hdf = leadfield_hdf

    @property
    def leadfield_path(self):
        return self._tdcs_opt_obj.leadfield_path

    @leadfield_path.setter
    def leadfield_path(self, leadfield_path):
        self._tdcs_opt_obj.leadfield_path = leadfield_path

    @property
    def mesh_path(self):
        return self._tdcs_opt_obj.mesh_path

    @mesh_path.setter
    def mesh_path(self, mesh_path):
        self._tdcs_opt_obj.mesh_path = mesh_path

    @property
    def name(self):
        return self._tdcs_opt_obj.name

    @name.setter
    def name(self, name):
        self._tdcs_opt_obj.name = name

    @property
    def leadfield(self):
        ''' Reads the leadfield from the HDF5 file'''
        self._tdcs_opt_obj.leadfield_hdf = self.leadfield_hdf
        return self._tdcs_opt_obj.leadfield

    @leadfield.setter
    def leadfield(self, leadfield):
        self._tdcs_opt_obj.leadfield = leadfield

    @property
    def mesh(self):
        self._tdcs_opt_obj.leadfield_hdf = self.leadfield_hdf
        return self._tdcs_opt_obj.mesh

    @mesh.setter
    def mesh(self, mesh):
        self._tdcs_opt_obj.mesh = mesh

    @property
    def field_name(self):
        self._tdcs_opt_obj.leadfield_hdf = self.leadfield_hdf
        return self._tdcs_opt_obj._field_name


    @field_name.setter
    def field_name(self, field_name):
        self._tdcs_opt_obj._field_name = field_name

    @property
    def field_units(self):
        self._tdcs_opt_obj.leadfield_hdf = self.leadfield_hdf
        return self._tdcs_opt_obj._field_units

    def to_mat(self):
        """ Makes a dictionary for saving a matlab structure with scipy.io.savemat()

        Returns
        --------------------
        dict
            Dictionaty for usage with scipy.io.savemat
        """
        mat = {}
        mat['type'] = 'TDCSDistributedOptimize'
        mat['leadfield_hdf'] = remove_None(self.leadfield_hdf)
        mat['max_total_current'] = remove_None(self.max_total_current)
        mat['max_individual_current'] = remove_None(self.max_individual_current)
        mat['max_active_electrodes'] = remove_None(self.max_active_electrodes)
        mat['open_in_gmsh'] = remove_None(self.open_in_gmsh)
        mat['name'] = remove_None(self.name)
        mat['target_image'] = remove_None(self.target_image)
        mat['mni_space'] = remove_None(self.mni_space)
        mat['subpath'] = remove_None(self.subpath)
        mat['intensity'] = remove_None(self.intensity)
        mat['min_img_value'] = remove_None(self.min_img_value)

        return mat

    @classmethod
    def read_mat_struct(cls, mat):
        '''Reads a .mat structure

        Parameters
        -----------
        mat: dict
            Dictionary from scipy.io.loadmat

        Returns
        ----------
        p: TDCSoptimize
            TDCSoptimize structure
        '''
        t = cls()
        leadfield_hdf = try_to_read_matlab_field(
            mat, 'leadfield_hdf', str, t.leadfield_hdf)
        max_total_current = try_to_read_matlab_field(
            mat, 'max_total_current', float, t.max_total_current)
        max_individual_current = try_to_read_matlab_field(
            mat, 'max_individual_current', float, t.max_individual_current)
        max_active_electrodes = try_to_read_matlab_field(
            mat, 'max_active_electrodes', int, t.max_active_electrodes)
        open_in_gmsh = try_to_read_matlab_field(
            mat, 'open_in_gmsh', bool, t.open_in_gmsh)
        name = try_to_read_matlab_field(
            mat, 'name', str, t.name)
        target_image = try_to_read_matlab_field(
            mat, 'target_image', str, t.target_image)
        mni_space = try_to_read_matlab_field(
            mat, 'mni_space', bool, t.mni_space)
        subpath = try_to_read_matlab_field(
            mat, 'subpath', str, t.subpath)
        intensity = try_to_read_matlab_field(
            mat, 'intensity', float, t.intensity)
        min_img_value = try_to_read_matlab_field(
            mat, 'min_img_value', float, t.min_img_value)

        return cls(
            leadfield_hdf=leadfield_hdf,
            max_total_current=max_total_current,
            max_individual_current=max_individual_current,
            max_active_electrodes=max_active_electrodes,
            name=name,
            target_image=target_image,
            mni_space=mni_space,
            subpath=subpath,
            intensity=intensity,
            min_img_value=min_img_value,
            open_in_gmsh=open_in_gmsh
        )

    def _target_distribution(self):
        ''' Gets the y and W fields, by interpolating the target_image

        Based on Eq. 1 from
        Ruffini et al. "Optimization of multifocal transcranial current
        stimulation for weighted cortical pattern targeting from realistic modeling of
        electric fields", NeuroImage, 2014
        '''
        assert self.mesh is not None, 'Please set a mesh'
        assert self.min_img_value >= 0, 'min_img_value must be >= 0'
        assert self.intensity is not None, 'intensity not set'
        # load image
        if isinstance(self.target_image, str):
            img = nibabel.load(self.target_image)
            vol = np.array(img.dataobj)
            affine = img.affine
        else:
            vol, affine = self.target_image
        vol=vol.squeeze() # fix when image is "4D", i.e. NxMxKx1
        if vol.ndim != 3:
            raise ValueError('Target image has to be 3D')
        vol[np.isnan(vol)] = 0.0

        # if in MNI space, tranfrom coordinates
        if self.mni_space:
            if self.subpath is None:
                raise ValueError('subpath not set!')
            nodes_mni = transformations.subject2mni_coords(
                self.mesh.nodes[:], self.subpath
            )
            orig_nodes = np.copy(self.mesh.nodes[:])
            self.mesh.nodes.node_coord = nodes_mni
        # Interpolate
        if self.lf_type == 'node':
            field = mesh_io.NodeData.from_data_grid(self.mesh, vol, affine)
        elif self.lf_type == 'element':
            field = mesh_io.ElementData.from_data_grid(self.mesh, vol, affine)
        field = np.float64(field[:])

        # setting values in eyes to zero
        if np.any(self.mesh.elm.tag1 == ElementTags.EYE_BALLS_TH_SURFACE):
            logger.info('setting target values in eyes to zero')
            if self.lf_type == 'node':
                eye_nodes=np.unique(self.mesh.elm.node_number_list[self.mesh.elm.tag1 == ElementTags.EYE_BALLS_TH_SURFACE,:])
                eye_nodes = eye_nodes[eye_nodes>0]
                field[eye_nodes-1] = 0.0 # node indices in mesh are 1-based
            elif self.lf_type == 'element':
                field[self.mesh.elm.tag1 == ElementTags.EYE_BALLS_TH_SURFACE] = 0.0

        if self.mni_space:
            self.mesh.nodes.node_coord = orig_nodes

        W = np.abs(field)
        W[np.abs(field) < self.min_img_value] = self.min_img_value
        y = field[:].copy()
        y[np.abs(field) < self.min_img_value] = 0
        y *= self.intensity

        if np.all(np.abs(field) < self.min_img_value):
            raise ValueError('Target image values are below min_img_value!')
        return y, W

    def normal_directions(self):
        assert self.mesh is not None, 'Please set a mesh'
        assert self.lf_type is not None, 'Please set a lf_type'

        if 4 in self.mesh.elm.elm_type:
            raise ValueError("Can't define a normal direction for volumetric data!")

        if self.lf_type == 'node':
            normals = self.mesh.nodes_normals()[:]
        elif self.lf_type == 'element':
            normals = self.mesh.triangle_normals()[:]

        return -normals

    def field(self, currents):
        ''' Outputs the electric fields caused by the current combination

        Parameters
        -----------
        currents: N_elec x 1 ndarray
            Currents going through each electrode, in A. Usually from the optimize
            method. The sum should be approximately zero

        Returns
        ----------
        E: simnibs.mesh.NodeData or simnibs.mesh.ElementData
            NodeData or ElementData with the field caused by the currents
        '''
        return self._tdcs_opt_obj.field(currents)

    def field_mesh(self, currents):
        ''' Creates showing the targets and the field
        Parameters
        -------------
        currents: N_elec x 1 ndarray
            Currents going through each electrode, in A. Usually from the optimize
            method. The sum should be approximately zero

        Returns
        ---------
        results: simnibs.msh.mesh_io.Msh
            Mesh file
        '''
        e_field = self.field(currents)
        e_magn_field = e_field.norm()
        normals = self.normal_directions()
        e_normal_field = np.sum(e_field[:]*normals, axis=1)
        target_map, W = self._target_distribution()
        erni = (target_map - W*e_normal_field) ** 2 - target_map ** 2
        erni *= len(target_map)/np.sum(W)

        m = copy.deepcopy(self.mesh)
        if self.lf_type == 'node':
            add_field = m.add_node_field
        elif self.lf_type == 'element':
            add_field = m.add_element_field

        add_field(e_field, e_field.field_name)
        add_field(e_magn_field, e_magn_field.field_name)
        add_field(e_normal_field, 'normal' + e_field.field_name)
        add_field(target_map, 'target_map')
        add_field(erni, 'ERNI')
        return m

    def optimize(self, fn_out_mesh=None, fn_out_csv=None):
        ''' Runs the optimization problem

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
        '''
        assert self.leadfield is not None, 'Leadfield not defined'
        assert self.mesh is not None, 'Mesh not defined'
        if self.max_active_electrodes is not None:
            assert self.max_active_electrodes > 1, \
                'The maximum number of active electrodes should be at least 2'

        if self.max_total_current is None:
            logger.warning('Maximum total current not set!')
            max_total_current = 1e3
        else:
            assert self.max_total_current > 0
            max_total_current = self.max_total_current

        if self.max_individual_current is None:
            max_individual_current = max_total_current

        else:
            assert self.max_individual_current > 0
            max_individual_current = self.max_individual_current

        assert self.min_img_value is not None, 'min_img_value not set'
        assert self.intensity is not None, 'intensity not set'

        y, W = self._target_distribution()
        normals = self.normal_directions()
        weights = np.sqrt(self._tdcs_opt_obj.get_weights())

        if self.max_active_electrodes is None:
            opt_problem = optimization_methods.TESDistributed(
                W[None, :, None] * self.leadfield,
                y[:, None]*normals, weights[:, None]*normals,
                max_total_current,
                max_individual_current
            )
        else:
            opt_problem = optimization_methods.TESDistributedElecConstrained(
                self.max_active_electrodes,
                W[None, :, None] * self.leadfield,
                y[:, None]*normals, weights[:, None]*normals,
                max_total_current,
                max_individual_current
            )

        currents = opt_problem.solve()

        logger.log(25, '\n' + self.summary(currents))

        if fn_out_mesh is not None:
            fn_out_mesh = os.path.abspath(fn_out_mesh)
            m = self.field_mesh(currents)
            m.write(fn_out_mesh)
            v = m.view()
            ## Configure view
            v.Mesh.SurfaceFaces = 0
            v.View[ElementTags.GM].Visible = 1
            # Electrode geo file
            el_geo_fn = os.path.splitext(fn_out_mesh)[0] + '_el_currents.geo'
            self._tdcs_opt_obj.electrode_geo(el_geo_fn, currents)
            v.add_merge(el_geo_fn)
            max_c = np.max(np.abs(currents))
            v.add_view(Visible=1, RangeType=2,
                       ColorTable=gmsh_view._coolwarm_cm(),
                       CustomMax=max_c, CustomMin=-max_c)
            v.write_opt(fn_out_mesh)
            if self.open_in_gmsh:
                mesh_io.open_in_gmsh(fn_out_mesh, True)


        if fn_out_csv is not None:
            self._tdcs_opt_obj.write_currents_csv(currents, fn_out_csv)

        return currents

    def __str__(self):
        s = 'Optimization set-up\n'
        s += '===========================\n'
        s += 'Leadfield file: {0}\n'.format(self.leadfield_hdf)
        s += 'Max. total current: {0} (A)\n'.format(self.max_total_current)
        s += 'Max. individual current: {0} (A)\n'.format(self.max_individual_current)
        s += 'Max. active electrodes: {0}\n'.format(self.max_active_electrodes)
        s += 'Name: {0}\n'.format(self.name)
        s += '----------------------\n'
        s += 'Target image: {0}\n'.format(self.target_image)
        s += 'MNI space: {0}\n'.format(self.mni_space)
        s += 'Min. image value: {0}\n'.format(self.min_img_value)
        s += 'Target intensity: {0}\n'.format(self.intensity)
        return s


    def summary(self, currents):
        ''' Returns a string with a summary of the optimization

        Parameters
        ------------
        field: ElementData or NodeData
            Field of interest

        Returns
        ------------
        summary: str
            Summary of field
        '''
        s = self._tdcs_opt_obj.summary(currents)
        # Calculate erri
        field = self.field(currents)[:]
        normals = self.normal_directions()
        field_normal = np.sum(field * normals, axis=1)
        y, W = self._target_distribution()
        erri =  np.sum((y - field_normal * W)**2 - y**2)
        erri *= len(y) /np.sum(W**2)
        # add Erri to messaga
        s += f'Error Relative to Non Intervention (ERNI): {erri:.2e}\n'
        return s


    def run(self, cpus=1):
        ''' Interface to use with the run_simnibs function

        Parameters
        ---------------
        cpus: int (optional)
            Does not do anything, it is just here for the common interface with the
            simulation's run function
        '''
        return TDCSoptimize.run(self)


def _save_TDCStarget_mat(target):
    target_dt = np.dtype(
        [('type', 'O'),
         ('indexes', 'O'), ('directions', 'O'),
         ('positions', 'O'), ('intensity', 'O'),
         ('max_angle', 'O'), ('radius', 'O'),
         ('tissues', 'O')])

    target_mat = np.empty(len(target), dtype=target_dt)

    for i, t in enumerate(target):
        target_mat[i] = np.array([
            ('TDCStarget',
             remove_None(t.indexes),
             remove_None(t.directions),
             remove_None(t.positions),
             remove_None(t.intensity),
             remove_None(t.max_angle),
             remove_None(t.radius),
             remove_None(t.tissues))],
             dtype=target_dt)

    return target_mat


def _save_TDCSavoid_mat(avoid):
    avoid_dt = np.dtype(
        [('type', 'O'),
         ('indexes', 'O'), ('positions', 'O'),
         ('weight', 'O'),('radius', 'O'),
         ('tissues', 'O')])

    avoid_mat = np.empty(len(avoid), dtype=avoid_dt)

    for i, t in enumerate(avoid):
        avoid_mat[i] = np.array([
            ('TDCSavoid',
             remove_None(t.indexes),
             remove_None(t.positions),
             remove_None(t.weight),
             remove_None(t.radius),
             remove_None(t.tissues))],
             dtype=avoid_dt)

    return avoid_mat



def _find_indexes(mesh, lf_type, indexes=None, positions=None, tissues=None, radius=0.):
    ''' Looks into the mesh to find either
        1. nodes/elements withn a given radius of a set of points (defined as positions)
        and in the specified tissues. The fist step will be to find the closest
        node/element
        2. Specific indexes
    Returns the indices of the nodes/elements in the mesh as well as a mapping saying
    from which of the oridinal points the new points were acquired'''



    if (positions is not None) == (indexes is not None): # negative XOR operation
        raise ValueError('Please define either positions or indexes')

    if indexes is not None:
        indexes = np.atleast_1d(indexes)
        return indexes, np.arange(len(indexes))

    if lf_type == 'node':
        if tissues is not None:
            mesh_indexes = mesh.elm.nodes_with_tag(tissues)
        else:
            mesh_indexes = mesh.nodes.node_number

        mesh_pos = mesh.nodes[mesh_indexes]

    elif lf_type == 'element':
        if tissues is not None:
            mesh_indexes = mesh.elm.elm_number[np.isin(mesh.elm.tag1, tissues)]
        else:
            mesh_indexes = mesh.elm.elm_number

        mesh_pos = mesh.elements_baricenters()[mesh_indexes]

    else:
        raise ValueError('lf_type must be either "node" or "element"')

    assert radius >= 0., 'radius should be >= 0'
    assert len(mesh_pos) > 0, 'Could not find any elements or nodes with given tags'
    kdtree = scipy.spatial.cKDTree(mesh_pos)
    pos_projected, indexes = kdtree.query(positions)
    indexes = np.atleast_1d(indexes)
    if radius > 1e-9:
        in_radius = kdtree.query_ball_point(mesh_pos[indexes], radius)
        original = np.concatenate([(i,)*len(ir) for i, ir in enumerate(in_radius)])
        in_radius, uq_idx = np.unique(np.concatenate(in_radius), return_index=True)
        return mesh_indexes[in_radius], original[uq_idx]
    else:
        return mesh_indexes[indexes],  np.arange(len(indexes))

def _find_directions(mesh, lf_type, directions, indexes, mapping=None):
    if directions is None:
        return None
    if directions == 'normal':
        if 4 in np.unique(mesh.elm.elm_type):
            raise ValueError("Can't define a normal direction for volumetric data!")
        if lf_type == 'node':
            directions = -mesh.nodes_normals()[indexes]
        elif lf_type == 'element':
            directions = -mesh.triangle_normals()[indexes]
        return directions
    else:
        directions = np.atleast_2d(directions)
        if directions.shape[1] != 3:
            raise ValueError(
                "directions must be the string 'normal' or a Nx3 array"
            )
        if mapping is None:
            if len(directions) == len(indexes):
                mapping = np.arange(len(indexes))
            else:
                raise ValueError('Different number of indexes and directions and no '
                                 'mapping defined')
        elif len(directions) == 1:
            mapping = np.zeros(len(indexes), dtype=int)

        directions = directions/np.linalg.norm(directions, axis=1)[:, None]
        return directions[mapping]


class TESoptimize():
    """
    Defines a TES optimization problem using node-wise current sources.

    Parameters
    --------------
    electrode : Electrode Object
        Electrode object containing ElectrodeArray instances
        (see /simulation/array_layout.py for pre-implemented examples)
    init_pos : str or list of str or list of str and np.ndarray of float [3]
        Initial positions of movable Electrode arrays (for each movable array)
    fn_eeg_cap : str, optional, default: 'EEG10-10_UI_Jurak_2007.csv'
        Filename of EEG cap to use for initial position (without path)
        - 'EEG10-10_UI_Jurak_2007.csv'
        - 'easycap_BC_TMS64_X21.csv'
        - 'EEG10-10_Cutini_2011.csv'
        - 'EEG10-10_Neuroelectrics.csv'
        - 'EEG10-20_extended_SPM12.csv'
        - 'EEG10-20_Okamoto_2004.csv'
    roi : list of RegionOfInterest class instances
        Region of interest(s) the field is evaluated in.
    anisotropy_type : str
        Specify type of anisotropy for simulation ('scalar', 'vn' or 'mc')
    weights : np.array of float [n_roi]
        Weights for optimizer for ROI specific goal function weighting
    min_electrode_distance : float, optional, default: None
        Minimum electrode distance to ensure during optimization (in mm).
    constrain_electrode_locations : bool, optional, default: False
        Constrains the possible locations of freely movable electrode arrays. Recommended for TTF optimizations,
        where two pairs of large electrode arrays are optimized. If True, parameter bounds for the optimization
        will be specified restricting the array locations to be frontal, parietal or occipital.
    overlap_factor : float, optional, default: 1
        Factor of overlap of allowed lambda regions to place electrodes. (1 corresponds to neatless regions,
        <1 regions have a gap in between, >1 regions are overlapping)
    plot : bool, optional, default: False
        Plot configurations in output folder for visualization and control
    polish : bool, optional, default: True
        If True (default), then scipy.optimize.minimize with the L-BFGS-B method is used to polish the best
        population member at the end, which can improve the minimization.
    optimize_init_vals : bool, optional, default: True
        If True, find initial values for optimization, guaranteeing a valid solution. If False, initial values
        are the center between bounds.
    e_postproc : str, optional, default: "norm"
        Specifies how the raw e-field in the ROI (Ex, Ey, Ez) is post-processed.
        - "norm": electric field magnitude (default)
        - "normal": determine normal component (required surface normals in dirvec)
        - "tangential": determine tangential component (required surface normals in dirvec)
        - "max_TI": maximum envelope for temporal interference fields
        - "dir_TI": directional sensitive maximum envelope for temporal interference fields
    optimizer : str, optional, default: "differential_evolution"
        Optimization algorithm
    goal : list of str [n_roi], or FunctionType, optional, default: ["mean"]
        Implemented or user provided goal functions:
        - "mean": maximize mean e-field in ROI
        - "max": maximize 99.9 percentile of electric field in ROI
        - "focality": Maximize focality  (goal: sensitivity = specificity = 1)
        - "focality_inv": Maximize inverse focality (goal: sensitivity(ROI) = 1, sensitivity(nonROI) = 1)
        - user provided function taking e-field as an input which is  a list of list of np.ndarrays of float
          [n_channel_stim][n_roi] containing np.array with e-field
    track_focality : bool, optional, default: False
        Tracks focality for each goal function value (requires ROI and non-ROI definition)
    run_final_electrode_simulation : bool, optional, default: True
        Runs final simulation with optimized parameters using real electrode model including remeshing.

    Attributes
    --------------
    electrode : Electrode Object
        Electrode object containing ElectrodeArray instances
        (see /simulation/array_layout.py for pre-implemented examples)
    ellipsoid : Ellipsoid Object
        Best fitting ellipsoid to valid skin reagion (used for coordinate system definition)
    msh_nodes_areas : np.ndarray of float [n_nodes_msh]
        Areas of nodes
    node_idx_msh : np.ndarray of int [n_nodes_skin]
        Indices of skin surface nodes in global msh
    """

    def __init__(self, matlab_struct=None):
        """ Initialized TESoptimize class instance """
        # folders and I/O
        self.date = time.strftime("%Y-%m-%d %H:%M:%S")
        self.time_str = time.strftime("%Y%m%d-%H%M%S")
        self.output_folder = None
        self.plot_folder = None
        self.plot = False
        self.fn_final_sim = []
        self.fn_results_hdf5 = None
        self.logger = None
        self.prepared = False

        # headmodel
        self.fn_eeg_cap = 'EEG10-20_Okamoto_2004.csv'
        self.fn_mesh = None
        self.subpath = None
        self.mesh = None
        self.mesh_relabel = None
        self.mesh_nodes_areas = None
        self.ff_templates = Templates()
        self.ff_subject = None
        self.skin_surface = None
        self.node_idx_msh = None
        self.ellipsoid = Ellipsoid()
        self.fn_electrode_mask = self.ff_templates.mni_volume_upper_head_mask

        # roi
        self.roi = []
        self.n_roi = None

        # electrode
        self.electrode = []
        self.electrode_pos = None
        self.electrode_pos_opt = None
        self.min_electrode_distance = 1e-3
        self.n_channel_stim = None
        self.n_iter_dirichlet_correction = None
        self.n_ele_free = None
        self.init_pos = None
        self.init_pos_subject_coords = None

        # goal function
        self.goal = None
        self.goal_dir = None
        self.e_postproc = None
        self.threshold = None
        self.optimizer = "differential_evolution"
        self.weights = None
        self.track_focality = False
        self.constrain_electrode_locations = False
        self.overlap_factor = None
        self.polish = False
        self.n_test = 0  # number of tries to place the electrodes
        self.n_sim = 0   # number of final simulations carried out (only valid electrode positions)
        self.optimize_init_vals = True
        self.bounds = None
        self.x0 = None

        # track goal fun value (in ROI 0) and focality measures for later analysis
        self.goal_fun_value = None
        self.AUC = None
        self.integral_focality = None

        # set default options for optimizer
        self.optimizer_options = None       # passed by user
        self.optimizer_options_std = {"bounds": self.bounds,
                                      "init_vals": self.x0,
                                      "vol_tol": None,
                                      "len_tol": 1. / 3600000000.,
                                      "f_min_rtol": 1e-12,
                                      "maxiter": 1000,
                                      "polish": True,
                                      "disp": True,
                                      "recombination": 0.7,     # differential evolution
                                      "mutation": (0.01, 0.5),  # differential evolution
                                      "popsize": 13,            # differential evolution
                                      "tol": 0.1,               # differential evolution
                                      "locally_biased": False}

        # FEM
        self.run_final_electrode_simulation = True
        self.dirichlet_node = None
        self.dataType = None
        self.anisotropy_type = "scalar"
        self.solver_options = "pardiso"
        self.ofem = None

        if matlab_struct:
            self.read_mat_struct(matlab_struct)

    def _prepare(self):
        """
        Prepares TESoptimize
        """
        # check variable assignments
        ################################################################################################################
        if self.output_folder is None:
            raise ValueError("Please define TESoptimize.output_folder !")

        if self.subpath is None:
            raise ValueError("Please define TESoptime.subpath: m2m_* folder containing the headmodel.")

        if self.roi is None:
            raise ValueError("Please define TESoptime.roi using the simnibs.RegionOfInterest class.")

        if self.electrode is None:
            raise ValueError("Please define TESoptime.electrode using the simnibs.ElectrodeArrayPair or "
                             "simnibs.CircularArray class.")

        if self.goal is None:
            raise ValueError("Please define type of goal function in TESoptimize.goal"
                             " ('mean', 'focality', 'focality_inv')")

        # setup output folders, logging and IO
        ################################################################################################################
        self.plot_folder = os.path.join(self.output_folder, "plots")
        self.fn_results_hdf5 = os.path.join(self.output_folder, "opt.hdf5")

        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

        if not os.path.exists(self.plot_folder):
            os.makedirs(self.plot_folder)

        # setup logger
        self.logger = setup_logger(os.path.join(self.output_folder,
                                                "simnibs_simulation_" + time.strftime("%Y%m%d-%H%M%S")))
        self.logger.log(20, "Setting up output folders, logging and IO ...")

        # setup headmodel
        ################################################################################################################
        self.logger.log(20, "Setting up headmodel ...")

        # get subject specific filenames
        self.ff_subject = SubjectFiles(subpath=self.subpath)

        # read mesh or store in self
        self.fn_mesh = self.ff_subject.fnamehead
        self.mesh = mesh_io.read_msh(self.fn_mesh)

        # Calculate node areas for whole mesh
        self.mesh_nodes_areas = self.mesh.nodes_areas()

        # relabel internal air
        self.mesh_relabel = relabel_internal_air(m=self.mesh,
                                                 subpath=os.path.split(self.mesh.fn)[0],
                                                 label_skin=1005,
                                                 label_new=1099,
                                                 label_internal_air=501)

        # make final skin surface including some additional distance
        self.skin_surface = surface.Surface(mesh=self.mesh_relabel, labels=1005)
        self.skin_surface = valid_skin_region(skin_surface=self.skin_surface,
                                              fn_electrode_mask=self.fn_electrode_mask,
                                              mesh=self.mesh_relabel,
                                              additional_distance=0)

        # get mapping between skin_surface node indices and global mesh nodes
        self.node_idx_msh = np.where(np.isin(self.mesh.nodes.node_coord, self.skin_surface.nodes).all(axis=1))[0]

        # fit optimal ellipsoid to valid skin points
        self.ellipsoid.fit(points=self.skin_surface.nodes)

        # plot skin surface and ellipsoid
        if self.plot:
            try:
                import pynibs
                pynibs.write_geo_hdf5_surf(out_fn=os.path.join(self.plot_folder, f"skin_surface_geo.hdf5"),
                                           points=self.skin_surface.nodes,
                                           con=self.skin_surface.tr_nodes,
                                           replace=True,
                                           hdf5_path='/mesh')

                pynibs.write_data_hdf5_surf(data=[np.zeros(self.skin_surface.tr_nodes.shape[0])],
                                            data_names=["domain"],
                                            data_hdf_fn_out=os.path.join(self.plot_folder, f"skin_surface_data.hdf5"),
                                            geo_hdf_fn=os.path.join(self.plot_folder, f"skin_surface_geo.hdf5"),
                                            replace=True)
            except:
                pass

            # save fitted ellipsoid
            beta = np.linspace(-np.pi / 2, np.pi / 2, 180)
            lam = np.linspace(0, 2 * np.pi, 360)
            coords_sphere_jac = np.array(np.meshgrid(beta, lam)).T.reshape(-1, 2)
            eli_coords_jac = self.ellipsoid.jacobi2cartesian(coords=coords_sphere_jac, return_normal=False)
            np.savetxt(os.path.join(self.output_folder, "plots", "fitted_ellipsoid.txt"), eli_coords_jac)

        # setup ROI
        ################################################################################################################
        self.logger.log(20, "Setting up ROI ...")

        if type(self.roi) is not list:
            self.roi = [self.roi]

        # initialize ROIs if not done already
        for i in range(len(self.roi)):
            if type(self.roi[i]) == RegionOfInterestInitializer:
                self.roi[i].mesh = self.mesh
                self.roi[i] = self.roi[i].initialize()

        self.n_roi = len(self.roi)

        # setup electrode
        ################################################################################################################
        self.logger.log(20, "Setting up electrodes ...")

        if type(self.electrode) is not list:
            self.electrode = [self.electrode]

        # initialize electrodes if not done already
        for i in range(len(self.electrode)):
            if type(self.electrode[i]) == ElectrodeInitializer:
                self.electrode[i] = self.electrode[i].initialize()

        # number of independent stimulation channels
        self.n_channel_stim = len(self.electrode)

        # initialize lists with number of dirichlet correction iterations for convergence analysis
        self.n_iter_dirichlet_correction = [[] for _ in range(self.n_channel_stim)]

        # list containing the number of freely movable arrays for each channel [i_channel_stim]
        self.n_ele_free = [len(ele.electrode_arrays) for ele in self.electrode]

        # list containing beta, lambda, alpha for each freely movable array and for each stimulation channel
        self.electrode_pos = [[0 for _ in range(n_ele_free)] for n_ele_free in self.n_ele_free]

        for i_channel_stim in range(self.n_channel_stim):
            for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
                if _electrode_array.optimize_alpha:
                    self.electrode_pos[i_channel_stim][i_array] = np.zeros(3)
                else:
                    self.electrode_pos[i_channel_stim][i_array] = np.zeros(2)

        # parameter bounds for optimizer (constrain if desired)
        self.bounds = self.get_bounds(constrain_electrode_locations=self.constrain_electrode_locations,
                                      overlap_factor=self.overlap_factor)

        # determine initial values
        self.x0 = self.get_init_vals(bounds=self.bounds, optimize=self.optimize_init_vals)

        # compile node arrays
        for _electrode in self.electrode:
            _electrode.compile_node_arrays()

        # plot electrodes
        if self.plot:
            for i_channel_stim in range(self.n_channel_stim):
                for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
                    _electrode_array.plot(show=False, fn_plot=os.path.join(
                        self.output_folder, "plots", f"electrode_stim_{i_channel_stim}_array_{i_array}.png"))

        # setup optimization
        ################################################################################################################
        self.logger.log(20, "Setting up optimization algorithm ...")

        # equal ROI weighting if None is provided
        if type(self.weights) is float or type(self.weights) is int:
            self.weights = None

        if self.weights is None:
            self.weights = np.ones(len(self.roi)) / len(self.roi)
        elif type(self.weights) is list:
            self.weights = np.array(self.weights)

        if type(self.goal) is not list:
            self.goal = [self.goal]

        if "focality" in self.goal and len(self.goal) != len(self.roi):
            self.goal = ["focality"] * len(self.roi)

        self.goal_dir = []

        for i_roi in range(len(self.roi)):
            if "normal" in self.e_postproc or "tangential" in self.e_postproc:
                self.goal_dir.append(self.roi[i_roi].triangles_normals)
            else:
                self.goal_dir.append(None)

        if ("focality" in self.goal or "focality_inv" in self.goal) and self.threshold is None:
            raise ValueError("Please define TESoptimze.threshold for focality optimization!")

        if (type(self.threshold) != list and type(self.threshold) != np.ndarray) and self.threshold is not None:
            self.threshold = [self.threshold]

        if type(self.e_postproc) is not list:
            self.e_postproc = [self.e_postproc]

        if len(self.e_postproc) != len(self.roi):
            self.e_postproc = self.e_postproc * len(self.roi)

        if self.track_focality and len(self.roi) != 2:
            raise ValueError("Focality can not be computed and tracked with only one ROI (non-ROI missing).")

        if not ((isinstance(self.goal[0], types.FunctionType) and len(
                self.goal) == 1) or "focality" in self.goal or "focality_inv" in self.goal):
            assert len(self.goal) == len(self.roi), "Please provide a goal function for each ROI."
            assert len(self.weights) == len(self.roi), "Number of weights has to match the number ROIs"

        if "focality" in self.goal and len(self.roi) != 2:
            raise ValueError("For focality optimization please provide ROI and non ROI region (in this order).")

        # track goal fun value (in ROI 0) and focality measures for later analysis
        self.goal_fun_value = [[] for _ in range(self.n_channel_stim)]
        self.AUC = [[] for _ in range(self.n_channel_stim)]
        self.integral_focality = [[] for _ in range(self.n_channel_stim)]

        # direct and shgo optimizer do not take init vals
        if self.optimizer in ["direct", "shgo"]:
            self.optimize_init_vals = False

        # define gpc parameters for current estimator
        if self.electrode[0].current_estimator is not None:
            if self.electrode[0].current_estimator.method == "gpc":
                min_idx = 0
                max_idx = 0

                for i_channel_stim in range(self.n_channel_stim):
                    for _electrode_array in self.electrode[i_channel_stim].electrode_arrays:
                        if _electrode_array.optimize_alpha:
                            max_idx += 3
                        else:
                            max_idx += 2

                    self.electrode[i_channel_stim].current_estimator.set_gpc_parameters(
                        lb=self.bounds.lb[min_idx:max_idx],
                        ub=self.bounds.ub[min_idx:max_idx])

                    min_idx = max_idx

        # set default options for optimizer
        self.optimizer_options_std["bounds"] = self.bounds
        self.optimizer_options_std["init_vals"] = self.x0
        self.optimizer_options_std["vol_tol"] = 1. / 3600000000. * 3 * np.sum(self.n_ele_free)

        # insert user specific options
        if self.optimizer_options is not None:
            for key in self.optimizer_options:
                self.optimizer_options_std[key] = self.optimizer_options[key]

        # setup FEM
        ################################################################################################################
        # set dirichlet node to closest node of center of gravity of head model (indexing starting with 1)
        self.dirichlet_node = get_dirichlet_node_index_cog(mesh=self.mesh, roi=self.roi)

        # always compute e-field components (Ex, Ey, Ez), it will be postprocessed later according to self.e_postproc
        self.dataType = [1] * len(self.roi)

        # prepare FEM
        self.ofem = OnlineFEM(mesh=self.mesh,
                              electrode=self.electrode,
                              method="TES",
                              roi=self.roi,
                              anisotropy_type=self.anisotropy_type,
                              solver_options=self.solver_options,
                              fn_results=self.fn_results_hdf5,
                              useElements=True,
                              dataType=self.dataType,
                              dirichlet_node=self.dirichlet_node)

        # log summary
        ################################################################################################################
        self.logger.log(25, f"=" * 100)
        self.logger.log(25, f"headmodel:                        {self.mesh.fn}")
        self.logger.log(25, f"n_roi:                            {self.n_roi}")
        self.logger.log(25, f"anisotropy type:                  {self.ofem.anisotropy_type}")
        self.logger.log(25, f"n_channel_stim:                   {self.n_channel_stim}")
        self.logger.log(25, f"fn_eeg_cap:                       {self.fn_eeg_cap}")
        self.logger.log(25, f"fn_electrode_mask:                {self.fn_electrode_mask}")
        self.logger.log(25, f"FEM solver options:               {self.ofem.solver_options}")
        self.logger.log(25, f"dirichlet_correction:             {self.electrode[0].dirichlet_correction}")
        self.logger.log(25, f"dirichlet_correction_detailed:    {self.electrode[0].dirichlet_correction_detailed}")
        self.logger.log(25, f"current_outlier_correction:       {self.electrode[0].current_outlier_correction}")
        self.logger.log(25, f"optimizer:                        {self.optimizer}")
        self.logger.log(25, f"goal:                             {self.goal}")
        self.logger.log(25, f"e_postproc:                       {self.e_postproc}")
        self.logger.log(25, f"threshold:                        {self.threshold}")
        self.logger.log(25, f"weights:                          {self.weights}")
        self.logger.log(25, f"output_folder:                    {self.output_folder}")
        self.logger.log(25, f"fn_results_hdf5:                  {self.fn_results_hdf5}")

        if self.optimizer_options is not None:
            for key in self.optimizer_options:
                if key != "bounds":
                    self.logger.log(25, f"{key}:                {self.optimizer_options[key]}")

        for i_channel_stim in range(self.n_channel_stim):
            self.logger.log(25, f"Stimulation: {i_channel_stim} (n_ele_free: {self.n_ele_free[i_channel_stim]})")
            self.logger.log(25, f"---------------------------------------------")

            for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
                self.logger.log(25, f"Electrode array [i_channel_stim][i_array]: [{i_channel_stim}][{i_array}]")
                self.logger.log(25, f"\tn_ele: {_electrode_array.n_ele}")
                # self.logger.log(25, f"\tinit_pos: {self.init_pos[i_channel_stim][i_array]}")
                # self.logger.log(25, f"\tcenter: {_electrode_array.center}")
                # self.logger.log(25, f"\tradius: {_electrode_array.radius}")
                # self.logger.log(25, f"\tlength_x: {_electrode_array.length_x}")
                # self.logger.log(25, f"\tlength_y: {_electrode_array.length_y}")

        self.logger.log(25, f"=" * 100)

        self.prepared = True

    def add_electrode(self, electrode=None):
        """
        Adds an electrode to the current TESoptimize

        Parameters
        ----------
        electrode: ElectrodeArrayPair or CircularArray class instance, optional, default=None
            Electrode structure.

        Returns
        -------
        electrode: ElectrodeInitializer
            ElectrodeInitializer structure added to TESoptimize
        """
        if electrode is None:
            electrode = ElectrodeInitializer()
        self.electrode.append(electrode)
        return electrode

    def add_roi(self, roi=None):
        """
        Adds an ROI to the current TESoptimize

        Parameters
        ----------
        roi: RegionOfInterest object, optional, default=None
            ROI structure.

        Returns
        -------
        roi: RegionOfInterestInitializer
            RegionOfInterestInitializer structure added to TESoptimize
        """
        if roi is None:
            roi = RegionOfInterestInitializer()
        self.roi.append(roi)
        return roi

    def sim_struct2mat(self):
        """ Convert sim_struct to mat """
        pass

        # cond
        # mat = SimuList.cond_mat_struct(self.ofem)
        mat = dict()

        # folders and I/O
        mat['date'] = remove_None(self.date)
        mat['output_folder'] = remove_None(self.output_folder)
        mat['plot_folder'] = remove_None(self.plot_folder)
        mat['plot'] = remove_None(self.plot)
        mat['fn_final_sim'] = remove_None(self.fn_final_sim)
        mat['fn_results_hdf5'] = remove_None(self.fn_results_hdf5)
        mat['prepared'] = remove_None(self.prepared)

        # headmodel
        mat['fn_eeg_cap'] = remove_None(self.fn_eeg_cap)
        mat['fn_mesh'] = remove_None(self.fn_mesh)
        mat['fn_electrode_mask'] = remove_None(self.fn_electrode_mask)

        # roi
        mat['n_roi'] = remove_None(self.n_roi)

        roi_dict = dict()

        for i_roi, _roi in enumerate(self.roi):
            roi_dict[f"roi_{i_roi}"] = dict()
            roi_dict[f"roi_{i_roi}"]["center"] = remove_None(_roi.center)
            roi_dict[f"roi_{i_roi}"]["nodes"] = remove_None(_roi.nodes)
            roi_dict[f"roi_{i_roi}"]["con"] = remove_None(_roi.con)
            roi_dict[f"roi_{i_roi}"]["domains"] = remove_None(_roi.domains)
            roi_dict[f"roi_{i_roi}"]["type"] = remove_None(_roi.type)
            roi_dict[f"roi_{i_roi}"]["roi_sphere_center_mni"] = remove_None(_roi.roi_sphere_center_mni)
            roi_dict[f"roi_{i_roi}"]["roi_sphere_center_subject"] = remove_None(_roi.roi_sphere_center_subject)
            roi_dict[f"roi_{i_roi}"]["roi_sphere_radius"] = remove_None(_roi.roi_sphere_radius)

        mat['roi_dict'] = roi_dict

        # electrode
        mat['electrode_pos'] = remove_None(self.electrode_pos)
        mat['electrode_pos_opt'] = remove_None(self.electrode_pos_opt)
        mat['min_electrode_distance'] = remove_None(self.min_electrode_distance)
        mat['n_channel_stim'] = remove_None(self.n_channel_stim)
        mat['n_iter_dirichlet_correction'] = remove_None(self.n_iter_dirichlet_correction)
        mat['n_ele_free'] = remove_None(self.n_ele_free)
        mat['init_pos'] = remove_None(self.init_pos)
        mat['init_pos_subject_coords'] = remove_None(self.init_pos_subject_coords)

        electrode_dict = dict()

        for i_ele, _electrode in enumerate(self.electrode):
            electrode_dict[f"ele_{i_ele}"] = dict()
            electrode_dict[f"ele_{i_ele}"]["type"] = remove_None(_electrode.__class__.__name__)
            electrode_dict[f"ele_{i_ele}"]["dirichlet_correction"] = remove_None(_electrode.dirichlet_correction)
            electrode_dict[f"ele_{i_ele}"]["dirichlet_correction_detailed"] = remove_None(_electrode.dirichlet_correction_detailed)
            electrode_dict[f"ele_{i_ele}"]["current_estimator_method"] = remove_None(_electrode.current_estimator_method)
            electrode_dict[f"ele_{i_ele}"]["current"] = remove_None(_electrode.current)

            if _electrode.__class__.__name__ == "CircularArray":
                electrode_dict[f"ele_{i_ele}"]["radius_inner"] = remove_None(_electrode.radius_inner)
                electrode_dict[f"ele_{i_ele}"]["radius_inner_bounds"] = remove_None(_electrode.radius_inner_bounds)
                electrode_dict[f"ele_{i_ele}"]["radius_outer"] = remove_None(_electrode.radius_outer)
                electrode_dict[f"ele_{i_ele}"]["radius_outer_bounds"] = remove_None(_electrode.radius_outer_bounds)
                electrode_dict[f"ele_{i_ele}"]["distance"] = remove_None(_electrode.distance)
                electrode_dict[f"ele_{i_ele}"]["distance_bounds"] = remove_None(_electrode.distance_bounds)
                electrode_dict[f"ele_{i_ele}"]["n_outer"] = remove_None(_electrode.n_outer)
                electrode_dict[f"ele_{i_ele}"]["n_outer_bounds"] = remove_None(_electrode.n_outer_bounds)

            elif _electrode.__class__.__name__ == "ElectrodeArrayPair":
                electrode_dict[f"ele_{i_ele}"]["center"] = remove_None(_electrode.center)
                electrode_dict[f"ele_{i_ele}"]["radius"] = remove_None(_electrode.radius)
                electrode_dict[f"ele_{i_ele}"]["radius_bounds"] = remove_None(_electrode.radius_bounds)
                electrode_dict[f"ele_{i_ele}"]["length_x"] = remove_None(_electrode.length_x)
                electrode_dict[f"ele_{i_ele}"]["length_x_bounds"] = remove_None(_electrode.length_x_bounds)
                electrode_dict[f"ele_{i_ele}"]["length_y"] = remove_None(_electrode.length_y)
                electrode_dict[f"ele_{i_ele}"]["length_y_bounds"] = remove_None(_electrode.length_y_bounds)

        mat['electrode_dict'] = electrode_dict

        # goal function
        mat['goal'] = remove_None(self.goal)

        if np.array([True for _t in self.goal_dir if _t is None]).any():
            mat['goal_dir'] = [''] * len(self.roi)
        else:
            mat['goal_dir'] = remove_None(self.goal_dir)

        mat['e_postproc'] = remove_None(self.e_postproc)
        mat['threshold'] = remove_None(self.threshold)
        mat['optimizer'] = remove_None(self.optimizer)
        mat['weights'] = remove_None(self.weights)
        mat['track_focality'] = remove_None(self.track_focality)
        mat['constrain_electrode_locations'] = remove_None(self.constrain_electrode_locations)
        mat['overlap_factor'] = remove_None(self.overlap_factor)
        mat['polish'] = remove_None(self.polish)
        mat['n_test'] = remove_None(self.n_test)
        mat['n_sim'] = remove_None(self.n_sim)
        mat['optimize_init_vals'] = remove_None(self.optimize_init_vals)
        mat['bounds'] = remove_None(self.bounds)
        mat['x0'] = remove_None(self.x0)
        mat['goal_fun_value'] = remove_None(self.goal_fun_value)
        mat['AUC'] = remove_None(self.AUC)
        mat['integral_focality'] = remove_None(self.integral_focality)
        mat['optimizer_options'] = remove_None(self.optimizer_options)
        mat['optimizer_options_std'] = remove_None(self.optimizer_options_std)

        # FEM
        mat['run_final_electrode_simulation'] = remove_None(self.run_final_electrode_simulation)
        mat['dirichlet_node'] = remove_None(self.dirichlet_node)
        mat['dataType'] = remove_None(self.dataType)
        mat['anisotropy_type'] = remove_None(self.anisotropy_type)
        mat['solver_options'] = remove_None(self.solver_options)

        return mat

    def read_mat_element(self, mat, tag):
        """
        Read element from matlab structure and return content

        Paramters
        ---------
        mat : dict
            Dictionary read from scipy.io.loadmat(fn_mat, simplify_cells=True)
        key : str
            Field tag to read

        Returns
        -------
        out : list, dict, np.ndarray
            Content of field
        """

        if type(mat[tag]) is not int and type(mat[tag]) is not float and \
            len(mat[tag]) == 0 and type(mat[tag]) == np.ndarray:
                return None
        else:
            return mat[tag]

    @classmethod
    def read_mat_struct(self, mat):
        """
        Reads parameters from matlab structure

        Parameters
        ----------
        mat: str or scipy.io.loadmat
            Filename of .mat file or loaded matlab structure
        """
        self = self()

        if type(mat) is str:
            mat = scipy.io.loadmat(mat, simplify_cells=True)

        # cond
        # SimuList.read_cond_mat_struct(self, mat)

        # folders and I/O
        self.date = self.read_mat_element(mat, 'date')
        self.output_folder = self.read_mat_element(mat, 'output_folder')
        self.plot_folder = self.read_mat_element(mat, 'plot_folder')
        self.plot = self.read_mat_element(mat, 'plot')
        self.fn_final_sim = self.read_mat_element(mat, 'fn_final_sim')
        self.fn_results_hdf5 = self.read_mat_element(mat, 'fn_results_hdf5')
        self.prepared = self.read_mat_element(mat, 'prepared')

        # headmodel
        self.fn_eeg_cap = self.read_mat_element(mat, 'fn_eeg_cap')
        self.fn_mesh = self.read_mat_element(mat, 'fn_mesh')
        self.fn_electrode_mask = self.read_mat_element(mat, 'fn_electrode_mask')
        self.mesh = mesh_io.read_msh(self.fn_mesh)

        # roi
        self.n_roi =self.read_mat_element(mat, 'n_roi')
        roi_dict = self.read_mat_element(mat, 'roi_dict')

        self.roi = []
        for i_roi in range(self.n_roi):
            roi_ = RegionOfInterestInitializer()
            roi_.type = roi_dict[f"roi_{i_roi}"]["type"]
            roi_.roi_sphere_radius = roi_dict[f"roi_{i_roi}"]["roi_sphere_radius"]

            if len(roi_dict[f"roi_{i_roi}"]["roi_sphere_center_mni"]) == 0:
                roi_.roi_sphere_center_mni = None
            else:
                roi_.roi_sphere_center_mni = roi_dict[f"roi_{i_roi}"]["roi_sphere_center_mni"]

            if len(roi_dict[f"roi_{i_roi}"]["roi_sphere_center_subject"]) == 0:
                roi_.roi_sphere_center_subject = None
            else:
                roi_.roi_sphere_center_subject = roi_dict[f"roi_{i_roi}"]["roi_sphere_center_subject"]

            if len(roi_dict[f"roi_{i_roi}"]["center"]) == 0:
                roi_.center = None
            else:
                roi_.center = roi_dict[f"roi_{i_roi}"]["center"]

            if len(roi_dict[f"roi_{i_roi}"]["nodes"]) == 0:
                roi_.nodes = None
            else:
                roi_.nodes = roi_dict[f"roi_{i_roi}"]["nodes"]

            if len(roi_dict[f"roi_{i_roi}"]["con"]) == 0:
                roi_.con = None
            else:
                roi_.con = roi_dict[f"roi_{i_roi}"]["con"]

            if len(roi_dict[f"roi_{i_roi}"]["domains"]) == 0:
                roi_.domains = None
            else:
                roi_.domains = roi_dict[f"roi_{i_roi}"]["domains"]

            roi_.mesh = self.mesh

            # reinitialize ROI
            self.roi.append(roi_)

        # electrode
        self.electrode_pos = self.read_mat_element(mat, 'electrode_pos')
        self.electrode_pos_opt = self.read_mat_element(mat, 'electrode_pos_opt')
        self.min_electrode_distance = self.read_mat_element(mat, 'min_electrode_distance')
        self.n_channel_stim = self.read_mat_element(mat, 'n_channel_stim')
        self.n_iter_dirichlet_correction = self.read_mat_element(mat, 'n_iter_dirichlet_correction')
        self.n_ele_free = self.read_mat_element(mat, 'n_ele_free')
        self.init_pos = self.read_mat_element(mat, 'init_pos')
        self.init_pos_subject_coords = self.read_mat_element(mat, 'init_pos_subject_coords')

        electrode_dict = mat['electrode_dict']
        self.electrode = []

        for i_ele in range(self.n_channel_stim):
            electrode_ = ElectrodeInitializer()
            electrode_.type = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "type")
            electrode_.dirichlet_correction = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "dirichlet_correction")
            electrode_.dirichlet_correction_detailed = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "dirichlet_correction_detailed")
            electrode_.current_estimator_method = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "current_estimator_method")
            electrode_.current = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "current")

            if electrode_.type == "CircularArray":
                electrode_.radius_inner = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "radius_inner")
                electrode_.radius_inner_bounds = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "radius_inner_bounds")
                electrode_.radius_outer = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "radius_outer")
                electrode_.radius_outer_bounds = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "radius_outer_bounds")
                electrode_.distance = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "distance")
                electrode_.distance_bounds = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "distance_bounds")
                electrode_.n_outer = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "n_outer")
                electrode_.n_outer_bounds = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "n_outer_bounds")

                if np.all(electrode_.radius_inner_bounds == electrode_.radius_inner_bounds[0]):
                    electrode_.radius_inner_bounds = None

                if np.all(electrode_.radius_outer_bounds == electrode_.radius_outer_bounds[0]):
                    electrode_.radius_outer_bounds = None

                if np.all(electrode_.distance_bounds == electrode_.distance_bounds[0]):
                    electrode_.distance_bounds = None

                if np.all(electrode_.n_outer_bounds == electrode_.n_outer_bounds[0]):
                    electrode_.n_outer_bounds = None

            elif electrode_.type == "ElectrodeArrayPair":
                electrode_.center = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "center")
                electrode_.radius = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "radius")
                electrode_.radius_bounds = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "radius_bounds")
                electrode_.length_x = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "length_x")
                electrode_.length_x_bounds = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "length_x_bounds")
                electrode_.length_y = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "length_y")
                electrode_.length_y_bounds = self.read_mat_element(electrode_dict[f"ele_{i_ele}"], "length_y_bounds")

                if np.all(electrode_.radius_bounds == electrode_.radius_bounds[0]):
                    electrode_.radius_bounds = None

                if np.all(electrode_.length_x_bounds == electrode_.length_x_bounds[0]):
                    electrode_.length_x_bounds = None

                if np.all(electrode_.length_y_bounds == electrode_.length_y_bounds[0]):
                    electrode_.length_y_bounds = None

            # reinitialize ROI
            self.electrode.append(electrode_.initialize())

        # goal function
        self.goal = self.read_mat_element(mat, 'goal')
        self.goal_dir = self.read_mat_element(mat, 'goal_dir')
        self.e_postproc = self.read_mat_element(mat, 'e_postproc')
        self.threshold = self.read_mat_element(mat, 'threshold')
        self.optimizer = self.read_mat_element(mat, 'optimizer')
        self.weights = self.read_mat_element(mat, 'weights')
        self.track_focality = self.read_mat_element(mat, 'track_focality')
        self.constrain_electrode_locations = self.read_mat_element(mat, 'constrain_electrode_locations')
        self.overlap_factor = self.read_mat_element(mat, 'overlap_factor')
        self.polish = self.read_mat_element(mat, 'polish')
        self.n_test = self.read_mat_element(mat, 'n_test')
        self.n_sim = self.read_mat_element(mat, 'n_sim')
        self.optimize_init_vals = self.read_mat_element(mat, 'optimize_init_vals')
        bounds_ = self.read_mat_element(mat, 'bounds')
        self.bounds = Bounds(lb=bounds_["lb"].flatten(), ub=bounds_["ub"].flatten())
        self.x0 = self.read_mat_element(mat, 'x0')
        self.goal_fun_value = self.read_mat_element(mat, 'goal_fun_value')
        self.AUC = self.read_mat_element(mat, 'AUC')
        self.integral_focality = self.read_mat_element(mat, 'integral_focality')
        self.optimizer_options = self.read_mat_element(mat, 'optimizer_options')
        self.optimizer_options_std = self.read_mat_element(mat, 'optimizer_options_std')

        # FEM
        self.run_final_electrode_simulation = self.read_mat_element(mat, 'run_final_electrode_simulation')
        self.dirichlet_node = self.read_mat_element(mat, 'dirichlet_node')
        self.dataType = self.read_mat_element(mat, 'dataType')
        self.anisotropy_type = self.read_mat_element(mat, 'anisotropy_type')
        self.solver_options = self.read_mat_element(mat, 'solver_options')
        self.prepared = False

        return self

    def run(self, cpus=None, allow_multiple_runs=False, save_mat=True, return_n_max=1):
        """
        Runs the tes optimization

        Parameters
        ----------
        cpus : int, optional, default: None
            Number of CPU cores to use (not used here because of pardiso solver)
        allow_multiple_runs: bool, optional, default: False
            Whether to allow multiple runs in one folder. (not implemented yet)
        save_mat: bool, optional, default: True
            Save the ".mat" file of this structure
        return_n_max: int, optional, default: 1
            Return n-th best solutions. (not implemented yet)

        Returns
        --------
        <files>: Results files (.hdf5) in self.output_folder.
        """

        # prepare optimization
        if not self.prepared:
            self._prepare()

        # save structure in .mat format
        if save_mat:
            mat = self.sim_struct2mat()
            scipy.io.savemat(os.path.join(self.output_folder, 'simnibs_simulation_{0}.mat'.format(self.time_str)), mat)

        # run optimization
        ################################################################################################################
        start = time.time()
        if self.optimizer == "direct":
            result = direct(self.goal_fun,
                            bounds=self.optimizer_options_std["bounds"],
                            vol_tol=self.optimizer_options_std["vol_tol"],
                            len_tol=self.optimizer_options_std["len_tol"],
                            f_min_rtol=self.optimizer_options_std["f_min_rtol"],
                            maxiter=self.optimizer_options_std["maxiter"],
                            locally_biased=self.optimizer_options_std["locally_biased"])

        elif self.optimizer == "Nelder-Mead":
            result = minimize(self.goal_fun, self.optimizer_options_std["init_vals"],
                              method='Nelder-Mead',
                              bounds=self.optimizer_options_std["bounds"],
                              options={"disp": self.optimizer_options_std["disp"]})

        elif self.optimizer == "differential_evolution":
            result = differential_evolution(self.goal_fun,
                                            x0=self.optimizer_options_std["init_vals"],
                                            strategy='best1bin',
                                            recombination=self.optimizer_options_std["recombination"],
                                            mutation=self.optimizer_options_std["mutation"],
                                            tol=self.optimizer_options_std["tol"],
                                            maxiter=self.optimizer_options_std["maxiter"],
                                            popsize=self.optimizer_options_std["popsize"],
                                            bounds=self.optimizer_options_std["bounds"],
                                            disp=self.optimizer_options_std["disp"],
                                            polish=False)  # we will decide if to polish afterwards

        elif self.optimizer == "shgo":
            result = shgo(self.goal_fun,
                          bounds=self.optimizer_options_std["bounds"],
                          options={"disp": self.optimizer_options_std["disp"]})

        elif self.optimizer == "basinhopping":
            result = basinhopping(self.goal_fun,
                                  x0=self.optimizer_options_std["init_vals"],
                                  disp=self.optimizer_options_std["disp"])
        else:
            raise NotImplementedError(f"Specified optimization method: '{self.optimizer}' not implemented.")

        self.logger.log(20, f"Optimization finished! Best electrode position: {result.x}")
        fopt_before_polish = result.fun
        stop = time.time()
        t_optimize = stop - start

        # polish optimization
        ################################################################################################################
        if self.polish:
            self.logger.log(20, f"Polishing optimization results!")
            result = minimize(self.goal_fun,
                              x0=result.x,
                              method='L-BFGS-B',
                              bounds=self.optimizer_options_std["bounds"],
                              jac='2-point',
                              options={"finite_diff_rel_step": 0.01})
            self.logger.log(20, f"Optimization finished! Best electrode position: {result.x}")

        # transform electrode pos from array to list of list
        self.electrode_pos_opt = self.get_electrode_pos_from_array(result.x)

        fopt = result.fun
        nfev = result.nfev

        # plot final solution and electrode position (with node-wise dirichlet correction)
        ################################################################################################################
        for _electrode in self.electrode:
            _electrode.dirichlet_correction = True
            _electrode.dirichlet_correction_detailed = True

        # compute best e-field again, plot field and electrode position
        e = self.update_field(electrode_pos=self.electrode_pos_opt, plot=True)

        # postprocess e-field
        e_pp = [[0 for _ in range(self.n_roi)] for _ in range(self.n_channel_stim)]
        e_plot = [[] for _ in range(self.n_roi)]
        e_plot_label = [[] for _ in range(self.n_roi)]

        if np.array(["TI" in _t for _t in self.e_postproc]).any():
            for i_roi in range(self.n_roi):
                e_pp[0][i_roi] = postprocess_e(e=e[0][i_roi],
                                               e2=e[1][i_roi],
                                               dirvec=self.goal_dir[i_roi],
                                               type=self.e_postproc[i_roi])
                e_pp[1][i_roi] = e_pp[0][i_roi]

                e_plot[i_roi].append(e[0][i_roi])
                e_plot[i_roi].append(e[0][i_roi])
                e_plot[i_roi].append(e_pp[0][i_roi])
                e_plot_label[i_roi].append(f"e_stim_0")
                e_plot_label[i_roi].append(f"e_stim_1")
                e_plot_label[i_roi].append(f"e_pp")

                # plot field
                if self.plot:
                    fn_out = os.path.join(self.plot_folder, f"e_roi_{i_roi}")
                    plot_roi_field(e=e_plot[i_roi], roi=self.roi[i_roi], e_label=e_plot_label[i_roi], fn_out=fn_out)
        else:
            for i_roi in range(self.n_roi):
                for i_channel_stim in range(self.n_channel_stim):
                    e_pp[i_channel_stim][i_roi] = postprocess_e(e=e[i_channel_stim][i_roi],
                                                                e2=None,
                                                                dirvec=self.goal_dir[i_roi],
                                                                type=self.e_postproc[i_roi])
                    e_plot[i_roi].append(e[i_channel_stim][i_roi])
                    e_plot[i_roi].append(e_pp[i_channel_stim][i_roi])
                    e_plot_label[i_roi].append(f"e_stim_{i_channel_stim}")
                    e_plot_label[i_roi].append(f"e_pp_stim_{i_channel_stim}")

                if self.plot:
                    fn_out = os.path.join(self.plot_folder, f"e_roi_{i_roi}")
                    plot_roi_field(e=e_plot[i_roi], roi=self.roi[i_roi], e_label=e_plot_label[i_roi], fn_out=fn_out)

        # run final simulation with real electrode including remeshing
        ################################################################################################################
        if self.run_final_electrode_simulation:
            for i_channel_stim in range(self.n_channel_stim):
                s = create_tdcs_session_from_array(ElectrodeArray=self.electrode[i_channel_stim],
                                                   fnamehead=self.mesh.fn,
                                                   pathfem=os.path.join(self.output_folder,
                                                                        f"final_sim_{i_channel_stim}"))
                self.fn_final_sim.append(run_simnibs(s)[0])

        # print optimization summary
        save_optimization_results(fname=os.path.join(self.output_folder, "summary"),
                                  optimizer=self.optimizer,
                                  optimizer_options=self.optimizer_options_std,
                                  fopt=fopt,
                                  fopt_before_polish=fopt_before_polish,
                                  popt=self.electrode_pos_opt,
                                  nfev=nfev,
                                  e=e,
                                  e_pp=e_pp,
                                  time=t_optimize,
                                  msh=self.mesh,
                                  electrode=self.electrode,
                                  goal=self.goal,
                                  n_test=self.n_test,
                                  n_sim=self.n_sim,
                                  n_iter_dirichlet_correction=self.n_iter_dirichlet_correction,
                                  goal_fun_value=self.goal_fun_value,
                                  AUC=self.AUC,
                                  integral_focality=self.integral_focality)

    def goal_fun(self, parameters):
        """
        Run function for optimization algorithms.

        Parameters
        ----------
        parameters : np.ndarray of float [n_channel_stim * n_free_arrays * 3]
            Electrodes positions in spherical coordinates (theta, phi, alpha) for each freely movable array.
            e.g.: np.array([theta_stim_1_1, phi_stim_1_1, alpha_stim_1_1, theta_stim_1_2, phi_stim_1_2, alpha_2, ...])

        Returns
        -------
        y : float
            Goal function value
        """
        self.n_test += 1
        parameters_str = f"Parameters: {parameters}"
        self.logger.log(20, parameters_str)

        # transform electrode pos from array to list of list
        self.electrode_pos = self.get_electrode_pos_from_array(parameters)

        # # extract electrode positions from parameters
        # self.electrode_pos = [[] for _ in range(self.n_channel_stim)]
        #
        # i_para = 0
        # for i_channel_stim in range(self.n_channel_stim):
        #     for i_ele_free in range(self.n_ele_free[i_channel_stim]):
        #         if self.electrode[i_channel_stim].electrode_arrays[i_ele_free].optimize_alpha:
        #             i_para_increment = 3
        #         else:
        #             i_para_increment = 2
        #         self.electrode_pos[i_channel_stim].append(parameters[i_para:(i_para + i_para_increment)])
        #         i_para += i_para_increment
        #
        # # extract geometrical electrode parameters from optimal parameters and update electrode
        # for i_channel_stim in range(self.n_channel_stim):
        #     if self.electrode[i_channel_stim].any_free_geometry:
        #         n_free_parameters = np.sum(self.electrode[i_channel_stim].free_geometry)
        #         self.electrode[i_channel_stim].set_geometrical_parameters_optimization(
        #             parameters[i_para:(i_para + n_free_parameters)])
        #         i_para += n_free_parameters

        # update field, returns list of list e[n_channel_stim][n_roi] (None if position is not applicable)
        e = self.update_field(electrode_pos=self.electrode_pos, plot=False)

        # post-process raw electric field (components Ex, Ey, Ez)
        if e is None:
            e_pp = None
        else:
            e_pp = [[0 for _ in range(self.n_roi)] for _ in range(self.n_channel_stim)]
            if np.array(["TI" in _t for _t in self.e_postproc]).any():
                for i_roi in range(self.n_roi):
                    e_pp[0][i_roi] = postprocess_e(e=e[0][i_roi],
                                                   e2=e[1][i_roi],
                                                   dirvec=self.goal_dir[i_roi],
                                                   type=self.e_postproc[i_roi])
                    e_pp[1][i_roi] = e_pp[0][i_roi]
            else:
                for i_channel_stim in range(self.n_channel_stim):
                    for i_roi in range(self.n_roi):
                        e_pp[i_channel_stim][i_roi] = postprocess_e(e=e[i_channel_stim][i_roi],
                                                                    e2=None,
                                                                    dirvec=self.goal_dir[i_roi],
                                                                    type=self.e_postproc[i_roi])

        # compute goal function value
        goal_fun_value = self.compute_goal(e_pp)

        self.logger.log(20, f"Goal ({self.goal}): {goal_fun_value:.3f} (n_sim: {self.n_sim}, n_test: {self.n_test})")
        self.logger.log(20, "-" * len(parameters_str))

        return goal_fun_value

    def compute_goal(self, e):
        """
        Computes goal function value from electric field.

        Parameters
        ----------
        e: list of list of np.ndarrays of float [n_channel_stim][n_roi][n_roi_nodes]
            Post-processed electric fields from simulated simulation conditions and ROIs.

        Returns
        -------
        goal_fun_value : float
            Accumulated goal function value. The average is taken over all stimulation conditions and the weighted
            average is taken according to self.weights over the different goal functions of the ROIs.
        """
        # calculate goal function value for every ROI
        y = np.zeros((self.n_channel_stim, self.n_roi))  # shape: [n_channel_stim x n_roi]

        if e is None:
            self.logger.log(20, f"Goal ({self.goal}): 2.0")
            return 2.0

        # user provided goal function
        ################################################################################################################
        elif isinstance(self.goal[0], types.FunctionType):
            y_weighted_sum = self.goal[0](e)

        # implemented goal functions
        else:
            # focality based goal functions
            ############################################################################################################
            if "focality" in self.goal:
                # TI focality (total field was previously calculated by the 2 channels, no loop over channel_stim here)
                if np.array(["TI" in _t for _t in self.e_postproc]).any():
                    y[:, :] = -100 * (np.sqrt(2) - ROC(e1=e[0][0],  # e-field in ROI
                                                       e2=e[0][1],  # e-field in non-ROI
                                                       threshold=self.threshold,
                                                       focal=True))

                # General focality (can be different for each channel for e.g. TTF)
                else:
                    for i_channel_stim in range(self.n_channel_stim):
                        y[i_channel_stim, :] = -100 * (np.sqrt(2) - ROC(e1=e[i_channel_stim][0],  # e-field in ROI
                                                                        e2=e[i_channel_stim][1],  # e-field in non-ROI
                                                                        threshold=self.threshold,
                                                                        focal=True))

            elif "focality_inv" in self.goal:
                # TI focality (total field was previously calculated by the 2 channels, no loop over channel_stim here)
                if np.array(["TI" in _t for _t in self.e_postproc]).any():
                    y[:, :] = -100 * ROC(e1=e[0][0],  # e-field in ROI
                                         e2=e[0][1],  # e-field in non-ROI
                                         threshold=self.threshold,
                                         focal=False)

                # General focality (can be different for each channel for e.g. TTF)
                else:
                    for i_channel_stim in range(self.n_channel_stim):
                        y[i_channel_stim, :] = -100 * ROC(e1=e[i_channel_stim][0],  # e-field in ROI
                                                          e2=e[i_channel_stim][1],  # e-field in non-ROI
                                                          threshold=self.threshold,
                                                          focal=False)

            # Mean/Max/ etc. based goal functions
            ############################################################################################################
            else:
                for i_roi in range(self.n_roi):
                    for i_channel_stim in range(self.n_channel_stim):
                        if e[i_channel_stim][i_roi] is None:
                            self.logger.log(20, f"Goal ({self.goal}): 2.0 (one e-field was None)")
                            return 2.0
                        else:
                            # mean electric field in the roi
                            if self.goal[i_roi] == "mean":
                                y[i_channel_stim, i_roi] = -np.mean(e[i_channel_stim][i_roi])
                            # negative mean electric field in the roi (for e.g. normal component)
                            elif self.goal[i_roi] == "neg_mean":
                                y[i_channel_stim, i_roi] = np.mean(e[i_channel_stim][i_roi])
                            # mean of absolute value of electric field in the roi (for e.g. normal component)
                            elif self.goal[i_roi] == "mean_abs":
                                y[i_channel_stim, i_roi] = -np.mean(np.abs(e[i_channel_stim][i_roi]))
                            # max electric field in the roi (percentile)
                            elif self.goal[i_roi] == "max":
                                y[i_channel_stim, i_roi] = -np.percentile(e[i_channel_stim][i_roi], 99.9)
                            # negative max electric field in the roi (percentile)
                            elif self.goal[i_roi] == "neg_max":
                                y[i_channel_stim, i_roi] = np.percentile(e[i_channel_stim][i_roi], 99.9)
                            # max of absolute value of electric field in the roi (percentile)
                            elif self.goal[i_roi] == "max_abs":
                                y[i_channel_stim, i_roi] = -np.percentile(np.abs(e[i_channel_stim][i_roi]), 99.9)

            # if desired, track focality measures and goal function values
            for i_channel_stim in range(self.n_channel_stim):
                if self.track_focality:
                    if np.array(["max_TI" in _t for _t in self.e_postproc]).any():
                        e1 = get_maxTI(E1_org=e[0][0], E2_org=e[1][0])
                        e2 = get_maxTI(E1_org=e[0][1], E2_org=e[1][1])
                    elif np.array(["dir_TI" in _t for _t in self.e_postproc]).any():
                        e1 = get_dirTI(E1=e[0][0], E2=e[1][0], dirvec_org=self.goal_dir)
                        e2 = get_dirTI(E1=e[0][1], E2=e[1][1], dirvec_org=self.goal_dir)
                    else:
                        e1 = e[i_channel_stim][0]
                        e2 = e[i_channel_stim][1]

                    # compute integral focality
                    self.integral_focality[i_channel_stim].append(
                        integral_focality(e1=e1,
                                          e2=e2,
                                          v1=self.roi[0].vol,
                                          v2=self.roi[1].vol))

                    # compute auc
                    self.AUC[i_channel_stim].append(AUC(e1=e1,
                                                        e2=e2))
                else:
                    self.AUC[i_channel_stim].append(0)
                    self.integral_focality[i_channel_stim].append(0)

                # goal fun value in roi 0
                self.goal_fun_value[i_channel_stim].append(np.mean(y[i_channel_stim, 0]))

            # average over all stimulations (channel_stim)
            ############################################################################################################
            y = np.mean(y, axis=0)

            # weight and sum the goal function values of the ROIs
            ############################################################################################################
            y_weighted_sum = np.sum(y * self.weights)

        self.n_sim += 1

        return y_weighted_sum

    def get_bounds(self, constrain_electrode_locations, overlap_factor=1.):
        """
        Get boundaries of freely movable electrode arrays for optimizer.

        Parameters
        ----------
        constrain_electrode_locations : bool
            Constrains the possible locations of freely movable electrode arrays. Recommended for TTF optimizations,
            where two pairs of large electrode arrays are optimized. If True, parameter bounds for the optimization
            will be specified restricting the array locations to be frontal, parietal or occipital.
        overlap_factor : float, optional, default: 1.
            Factor of overlap of allowed lambda regions to place electrodes. (1 corresponds to neatless regions,
            <1 regions have a gap in between, >1 regions are overlapping)

        Returns
        -------
        bounds : scipy.optimize.Bounds instance [n_ele_free]
            Boundaries of freely movable electrode arrays tuple of length [n_ele_free] with lower bounds (lb) and
            upper bounds (ub) of beta, lambda, alpha.
        """

        if constrain_electrode_locations:
            # read fiducials from subject data
            fn_fiducials = os.path.join(self.ff_subject.eeg_cap_folder, "Fiducials.csv")
            with open(fn_fiducials, newline='') as f:
                reader = csv.reader(f)
                for row in reader:
                    if "Fiducial" in row and "Nz" in row:
                        # Nz (Nasion), Iz (Inion), LPHA (left ear), RPA (right ear)
                        Nz = np.array([row[1], row[2], row[3]]).astype(float)

            # project nasion to skin surface and determine normal vector
            con_skin = self.mesh.elm.node_number_list[self.mesh.elm.tag1 == 1005,][:, :3] - 1
            tri_skin_center = np.mean(self.mesh.nodes.node_coord[con_skin,], axis=1)
            idx_min = np.argmin(np.linalg.norm(tri_skin_center - Nz, axis=1))
            p1_tri = self.mesh.nodes.node_coord[con_skin[:, 0], :]
            p2_tri = self.mesh.nodes.node_coord[con_skin[:, 1], :]
            p3_tri = self.mesh.nodes.node_coord[con_skin[:, 2], :]
            tri_normal = np.cross(p2_tri - p1_tri, p3_tri - p1_tri)
            tri_normal /= np.linalg.norm(tri_normal, axis=1)[:, np.newaxis]
            Nz_normal = tri_normal[idx_min]

            # project nasion from skin surface to ellipsoid
            Nz_eli = subject2ellipsoid(coords=Nz, normals=Nz_normal, ellipsoid=self.ellipsoid)
            Nz_cart = self.ellipsoid.ellipsoid2cartesian(coords=Nz_eli, norm=False)
            Nz_jacobi = self.ellipsoid.cartesian2jacobi(coords=Nz_cart, norm=False)

            # go a small step into positive z-direction
            Nz_cart_test = self.ellipsoid.jacobi2cartesian(coords=Nz_jacobi + np.array([0, 1e-2]), norm=False)

            if Nz_cart_test[0, 2] > Nz_cart[0, 2]:
                lambda_sign = 1
            else:
                lambda_sign = -1

            # if Nz_jacobi[0, 1] + np.pi / 2 > np.pi:
            #     lambda_sign = -1
            # else:
            #     lambda_sign = 1
            beta_min = np.array([-np.pi / 4, -np.pi / 4, -np.pi / 2, np.pi / 4])
            beta_max = np.array([np.pi / 4, np.pi / 4, -np.pi / 4, np.pi / 2])
            lambda_min = np.array(
                [Nz_jacobi[0, 1], Nz_jacobi[0, 1] + lambda_sign * (np.pi / 2 + np.pi / 16), -np.pi, -np.pi])
            lambda_max = np.array([Nz_jacobi[0, 1] + lambda_sign * (np.pi / 2 + np.pi / 16),
                                   Nz_jacobi[0, 1] + lambda_sign * 3 * np.pi / 2 - np.pi / 8, np.pi, np.pi])
            alpha_min = -np.pi * np.ones(np.sum(self.n_ele_free))
            alpha_max = np.pi * np.ones(np.sum(self.n_ele_free))

            # rearrange reg_lam_center such that electrode array pairs for stim on opposite sites are coming next
            # to each other in order
            lb = np.ravel([beta_min, lambda_min, alpha_min], 'F')
            ub = np.ravel([beta_max, lambda_max, alpha_max], 'F')

            # sort bounds (min, max)
            lbub = np.sort(np.vstack((lb, ub)), axis=0)
            lb = lbub[0, :]
            ub = lbub[1, :]

            # write txt files with some points for visualization
            beta_region_0 = np.linspace(lb[0], ub[0], 100)
            beta_region_1 = np.linspace(lb[3], ub[3], 100)
            beta_region_2 = np.linspace(lb[6], ub[6], 100)
            beta_region_3 = np.linspace(lb[9], ub[9], 100)
            lam_region_0 = np.linspace(lb[1], ub[1], 100)
            lam_region_1 = np.linspace(lb[4], ub[4], 100)
            lam_region_2 = np.linspace(lb[7], ub[7], 100)
            lam_region_3 = np.linspace(lb[10], ub[10], 100)

            coords_region_0_jac = np.array(np.meshgrid(beta_region_0, lam_region_0)).T.reshape(-1, 2)
            coords_region_1_jac = np.array(np.meshgrid(beta_region_1, lam_region_1)).T.reshape(-1, 2)
            coords_region_2_jac = np.array(np.meshgrid(beta_region_2, lam_region_2)).T.reshape(-1, 2)
            coords_region_3_jac = np.array(np.meshgrid(beta_region_3, lam_region_3)).T.reshape(-1, 2)

            coords_region_0 = self.ellipsoid.jacobi2cartesian(coords=coords_region_0_jac, return_normal=False)
            coords_region_1 = self.ellipsoid.jacobi2cartesian(coords=coords_region_1_jac, return_normal=False)
            coords_region_2 = self.ellipsoid.jacobi2cartesian(coords=coords_region_2_jac, return_normal=False)
            coords_region_3 = self.ellipsoid.jacobi2cartesian(coords=coords_region_3_jac, return_normal=False)

            np.savetxt(os.path.join(self.plot_folder, "coords_region_0.txt"), coords_region_0)
            np.savetxt(os.path.join(self.plot_folder, "coords_region_1.txt"), coords_region_1)
            np.savetxt(os.path.join(self.plot_folder, "coords_region_2.txt"), coords_region_2)
            np.savetxt(os.path.join(self.plot_folder, "coords_region_3.txt"), coords_region_3)

        else:
            # beta, lambda, alpha
            # order: [stim_1_array_1, stim_1_array_2, stim_2_array_1, stim_2_array_2, ... ]
            lb = np.array([-np.pi / 2, -np.pi, -np.pi] * np.sum(self.n_ele_free))
            ub = np.array([np.pi / 2, np.pi, np.pi] * np.sum(self.n_ele_free))

        # check if optimize_alpha is set for movable electrode_arrays, if not, remove alpha from the bounds
        i_para = 2
        idx_alpha_remove = []
        for i_channel_stim in range(self.n_channel_stim):
            for _electrode_array in self.electrode[i_channel_stim].electrode_arrays:
                if not _electrode_array.optimize_alpha:
                    idx_alpha_remove.append(i_para)
                i_para += 3
        lb = np.delete(lb, idx_alpha_remove)
        ub = np.delete(ub, idx_alpha_remove)

        # TODO: think this works only for one channel_stim right now (HDTES), test with 2 channel stim and adapt
        # add bounds of geometry parameters of electrode (if any)
        for i_channel_stim in range(self.n_channel_stim):
            if self.electrode[i_channel_stim].any_free_geometry:
                lb = np.append(lb, self.electrode[i_channel_stim].geo_para_bounds[
                    self.electrode[i_channel_stim].free_geometry, 0])
                ub = np.append(ub, self.electrode[i_channel_stim].geo_para_bounds[
                    self.electrode[i_channel_stim].free_geometry, 1])

        bounds = Bounds(lb=lb, ub=ub)

        return bounds

    def get_electrode_pos_from_array(self, electrode_pos_array):
        """
        Transforms an electrode_pos in array (1D) format to a list [n_channel_stim] of list [n_electrode_arrays] of
        numpy arrays (2 or 3, i.e. with or without alpha optimization).
        Changes electrode geometry in place in self.electrode in case of geometrical optimization.

        Parameters:
        -----------
        electrode_pos_array : np.array of float
            Electrode position in array format

        Returns
        -------
        electrode_pos : list of list of np.ndarray of float [2 or 3] of length [n_channel_stim][n_ele_free]
            Spherical coordinates (beta, lambda) and orientation angle (alpha) for each electrode array.
                      electrode array 1                        electrode array 2
            [ np.array([beta_1, lambda_1, alpha_1]),   np.array([beta_2, lambda_2, alpha_2]) ]
        """

        # extract electrode positions from optimal parameters
        electrode_pos = [[] for _ in range(self.n_channel_stim)]

        i_para = 0
        for i_channel_stim in range(self.n_channel_stim):
            for i_ele_free in range(self.n_ele_free[i_channel_stim]):
                if self.electrode[i_channel_stim].electrode_arrays[i_ele_free].optimize_alpha:
                    i_para_increment = 3
                else:
                    i_para_increment = 2
                electrode_pos[i_channel_stim].append(electrode_pos_array[i_para:(i_para + i_para_increment)])
                i_para += i_para_increment

        # TODO: same here I think it only works for one channel_stim
        # extract geometrical electrode parameters from optimal parameters and update electrode
        for i_channel_stim in range(self.n_channel_stim):
            if self.electrode[i_channel_stim].any_free_geometry:
                n_free_parameters = np.sum(self.electrode[i_channel_stim].free_geometry)
                self.electrode[i_channel_stim].set_geometrical_parameters_optimization(
                    electrode_pos_array[i_para:(i_para + n_free_parameters)])
                i_para += n_free_parameters

        return electrode_pos

    def get_init_vals(self, bounds, optimize=True):
        """
        Determine initial values for optimization, guaranteeing a valid electrode position.

        Parameters
        ----------
        bounds : Bounds instance
            Lower and upper bounds of optimization problem
        optimize : bool, optional, default: True
            If True, find initial values for optimization, guaranteeing a valid solution. If False, initial values
            are the center between bounds.

        Returns
        -------
        x0 : ndarray of float [n_para]
            Initial values
        """

        if optimize:
            self.logger.log(20, "Finding valid initial values for electrode position for optimization.")
            n_max = 5000
            n_para = len(bounds.lb)

            # make a list of possible parameter combinations within bounds
            para_test_grid = np.random.rand(n_max, n_para)
            para_test_grid = para_test_grid * (bounds.ub - bounds.lb) + bounds.lb
            para_test_grid_orig = copy.deepcopy(para_test_grid)

            for i in range(n_max):
                # transform electrode pos from array to list of list
                electrode_pos = self.get_electrode_pos_from_array(para_test_grid[i, :])

                # test position
                node_idx_dict = self.get_nodes_electrode(electrode_pos=electrode_pos)

                if type(node_idx_dict[0]) is str:
                    valid = False
                    electrode_pos_valid = node_idx_dict[1]

                    # write valid electrode positions into test grid to ease future iterations
                    i_para = 0
                    for i_channel_stim in range(self.n_channel_stim):
                        for i_ele_free in range(self.n_ele_free[i_channel_stim]):
                            if self.electrode[i_channel_stim].electrode_arrays[i_ele_free].optimize_alpha:
                                i_para_increment = 3
                            else:
                                i_para_increment = 2

                            if electrode_pos_valid[i_channel_stim][i_ele_free] is not None:
                                para_test_grid[i:, i_para:(i_para + i_para_increment)] = electrode_pos_valid[i_channel_stim][
                                    i_ele_free]
                            else:
                                para_test_grid[:, i_para:(i_para + i_para_increment)] = para_test_grid_orig[:, i_para:(i_para + i_para_increment)]
                            i_para += i_para_increment
                else:
                    valid = True

                self.logger.log(20, f"Testing position #{i + 1}: {para_test_grid[i, :]} -> {valid}")

                if not valid:
                    self.logger.log(20, f"> electrode_pos_valid: {node_idx_dict[1]}")

                if valid:
                    return para_test_grid[i, :]

        else:
            if self.n_channel_stim > 1 and type(self.init_pos) is str:
                raise AssertionError("Please provide a list of initial positions for each stimulation channel"
                                     "containing a list of initial positions for each freely movable array.")

            if self.n_channel_stim == 1 and type(self.init_pos) is str:
                self.init_pos = [[self.init_pos]]

            if self.n_channel_stim == 1 and type(self.init_pos) is np.ndarray:
                self.init_pos = [[self.init_pos]]

            self.init_pos_subject_coords = [[] for _ in range(self.n_channel_stim)]

            init_pos_list = [["Fz", "Pz", "P3", "F4"],  # init defaults for first i_channel_stim
                             ["C3", "C4", "F3", "P4"]]  # init defaults for second i_channel_stim

            # set initial positions of electrodes if nothing is provided
            assert self.n_channel_stim <= len(init_pos_list), "Please provide initial electrode positions."

            if self.init_pos is None:
                self.init_pos = [0 for _ in range(self.n_channel_stim)]
                for i_channel_stim in range(self.n_channel_stim):
                    if self.n_ele_free[i_channel_stim] > len(init_pos_list[i_channel_stim]):
                        raise NotImplementedError(
                            "Please specify initial coordinates or EEG electrode positions for each"
                            "freely movable electrode array (init_pos)!")
                    self.init_pos[i_channel_stim] = [init_pos_list[i_channel_stim][i]
                                                     for i in range(self.n_ele_free[i_channel_stim])]

            # get subject coordinates of initial positions
            x0 = np.array([])
            for i_channel_stim in range(self.n_channel_stim):
                # user provided EEG electrode position as str (e.g. "C3", ...)
                if type(self.init_pos[i_channel_stim][0]) is str:
                    for eeg_pos in self.init_pos[i_channel_stim]:
                        tmp = ELECTRODE()
                        tmp.centre = eeg_pos
                        tmp.substitute_positions_from_cap(
                            cap=self.ff_subject.get_eeg_cap(cap_name=self.fn_eeg_cap))
                        self.init_pos_subject_coords[i_channel_stim].append(tmp.centre)
                # user provided coordinates in subject space as np.array
                else:
                    self.init_pos_subject_coords[i_channel_stim] = self.init_pos[i_channel_stim]

                # transform initial positions from subject to ellipsoid space
                for i_ele_free, coords in enumerate(self.init_pos_subject_coords[i_channel_stim]):
                    # get closest point idx on subject surface
                    point_idx = np.argmin(np.linalg.norm(coords - self.skin_surface.nodes, axis=1))

                    # electrode positon in ellipsoid space (jacobi coordinates)
                    self.electrode_pos[i_channel_stim][i_ele_free][:2] = self.ellipsoid.cartesian2jacobi(
                        coords=self.ellipsoid.ellipsoid2cartesian(
                            coords=subject2ellipsoid(
                                coords=self.skin_surface.nodes[point_idx, :],
                                normals=self.skin_surface.nodes_normals[point_idx, :],
                                ellipsoid=self.ellipsoid)))

                    # set initial orientation alpha to zero
                    if len(self.electrode_pos[i_channel_stim][i_ele_free]) > 2:
                        self.electrode_pos[i_channel_stim][i_ele_free][2] = 0.

                    # append position to initial values
                    x0 = np.append(x0, self.electrode_pos[i_channel_stim][i_ele_free])

                    # add geometric parameters if applicable
                    x0 = np.append(x0, self.electrode[i_channel_stim].geo_para_mean[self.electrode[i_channel_stim].free_geometry])

        # self.x0 = np.mean(np.vstack((bounds.lb, bounds.ub)), axis=0)
        return x0

    def get_nodes_electrode(self, electrode_pos):
        """
        Assigns the skin points of the electrodes in electrode array and writes the points in
        electrode_array.electrodes[i].nodes and electrode_array.electrodes[i].node_area.
        Estimate optimal electrode currents based on previous simulations.

        Parameters
        ----------
        electrode_pos : list of list of np.ndarray of float [2 or 3] of length [n_channel_stim][n_ele_free]
            Spherical coordinates (beta, lambda) and orientation angle (alpha) for each electrode array.
                      electrode array 1                        electrode array 2
            [ np.array([beta_1, lambda_1, alpha_1]),   np.array([beta_2, lambda_2, alpha_2]) ]

        Returns
        -------
        node_idx_dict : list of dict
            List [n_channel_stim] containing dicts with electrode channel IDs as keys and node indices.
        """

        electrode_coords_subject = [0 for _ in range(self.n_channel_stim)]
        electrode_pos_valid = [[None for _ in range(len(self.electrode[i_channel_stim].electrode_arrays))] for
                               i_channel_stim in range(self.n_channel_stim)]
        node_idx_dict = [dict() for _ in range(self.n_channel_stim)]
        node_coords_list = [[] for _ in range(np.sum(self.n_ele_free))]
        i_array_global = 0

        n = []
        cx = []
        cy = []

        for i_channel_stim in range(self.n_channel_stim):
            # collect all parameters
            start = np.zeros((self.n_ele_free[i_channel_stim], 3))
            a = np.zeros((self.n_ele_free[i_channel_stim], 3))
            b = np.zeros((self.n_ele_free[i_channel_stim], 3))
            cx.append(np.zeros((self.n_ele_free[i_channel_stim], 3)))
            cy.append(np.zeros((self.n_ele_free[i_channel_stim], 3)))
            n_tmp = np.zeros((self.n_ele_free[i_channel_stim], 3))
            start_shifted_ = np.zeros((len(electrode_pos[i_channel_stim]), 3))
            distance = []
            alpha = []

            for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
                start[i_array, :] = self.ellipsoid.jacobi2cartesian(coords=electrode_pos[i_channel_stim][i_array][:2])

                c0, n_tmp[i_array, :] = self.ellipsoid.jacobi2cartesian(
                    coords=electrode_pos[i_channel_stim][i_array][:2], return_normal=True)
                a[i_array, :] = self.ellipsoid.jacobi2cartesian(coords=np.array(
                    [electrode_pos[i_channel_stim][i_array][0] - 1e-2, electrode_pos[i_channel_stim][i_array][1]])) - c0
                b[i_array, :] = self.ellipsoid.jacobi2cartesian(coords=np.array(
                    [electrode_pos[i_channel_stim][i_array][0], electrode_pos[i_channel_stim][i_array][1] - 1e-2])) - c0
                a[i_array, :] /= np.linalg.norm(a[i_array, :])
                b[i_array, :] /= np.linalg.norm(b[i_array, :])

                if len(electrode_pos[i_channel_stim][i_array]) > 2:
                    start_shifted_[i_array, :] = c0 + (
                                1e-3 * ((a[i_array, :]) * np.cos(electrode_pos[i_channel_stim][i_array][2]) +
                                        (b[i_array, :]) * np.sin(electrode_pos[i_channel_stim][i_array][2])))
                else:
                    start_shifted_[i_array, :] = c0 + 1e-3 * a[i_array, :]

                cy[i_channel_stim][i_array, :] = start_shifted_[i_array, :] - start[i_array, :]
                cy[i_channel_stim][i_array, :] /= np.linalg.norm(cy[i_channel_stim][i_array, :])
                cx[i_channel_stim][i_array, :] = np.cross(cy[i_channel_stim][i_array, :], -n_tmp[i_array, :])
                cx[i_channel_stim][i_array, :] /= np.linalg.norm(cx[i_channel_stim][i_array, :])

                distance.append(_electrode_array.distance)
                if _electrode_array.optimize_alpha:
                    alpha.append(electrode_pos[i_channel_stim][i_array][2] + _electrode_array.angle)
                else:
                    alpha.append(_electrode_array.angle)

            distance = np.array(distance).flatten()
            alpha = np.array(alpha).flatten()

            for i_a, _alpha in enumerate(alpha):
                if _alpha > np.pi:
                    alpha[i_a] = _alpha - 2 * np.pi
                elif _alpha < -np.pi:
                    alpha[i_a] = _alpha + 2 * np.pi

            start = np.vstack([np.tile(start[i_array, :], (_electrode_array.n_ele, 1))
                               for i_array, _electrode_array in
                               enumerate(self.electrode[i_channel_stim].electrode_arrays)])

            electrode_array_idx = np.hstack([i_array * np.ones(_electrode_array.n_ele)
                                             for i_array, _electrode_array in
                                             enumerate(self.electrode[i_channel_stim].electrode_arrays)])

            # determine electrode center on ellipsoid
            if not (distance == 0.).all():
                electrode_coords_eli_cart = self.ellipsoid.get_geodesic_destination(start=start,
                                                                                    distance=distance,
                                                                                    alpha=alpha,
                                                                                    n_steps=400)
            else:
                electrode_coords_eli_cart = start

            n.append(self.ellipsoid.get_normal(coords=electrode_coords_eli_cart))

            # transform to ellipsoidal coordinates
            electrode_coords_eli_eli = self.ellipsoid.cartesian2ellipsoid(coords=electrode_coords_eli_cart)

            # project coordinates to subject
            tmp_arrays = []
            i_ele = 0
            for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
                ele_idx, tmp = ellipsoid2subject(coords=electrode_coords_eli_eli[electrode_array_idx == i_array, :],
                                                 ellipsoid=self.ellipsoid,
                                                 surface=self.skin_surface)
                tmp_arrays.append(tmp)

                if len(ele_idx) != len(alpha[electrode_array_idx == i_array]):
                    return "Electrode position: invalid (not all electrodes in valid skin region)", electrode_pos_valid
                else:
                    electrode_pos_valid[i_channel_stim][i_array] = electrode_pos[i_channel_stim][i_array]
                    # print("Electrode position: invalid (not all electrodes in valid skin region)")

                electrode_coords_subject[i_channel_stim] = np.vstack(tmp_arrays)

                # loop over electrodes and determine node indices
                for _electrode in _electrode_array.electrodes:
                    if _electrode.type == "spherical":
                        # mask with a sphere
                        mask = np.linalg.norm(
                            self.skin_surface.nodes - electrode_coords_subject[i_channel_stim][i_ele, :],
                            axis=1) < _electrode.radius

                        # save position of electrode in subject space to posmat field
                        _electrode.posmat[:3, 3] = electrode_coords_subject[i_channel_stim][i_ele, :]

                    elif _electrode.type == "rectangular":
                        cx_local = np.cross(n[i_channel_stim][i_ele, :], cy[i_channel_stim][i_array, :])

                        # rotate skin nodes to normalized electrode space
                        rotmat = np.array([[cx_local[0], cy[i_channel_stim][i_array, 0], n[i_channel_stim][i_ele, 0]],
                                           [cx_local[1], cy[i_channel_stim][i_array, 1], n[i_channel_stim][i_ele, 1]],
                                           [cx_local[2], cy[i_channel_stim][i_array, 2], n[i_channel_stim][i_ele, 2]]])
                        center = np.array([electrode_coords_subject[i_channel_stim][i_ele, 0],
                                           electrode_coords_subject[i_channel_stim][i_ele, 1],
                                           electrode_coords_subject[i_channel_stim][i_ele, 2]])

                        # save position of electrode in subject space to posmat field
                        _electrode.posmat = np.vstack(
                            (np.hstack((rotmat, center[:, np.newaxis])), np.array([0, 0, 0, 1])))

                        skin_nodes_rotated = (self.skin_surface.nodes - center) @ rotmat

                        # mask with a box
                        mask_x = np.logical_and(skin_nodes_rotated[:, 0] > -_electrode.length_x / 2,
                                                skin_nodes_rotated[:, 0] < +_electrode.length_x / 2)
                        mask_y = np.logical_and(skin_nodes_rotated[:, 1] > -_electrode.length_y / 2,
                                                skin_nodes_rotated[:, 1] < +_electrode.length_y / 2)
                        mask_z = np.logical_and(skin_nodes_rotated[:, 2] > -30,
                                                skin_nodes_rotated[:, 2] < +30)
                        mask = np.logical_and(np.logical_and(mask_x, mask_y), mask_z)
                    else:
                        raise AssertionError("Electrodes have to be either 'spherical' or 'rectangular'")

                    # node areas
                    _electrode.node_area = self.skin_surface.nodes_areas[mask]

                    # total effective area of all nodes
                    _electrode.area_skin = _electrode.node_area_total

                    # electrode position is invalid if it overlaps with invalid skin region and area is not "complete"
                    if _electrode.area_skin < 0.90 * _electrode.area:
                        electrode_pos_valid[i_channel_stim][i_array] = None
                        # print("Electrode position: invalid (partly overlaps with invalid skin region)")
                        return "Electrode position: invalid (partly overlaps with invalid skin region)", electrode_pos_valid

                    # save node indices (referring to global mesh)
                    _electrode.node_idx = self.node_idx_msh[mask]

                    # save node coords (refering to global mesh)
                    _electrode.node_coords = self.skin_surface.nodes[mask]

                    # save number of nodes assigned to this electrode
                    # _electrode.n_nodes = len(_electrode.node_idx)

                    node_coords_list[i_array_global].append(_electrode.node_coords)

                    # group node indices of same channel IDs
                    if _electrode.channel_id in node_idx_dict[i_channel_stim].keys():
                        node_idx_dict[i_channel_stim][_electrode.channel_id] = np.append(
                            node_idx_dict[i_channel_stim][_electrode.channel_id], _electrode.node_idx)
                    else:
                        node_idx_dict[i_channel_stim][_electrode.channel_id] = _electrode.node_idx

                    i_ele += 1

                electrode_pos_valid[i_channel_stim][i_array] = electrode_pos[i_channel_stim][i_array]

                # gather all electrode node coords of freely movable arrays
                node_coords_list[i_array_global] = np.vstack(node_coords_list[i_array_global])
                i_array_global += 1

        # check if electrode distance is sufficient
        invalid = False
        i_array_global_lst = np.hstack(
            [np.arange(self.n_ele_free[i_channel_stim]) for i_channel_stim in range(self.n_channel_stim)]).astype(int)
        i_channel_stim_global_lst = np.hstack(
            [i_channel_stim * np.ones(self.n_ele_free[i_channel_stim]) for i_channel_stim in
             range(self.n_channel_stim)]).astype(int)
        if self.min_electrode_distance is not None and self.min_electrode_distance >= 0:
            i_array_test_start = 1
            # start with first array and test if all node coords are too close to other arrays
            for i_array_global in range(np.sum(self.n_ele_free)):
                for node_coord in node_coords_list[i_array_global]:
                    for i_array_test in range(i_array_test_start, np.sum(self.n_ele_free)):
                        # calculate euclidean distance between node coords
                        min_dist = np.min(np.linalg.norm(node_coords_list[i_array_test] - node_coord, axis=1))
                        # stop testing if an electrode is too close
                        if min_dist <= self.min_electrode_distance:
                            # remove tested array from valid list
                            electrode_pos_valid[i_channel_stim_global_lst[i_array_test]][
                                i_array_global_lst[i_array_test]] = None
                            # print("Electrode position: invalid (minimal distance between electrodes too small)")
                            invalid = True

                i_array_test_start += 1

        if invalid:
            return "Electrode position: invalid (minimal distance between electrodes too small)", electrode_pos_valid

        # save electrode_pos in ElectrodeArray instances
        for i_channel_stim in range(self.n_channel_stim):
            for i_array, _electrode_array in enumerate(self.electrode[i_channel_stim].electrode_arrays):
                _electrode_array.electrode_pos = electrode_pos[i_channel_stim][i_array]

        # estimate optimal electrode currents from previous simulations
        for i_channel_stim in range(self.n_channel_stim):
            if self.electrode[i_channel_stim].current_estimator is not None:

                # estimate optimal currents electrode wise
                currents_estimate = self.electrode[i_channel_stim].estimate_currents(electrode_pos[i_channel_stim])

                # write currents in electrodes
                if currents_estimate is not None:
                    currents_estimate = currents_estimate.flatten()
                    for _electrode_array in self.electrode[i_channel_stim].electrode_arrays:
                        for _ele in _electrode_array.electrodes:
                            mask_estimator = (self.electrode[i_channel_stim].current_estimator.ele_id == _ele.ele_id) * \
                                             (self.electrode[
                                                  i_channel_stim].current_estimator.channel_id == _ele.channel_id)
                            _ele.ele_current = currents_estimate[mask_estimator]
            else:
                # reset to original currents
                for _electrode_array in self.electrode[i_channel_stim].electrode_arrays:
                    for _ele in _electrode_array.electrodes:
                        _ele.ele_current = _ele.ele_current_init

                self.electrode[i_channel_stim].compile_node_arrays()

        # compile node arrays
        for i_channel_stim in range(self.n_channel_stim):
            self.electrode[i_channel_stim].compile_node_arrays()

        return node_idx_dict

    def update_field(self, electrode_pos, plot=False):
        """
        Calculate the E field for given electrode positions.

        Parameters
        ----------
        electrode_pos : list of list of np.ndarray of float [3] of length [n_channel_stim][n_ele_free][3]
            Spherical coordinates (beta, lambda) and orientation angle (alpha) for each electrode array.
                      electrode array 1                        electrode array 2
            [ np.array([beta_1, lambda_1, alpha_1]),   np.array([beta_2, lambda_2, alpha_2]) ]
        plot : bool, optional, default: False
            Save data to plot e-field and electrode positions

        Returns
        -------
        e : list of list of np.ndarray [n_channel_stim][n_roi]
            Electric field for different stimulations in ROI(s).
        """

        e = [[0 for _ in range(self.n_roi)] for _ in range(self.n_channel_stim)]

        # assign surface nodes to electrode positions and estimate optimal currents
        # start = time.time()
        node_idx_dict = self.get_nodes_electrode(electrode_pos=electrode_pos)
        # stop = time.time()
        # print(f"Time: get_nodes_electrode: {stop-start}")

        # perform one electric field calculation for every stimulation condition (one at a time is on)
        for i_channel_stim in range(self.n_channel_stim):
            if type(node_idx_dict[0]) is str:
                self.logger.log(20, node_idx_dict[0])
                return None
            self.logger.log(20, f"Electrode position for stimulation {i_channel_stim}: valid")

            # set RHS
            b = self.ofem.set_rhs(electrode=self.electrode[i_channel_stim])

            # solve system
            if self.electrode[i_channel_stim].dirichlet_correction:
                if plot:
                    fn_electrode_txt = os.path.join(self.plot_folder,
                                                    f"electrode_coords_nodes_subject_{i_channel_stim}.txt")
                else:
                    fn_electrode_txt = None

                v = self.ofem.solve_dirichlet_correction(electrode=self.electrode[i_channel_stim],
                                                         fn_electrode_txt=fn_electrode_txt)

                # store number of dirichlet iterations for convergence analysis
                self.n_iter_dirichlet_correction[i_channel_stim].append(self.ofem.n_iter_dirichlet_correction)
            else:
                v = self.ofem.solve(b)

                if plot:
                    fn_electrode_txt = os.path.join(self.plot_folder,
                                                    f"electrode_coords_nodes_subject_{i_channel_stim}.txt")
                    np.savetxt(fn_electrode_txt,
                               np.hstack((self.electrode[i_channel_stim].node_coords,
                                          self.electrode[i_channel_stim].node_current[:, np.newaxis])))

            # Determine electric field in ROIs
            # start = time.time()
            for i_roi, r in enumerate(self.roi):
                if v is None:
                    e[i_channel_stim][i_roi] = None
                    self.logger.log(20, "Warning! Simulation failed! Returning e-field: None!")
                else:
                    e[i_channel_stim][i_roi] = r.calc_fields(v, dataType=self.dataType[i_roi])
            # stop = time.time()
            # print(f"Time: calc fields: {stop - start}")

        return e


def save_optimization_results(fname, optimizer, optimizer_options, fopt, fopt_before_polish, popt, nfev, e, e_pp, time,
                              msh, electrode, goal, n_test=None, n_sim=None, n_iter_dirichlet_correction=None,
                              goal_fun_value=None, AUC=None, integral_focality=None):
    """
    Saves optimization settings and results in an <fname>.hdf5 file and prints a summary in a <fname>.txt file.

    Parameters
    ----------
    fname : str
        Filename of the summary.txt file.
    optimizer : str
        Name of optimization method.
    optimizer_options : dict
        Dictionary containing the optimization setting.
    fopt : float
        Objective function value in optimum.
    popt : list of list of np.ndarray of float [n_channel_stim][n_free_electrodes]
        List containing the optimal parameters for each freely movable electrode array
        [np.array([beta_1, lambda_1, alpha_1]), np.array([beta_2, lambda_2, alpha_2]), ...]
    nfev : int
        Number of function evaluations during optimization
    e : list of np.ndarray [n_channel_stim][n_roi]
        List of list containing np.ndarrays of the (raw) electric fields in the ROIs (Ex, Ey, Ez)
    e_pp : list of np.ndarray [n_channel_stim][n_roi]
        List of list containing np.ndarrays of the postprocessed electric fields in the ROIs (norm or normal etc...)
    time : float
        Runtime of optimization in s
    electrode : list of ElectrodeArray objects [n_channel_stim]
        List of ElectrodeArray objects for every stimulation
    goal : str
        Goal function definition
    n_test : int, optional, default: None
        Number of runs to place the electrodes
    n_sim : int, optional, default: None
        Number of actual FEM simulations (for valid electrode placements only)
    n_iter_dirichlet_correction : list of int [n_channel_stim][n_iter], optional, default: None
        Number of iterations required to determine optimal currents in case of dirichlet correction
        for each call of "solve_dirichlet_correction"
    goal_fun_value : list of list of float [n_channel_stim][n_opt_runs]
        Goal function values of all stimulation conditions during optimization of ROI 0
    AUC : list of list of float [n_channel_stim][n_opt_runs]
        Area under curve focality measure for all stimulation conditions.
    integral_focality : list of list of float [n_channel_stim][n_opt_runs]
        Integral focality measure for all stimulation conditions.
    """

    def sep(x):
        if x >= 0:
            sep = " "
        else:
            sep = ""

        return sep

    fname_txt = fname + ".txt"
    fname_hdf5 = fname + ".hdf5"

    # print summary in <fname>.txt file
    ####################################################################################################################
    d = datetime.datetime.now()

    with open(fname_txt, 'w') as f:
        f.write(f"Optimization summary:\n")
        f.write(f"===================================================================\n")
        f.write(f"Date: {d.year}-{d.month}-{d.day}, {d.hour}:{d.minute}:{d.second}\n")
        f.write(f"Simulation time: {str(datetime.timedelta(seconds=time))[:-7]}\n")
        f.write(f"headmodel: {msh.fn}\n")
        f.write(f"Goal: {goal}\n")
        f.write(f"Number of Channels: {len(e)}\n")
        f.write(f"Number of ROIs: {len(e[0])}\n")
        f.write(f"\n")
        f.write(f"Electrode coordinates:\n")
        f.write(f"===================================================================\n")
        f.write(f"Ellipsoid space (Jacobian coordinates):\n")
        f.write(f"---------------------------------------\n")

        for i_stim in range(len(popt)):
            f.write(f"Stimulation {i_stim}:\n")
            for i, p in enumerate(popt[i_stim]):
                f.write(f"Array {i}:\n")
                f.write(f"\tbeta:   {sep(p[0])}{p[0]:.3f}\n")
                f.write(f"\tlambda: {sep(p[1])}{p[1]:.3f}\n")
                if len(p) > 2:
                    f.write(f"\talpha:  {sep(p[2])}{p[2]:.3f}\n")
                else:
                    f.write(f"\talpha:  {sep(0.000)}{0.000}\n")

                if type(electrode[i_stim]) is CircularArray:
                    f.write(f"\tradius_inner:  {sep(electrode[i_stim].radius_inner)}{electrode[i_stim].radius_inner}\n")
                    f.write(f"\tradius_outer:  {sep(electrode[i_stim].radius_outer)}{electrode[i_stim].radius_outer}\n")
                    f.write(f"\tdistance:  {sep(electrode[i_stim].distance)}{electrode[i_stim].distance}\n")
                    f.write(f"\tn_outer:  {sep(electrode[i_stim].n_outer)}{electrode[i_stim].n_outer}\n")

        f.write(f"\n")
        f.write(f"Subject space (Cartesian coordinates):\n")
        f.write(f"--------------------------------------\n")
        for i_stim in range(len(popt)):
            f.write(f"Stimulation {i_stim}:\n")
            for i_array, _electrode_array in enumerate(electrode[i_stim].electrode_arrays):
                f.write(f"Array {i_array}:\n")
                for i_electrode, _electrode in enumerate(_electrode_array.electrodes):
                    f.write(f"\tElectrode {i_electrode} ({_electrode.type}):\n")
                    for i_row in range(4):
                        f.write("\t\t" + sep(_electrode.posmat[i_row, 0]) + f"{_electrode.posmat[i_row, 0]:.3f}, " +
                                sep(_electrode.posmat[i_row, 1]) + f"{_electrode.posmat[i_row, 1]:.3f}, " +
                                sep(_electrode.posmat[i_row, 2]) + f"{_electrode.posmat[i_row, 2]:.3f}, " +
                                sep(_electrode.posmat[i_row, 3]) + f"{_electrode.posmat[i_row, 3]:.3f}\n")
        f.write("\n")
        f.write("Optimization method:\n")
        f.write("===================================================================\n")
        f.write(f"Optimizer: {optimizer}\n")
        f.write(f"Settings:\n")
        if optimizer_options is not None:
            for key in optimizer_options:
                if type(optimizer_options[key]) is Bounds:
                    f.write(f"\tlb: {optimizer_options[key].lb}\n")
                    f.write(f"\tub: {optimizer_options[key].ub}\n")
                else:
                    f.write(f"\t{key}: {optimizer_options[key]}\n")
        else:
            f.write("None\n")
        f.write(f"nfev: {nfev}\n")
        f.write(f"fopt: {fopt}\n")
        f.write(f"fopt_before_polish: {fopt_before_polish}\n")

        if n_sim is not None:
            f.write(f"n_sim: {n_sim}\n")

        if n_test is not None:
            f.write(f"n_test: {n_test}\n")

    # save optimization settings and results in <fname>.hdf5 file
    ####################################################################################################################
    with h5py.File(fname_hdf5, "w") as f:
        # general info
        f.create_dataset(data=msh.fn, name="fnamehead")
        f.create_dataset(data=f"{d.year}-{d.month}-{d.day}, {d.hour}:{d.minute}:{d.second}", name="date")
        f.create_dataset(data=len(e), name="n_channel")
        f.create_dataset(data=len(e[0]), name="n_roi")

        # optimizer
        f.create_dataset(data=optimizer, name="optimizer/optimizer")
        f.create_dataset(data=fopt, name="optimizer/fopt")
        f.create_dataset(data=fopt_before_polish, name="optimizer/fopt_before_polish")
        f.create_dataset(data=nfev, name="optimizer/nfev")
        f.create_dataset(data=time, name="optimizer/time")
        f.create_dataset(data=goal, name="optimizer/goal")

        if goal_fun_value is not None:
            f.create_dataset(data=np.array(goal_fun_value), name="optimizer/goal_fun_value")

        if AUC is not None:
            f.create_dataset(data=np.array(AUC), name="optimizer/AUC")

        if integral_focality is not None:
            f.create_dataset(data=np.array(integral_focality), name="optimizer/integral_focality")

        if n_iter_dirichlet_correction is not None:
            for i_channel_stim, n_iter in enumerate(n_iter_dirichlet_correction):
                f.create_dataset(data=n_iter, name=f"optimizer/n_iter_dirichlet_correction/channel_{i_channel_stim}")

        if n_sim is not None:
            f.create_dataset(data=n_sim, name="optimizer/n_sim")

        if n_test is not None:
            f.create_dataset(data=n_test, name="optimizer/n_test")

        for key in optimizer_options:
            if type(optimizer_options[key]) is Bounds:
                f.create_dataset(data=optimizer_options[key].lb, name=f"optimizer/optimizer_options/lb")
                f.create_dataset(data=optimizer_options[key].ub, name=f"optimizer/optimizer_options/ub")
            else:
                f.create_dataset(data=optimizer_options[key], name=f"optimizer/optimizer_options/{key}")

        # electrodes
        for i_stim in range(len(electrode)):
            f.create_dataset(data=electrode[i_stim].center, name=f"electrode/channel_{i_stim}/center")
            f.create_dataset(data=electrode[i_stim].radius, name=f"electrode/channel_{i_stim}/radius")
            f.create_dataset(data=electrode[i_stim].length_x, name=f"electrode/channel_{i_stim}/length_x")
            f.create_dataset(data=electrode[i_stim].length_y, name=f"electrode/channel_{i_stim}/length_y")
            f.create_dataset(data=electrode[i_stim].current, name=f"electrode/channel_{i_stim}/current")

            if type(electrode[i_stim]) is CircularArray:
                f.create_dataset(data=electrode[i_stim].radius_inner, name=f"electrode/channel_{i_stim}/radius_inner")
                f.create_dataset(data=electrode[i_stim].radius_outer, name=f"electrode/channel_{i_stim}/radius_outer")
                f.create_dataset(data=electrode[i_stim].distance, name=f"electrode/channel_{i_stim}/distance")
                f.create_dataset(data=electrode[i_stim].n_outer, name=f"electrode/channel_{i_stim}/n_outer")

            for i_array, _electrode_array in enumerate(electrode[i_stim].electrode_arrays):
                f.create_dataset(data=popt[i_stim][i_array][0],
                                 name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/beta")
                f.create_dataset(data=popt[i_stim][i_array][1],
                                 name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/lambda")
                if len(popt[i_stim][i_array]) > 2:
                    f.create_dataset(data=popt[i_stim][i_array][2],
                                     name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/alpha")
                else:
                    f.create_dataset(data=0.0,
                                     name=f"electrode/channel_{i_stim}/popt/electrode_array_{i_array}/alpha")

                for i_electrode, _electrode in enumerate(_electrode_array.electrodes):
                    f.create_dataset(data=_electrode.posmat,
                                     name=f"electrode/channel_{i_stim}/posmat/electrode_array_{i_array}/electrode_{i_electrode}/posmat")

        # electric field in ROIs
        for i_stim in range(len(e)):
            for i_roi in range(len(e[i_stim])):
                f.create_dataset(data=e[i_stim][i_roi], name=f"e/channel_{i_stim}/e_roi_{i_roi}")
                f.create_dataset(data=e_pp[i_stim][i_roi], name=f"e_pp/channel_{i_stim}/e_roi_{i_roi}")


def valid_skin_region(skin_surface, mesh, fn_electrode_mask, additional_distance=0.):
    """
    Determine the nodes of the scalp surface where the electrode can be applied (not ears and face etc.)

    Parameters
    ----------
    skin_surface : Surface object
        Surface of the mesh (mesh_tools/surface.py)
    mesh : Msh object
        Mesh object created by SimNIBS (mesh_tools/mesh_io.py)
    additional_distance : float, optional, default: 0
        Additional distance in anterior part to put between original MNI template registration
    """
    nodes_all = copy.deepcopy(skin_surface.nodes)
    tr_nodes_all = copy.deepcopy(skin_surface.tr_nodes)
    # load mask of valid electrode positions (in MNI space)
    mask_img = nibabel.load(fn_electrode_mask)
    mask_img_data = mask_img.get_fdata()

    # add a certain distance to mask out closer to the eyes
    for i_border in range(mask_img_data.shape[0]):
        if mask_img_data[:, :, i_border].all():
            break

    mask_img_data[:, :, (i_border - additional_distance):i_border] = 1

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

    unique_points = np.unique(
        skin_surface.tr_nodes[skin_surface.mask_valid_nodes[skin_surface.tr_nodes].all(axis=1), :])
    for point in unique_points:
        idx_where = np.where(skin_surface.tr_nodes == point)
        skin_surface.mask_valid_tr[idx_where[0], idx_where[1]] = True
    skin_surface.mask_valid_tr = skin_surface.mask_valid_tr.all(axis=1)

    # determine connectivity list of valid skin region (creates new node and connectivity list)
    skin_surface.nodes, skin_surface.tr_nodes = create_new_connectivity_list_point_mask(
        points=skin_surface.nodes,
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
            nodes_idx_of_domain = np.unique(
                np.append(nodes_idx_of_domain, skin_surface.tr_nodes[tri_idx_of_domain, :])).astype(int)
            tri_idx_of_domain = np.isin(skin_surface.tr_nodes, nodes_idx_of_domain).any(axis=1)
            n_current = np.sum(tri_idx_of_domain)
            # print(f"domain: {domain}, n_current: {n_current}")

        tri_domain[tri_idx_of_domain] = domain
        point_domain[nodes_idx_of_domain] = domain
        domain += 1

    domain_idx_main = np.argmax([np.sum(point_domain == d) for d in range(domain)])

    skin_surface.nodes, skin_surface.tr_nodes = create_new_connectivity_list_point_mask(
        points=skin_surface.nodes,
        con=skin_surface.tr_nodes,
        point_mask=point_domain == domain_idx_main)

    # update masks
    skin_surface.mask_valid_nodes = np.zeros(nodes_all.shape[0]).astype(bool)
    for i_p, p in enumerate(nodes_all):
        if p in skin_surface.nodes:
            skin_surface.mask_valid_nodes[i_p] = True

    skin_surface.mask_valid_tr = np.zeros(tr_nodes_all.shape[0]).astype(bool)
    for i_t, t in enumerate(tr_nodes_all):
        if t in skin_surface.tr_nodes:
            skin_surface.mask_valid_tr[i_t] = True

    skin_surface.nodes_areas = skin_surface.nodes_areas[skin_surface.mask_valid_nodes]
    skin_surface.nodes_normals = skin_surface.nodes_normals[skin_surface.mask_valid_nodes, :]
    skin_surface.surf2msh_nodes = skin_surface.surf2msh_nodes[skin_surface.mask_valid_nodes]
    skin_surface.surf2msh_triangles = skin_surface.surf2msh_triangles[skin_surface.mask_valid_tr]
    skin_surface.tr_areas = skin_surface.tr_areas[skin_surface.mask_valid_tr]
    skin_surface.tr_centers = skin_surface.tr_centers[skin_surface.mask_valid_tr, :]
    skin_surface.tr_normals = skin_surface.tr_normals[skin_surface.mask_valid_tr, :]

    return skin_surface


def relabel_internal_air(m, subpath, label_skin=1005, label_new=1099, label_internal_air=501):
    """
    Relabels skin in internal air cavities to something else; relevant for charm meshes

    Parameters
    ----------
    m : Msh object
        Mesh object with internal air
    subpath :
        Path to subject m2m folder
    label_skin : int
        Original skin label
    label_new : int
        New skin label
    label_internal_air : int
        New label of internal air

    Returns
    -------
    m : Msh object
        Mesh with relabeled internal air
    """
    subject_files = SubjectFiles(subpath=subpath)

    # relabel internal skin to some other label
    label_nifti = nibabel.load(subject_files.labeling)
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


def get_array_direction(electrode_pos, ellipsoid):
    """
    Determine electrode array direction given in ellipsoidal coordinates [theta, phi, alpha] and return direction
    angle alpha in Jacobian coordinates.

    Parameters
    ----------
    electrode_pos : np.ndarray of float [3]
        Spherical coordinates (theta, phi) and orientation angle (alpha) for electrode array (in this order)
    ellipsoid : Ellipsoid class instance
        Ellipsoid

    Returns
    -------
    alpha_jac : float
        Angle of array with respect to constant lambda (Jacobian coordinates)
    """
    # create cartesian vector in direction of electrode array in ellipsoidal coordinates (c1 -> c4)
    c1 = ellipsoid.ellipsoid2cartesian(coords=electrode_pos[:2])
    c2 = ellipsoid.ellipsoid2cartesian(coords=np.array([electrode_pos[0] - 1e-3, electrode_pos[1]]))
    c3 = ellipsoid.ellipsoid2cartesian(coords=np.array([electrode_pos[0], electrode_pos[1] - 1e-3]))
    c4 = c1 + (1e-3 * ((c3 - c1) * np.sin(electrode_pos[2]) + (c2 - c1) * np.cos(electrode_pos[2])))

    # create vector in direction of constant lambda
    l1 = ellipsoid.cartesian2jacobi(coords=c1)
    l2 = ellipsoid.jacobi2cartesian(coords=np.array([l1[0, 0] + 1e-3, l1[0, 1]]))

    l2c1 = ((l2 - c1) / np.linalg.norm(l2 - c1)).flatten()
    c4c1 = ((c4 - c1) / np.linalg.norm(c4 - c1)).flatten()
    angle_jac = np.arccos(np.dot((c4c1), (l2c1)))

    return angle_jac


def setup_logger(logname, filemode='w', format='[ %(name)s ] %(levelname)s: %(message)s',
                 datefmt='%H:%M:%S'):
    """
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


def plot_roi_field(e, roi, fn_out, e_label=None):
    """
    Creates plot files for paraview to visualize electric fields.
    Creates file triple of *_geo.hdf5,  *_data.hdf5, and *_data.xdmf, which can be loaded with Paraview.

    Parameters
    ----------
    e : np.ndarray of float [n_roi_ele, ] or list of np.ndarray of float [n_e-fields]
        Electric fields to visualize. Multiple fields can be passed in a list. Each field has to match the ROI.
    roi : RegionOfInterest class instance
        RegionOfInterest Object the data is associated with.
    fn_out : str
        Prefix of output file name, will append *_geo.hdf5, *_data.hdf5, and *_data.xdmf
    e_label : str or list of str [n_e-fields]
        Data names
    """
    if type(e) is not list:
        e = [e]

    if e_label is None:
        e_label = [f"data_{str(i)}" for i in range(len(e))]

    if type(e_label) is not list:
        e_label = [e_label]

    # if we have a connectivity (surface):
    if roi.con is not None:
        import pynibs

        fn_geo = fn_out + "_geo.hdf5"
        fn_data = fn_out + "_data.hdf5"

        # surface plot
        if roi.con.shape[1] == 3:
            pynibs.write_geo_hdf5_surf(out_fn=fn_geo, points=roi.nodes, con=roi.con, replace=True, hdf5_path='/mesh')
            pynibs.write_data_hdf5_surf(data=e, data_names=e_label, data_hdf_fn_out=fn_data, geo_hdf_fn=fn_geo,
                                        replace=True)

        # volume plot
        else:
            pynibs.write_geo_hdf5_vol(out_fn=fn_geo, points=roi.nodes, con=roi.con, replace=True, hdf5_path='/mesh')
            pynibs.write_data_hdf5_vol(data=e, data_names=e_label, data_hdf_fn_out=fn_data, geo_hdf_fn=fn_geo,
                                       replace=True)

    else:
        # if we just have points and data w/o connectivity information:
        e = np.hstack(e)
        np.savetxt(fn_out + "_data.txt", np.hstack((roi.center, e)))

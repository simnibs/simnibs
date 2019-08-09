import os
from collections import OrderedDict
import time
import glob
import functools

import numpy as np
import h5py

from . import optimize_tms
from ..simulation import cond
from ..simulation import fem
from ..simulation.sim_struct import SESSION, TMSLIST, SimuList, save_matlab_sim_struct
from ..msh import mesh_io
from ..utils.simnibs_logger import logger
from ..utils.file_finder import SubjectFiles
from ..utils.matlab_read import try_to_read_matlab_field, remove_None


class TMSOptimize():
    '''
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
    pos(optional): (3x1) array or None
        Position in scalp to use as a reference for the search space. By dafault, will
        project the target to the scalp and use it as the "pos" variable
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
        Spatial resolution of search area, in mm. Default: 2
    search_angles (optional): (2x1)float
        Limit angles to use in search, in degrees. If pos_ydir is set, will default to +/- 60 degrees
        around the original y axis. if not, will search in 360 degrees
    angle_resolution (optional): float
        Resolution to use for the angles, in degrees. Default: 20
    open_in_gmsh(optional): bool
        Wether to open the results in gmsh
    '''
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
        self.name = ''  # This is here only for leagacy reasons, it doesnt do anything
        # Optimization stuff
        self.target = None
        self.tissues = [2]
        self.target_size = 5
        self.pos = []
        self.pos_ydir = []
        self.distance = 4.
        self.didt = 1e6
        self.search_radius = 20
        self.spatial_resolution = 2
        self.search_angles = []
        self.angle_resolution = 20

        self.open_in_gmsh = True

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
        SESSION._set_logger(self)

    def _finish_logger(self):
        SESSION._finish_logger(self)

    def _prepare(self):
        """Prepares Leadfield for simulations
        relative paths are made absolute,
        empty fields are set to default values,
        check if required fields exist
        """
        self.fnamehead = os.path.abspath(os.path.expanduser(self.fnamehead))
        if not os.path.isfile(self.fnamehead):
            raise IOError('Cannot locate head mesh file: %s' % self.fnamehead)

        sub_files = SubjectFiles(self.fnamehead, self.subpath)
        self.fnamehead = sub_files.fnamehead
        self.subpath = sub_files.subpath

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
        assert self.target is not None and len(self.target) == 3
        assert self.search_radius > 0
        assert self.spatial_resolution > 0
        assert self.angle_resolution > 0
        # fix a few variables
        if len(self.pos_ydir) == 0:
            self.pos_ydir = None
        if len(self.pos) == 0:
            self.pos = np.copy(self.target)
        if len(self.search_angles) == 0:
            self.search_angles = None

    def sim_struct2mat(self):
        mat = SimuList.cond_mat_struct(self)
        mat['type'] = self.type
        mat['date'] = remove_None(self.date)
        mat['subpath'] = remove_None(self.subpath)
        mat['fnamehead'] = remove_None(self.fnamehead)
        mat['pathfem'] = remove_None(self.pathfem)
        mat['fname_tensor'] = remove_None(self.fname_tensor)
        mat['tissues'] = remove_None(self.tissues)

        mat['target'] = remove_None(self.target)
        mat['target_size'] = remove_None(self.target_size)
        mat['pos'] = remove_None(self.pos)
        mat['pos_ydir'] = remove_None(self.pos_ydir)
        mat['distance'] = remove_None(self.distance)
        mat['didt'] = remove_None(self.didt)
        mat['search_radius'] = remove_None(self.search_radius)
        mat['spatial_resolution'] = remove_None(self.spatial_resolution)
        mat['search_angles'] = remove_None(self.search_angles)
        mat['angle_resolution'] = remove_None(self.angle_resolution)
        mat['open_in_gmsh'] = remove_None(self.open_in_gmsh)

        return mat

    def read_mat_struct(self, mat):
        """ Reads form matlab structure
        Parameters
        ------------------
        mat: scipy.io.loadmat
            Loaded matlab structure
        """
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
        self.fname_tensor = try_to_read_matlab_field(
            mat, 'fname_tensor', str, self.fname_tensor
        )
        self.target = try_to_read_matlab_field(
            mat, 'target', list, self.target
        )
        self.target_size = try_to_read_matlab_field(
            mat, 'target_size', float, self.target_size
        )
        self.pos = try_to_read_matlab_field(
            mat, 'pos', list, self.pos
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
        self.search_angles = try_to_read_matlab_field(
            mat, 'search_angles', list, self.search_angles
        )
        self.angle_resolution = try_to_read_matlab_field(
            mat, 'angle_resolution', float, self.angle_resolution
        )
        self.open_in_gmsh = try_to_read_matlab_field(
            mat, 'open_in_gmsh', bool, self.open_in_gmsh
        )

    def run(self, cpus=1, allow_multiple_runs=False, save_mat=True):
        ''' Runs the tms optimization

        Parameters
        -----------
        cpus: int (optional)
            Number of cpus to use. Not nescessaraly will use all cpus. Default: 1
        allow_multiple_runs: bool (optinal)
            Wether to allow multiple runs in one folder. Default: False
        save_mat: bool (optional)
            Whether to save the ".mat" file of this structure
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

        self._prepare()
        if save_mat:
            save_matlab_sim_struct(
                self,
                os.path.join(
                    dir_name,
                    'simnibs_simulation_{0}.mat'.format(self.time_str)))
        logger.info(str(self))
        pos_matrices = optimize_tms.get_opt_grid(
            self.mesh, self.pos,
            handle_direction_ref=self.pos_ydir,
            distance=self.distance, radius=self.search_radius,
            resolution_pos=self.spatial_resolution,
            resolution_angle=self.angle_resolution,
            angle_limits=self.search_angles
        )
        cond = SimuList.cond2elmdata(self)
        didt_list = [self.didt for i in pos_matrices]

        # Define target region
        target_region = optimize_tms.define_target_region(
            self.mesh,
            self.target,
            self.target_size,
            self.tissues
        )
        if len(target_region) == 0:
            raise ValueError('Did not find any elements within the defined target region')

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
        # Write out the grid
        optimize_tms.plot_matsimnibs_list(
            pos_matrices,
            np.ones(len(pos_matrices)),
            "Grid",
            os.path.join(self.pathfem, 'coil_positions.geo')
        )
        v.add_merge(os.path.join(self.pathfem, 'coil_positions.geo'))
        v.add_view(VectorType=4, CenterGlyphs=0, Visible=1)
        v.write_opt(fn_target)
        if self.open_in_gmsh:
            mesh_io.open_in_gmsh(fn_target, True)

        # Run simulations
        fn_hdf5 = os.path.join(self.pathfem, self._name_hdf5())
        if os.path.isfile(fn_hdf5):
            os.remove(fn_hdf5)
        dataset = 'tms_optimization/E_norm'
        volumes = self.mesh.elements_volumes_and_areas()[target_region]
        # Define postporcessing to calculate average field norm
        def postprocessing(E, target_region, volumes):
            return np.average(
                np.linalg.norm(E[target_region - 1], axis=1),
                weights=volumes
            )
        postpro = functools.partial(
            postprocessing,
            target_region=target_region,
            volumes=volumes
        )
        fem.tms_many_simulations(
            self.mesh, cond,
            self.fnamecoil,
            pos_matrices, didt_list,
            fn_hdf5, dataset,
            post_pro=postpro,
            n_workers=cpus
        )
        # Read the field norms
        with h5py.File(fn_hdf5, 'a') as f:
            normE = f[dataset][:]
        os.remove(fn_hdf5)
        # Update the .geo file with the normE value
        optimize_tms.plot_matsimnibs_list(
            pos_matrices,
            normE,
            "normE at target",
            os.path.join(self.pathfem, 'coil_positions.geo')
        )
        v.add_merge(os.path.join(self.pathfem, 'coil_positions.geo'))
        v.add_view(VectorType=4, CenterGlyphs=0)
        v.write_opt(fn_target)

        # Run one extra simulation with the best position
        logger.info('Re-running best position')
        fn_out = fn_hdf5[:-5] + '.msh'
        fn_geo = fn_hdf5[:-5] + '_coil_pos.geo'
        fem.tms_coil(
            self.mesh, cond, self.fnamecoil, 'eEjJ',
            [pos_matrices[np.argmax(normE)]],
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
            visible_fields=['normE']
        )
        v.View[-1].CustomMax = 1
        v.View[-1].CustomMin = 0
        v.add_merge(fn_geo)
        v.write_opt(fn_out)
        if self.open_in_gmsh:
            mesh_io.open_in_gmsh(fn_out, True)
        # Another output with coil positions, possibly also transformed

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

    def __str__(self):
        string = 'Subject Folder: %s\n' % self.subpath
        string += 'Mesh file name: %s\n' % self.fnamehead
        string += 'Coil file: %s\n' % self.fnamecoil
        string += 'Target: %s\n' % self.target
        string += 'Reference position: %s\n' % self.pos
        string += 'Reference y: %s\n' % self.pos_ydir
        string += 'Coil distance: %s\n' % self.distance
        string += 'Search radius: %s\n' % self.search_radius
        string += 'Spatial resolution: %s\n' % self.spatial_resolution
        string += 'Angle resolution: %s' % self.angle_resolution
        return string

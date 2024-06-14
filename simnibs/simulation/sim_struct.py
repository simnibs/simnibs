# -*- coding: utf-8 -*-\
'''
    Structures for Simnibs.py
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.
    Copyright (C) 2018 Andre Antunes, Guilherme B Saturnino, Axel Thielscher

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
import os
import time
import copy
import glob
import gc
import logging
import functools
from typing import Union

import numpy as np
import scipy.io
from scipy.spatial.transform import Rotation as R
import nibabel
import h5py

import simnibs.utils.cond_utils
from simnibs.utils.cond_utils import COND
from simnibs.utils.mesh_element_properties import ElementTags

from ..mesh_tools import mesh_io, gmsh_view
from ..utils import transformations
from ..utils import simnibs_logger
from ..utils import file_finder
from ..utils.simnibs_logger import logger
from ..utils.file_finder import SubjectFiles
from ..utils.matlab_read import try_to_read_matlab_field, remove_None
from ..utils.csv_reader import read_csv_positions, _get_eeg_positions
from ..utils.transformations import project_points_on_surface
from . import fem
from . import electrode_placement
from .. import  __version__

class SESSION(object):
    """Class that defines a set of simnibs simulations

    Attributes
    ------------------------------
    date: str
        Date and time when the session struct was initiated
    poslist: list of sim_struct.SimuList()
        simulation set-up
    subpath: str
        path to m2m folder
    fnamehead: str
        file name of mesh
    pathfem: str
        path where the simulation should be saved
    fname_tensor: str
        name of DTI tensor file
    open_in_gmsh: bool
        Whether to open simulation in GMSH after finished
    map_to_surf: bool
        Whether to map fields to middle GM surface
    map_to_fsavg: bool
        Whether to map fields to FSAverage template
    map_to_vol: bool
        Whether to map fields to a NifTI volume
    map_to_MNI: bool
        Whether to map fields to the MNI template
    tissues_in_niftis: str or list of tissue tags, optional
        determines for which tissues the fields will be written
        to the NifTI volumes (for map_to_vol and map_to_MNI);
        either 'all' or a list of tags (standard: 2 for GM)
    fiducials: FIDUCIALS
        Structure with fiducil points
    fields: str
        Fields to be calculated for the simulations
    eeg_cap: str
        Name of eeg cap (in subject space)

    Parameters
    ------------------------
    matlab_struct: (optional) scipy.io.loadmat()
        matlab structure

    """

    def __init__(self, matlab_struct=None):
        # : Date when the session was initiated
        self.date = time.strftime("%Y-%m-%d %H:%M:%S")
        self.time_str = time.strftime("%Y%m%d-%H%M%S")
        # self.volfn = ''
        # self.vol = VOLUME()
        self.poslists = []
        self.fnamehead = None
        self.subpath = None
        self.pathfem = None
        self.fname_tensor = None
        self.open_in_gmsh = True
        self.map_to_surf = False
        self.map_to_fsavg = False
        self.map_to_vol = False
        self.map_to_MNI = False
        self.tissues_in_niftis = None
        self.fiducials = FIDUCIALS()
        self.fields = 'eE'
        self.eeg_cap = None
        self._prepared = False
        self._log_handlers = []

        if matlab_struct:
            self.read_mat_struct(matlab_struct)

    @property
    def type(self):
        return self.__class__.__name__

    def add_poslist(self, pl):
        """ Adds a SimList object to the poslist variable

        Parameters:
        ----------------
        pl: sim_struct.SimuList
            SimuList object
        """
        if not isinstance(pl, SimuList):
            raise TypeError(
                'Elements in poslist must be subclasses of SimuList')
        self.poslists.append(pl)

    def remove_poslist(self, number):
        """Removes the specified poslist

        Parameters:
        ------------------------------
        number: int
        indice of postist to be removed
        """
        del self.poslists[number]

    def clear_poslist(self):
        """ Removes all poslists
        """
        self.poslists = []

    def _prepare(self):
        """Prepares session for simulations
        relative paths are made absolute,
        empty fields are set to default values,
        check if required fields exist
        """
        # if self._prepared:
        #    raise Exception('Re-using a Python SESSION '
        #                    'structure can cause bugs!'
        #                    ' Please initialize a new SESSION')

        if self.fnamehead:
            self.fnamehead = os.path.abspath(
                os.path.expanduser(self.fnamehead))
            if not os.path.isfile(self.fnamehead):
                raise IOError('Cannot locate head mesh file: %s' %
                              self.fnamehead)
        else:
            self.subpath = os.path.abspath(os.path.expanduser(self.subpath))
            if not os.path.isdir(self.subpath):
                raise IOError(
                    'Cannot locate subjects m2m folder: %s' % self.subpath)

        sub_files = SubjectFiles(self.fnamehead, self.subpath)
        self.fnamehead = sub_files.fnamehead
        self.subpath = sub_files.subpath

        if not os.path.isdir(self.subpath):
            logger.warning('Cannot locate subjects m2m folder')
            logger.warning('some postprocessing options might fail')
            self.subpath = None

        if not self.fname_tensor:
            self.fname_tensor = sub_files.tensor_file

        if not self.eeg_cap:
            self.eeg_cap = sub_files.eeg_cap_1010

        if not self.tissues_in_niftis:
            self.tissues_in_niftis = [2]

        logger.info('Head Mesh: {0}'.format(self.fnamehead))
        logger.info('Subject Path: {0}'.format(self.subpath))
        self.pathfem = os.path.abspath(os.path.expanduser(self.pathfem))
        logger.info('Simulation Folder: {0}'.format(self.pathfem))

        if os.path.isfile(self.fnamehead):
            mesh = mesh_io.read_msh(self.fnamehead)
            mesh.fix_surface_labels()
            for PL in self.poslists:
                PL.postprocess = self.fields
                PL.fn_tensor_nifti = self.fname_tensor
                PL.eeg_cap = self.eeg_cap
                PL._prepare()
                if not PL.mesh:
                    PL.mesh = mesh
        else:
            raise IOError(
                'Could not find head mesh file: {0}'.format(self.fnamehead))

        self._prepared = True

    def run(self, cpus=1, allow_multiple_runs=False, save_mat=True):
        """ Run simulations in the current session

        Parameters
        -----------
        cpus: int (optional)
            Number of cpus to use. Not necessarily will use all cpus. Default: 1
        allow_multiple_runs: bool (optinal)
            Whether to allow multiple runs in one folder. Default: False
        save_mat: bool (optional)
            Whether to save the ".mat" file of this structure

        Returns
        ---------
        Writes the simulations

        """
        self._set_logger()
        self._prepare()
        dir_name = os.path.abspath(os.path.expanduser(self.pathfem))
        final_names = []
        final_names_geo = []

        if os.path.isdir(dir_name):
            g = glob.glob(os.path.join(dir_name, 'simnibs_simulation*.mat'))
            if len(g) > 0 and not allow_multiple_runs:
                raise IOError('Found already existing simulation results in directory.'
                              ' Please run the simulation in a new directory or delete'
                              ' the simnibs_simulation*.mat files from the folder : {0}'.format(dir_name))
            logger.info(
                'Running simulations in the directory: {0}'.format(dir_name))
        else:
            logger.info(
                'Running simulations on new directory: {0}'.format(dir_name))
            os.makedirs(dir_name)

        if save_mat:
            save_matlab_sim_struct(
                self,
                os.path.join(
                    dir_name,
                    'simnibs_simulation_{0}.mat'.format(self.time_str)))

        name = os.path.split(self.fnamehead)[1]
        name = os.path.splitext(name)[0]
        for i, PL in enumerate(self.poslists):
            logger.info('Running Poslist Number: {0}'.format(i + 1))
            if PL.name:
                simu_name = os.path.join(dir_name, PL.name)
            else:
                if PL.type == 'TMSLIST':
                    simu_name = os.path.join(
                        dir_name, '{0}_TMS_{1}'.format(name, i + 1))
                elif PL.type == 'TDCSLIST':
                    simu_name = os.path.join(
                        dir_name, '{0}_TDCS_{1}'.format(name, i + 1))
                else:
                    simu_name = os.path.join(dir_name, '{0}'.format(i + 1))
            fn, fn_geo = PL.run_simulation(
                simu_name, cpus=cpus, view=self.open_in_gmsh)
            PL.mesh = None
            final_names += fn
            final_names_geo += fn_geo
            logger.info('Finished Running Poslist Number: {0}'.format(i + 1))
            logger.info('Result Files:\n{0}'.format('\n'.join(fn)))
            gc.collect()

        folders = []
        if self.map_to_surf or self.map_to_fsavg or self.map_to_vol or self.map_to_MNI:
            if not self.subpath:
                raise IOError('Cannot run postprocessing: subpath not set')
            elif not os.path.isdir(self.subpath):
                raise IOError('Cannot run postprocessing: {0} is not a '
                              'directory'.format(self.subpath))

        if self.map_to_surf or self.map_to_fsavg:
            logger.info('Interpolating to the middle of Gray Matter')
            out_folder = os.path.join(dir_name, 'subject_overlays')
            folders += [out_folder]
            if self.map_to_fsavg:
                out_fsavg = os.path.join(dir_name, 'fsavg_overlays')
                folders += [out_fsavg]
            else:
                out_fsavg = None
            for f, f_geo in zip(final_names, final_names_geo):
                if f.endswith('.msh'):
                    transformations.middle_gm_interpolation(
                        f, self.subpath, out_folder,
                        out_fsaverage=out_fsavg, depth=0.5,
                        open_in_gmsh=self.open_in_gmsh, f_geo=f_geo)

        keep_tissues = None
        if type(self.tissues_in_niftis) == list:
            keep_tissues = self.tissues_in_niftis

        if self.map_to_vol:
            logger.info('Mapping to volume')
            out_folder = os.path.join(dir_name, 'subject_volumes')
            folders += [out_folder]
            if not os.path.isdir(out_folder):
                os.mkdir(out_folder)
            for f in final_names:
                if f.endswith('.msh'):
                    name = os.path.split(f)[1]
                    name = os.path.splitext(name)[0] + '.nii.gz'
                    name = os.path.join(out_folder, name)
                    transformations.interpolate_to_volume(
                        f, self.subpath, name, keep_tissues=keep_tissues)

        if self.map_to_MNI:
            logger.info('Mapping to MNI space')
            out_folder = os.path.join(dir_name, 'mni_volumes')
            folders += [out_folder]
            if not os.path.isdir(out_folder):
                os.mkdir(out_folder)

            for f in final_names:
                if f.endswith('.msh'):
                    name = os.path.split(f)[1]
                    name = os.path.splitext(name)[0] + '.nii.gz'
                    name = os.path.join(out_folder, name)
                    transformations.warp_volume(
                        f, self.subpath, name, keep_tissues=keep_tissues)

        logger.info('=====================================')
        logger.info('SimNIBS finished running simulations')
        logger.info('Simulation Result Meshes:')
        [logger.info(f) for f in final_names]
        if len(folders) > 0:
            logger.info('Postprocessing Folders:')
            [logger.info(f) for f in folders]
        logger.info('=====================================')
        self._finish_logger()
        return final_names

    def read_mat_struct(self, mat):
        """ Reads form matlab structure
        Parameters
        ------------------
        mat: scipy.io.loadmat or str
            Loaded matlab structure or filename
        """
        if type(mat) == str:
            mat = scipy.io.loadmat(mat)

        self.date = try_to_read_matlab_field(mat, 'date', str, self.date)
        # self.volfn = try_to_read_matlab_field(mat, 'volfn', str, self.volfn)
        # self.vol = try_to_read_matlab_field(mat, 'volfn', VOLUME, VOLUME())
        self.subpath = try_to_read_matlab_field(mat, 'subpath', str,
                                                self.subpath)
        self.fnamehead = try_to_read_matlab_field(
            mat, 'fnamehead', str, self.fnamehead)
        self.pathfem = try_to_read_matlab_field(
            mat, 'pathfem', str, self.pathfem)
        self.fname_tensor = try_to_read_matlab_field(mat, 'fname_tensor', str,
                                                     self.fname_tensor)
        self.eeg_cap = try_to_read_matlab_field(
            mat, 'eeg_cap', str, self.eeg_cap)

        self.open_in_gmsh = try_to_read_matlab_field(
            mat, 'open_in_gmsh', bool, self.open_in_gmsh)
        self.map_to_vol = try_to_read_matlab_field(
            mat, 'map_to_vol', bool, self.map_to_vol)
        self.map_to_MNI = try_to_read_matlab_field(
            mat, 'map_to_MNI', bool, self.map_to_MNI)
        self.map_to_surf = try_to_read_matlab_field(
            mat, 'map_to_surf', bool, self.map_to_surf)
        self.map_to_fsavg = try_to_read_matlab_field(
            mat, 'map_to_fsavg', bool, self.map_to_fsavg)
        self.tissues_in_niftis = try_to_read_matlab_field(
            mat, 'tissues_in_niftis', list, self.tissues_in_niftis)
        if (self.tissues_in_niftis is not None
            and type(self.tissues_in_niftis[0]) == str
                and self.tissues_in_niftis[0] == 'a'):
            self.tissues_in_niftis = 'all'

        self.fields = try_to_read_matlab_field(
            mat, 'fields', str, self.fields)

        self.fiducials.read_mat_struct(mat)
        if len(mat['poslist']) > 0:
            for PL in mat['poslist'][0]:
                if PL['type'][0] == 'TMSLIST':
                    self.add_poslist(TMSLIST(PL[0][0]))
                elif PL['type'][0] == 'TDCSLIST':
                    self.add_poslist(TDCSLIST(PL[0][0]))
                else:
                    raise IOError(
                        "poslist type is not of type TMSLIST or TDCSLIST")

    def sim_struct2mat(self):
        """ Makes a dictionary for saving a matlab structure with scipy.io.savemat()

        Returns
        --------------------
        dict
            Dictionaty for usage with scipy.io.savemat
        """
        mat = {}
        mat['type'] = 'SESSION'
        mat['date'] = remove_None(self.date)
        # mat['volfn'] = remove_None(self.volfn)
        mat['subpath'] = remove_None(self.subpath)
        mat['eeg_cap'] = remove_None(self.eeg_cap)
        mat['fnamehead'] = remove_None(self.fnamehead)
        mat['pathfem'] = remove_None(self.pathfem)
        mat['fname_tensor'] = remove_None(self.fname_tensor)
        # mat['vol'] = self.vol.sim_struct2mat()
        mat['map_to_vol'] = remove_None(self.map_to_vol)
        mat['map_to_MNI'] = remove_None(self.map_to_MNI)
        mat['map_to_fsavg'] = remove_None(self.map_to_fsavg)
        mat['map_to_surf'] = remove_None(self.map_to_surf)
        mat['tissues_in_niftis'] = remove_None(self.tissues_in_niftis)
        mat['fields'] = remove_None(self.fields)
        mat['fiducials'] = self.fiducials.sim_struct2mat()
        mat['poslist'] = []
        for PL in self.poslists:
            mat['poslist'].append(PL.sim_struct2mat())
            # pass

        return mat

    def add_tdcslist(self, tdcslist=None):
        ''' Appends a TDCSLIST to the current SESSION

        Parameters:
        ------------
        tdcslist: TDCSLIST (optional)
            tdcslist to be added. (Default: empty TDCSLIST)

        Returns
        -------
        tdcslist: TDCSLIST
            the tdcslist added to this SESSION
        '''
        if tdcslist is None:
            tdcslist = TDCSLIST()

        self.poslists.append(tdcslist)
        return tdcslist

    def add_tmslist(self, tmslist=None):
        ''' Appends a TMSLIST to the current SESSION

        Parameters:
        ------------
        tmslist: TMSLIST (optional)
            tmslist to be added. (Default: empty TMSLIST)

        Returns
        -------
        tmslist: TMSLIST
            the tmslist added to this SESSION
        '''
        if tmslist is None:
            tmslist = TMSLIST()

        self.poslists.append(tmslist)
        return tmslist

    def _set_logger(self, fname_prefix='simnibs_simulation', summary=True):
        """
        Set-up loggger to write to a file

        Parameters
        ----------
        fname_prefix: str, optional
            Prefix of log-file. Defaults to 'simnibs_simulation'.
        summary: bool, optional
            Create summary file 'fields_summary.txt'. Default: True.
        """
        if not os.path.isdir(self.pathfem):
            os.mkdir(self.pathfem)
        log_fn = os.path.join(
            self.pathfem,
            fname_prefix + '_{0}.log'.format(self.time_str))
        fh = logging.FileHandler(log_fn, mode='w')
        formatter = logging.Formatter(
            f'[ %(name)s {__version__} - %(asctime)s - %(process)d ]%(levelname)s: %(message)s')
        fh.setFormatter(formatter)
        fh.setLevel(logging.DEBUG)
        logger = logging.getLogger("simnibs")
        logger.addHandler(fh)
        self._log_handlers += [fh]

        if summary:
            fn_summary = os.path.join(self.pathfem, 'fields_summary.txt')
            fh_s = logging.FileHandler(fn_summary, mode='w')
            fh_s.setFormatter(logging.Formatter('%(message)s'))
            fh_s.setLevel(25)
            logger.addHandler(fh_s)
            self._log_handlers += [fh_s]
        simnibs_logger.register_excepthook(logger)

    def _finish_logger(self):
        logger = logging.getLogger("simnibs")
        [logger.removeHandler(lh) for lh in self._log_handlers]
        self._log_handlers = []
        simnibs_logger.unregister_excepthook()

    def __str__(self):
        string = 'Subject Folder: %s\n' % self.subpath
        string += 'Mesh file name: %s\n' % self.fnamehead
        string += 'Date: %s\n' % self.date
        string += 'Number of Poslists:%d' % len(self.poslists)
        return string

    def __eq__(self, other):
        try:
            return self.__dict__ == other.__dict__
        except AttributeError:
            return False


class FIDUCIALS(object):
    def __init__(self, Nz=[], Iz=[], LPA=[], RPA=[]):
        self.Nz = Nz
        self.Iz = Iz
        self.LPA = LPA
        self.RPA = RPA

    def read_mat_struct(self, mat):
        """ Reads form matlab structure
        Parameters
        ------------------
        mat: scipy.io.loadmat
            Loaded matlab structure
        """
        self.Nz = try_to_read_matlab_field(mat, 'Nz', list, self.Nz)
        self.Iz = try_to_read_matlab_field(mat, 'Iz', list, self.Iz)
        self.LPA = try_to_read_matlab_field(mat, 'LPA', list, self.LPA)
        self.RPA = try_to_read_matlab_field(mat, 'RPA', list, self.RPA)

    def sim_struct2mat(self):
        """ Makes a dictionary for saving a matlab structure with scipy.io.savemat()

        Returns
        --------------------
        dict
            Dictionaty for usage with scipy.io.savemat
        """
        mat = {}
        mat['type'] = 'FIDUCIALS'
        mat['Nz'] = remove_None(self.Nz)
        mat['Iz'] = remove_None(self.Iz)
        mat['LPA'] = remove_None(self.LPA)
        mat['RPA'] = remove_None(self.RPA)
        return mat

    def from_csv(self, fn_csv):
        ''' Reads a csv file and loads the fiducials defined in it

        ----------
        fn_csv: str
            CSV file with the fields
                 Type, pos_x, pos_y, pos_z, name, whatever
            Type must be Fiducial, and name Nz, Iz, LPA, RPA
        '''
        type_, coordinates, _, name, _, _ = read_csv_positions(fn_csv)
        for t, c, n in zip(type_, coordinates, name):
            if t == 'Fiducial':
                if n in ['Nz', 'Iz', 'LPA', 'RPA']:
                    self.__dict__[n] = c
                else:
                    logger.warning(
                        'Unrecognized Fiducial: {0} '
                        'Acceptable fiducuals are: {1}'
                        .format(n, ['Nz', 'Iz', 'LPA', 'RPA']))


class SimuList(object):
    """ Parent class

    Attributes:
    ----------------------------
    cond: list
        list of COND structures with conductivity information
    mesh: simnibs.mesh_io.Msh
        Mesh where the simulations will be performed
    anisotropy_type: property, can be 'scalar', 'vn' or 'mc'
        type of anisotropy for simulation
    fn_tensor_nifti: str
        file name of nifti with tensor information
    postprocess: property
        fields to be calculated. valid fields are: 'v' , 'E', 'e', 'J', 'j', 'g', 's', 'D', 'q'
    anisotropy_vol: ndarray (optional)
        Volume with anisotropy information (lower priority over fn_tensor_nifti)
    anisotropy_affine: ndarray (optional)
        4x4 affine matrix describing the transformation from the regular grid to the mesh
        space (lower priority over fn_tensor_nifti)
    anisotropic_tissues: list (optional)
        List with tissues with anisotropic conductivities
    eeg_cap: str (optional)
        Name of csv file with EEG positions
    aniso_maxratio: float (optional)
        Maximum ration between the largest and smallest eigenvalue in a conductivity
        tensor.
    aniso_cond: float (optional)
        Maximum eigenvalue of a conductivity tensor.
    solver_options: str (optional)
        Options for the FEM solver
    """

    def __init__(self, mesh=None):
        # list of conductivities (using COND class)
        self.cond = simnibs.utils.cond_utils.standard_cond()
        self.mesh = mesh  # The mesh where the simulation will be performed
        self.fn_tensor_nifti = None  # File name with anisotropy information
        # The 2 variables bellow are set when the _get_vol() method is called
        # If set, they have priority over fn_tensor_nifti
        self.anisotropy_vol = None  # 4-d data with anisotropy information
        self.anisotropy_affine = None  # 4x4 affine transformation from the regular grid
        # if an anisotropic conductivity is to be used,
        self.anisotropic_tissues = [1, 2]
        self.suppl = []
        self.name = None  # Name to be given by simulations
        self.eeg_cap = None
        self.aniso_maxratio = 10
        self.aniso_maxcond = 2
        self.solver_options = None
        self._anisotropy_type = 'scalar'
        self._postprocess = ['e', 'E', 'j', 'J']

    @property
    def type(self):
        return self.__class__.__name__

    @property
    def anisotropy_type(self):
        return self._anisotropy_type

    @anisotropy_type.setter
    def anisotropy_type(self, value):
        if value not in ['scalar', 'dir', 'vn', 'mc']:
            raise ValueError("Invalid anisotroy type: {0}".format(value))
        else:
            self._anisotropy_type = value

    @property
    def postprocess(self):
        return self._postprocess

    @postprocess.setter
    def postprocess(self, value):
        valid_postptocessing = set(
            ['v', 'E', 'e', 'J', 'j', 'g', 's', 'D'])
        if not set(value).issubset(valid_postptocessing):
            raise ValueError('postprocessing operation: {0} \n'
                             'Not Supported, Supported operations are: '
                             '{1}'.format(list(set(value) - valid_postptocessing),
                                          list(valid_postptocessing)))
        else:
            self._postprocess = value

    @property
    def conductivity(self):
        c = copy.copy(self.cond)
        c.insert(0, None)
        return c

    def check_conductivities(self):
        if not self.anisotropy_type:
            self.anisotropy_type = 'scalar'

        if self.anisotropy_type not in ['vn', 'mc', 'dir', 'scalar']:
            raise ValueError(
                'Unrecognized anisotropy type:{0}'.format(self.anisotropy_type))

        if self.anisotropy_type != 'scalar' and not self.fn_tensor_nifti:
            raise IOError(
                'Cannot perform simulation with anisotropic conductivity. '
                'Invalid tensor file')

    def cond_mat_struct(self):
        """Returns a mat structure for the conductivity list

        Returns
        --------------------------
        dict
            Dictionary for scipy.io.savemat
        """
        mat_cond = {}

        # cond data type
        cond_dt = np.dtype([('type', 'O'),
                            ('name', 'O'), ('value', 'O'), ('descrip', 'O'),
                            ('distribution_type', 'O'),
                            ('distribution_parameters', 'O')])
        cond_mat = np.empty((0,), dtype=cond_dt)

        for c in self.cond:
            c = copy.deepcopy(c)
            if not c.value:  # for data type reasons
                c.value = []
            if not c.name:  # for data type reasons
                c.name = []
            if not c.distribution_type:
                c._distribution_type = []
            cond_mat = np.append(cond_mat, np.array(
                [('COND', c.name, c.value, c.descrip,
                  c.distribution_type, c.distribution_parameters)],
                dtype=cond_dt))

        mat_cond['cond'] = cond_mat
        mat_cond['anisotropy_type'] = remove_None(self.anisotropy_type)
        mat_cond['fn_tensor_nifti'] = remove_None(self.fn_tensor_nifti)
        mat_cond['aniso_maxratio'] = remove_None(self.aniso_maxratio)
        mat_cond['aniso_maxcond'] = remove_None(self.aniso_maxcond)
        # Not really related to the conductivity
        mat_cond['name'] = remove_None(self.name)
        mat_cond['solver_options'] = remove_None(self.name)

        return mat_cond

    def read_cond_mat_struct(self, mat_struct):
        """ Reads the conductivity part of the matlab structure
        sets self.cond and set.anisotropy_type

        Parameters
        ------------------------
        mat_struct: matlab structure, as loaded by scipy.io.loadmat()

        """
        self.cond = []
        try:
            for c in mat_struct['cond'][0]:
                self.cond.append(COND(c))
        except (KeyError, ValueError):
            pass

        self.anisotropy_type = try_to_read_matlab_field(
            mat_struct, 'anisotropy_type', str, self.anisotropy_type)
        self.fn_tensor_nifti = try_to_read_matlab_field(
            mat_struct, 'fn_tensor_nifti', str, self.fn_tensor_nifti)
        self.name = try_to_read_matlab_field(
            mat_struct, 'name', str, self.name)
        self.aniso_maxcond = try_to_read_matlab_field(
            mat_struct, 'aniso_maxcond', float, self.aniso_maxcond)
        self.aniso_maxratio = try_to_read_matlab_field(
            mat_struct, 'aniso_maxratio', float, self.aniso_maxratio)
        self.solver_options = try_to_read_matlab_field(
            mat_struct, 'solver_options', str, self.solver_options)

    def compare_conductivities(self, other):
        if self.anisotropy_type != other.anisotropy_type:
            return False

        if len(self.cond) != len(other.cond):
            return False

        for ii in range(len(self.cond)):
            if not self.cond[ii] == other.cond[ii]:
                return False
        return True

    def __eq__(self, other):
        if not isinstance(other, SimuList):
            return False

        return self.compare_conductivities(other)

    def cond2elmdata(self, mesh=None, excentricity_scale=None, logger_level=20):
        ''' Transforms a conductivity list to an ElementData field

        Parameters
        -----------
        mesh: simnibs.mesh_io.Msh (optional
            Mesh where the conductivities will be applied. Default: self.mesh

        excentricity_scale: float (optional)
            Scales the excentricity of conductivity tensors. Used in gPC simulations. Default: do not scale
            excentricities

        Returns
        --------
        conductivity: mesh_io.msh.ElementData()
            ElementData structure with conductivity information for each tetrahedron
        '''
        if mesh is None:
            mesh = self.mesh

        if mesh is None:
            raise ValueError('The mesh for this simulation is not set')

        cond_list = [c.value for c in self.cond]
        level = logger_level

        if self.anisotropy_type == 'scalar':
            logger.log(level, 'Using isotropic conductivities')
            return simnibs.utils.cond_utils.cond2elmdata(mesh, cond_list)

        elif self.anisotropy_type == 'dir':
            logger.log(
                level,
                'Using anisotropic direct conductivities based on the file:'
                ' {0}'.format(self.fn_tensor_nifti))
            image, affine = self._get_vol_info()
            return simnibs.utils.cond_utils.cond2elmdata(mesh, cond_list,
                                     anisotropy_volume=image,
                                     affine=affine,
                                     aniso_tissues=self.anisotropic_tissues,
                                     max_cond=self.aniso_maxcond,
                                     max_ratio=self.aniso_maxratio,
                                     excentricity_scaling=excentricity_scale)

        elif self.anisotropy_type == 'vn':
            logger.log(level, 'Using anisotropic volume normalized conductivities based on the file:'
                              ' {0}'.format(self.fn_tensor_nifti))
            image, affine = self._get_vol_info()
            return simnibs.utils.cond_utils.cond2elmdata(mesh, cond_list,
                                     anisotropy_volume=image,
                                     affine=affine,
                                     aniso_tissues=self.anisotropic_tissues,
                                     normalize=True,
                                     max_cond=self.aniso_maxcond,
                                     max_ratio=self.aniso_maxratio,
                                     excentricity_scaling=excentricity_scale)

        elif self.anisotropy_type == 'mc':
            logger.log(level, 'Using isotropic mean conductivities based on the file:'
                              ' {0}'.format(self.fn_tensor_nifti))
            image, affine = self._get_vol_info()
            return simnibs.utils.cond_utils.cond2elmdata(mesh, cond_list,
                                     anisotropy_volume=image,
                                     affine=affine,
                                     aniso_tissues=self.anisotropic_tissues,
                                     max_cond=self.aniso_maxcond,
                                     max_ratio=self.aniso_maxratio,
                                     excentricity_scaling=0.)

        else:
            raise ValueError('Invalid anisotropy_type: {0}'
                             ' valid types are: "scalar", "mc", "dir", "vn" '
                             ''.format(self.anisotropy_type))

    def _get_vol_info(self):
        if self.anisotropy_vol is not None:
            if self.anisotropy_affine is not None:
                return self.anisotropy_vol, self.anisotropy_affine

        if not self.fn_tensor_nifti:
            raise ValueError('could not get anisotropy information: '
                             'fn_tensor_nifti not set')

        fn_nifti = \
            os.path.abspath(os.path.expanduser(self.fn_tensor_nifti))
        if not os.path.isfile(fn_nifti):
            raise ValueError(
                'Could not find file \'{0}\' to get anisotropy '
                'information'.format(self.fn_tensor_nifti))

        # Load the nifti and interpolate the conductivieis
        image = nibabel.load(fn_nifti)
        affine = image.affine
        return image.dataobj, affine

    def _write_conductivity_to_hdf5(self, fn_hdf5, path='cond/'):
        """
        Parameters
        -----------
        fn_hdf5: str
            file name of hdf5 file
        path: str
            path in the hdf5 file where the conductivity information should be saved
        """
        with h5py.File(fn_hdf5, 'a') as f:
            try:
                g = f.create_group(path)
            except ValueError:
                g = f[path]
            value_array = np.array([c.value for c in self.cond], dtype=float)
            g.create_dataset('values', data=value_array)
            g.create_dataset('names',
                             data=np.array([c.name for c in self.cond],
                                           dtype=np.string_))
            g.create_dataset('distribution_types',
                             data=np.array(
                                 [c.distribution_type for c in self.cond],
                                 dtype=np.string_))

            distribution_parameters = np.nan * \
                np.zeros((len(self.cond), 4), dtype=float)
            for i, c in enumerate(self.cond):
                distribution_parameters[i, :len(
                    c.distribution_parameters)] = c.distribution_parameters
            g.create_dataset('distribution_parameters',
                             data=distribution_parameters)
            g.attrs['anisotropic_tissues'] = self.anisotropic_tissues
            g.attrs['anisotropy_type'] = self.anisotropy_type
            aniso = False
            if self.anisotropy_type in ['dir', 'vn', 'mc']:
                image, affine = self._get_vol_info()
                g.create_dataset('anisotropy_vol', data=image)
                g.create_dataset('anisotropy_affine', data=affine)

    def _get_conductivity_from_hdf5(self, fn_hdf5, path='cond/'):
        '''get conductivity information from HDF5 file '''
        with h5py.File(fn_hdf5, 'r') as f:
            try:
                g = f[path]
            except:
                raise IOError('Could not find the group {0} '
                              'in the HDF5 file'.format(path))
            self.anisotropy_type = g.attrs['anisotropy_type']
            self.anisotropic_tissues = g.attrs['anisotropic_tissues']
            value_array = g['values'][:]
            name_array = g['names'][:]
            dist_array = g['distribution_types'][:]
            dist_p_array = g['distribution_parameters'][:]
            self.cond = [COND()
                         for i in range(max(len(value_array), len(name_array)))]
            for i, c in enumerate(self.cond):
                if not np.isnan(value_array[i]):
                    self.cond[i].value = value_array[i]
                if name_array[i] != b'None':
                    self.cond[i].name = name_array[i].decode()
                if dist_array[i] != b'None':
                    self.cond[i].distribution_type = dist_array[i].decode()
                self.cond[i].distribution_parameters = \
                    dist_p_array[i][~np.isnan(dist_p_array[i])].tolist()

            try:
                self.anisotropy_affine = g['anisotropy_affine'][:]
            except KeyError:
                self.anisotropy_affine = None
            try:
                self.anisotropy_vol = g['anisotropy_vol'][:]
            except KeyError:
                self.anisotropy_vol = None

    def _scalp_geo(self, m, fn_out, scalp_idx : int = ElementTags.SCALP_TH_SURFACE):
        ''' write out scalp surface as geo file '''
        idx = (m.elm.tag1 == scalp_idx) & (m.elm.elm_type == 2)
        mesh_io.write_geo_triangles(m.elm[idx, :3]-1,
                                    m.nodes.node_coord, fn_out,
                                    name='scalp', mode='ba')


class TMSLIST(SimuList):
    """List of TMS coil position

    Note: Children of SimuList class

    Parameters
    ----------
    matlab_struct(optional): scipy.io.loadmat structure
        matlab structure defining the posist

    Attributes
    ----------
    fnamecoil: str
        Relative path of coil file.
    pos: list of simnibs.simulation.sim_struct.POSITION() structures
        Definition of coil positions.
    """

    def __init__(self, matlab_struct=None):
        SimuList.__init__(self)
        self.fnamecoil = ''
        self.pos = []
        self.postprocess = ['E', 'e', 'J', 'j']

        if matlab_struct is not None:
            self.read_mat_struct(matlab_struct)

    def _prepare(self):
        """Prepares structures for simulations
        Changes anisotropy_type and fnamecoil, _prepares pos structures
        """
        self.check_conductivities()
        self.resolve_fnamecoil()
        for pos in self.pos:
            pos.eeg_cap = self.eeg_cap
            pos._prepare()

    def read_mat_struct(self, PL):
        """ Reads matlab poslist structure

        Changes all the fields in poslist, as well as poscoil

        Parameters:
        ------------------------------
        PL: scipy.io.loadmat object or str
            Output of scipy.io.loadmat() or path to scipy.io.savemat() file.
        """
        if type(PL) == str:
            PL = scipy.io.loadmat(PL)
        self.read_cond_mat_struct(PL)

        try:
            self.fnamecoil = str(PL['fnamecoil'][0])
        except:
            self.fnamecoil = ''

        if len(PL['pos']) > 0:
            for pos in PL['pos'][0]:
                self.pos.append(POSITION(pos))

        try:
            self.suppl = PL['suppl'][0]
        except:
            pass

    def sim_struct2mat(self):
        """ Dictionaty for saving as a matlab structure with scipy.io.savemat

        Returns
        ----------------------------
        dict
            Dictionary with poslist parameters for saving into a matlab structure

        """
        mat_poslist = self.cond_mat_struct()

        mat_poslist['type'] = 'TMSLIST'
        mat_poslist['fnamecoil'] = remove_None(self.fnamecoil)
        mat_poslist['suppl'] = remove_None(self.suppl)

        # pos data type
        pos_dt = np.dtype([('type', 'O'),
                           ('name', 'O'), ('date', 'O'),
                           ('matsimnibs', 'O'), ('didt', 'O'),
                           ('fnamefem', 'O'), ('centre', 'O'),
                           ('pos_ydir', 'O'), ('distance', 'O')])

        pos_mat = np.empty((0,), dtype=pos_dt)
        for pos in self.pos:
            pos_array = np.array([('POSITION', remove_None(pos.name),
                                   remove_None(pos.date),
                                   remove_None(pos.matsimnibs),
                                   remove_None(pos.didt),
                                   remove_None(pos.fnamefem),
                                   remove_None(pos.centre),
                                   remove_None(pos.pos_ydir),
                                   remove_None(pos.distance))],
                                 dtype=pos_dt)
            pos_mat = np.append(pos_mat, pos_array)

        mat_poslist['pos'] = pos_mat

        return mat_poslist

    def add_positions_from_csv(self, fn_csv):
        ''' Reads a csv file and adds the positions defined to the tmslist

        Parameters
        ----------
        fn_csv: str
            CSV file with the fields
                Type, pos_x, pos_y, pos_z, ez_x, ez_y, ez_z, ey_x, ey_y, ey_z, dist, name, ...
                "Type" needs to be CoilPos. The positions are in subject space. The
                transformations module can transfrom from MNI to subject space
        '''
        type_, coordinates, extra, name, _, _ = read_csv_positions(fn_csv)
        for t, c, e, n in zip(type_, coordinates, extra, name):
            if t == 'CoilPos':
                p = POSITION()
                vz = e[:3]
                vy = e[3:6]
                vx = np.cross(vy, vz)
                mat = np.eye(4)
                mat[:3, 0] = vx
                mat[:3, 1] = vy
                mat[:3, 2] = vz
                mat[:3, 3] = c
                p.matsimnibs = mat.tolist()
                p.name = n
                self.pos.append(p)

    def resolve_fnamecoil(self):
        try:
            fnamecoil = os.path.expanduser(self.fnamecoil)
        except TypeError:
            print(self.fnamecoil)
            raise TypeError('fnamecoil is not a string')
        if os.path.isfile(fnamecoil):
            self.fnamecoil = fnamecoil
        else:
            fnamecoil = os.path.join(file_finder.coil_models, self.fnamecoil)
            if os.path.isfile(fnamecoil):
                self.fnamecoil = fnamecoil
            else:
                raise IOError(
                    'Could not find coil file: {0}'.format(self.fnamecoil))

    def run_simulation(self, fn_simu, cpus=1, view=True):
        """ Runs the TMS simulations. One per "position" object defined in the ".pos"
        field of the current structure.

        Parameters
        ----------
        fn_simu: str
            simulation path and name
        cpus: int (Optional)
            Number of parallel processes to run. Default: 1
        view: bool
            Whether to open the result in Gmsh
        Returns
        ---------
        final_name: list
          List with the names of the output files.
          We return a list for consistency with the TMS version
        geo_names: list
            List with the names of the geo-files that accompany the output files.
        """
        add_scalp_to_geo = True
        if len(self.pos) == 0:
            raise ValueError(
                'There are no positions defined for this poslist!')
        fn_simu = os.path.abspath(os.path.expanduser(fn_simu))
        assert isinstance(self.mesh, mesh_io.Msh), \
            'mesh property not set appropriately!'
        # gPC
        for c in self.cond:
            if c.distribution_type:
                logger.warning(
                    'Distribution value for conductivity found, starting gPC')
                fn = self.run_gpc(fn_simu, cpus=cpus)
                return fn, [None]*len(fn)

        logger.info('Began to run TMS simulations')
        logger.info(f'Coil file: {self.fnamecoil}')
        self._prepare()
        # Create direcory
        path, basename = os.path.split(fn_simu)
        if not os.path.isdir(path) and path != '':
            os.mkdir(path)

        # Prepare arguments for tms_coil
        cond = self.cond2elmdata()
        matsimnibs_list = [p.calc_matsimnibs(self.mesh) for p in self.pos]
        didt_list = [p.didt for p in self.pos]
        dirname = os.path.dirname(self.fnamecoil)
        #fname = os.path.splitext(os.path.basename(self.fnamecoil))[0]
        #fname = fname.split('.nii')[0]+'.stl'
        #fn_stl = os.path.join(dirname,fname)
        #if not os.path.isfile(fn_stl):
        #    fn_stl = None

        # Output names
        coil_name = os.path.splitext(os.path.basename(self.fnamecoil))[0]
        if coil_name.endswith('.nii'):
            coil_name = coil_name[:-4] + '_nii'
        fn_simu = ["{0}-{1:0=4d}_{2}_".format(fn_simu, i + 1, coil_name)
                   for i in range(len(self.pos))]
        output_names = [f + self.anisotropy_type + '.msh' for f in fn_simu]
        geo_names = [f + 'coil_pos.geo' for f in fn_simu]

        # call tms_coil
        fem.tms_coil(self.mesh, cond, self.cond, self.fnamecoil, self.postprocess,
                     matsimnibs_list, didt_list, output_names, geo_names,
                     solver_options=self.solver_options, n_workers=cpus)


        logger.info('Creating visualizations')
        summary = ''
        for p, n, s in zip(self.pos, output_names, fn_simu):
            p.fnamefem = n
            m = mesh_io.read_msh(n)

            if view:
                mesh_io.open_in_gmsh(n, True)

            summary += f'\n{os.path.split(s)[1][:-1]}\n'
            summary += len(os.path.split(s)[1][:-1]) * '=' + '\n'
            summary += 'Gray Matter\n\n'
            summary += m.fields_summary(roi=2)

        logger.log(25, summary)

        del cond
        gc.collect()
        return output_names, geo_names

    def run_gpc(self, fn_simu, cpus=1, tissues=[ElementTags.GM], eps=1e-2):
        from .gpc import run_tms_gpc
        gPC_regression = run_tms_gpc(
            self, fn_simu, cpus=cpus, tissues=tissues, eps=eps)
        return gPC_regression

    def add_position(self, position=None):
        """ Adds a position to the current TMSLIST

        Parameters
        -----
        position: POSITION (Optional)
            Position structure defining the coil position (Default: empty POSITION())

        Returns
        ------
        position: POSITION
            POSITION structure defining the coil position
        """
        if position is None:
            position = POSITION()

        self.pos.append(position)
        return position

    def __str__(self):
        string = "type: {0} \n" \
                 " fnamecoil: {1}, \n" \
                 " nr coil positions: {2} \n" \
                 " anisotropy: {3}" \
                 "".format(self.type,
                           self.fnamecoil,
                           len(self.pos),
                           self.anisotropy_type)
        return string

    def __eq__(self, other):
        if not isinstance(other, TMSLIST):
            return False

        if self.type != other.type:
            return False

        if self.anisotropy_type != other.anisotropy_type:
            return False

        if len(self.cond) != len(other.cond):
            return False

        for ii in range(len(self.cond)):
            if not self.cond[ii] == other.cond[ii]:
                return False

        if len(self.pos) != len(other.pos):
            return False

        for ii in range(len(self.pos)):
            if not self.pos[ii] == other.pos[ii]:
                return False
        return True


class POSITION(object):
    """ TMS coil position, orientation, name and dI/dt

    Parameters:
    -----------------------------------------
    matlab_struct(optional): output from scipy.io.loadmat
        Matlab structure defining the position

    Attributes:
    ---------------------------------------
    name: str
        name of position
    date: str
        date when stimulation was done
    matsimnibs: list with 16 numbers
        4x4 matrix defining coil center and orientation
        in simnibs coordinate system.
        HAS PREFERENCE OVER (centre, pos_y, distance)
    dIdt: Union[float, list[float]]
        Change of current in coil, in A/s
    fnamefem: str
        Name of simulation output
    ---
    The 3 variables bellow offer an alternative way to set-up a simulation.
    matsimnibs has preference over them.
    ---
    centre: list or str
        Center of the coil. Will be projected in the head surface.
        if a string, also define an eeg_cap
    pos_ydir: list of str
        Reference position for the prolongation of the coil handle.
        Will be projected in the head surface. if a string, also define an eeg_cap
    distance: float
        Distance to head
    """

    def __init__(self, matlab_struct=None):
        self.name = ''
        self.date = None
        self.matsimnibs = None
        self.didt: Union[float, list[float]] = 1e6  # corresponding to 1 A/us
        self.fnamefem = ''
        self.centre = None
        self.pos_ydir = None
        self.distance = 4.
        self.eeg_cap = None

        if matlab_struct is not None:
            self.read_mat_struct(matlab_struct)

    def _prepare(self):
        """ Prepares structure for simulation
        """
        return

    # read TMS coil position structure
    def read_mat_struct(self, p):
        """ Reads matlab structure

        Parameters
        --------------------------
        p: scipy.io.loadmat
            strucure as loaded by scipy
        """
        self.name = try_to_read_matlab_field(p, 'name', str, self.name)
        self.date = try_to_read_matlab_field(p, 'date', str, self.date)
        self.didt = try_to_read_matlab_field(p, 'didt', float, self.didt)
        self.fnamefem = try_to_read_matlab_field(
            p, 'fnamefem', str, self.fnamefem)
        try:
            self.matsimnibs = (p['matsimnibs']).tolist()
        except:
            pass
        self.centre = try_to_read_matlab_field(p, 'centre', list, self.centre)
        self.centre = try_to_read_matlab_field(p, 'center', list, self.centre)
        self.pos_ydir = try_to_read_matlab_field(
            p, 'pos_ydir', list, self.pos_ydir)
        self.distance = try_to_read_matlab_field(
            p, 'distance', float, self.distance)

        # Parse string values for centre and pos_ydir
        if self.centre and isinstance(self.centre[0], str):
            self.centre = ''.join(self.centre)

        if self.pos_ydir and isinstance(self.pos_ydir[0], str):
            self.pos_ydir = ''.join(self.pos_ydir)

    def to_dict(self) -> dict:
        """ Makes a dictionary storing all settings as key value pairs

        Returns
        --------------------
        dict
            Dictionary containing settings as key value pairs
        """
        # Generate dict from instance variables (excluding variables starting with _ or __)
        settings = {
            key:value for key, value in self.__dict__.items() 
            if not key.startswith('__')
            and not key.startswith('_')
            and not callable(value) 
            and not callable(getattr(value, "__get__", None))
            and value is not None
        }

        return settings
    
    def from_dict(self, settings: dict) -> "POSITION":
        """ Reads parameters from a dict

        Parameters
        ----------
        settings: dict
            Dictionary containing parameter as key value pairs

        Returns
        -------
        POSITION
            Self with applied settings
        """
        for key, value in self.__dict__.items():
            if key.startswith('__') or key.startswith('_') or callable(value) or callable(getattr(value, "__get__", None)):
                continue
            setattr(self, key, settings.get(key, value))

        return self

    def substitute_positions_from_cap(self, cap=None):
        if cap is None:
            cap = self.eeg_cap
        self.centre = _substitute_el(self.centre, cap)
        self.pos_ydir = _substitute_el(self.pos_ydir, cap)

    def matsimnibs_is_defined(self):
        if isinstance(self.matsimnibs, np.ndarray):
            if self.matsimnibs.ndim == 2 and \
                    self.matsimnibs.shape == (4, 4):
                return True
        elif self.matsimnibs and \
                np.array(self.matsimnibs).ndim == 2 and \
                np.array(self.matsimnibs).shape == (4, 4):
            return True
        else:
            return False

    def calc_matsimnibs(self, msh, cap=None, msh_surf=None):
        if type(msh) == str:
            msh = mesh_io.read_msh(msh)
        if cap is None:
            cap = self.eeg_cap
        if self.matsimnibs_is_defined():
            self.matsimnibs = np.array(self.matsimnibs)
        else:
            logger.info(
                'Calculating Coil position from (centre, pos_y, distance)')
            if not isinstance(self.centre, np.ndarray) and not self.centre:
                raise ValueError('Coil centre not set!')
            if not isinstance(self.pos_ydir, np.ndarray) and not self.pos_ydir:
                raise ValueError('Coil pos_ydir not set!')
            if not self.distance:
                raise ValueError('Coil distance not set!')
            self.substitute_positions_from_cap(cap=cap)
            if len(self.centre) != 3:
                raise ValueError('Coil centre must be 3-dimensional')
            if len(self.pos_ydir) != 3:
                raise ValueError('Coil pos_ydir must be 3-dimensional')
            self.matsimnibs = msh.calc_matsimnibs(
                self.centre, self.pos_ydir, self.distance, msh_surf=msh_surf)

        logger.info('matsimnibs: \n{0}'.format(self.matsimnibs))
        if ElementTags.GM_TH_SURFACE in msh.elm.tag1:
            cc_distance = np.min(np.linalg.norm(
                msh.elements_baricenters(
                )[msh.elm.tag1 == ElementTags.GM_TH_SURFACE] - self.matsimnibs[:3, 3],
                axis=1)
            )
            logger.info(f'coil-cortex distance: {cc_distance:.2f}mm')
        return self.matsimnibs

    def __eq__(self, other):
        if self.name != other.name or self.date != other.date or \
                self.matsimnibs != other.matsimnibs or self.didt != other.didt:
            return False

        else:
            return True

    def __str__(self):
        s = 'Coil Position Matrix: {0}\n'.format(self.matsimnibs)
        s += 'dIdt: {0}\n'.format(self.didt)
        if not self.matsimnibs_is_defined():
            s += 'Centre: {0}\n'.format(self.centre)
            s += 'pos_ydir: {0}\n'.format(self.pos_ydir)
            s += 'distance: {0}\n'.format(self.distance)
        return s


class TDCSLIST(SimuList):
    """Structure that defines a tDCS simulation

    Parameters
    ----------------------------------------
    matlab_struct(optinal): dict
        matlab structure as red with scipy.io, not compressed

    Attributes
    ------------------------------------------
    currents: list of floats
        current in each channel
    electrode: list of sim_struct.ELECTRODE structures
        electrodes
    cond: list
        list of COND structures with conductivity information
    anisotropy_type: property, can be 'scalar', 'vn' or 'mc'
        type of anisotropy for simulation
    postprocess: property
        fields to be calculated. valid fields are: 'v' , 'E', 'e', 'J', 'j', 'g', 's', 'D', 'q'
    """

    def __init__(self, matlab_struct=None):
        SimuList.__init__(self)
        # currents in A (not mA!; given per stimulator channel)
        self.currents = []
        self.electrode = []
        self.fnamefem = ''
        self.postprocess = 'eEjJ'

        # internal to simnibs
        self.tdcs_msh_name = None
        self.tdcs_msh_electrode_name = None

        if matlab_struct is not None:
            self.read_mat_struct(matlab_struct)

    @property
    def unique_channels(self):
        return list(set([el.channelnr for el in self.electrode]))

    def _prepare(self):
        if None in self.unique_channels:
            raise ValueError('Found a None in an Electrode Channel number',
                             'please connect all electrodes to a channel')

        if len(self.unique_channels) != len(self.currents):
            raise ValueError("Number of channels should correspond to" +
                             "the size of the currents array:\n" +
                             "unique channels:" +
                             str(self.unique_channels) + " "
                             "Currents:" + str(self.currents))

        for i in self.unique_channels:
            while len(self.cond) < (ElementTags.SALINE_START - 1) + 1 + i:
                self.cond.append(COND())
            if self.cond[(ElementTags.ELECTRODE_RUBBER_START - 1) + i].value is None:
                self.cond[(ElementTags.ELECTRODE_RUBBER_START - 1) + i].name = str(i) + ' el'
                self.cond[(ElementTags.ELECTRODE_RUBBER_START - 1) + i].value = self.cond[(ElementTags.ELECTRODE_RUBBER_START - 1)].value
            if self.cond[(ElementTags.SALINE_START - 1) + i].value is None:
                self.cond[(ElementTags.SALINE_START - 1) + i].name = str(i) + ' gel_sponge'
                self.cond[(ElementTags.SALINE_START - 1) + i].value = self.cond[(ElementTags.SALINE_START - 1)].value

        self.check_conductivities()
        if not np.isclose(np.sum(self.currents), 0):
            raise ValueError('Sum of currents should be zero!')

        if np.allclose(self.currents, 0):
            logger.warning('All current values are set to zero!')

        if len(self.unique_channels) != len(self.currents):
            raise ValueError('Number of unique channels is not equal to the number of'
                             'current values')

        for electrode in self.electrode:
            electrode.eeg_cap = self.eeg_cap

    def read_mat_struct(self, PL):
        self.read_cond_mat_struct(PL)
        self.currents = try_to_read_matlab_field(
            PL, 'currents', list, self.currents)
        self.fnamefem = try_to_read_matlab_field(
            PL, 'fnamefem', str, self.fnamefem)

        if len(PL['electrode']) > 0:
            for el in PL['electrode'][0]:
                self.electrode.append(ELECTRODE(el))

    def sim_struct2mat(self):
        mat_poslist = self.cond_mat_struct()
        mat_poslist['type'] = 'TDCSLIST'
        mat_poslist['currents'] = remove_None(self.currents)
        mat_poslist['fnamefem'] = remove_None(self.fnamefem)
        mat_poslist['electrode'] = save_electrode_mat(self.electrode)
        return mat_poslist

    def add_electrode(self, electrode=None):
        ''' Adds an electrode to the current TDCSLIST

        Parameters
        -----
        electrode: ELECTRODE (Optional)
            Electrode structure. (Default: empty ELECTRODE)

        Returns
        ------
        electrode: ELECTRODE
            electrode structure added to this TDCSLIST
        '''
        if electrode is None:
            electrode = ELECTRODE()
        self.electrode.append(electrode)
        return electrode

    def expand_to_center_surround(self, subpath, radius_surround=50, N=4,
                                  pos_dir_1stsurround=None, multichannel=False,
                                  phis_surround=None, el_surround=None):
        """
        Generate a center-surround montage (standard: 4x1) from a TDCSLIST.

        The TDCSLIST has to contain only the center electrode. Copies of this
        electrode are then placed in a circle around the centre

        Parameters
        ----------
        subpath : string
            m2m_folder of the subject
        radius_surround : float or array (N,), optional
            Distance (centre-to-centre) between the centre and surround
            electrodes. The default is 50. Either a single number (same
            radius for all electrodes) or an array with N numbers
            (N: number of surround electrodes)
        N : int, optional
            Number of surround electrodes. The default is 4.
        pos_dir_1stsurround : array (3,) or string, optional
            A position indicating the direction from center_pos to
            the position of the first surround electrode. The default is None.
        multichannel : Boolean, optional
            When set to True, a multichannel stimulator with each suround channel
            receiving 1/N-th of the of the center channel will be simulated
            (standard: False, i.e. all surround electrodes connected to the
             same return channel).
        phis_surround : array (N,), optional
            Angles in degree at which the electrodes will be placed relative to the
            direction defined by pos_dir_1stsurround. The default is None, in which
            case the electrodes will be placed at [0, 1/N*360, ..., (N-1)/N*360]
            degrees.
        el_surround : sim_struct.ELECTRODE, optional
            when supplied, this electrode definition is used for the surround
            electrodes instead of replicating the centre electrode

        """
        if len(self.electrode) != 1:
            raise ValueError(
                'The TDCSLIST has to contain exactly one ELECTRODE.')
        if not os.path.isdir(subpath):
            raise IOError('Could not find m2m-folder: {0}'.format(subpath))

        C = self.electrode[0]
        C.channelnr = 1  # Connect center to channel 1
        if not len(C.name):
            C.name = 'centre'

        # set surround channels and current strengths
        if type(self.currents) == float:
            C_current = self.currents
        else:
            C_current = self.currents[0]

        if multichannel:
            self.currents = -C_current/N*np.ones(N+1)
            self.currents[0] = C_current
            Channel_surround = np.arange(2, N+2)
        else:
            self.currents = [C_current, -C_current]
            Channel_surround = 2*np.ones(N, dtype=int)

        # get centers of surround electrodes
        ff = SubjectFiles(subpath=subpath)
        P_surround = get_surround_pos(C.centre, ff.fnamehead,
                                      radius_surround=radius_surround, N=N,
                                      pos_dir_1stsurround=pos_dir_1stsurround,
                                      phis_surround=phis_surround)

        if el_surround is not None:
            C = el_surround

        # get direction vector
        ydir = []
        if len(C.pos_ydir):
            tmp = copy.deepcopy(C)
            tmp.substitute_positions_from_cap(ff.get_eeg_cap())
            ydir = tmp.pos_ydir - tmp.centre

        # add surround electrodes to TDCSLIST
        for i in range(N):
            self.electrode.append(copy.deepcopy(C))
            El = self.electrode[-1]
            El.centre = P_surround[i]
            El.channelnr = Channel_surround[i]
            El.name = 'surround '+str(i+1)
            if len(ydir):
                El.pos_ydir = El.centre + ydir
        return

    def _place_electrodes(self):
        """ Add the defined electrodes to a mesh

        Parameters:
        ------------
        fn_out: str
            name of output file
        """

        w_elec = copy.deepcopy(self.mesh)
        w_elec.fix_tr_node_ordering()
        electrode_surfaces = [None for i in range(len(self.electrode))]
        for i, el in enumerate(self.electrode):
            logger.info('Placing Electrode:\n{0}'.format(str(el)))
            w_elec, n = el.add_electrode_to_mesh(w_elec)
            electrode_surfaces[i] = n

        w_elec.fix_th_node_ordering()
        w_elec.fix_tr_node_ordering()

        gc.collect()

        return w_elec, electrode_surfaces

    def run_simulation(self, fn_simu, cpus=1, view=True):
        """ Runs the tDCS simulation defined by this structure

        Parameters
        ----------
        fn_simu: str
            simulation path and name
        cpus: int (Optional)
            Number of parallel processes to run. Default: 1
        view: bool
            Whether to open the simulation result in Gmsh
        Returns
        ---------
        final_name: list
          List with one element: the name of the output file.
          We return a list for consistency with the TMS version
        geo_name: list
            List with the name of the geo-file that accompanies the output file.
        """
        add_scalp_to_geo = True
        fn_simu = os.path.abspath(os.path.expanduser(fn_simu))
        if not self.mesh:
            raise ValueError('The mesh for this simulation is not set')

        for c in self.cond:
            if c.distribution_type:
                logger.warning(
                    'Distribution value for conductivity found, starting gPC')
                return self.run_gpc(fn_simu, cpus=cpus), [None]

        logger.info('Began to run tDCS simulation')
        logger.info('Channels: {0}'.format(self.unique_channels))
        logger.info('Currents (A): {0}'.format(self.currents))
        self._prepare()
        path, basename = os.path.split(fn_simu)

        if not os.path.isdir(path) and path != '':
            os.mkdir(path)

        fn_no_extension, extension = os.path.splitext(fn_simu)

        mesh_elec, electrode_surfaces = self._place_electrodes()
        cond = self.cond2elmdata(mesh_elec)
        v = fem.tdcs(mesh_elec, cond, self.currents,
                     np.unique(electrode_surfaces),
                     solver_options=self.solver_options,
                     n_workers=cpus)
        m = fem.calc_fields(v, self.postprocess, cond=cond)
        final_name = fn_simu + '_' + self.anisotropy_type + '.msh'
        mesh_io.write_msh(m, final_name)
        self.fnamefem = final_name

        logger.info('Creating visualizations')
        # .geo-file
        el_geo_fn = fn_simu + '_el_currents.geo'
        self._electrode_current_geo(m, el_geo_fn)
        if add_scalp_to_geo:
            self._scalp_geo(m, el_geo_fn)
        # .opt-file
        v = m.view(
            visible_tags=_surf_preferences(m),
            visible_fields=_field_preferences(self.postprocess),
            cond_list=self.cond)
        v.add_merge(el_geo_fn)
        v.add_view(ColormapNumber=10, ColormapAlpha=.5,
                   Visible=1)  # el_currents
        if add_scalp_to_geo:
            v.add_view(ColormapNumber=8, ColormapAlpha=.3,
                       Visible=0, ShowScale=0)
        v.write_opt(final_name)

        if view:
            mesh_io.open_in_gmsh(final_name, True)

        summary = f'\n{os.path.split(fn_simu)[1]}\n'
        summary += len(os.path.split(fn_simu)[1]) * '=' + '\n'
        summary += 'Gray Matter\n\n'
        summary += m.fields_summary(roi=2)
        logger.log(25, summary)

        del m

        gc.collect()
        return [final_name], [el_geo_fn]

    def run_gpc(self, fn_simu, cpus=1, tissues=[ElementTags.GM], eps=1e-2):
        from .gpc import run_tcs_gpc
        gPC_regression = run_tcs_gpc(
            self, fn_simu, cpus=cpus, tissues=tissues, eps=eps)
        return gPC_regression

    def __str__(self):
        string = "Numer of electrodes: {0}\n".format(len(self.electrode)) + \
                 "Currents: {0}".format(self.currents)
        return string

    def __eq__(self, other):
        if self.type != other.type:
            return False

        value = self.compare_conductivities(other)

        if len(self.electrode) != len(other.electrode):
            return False

        for elec1, elec2 in zip(self.electrode, other.electrode):
            value *= elec1 == elec2

        value *= self.currents == other.currents

        return value

    def _electrode_current_geo(self, m, fn_out):
        triangles = []
        values = []
        for t, c in zip(self.unique_channels, self.currents):
            triangles.append(m.elm[m.elm.tag1 == ElementTags.SALINE_TH_SURFACE_START + t, :3])
            values.append(c * np.ones(len(triangles[-1])))

        triangles = np.concatenate(triangles, axis=0)
        values = np.concatenate(values, axis=0)
        mesh_io.write_geo_triangles(
            triangles - 1, m.nodes.node_coord,
            fn_out, values, 'electrode_currents')


class ELECTRODE(object):
    """ Class defining tDCS electrodes

    Attributes:
    ------------
    name: str
        Electrode name
    definition: 'plane' or 'conf'
        Where the electrode is defined: in a plane or in confom (suject) space
    shape: 'rect', 'ellipse' or 'custom'
        shape of electrode (plane)
    centre: list
        centre of electrode (in subject space) (plane)
    pos_ydir: list
        position in the head to define Y direction in the electrode. The Y direction is
        from the center to pos_ydir (plane)
    dimensions: list
        dimensions (x, y) of the electrode
    vertices: nx2 ndarray
        List of vertice positions in a plane, used to define custom electrodes(plane)
    thickness: int or list
        List of thickness. The number of elements specifies the type of electrode
        (conf+plane)
        1 element: simple electrode
        2 elements: electrode + gel
        3 elements: electrode in sponge
    dimensions_sponge: list
        dimensions (x, y) of the electrode

    """

    def __init__(self, matlab_struct=None):
        self.name = ''
        # how the electrode's position is evaluated: plane (2D) or conf
        # (vertices - 3D)
        self.definition = 'plane'

        # for definition = 'plane'
        self._shape = ''  # 'rect' or 'ellipse' or 'custom' for the general polygon
        self.centre = None  # centre of the 2D plane in conf coordinates
        self.dimensions = None  # for rectangle / ellipse
        self.pos_ydir = []  # second position used to define electrode orientation

        # for definition = 'plane' and 'conf'
        self.vertices = []  # array of vertices (nx2) for S.definition=='plane'
        # or (nx3) for S.definition=='conf'
        # common fields:
        # can have up to 3 arguments. 1st argument is the lower part of the
        # electrode
        self.thickness = []
        # 2nd arrgument is the rubber; 3rd is the upper
        # layer
        self.channelnr = None
        self.holes = []
        self.plug = []
        self.dimensions_sponge = None
        self.eeg_cap = None  # is overwitten by TDCSLIST._prepare()

        if matlab_struct is not None:
            self.read_mat_struct(matlab_struct)

    @property
    def type(self):
        return self.__class__.__name__

    @property
    def shape(self):
        return self._shape

    @shape.setter
    def shape(self, s):
        if s in ['rectangle', 'rect', 'ellipse', 'custom', '']:
            self._shape = s
        else:
            raise ValueError('Electrode shape must be'
                             '\'rect\', \'ellipse\', \'custom\' or \'\'')

    def _prepare(self):
        self.thickness = np.atleast_1d(self.thickness)
        if len(self.thickness) == 0:
            raise ValueError("Electrode thickness not defined!")

        if self.channelnr is None:
            logger.warning('Please connect the electrode to a channel')

    def read_mat_struct(self, el):
        self.name = try_to_read_matlab_field(el, 'name', str, self.name)
        self.definition = try_to_read_matlab_field(
            el, 'definition', str, self.definition)
        self.shape = try_to_read_matlab_field(el, 'shape', str, self.shape)
        self.dimensions = try_to_read_matlab_field(el, 'dimensions', list,
                                                   self.dimensions)
        self.centre = try_to_read_matlab_field(el, 'centre', list, self.centre)
        self.centre = try_to_read_matlab_field(el, 'center', list, self.centre)
        self.pos_ydir = try_to_read_matlab_field(
            el, 'pos_ydir', list, self.pos_ydir)
        self.thickness = try_to_read_matlab_field(
            el, 'thickness', list, self.thickness)
        self.channelnr = try_to_read_matlab_field(
            el, 'channelnr', int, self.channelnr)
        self.dimensions_sponge = try_to_read_matlab_field(el, 'dimensions_sponge', list,
                                                          self.dimensions_sponge)

        if self.dimensions is not None and len(self.dimensions) == 1:
            self.dimensions *= 2

        if self.dimensions_sponge is not None and len(self.dimensions_sponge) == 1:
            self.dimensions_sponge *= 2

        try:
            self.vertices = el['vertices'].tolist()
        except:
            # array of vertices (nx2) for S.definition=='plane'
            self.vertices = []

        try:
            el['holes']
        except ValueError:
            pass
        else:
            if len(el['holes']) > 0:
                for h in el['holes'][0]:
                    self.holes.append(ELECTRODE(h))

        try:
            el['plug']
        except ValueError:
            pass
        else:
            if len(el['plug']) > 0:
                for p in el['plug'][0]:
                    self.plug.append(ELECTRODE(p))

        if not self.definition or self.definition == '[]':  # simplify matlab syntax
            self.definition = 'plane'

        # Parse string values for centre and pos_ydir
        if self.centre and isinstance(self.centre[0], str):
            self.centre = ''.join(self.centre)

        if self.pos_ydir and isinstance(self.pos_ydir[0], str):
            self.pos_ydir = ''.join(self.pos_ydir)

    def add_electrode_to_mesh(self, mesh):
        """ Uses information in the structure in order to place an electrode

        Parameters:
        ------------
        mesh: simnibs.msh.mesh_io.Msh
            Mesh where the electrode is to be placed

        Retutns:
        -----------
        msh_elec: simnibs.msh.mesh_io.Msh
            Mesh with electrode
        tag: int
            Tag of electrode surface
        """
        self._prepare()
        m, t = electrode_placement.put_electrode_on_mesh(
            self, mesh, ElementTags.ELECTRODE_RUBBER_START + self.channelnr)
        return m, t

    def add_hole(self, hole=None):
        ''' Adds a hole to the current Electrode

        Parameters
        -----
        hole: ELECTRODE (Optional)
            Electrode structure defining the hole (Default: empty ELECTRODE)

        Returns
        ------
        hole: ELECTRODE
            electrode structure defining the hole
        '''
        if hole is None:
            hole = ELECTRODE()
        self.holes.append(hole)
        return hole

    def add_plug(self, plug=None):
        ''' Adds a plug to the current Electrode

        Parameters
        -----
        plug: ELECTRODE (Optional)
            Electrode structure defining the plug (Default: empty ELECTRODE)

        Returns
        ------
        plug: ELECTRODE
            electrode structure defining the plug
        '''
        if plug is None:
            plug = ELECTRODE()
        self.plug.append(plug)
        return plug

    def __str__(self):
        string = "definition: {0}\n".format(self.definition)
        if self.definition == 'plane':
            string += "shape: {0}\n".format(self.shape)
            string += "centre: {0}\n".format(self.centre)
            string += "pos_ydir: {0}\n".format(self.pos_ydir)
            string += "dimensions: {0}\n".format(self.dimensions)
        string += "thickness:{0}\n".format(self.thickness)
        if self.definition == 'conf' or self.shape == 'custom':
            string += "vertices: {0}\n".format(self.vertices)
        string += "channelnr: {0}\n".format(self.channelnr)
        string += "number of holes: {0}\n".format(len(self.holes))
        return string

    def substitute_positions_from_cap(self, cap=None):
        if cap is None:
            cap = self.eeg_cap
        for h in self.holes:
            h.eeg_cap = self.eeg_cap
            h.substitute_positions_from_cap(cap=cap)
        for p in self.plug:
            p.eeg_cap = self.eeg_cap
            p.substitute_positions_from_cap(cap=cap)

        self.centre = _substitute_el(self.centre, cap)
        self.pos_ydir = _substitute_el(self.pos_ydir, cap)

    def __eq__(self, other):
        try:
            return self.__dict__ == other.__dict__

        except AttributeError:
            return False


# class VOLUME:
#     ''' This is Axel's stuff, it doesn't do anything '''

#     def __init__(self, matlab_struct=None):
#         self.org = []  # used for parsing neuronavigation data; not stored permanently; optional
#         self.fname = ''  # string; points towards neuronavigation file specifying details of structural MRI; optional
#         # string; file-type of the T1 used by the sim_struct-system ('NIFTI' or 'DICOM')
#         self.ftype = ''
#         self.manufacturer = 'unknown'  # string; currently, only 'LOCALITE' is supported
#         # list of files of the T1 (one file for NIFTI; many for DICOM)
#         self.volfiles = []
#         self.img = []  # used to temporarily store the T1; can be deleted after coregistration to simnibs T1
#         self.voxsize = []  # voxel size of the T1
#         self.dim = []  # dimensions of the T1
#         self.m_qform = []  # qform of the T1
#         self.fname_conf = ''  # path and filename of the simnibs T1 of the subject
#         # 4x4 transformation matrix from sim_struct T1 to simnibs T1 (for mm-to-mm mapping of real world coordinates)
#         self.m_toconform = []

#         if matlab_struct:
#             self.read_mat_struct(matlab_struct)

#     def read_mat_struct(self, v):
#         self.org = try_to_read_matlab_field(v, 'org', list, self.org)
#         self.fname = try_to_read_matlab_field(v, 'fname', str, self.fname)
#         self.ftype = try_to_read_matlab_field(v, 'ftype', str, self.ftype)
#         self.manufacturer = try_to_read_matlab_field(v, 'manufacturer', str,
#                                                      self.manufacturer)
#         self.volfiles = try_to_read_matlab_field(
#             v, 'volfiles', list, self.volfiles)
#         self.img = try_to_read_matlab_field(v, 'img', list, self.img)
#         self.voxsize = try_to_read_matlab_field(
#             v, 'voxsize', list, self.voxsize)
#         self.dim = try_to_read_matlab_field(v, 'dim', list, self.dim)
#         self.fname_conf = try_to_read_matlab_field(
#             v, 'fname_conf', str, self.fname_conf)
#         try:
#             self.m_qform = v['m_qform'].tolist()
#         except:
#             pass
#         try:
#             self.m_toconform = v['m_toconform'].tolist()
#         except:
#             pass

#     def sim_struct2mat(self):
#         mat_vol = {}
#         mat_vol['org'] = remove_None(self.org)
#         mat_vol['fname'] = remove_None(self.fname)
#         mat_vol['ftype'] = remove_None(self.ftype)
#         mat_vol['manufacturer'] = remove_None(self.manufacturer)
#         mat_vol['volfiles'] = remove_None(self.volfiles)
#         mat_vol['img'] = remove_None(self.img)
#         mat_vol['voxsize'] = remove_None(self.voxsize)
#         mat_vol['dim'] = remove_None(self.dim)
#         mat_vol['m_qform'] = remove_None(self.m_qform)
#         mat_vol['m_toconform'] = remove_None(self.m_toconform)
#         mat_vol['fname_conf'] = remove_None(self.fname_conf)
#         return mat_vol


class LEADFIELD():
    ''' Parent class for leadfields, not meant to be called directly

    Attributes:
    ------------------------------
    fnamehead: str
        file name of mesh
    pathfem: str
        path where the leadfield should be saved
    subpath: str
        path to m2m folder. Dafault: discovered from fnamehead
    fname_tensor: str
        name of DTI tensor file. Only used in simulations with anisotropic conductivity.
        Default: discovered from subpath.
    mesh: mesh_io.Msh
        Mesh structure. Automatically loaded from fnamehead
    field: 'E' or 'J'
        Field in the leadfield (Electric field or current density). Defult: 'J'
    tissues: list
        List of tags in the mesh corresponding to the region of interest. Default: [2]
    interpolation : None | str | list | tuple
        Interpolate solution. If None, no interpolation is done. If
        'middle gm', interpolate solution to the middle of gray matter based
        on the surfaces generated by CAT/CHARM. If list or tuple, a sequence of
        filenames of surface meshes to interpolate the solution to (default =
        'middle gm').
    interpolation_subsampling : None | int
        Specify whether to interpolate to the full resolution gray matter
        surfaces (e.g., [lh/rh].central.gii) or a subsampled version of these
        (e.g., [lh/rh].central.10000.gii). The subsampled surfaces must already
        exist. Only used if 'interpolation' is 'middle gm' (default = None).
    interpolation_tissue : tuple | list
        If 'interpolation' is specified, this determines the tissue
        compartment(s) from which to interpolate the field. Interpolation
        tissues must be volumes, i.e., < 1000. (Default = [2])
    cond: list
        list of COND structures with conductivity information
    anisotropy_type: property, can be 'scalar', 'vn' or 'mc'
        type of anisotropy for simulation
    solver_options (optional): str
        Options for the FEM solver. Default: CG+AMG
    Parameters
    ------------------------
    matlab_struct: (optional) scipy.io.loadmat()
        matlab structure
    '''

    def __init__(self, matlab_struct=None):
        # : Date when the session was initiated
        self.date = time.strftime("%Y-%m-%d %H:%M:%S")
        self.time_str = time.strftime("%Y%m%d-%H%M%S")
        self.fnamehead = None
        self.subpath = None
        self.pathfem = None
        self.fname_tensor = None
        self.mesh = None
        self.field = 'E'
        self.tissues = [ElementTags.GM]
        self.interpolation = 'middle gm'
        self.interpolation_subsampling = None
        self.interpolation_tissue = [ElementTags.GM]
        self.cond = simnibs.utils.cond_utils.standard_cond()
        self.anisotropy_type = 'scalar'
        self.aniso_maxratio = 10
        self.aniso_maxcond = 2
        # The 2 variables bellow are set when the _get_vol() method is called
        # If set, they have priority over fn_tensor_nifti
        self.anisotropy_vol = None  # 4-d data with anisotropy information
        self.anisotropy_affine = None  # 4x4 affine transformation from the regular grid
        # if an anisotropic conductivity is to be used,
        self.anisotropic_tissues = [ElementTags.WM, ElementTags.GM]

        self.name = ''  # This is here only for leagacy reasons, it doesnt do anything
        self._log_handlers = []

        self.solver_options = ''
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

    def read_mat_struct(self, mat):
        """ Reads form matlab structure
        Parameters
        ------------------
        mat: scipy.io.loadmat
            Loaded matlab structure
        """
        SimuList.read_cond_mat_struct(self, mat)
        self.date = try_to_read_matlab_field(mat, 'date', str, self.date)
        self.subpath = try_to_read_matlab_field(
            mat, 'subpath', str, self.subpath)
        self.fnamehead = try_to_read_matlab_field(
            mat, 'fnamehead', str, self.fnamehead)
        self.pathfem = try_to_read_matlab_field(
            mat, 'pathfem', str, self.pathfem)
        self.field = try_to_read_matlab_field(mat, 'field', str, self.field)
        self.fname_tensor = try_to_read_matlab_field(
            mat, 'fname_tensor', str, self.fname_tensor)
        # interpolation takes different kinds of arguments so loading it is a
        # little more complex
        if len(mat['interpolation']) == 0:
            self.interpolation = None
        elif isinstance(mat['interpolation'][0], str):
            self.interpolation = mat['interpolation'][0]
        elif isinstance(mat['interpolation'][0], np.ndarray):
            self.interpolation = [s[0] for s in mat['interpolation'][0]]
        self.interpolation_subsampling = try_to_read_matlab_field(
            mat, 'interpolation_subsampling', int,
            self.interpolation_subsampling)
        self.interpolation_tissue = try_to_read_matlab_field(
            mat, 'interpolation_tissue', list, self.interpolation_tissue)
        self.tissues = try_to_read_matlab_field(
            mat, 'tissues', list, self.tissues)
        self.solver_options = try_to_read_matlab_field(mat, 'solver_options', str,
                                                       self.solver_options)

    def sim_struct2mat(self):
        mat = SimuList.cond_mat_struct(self)
        mat['type'] = self.type
        mat['date'] = remove_None(self.date)
        mat['subpath'] = remove_None(self.subpath)
        mat['fnamehead'] = remove_None(self.fnamehead)
        mat['pathfem'] = remove_None(self.pathfem)
        mat['field'] = remove_None(self.field)
        mat['fname_tensor'] = remove_None(self.fname_tensor)
        mat['interpolation'] = remove_None(self.interpolation)
        mat['tissues'] = remove_None(self.tissues)
        mat['solver_options'] = remove_None(self.solver_options)
        return mat

    def run(self, **kwargs):
        raise NotImplementedError()


class TDCSLEADFIELD(LEADFIELD):
    """Class that defines a set of simnibs simulations

    Attributes:
    ------------------------------
    fnamehead: str
        file name of mesh
    subpath: str
        path to m2m folder
    pathfem: str
        path where the leadfield should be saved
    tissues: list
        List of tags in the mesh corresponding to the region of interest in addition to
        what is defined in 'interpolation'. Default: [1006] (eye surfaces)
    interpolation : None | str | list | tuple
        Interpolate solution. If None, no interpolation is done. If
        'middle gm', interpolate solution to the middle of gray matter based
        on the surfaces generated by CAT/CHARM. If list or tuple, a sequence of
        filenames of surface meshes to interpolate the solution to (default =
        'middle gm').
    interpolation_subsampling : None | int
        Specify whether to interpolate to the full resolution gray matter
        surfaces (e.g., [lh/rh].central.gii) or a subsampled version of these
        (e.g., [lh/rh].central.10000.gii). The subsampled surfaces must already
        exist. Only used if 'interpolation' is 'middle gm' (default = None).
    interpolation_tissue : tuple | list
        If 'interpolation' is specified, this determines the tissue
        compartment(s) from which to interpolate the field. Interpolation
        tissues must be volumes, i.e., < 1000. (Default = [2])
    cond: list
        list of COND structures with conductivity information
    anisotropy_type: property, can be 'scalar', 'vn' or 'mc'
        type of anisotropy for simulation. Default: 'scalar'
    eeg_cap: str or None
        Name of eeg cap (in subject space). by default, will look for the 10-10 cap.
    electrode: ELECTRODE object or list of ELECTRODE objects
        Electrode to be used during the simulations. Default: circular electrode of 10mm
        diameter and 4mm thickness.
        If a single electrode, will be used as a template for all the electrode positions
        in the eeg cap. If a list of electrodes, each electrode should correspond to an
        entry in the cap. If is a list and eeg_cap is set to None, will place the
        electrodes based on their centre and pos_ydir properties.
    Parameters
    ------------------------
    matlab_struct: (optional) scipy.io.loadmat()
        matlab structure

    """

    def __init__(self, matlab_struct=None):
        super().__init__()
        self.eeg_cap = 'EEG10-10_UI_Jurak_2007.csv'
        self.pathfem = 'tdcs_leadfield/'
        self.tissues = [ElementTags.EYE_BALLS_TH_SURFACE]
        self.electrode = ELECTRODE()
        self.electrode.shape = 'ellipse'
        self.electrode.dimensions = [10, 10]
        self.electrode.thickness = [4]

        if matlab_struct:
            self.read_mat_struct(matlab_struct)

    @property
    def unique_channels(self):
        return list(set([el.channelnr for el in self.electrode]))

    def _prepare(self):
        """Prepares Leadfield for simulations
        relative paths are made absolute,
        empty fields are set to default values,
        check if required fields exist
        """
        LEADFIELD._prepare(self)
        sub_files = SubjectFiles(self.fnamehead, self.subpath)
        if self.eeg_cap is not None and not os.path.isfile(self.eeg_cap):
            self.eeg_cap = sub_files.get_eeg_cap(self.eeg_cap)

        logger.info('EEG Cap: {0}'.format(self.eeg_cap))

    def _add_electrodes_from_cap(self):
        ''' Reads a csv file and adds the electrodes defined to the tdcslist
        How this function behaves depends on how the "electrode" attribute is defined
        ELECTRODE: The electrode structured is copyed over. The electrode centre and
             y_dir are changed according to the information in the CSV
        list: must have the same number of electrodes defined in self.eeg_cap.
            Keeps the electrode geometry for each electrode
        '''
        if not os.path.isfile(self.eeg_cap):
            raise ValueError(
                'Could not find EEG cap file: {0}'.format(self.eeg_cap))
        type_, coordinates, extra, name, _, _ = \
            read_csv_positions(self.eeg_cap)
        count_csv = len([t for t in type_ if t in
                         ['Electrode', 'ReferenceElectrode']])

        # Create a dummy electrode object if set to None
        if self.electrode in [None, 'none', '']:
            self.electrode = ELECTRODE()
            self.electrode.shape = ''
        try:
            count_struct = len(self.electrode)
        except TypeError:
            self.electrode = \
                [copy.deepcopy(self.electrode) for i in range(count_csv)]
            count_struct = len(self.electrode)

        if count_struct != count_csv:
            raise IOError(
                'The number of electrodes in the structure is'
                ' not 0, 1 or the same number as in the CSV file'
            )
        i = 0
        ref_idx = None
        for t, c, e, n in zip(type_, coordinates, extra, name):
            if t in ['Electrode', 'ReferenceElectrode']:
                el = self.electrode[i]
                el.centre = c
                el.name = n
                if e is not None:
                    el.pos_ydir = e
                else:
                    el.pos_ydir = []
                if t == 'ReferenceElectrode':
                    ref_idx = i
                i += 1

        if ref_idx is not None:
            self.electrode[0], self.electrode[ref_idx] = \
                self.electrode[ref_idx], self.electrode[0]

        for i, el in enumerate(self.electrode):
            el.channelnr = i + 1

    def _add_el_conductivities(self):
        if None in self.unique_channels:
            raise ValueError('Found a None in an Electrode Channel number',
                             'please connect all electrodes to a channel')

        for i in self.unique_channels:
            while len(self.cond) < (ElementTags.SALINE_START - 1) + 1 + i:
                self.cond.append(COND())
            if self.cond[(ElementTags.ELECTRODE_RUBBER_START - 1) + i].value is None:
                self.cond[(ElementTags.ELECTRODE_RUBBER_START - 1) + i].name = 'el' + str(i)
                self.cond[(ElementTags.ELECTRODE_RUBBER_START - 1) + i].value = self.cond[(ElementTags.ELECTRODE_RUBBER_START - 1)].value
            if self.cond[(ElementTags.SALINE_START - 1) + i].value is None:
                self.cond[(ElementTags.SALINE_START - 1) + i].name = 'gel_sponge' + str(i + 1)
                self.cond[(ElementTags.SALINE_START - 1) + i].value = self.cond[(ElementTags.SALINE_START - 1)].value

    def _place_electrodes(self):
        """ Add the defined electrodes to a mesh """
        return TDCSLIST._place_electrodes(self)

    def _lf_name(self):
        try:
            subid = os.path.splitext(os.path.basename(self.fnamehead))[0]
            subid += '_'
        except TypeError:
            subid = ''
        try:
            eeg_cap = os.path.splitext(os.path.basename(self.eeg_cap))[0]
            eeg_cap = '_' + eeg_cap
        except TypeError:
            eeg_cap = ''
        name = '{0}leadfield{1}.hdf5'.format(subid, eeg_cap)
        return name

    def _el_name(self):
        try:
            subid = os.path.splitext(os.path.basename(self.fnamehead))[0]
            subid += '_'
        except TypeError:
            subid = ''
        try:
            eeg_cap = os.path.splitext(os.path.basename(self.eeg_cap))[0]
            eeg_cap = '_' + eeg_cap
        except TypeError:
            eeg_cap = ''
        name = '{0}electrodes{1}.msh'.format(subid, eeg_cap)
        return name

    def _mesh_roi_name(self):
        try:
            subid = os.path.splitext(os.path.basename(self.fnamehead))[0]
            subid += '_'
        except TypeError:
            subid = ''

        name = '{0}ROI.msh'.format(subid)
        return name

    def run(self, cpus=1, allow_multiple_runs=False, save_mat=True):
        ''' Runs the calculations for the leadfield

        Parameters
        -----------
        cpus: int (optional)
            Number of cpus to use. Not nescessaraly will use all cpus. Default: 1
        allow_multiple_runs: bool (optinal)
            Wether to allow multiple runs in one folder. Default: False
        save_mat: bool (optional)
            Whether to save the ".mat" file of this structure

        Returns
        ---------
        Writes the simulations

        '''
        def contains_volume(x): return any([ElementTags.TH_START <= i <= ElementTags.TH_END for i in x])
        def contains_surface(x): return any([ElementTags.TH_SURFACE_END >= i >= ElementTags.TH_SURFACE_START for i in x])

        assert isinstance(self.tissues, list)
        assert isinstance(self.interpolation, (type(None), str, list, tuple))
        assert isinstance(self.interpolation_tissue, (list, tuple))

        if self.interpolation:
            assert self.interpolation_tissue, "Interpolation tissues must be specified when interpolating solution"
            assert not contains_surface(
                self.interpolation_tissue), "Interpolation tissues must be volumes"
            if self.tissues:
                assert not contains_volume(
                    self.tissues), "Mixing volume ROIs and surface interpolation is not allowed"
        else:
            assert self.tissues, "Either 'interpolation' or 'tissues' must be specified"
            assert contains_volume(self.tissues) ^ contains_surface(
                self.tissues), "Tissue ROIs must be either surfaces OR volumes"

        if self.interpolation:
            # Leadfield is INTERPOLATED from ROI
            roi = self.interpolation_tissue + self.tissues
        else:
            # Leadfield is SAMPLED in ROI
            roi = self.tissues

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
                    'simnibs_simulation_{0}.mat'.format(self.time_str))
            )
        # For simulations without electrodes
        has_electrodes = self.electrode is not None
        # Set electrode positions
        if self.eeg_cap is not None:
            if os.path.isfile(self.eeg_cap):
                self._add_electrodes_from_cap()
            else:
                raise IOError(
                    'Could not find EEG cap file: {0}'.format(self.eeg_cap))
        try:
            len(self.electrode)
        except TypeError:
            raise ValueError(
                'Please define either an EEG cap or a list of electrodes')

        for el in self.electrode:
            if len(el.thickness) == 3:
                raise ValueError('Can not run leadfield on sponge electrodes')

        # Handle the ROI
        if np.any(np.array(self.tissues) >= ElementTags.TH_SURFACE_START) and np.any(np.array(self.tissues) <= ElementTags.TH_END):
            raise ValueError('Mixing Volumes and Surfaces in ROI!')

        if self.field != 'E' and self.field != 'J':
            raise ValueError("field parameter should be E or J. "
                             "found: {0}".format(self.field))

        self._add_el_conductivities()

        # Get names for leadfield and file of head with cap
        fn_hdf5 = os.path.join(dir_name, self._lf_name())
        fn_el = os.path.join(dir_name, self._el_name())
        if has_electrodes:
            # Place electrodes
            logger.info('Placing Electrodes')
            w_elec, electrode_surfaces = self._place_electrodes()
            mesh_io.write_msh(w_elec, fn_el)
            scalp_electrodes = w_elec.crop_mesh([ElementTags.SCALP_TH_SURFACE] + electrode_surfaces)
            scalp_electrodes.write_hdf5(fn_hdf5, 'mesh_electrodes/')
            input_type = 'tag'
            current = 1.0
            weigh_by_area = True
        else:
            # Current is injected at the nodes of the triangle on which an
            # electrode is projected weighted by the resulting barycentric
            # coordinates
            logger.info("Using point electrodes")
            w_elec = self.mesh
            skin_faces = w_elec.elm[w_elec.elm.tag1 == ElementTags.SCALP_TH_SURFACE, :3]-1
            _, v_outside, _, _ = w_elec.partition_skin_surface()
            v_outside -= 1

            pts = np.array([e.centre for e in self.electrode])
            surf = dict(points=w_elec.nodes.node_coord, tris=skin_faces)
            pttris = transformations._get_nearest_triangles_on_surface(pts, surf, 3, v_outside)
            tris, weights, projs, _ = transformations._project_points_to_surface(
                pts,
                surf,
                pttris
            )
            electrode_surfaces = surf["tris"][tris] + 1 # back to 1 indexing...
            current = weights[1:]
            input_type = "nodes"
            weigh_by_area = False

            # For debugging

            # data = np.zeros(w_elec.nodes.nr)
            # data[subset] = 1
            # w_elec.nodedata.append(mesh_io.NodeData(data, "subset"))
            # w_elec.write(self.fnamehead[:-4] + "_skin_outer_annot.msh")
            mesh_io.write_geo_spheres(pts, fn_el[:-4] + '.geo')
            mesh_io.write_geo_spheres(projs, fn_el[:-4] + '_proj.geo')

        # Write roi, scalp and electrode surfaces hdf5
        # We add WM and CSF to make the mesh convex.
        roi_fill = roi + [1, 3]

        # Write roi, scalp and electrode surfaces hdf5
        roi_msh = w_elec.crop_mesh(roi_fill)

        # 'roi':  the GM and eye surfaces ([self.tissues + [2])
        in_roi = np.in1d(roi_msh.elm.tag1, roi)

        # th_indices': the volume in the 'roi_msh'
        th_indices = roi_msh.elm.elm_number[in_roi]

        if self.interpolation is not None:
            interp_to = []
            interp_id = []
            if isinstance(self.interpolation, str):
                if self.interpolation == 'middle gm':
                    sub_files = SubjectFiles(self.fnamehead, self.subpath)
                    for hemi in sub_files.hemispheres:
                        interp_to.append(
                            mesh_io.read(
                                sub_files.get_surface(
                                    hemi,
                                    'central',
                                    self.interpolation_subsampling)
                                    )
                        )
                        interp_id.append(hemi)
                else:
                    raise ValueError(
                        'Invalid string argument to "interpolation".')
            elif isinstance(self.interpolation, (tuple, list)):
                for f in self.interpolation:
                    if os.path.isfile(f):
                        interp_to.append(mesh_io.read(f))
                        interp_id.append(f)
                    else:
                        logger.warning('{} does not exist.'.format(f))
            else:
                raise ValueError(
                    '"interpolation" must be either None, str, list, or tuple.')

            # Join meshes
            for i, m in enumerate(interp_to):
                m.elm.tag1 = (i+1)*np.ones(m.elm.nr)
                m.elm.tag2 = m.elm.tag1.copy()
            mesh_lf = interp_to[0]
            for m in interp_to[1:]:
                mesh_lf = mesh_lf.join_mesh(m)

            if len(self.tissues) > 0:
                try:
                    mesh_lf = mesh_lf.join_mesh(
                        roi_msh.crop_mesh(self.tissues))
                except ValueError:
                    logger.warning(
                        f'Could not find tissues number {self.tissues}')

            # Create interpolation matrix
            M = roi_msh.interp_matrix(
                mesh_lf.nodes.node_coord,
                out_fill='nearest',
                th_indices=th_indices,
                element_wise=True
            )

            # Define postprocessing operation

            def post(out_field, M):
                return M.dot(out_field)
            post_pro = functools.partial(post, M=M)
        else:

            # Write roi, scalp and electrode surfaces hdf5
            roi_msh = w_elec.crop_mesh(roi)

            mesh_lf = roi_msh
            post_pro = None

        fn_roi = os.path.join(dir_name, self._mesh_roi_name())
        mesh_lf.write(fn_roi)
        mesh_lf.write_hdf5(fn_hdf5, 'mesh_leadfield/')

        # Write information about number of nodes/elements in each surface
        # nix : first node of each surface
        # eix : first element of each surface
        # as found in f['mesh_leadfield']['nodes'] and
        # f['mesh_leadfield']['elm']
        if self.interpolation:
            with h5py.File(fn_hdf5, 'a') as f:
                nix = np.cumsum([0]+[m.nodes.nr for m in interp_to[:-1]])
                eix = np.cumsum([0]+[m.elm.nr for m in interp_to[:-1]])
                f['mesh_leadfield'].attrs['node_ix'] = nix
                f['mesh_leadfield'].attrs['elm_ix'] = eix
                f['mesh_leadfield'].attrs['interp_id'] = interp_id

        # Run Leadfield
        dset = 'mesh_leadfield/leadfields/tdcs_leadfield'

        logger.info('Running Leadfield')
        c = SimuList.cond2elmdata(self, w_elec)
        fem.tdcs_leadfield(
            w_elec, c, electrode_surfaces, fn_hdf5, dset,
            current=current, roi=roi,
            post_pro=post_pro, field=self.field,
            solver_options=self.solver_options,
            n_workers=cpus,
            input_type=input_type,
            weigh_by_area=weigh_by_area,
        )

        with h5py.File(fn_hdf5, 'a') as f:
            f[dset].attrs['electrode_names'] = [el.name for el in self.electrode]
            f[dset].attrs['reference_electrode'] = self.electrode[0].name
            f[dset].attrs['electrode_pos'] = [el.centre for el in self.electrode]
            f[dset].attrs['electrode_cap'] = self.eeg_cap or "none"
            f[dset].attrs['electrode_tags'] = electrode_surfaces
            f[dset].attrs['tissues'] = self.tissues
            f[dset].attrs['field'] = self.field
            f[dset].attrs['current'] = '1A'
            if self.field == 'E':
                f[dset].attrs['units'] = 'V/m'
            elif self.field == 'J':
                f[dset].attrs['units'] = 'A/m2'
            else:
                f[dset].attrs['units'] = 'Au'
            if self.interpolation is None:
                f[dset].attrs['d_type'] = 'element_data'
                f[dset].attrs['interpolation'] = 'None'
            else:
                f[dset].attrs['d_type'] = 'node_data'
                f[dset].attrs['interpolation'] = \
                    'middle gm' if self.interpolation == 'middle gm' \
                        else 'True'
            f[dset].attrs['interpolation_subsampling'] = \
                self.interpolation_subsampling if \
                    self.interpolation_subsampling is not None else 'None'

        logger.info('=====================================')
        logger.info('SimNIBS finished running leadfield')
        logger.info('Final HDF5 file:')
        logger.info(fn_hdf5)
        logger.info('=====================================')
        self._finish_logger()

    def read_mat_struct(self, mat):
        """ Reads form matlab structure
        Parameters
        ------------------
        mat: scipy.io.loadmat
            Loaded matlab structure
        """
        LEADFIELD.read_mat_struct(self, mat)
        SimuList.read_cond_mat_struct(self, mat)
        self.eeg_cap = try_to_read_matlab_field(
            mat, 'eeg_cap', str, self.eeg_cap)

        if len(mat['electrode']) > 0:
            if type(mat['electrode'][0]) is np.str_ and mat['electrode'][0] == 'none':
                self.electrode = None
            else:
                self.electrode = []
                for el in mat['electrode'][0]:
                    self.electrode.append(ELECTRODE(el))
        if self.electrode is not None and len(self.electrode) == 1:
            self.electrode = self.electrode[0]

    def sim_struct2mat(self):
        mat = LEADFIELD.sim_struct2mat(self)
        if self.electrode is None:
            mat['electrode'] = 'none'
        else:
            try:
                len(self.electrode)
                electrode = self.electrode
            except TypeError:
                electrode = [self.electrode]
            mat['electrode'] = save_electrode_mat(electrode)
        mat['eeg_cap'] = remove_None(self.eeg_cap)
        return mat


"""
    EXPORT FUNCTIONS
"""


def save_matlab_sim_struct(struct, fn):
    mat = struct.sim_struct2mat()
    scipy.io.savemat(fn, mat)


def save_electrode_mat(electrode_list):
    elec_dt = np.dtype([('type', 'O'),
                        ('name', 'O'), ('definition', 'O'),
                        ('shape', 'O'), ('centre', 'O'),
                        ('dimensions', 'O'),
                        ('pos_ydir', 'O'), ('vertices', 'O'),
                        ('thickness', 'O'),
                        ('channelnr', 'O'), ('holes', 'O'),
                        ('plug', 'O'), ('dimensions_sponge', 'O')])

    elec_mat = np.empty((0,), dtype=elec_dt)

    for elec in electrode_list:
        holes = save_electrode_mat(elec.holes)
        plug = save_electrode_mat(elec.plug)
        elec_array = np.array([('ELECTRODE',
                                remove_None(elec.name),
                                remove_None(elec.definition),
                                remove_None(elec.shape),
                                remove_None(elec.centre),
                                remove_None(elec.dimensions),
                                remove_None(elec.pos_ydir),
                                remove_None(elec.vertices),
                                remove_None(elec.thickness),
                                remove_None(elec.channelnr),
                                remove_None(holes),
                                remove_None(plug),
                                remove_None(elec.dimensions_sponge))],
                              dtype=elec_dt)
        elec_mat = np.append(elec_mat, elec_array)

    return elec_mat


def _substitute_el(pos, eeg_cap):
    if isinstance(pos, str):
        if eeg_cap:
            eeg_pos = _get_eeg_positions(eeg_cap)
            try:
                pos = eeg_pos[pos]
            except:
                raise ValueError(
                    'Could not find position {0} in cap {1}'.format(
                        pos, eeg_cap))
        else:
            raise ValueError(
                'Tried to read position: {0} but eeg_cap is not set'.format(pos))

    return pos


def _field_preferences(postprocess):
    if 'e' in postprocess:
        return ['magnE']
    elif 'j' in postprocess:
        return ['magnJ']
    elif 'E' in postprocess:
        return ['E']
    elif 'J' in postprocess:
        return ['J']
    else:
        return 'all'


def _surf_preferences(mesh):
    if ElementTags.GM_TH_SURFACE in mesh.elm.tag1:
        return [ElementTags.GM_TH_SURFACE.value]
    elif ElementTags.GM in mesh.elm.tag1:
        return [ElementTags.GM.value]
    else:
        return None


def _volume_preferences(mesh):
    if ElementTags.GM in mesh.elm.tag1:
        return [ElementTags.GM.value]
    else:
        return None


def _show_for_debugging(m, sph_centre, radius, P_centre, P_surround, surround_fit, M_sph):
    """Show some results of get_surround_pos() in gmsh for debugging."""
    import tempfile

    fn_geo = tempfile.NamedTemporaryFile(suffix='.geo').name
    mesh_io.write_geo_spheres(sph_centre.reshape((1, 3)),
                              fn_geo, name='center', mode='bw')
    mesh_io.write_geo_spheres((M_sph @ [radius, 0, 0, 1])[:3].reshape((1, 3)),
                              fn_geo, name='x', mode='ba')
    mesh_io.write_geo_spheres((M_sph @ [0, radius, 0, 1])[:3].reshape((1, 3)),
                              fn_geo, name='y', mode='ba')
    mesh_io.write_geo_spheres((M_sph @ [0, 0, radius, 1])[:3].reshape((1, 3)),
                              fn_geo, name='z', mode='ba')

    TH_DEBUG = np.arange(-1.0, 1.01, 0.1)*np.pi
    PHI_DEBUG = np.arange(0., 1.01, 0.05)*2*np.pi
    TH_DEBUG, PHI_DEBUG = np.meshgrid(TH_DEBUG, PHI_DEBUG)
    TH_DEBUG = TH_DEBUG.flatten()
    PHI_DEBUG = PHI_DEBUG.flatten()
    R_DEBUG = radius*np.ones_like(TH_DEBUG)

    pts = sph2cart(PHI_DEBUG, TH_DEBUG, R_DEBUG)
    pts = np.vstack((pts, np.ones((1, pts.shape[1]))))
    mesh_io.write_geo_spheres((M_sph @ pts)[:3, :].T, fn_geo,
                              name='sphere', mode='ba')

    mesh_io.write_geo_spheres(P_centre.reshape((1, 3)),
                              fn_geo, name='centre', mode='ba')
    for i in range(len(P_surround)):
        mesh_io.write_geo_spheres(P_surround[i].reshape((1, 3)),
                                  fn_geo, name='surr '+str(i), mode='ba')

    N_pts = 50
    for i in range(len(surround_fit)):
        tmp_centre = surround_fit[i][0]
        tmp_r = surround_fit[i][1]
        tmp_theta = surround_fit[i][2]
        tmp_theta_z0 = surround_fit[i][3]
        tmp_M = surround_fit[i][4]

        tmp_arc = np.vstack((
            tmp_r*np.sin(tmp_theta_z0 + (tmp_theta-tmp_theta_z0)
                         * np.arange(N_pts)/(N_pts-1)) + tmp_centre[0],
            np.zeros((1, N_pts)),
            tmp_r*np.cos(tmp_theta_z0 + (tmp_theta-tmp_theta_z0)
                         * np.arange(N_pts)/(N_pts-1)) + tmp_centre[1],
            np.ones((1, N_pts))
        ))
        tmp_arc = (tmp_M @ tmp_arc)[:3].T
        mesh_io.write_geo_spheres(
            tmp_arc, fn_geo, name='arc '+str(i), mode='ba')

    vis = mesh_io.gmsh_view.Visualization(m)
    vis.add_merge(fn_geo)
    vis.show()
    os.remove(fn_geo)


def sph2cart(az, el, r):  # phi, theta, radius
    """Conversion from spherical to cartesian coordinates."""
    rcos_theta = r * np.cos(el)
    pts = np.zeros(((3,) + rcos_theta.shape))
    pts[0, :] = rcos_theta * np.cos(az)
    pts[1, :] = rcos_theta * np.sin(az)
    pts[2, :] = r * np.sin(el)
    return pts


def sphereFit(pts, bounds=None):
    """
    Fit a circle or sphere to a point cloud.

    returns the radius and center points of the best fit sphere
    (adapted from https://jekel.me/2015/Least-Squares-Sphere-Fit/)

    Parameters
    ----------
    pts : array (Nx2) or (Nx3)
        point cloud

    Returns
    -------
    R : float64
        radius
    centre: ndarray (2,) or (3,)
        centre position

    """
    A = np.hstack((2*pts, np.ones((pts.shape[0], 1))))
    f = np.sum(pts**2, 1)
    C, residuals, rank, singval = np.linalg.lstsq(A, f, rcond=None)

    dim = pts.shape[1]
    R = np.sqrt(np.sum(C[0:dim]**2) + C[dim])
    return R, C[0:dim]


def get_surround_pos(center_pos, fnamehead, radius_surround=50, N=4,
                     pos_dir_1stsurround=None, phis_surround=None,
                     tissue_idx=ElementTags.SCALP_TH_SURFACE, DEBUG=False):
    """
    Determine the positions of surround electrodes.

    Parameters
    ----------
    center_pos : array (3,) or string
        Center position of the central electrode.
    fnamehead : string
        Filename of head mesh.
    radius_surround : float or array (N,), optional
        Distance (centre-to-centre) between the centre and surround
        electrodes. The default is 50. Either a single number (same
        radius for all electrodes) or an array with N numbers
        (N: number of surround electrodes)
    N : int, optional
        Number of surround electrodes. The default is 4.
    pos_dir_1stsurround : array (3,) or string, optional
        A position indicating the direction from center_pos to
        the position of the first surround electrode. The default is None.
    phis_surround : array (N,), optional
        Angles in degree at which the electrodes will be place relative to the
        direction defined by pos_dir_1stsurround. The default is None, in which
        case the electrodes will be placed at [0, 1/N*360, ..., (N-1)/N*360]
        degrees.
    tissue_idx : int, optional
        Index of the tissue on which the surround positions will be planned on.
        (standard: 1005 for skin)
    DEBUG : Boolean, optional
        When set to True, a visualization in gmsh will open for control
        (standard: False)

    Returns
    -------
    P_surround : list of arrays (3,)
        List of the centre positions of the surround electrodes.

    """

    # replace electrode name with position if needed
    # and get skin ROI around centre position
    ff = SubjectFiles(fnamehead=fnamehead)
    tmp = ELECTRODE()
    tmp.centre = center_pos
    tmp.substitute_positions_from_cap(ff.get_eeg_cap())

    m = mesh_io.read_msh(fnamehead)
    if tissue_idx < ElementTags.TH_SURFACE_START:
        tissue_idx += ElementTags.TH_SURFACE_START
    idx = (m.elm.elm_type == 2) & (
        (m.elm.tag1 == tissue_idx) | (m.elm.tag1 == tissue_idx - ElementTags.TH_SURFACE_START))
    m = m.crop_mesh(elements=m.elm.elm_number[idx])
    #P_centre = m.find_closest_element(tmp.centre)
    P_centre = project_points_on_surface(m, tmp.centre)
    idx = np.sum((m.nodes[:] - P_centre)**2,
                 1) <= (np.max(radius_surround)+10)**2
    m = m.crop_mesh(nodes=m.nodes.node_number[idx])
    idx = m.elm.connected_components()
    m = m.crop_mesh(elements=max(idx, key=np.size))

    # fit sphere to skin ROI to build local coordinate system
    #   origin: sphere center
    #   x-axis: direction of first surround
    #   z-axis: from sphere center to centre electrode
    r_sph, sph_centre = sphereFit(m.nodes[:])

    M_sph = np.eye(4)
    M_sph[:3, 3] = sph_centre
    tmp = P_centre - sph_centre
    M_sph[:3, 2] = tmp/np.linalg.norm(tmp)
    # direction of first surround
    if pos_dir_1stsurround is not None:
        # replace electrode name with position if needed
        tmp = ELECTRODE()
        tmp.centre = pos_dir_1stsurround
        tmp.substitute_positions_from_cap(ff.get_eeg_cap())
        tmp = tmp.centre - P_centre  # this is not orthogonal to Z
    else:
        # get a vector orthogonal to z-axis
        tmp = np.cross(M_sph[:3, 2], np.eye(3))
        tmp = tmp[:, np.argmax(np.linalg.norm(tmp, axis=1))]
    M_sph[:3, 1] = np.cross(M_sph[:3, 2], tmp)
    M_sph[:3, 1] /= np.linalg.norm(M_sph[:3, 1])
    M_sph[:3, 0] = np.cross(M_sph[:3, 1], M_sph[:3, 2])

    # fit arcs to the skin to get the distances accurate
    if phis_surround is not None:
        if len(phis_surround) != N:
            raise ValueError('exactly N angles are required')
        phis_surround = np.asarray(phis_surround)/180*np.pi  # convert to rad
    else:
        phis_surround = np.arange(N)/N*2*np.pi

    radius_surround = np.array(radius_surround)
    if radius_surround.size == 1:
        radius_surround = np.tile(radius_surround, N)

    N_pts = 50
    P_surround = []
    surround_fit = []
    for phi in range(N):
        theta_on_sph = radius_surround[phi]/r_sph
        arc = np.vstack((r_sph*np.sin(theta_on_sph*np.arange(N_pts)/(N_pts-1)),
                         np.zeros((1, N_pts)),
                         r_sph*np.cos(theta_on_sph *
                                      np.arange(N_pts)/(N_pts-1)),
                         np.ones((1, N_pts))))

        M_rot = np.eye(4)
        M_rot[:3, :3] = R.from_euler('z', phis_surround[phi]).as_matrix()
        M_to_world = M_sph @ M_rot
        M_from_world = np.linalg.inv(M_to_world)

        # project skin points into XZ-plane that contains the arc
        Pts = (M_to_world @ arc).T
        #Pts[:, :3] = m.find_closest_element(Pts[:, :3])
        Pts[:, :3] = project_points_on_surface(m, Pts[:, :3])
        Pts = M_from_world @ Pts.T

        # fit individual arc
        r_arc, arc_centre = sphereFit(Pts[(0, 2), :].T)

        if np.abs(arc_centre[0]) > r_arc:
            # z-axis does not intersect with circle
            # --> use initial sphere instead
            r_arc = r_sph
            arc_centre *= 0

        theta_z0_on_arc = -np.arcsin(arc_centre[0]/r_arc)
        if arc_centre[1] < np.mean(Pts[2, :]):
            theta_on_arc = radius_surround[phi]/r_arc + theta_z0_on_arc
        else:
            # best fitting arc has opposite curvature compared
            # to initial sphere
            theta_z0_on_arc = - theta_z0_on_arc + np.pi
            theta_on_arc = theta_z0_on_arc - radius_surround[phi]/r_arc

        # get centre of surround electrode
        tmp = np.array((r_arc*np.sin(theta_on_arc) + arc_centre[0],
                        0.,
                        r_arc*np.cos(theta_on_arc) + arc_centre[1],
                        1.))
        #P_surround.append(m.find_closest_element((M_to_world @ tmp).T[:3]))
        P_surround.append(project_points_on_surface(m, (M_to_world @ tmp).T[:3]))

        if DEBUG:
            surround_fit.append(
                (arc_centre, r_arc, theta_on_arc, theta_z0_on_arc, M_to_world))
    if DEBUG:
        print('achieved distances:')
        print(np.linalg.norm(np.array(P_surround)-P_centre, axis=1))
        _show_for_debugging(m, sph_centre, r_sph, P_centre,
                            P_surround, surround_fit, M_sph)

    return P_surround

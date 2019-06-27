#!/usr/bin/python2.7 -u
# -*- coding: utf-8 -*-\
'''
Global Polynomial Chaos things for SimNIBS
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.
    Copyright (C) 2017, 2018 Konstantin Weise, Guilherme B Saturnino

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
from __future__ import print_function
import os

import h5py
import numpy as np
import copy
from collections import OrderedDict


from .. import pygpc
from ..msh import mesh_io
from .sim_struct import SimuList
from . import fem
from . import coil_numpy as coil
from ..utils.simnibs_logger import logger


def write_data_hdf5(data, data_name, hdf5_fn, path='data/'):
    ''' Saves a field in an hdf5 file

    Parameters:
    ---------------
    data: np.array
        data to be saved
    data_name: str
        name of data
    path: str (optional)
        path inside hdf5 file (default: data/)
    '''
    with h5py.File(hdf5_fn, 'a') as f:
        f.create_dataset(path + data_name, data=data)


def read_data_hdf5(data_name, hdf5_fn, path='data/'):
    ''' Saves a field in an hdf5 file

    Parameters:
    ---------------
    data_name: str
        name of data
    hdf5_fn: str
        name of hdf5 file
    path: str (optional)
        path inside hdf5 file (default: data/)
    '''
    with h5py.File(hdf5_fn, 'r') as f:
        return np.squeeze(f[path + data_name])


class gPC_regression(pygpc.RegularizedRegression):
    ''' Class defining the parameters of a gPC regression

    Inherits from pygpc's "reg" class

    Attributes:
    ------------------
    pdftype: list
        List of strings with PDF types. should be 'beta' or 'normal'
        uniform distributions are beta distributed with both parameters = 1
    pdfshape: list of lists
        Parameters for pdf distributions [[p1, p2, mean1, ...], [q1, q2, var1, ...]]
    limits: list of lists
        Limits for pdfs [[min1, min2, None, ..], [max1, max2, None, ...]]
        Normal ditributions do not have limits, therefore both are set to None
    poly_idx: np.ndarray
        Indices (i, j, k ,...) for each polynomial
    order: int
        Raximum individual expansion order list [1 x DIM]
        generates individual polynomials also if maximum expansion order
        in order_max is exceeded
    order_max: int
        maximum expansion order (sum of all exponents) [1]
        the maximum expansion order considers the sum of the
        orders of combined polynomials only
    interaction_order: int
        number of random variables, which can interact with each other [1]
        all polynomials are ignored, which have an interaction order greater than the specified
    grid: pygpc.grid
        grid object gnerated in pygpc.grid, including grid.coords and grid.coords_norm
    data_file: str
        Path to file with regression data
    Parameters
    ----------------
    random_vars: list
        List of random variables. Integers is conductivities, else strings
    pdftype: list
        List of strings with PDF types. should be 'beta' or 'normal'
        uniform distributions are beta distributed with both parameters = 1
    pdfshape: list of lists
        Parameters for pdf distributions [[p1, p2, mean1, ...], [q1, q2, var1, ...]]
    limits: list of lists
        Limits for pdfs [[min1, min2, None, ..], [max1, max2, None, ...]]
        Normal ditributions do not have limits, therefore both are set to None
    poly_idx: np.ndarray
        Indices (i, j, k ,...) for each polynomial
    coords_norm: np.ndarray
        List of sampled points in normalized space
    sim_type: {'TMS', 'TCS'}
        Type of siulation 
    data_file: str, optional
        Path to file with raw data
    '''
    def __init__(self, random_vars, pdftype, pdfshape, limits, poly_idx, coords_norm, sim_type, data_file=None):
        if isinstance(pdftype, str):
            raise ValueError('pdftype must be a list of strings')

        if sim_type in ['TMS', 'TCS']:
            self._sim_type = sim_type
        else:
            raise ValueError('invalid sim_type')

        self.random_vars = random_vars
        # initialize grid class
        grid = pygpc.randomgrid(pdftype=pdftype,
                                gridshape=pdfshape,
                                limits=limits,
                                N=0)

        # overwrite grid information
        grid.coords = pygpc.denorm(coords_norm, pdftype, pdfshape, limits)
        grid.coords_norm = coords_norm

        # initialize reg class
        super().__init__(pdftype, pdfshape, limits, [0]*len(pdftype), 0, len(pdftype),
                         grid=grid)
        # enrich polynomial basis
        self.enrich_polynomial_basis(poly_idx, form_A=False)

        # set up data (HDF5) file
        self.data_file = data_file

    @property
    def sim_type(self):
        '''Simulation type, TMS or TCS'''
        return self._sim_type

    @property
    def mesh_file(self):
        ''' Name of mesh file '''
        return os.path.splitext(self.data_file)[0] + '.msh'

    def postprocessing(self, postprocessing_type, order_sobol_max=1):
        ''' Postprocessing for TMS gPC operation
        Makes the operation for the postprodesssing and saves in in the data_file

        Parameters
        -----------------------
        postprocessing_type: str
            A combination of 'v'(potential), 'E'(electric field vector),
            'e'(electric field norm), 'J'(current density vector) and
            'j'(current density norm) eg: 'eEJ'

        order_sobol_max: int (Optional)
            Maximum number of components for Sobol coefficients. Default: 1
        Returns
        -----------------------
        Writes mean, std, expansion coefficiets, sobol coefficiets, sensitivity and
        global sensitivity to file
        '''

        if self.data_file is None:
            raise ValueError('Please set a data_file before running the postprocessing')
        if not os.path.exists(self.data_file):
            raise ValueError('Could not find the data_file: {0}'.format(self.data_file))
        # make the vector into a list
        postprocessing_type = [p for p in postprocessing_type]
        for p in postprocessing_type:
            if p not in ['v', 'E', 'e', 'J', 'j']:
                raise ValueError(
                    'Unrecognized postprocessing type: {0}'.format(p))

        # load the cropped mesh
        msh = mesh_io.Msh.read_hdf5(self.data_file, 'mesh_roi/')

        lst = SimuList()
        lst._get_conductivity_from_hdf5(self.data_file)
        lst.mesh = msh

        potentials = read_data_hdf5('v_samples', self.data_file,
                                    'mesh_roi/data_matrices/')

        nr_simu = potentials.shape[0]

        if 'v' in postprocessing_type:
            self._postprocessing_core(potentials.T, 'v', False)
            postprocessing_type.remove('v')

        # This will NOT work when changing positions
        if self.sim_type == 'TMS':
            dAdt = read_data_hdf5('dAdt', self.data_file, 'mesh_roi/elmdata/')
            dAdt = mesh_io.ElementData(dAdt, mesh=msh)
        else:
            dAdt = None

        # See which fields have already been calculated
        fields_dict = dict.fromkeys(postprocessing_type)
        to_calc = []
        field_name_dict = {'E': 'E', 'e': 'normE', 'J': 'J', 'j': 'normJ'}
        for postprocess in postprocessing_type:
            try:
                fields_dict[postprocess] = read_data_hdf5(
                    field_name_dict[postprocess] +
                    '_samples', self.data_file, 'mesh_roi/data_matrices/')
            except:
                if postprocess in ['e', 'j']:
                    fields_dict[postprocess] = np.nan * np.ones((nr_simu,msh.elm.nr), dtype=float)
                elif postprocess in ['E', 'J']:
                    fields_dict[postprocess] = np.nan * np.ones((nr_simu, msh.elm.nr, 3), dtype=float)
                else:
                    raise ValueError('Unrecognized postprocessing option: ' + postprocess)
                to_calc.append(postprocess)

        # Calculate the ones which have not been calculated
        # This step is slow! I've opted for greater flexibility and reliability rather
        # than just speed. it uses calc_fields
        if len(to_calc) > 0:
            logger.info('Calculating fields: {0}'.format(''.join(to_calc)))
            for i, phi, gr in zip(range(nr_simu), potentials, self.grid.coords):
                # set the conductivities right
                for rv, g in zip(self.random_vars, gr):
                    if isinstance(rv, int):
                        lst.cond[rv - 1].value = g
                    elmdata = lst.cond2elmdata(logger_level=10)
                # sets the potential
                pot = mesh_io.NodeData(phi, mesh=msh)
                # calclate the remaining fields
                # We can use 'E' to calculate the remaining fields, if it has been
                # defined
                if 'E' in fields_dict.keys() and not np.all(np.isnan(fields_dict['E'][i])):
                    m = fem.calc_fields(pot, to_calc, cond=elmdata, dadt=dAdt,
                                        E=fields_dict['E'][i])
                #Else we also need to calculate E
                else:
                    m = fem.calc_fields(pot, to_calc, cond=elmdata, dadt=dAdt)

                for p in to_calc:
                    fields_dict[p][i] = m.field[field_name_dict[p]].value

        for p, f in fields_dict.items():
            logger.info('Expanding field: {0}'.format(field_name_dict[p]))
            if f is 'v':
                dtype = 'nodedata'
            else:
                dtype = 'elmdata'
            self._postprocessing_core(
                f, 'mesh_roi/' + dtype + '/',
                field_name_dict[p],
                order_sobol_max)


    def _postprocessing_core(self, data, path, name, order_sobol_max):
        ''' Convinience function to calculate postprocessing output '''
        data_dims = data.shape[1:]
        if data.ndim == 3:
            data = data.reshape(data.shape[0], -1)
        coeffs, error = self.expand(data)
        logger.info('Estimated error: {0:1e} for field {1}'.format(error, name))

        # Mean and std
        mean = self.mean(coeffs)
        std = self.std(coeffs)
        # Sobol
        sobol, sobol_idx = self.sobol(coeffs)
        # Filter sobol
        sobol_order = np.array([len(idx) for idx in sobol_idx])
        sobol = sobol[sobol_order <= order_sobol_max]
        sobol_idx = [s_idx for s_idx in sobol_idx if len(s_idx) <= order_sobol_max]
        # Sensitivity
        sens = self.globalsens(coeffs)

        # Put everything in Data structures
        mean = mesh_io.Data(mean.reshape(*data_dims), name=name + '_mean')
        mean.write_hdf5(self.data_file, path=path)
        std = mesh_io.Data(std.reshape(*data_dims), name=name + '_std')
        std.write_hdf5(self.data_file, path=path)
        for s, s_idx in zip(sobol, sobol_idx):
            s_name = name + '_sobol_' + \
                    '_'.join([str(self.random_vars[s_i]) for s_i in s_idx])
            sobol_dat = mesh_io.Data(s.reshape(*data_dims), name=s_name)
            sobol_dat.write_hdf5(self.data_file, path=path)

        for i, s in enumerate(sens):
            s_name = name + '_sensitivity_' + str(self.random_vars[i])
            sens_dat = mesh_io.Data(s.reshape(*data_dims), name=s_name)
            sens_dat.write_hdf5(self.data_file, path=path)

    def save_hdf5(self, data_file=None):
        """ Saves the gPC information in an hdf5 file
        pdftype, pdfshape, limits, poly_idx and coords_norm are saved
        in "gpc_object" in the data file

        Parameters:
        -----------------
        data_file (optional): str
            name of data file. default: self.data_file

        """
        if data_file is None:
            data_file = self.data_file
        if data_file is None:
            raise ValueError('Please specify a data file')
        with h5py.File(data_file, 'a') as f:
            f.attrs['type'] = self.sim_type
            if 'gpc_object' in f.keys():
                del f['gpc_object']
            f.create_dataset('gpc_object/random_vars',
                             data=np.array(self.random_vars, dtype=np.string_))
            f.create_dataset('gpc_object/pdftype',
                             data=np.array(self.pdftype, dtype='S10'))
            f.create_dataset('gpc_object/pdfshape',
                             data=np.array(self.pdfshape))
            lim = [[], []]
            lim[0] = [l if l is not None else -1e10 for l in self.limits[0]]
            lim[1] = [l if l is not None else 1e10 for l in self.limits[1]]
            f.create_dataset('gpc_object/limits',
                             data=np.array(lim))
            f.create_dataset('gpc_object/poly_idx',
                             data=np.array(self.poly_idx))
            f.create_dataset('gpc_object/grid/coords_norm',
                             data=np.array(self.grid.coords_norm))

    @classmethod
    def read_hdf5(cls, fn_hdf5):
        """ Reads gPC information from hdf5 file
        Information must have the same format as in gPC_regression.save_hdf5

        Parameteres:
        -----------------
        fn_hdf5: str
            Name of hdf5 file

        Returns:
        -----------------
        gPC_regression: gPC_regression
            regression object
        """
        with h5py.File(fn_hdf5, 'r') as f:
            sim_type = f.attrs['type']
            random_vars = f['gpc_object/random_vars'][()].tolist()
            processed = []
            for rv in random_vars:
                try:
                    processed.append(int(rv.decode()))
                except(AttributeError, ValueError):
                    processed.append(rv.decode)
            random_vars = processed
            pdftype = f['gpc_object/pdftype'][()].tolist()
            pdftype = [s.decode() for s in pdftype]
            pdfshape = f['gpc_object/pdfshape'][()].tolist()
            limits = f['gpc_object/limits'][()].tolist()
            limits[0] = [None if np.isclose(l, -1e10) else l for l in limits[0]]
            limits[1] = [None if np.isclose(l, 1e10) else l for l in limits[1]]
            poly_idx = f['gpc_object/poly_idx'][()]
            coords_norm = f['gpc_object/grid/coords_norm'][()]

        return cls(random_vars, pdftype, pdfshape, limits,
                   poly_idx, coords_norm, sim_type, data_file=fn_hdf5)

    def visualize(self):
        ''' Creates a mesh file for visualization

        Returns:
        --------
        writes a mesh file in the same folder as the hdf5 file
        msh: simnibs.msh.mesh_io
            mesh with gpc fields
        '''
        msh = mesh_io.Msh.read_hdf5(self.data_file, path='mesh_roi/')
        mesh_io.write_msh(msh, self.mesh_file)
        return msh

    def expand_quantity(self, func, field='E'):
        ''' Expand an arbitrary quantity

        Parameters
        --------......
        func: function
            Function which takes up a single argument and returns a single number or a
            vector to be expanded by gpc.
            The arguments corresponds the the electric field in the format
            [N_simulations x N_roi x 3]
        field: {v, e, E, J, j}
            field to be passed as an argument to the function

        Returns
        --------
        coeffs: ndarray
            gPC polynomial coefficients. You can then call the attributes mean, std,
            sobol, and globalsens to calculate the mean, standard deviation, sobol coefficients
            and sensitivity of you quantity of interest
        '''
        if field != 'E':
            raise NotImplementedError('For now, can only expand E')
        E = read_data_hdf5('E_samples', self.data_file, 'mesh_roi/data_matrices/')
        f = func(E)
        coeffs, cv = self.expand(f)
        logger.info('CV value: {0:1e}'.format(cv))
        return coeffs


def prep_gpc(simlist):
    cond = simlist.cond
    random_vars = []
    pdf_types = []
    parameters = []

    for i, c in enumerate(cond):
        if c.distribution_type is not None:
            random_vars.append(i + 1)
            pdf_types.append(c.distribution_type)
            parameters.append(c.distribution_parameters)

    if len(random_vars) == 0:
        raise ValueError('No random variables found for simulation') 

    # assign variables in a way gpc understands
    limits = [[], []]
    pdfshape = [[], []]
    for pdf, pars in zip(pdf_types, parameters):
        if pdf == 'uniform':
            if len(pars) != 2:
                raise ValueError('uniform random variables must have 2 parameters')
            pdfshape[0].append(1)
            pdfshape[1].append(1)
            limits[0].append(pars[0])
            limits[1].append(pars[1])
        elif pdf == 'beta':
            if len(pars) != 4:
                raise ValueError('Beta random variables must have 4 parameters')
            pdfshape[0].append(pars[0])
            pdfshape[1].append(pars[1])
            limits[0].append(pars[2])
            limits[1].append(pars[3])
        elif pdf == 'normal':
            if len(pars) != 2:
                raise ValueError('Normal random variables must have 2 parameters')
            pdfshape[0].append(pars[0])
            pdfshape[1].append(pars[1])
            limits[0].append(None)
            limits[1].append(None)
        else:
            raise ValueError('Invalid distribution_type: {0}'.format(pdf))

    # change "uniform" to "beta"
    pdf_type = ['beta' if 'uniform' == p else p for p in pdf_types]
    return random_vars, pdf_type, pdfshape, limits


def run_tms_gpc(poslist, fn_simu, cpus=1, tissues=[2], eps=1e-2,
                max_iter=1000, min_iter=2, data_poly_ratio=2):
    ''' Runs one TMS gPC for each position in the current TMSLIST

    Parameters
    ------------
    poslist: simnibs.sim_struct.TDCSLIST
        TDCSLIST structure defining the simulation
    fn_simu: str
        Output name
    cpus: int (optional)
        Number of CPUs to use
    tissues: list (Optional)
        List of tissue tags where to evaluate the electric field. Default: [2]
    eps: float (optional)
        Tolerance for gPC expansions. Default:1e-2
    max_iter: int (optinal)
        Maximum number of adaptive gPC expansion. Defaut:1000
    min_iter: int (optinal)
        Minimum number of adaptive gPC expansion interations. Defaut:2
    data_poly_ratio: int
        Ratio of number of new simulation per new polynomial. Default:2

    Returns
    --------
    fns: list
        List of mesh file names
    '''
    poslist._prepare()
    fn_simu = os.path.abspath(os.path.expanduser(fn_simu))
    if cpus > 1:
        logger.warning("Can't run GPC with multiprocessing (for now)")
        logger.warning("Setting cpus=1")
        cpus = 1

    logger.info('Running a gPC expansion with tolerance: {0:1e}'.format(eps))
    # run simulations
    random_vars, pdf_type, pdfshape, limits = prep_gpc(poslist)
    path, basename = os.path.split(fn_simu)
    if not os.path.isdir(path) and path != '':
        os.mkdir(path)

    fns = []
    for i, p in enumerate(poslist.pos):
        fn_hdf5 = fn_simu+'_{0:0=4d}_gpc.hdf5'.format(i + 1)
        if os.path.isfile(fn_hdf5):
            raise IOError('Output file ' + fn_hdf5 + ' already exists')
        matsimnibs = p.calc_matsimnibs(poslist.mesh)
        sampler = TMSgPCSampler(
            poslist.mesh, poslist, fn_hdf5,
            poslist.fnamecoil, matsimnibs, p.didt,
            roi=tissues)
        sampler.create_hdf5()
        reg, phi = pygpc.adaptive.run_reg_adaptive_grid(
            pdf_type, pdfshape, limits,
            sampler.run_simulation,
            data_poly_ratio=data_poly_ratio,
            max_iter=max_iter,
            eps=eps,
            n_cpus=cpus,
            print_function=logger.info,
            min_iter=min_iter)
        gpc_reg = gPC_regression(random_vars,
                                 pdf_type, pdfshape, limits, reg.poly_idx,
                                 reg.grid.coords_norm, 'TMS',
                                 data_file=fn_hdf5)
        gpc_reg.save_hdf5()
        gpc_reg.postprocessing(poslist.postprocess)
        gpc_reg.visualize()
        fns.append(gpc_reg.mesh_file)

    return fns


def run_tcs_gpc(poslist, fn_simu, cpus=1, tissues=[2], eps=1e-2,
                max_iter=1000, min_iter=2, data_poly_ratio=2):
    ''' Runs a tDCS gPC expansion

    Parameters
    ------------
    poslist: simnibs.sim_struct.TDCSLIST
        TDCSLIST structure defining the simulation
    fn_simu: str
        Output name
    cpus: int (optional)
        Number of CPUs to use
    tissues: list (Optional)
        List of tissue tags where to evaluate the electric field. Default: [2]
    eps: float (optional)
        Tolerance for gPC expansions. Default:1e-2
    max_iter: int (optinal)
        Maximum number of adaptive gPC expansion iterations. Defaut:1000
    min_iter: int (optinal)
        Minimum number of adaptive gPC expansion interations. Defaut:2
    data_poly_ratio(optional): int
        Ratio of number of new simulation per new polynomial. Default:2

    Returns
    --------
    fns: list
        List of mesh file names
    '''
    poslist._prepare()
    fn_simu = os.path.abspath(os.path.expanduser(fn_simu))
    if cpus > 1:
        logger.warning("Can't run GPC with multiprocessing (for now)")
        logger.warning("Setting cpus=1")
        cpus = 1
    logger.info('Running a gPC expansion with tolerance: {0:1e}'.format(eps))
    fn_hdf5 = fn_simu+'_gpc.hdf5'
    if os.path.isfile(fn_hdf5):
        raise IOError('Output file ' + fn_hdf5 + ' already exists')
    random_vars, pdf_type, pdfshape, limits = prep_gpc(poslist)

    path, basename = os.path.split(fn_simu)
    if not os.path.isdir(path) and path != '':
        os.mkdir(path)
    # place electrodes
    fn_no_extension, extension = os.path.splitext(fn_simu)
    m, electrode_surfaces = poslist._place_electrodes()
    mesh_io.write_msh(m, fn_simu + '_electrodes.msh')
    sampler = TDCSgPCSampler(
        m, poslist, fn_simu+'_gpc.hdf5',
        electrode_surfaces, poslist.currents,
        roi=tissues)
    sampler.create_hdf5()
    reg, phi = pygpc.adaptive.run_reg_adaptive_grid(
        pdf_type, pdfshape, limits,
        sampler.run_simulation,
        data_poly_ratio=data_poly_ratio,
        max_iter=max_iter,
        eps=eps,
        n_cpus=cpus,
        print_function=logger.info,
        min_iter=min_iter)
    gpc_reg = gPC_regression(random_vars,
                             pdf_type, pdfshape, limits, reg.poly_idx,
                             reg.grid.coords_norm, 'TCS',
                             data_file=fn_hdf5)
    gpc_reg.save_hdf5()
    gpc_reg.postprocessing(poslist.postprocess)
    gpc_reg.visualize()

    return [gpc_reg.mesh_file]

class gPCSampler(object):
    ''' Object used by pygpc to sample

    Attributes
    -----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh structure
    roi: list of integers
        List of tags defining the ROI
    mesh_roi: simnibs.msh.mesh_io.Msh
        Mesh of the ROI only
    fn_hdf5: string
        Name of hdf5 file with simulation output
    poslist: simnibs.simulation.sim_struct.SimuList
        Structure where the conductivity is defined
    identifiers: list
        List of random variable identifiers
    qoi_function: OrderedDict
        dictionaty with functions for each QOI.
        The first QOI will be passed to the gPC algorithm.


    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh structure
    poslist: simnibs.simulation.sim_struct.SimuList
        Structure where the conductivity is defined
    fn_hdf5: string
        Name of hdf5 file with simulation output
    roi: list of integers (Optional)
        List of tags defining the ROI. Default: [2]
    '''
    def __init__(self, m, poslist, fn_hdf5, roi=[2]):
        self.mesh = m
        self.roi = roi
        if self.roi is not None:
            self.mesh_roi = m.crop_mesh(roi)
        else:
            self.mesh_roi = m
            self.roi = np.unique(m.elm.tag1).tolist()
        self.fn_hdf5 = fn_hdf5
        self.poslist = poslist
        self._gpc_vars = prep_gpc(poslist)
        self.identifiers = self._gpc_vars[0]
        self.qoi_function = OrderedDict([('E', self._calc_E)])

    def create_hdf5(self):
        '''Creates an HDF5 file to store the data '''
        # if the hdf5 file does not exist, create it
        file_exists = os.path.exists(self.fn_hdf5)
        if file_exists:
            raise IOError('Cannot create hdf5 file: {0} '
                          'it already exists!'.format(self.fn_hdf5))

        self.mesh.write_hdf5(self.fn_hdf5, 'mesh/')
        self.mesh_roi.write_hdf5(self.fn_hdf5, 'mesh_roi/')
        self.poslist._write_conductivity_to_hdf5(self.fn_hdf5)
        with h5py.File(self.fn_hdf5, 'a') as f:
            f.create_dataset('roi', data=np.atleast_1d(np.array(self.roi, dtype=int)))

    @classmethod
    def load_hdf5(cls, fn_hdf5):
        '''Loads structure from hdf5 file '''
        mesh = mesh_io.Msh.read_hdf5(fn_hdf5, 'mesh/')
        poslist = SimuList()
        poslist._get_conductivity_from_hdf5(fn_hdf5)
        with h5py.File(fn_hdf5, 'r') as f:
            roi = f['roi'][()].tolist()
        return cls(mesh, poslist, fn_hdf5, roi)

    def record_data_matrix(self, data, name, group):
        ''' Appends or create data to the HDF5 file 

        Parameters:
        -------------
        data: np.ndarray
            Data to be appended. Will be appended along the first dimension
        name: str
            Name of data seet
        group: str
            Group where to place data set
        '''
        data = np.array(data).squeeze()
        data = np.atleast_1d(data)
        with h5py.File(self.fn_hdf5, 'a') as f:
            try:
                g = f.create_group(group)
            except:
                g = f[group]
            if name not in g.keys():
                g.create_dataset(name,
                                 shape=(0, ) + data.shape,
                                 maxshape=(None, ) + data.shape,
                                 dtype=data.dtype,
                                 chunks=(1, ) + data.shape)

            dset = g[name]
            dset.resize((dset.shape[0] + 1, ) + data.shape)
            dset[-1, ...] = data

    def run_simulation(self, random_vars):
        raise NotImplementedError('This method is to be implemented in a subclass!')

    def _calc_E(self, v, dAdt=None):
        grad = v.gradient()
        grad.assign_triangle_values()
        E = -grad.value * 1e3
        if dAdt is not None:
            E -= dAdt.value
        return E

    def _update_poslist(self, random_vars):
        poslist = copy.deepcopy(self.poslist)
        for i, iden in enumerate(self.identifiers):
            if type(iden) == int:
                poslist.cond[iden-1].value = random_vars[i]
        return poslist

    def run_N_random_simulations(self, N):
        grid = pygpc.randomgrid(pdftype=self._gpc_vars[1],
                                gridshape=self._gpc_vars[2],
                                limits=self._gpc_vars[3],
                                N=N)
        for i, x in enumerate(grid.coords):
            logger.info('Running simulation {0} out of {1}'.format(i + 1, N))
            self.run_simulation(x)


class TDCSgPCSampler(gPCSampler):
    ''' Object used by pygpc to sample a tDCS problem

    Attributes
    -----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh structure
    roi: list of integers
        List of tags defining the ROI
    mesh_roi: simnibs.msh.mesh_io.Msh
        Mesh of the ROI only
    fn_hdf5: string
        Name of hdf5 file with simulation output
    poslist: simnibs.simulation.sim_struct.SimuList
        Structure where the conductivity is defined
    gpc_vars: list
        List with variables for inputing to pygpc
    identifiers: list
        List of random variable identifiers
    qoi_function: OrderedDict
        dictionaty with functions for each QOI.
        The first QOI will be passed to the gPC algorithm, flattened.
        The QOI function should take only one argument (the potential v)
    el_tags: list of integers
        List of integers with surface tags of electrodes
    el_currents: list of floats
        List with current values for each electrode (in A)

    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh structure
    poslist: simnibs.simulation.sim_struct.SimuList
        Structure where the conductivity is defined
    fn_hdf5: string
        Name of hdf5 file with simulation output
    el_tags: list of integers
        List of integers with surface tags of electrodes
    el_currents: list of floats
        List with current values for each electrode (in A)
    roi: list of integers (Optional)
        List of tags defining the ROI. Default: [2]

    '''

    def __init__(self, mesh, poslist, fn_hdf5, el_tags, el_currents, roi=[2]):
        super(TDCSgPCSampler, self).__init__(
            mesh, poslist, fn_hdf5, roi=roi)
        self.el_tags = el_tags
        self.el_currents = el_currents

    def create_hdf5(self):
        super(TDCSgPCSampler, self).create_hdf5()
        with h5py.File(self.fn_hdf5, 'a') as f:
            f.create_dataset('el_tags', data=np.array(self.el_tags, dtype=int))
            f.create_dataset('el_currents', data=np.array(self.el_currents,
                                                          dtype=float))

    @classmethod
    def load_hdf5(cls, fn_hdf5):
        s = gPCSampler.load_hdf5(fn_hdf5)
        with h5py.File(fn_hdf5, 'r') as f:
            el_tags = f['el_tags'][()].tolist()
            el_currents = f['el_currents'][()].tolist()
        return cls(
            s.mesh, s.poslist, s.fn_hdf5, el_tags, el_currents,
            roi=s.roi)

    def run_simulation(self, random_vars):
        poslist = self._update_poslist(random_vars)
        cond = poslist.cond2elmdata(self.mesh)
        v = fem.tdcs(
            self.mesh, cond, self.el_currents,
            self.el_tags, units='mm')


        self.mesh.nodedata = [v]
        cropped = self.mesh.crop_mesh(self.roi)
        v_c = cropped.nodedata[0]
        self.mesh.nodedata = []

        qois = []
        for qoi_name, qoi_f in self.qoi_function.items():
            qois.append(qoi_f(v_c))

        self.record_data_matrix(random_vars, 'random_var_samples', '/')
        self.record_data_matrix(v.value, 'v_samples', 'mesh/data_matrices')
        self.record_data_matrix(v_c.value, 'v_samples', 'mesh_roi/data_matrices')
        for qoi_name, qoi_f in self.qoi_function.items():
            self.record_data_matrix(
                qois[-1], qoi_name + '_samples', 'mesh_roi/data_matrices')

        del cropped
        del cond
        del v
        del v_c

        return np.atleast_1d(qois[0]).reshape(-1)


class TMSgPCSampler(gPCSampler):
    ''' Object used by pygpc to sample a TMS problem

    Attributes
    -----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh structure
    roi: list of integers
        List of tags defining the ROI
    mesh_roi: simnibs.msh.mesh_io.Msh
        Mesh of the ROI only
    fn_hdf5: string
        Name of hdf5 file with simulation output
    poslist: simnibs.simulation.sim_struct.SimuList
        Structure where the conductivity is defined
    gpc_vars: list
        List with variables for inputing to pygpc
    identifiers: list
        List of random variable identifiers
    qoi_function: OrderedDict
        dictionaty with functions for each QOI.
        The first QOI will be passed to the gPC algorithm.
        The QOI function should take only 2 arguments: (the potential v and dAdt)
    matsimnibs: np.ndarray
        Matrix defining coil position
    dIdt: float
        Current intensity in coil

    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh structure
    poslist: simnibs.simulation.sim_struct.SimuList
        Structure where the conductivity is defined
    fn_hdf5: string
        Name of hdf5 file with simulation output
    pos: simnibs.simulation.sim_struct.POS
       Coil position definition
    roi: list of integers (Optional)
        List of tags defining the ROI. Default: [2]

    '''

    def __init__(self, mesh, poslist, fn_hdf5, fnamecoil, matsimnibs, didt, roi=[2]):
        super(TMSgPCSampler, self).__init__(
            mesh, poslist, fn_hdf5, roi=roi)
        self.matsimnibs = np.ascontiguousarray(matsimnibs)
        self.didt = didt
        self.fnamecoil = fnamecoil
        self.constant_dAdt = True

    def create_hdf5(self):
        super(TMSgPCSampler, self).create_hdf5()
        with h5py.File(self.fn_hdf5, 'a') as f:
            f.create_dataset('matsimnibs', data=self.matsimnibs)
            f.create_dataset('didt', data=np.array(self.didt, dtype=float))
            f.create_dataset('fnamecoil', data=np.array(self.fnamecoil, dtype=np.string_))

    @classmethod
    def load_hdf5(cls, fn_hdf5):
        s = gPCSampler.load_hdf5(fn_hdf5)
        with h5py.File(fn_hdf5, 'r') as f:
            matsimnibs = f['matsimnibs'][()]
            fnamecoil = f['fnamecoil'][()].decode()
            didt = f['didt'][()]

        return cls(
            s.mesh, s.poslist, s.fn_hdf5, fnamecoil,
            matsimnibs, didt, roi=s.roi)

    def run_simulation(self, random_vars):
        poslist = self._update_poslist(random_vars)
        cond = poslist.cond2elmdata(self.mesh)
        if self.constant_dAdt:
            try:
                dAdt = self.dAdt
                dAdt_roi = self.dAdt_roi
            except AttributeError:
                dAdt = coil.set_up_tms_dAdt(
                    self.mesh,
                    self.fnamecoil,
                    self.matsimnibs,
                    didt=self.didt,
                    fn_geo=self.fn_hdf5[:-5]+'_coil.geo')
                if isinstance(dAdt, mesh_io.NodeData):
                    dAdt = dAdt.node_data2elm_data()
                dAdt.field_name = 'dAdt'
                dAdt.write_hdf5(self.fn_hdf5, 'mesh/elmdata/')
                self.dAdt = dAdt
                self.mesh.elmdata = [dAdt]
                cropped = self.mesh.crop_mesh(self.roi)
                dAdt_roi = cropped.elmdata[0]
                dAdt_roi.write_hdf5(self.fn_hdf5, 'mesh_roi/elmdata/')
                self.mesh.elmdata = []
                self.dAdt_roi = dAdt_roi
        else:
            raise NotImplementedError

        v = fem.tms_dadt(self.mesh, cond, dAdt)
        self.mesh.nodedata = [v]
        cropped = self.mesh.crop_mesh(self.roi)
        v_c = cropped.nodedata[0]
        self.mesh.nodedata = []

        qois = []
        for qoi_name, qoi_f in self.qoi_function.items():
            qois.append(qoi_f(v_c, dAdt_roi))

        self.record_data_matrix(random_vars, 'random_var_samples', '/')
        self.record_data_matrix(v.value, 'v_samples', 'mesh/data_matrices')
        self.record_data_matrix(v_c.value, 'v_samples',
                                'mesh_roi/data_matrices')
        for qoi_name, qoi_f in self.qoi_function.items():
            self.record_data_matrix(
                qois[-1], qoi_name + '_samples', 'mesh_roi/data_matrices')

        del cropped
        del cond
        del v

        return np.atleast_1d(qois[0]).reshape(-1)

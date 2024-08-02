# -*- coding: utf-8 -*-\
'''
    Functions for assembling and solving FEM systems
'''

import gc
import multiprocessing
import time
import copy
import warnings
import atexit
import h5py
import logging
import numpy as np
import scipy.sparse as sparse
from simnibs.mesh_tools import gmsh_view
from simnibs.simulation.tms_coil.tms_coil import TmsCoil

from simnibs.utils.mesh_element_properties import ElementTags

from ..mesh_tools import mesh_io
from ..utils import cond_utils as cond_lib
import mumps
from . import pardiso
from . import petsc_solver
from ..utils.simnibs_logger import logger

DEFAULT_SOLVER_OPTIONS = \
        '-ksp_type cg ' \
        '-ksp_rtol 1e-10 -pc_type hypre ' \
        '-pc_hypre_type boomeramg ' \
        '-pc_hypre_boomeramg_coarsen_type HMIS'

'''
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2018-2019  Guilherme B Saturnino


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

_PETSC_IS_INITILIZED = False
def _initialize_petsc():
    global _PETSC_IS_INITILIZED
    if not _PETSC_IS_INITILIZED:
        # Initialize PETSc
        petsc_solver.petsc_initialize()
        # Register the finalizer for PETSc
        atexit.register(petsc_solver.petsc_finalize)
        _PETSC_IS_INITILIZED = True

class MUMPS_Solver:
    def __init__(self, A=None, isSymmetric=True, log_level=20):
        self.log_level = log_level
        start = time.time()
        self.ctx = mumps.Context()
        self.ctx.set_matrix(A, symmetric=isSymmetric)
        logger.log(self.log_level, f'{time.time()-start:.2f} seconds to init solver')
        start = time.time()
        self.ctx.analyze()
        logger.log(self.log_level, f'{time.time()-start:.2f} seconds to analyze matrix')
        start = time.time()
        self.ctx.factor()
        logger.log(self.log_level, f'{time.time()-start:.2f} seconds to factorize matrix')
        
    def solve(self, b):
        start = time.time()
        x = self.ctx._solve_dense(b)
        logger.log(self.log_level, f'{time.time()-start:.2f} seconds to solve system')
        return x

        
def calc_fields(potentials, fields, cond=None, dadt=None, units='mm', E=None):
    ''' Given a mesh and the electric potentials at the nodes,
    calculates the fields

    Parameters
    ----------
    potentials: simnibs.msh.mesh_io.NodeData
        NodeData field with potentials.
        Attention: the mesh property should be set
    fields: Any combination of 'vEeJjsDg'
        Fields to output
        v: electric potential at the nodes
        E: Electric field at the elements
        e: Electric field magnitude at the elements
        J: Current density at the elements
        j: Current density magnitude at the elements
        s: Conductivity at the elements
        D: dA/dt at the nodes
        g: gradient of the potential at the elements
    cond: simnibs.mesh.mesh_io.ElementData (optional)
        Conductivity at the elements, used to calculate J, j and s.
        Might be a scalar or a tensor.
    dadt: simnibs.msh.mesh_io.NodeData (optional)
        dA/dt at the nodes for TMS simulations
    units: {'mm' or 'm'} (optional)
        Mesh units, either milimiters (mm) or meters (m). Default: mm
    E: np.ndarray or simnibs.msh.mesh_io.ElementData
        Electric field, if it has been already calculated

    Returns
    -------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh object with the calculated fields
    '''
    if units == 'mm':
        scaling_factor = 1e3
    elif units == 'm':
        scaling_factor = 1
    else:
        raise ValueError('Invalid units: {0}'.format(units))

    mesh = potentials.mesh
    if mesh is None:
        raise ValueError('potential does not have the mesh property set')

    assert mesh.nodes.nr == potentials.nr, \
        ('The number of nodes in the mesh and of data points in the'
         ' potential does not match')

    if cond is not None:
        assert mesh.elm.nr == cond.nr, \
            ('The number of elements in the mesh and of data points in the'
             ' conductivity field does not match')

    out_mesh = copy.deepcopy(mesh)
    out_mesh.elmdata = []
    out_mesh.nodedata = []
    if 'v' in fields:
        out_mesh.nodedata.append(
            mesh_io.NodeData(
                potentials.value,
                name='v',
                mesh=out_mesh))
    if 'D' in fields:
        if dadt is None:
            warnings.warn('Cannot calculate D field: needs dadt input')
        elif isinstance(dadt, mesh_io.NodeData):
            out_mesh.nodedata.append(
                mesh_io.NodeData(
                    dadt.value,
                    name='D',
                    mesh=out_mesh))
        else:
            out_mesh.elmdata.append(
                mesh_io.ElementData(
                    dadt.value,
                    name='D',
                    mesh=out_mesh))

    if any(f in ['E', 'e', 'J', 'j', 'g', 's'] for f in fields):
        if 'g' in fields or E is None:
            grad = potentials.gradient() * scaling_factor
            grad.assign_triangle_values()
            grad.field_name = 'g'
            grad.mesh = out_mesh

        if 'g' in fields:
            out_mesh.elmdata.append(grad)

        if E is None:
            if dadt is not None:
                if isinstance(dadt, mesh_io.NodeData):
                    dadt_elmdata = dadt.node_data2elm_data()
                else:
                    dadt_elmdata = dadt
                dadt_elmdata.assign_triangle_values()
                E = mesh_io.ElementData(
                    -grad.value - dadt_elmdata.value,
                    name='E',
                    mesh=out_mesh)
            else:
                E = mesh_io.ElementData(
                    -grad.value,
                    name='E',
                    mesh=out_mesh)
        else:
            if not isinstance(E, mesh_io.ElementData):
                E = mesh_io.ElementData(E, name='E', mesh=out_mesh)
            if E.nr != out_mesh.elm.nr:
                raise ValueError(
                    'Provided E does not have the same number of'
                    ' samples as the mesh!')
            if E.nr_comp != 3:
                raise ValueError('Provided E does not have 3 components!')

        if 'E' in fields:
            out_mesh.elmdata.append(E)
        if 'e' in fields:
            e = np.linalg.norm(E.value, axis=1)
            out_mesh.elmdata.append(
                mesh_io.ElementData(
                    e, name='magnE', mesh=out_mesh))

        if any(f in ['J', 'j', 's'] for f in fields):
            if cond is None:
                raise ValueError(
                    'Cannot calculate J, j os s field: No conductivity input')
            cond.assign_triangle_values()
            if 's' in fields:
                cond.field_name = 'conductivity'
                cond.mesh = out_mesh
                if cond.nr_comp == 9:
                    out_mesh.elmdata += cond_lib.visualize_tensor(cond, out_mesh)
                else:
                    out_mesh.elmdata.append(cond)

            J = mesh_io.ElementData(calc_J(E, cond),
                                 name='J', mesh=out_mesh)

            if 'J' in fields:
                out_mesh.elmdata.append(J)
            if 'j' in fields:
                j = np.linalg.norm(J.value, axis=1)
                out_mesh.elmdata.append(
                    mesh_io.ElementData(
                        j, name='magnJ', mesh=mesh))

    return out_mesh


def calc_J(E, cond):
    '''Calculates J

    Parameters
    ----------
    E: ndarray of mesh_io.Data
        Electric field. Nx3 vector

    cond: ndarray or mesh_io.Data
        Conductivity. A Nx1 (scalar) or Nx9 (tensor) vector.

    Returns
    -------
    J: ndarray
        Current density
    '''
    if isinstance(E, mesh_io.Data):
        E = E.value
    if isinstance(cond, mesh_io.Data):
        cond = cond.value
    if cond.ndim == 1 or cond.shape[1] == 1:
        J = E*cond[:, None]
    elif cond.shape[1] == 9:
        J = np.einsum('ikj, ik -> ij', cond.reshape(-1, 3, 3), E)
    else:
        raise ValueError('Conductivity should be a Nx1 or an Nx9 vector')
    return J


class dofMap(object):
    '''Dictionary mapping degrees of freedom to 1-N indices '''
    def __init__(self, inverse):
        self.from_inverse_map(inverse)

    def from_inverse_map(self, inverse):
        self.inverse = np.array(inverse, dtype=int)
        self._map = -99999 * np.ones(np.max(inverse) + 1, dtype=int)
        self._map[inverse] = np.arange(len(inverse), dtype=int)
        self.nr = len(inverse)

    def order_like(self, other_dof, array=None):
        sort = self.inverse.argsort()
        pos = np.searchsorted(self.inverse[sort], other_dof.inverse)
        indices = sort[pos]
        if array is not None:
            return dofMap(self.inverse[indices]), array[indices]
        else:
            return dofMap(self.inverse[indices])

    def __getitem__(self, index):
        ret = self._map.__getitem__(index)
        if np.any(ret == -99999):
            raise IndexError(
                'Index out of range')
        return ret

    def __eq__(self, other):
        return np.all(self._map == other._map) and np.all(self.inverse == other.inverse)


class DirichletBC(object):
    '''Class Defining  dirichlet boundary conditions

    Attributes
    ----------
    nodes: list
        List of nodes there the BC should be applied

    values: list
        Value at each node

    Parameters
    ----------
    nodes: list
        List of nodes there the BC should be applied

    values: list
        Value at each node
    '''
    def __init__(self, nodes, values):
        assert len(nodes) == len(values), 'There should be one value for each node'
        self.nodes = nodes
        self.values = values

    def apply(self, A, b, dof_map):
        ''' Applies the dirichlet BC to the system

        Parameters
        ----------
        A: scipy.sparse.csr
            Sparse matrix
        b: numpy array or None
            Righ-hand side. if None, it will return None
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b

        Returns
        -------
        A: scipy.sparse.csr
            Sparse matrix, modified
        b: numpy array
            Righ-hand side, modified
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b, modified
        '''
        if np.any(~np.in1d(self.nodes, dof_map.inverse)):
            raise ValueError('BC node indices not found in dof_map')
        stay = np.ones(A.shape[0], dtype=bool)
        stay[dof_map[self.nodes]] = False
        if b is not None:
            b, _ = self.apply_to_rhs(A, b, dof_map)
        A, dof_map = self.apply_to_matrix(A, dof_map)
        return A, b, dof_map

    def apply_to_rhs(self, A, b, dof_map):
        ''' Applies the dirichlet BC to the system

        Parameters
        ----------
        A: scipy.sparse.csr
            Sparse matrix
        b: numpy array or None
            Righ-hand side. if None, it will return None
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b

        Returns
        -------
        b: numpy array
            Righ-hand side, modified
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b, modified
        '''
        if np.any(~np.in1d(self.nodes, dof_map.inverse)):
            raise ValueError('BC node indices not found in dof_map')
        stay = np.ones(b.shape[0], dtype=bool)
        stay[dof_map[self.nodes]] = False
        b = np.atleast_2d(b)
        if b.shape[0] < b.shape[1]:
            b = b.T
        A = A.tocsc()
        s = A[:, dof_map[self.nodes]].dot(self.values)
        s = np.atleast_2d(s)
        if s.shape[0] < s.shape[1]:
            s = s.T
        b = b - s
        b = b[stay]
        dof_map = dofMap(dof_map.inverse[stay])
        return b, dof_map

    def apply_to_matrix(self, A, dof_map):
        ''' Applies the dirichlet BC to the system

        Parameters
        ----------
        A: scipy.sparse.csr
            Sparse matrix
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b

        Returns
        -------
        A: numpy array
            System matrix, modified
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b, modified
        '''
        if np.any(~np.in1d(self.nodes, dof_map.inverse)):
            raise ValueError('BC node indices not found in dof_map')
        stay = np.ones(A.shape[0], dtype=bool)
        stay[dof_map[self.nodes]] = False
        A = A.tocsr()
        # Remove rows
        for n in dof_map[self.nodes]:
            A.data[A.indptr[n]:A.indptr[n+1]] = 0
        A = A[stay, :]
        A.eliminate_zeros()
        A = A.tocsc()
        # Remove columns
        for n in dof_map[self.nodes]:
            A.data[A.indptr[n]:A.indptr[n+1]] = 0
        A = A[:, stay]
        A.eliminate_zeros()
        dof_map = dofMap(dof_map.inverse[stay])
        return A, dof_map


    def apply_to_vector(self, v, dof_map):
        ''' Apply to an lhs vector by just removing the entries
        '''
        stay = np.ones(dof_map.nr, dtype=bool)
        stay[dof_map[self.nodes]] = False
        return v[stay], dofMap(dof_map.inverse[stay])

    def apply_to_solution(self, x, dof_map):
        ''' Applies the dirichlet BC to a solution

        Parameters
        ----------
        x: numpy array
            Righ-hand side
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b

        Returns
        -------
        x: numpy array
            Righ-hand side, modified
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b, modified
        '''
        if np.any(np.in1d(self.nodes, dof_map.inverse)):
            raise ValueError('Found DOFs already defined')
        dof_inverse = np.hstack((dof_map.inverse, self.nodes))
        x = np.atleast_2d(x)
        if x.shape[0] < x.shape[1]:
            x = x.T
        x = np.vstack((x, np.tile(self.values, (x.shape[1], 1)).T))
        return x, dofMap(dof_inverse)

    @classmethod
    def join(cls, list_of_bcs):
        ''' Join many BCs into one

        Parameters
        ----------
        list_of_bcs: list
            List of DirichletBC objects

        Returns
        -------
        bc: DirichletBC
            DirichletBC corresponding to a union of the list
        '''
        nodes = np.hstack([bc.nodes for bc in list_of_bcs])
        values = np.hstack([bc.values for bc in list_of_bcs])
        return cls(nodes, values)


def set_ground_at_nodes(mesh, nodes=None):
    """Set a Dirichlet boundary condition with a value of 0 at the specified
    nodes. If `nodes` is None, the single lowest node in the mesh (smallest z
    coordinate) is used. In the case of TMS and direct approach this is used to
    ensure uniqueness of the solution.

    Add a dirichlet BC to the single lowest node to make the problem SPD.

    Parameters
    ----------
    nodes : int | list | ndarray
        The nodes to use as ground.

    Returns
    -------
    bc : DirichletBC
        The Dirichlet boundary condition object.
    """
    if nodes is None:
        nodes = mesh.nodes.node_number[mesh.nodes.node_coord[:, 2].argmin()]
    nodes = np.atleast_1d(nodes)
    return DirichletBC(nodes, np.zeros_like(nodes, dtype=float))


class FEMSystem(object):
    ''' Class defining the equation system to be used in FEM calculations

    Parameters
    ----------
    mesh: mesh_io.Msh
       Mesh object
    cond: mesh_io.ElementData
        conductivity information
    dirichlet: list of DirichletBC (optional)
        Dirichlet boundary conditions
    units: {'mm' or 'm'} (optional)
        Units of the mesh nodes. Default: mm
    store_G: bool (optional)
        Wether to store the gradient matrix. Default: False
    solver_options: str
        Options to be used by the solver. Default: DEFAULT_SOLVER_OPTIONS
    solver_loglevel: int (optional)
        sets the log level of the standard logging messages of the solvers.
        Default: logging.INFO 

    Attributes
    ----------
    mesh: mesh_io.Msh
       Mesh object
    cond: mesh_io.ElementData
        conductivity information
    dirichlet: DirichletBC (optional)
        Dirichlet boundary condition
    units: {'mm' or 'm'}
        Units of the mesh nodes
    A: scipy.sparse.csr_matrix
        Sparse system matric
    dof_map: dofMap
        Mapping between rows/columns of A and DOFs

    Notes
    -----
    Once created, do NOT change the attributes of this class.

    '''
    def __init__(self, mesh, cond, dirichlet=None, units='mm', store_G=False,
                 solver_options=None, solver_loglevel=logging.INFO):
        if units in ['mm', 'm']:
            self.units = units
        else:
            raise ValueError('Invalid unit: {0}'.format(units))
        self.solver_loglevel = solver_loglevel
        self._mesh = mesh
        if isinstance(cond, mesh_io.ElementData):
            cond = cond.value.squeeze()
            if cond.ndim == 2:
                cond = cond.reshape(-1, 3, 3)
        if self.mesh.elm.nr != len(cond):
            raise ValueError('Please define one conductivity for each element')
        self._cond = cond
        self._dirichlet = dirichlet
        self._dof_map = dofMap(mesh.nodes.node_number)
        self._A = None
        self._solver = None
        self._G = None # Gradient operator
        self._D = None # Gradient matrix
        if solver_options in [None, '']:
            self._solver_options = DEFAULT_SOLVER_OPTIONS
        else:
            self._solver_options = solver_options
        self.assemble_fem_matrix(store_G=store_G)

    @property
    def mesh(self):
        return self._mesh

    @property
    def cond(self):
        return self._cond

    @property
    def dirichlet(self):
        return self._dirichlet

    @property
    def A(self):
        return self._A

    @property
    def dof_map(self):
        return self._dof_map

    def assemble_fem_matrix(self, store_G=False):
        ''' Assembly of the l.h.s matrix A. !Only works with symmetric matrices!
        Based in the OptVS algorithm in Cuvelier et. al. 2016 '''
        logger.info('Assembling FEM Matrix')
        start = time.time()
        msh = self.mesh
        cond = self.cond[msh.elm.elm_type == 4]
        th_nodes = msh.elm.node_number_list[msh.elm.elm_type == 4]
        G = _gradient_operator(msh)
        if store_G:
            self._G = G  # stores the operator in case we need it later (TMS)
        vols = _vol(msh)
        dof_map = self.dof_map
        self._A = _assemble_matrix(vols, G, th_nodes, cond, dof_map,
                                   units=self.units)
        if np.any(np.diff(self.A.indptr) == 0):
            raise ValueError('Found a column of zeros in the stiffness matrix'
                             ' disconected nodes?')

        time_assemble = time.time() - start
        logger.info(
            '{0:.2f}s to assemble FEM matrix'.format(time_assemble))

    def prepare_solver(self):
        '''Prepares the object to solve FEM systems

        Notes
        -----
        After running this method, do NOT change any attributes of the class!
        '''
        logger.info(f'Using solver options: {self._solver_options}')
        A = sparse.csc_matrix(self.A, copy=True)
        A.sort_indices()
        dof_map = copy.deepcopy(self.dof_map)
        if self.dirichlet is not None:
            A, dof_map = self.dirichlet.apply_to_matrix(A, dof_map)

        if self._solver_options == 'pardiso':
            self._solver = pardiso.Solver(A, log_level=self.solver_loglevel)
        elif self._solver_options == 'mumps':
            self._solver = MUMPS_Solver(A, isSymmetric=True, log_level=self.solver_loglevel)
        else:
            _initialize_petsc()
            self._A_reduced = A  # We need to save this as PETSc does not copy the vectors
            self._solver = petsc_solver.Solver(self._solver_options, A, log_level=self.solver_loglevel)

    def solve(self, b=None):
        ''' Solves the FEM system

        Parameters
        ----------
        b: np.ndarray (Optional):
            Right-hand side. If not set, will assume zeros

        Returns
        -------
        x: ndarray
            array with solution

        Notes
        -----
        After running this method, do NOT change any attributes of the class!
        '''
        logger.debug('Solving FEM System')
        if b is None:
            b = np.zeros(self.dof_map.nr, dtype=float)
        else:
            b = np.copy(b)

        if self._solver is None:
            self.prepare_solver()

        # We also need the A matrix here because the DOFs change
        dof_map = copy.deepcopy(self.dof_map)
        if self.dirichlet is not None:
            b, dof_map = self.dirichlet.apply_to_rhs(self.A, b, dof_map)

        x = self._solver.solve(b)

        if self.dirichlet is not None:
            x, dof_map = self.dirichlet.apply_to_solution(x, dof_map)
        dof_map, x = dof_map.order_like(self.dof_map, array=x)

        return np.squeeze(x)


    def calc_gradient(self, v):
        ''' Calculates gradients

        Parameters
        ----------
        v: np.ndarray
            Array with fields at the nodes. Can be 1d or 2d (n_nodes x n)

        Returns
        -------
        grad: np.ndarray
            Array with gradients at the tetrahedra. Can be 2d if v in 1d or 3d (n_th x 3
            x n), if v is 2d.
        '''
        if self._G is None:
            G = _gradient_operator(self.mesh)
        else:
            G = self._G
        if self._D is None:
            self._D = grad_matrix(self.mesh, G)
        grad = self._D.dot(v)
        if v.ndim == 1:
            return grad.reshape(-1, 3)
        elif v.ndim == 2:
            return grad.reshape(-1, 3, v.shape[1])


# Classes for specific types of FEM systems

class TMSFEM(FEMSystem):
    def __init__(
            self,
            mesh,
            cond,
            solver_options=None,
            units='mm',
            store_G=True,
            solver_loglevel=logging.INFO
        ):
        '''Set up a TMS problem.

        Parameters
        ----------
        mesh: simnibs.mesh_io.msh.Msh
            Mesh structure.
        cond: ndarray or simnibs.mesh_io.msh.ElementData
            Conductivity of each element.
        solver_options: str
            Options to be used by the solver. Default: DEFAULT_SOLVER_OPTIONS
        '''
        dirichlet_bc = set_ground_at_nodes(mesh)
        super().__init__(mesh, cond, dirichlet_bc, units, store_G, solver_options, solver_loglevel)

    def assemble_rhs(self, dadt):
        '''Assemble the right-hand side for a TMS simulation.

        Parameters
        ----------
        dadt: NodeData or ElementData
            dA/dt field at each node or element

        Returns
        -------
        b: np.array
            Right-hand side

        References
        ----------
        Gomez, Luis J., Moritz Dannhauer, and Angel V. Peterchev. "Fast
            computational optimization of TMS coil placement for individualized
            electric field targeting." bioRxiv (2020).
        '''
        msh = self.mesh
        cond = self.cond[msh.elm.elm_type == 4]
        if isinstance(dadt, mesh_io.NodeData):
            dadt = dadt.node_data2elm_data()
        dadt = dadt.value
        dadt = dadt[msh.elm.elm_type == 4]
        G = _gradient_operator(msh) if self._G is None else self._G
        vols = _vol(msh)
        th_nodes = msh.elm.node_number_list[msh.elm.elm_type == 4]
        # integrate in each node of each element, the value for repeated nodes will be summed
        # together later
        elm_node_integral = np.zeros((len(th_nodes), 4), dtype=np.float64)
        if cond.ndim == 1:
            sigma_dadt = cond[:, None] * dadt
        elif cond.ndim == 3:
            sigma_dadt = np.einsum('aij, aj -> ai', cond, dadt)
        else:
            raise ValueError('Invalid cond array')

        for i in range(4):
            elm_node_integral[:, i] = \
                        -vols * (sigma_dadt * G[:, i, :]).sum(axis=1)

        if self.units == 'mm':
            elm_node_integral *= 1e-6

        b = np.bincount(self.dof_map[th_nodes.reshape(-1)],
                        elm_node_integral.reshape(-1))
        #self.b = np.atleast_2d(self.b).T
        return b


class TDCSFEMDirichlet(FEMSystem):
    def __init__(
            self,
            mesh,
            cond,
            electrodes,
            potentials,
            # input_type="tag", # currently only supports "tag"!
            solver_options=None,
            units='mm',
            store_G=False,
            solver_loglevel=logging.INFO
        ):
        '''Set up a TDCS problem using Dirichlet boundary conditions in all
        electrodes.

        Parameters
        ----------
        mesh: simnibs.mesh_io.msh.Msh
            Mesh structure.
        cond: ndarray or simnibs.mesh_io.msh.ElementData
            Conductivity of each element.
        electrode_tags: list
            list of the surfaces where the dirichlet BC is to be applied.
        potentials: list
            list of the potentials each surface is to be set.
        solver_options: str
            Options to be used by the solver. Default: DEFAULT_SOLVER_OPTIONS
        '''
        self.electrodes = electrodes
        self.potentials = potentials
        # self.input_type = input_type

        dirichlet_bc = self._init_dirichlet_bcs(mesh)
        super().__init__(mesh, cond, dirichlet_bc, units, store_G, solver_options, solver_loglevel)

    def _init_dirichlet_bcs(self, mesh):
        """Set Dirichlet boundary conditions on all electrodes."""

        assert len(self.electrodes) == len(self.potentials)
        bcs = []
        for t, p in zip(self.electrodes, self.potentials):
            elements_in_surface = (mesh.elm.tag1 == t) * (mesh.elm.elm_type == 2)
            if np.sum(elements_in_surface) == 0:
                raise ValueError('Did not find any surface with tag: {0}'.format(t))
            n = np.unique(mesh.elm.node_number_list[elements_in_surface, :3])
            bcs.append(DirichletBC(n, p * np.ones_like(n, dtype=float)))
        return DirichletBC.join(bcs)            # =====


class TDCSFEMNeumann(FEMSystem):
    def __init__(
            self,
            mesh,
            cond,
            ground_electrode,
            input_type="tag",
            weigh_by_area=True,
            solver_options=None,
            units='mm',
            store_G=False,
            solver_loglevel=logging.INFO
        ):
        '''Set up a TDCS problem using Dirichlet boundary conditions in the
        ground electrode and Neumann boundary conditions in the other
        electrodes.

        Parameters
        ----------
        mesh: simnibs.mesh_io.msh.Msh
            Mesh structure
        cond: ndarray or simnibs.mesh_io.msh.ElementData
            Conductivity of each element
        ground_electrode: int
            Tag of the ground electrode surface
        solver_options: str (optional)
            Options to be used by the solver. Default: DEFAULT_SOLVER_OPTIONS
        input_type: 'tag' or "nodes" (optional)
            Input can be either the tag of the electrode surface (default) or a
            list of nodes
        '''
        assert input_type in {"tag", "nodes"}

        self.ground_electrode = ground_electrode
        self.input_type = input_type
        if input_type == "tag" and not weigh_by_area:
            warnings.warn("Parameters `weigh_by_area=False` and `input_type='tag' are incompatible. Forcing weigh_by_area=True!")
        self.weigh_by_area = weigh_by_area | (input_type == "tag")
        self.areas = mesh.nodes_areas() if self.weigh_by_area else None

        dirichlet_bc = self._init_dirichlet_bc(mesh)
        super().__init__(mesh, cond, dirichlet_bc, units, store_G, solver_options, solver_loglevel)

    def _init_dirichlet_bc(self, mesh):
        """Set Dirichlet boundary condition on the ground electrode only."""

        # The first surface is set to a DirichletBC
        if self.input_type == 'tag':
            # Find the nodes in the tag
            elements_in_surface = (
                (mesh.elm.tag1 == self.ground_electrode) *
                (mesh.elm.elm_type == 2)
            )
            if np.sum(elements_in_surface) == 0:
                raise ValueError(
                    'Did not find any surface with tag: {0}'.format(self.ground_electrode))
            n = np.unique(mesh.elm.node_number_list[elements_in_surface, :3])
        elif self.input_type == "nodes":
            n = np.atleast_1d(self.ground_electrode)
        else:
            raise ValueError('Invalid value for `input_type`')
        return set_ground_at_nodes(mesh, n)

    def assemble_rhs(self, electrodes, currents):
        '''Assemble the right-hand side for a TDCS simulation with Neumann
        boundary conditions.

        Parameters
        ----------
        electrodes: list
            list of the surface tags or nodes where the currents to be applied.
            WARNING: should NOT include the ground electrode
        currents: list
            list of the currents in each surface
            WARNING: should NOT include the ground electrode

        Returns
        -------
        b: np.ndarray
            Right-hand-side of FEM system
        '''
        # if self.input_type == "nodes":
        #     electrodes = np.atleast_2d(electrodes)
        assert len(electrodes) == len(currents)
        b = np.zeros(self.dof_map.nr, dtype=np.float64)
        for e, c in zip(electrodes, currents):
            if self.input_type == 'tag':
                # Find the nodes in the tag
                nodes = np.unique(
                    self.mesh.elm[
                        (self.mesh.elm.tag1 == e) *
                        (self.mesh.elm.elm_type == 2),
                        :3
                    ])
            elif self.input_type == "nodes":
                nodes = e
            else:
                raise ValueError('Invalid value for `input_type`')
            b += self._rhs_node(nodes, c)
        return b

    def _rhs_node(self, nodes, current):
        '''Assemble the Neumann RHS on a set of nodes.'''
        b = np.zeros(self.dof_map.nr, dtype=np.float64)
        ix = self.dof_map[nodes]
        b[ix] = current
        if self.weigh_by_area:  # self.areas is not None
            b[ix] *= self.areas[nodes] / self.areas[nodes].sum()
        return b


class DipoleFEM(FEMSystem):
    def __init__(
            self,
            mesh,
            cond,
            solver_options=None,
            units='mm',
            store_G=True,
            solver_loglevel=logging.INFO
        ):
        '''Set up an electric dipole simulation using the selected source
        model (i.e., the "direct" approach).

        Parameters
        ----------
        mesh: simnibs.mesh_io.msh.Msh
            Mesh structure
        cond: ndarray or simnibs.mesh_io.msh.ElementData
            Conductivity of each element
        solver_options: str (optional)
            Options to be used by the solver. Default: DEFAULT_SOLVER_OPTIONS
        '''
        dirichlet_bc = set_ground_at_nodes(mesh)
        super().__init__(mesh, cond, dirichlet_bc, units, store_G, solver_options, solver_loglevel)

    # def assemble_rhs(self, primary_j, source_model):
    def assemble_rhs(self, dip_pos, dip_mom, source_model):
        '''Assemble the right-hand side of the system equation using the
        specified source model. Supported source models

        * Partial Integration
            Sources are placed on nodes of the tetrahedron in which a dipole
            resides.
            Loads are calculated by computing the similarity of the dipole
            moment with the gradient in a particular element and weighing each
            node accordingly (*integrating* the contributions from each dot
            product). (Node positions are related to the basis vectors of the
            element, e.g., e1=v1-v0, which are related to the gradients.)

            This is insensitive to the position of a dipole within a given
            volume element!

        * St. Venant
            Sources are placed on the node closest to the dipole as well as
            first degree neighbors. A linear system of equations is solved
            such that the sum of the monopolar moments is 0, the dipole moment
            is equal to the specified moment, and the squared dipole moment is
            zero as well. Thus, there are 7 equations with n unknowns where n
            is the number of nodes on which loads are placed.

        Parameters
        ----------
        dip_pos: ndarray (n, 3)
            Dipole positions.
        dip_mom: ndarray (n, 3)
            Dipole moments in ampere-meter (Am). In case one wants to supply a
            dipole moment matching a particular (primary) current density, J,
            instead, use the following relationship between the dipole moment,
            p, and J, to calculate the proper value of p.
                p = int J dV
                p = J*V_i
                J = p/V_i
            Here V_i is the volume of the ith element.
        source_model: str
            Select `partial integration` or `st. venant`.

        Returns
        -------
        b: ndarray
            Right-hand side.

        References
        ----------
        Weinstein, David, Leonid Zhukov, and Chris Johnson. "Lead-field bases
            for electroencephalography source imaging." Annals of biomedical
            engineering 28.9 (2000): 1059-1065.
        '''
        dip_pos = np.atleast_2d(dip_pos).astype(float)
        dip_mom = np.atleast_2d(dip_mom).astype(float)
        assert dip_pos.shape == dip_mom.shape, "`dip_pos` and `dip_mom` must have the same dimensions"
        n_dip = dip_pos.shape[0]

        if self.units == "mm":
            dip_mom *= 1e3 # Am to Amm

        if source_model == "partial integration":
            # This is very similar to the TMS case. Guilherme's original
            # comment was
            #   The RHS is basically the same as the TMS one, but with the
            #   conductiviy is already incorporated in dipole_vectors

            # find element containing each source
            tetra_idx = self.mesh.find_tetrahedron_with_points(dip_pos, False)

            # convert 1 to 0 indexing!
            tetra_idx -= 1
            tetra_idx = np.atleast_1d(tetra_idx)
            nodes_idx = self.mesh.elm.node_number_list[tetra_idx]-1

            n_nodes_idx = list(map(len, nodes_idx))
            rows = np.concatenate(nodes_idx)
            cols = np.concatenate([i*np.ones(n_nodes_idx[i], dtype=int) for i in range(n_dip)], dtype=int)

            grad_op = _gradient_operator(self.mesh) if self._G is None else self._G
            # grad_op is only computed for tetrahedra so reindex `tetra_idx`
            reindexer = np.zeros(self.mesh.elm.nr, int)
            reindexer[self.mesh.elm.tetrahedra-1] = np.arange(len(self.mesh.elm.tetrahedra))
            tetra_idx = reindexer[tetra_idx]

            # compute the integrals
            data = np.ravel(grad_op[tetra_idx] @ dip_mom[..., None])
            b = sparse.csc_matrix((data, (rows, cols)), shape=(self.mesh.nodes.nr, n_dip))
        elif source_model == "st. venant":
            raise NotImplementedError("St. Venant implementation has not been validated yet...")
            # find closest mesh node to each source position
            _, src_idx = self.mesh.nodes.find_closest_node(dip_pos, return_index=True)
            # find all elements of closest nodes and keep only tetrahedra and
            # reindex to 0
            element_idx = [
                np.intersect1d(
                    self.mesh.elm.find_all_elements_with_node(i), self.mesh.elm.tetrahedra,
                    assume_unique=True,
                )-1 for i in src_idx
            ]
            # find all unique nodes of those elements
            node_indices = [np.unique(self.mesh.elm.node_number_list[e]-1) for e in element_idx]

            b = compute_st_venant_loads(dip_pos, dip_mom, node_indices, self.mesh.nodes.node_coord)
        else:
            raise ValueError

        return np.array(b.todense()).squeeze()


def compute_st_venant_loads(dip_pos, dip_mom, dip_node_indices, vertices):
    """Compute the loads for the St. Venant source model.

    Parameters
    ----------
    dip_pos : ndarray (n, 3)
        Dipole positions.
    dip_mom : ndarray (n, 3)
        Dipole moments.
    dip_node_indices : list (n, ) of ndarray (m, 3)
        Indices of the nodes surrounding each dipole (i.e., all nodes
        connected to the closest node of each dipole).

    Returns
    -------
    b : scipy.sparse.csc_matrix
        Sparse matrix of shape (len(vertices), len(dip_pos)) containing the RHS
        corresponding to each dipole in columns.
    """
    # constants
    aref = 20
    lambda_ = 1e-5
    r = 1

    ntotal = sum(map(len, dip_node_indices))
    rows = np.concatenate(dip_node_indices)
    cols = np.concatenate([i*np.ones(len(dip_node_indices[i]), dtype=int) for i in range(len(dip_node_indices))], dtype=int)
    data = np.zeros(ntotal)
    for i, (vert, d, ni) in enumerate(zip(dip_pos, dip_mom, dip_node_indices)):
        nn = len(ni)

        # vector from vertex to neighboring vertices
        v = vertices[ni] - np.atleast_2d(vert)
        v /= aref

        # system matrix

        # in fieldtrip the order is 1, x, x**2, 1, y, y**2, 1, z, z**2 like
        #
        #   X = np.ones((9, nn))
        #   X[1::3] = v.T
        #   X[2::3] = v.T**2
        #
        # not sure why they use 9 rows as we only have 7 eqs.

        # in simnibs the order is 1 x y z x**2 y**2 z**2
        X = np.ones((7, nn))
        X[1:4] = v.T
        X[4:] = v.T**2

        # moments (RHS)

        # in fieldtrip
        #   t = np.zeros(9)
        #   t[1::3] = d / aref

        t = np.zeros(7)
        t[1:4] = d/aref

        # 3*nn x nn
        W = np.zeros((3 * nn, nn))
        W[:nn] = np.diag(v[:,0])**r
        W[nn:2*nn] = np.diag(v[:,1])**r
        W[2*nn:] = np.diag(v[:,2])**r

        # Calculate the loads, q
        q = np.linalg.solve(X.T @ X + lambda_ * W.T @ W, X.T @ t)
        data[cols==i] = q
    return sparse.csc_matrix((data, (rows, cols)), shape=(vertices.shape[0], dip_pos.shape[0]))


def assemble_diagonal_mass_matrix(msh, units='mm'):
    ''' Assemble a Mass matrix by doing a first-order integration at the nodes
    Results in a diagonal matrix

    Parameters
    ----------
    msh: simnibs.msh.mesh_io.Msh
        Mesh structure
    units: {'m' or 'mm'}
        Units where the mesh is defined. The matrix will be scaled accordingly

    Returns
    -------
    M: scipy.sparse.csc_matrix:
        Diagonal matrix
    '''
    th_nodes = msh.elm.node_number_list[msh.elm.elm_type == 4]
    vols = _vol(msh)

    # I'm using csc for consistency
    dof_map = dofMap(msh.nodes.node_number)
    M = sparse.csc_matrix((dof_map.nr, dof_map.nr))
    for i in range(4):
        M += sparse.csc_matrix(
            (.25 * vols,
             (dof_map[th_nodes[:, i]],
              dof_map[th_nodes[:, i]])),
            shape=(dof_map.nr, dof_map.nr))

    M.sum_duplicates()
    if units == 'mm':
        M *= 1e-9

    return M


def _gradient_operator(msh, volume_tag=None):
    ''' G calculates the gradient of a function in each tetrahedra
    The way it works: The operator has 2 parts
    G = T^{-1}A
    A is a projection matrix
    A = [-1, 1, 0, 0]
        [-1, 0, 1, 0]
        [-1, 0, 0, 1]
    And T is the transfomation to baricentric coordinates
    '''
    if volume_tag is None:
        th = msh.nodes[msh.elm.node_number_list[msh.elm.elm_type == 4]]
    else:
        th = msh.nodes[msh.elm.node_number_list[(msh.elm.elm_type == 4) *
                                                (msh.elm.tag1 == volume_tag)]]
    A = np.hstack([-np.ones((3, 1)), np.eye(3)])
    G = np.linalg.solve(th[:, 1:4] - th[:, 0, None], A[None, :, :])
    G = np.transpose(G, (0, 2, 1))
    return G

def _assemble_matrix(vols, G, th_nodes, cond, dof_map, units='mm'):
    '''Based in the OptVS algorithm in Cuvelier et. al. 2016 '''
    A = sparse.csc_matrix((dof_map.nr, dof_map.nr), dtype=np.float64)
    if cond.ndim == 1:
        vGc = vols[:, None, None]*G*cond[:, None, None]
    elif cond.ndim == 3:
        vGc = vols[:, None, None]*np.einsum('aij, ajk -> aik', G, cond)
    else:
        raise ValueError('Invalid cond array')
    ''' Off-diagonal '''
    for i in range(4):
        for j in range(i+1, 4):
            Kg = (vGc[:, i, :]*G[:, j, :]).sum(axis=1)
            A += sparse.csc_matrix(
                (Kg, (dof_map[th_nodes[:, i]],
                      dof_map[th_nodes[:, j]])),
                shape=(dof_map.nr, dof_map.nr),
                dtype=np.float64)

    A += A.T
    ''' Diagonal'''
    for i in range(4):
        Kg = (vGc[:, i, :]*G[:, i, :]).sum(axis=1)
        A += sparse.csc_matrix(
            (Kg, (dof_map[th_nodes[:, i]],
                  dof_map[th_nodes[:, i]])),
            shape=(dof_map.nr, dof_map.nr),
            dtype=np.float64)

    if units == 'mm':
        A *= 1e-3  # * 1e6 from the gradiend operator, 1e-9 from the volume

    A.eliminate_zeros()
    return A


def grad_matrix(msh, G=None, split=False):
    ''' Matrix that calculates the gradients at the elements

    Parameters
    ----------
    msh: simnibs.msh.mesh_io
        Mesh structure
    G: sparse matrix (optional)
        G matrix to avoid re-calculations. If not set, it will be re-calculated
    split: bool (optional)
        If true, will return a list of sparse matrices, one for each component.
        Default: False

    Returns
    -------
    (if split=False, default):
    D: sparse matrix
        Matrix such that D.dot(x).reshape(-1, 3) is grad(x)
        The triangle values are also assigned

    (if split=True):
    D: list of sparse matrices
        list of sparse matrices such that D[i].dot(x)
        is the i-th component of grad(x)
        The triangle values are also assigned

    '''
    if G is None:
        G = _gradient_operator(msh)
    th = msh.elm.elm_number[msh.elm.elm_type == 4] - 1
    tr = msh.elm.elm_number[msh.elm.elm_type == 2] - 1
    cp = msh.find_corresponding_tetrahedra() - 1
    G_expanded = np.empty((msh.elm.nr, 4, 3), dtype=float)
    G_expanded[th] = G
    G_expanded[tr] = G_expanded[cp]
    G = G_expanded
    th_nodes = np.zeros((msh.elm.nr, 4), dtype=int)
    th_nodes[th] = msh.elm.node_number_list[msh.elm.elm_type == 4]
    th_nodes[tr] = th_nodes[cp]
    if not split:
        D = sparse.csc_matrix((3 * msh.elm.nr, msh.nodes.nr))
        for j in range(3):
            for i in range(4):
                D += sparse.csc_matrix(
                    (G[:, i, j], (j + 3 * np.arange(msh.elm.nr),
                     th_nodes[:, i] - 1)),
                    shape=D.shape)
    if split:
        D = []
        for j in range(3):
            D.append(sparse.csc_matrix((msh.elm.nr, msh.nodes.nr)))
            for i in range(4):
                D[-1] += sparse.csc_matrix(
                    (G[:, i, j], (np.arange(msh.elm.nr),
                     th_nodes[:, i] - 1)),
                    shape=D[-1].shape)

    return D


def _vol(msh, volume_tag=None):
    '''Volume of the tetrahedra '''
    if volume_tag is None:
        th = msh.nodes[msh.elm.node_number_list[msh.elm.elm_type == 4]]
    else:
        th = msh.nodes[msh.elm.node_number_list[(msh.elm.elm_type == 4) *
                                                (msh.elm.tag1 == volume_tag)]]
    return np.abs(np.linalg.det(th[:, 1:] - th[:, 0, None])) / 6.


def tdcs(mesh, cond, currents, electrode_surface_tags, n_workers=1, units='mm',
         solver_options=None):
    ''' Simulates a tDCS electric potential.

    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh file with geometry information.
    cond: simnibs.msh.mesh_io.ElementData
        An ElementData field with conductivity information.
    currents: list or ndarray
        A list of currents going though each electrode.
    electrode_surface_tags: list
        A list of the indices of the surfaces where the dirichlet BC is to be
        applied.

    Returns
    -------
    potential: simnibs.msh.mesh_io.NodeData
        Total electric potential
    '''
    assert len(currents) == len(electrode_surface_tags),\
        'there should be one channel for each current'

    surf_tags = np.unique(mesh.elm.tag1[mesh.elm.elm_type == 2])
    assert np.all(np.in1d(electrode_surface_tags, surf_tags)),\
        'Could not find all the electrode surface tags in the mesh'

    assert np.isclose(np.sum(currents), 0),\
        'Currents should sum to 0'

    ref_electrode = electrode_surface_tags[0]
    total_p = np.zeros(mesh.nodes.nr, dtype=float)

    n_workers = min(len(currents) - 1, n_workers)
    if n_workers == 1:
        for el_surf, el_c in zip(electrode_surface_tags[1:], currents[1:]):
            total_p += _sim_tdcs_pair(
                mesh, cond, ref_electrode, el_surf, el_c, units, solver_options)
    else:
        with multiprocessing.Pool(processes=n_workers) as pool:
            sims = []
            for el_surf, el_c in zip(electrode_surface_tags[1:], currents[1:]):
                sims.append(
                    pool.apply_async(
                        _sim_tdcs_pair,
                        (mesh, cond, ref_electrode, el_surf, el_c, units, solver_options)))
            for s in sims:
                total_p += s.get()
            pool.close()
            pool.join()

    return mesh_io.NodeData(total_p, 'v', mesh=mesh)


def _sim_tdcs_pair(mesh, cond, ref_electrode, el_surf, el_c, units, solver_options):
    logger.info('Simulating electrode pair {0} - {1}'.format(
        ref_electrode, el_surf))

    s = TDCSFEMDirichlet(mesh, cond,  [ref_electrode, el_surf], [0., 1.], solver_options)
    v = s.solve()

    v = mesh_io.NodeData(v, name='v', mesh=mesh)
    flux = np.array([
        _calc_flux_electrodes(v, cond,
                              [el_surf - 1000, el_surf - 600,
                               el_surf - 2000, el_surf - 1600],
                              units=units),
        _calc_flux_electrodes(v, cond,
                              [ref_electrode - 1000, ref_electrode - 600,
                               ref_electrode - 2000, ref_electrode - 1600],
                              units=units)])
    current = np.average(np.abs(flux))
    error = np.abs(np.abs(flux[0]) - np.abs(flux[1])) / current
    if error <= 0.1:
        logger.info('Estimated current calibration error: {0:.1%}'.format(error))
    else:
        logger.warning(f'The current calibration error exceeded 10%! Estimated error value: {error*100:.2f}%')
    del s
    gc.collect()
    return el_c / current * v.value


def _calc_flux_electrodes(v, cond, el_volume, scalp_tag=[ElementTags.SCALP, ElementTags.SCALP_TH_SURFACE], units='mm'):
    # Set-up a mesh with a mesh
    m = copy.deepcopy(v.mesh)
    m.nodedata = [v]
    m.elmdata = [cond]
    # Select mesh nodes wich are is in one electrode as well as the scalp
    # Triangles in scalp
    tr_scalp = np.in1d(m.elm.tag1, scalp_tag) * (m.elm.elm_type == 2)
    if not np.any(tr_scalp):
        raise ValueError('Could not find skin surface')
    tr_scalp_nodes = m.elm.node_number_list[tr_scalp, :3]
    tr_index = m.elm.elm_number[tr_scalp]

    # Tetrahehedra in electrode
    th_el = np.in1d(m.elm.tag1, el_volume) * (m.elm.elm_type == 4)
    if not np.any(th_el):
        raise ValueError('Could not find electrode volume')
    th_el_nodes = m.elm.node_number_list[th_el]
    nodes_el = np.unique(th_el_nodes)
    th_index = m.elm.elm_number[th_el]

    # Triangles in interface
    tr_interface = tr_index[
        np.all(np.isin(tr_scalp_nodes, nodes_el), axis=1)]
    if len(tr_interface) == 0:
        raise ValueError('Could not find skin-electrode interface')
    keep = np.hstack((th_index, tr_interface))

    # Now make a mesh with only tetrahedra and triangles in interface
    crop = m.crop_mesh(elements=keep)
    crop.elm.tag1 = np.ones_like(crop.elm.tag1)
    crop.elm.tag2 = np.ones_like(crop.elm.tag2)

    # Calculate J in the interface
    crop = calc_fields(crop.nodedata[0], 'J', crop.elmdata[0], units=units)

    # Calculate flux
    flux = crop.elmdata[0].calc_flux()
    if units == 'mm':
        flux *= 1e-6

    del m
    del crop
    return flux


def tms_dadt(mesh, cond, dAdt, solver_options=None):
    ''' Simulates a TMS electric potential from a dA/dt field.

    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh file with geometry information
    cond: simnibs.msh.mesh_io.ElementData
        An ElementData field with conductivity information
    dAdt: simnibs.msh.mesh_io.NodeData or simnibs.msh.mesh_io.ElementData
        dAdt information

    Returns
    -------
    v:  simnibs.msh.mesh_io.NodeData
        NodeData instance with potential at the nodes
    '''
    s = TMSFEM(mesh, cond, solver_options)
    b = s.assemble_rhs(dAdt)
    v = s.solve(b)
    
    del s, b
    gc.collect()
    return mesh_io.NodeData(v, name='v', mesh=mesh)


def tms_coil(mesh, cond, cond_list, fn_coil, fields, matsimnibs_list, didt_list,
             output_names, geo_names=None, solver_options=None, n_workers=1):
    '''Simulates TMS fields using a coil + matsimnibs + dIdt definition.

    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh structure
    cond: simnibs.msh.mesh_io.ElementData
        Conductivity field
    fields: str or list of str
        Fields to be calculated for each position
    fn_coil: string
        Relative path of coil file
    matsimnibs_list: list
        List of "matsimnibs" matrices, one per position
    didt_list: list
        List of dIdt values, one per position or a list of dIdt values, one per stimulator
    output_names: list of str
        List of output mesh file names, one per position
    geo_names: list of str
        List of output mesh file names, one per position
    solver_options: str
        Options for the solver
    n_workers: int
        Number of workers to use
    fn_stl: string
        Name of stl-file for coil visualization

    Returns
    -------
    Writes output meshes to the files specified in output_names
    '''
    assert len(matsimnibs_list) == len(didt_list)
    assert len(output_names) == len(didt_list)
    n_sims = len(matsimnibs_list)
    n_workers = min(n_sims, n_workers)

    if geo_names is None:
        geo_names = [None for i in range(n_sims)]

    S = TMSFEM(mesh, cond, solver_options)
    if n_workers == 1:
        _set_up_global_solver(S)
        for matsimnibs, didt, fn_out, fn_geo in zip(
                matsimnibs_list, didt_list, output_names, geo_names):
            _run_tms(
                mesh, cond, cond_list, fn_coil, fields,
                matsimnibs, didt, fn_out, fn_geo)
        _finalize_global_solver()
    else:
        with multiprocessing.Pool(processes=n_workers,
                                  initializer=_set_up_global_solver,
                                  initargs=(S,)) as pool:
            sims = []
            for matsimnibs, didt, fn_out, fn_geo in zip(
                    matsimnibs_list, didt_list, output_names, geo_names):
                sims.append(
                    pool.apply_async(
                        _run_tms,
                        (mesh, cond, cond_list, fn_coil, fields,
                         matsimnibs, didt, fn_out, fn_geo)))
            pool.close()
            pool.join()


def _set_up_global_solver(S):
    global tms_global_solver
    tms_global_solver = S


def _run_tms(mesh, cond, cond_list, fn_coil, fields, matsimnibs, didt, fn_out, fn_geo):
    global tms_global_solver
    logger.info('Calculating dA/dt field')
    start = time.time()

    dAdt = _get_da_dt_from_coil(fn_coil, mesh, didt, matsimnibs)

    #dAdt = coil_lib.set_up_tms(mesh, fn_coil, matsimnibs, didt,
    #                           fn_geo=fn_geo, fn_stl=fn_stl)
    logger.info(f'{time.time() - start:.2f}s to calculate dA/dt')
    b = tms_global_solver.assemble_rhs(dAdt)
    v = tms_global_solver.solve(b)

    v = mesh_io.NodeData(v, name='v', mesh=mesh)
    v.mesh = mesh
    out = calc_fields(v, fields, cond=cond, dadt=dAdt)
    mesh_io.write_msh(out, fn_out)

    if fn_geo is not None:
        logger.info('Creating visualizations')
        #summary = ''
        skin_mesh = mesh.crop_mesh(tags = [ElementTags.SCALP_TH_SURFACE])

        # write .opt-file
        v = out.view(
            visible_tags=[ElementTags.GM_TH_SURFACE.value],
            visible_fields=['magnE'],
            cond_list=cond_list)
        TmsCoil.from_file(fn_coil).append_simulation_visualization(v, fn_geo, skin_mesh, matsimnibs)

        mesh_io.write_geo_triangles(skin_mesh.elm.node_number_list - 1,
                                        skin_mesh.nodes.node_coord, fn_geo,
                                        name='scalp', mode='ba')
        v.add_view(ColormapNumber=8, ColormapAlpha=.3,
                Visible=0, ShowScale=0)  # scalp
        v.add_merge(fn_geo)
        v.write_opt(fn_out)

    #if view:
    #    mesh_io.open_in_gmsh(fn_out, True)
#
    #summary += f'\n{os.path.split(s)[1][:-1]}\n'
    #summary += len(os.path.split(s)[1][:-1]) * '=' + '\n'
    #summary += 'Gray Matter\n\n'
    #summary += m.fields_summary(roi=2)
#
    #logger.log(25, summary)
    
    del dAdt, v, b
    gc.collect()

def _get_da_dt_from_coil(fn_coil, mesh, didt, matsimnibs):
    tms_coil = TmsCoil.from_file(fn_coil)

    didt = np.atleast_1d(didt)
    if len(didt) == 1:
        for stimulator in tms_coil.get_elements_grouped_by_stimulators().keys():
            stimulator.di_dt = didt
    else:
        for stimulator, stimulator_didt in zip(tms_coil.get_elements_grouped_by_stimulators().keys(), didt):
            stimulator.di_dt = stimulator_didt
    return tms_coil.get_da_dt(mesh, matsimnibs).node_data2elm_data()

def _finalize_global_solver():
    global tms_global_solver
    del tms_global_solver
    gc.collect()


def tdcs_neumann(mesh, cond, currents, electrode_surface_tags):
    '''Simulates a tDCS electric potential using Neumann boundary conditions on
    the electrodes.

    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh file with geometry information
    cond: simnibs.msh.mesh_io.ElementData
        An ElementData field with conductivity information
    currents: list or ndarray
        A list of currents going though each electrode
    electrode_surface_tags: list
        A list of the indices of the surfaces where the dirichlet BC is to be applied

    Returns
    -------
    potential: simnibs.msh.mesh_io.NodeData
        Total electric potential
    '''
    assert len(electrode_surface_tags) == len(currents),\
        'Please define one current per electrode'
    assert np.isclose(np.sum(currents), 0.), 'Sum of currents must be zero'
    S = TDCSFEMNeumann(mesh, cond, electrode_surface_tags[0])
    b = S.assemble_rhs(electrode_surface_tags[1:], currents[1:])
    v = S.solve(b)
    
    del S, b
    gc.collect()
    return mesh_io.NodeData(v, name='v', mesh=mesh)


def tdcs_leadfield(mesh, cond, electrode_surface, fn_hdf5, dataset,
                   current=1., roi=None, post_pro=None, field='E',
                   solver_options=None, n_workers=1, input_type='tag',
                   weigh_by_area=True):
    '''Simulates tDCS fields using Neumann boundary conditions and writes the
    output electric fields to an HDF5 file.

    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh file with geometry information
    cond: simnibs.msh.mesh_io.ElementData
        An ElementData field with conductivity information
    electrode_surface: list
        If `input_type` is "tag", then this is a list of the tags of the
        electrode surfaces. If `input_type` is "nodes", then this is a list
        (or list of lists) of nodes. The first electrode is used as a
        reference.
    fn_hdf5: str
        Name of hdf5 where simulations will be saved
    dataset: str
        Name of dataset where data is to be saved
    current: float | iterable
        Specifies the current (in ampere) to use in each simulation. If float,
        this current will be used in all simulations. If iterable (where the
        length of the iterable is number of electrodes minus one as the
        reference is set to zero), it explicitly specifies the current in each
        simulation. (Default = 1).
    roi: list or None (optional)
        Regions of interest where the fields is to be saved.
        If set to None, will save the electric field in all tissues.
        Default: None
    field: 'E' or 'J' (optional)
        Which field to save (electric field E or current density J). Default: 'E'
    post_pro: callable (optional)
        callable f_post = post_pro(f), where f is an input field in the ROI and
        f_post is an Nx3 ndarray. The postprocessing result will be saved instead of the
        field
    solver_options: str (optional)
        Options to be used by the solver. Default: Hypre solver
    n_workers: int
        Number of workers to use
    input_type: 'tag'  or "nodes" (optional)
        Whether electrode_surface refers to surface tags (default) or nodes.
    weigh_by_area: bool
        Weigh current by node area. If `input_type == "tag"` this is ignored
        and area weighting is implied.

    Returns
    -------
    None
        Writes the field resulting from each simulation to a dataset called
        fn_dataset in an hdf5 file called fn_hdf5.

    Notes
    -----
    Possible combinations/uses of `electrode_surface` and `current` when
    `input_type="nodes"`.

    (1) One node per electrode (weigh_by_area has no effect)

        tdcs_leadfield(..., el=[1,2,3,4], current=1., input_type="nodes", ...)
        tdcs_leadfield(..., el=[1,2,3,4], current=[1.,2.,1.], input_type="nodes", ...)

    (2) Several nodes per electrode

    Set same current for all electrodes (or one value per electrode) and weigh
    according to area so that the total current per electrode is equal to the
    input current

        tdcs_leadfield(..., el=[[1,2,3],[4,5,6],[7,8,9]], current=1., input_type="nodes", ...)
        tdcs_leadfield(..., el=[[1,2,3],[4,5,6],[7,8,9]], current=[1., 2.], input_type="nodes", ...)

    Set the weighting eplicitly by specifying the current per node for all
    electrodes (disable weigh_by_area)

        tdcs_leadfield(...,
            el = [[1,2,3],[4,5,6],[7,8,9]],
            current = [[0.1,0.2,0.7],[0.2,0.3,0.5]],
            input_type = "nodes",
            weigh_by_area = False,
        ...)

    '''
    if field not in ('E', 'J'):
        raise ValueError(f"Field shoud be either 'E' or 'J' (got {field})")

    # Construct system and gradient matrix
    S = TDCSFEMNeumann(
        mesh,
        cond,
        electrode_surface[0],
        input_type,
        weigh_by_area,
        solver_options,
    )

    logger.info("Computing gradient matrix")
    D = grad_matrix(mesh, split=True)
    n_out = mesh.elm.nr
    # Separate out the part of the gradiend that is in the ROI
    if roi is not None:
        roi = np.in1d(mesh.elm.tag1, roi)
        D = [d.tocsc() for d in D]
        D = [d[roi] for d in D]
        n_out = np.sum(roi)
        cond_roi = cond.value[roi]

    # Figure out size of the postprocessing output
    if post_pro is not None:
        n_out = len(post_pro(np.zeros((n_out, 3))))

    # Create HDF5 dataset
    with h5py.File(fn_hdf5, 'a') as f:
        f.create_dataset(
            dataset,
            (len(electrode_surface) - 1, n_out, 3),
            dtype=float, compression="gzip")

    n_sims = len(electrode_surface) - 1
    currents = [current]*n_sims if isinstance(current, float) else current
    assert len(currents) == n_sims, f"Number of currents ({len(currents)}) do not correspond to the number of simulations ({n_sims})"

    # Run simulations (sequential)
    if n_workers == 1:
        for i, (el_tag, current) in enumerate(zip(electrode_surface[1:], currents)):
            logger.info('Running Simulation {0} out of {1}'.format(
                i+1, n_sims))
            b = S.assemble_rhs([el_tag], [current])
            v = S.solve(b)

            #TODO implement calibration error also for element/node defined electrodes
            # when input_type == "nodes"
            if input_type == "tag":
                # estimate calibration error
                ref_electrode = el_tag
                # other_electrodes = [x for x in electrode_surface if np.all(x!=ref_electrode)][0]  
                other_electrodes = np.array([x for x in electrode_surface if x!=ref_electrode])             

                v_ = mesh_io.NodeData(v, name='v', mesh=mesh)
                flux = np.array([
                    _calc_flux_electrodes(v_, cond,
                                        [other_electrodes - 1000, other_electrodes - 600,
                                        other_electrodes - 2000, other_electrodes - 1600],
                                        units='mm'),
                    _calc_flux_electrodes(v_, cond,
                                        [ref_electrode - 1000, ref_electrode - 600,
                                        ref_electrode - 2000, ref_electrode - 1600],
                                        units='mm')])
                current_ = np.average(np.abs(flux))
                error = np.abs(np.abs(flux[0]) - np.abs(flux[1])) / current_
                if error > 0.1:
                    logger.warning(f'The current calibration error exceeded 10%! Estimated error value: {error*100:.2f}%')

            E = np.vstack([-d.dot(v) for d in D]).T * 1e3
            if field == 'E':
                out_field = E
            elif field == 'J':
                out_field = calc_J(E, cond_roi)
            else:
                raise ValueError
            if post_pro is not None:
                out_field = post_pro(out_field)
            with h5py.File(fn_hdf5, 'a') as f:
                f[dataset][i] = out_field

        del S, b, v
        gc.collect()
        
    # Run simulations (parallel)
    else:
        # Lock has to be passed through inheritance
        S.lock = multiprocessing.Lock()
        with multiprocessing.Pool(processes=n_workers,
                                  initializer=_set_up_tdcs_global_solver,
                                  initargs=(S, n_sims, D, post_pro, cond_roi, field)) as pool:
            sims = []
            for i, (el_tag, current) in enumerate(zip(electrode_surface[1:], currents)):
                if input_type == "tag":
                    ref_electrode = el_tag
                    other_electrodes = np.array([x for x in electrode_surface if x!=ref_electrode])
                else:
                    ref_electrode = el_tag
                    other_electrodes = [x for x in electrode_surface if np.all(x!=ref_electrode)][0]
                sims.append(
                    pool.apply_async(
                        _run_tdcs_leadfield,
                        (i, [el_tag], [current], fn_hdf5, dataset, input_type, mesh, cond, ref_electrode, other_electrodes)))
            [s.get() for s in sims]
            pool.close()
            pool.join()


# ### Functions for running tDCS leadfields in parallel ####
def _set_up_tdcs_global_solver(S, n, D, post_pro, cond, field):
    global tdcs_global_solver
    global tdcs_global_nsims
    global tdcs_global_grad_matrix
    global tdcs_global_post_pro
    global tdcs_global_cond
    global tdcs_global_field
    tdcs_global_solver = S
    tdcs_global_nsims = n
    tdcs_global_grad_matrix = D
    tdcs_global_post_pro = post_pro
    tdcs_global_cond = cond
    tdcs_global_field = field


def _run_tdcs_leadfield(i, el_tags, currents, fn_hdf5, dataset, input_type, mesh, cond, ref_electrode, other_electrodes):
    global tdcs_global_solver
    global tdcs_global_nsims
    global tdcs_global_grad_matrix
    global tdcs_global_post_pro
    global tdcs_global_cond
    global tdcs_global_field
    logger.info('Running Simulation {0} out of {1}'.format(
        i+1, tdcs_global_nsims))
    b = tdcs_global_solver.assemble_rhs(el_tags, currents)
    v = tdcs_global_solver.solve(b)

    #TODO implement calibration error also for element/node defined electrodes
    # when input_type == "nodes"
    if input_type == "tag":
        v_ = mesh_io.NodeData(v, name='v', mesh=mesh)
        flux = np.array([
            _calc_flux_electrodes(v_, cond,
                                [other_electrodes - 1000, other_electrodes - 600,
                                other_electrodes - 2000, other_electrodes - 1600],
                                units='mm'),
            _calc_flux_electrodes(v_, cond,
                                [ref_electrode - 1000, ref_electrode - 600,
                                ref_electrode - 2000, ref_electrode - 1600],
                                units='mm')])
        current_ = np.average(np.abs(flux))
        error = np.abs(np.abs(flux[0]) - np.abs(flux[1])) / current_
        if error > 0.1:
            logger.warning(f'The current calibration error exceeded 10%! Estimated error value: {error*100:.2f}%')

    # Calculate E and postprocessing
    E = np.vstack([-d.dot(v) for d in tdcs_global_grad_matrix]).T * 1e3
    if tdcs_global_field == 'E':
        out_field = E
    elif tdcs_global_field == 'J':
        out_field = calc_J(E, tdcs_global_cond)
    else:
        raise ValueError

    if tdcs_global_post_pro is not None:
        out_field = tdcs_global_post_pro(out_field)
    # Write out
    tdcs_global_solver.lock.acquire()
    with h5py.File(fn_hdf5, 'a') as f:
        f[dataset][i] = out_field
    tdcs_global_solver.lock.release()
    
    del b, v
    gc.collect()


def _finalize_tdcs_global_solver():
    global tdcs_global_solver
    del tdcs_global_solver
    global tdcs_global_nsims
    del tdcs_global_nsims
    global tdcs_global_grad_matrix
    del tdcs_global_grad_matrix
    global tdcs_global_post_pro
    del tdcs_global_post_pro
    gc.collect()


def tms_many_simulations(
    mesh, cond, fn_coil, matsimnibs_list, didt_list,
    fn_hdf5, dataset, roi=None, field='E', post_pro=None,
    solver_options=None, n_workers=1):
    ''' Function for running a large amount of TMS simulations.

    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh structure
    cond: simnibs.msh.mesh_io.ElementData
        Conductivity field
    fn_coil: string
        Relative path of coil file
    matsimnibs_list: list
        List of "matsimnibs" matrices, one per position
    didt_list: list
        List of dIdt values, one per position
    fn_hdf5: str
        Name of hdf5 where simulations will be saved
    dataset: str
        Name of dataset where data is to be saved
    roi: list or None (optional)
        Regions of interest where the fields is to be saved.
        If set to None, will save the electric field in all tissues.
        Default: None
    field: str
        Which field/s to save, any combination of E, D, J, v. E.g. 'EDJ'. Default: 'E'.
        Note: When post-processing is used the field can only be either 'E' or 'J'
    post_pro: list of callables (optional)
        List of callables f_post = post_pro(f), where f is the requested field (or multiple field as a tuple of
        the ordering E,D,J,v) in the ROI and f_post is an Nx3 ndarray. The postprocessing result will be saved
        instead of the field/s.
    solver_options: str (optional)
        Options to be used by the solver. Default: Hypre solver
    n_workers: int
        Number of workers to use
    '''
    for f in field:
        if f not in 'EDJv':
            raise ValueError("Field must be one or more of 'E', 'D', 'J', 'v'")
    if len(matsimnibs_list) != len(didt_list):
        raise ValueError("matsimnibs_list and didt_list should have the same length")
    D = grad_matrix(mesh, split=True)
    S = TMSFEM(mesh, cond, solver_options)
    n_out = mesh.elm.nr
    # Separate out the part of the gradient that is in the ROI
    if roi is not None:
        roi = np.in1d(mesh.elm.tag1, roi)
        D = [d.tocsc() for d in D]
        D = [d[roi] for d in D]
        cond = cond.value[roi]
    else:
        roi = np.ones(mesh.elm.nr, dtype=bool)

    n_roi = np.sum(roi)
    # Figure out size of the postprocessing output
    if post_pro is not None:
        if len(field) != 1:
            raise ValueError("Post-processing work only with a single field.")
        if field not in 'EJ':
            raise ValueError("Field has to be E or J for post-processing.")
        n_out = np.array(post_pro(np.zeros((n_roi, 3)))).shape
    else:
        n_out = (n_roi, 3)

    n_sims = len(matsimnibs_list)
    # Create HDF5 dataset
    with h5py.File(fn_hdf5, 'a') as f:
        f.create_dataset(
            dataset,
            (n_sims,) + n_out,
            dtype=float, compression="gzip")

    # Run sequentially
    if n_workers == 1:
        for i, matsimnibs, didt in zip(range(n_sims), matsimnibs_list, didt_list):
            logger.info(
                f'Running Simulation {i+1} out of {n_sims}')
            dAdt = _get_da_dt_from_coil(fn_coil, mesh, didt, matsimnibs)
            # b = S.assemble_tms_rhs(dAdt)
            b = S.assemble_rhs(dAdt)
            v = S.solve(b)
            E = np.vstack([-d.dot(v) for d in D]).T * 1e3
            dAdt = dAdt[roi]
            E -= dAdt

            # build output fields
            out_field = []
            if 'E' in field:
                out_field.append(E)
            if 'D' in field:
                out_field.append(dAdt)
            if 'J' in field:
                out_field.append(calc_J(E, cond))
            if 'v' in field:
                out_field.append(v)
            out_field = tuple(out_field)

            # if only one field to output, un-tuple
            if len(out_field) == 1:
                out_field = out_field[0]
            if post_pro is not None:
                out_field = post_pro(out_field)
            with h5py.File(fn_hdf5, 'a') as f:
                f[dataset][i] = out_field

            del b
            gc.collect()
            
        del S
        gc.collect()

    # Run in parallel
    else:
        # Lock has to be passed through inheritance
        S.lock = multiprocessing.Lock()
        with multiprocessing.Pool(
                processes=n_workers,
                initializer=_set_up_tms_many_global_solver,
                initargs=(S, fn_coil, n_sims, D, post_pro, cond, field, roi)) as pool:
            sims = []
            for i, matsimnibs, didt in zip(
                    range(n_sims), matsimnibs_list, didt_list):
                sims.append(
                    pool.apply_async(
                        _run_tms_many_simulations,
                        (i, matsimnibs, didt,
                         fn_hdf5, dataset)))
            [s.get() for s in sims]
            pool.close()
            pool.join()


### Functions for running man TMS simulations in parallel ####
def _set_up_tms_many_global_solver(S, fn_coil, n, D, post_pro, cond, field, roi):
    global tms_many_global_solver
    global tms_many_global_fn_coil
    global tms_many_global_nsims
    global tms_many_global_grad_matrix
    global tms_many_global_post_pro
    global tms_many_global_cond
    global tms_many_global_field
    global tms_many_global_roi
    tms_many_global_solver = S
    tms_many_global_fn_coil = fn_coil
    tms_many_global_nsims = n
    tms_many_global_grad_matrix = D
    tms_many_global_post_pro = post_pro
    tms_many_global_cond = cond
    tms_many_global_field = field
    tms_many_global_roi = roi


def _run_tms_many_simulations(i, matsimnibs, didt, fn_hdf5, dataset):
    global tms_many_global_solver
    global tms_many_global_fn_coil
    global tms_many_global_nsims
    global tms_many_global_grad_matrix
    global tms_many_global_post_pro
    global tms_many_global_cond
    global tms_many_global_field
    global tms_many_global_roi
    logger.info('Running Simulation {0} out of {1}'.format(
        i+1, tms_many_global_nsims))
    # RHS
    dAdt = _get_da_dt_from_coil(
        tms_many_global_fn_coil,
        tms_many_global_solver.mesh,
        didt,
        matsimnibs
    )
    b = tms_many_global_solver.assemble_rhs(dAdt)
    # Simulate
    v = tms_many_global_solver.solve(b)
    # Calculate E and postprocessing
    E = np.vstack([-d.dot(v) for d in tms_many_global_grad_matrix]).T * 1e3
    E -= dAdt[tms_many_global_roi]

    # build output fields
    out_field = []
    if 'E' in tms_many_global_field:
        out_field.append(E)
    if 'D' in tms_many_global_field:
        out_field.append(dAdt)
    if 'J' in tms_many_global_field:
        out_field.append(calc_J(E, tms_many_global_cond))
    if 'v' in tms_many_global_field:
        out_field.append(v)
    out_field = tuple(out_field)

    # if only one field to output, un-tuple
    if len(out_field) == 1:
        out_field = out_field[0]

    if tms_many_global_post_pro is not None:
        out_field = tms_many_global_post_pro(out_field)
    # Write out
    tms_many_global_solver.lock.acquire()
    with h5py.File(fn_hdf5, 'a') as f:
        f[dataset][i] = out_field
    tms_many_global_solver.lock.release()
    
    del b
    gc.collect()



def _finalize_tms_many_simulations_global_solver():
    global tms_many_global_solver
    global tms_many_global_fn_coil
    global tms_many_global_nsims
    global tms_many_global_grad_matrix
    global tms_many_global_post_pro
    global tms_many_global_cond
    global tms_many_global_field
    global tms_many_global_roi

    del tms_many_global_solver
    del tms_many_global_fn_coil
    del tms_many_global_nsims
    del tms_many_global_grad_matrix
    del tms_many_global_post_pro
    del tms_many_global_cond
    del tms_many_global_field
    del tms_many_global_roi
    gc.collect()
### Finished function to run many TMS simulations in parallel ####

def electric_dipole(mesh, cond, dipole_positions, dipole_moments, source_model,
                    solver_options=None, units="mm"):
    ''' Electric dipole simulations using the partial integration method

    Parameters
    ----------
    mesh: simnibs.mesh_io.Msh
        Mesh file with geometry information
    cond: simnibs.msh.mesh_io.ElementData
        An ElementData field with conductivity information
    dipole_positions: Nx3 ndarray
        Positions of the dipoles. Each dipole will be a separate simulation
    dipole_moments: Nx3 ndarray
        Moment of each dipole in ampere-meter (Am).
    source_model: str
        Source model to use (partial integration, st. venant).
    solver_options: str (optional)
        Options for the sparse solver. Default: CG + AMG

    Returns
    -------
    v: np.ndarray of size Nxmesh.nodes.nr
        Electric potential caused by each dipole
    '''

    dipole_positions = np.atleast_2d(dipole_positions)
    dipole_moments = np.atleast_2d(dipole_moments)

    assert dipole_positions.shape[1] == 3, 'dipole_positions should be in Nx3 format!'
    assert dipole_moments.shape[1] == 3, 'dipole_moments should be in Nx3 format!'
    if dipole_positions.shape[0] != dipole_moments.shape[0]:
        raise ValueError('Different number of entries for '
                         'dipole_positions and dipole_moments')

    S = DipoleFEM(mesh, cond, solver_options, units)
    b = S.assemble_rhs(dipole_positions, dipole_moments, source_model)
    return np.atleast_2d(S.solve(b).T)


def get_dirichlet_node_index_cog(mesh, roi=None):
    """
    Get closest node to center of gravity of head model on lower 10% quantile in z-direction (neck direction)
    ensuring that it does not lie on a surface. (indexing starting with 1)

    Parameters
    ----------
    mesh : Msh object
        Mesh object
    roi : list of RegionOfInterest instances, optional, default: None
        List of ROI surfaces. The Dirichlet node will be set such it does not lay on it.

    Returns
    -------
    node_idx : int
        Index of the node with node indexing starting with 1
    """

    # center of whole head model but lower 25% quantile of z-axis (into neck direction)
    target_coords = np.mean(mesh.nodes.node_coord, axis=0)
    target_coords[2] = np.quantile(mesh.nodes.node_coord[:, 2], 0.1)

    # get list of node indices, which are closest to center of gravity of mesh
    node_idx = np.argsort(np.linalg.norm(
        mesh.nodes.node_coord - target_coords, axis=1))
    node_idx_triangles = np.unique(mesh.elm.node_number_list[mesh.elm.triangles-1, :][:, :-1]-1)

    if roi is not None:
        if type(roi) is not list:
            roi = [roi]

        for _roi in roi:
            node_idx_triangles = np.hstack((node_idx_triangles, _roi.node_index_list.flatten()))

        node_idx_triangles = np.unique(node_idx_triangles)

    # test and return first node, which is not lying on a surface
    for idx in node_idx:
        if idx not in node_idx_triangles:
            return idx + 1
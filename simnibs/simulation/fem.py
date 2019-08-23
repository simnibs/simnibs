# -*- coding: utf-8 -*-\
'''
    Functions for assembling and solving FEM systems
'''

import multiprocessing
import time
import copy
import warnings
import atexit
import h5py
import numpy as np
import scipy.sparse as sparse

from ..msh import mesh_io
from . import cond as cond_lib
from . import coil_numpy as coil_lib
from ..utils.simnibs_logger import logger
from ..cython_code import petsc_solver

PETSC_OPTIONS = '-ksp_type cg ' +\
                '-ksp_rtol 1e-10 -pc_type hypre ' +\
                '-pc_hypre_type boomeramg ' +\
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

# Initialize PETSc
petsc_solver.petsc_initialize()
# Register the finalizer for PETSc
atexit.register(petsc_solver.petsc_finalize)


def calc_fields(potentials, fields, cond=None, dadt=None, units='mm', E=None):
    ''' Given a mesh and the electric potentials at the nodes,
    calculates the fields

    Parameters
    ------------
    potentials: simnibs.msh.mesh_io.NodeData
        NodeData field with potentials.
        Attention: the mesh property should be set
    fields: Any combination of 'vEeJjsDg'
        Fields to output
        v: electric potential at the nodes
        E: Electric field at the elements
        e: Electric field norm at the elements
        J: Current density at the elements
        j: Current density norm at the elements
        s: Conductivity at the elements
        D: dA/dt at the nodes
        g: gradiet of the potential at the elements
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
    ---------
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
                    e, name='normE', mesh=out_mesh))

        if any(f in ['J', 'j', 's'] for f in fields):
            if cond is None:
                raise ValueError(
                    'Cannot calculate J, j os s field: No conductivity input')
            cond.assign_triangle_values()
            if 's' in fields:
                cond.field_name = 'conductivity'
                cond.mesh = out_mesh
                if cond.nr_comp == 9:
                    out_mesh.elmdata += cond_lib.TensorVisualization(cond, out_mesh)
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
                        j, name='normJ', mesh=mesh))

    return out_mesh


def calc_J(E, cond):
    '''Calculates J

    Parameters
    ------------
    E: ndarray of mesh_io.Data
        Electric field. Nx3 vector

    cond: ndarray or mesh_io.Data
        Conductivity. A Nx1 (scalar) or Nx9 (tensor) vector.

    Returns
    ---------
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
    ---------
    nodes: list
        List of nodes there the BC should be applied

    values: list
        Value at each node

    Parameters
    ---------
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
        Parameters:
        -------
        A: scipy.sparse.csr
            Sparse matrix
        b: numpy array or None
            Righ-hand side. if None, it will return None
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b
        Returns:
        ------
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
        Parameters:
        -------
        A: scipy.sparse.csr
            Sparse matrix
        b: numpy array or None
            Righ-hand side. if None, it will return None
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b
        Returns:
        ------
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

        Parameters:
        -------
        A: scipy.sparse.csr
            Sparse matrix
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b

        Returns:
        ------
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
        Parameters:
        -------
        x: numpy array
            Righ-hand side
        dof_map: dofMap
            Mapping of node indexes to rows and columns in A and b

        Returns:
        ------
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
        ------------
        list_of_bcs: list
            List of DirichletBC objects

        Returns
        ---------
        bc: DirichletBC
            DirichletBC corresponding to a union of the list
        '''
        nodes = np.hstack([bc.nodes for bc in list_of_bcs])
        values = np.hstack([bc.values for bc in list_of_bcs])
        return cls(nodes, values)


class FEMSystem(object):
    ''' Class defining the equation system to be used in FEM calculations

    Parameters
    ------------
    mesh: mesh_io.Msh
       Mesh object
    cond: mesh_io.ElementData
        conductivity information
    dirichlet: list of DirichletBC (optional)
        Dirichlet boundary conditions
    units: {'mm' or 'm'}
        Units of the mesh nodes

    Attributes
    ------------
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

    Note
    ------
    Once created, do NOT change the attributes of this class.

    '''
    def __init__(self, mesh, cond, dirichlet=None, units='mm', store_G=False):
        if units in ['mm', 'm']:
            self.units = units
        else:
            raise ValueError('Invalid unit: {0}'.format(units))
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
        self._petsc_solver = None
        self._G = None # Gradient operator
        self._D = None # Gradient matrix
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
        logger.info('Assembling FEM Matrx')
        start = time.time()
        msh = self.mesh
        cond = self.cond[msh.elm.elm_type == 4]
        th_nodes = msh.elm.node_number_list[msh.elm.elm_type == 4]
        G = _gradient_operator(msh)
        if store_G: self._G = G  # stores the operator in case we need it later (TMS)
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

    def prepare_solver(self, petsc_options=None):
        '''Prepares the object to solve FEM systems

        Parameters
        -------------
        petsc_options: str (optional)
            Options to be used by PETSc. Default:
            '-ksp_type cg -ksp_rtol 1e-10 -pc_type hypre -pc_hypre_type boomeramg '

        Note
        -------
        After running this method, do NOT change any attributes of the class!
        '''
        if petsc_options is None:
            petsc_options = PETSC_OPTIONS

        A = sparse.csc_matrix(self.A, copy=True)
        A.sort_indices()
        dof_map = copy.deepcopy(self.dof_map)
        if self.dirichlet is not None:
            A, dof_map = self.dirichlet.apply_to_matrix(A, dof_map)
        self._A_reduced = A  # We need to save this as PETSc does not copy the vectors
        self._petsc_solver = petsc_solver.Solver(petsc_options, A)

    def solve(self, b=None, petsc_options=None):
        ''' Solves the FEM system
        Calls PETSc to solve the FEM system

        Parameters
        ------------
        b: np.ndarray (Optional):
            Right-hand side. If not set, will assume zeros
        petsc_options: str (optional)
            Options to be used by PETSc. Default:
            '-ksp_type cg -ksp_rtol 1e-10 -pc_type hypre -pc_hypre_type boomeramg '
        Returns
        -----------
        x: ndarray
            array with solution
        Note
        -------
        After running this method, do NOT change any attributes of the class!
        '''
        logger.debug('Solving FEM System')

        if petsc_options is None:
            petsc_options = PETSC_OPTIONS

        logger.debug('PETSc Options: {0}'.format(petsc_options))

        if self._petsc_solver is None:
            self.prepare_solver(petsc_options)
            pass

        if b is None:
            b = np.zeros(self.dof_map.nr, dtype=float)
        else:
            b = np.copy(b)

        # We also need the A matrix here because the DOFs change
        dof_map = copy.deepcopy(self.dof_map)
        if self.dirichlet is not None:
            b, dof_map = self.dirichlet.apply_to_rhs(self.A, b, dof_map)

        x = self._petsc_solver.solve(b)

        if self.dirichlet is not None:
            x, dof_map = self.dirichlet.apply_to_solution(x, dof_map)
        dof_map, x = dof_map.order_like(self.dof_map, array=x)

        return np.squeeze(x)

    @classmethod
    def tms(cls, mesh, cond):
        '''Sets up a TMS problem

        Parameters:
        ------
        mesh: simnibs.mesh_io.msh.Msh
            Mesh structure
        cond: ndarray or simnibs.mesh_io.msh.ElementData
            Conductivity of each element
        Returns:
        ------
        S: FEMSystem
            FEMSystem structure with a ".solve()" command
        '''
        # Add a dirichlet BC to the single lowest node to make the problem SPD
        lowest_node = mesh.nodes.node_number[mesh.nodes.node_coord[:, 2].argmin()]
        bc = DirichletBC([lowest_node], [0])
        S = cls(mesh, cond, dirichlet=bc, store_G=True)
        return S

    def assemble_tms_rhs(self, dadt):
        ''' Assembles the right hand side for TMS.

        Parameters:
        -------------
        dadt: NodeData or ElementData
            dA/dt field at each node or element

        Returns:
        -----------
        b: np.array
            Right-hand side
        '''

        msh = self.mesh
        cond = self.cond[msh.elm.elm_type == 4]
        if isinstance(dadt, mesh_io.NodeData):
            dadt = dadt.node_data2elm_data()
        dadt = dadt.value
        dadt = dadt[msh.elm.elm_type == 4]
        if self._G is None:
            G = _gradient_operator(msh)
        else:
            G = self._G
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

    @classmethod
    def tdcs(cls, mesh, cond, electrode_tags, potentials):
        '''Sets up a tDCS problem using Dirichled boundary conditions

        Parameters:
        ------
        mesh: simnibs.mesh_io.msh.Msh
            Mesh structure
        cond: ndarray or simnibs.mesh_io.msh.ElementData
            Conductivity of each element
        electrode_tags: list
            list of the surfaces where the dirichlet BC is to be applied
        potentials: list
            list of the potentials each surface is to be set
        Returns:
        ------
        S: FEMSystem
            FEMSystem structure with a ".solve()" command
        '''
        assert len(electrode_tags) == len(potentials)
        bcs = []
        for t, p in zip(electrode_tags, potentials):
            elements_in_surface = (mesh.elm.tag1 == t) * (mesh.elm.elm_type == 2)
            if np.sum(elements_in_surface) == 0:
                raise ValueError('Did not find any surface with tag: {0}'.format(t))
            n = np.unique(mesh.elm.node_number_list[elements_in_surface, :3])
            bcs.append(DirichletBC(n, p * np.ones_like(n, dtype=float)))
        bc = DirichletBC.join(bcs)

        return cls(mesh, cond, dirichlet=bc)

    @classmethod
    def tdcs_neumann(cls, mesh, cond, ground_electrode):
        '''Sets up a tDCS problem using Dirichlet boundary conditions in the ground
        electrode and Neumann BC in the others

        Parameters:
        ------
        mesh: simnibs.mesh_io.msh.Msh
            Mesh structure
        cond: ndarray or simnibs.mesh_io.msh.ElementData
            Conductivity of each element
        ground_electrode: int
            Tag of the ground electrode surface
        Returns:
        ------
        S: FEMSystem
            FEMSystem structure with a ".solve()" command
        '''
        # The first surface is set to a DirichletBC
        elements_in_surface = (mesh.elm.tag1 == ground_electrode) * (mesh.elm.elm_type == 2)
        if np.sum(elements_in_surface) == 0:
            raise ValueError(
                'Did not find any surface with tag: {0}'.format(ground_electrode))
        n = np.unique(mesh.elm.node_number_list[elements_in_surface, :3])
        bcs = DirichletBC(n, np.zeros_like(n, dtype=float))
        S = cls(mesh, cond, dirichlet=bcs)
        return S

    def assemble_tdcs_neumann_rhs(self, electrode_tags, currents):
        ''' Assemble the RHS for a tDCS simulation with Neumann boundary conditions 

        Parameters:
        ---------------
        electrode_tags: list
            list of the surfaces where the currents to be applied.
            WARNING: should NOT include the ground electrode
        currents: list
            list of the currents in each surface
            WARNING: should NOT include the ground electrode

        Returns:
        -------------
        b: np.ndarray
            Right-hand-side of FEM system

        '''
        assert len(electrode_tags) == len(currents)
        b = np.zeros(self.dof_map.nr, dtype=np.float64)
        areas = self.mesh.elements_volumes_and_areas().value
        for t, c in zip(electrode_tags, currents):
            t = int(t)
            elements_in_surface = (self.mesh.elm.tag1 == t) * (self.mesh.elm.elm_type == 2)
            if np.sum(elements_in_surface) == 0:
                raise ValueError('Did not find any surface with tag: {0}'.format(t))
            total_area = areas[elements_in_surface].sum()
            # We should not have any conductivity here - doing the math, the conductivity
            # shows up naturally in the Neumann term, so we do not need to add anything to
            # here to the RHS
            flux = c / total_area * areas[elements_in_surface]
            tr_nodes = self.mesh.elm.node_number_list[elements_in_surface, :3]
            # The scaling factors in the area cancel each other out
            for i in range(3):
                b += np.bincount(self.dof_map[tr_nodes[:, i]], flux, self.dof_map.nr) / 3.0
        return b

    def calc_gradient(self, v):
        ''' Calculates gradients

        Parameters
        -----------
        v: np.ndarray
            Array with fields at the nodes. Can be 1d or 2d (n_nodes x n)

        Returns
        ----------
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


def assemble_diagonal_mass_matrix(msh, units='mm'):
    ''' Assemble a Mass matrix by doing a first-order integration at the nodes
    Results in a diagonal matrix

    Parameters
    ---------
    msh: simnibs.msh.mesh_io.Msh
        Mesh structure
    units: {'m' or 'mm'}
        Units where the mesh is defined. The matrix will be scaled accordingly

    Result
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
    ------------
    msh: simnibs.msh.mesh_io
        Mesh structure
    G: sparse matrix (optional)
        G matrix to avoid re-calculations. If not set, it will be re-calculated
    split: bool (optional)
        If true, will return a list of sparse matrices, one for each component.
        Default: False

    Returns
    ---------
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


def tdcs(mesh, cond, currents, electrode_surface_tags, n_workers=1, units='mm'):
    ''' Simulates a tDCS field using PETSc

    Parameters:
    ------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh file with geometry information
    cond: simnibs.msh.mesh_io.ElementData
        An ElementData field with conductivity information
    currents: list or ndarray
        A list of currents going though each electrode
    electrode_surface_tags: list
        A list of the indices of the surfaces where the dirichlet BC is to be applied

    Returns:
    ---------------
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
    total_p = np.zeros(mesh.nodes.nr, dtype=np.float)

    n_workers = min(len(currents) - 1, n_workers)
    if n_workers == 1:
        for el_surf, el_c in zip(electrode_surface_tags[1:], currents[1:]):
            total_p += _sim_tdcs_pair(
                mesh, cond, ref_electrode, el_surf, el_c, units)
    else:
        with multiprocessing.Pool(processes=n_workers) as pool:
            sims = []
            for el_surf, el_c in zip(electrode_surface_tags[1:], currents[1:]):
                sims.append(
                    pool.apply_async(
                        _sim_tdcs_pair,
                        (mesh, cond, ref_electrode, el_surf, el_c, units)))
            for s in sims:
                total_p += s.get()
            pool.close()
            pool.join()

    return mesh_io.NodeData(total_p, 'v', mesh=mesh)


def _sim_tdcs_pair(mesh, cond, ref_electrode, el_surf, el_c, units):
    logger.info('Simulating electrode pair {0} - {1}'.format(
        ref_electrode, el_surf))
    S = FEMSystem.tdcs(mesh, cond, [ref_electrode, el_surf], [0., 1.])
    v = S.solve()
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
    logger.info('Estimated current calibration error: {0:.1%}'.format(error))
    return el_c / current * v.value


def _calc_flux_electrodes(v, cond, el_volume, scalp_tag=[5, 1005], units='mm'):
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


def tms_dadt(mesh, cond, dAdt):
    ''' Simulates a TMS field using PETSc
    If strings are used, it will use already existing files. Otherwise, temporary files
    are created

    Parameters:
    ------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh file with geometry information
    cond: simnibs.msh.mesh_io.ElementData
        An ElementData field with conductivity information
    dAdt: simnibs.msh.mesh_io.NodeData or simnibs.msh.mesh_io.ElementData
        dAdt information
    Returns:
    -------
    v:  simnibs.msh.mesh_io.NodeData
        NodeData instance with potential at the nodes
    '''
    S = FEMSystem.tms(mesh, cond)
    b = S.assemble_tms_rhs(dAdt)
    v = S.solve(b)
    v = mesh_io.NodeData(v, name='v', mesh=mesh)
    return v


def tms_coil(mesh, cond, fn_coil, fields, matsimnibs_list, didt_list,
             output_names, geo_names=None, n_workers=1):
    ''' Simulates TMS fields using a coild + matsimnibs + dIdt definition

    Parameters
    ------------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh structure
    cond: simnibs.msh.mesh_io.ElementData
        Conductivity field
    fields: str
        Fields to be calculated for each position
    fn_coil: string
        Name of coil file
    matsimnibs_list: list
        List of "matsimnibs" matrices, one per position
    didt_list: list
        List of dIdt values, one per position
    output_names: list
        List of output mesh file names, one per position
    geo_names: list
        List of output mesh file names, one per position
    n_workers: int
        Number of workers to use

    Returns
    --------
    Writes output meshes to the files specified in output_names
    '''
    assert len(matsimnibs_list) == len(didt_list)
    assert len(output_names) == len(didt_list)
    n_sims = len(matsimnibs_list)
    n_workers = min(n_sims, n_workers)

    if geo_names is None:
        geo_names = [None for i in range(n_sims)]

    S = FEMSystem.tms(mesh, cond)
    if n_workers == 1:
        _set_up_global_solver(S)
        for matsimnibs, didt, fn_out, fn_geo in zip(
                matsimnibs_list, didt_list, output_names, geo_names):
            _run_tms(
                mesh, cond, fn_coil, fields,
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
                        (mesh, cond, fn_coil, fields,
                         matsimnibs, didt, fn_out, fn_geo)))
            pool.close()
            pool.join()

def _set_up_global_solver(S):
    global tms_global_solver
    tms_global_solver = S

def _run_tms(mesh, cond, fn_coil, fields, matsimnibs, didt, fn_out, fn_geo):
    global tms_global_solver
    logger.info('Calculating dA/dt field')
    dAdt = coil_lib.set_up_tms(mesh, fn_coil, matsimnibs, didt, fn_geo=fn_geo)
    b = tms_global_solver.assemble_tms_rhs(dAdt)
    v = tms_global_solver.solve(b)
    v = mesh_io.NodeData(v, name='v', mesh=mesh)
    v.mesh = mesh
    out = calc_fields(v, fields, cond=cond, dadt=dAdt)
    mesh_io.write_msh(out, fn_out)

def _finalize_global_solver():
    global tms_global_solver
    del tms_global_solver

def tdcs_neumann(mesh, cond, currents, electrode_surface_tags):
    ''' Simulates a tDCS field using PETSc and Neumann boundary conditions on the
    electrodes

    Parameters:
    ------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh file with geometry information
    cond: simnibs.msh.mesh_io.ElementData
        An ElementData field with conductivity information
    currents: list or ndarray
        A list of currents going though each electrode
    electrode_surface_tags: list
        A list of the indices of the surfaces where the dirichlet BC is to be applied

    Returns:
    ---------------
    potential: simnibs.msh.mesh_io.NodeData
        Total electric potential
    '''
    assert len(electrode_surface_tags) == len(currents),\
        'Please define one current per electrode'
    assert np.isclose(np.sum(currents), 0.), 'Sum of currents must be zero'

    S = FEMSystem.tdcs_neumann(mesh, cond, electrode_surface_tags[0])
    b = S.assemble_tdcs_neumann_rhs(electrode_surface_tags[1:], currents[1:])
    v = S.solve(b)
    v = mesh_io.NodeData(v, name='v', mesh=mesh)
    return v


def tdcs_leadfield(mesh, cond, electrode_surface_tags, fn_hdf5, dataset,
                   current=1., roi=None, post_pro=None, field='E', n_workers=1):
    ''' Simulates tDCS fields using Neumann boundary conditions and writes the output
    Electric fields to an HDF5 file

    Parameters
    -----------
    mesh: simnibs.msh.mesh_io.Msh
        Mesh file with geometry information
    cond: simnibs.msh.mesh_io.ElementData
        An ElementData field with conductivity information
    electrode_surface_tags: list
        A list of the indices of the surfaces. The first will be used as a reference
    fn_hdf5: str
        Name of hdf5 where simulations will be saved
    dataset: str
        Name of dataset where data is to be saved
    current: float (optional)
        Current to use in each simulation. Default: 1 A
    roi: list or None (optional)
        Regions of interest where the fields is to be saved.
        If set to None, will save the electric field in all tissues.
        Default: None
    field: 'E' or 'J' (optional)
        Which field to save (electric field E or current density J). Default: 'E'
    post_pro: callable (optional)
        callable f_post = post_pro(f), where f is an input field in the ROI and
        f_post is an Nx3 ndarray
    n_workers: int
        Number of workers to use
    Returns
    --------
    Writes the field resulting from each simulation to a dataset called
    fn_dataset in an hdf5 file called fn_hdf5
    '''
    if field != 'E' and field != 'J':
        raise ValueError("Field shoud be either 'E' or 'J'")
    # Construct system and gradient matrix
    S = FEMSystem.tdcs_neumann(mesh, cond, electrode_surface_tags[0])
    D = grad_matrix(mesh, split=True)
    n_out = mesh.elm.nr
    # Separate out the part of the gradiend that is in the ROI
    if roi is not None:
        roi = np.in1d(mesh.elm.tag1, roi)
        D = [d.tocsc() for d in D]
        D = [d[roi] for d in D]
        n_out = np.sum(roi)
        cond = cond.value[roi]

    # Figure out size of the postprocessing output
    if post_pro is not None:
        n_out = len(post_pro(np.zeros((n_out, 3))))

    # Create HDF5 dataset
    with h5py.File(fn_hdf5, 'a') as f:
        f.create_dataset(
            dataset,
            (len(electrode_surface_tags) - 1, n_out, 3),
            dtype=float, compression="gzip")

    n_sims = len(electrode_surface_tags) - 1
    # Run simulations (sequential)
    if n_workers == 1:
        for i, el_tag in enumerate(electrode_surface_tags[1:]):
            logger.info('Running Simulation {0} out of {1}'.format(
                i+1, n_sims))
            b = S.assemble_tdcs_neumann_rhs([el_tag], [current])
            v = S.solve(b)
            E = np.vstack([-d.dot(v) for d in D]).T * 1e3
            if field == 'E':
                out_field = E
            elif field == 'J':
                out_field = calc_J(E, cond)
            if post_pro is not None:
                out_field = post_pro(out_field)
            with h5py.File(fn_hdf5) as f:
                f[dataset][i] = out_field

    # Run simulations (parallel)
    else:
        # Lock has to be passed through inheritance
        S.lock = multiprocessing.Lock()
        with multiprocessing.Pool(processes=n_workers,
                                  initializer=_set_up_tdcs_global_solver,
                                  initargs=(S, n_sims, D, post_pro, cond, field)) as pool:
            sims = []
            for i, el_tag in enumerate(electrode_surface_tags[1:]):
                sims.append(
                    pool.apply_async(
                        _run_tdcs_leadfield,
                        (i, [el_tag], [current],
                         fn_hdf5, dataset)))
            [s.get() for s in sims]
            pool.close()
            pool.join()

#### Functions for running tDCS leadfields in parallel ####
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


def _run_tdcs_leadfield(i, el_tags, currents, fn_hdf5, dataset):
    global tdcs_global_solver
    global tdcs_global_nsims
    global tdcs_global_grad_matrix
    global tdcs_global_post_pro
    global tdcs_global_cond
    global tdcs_global_field
    logger.info('Running Simulation {0} out of {1}'.format(
        i+1, tdcs_global_nsims))
    # RHS
    b = tdcs_global_solver.assemble_tdcs_neumann_rhs(el_tags, currents)
    # Simulate
    v = tdcs_global_solver.solve(b)
    # Calculate E and postprocessing
    E = np.vstack([-d.dot(v) for d in tdcs_global_grad_matrix]).T * 1e3
    if tdcs_global_field == 'E':
        out_field = E
    elif tdcs_global_field == 'J':
        out_field = calc_J(E, tdcs_global_cond)

    if tdcs_global_post_pro is not None:
        out_field = tdcs_global_post_pro(out_field)
    # Write out
    tdcs_global_solver.lock.acquire()
    with h5py.File(fn_hdf5, 'a') as f:
        f[dataset][i] = out_field
    tdcs_global_solver.lock.release()


def _finalize_tdcs_global_solver():
    global tdcs_global_solver
    del tdcs_global_solver
    global tdcs_global_nsims
    del tdcs_global_nsims
    global tdcs_global_grad_matrix
    del tdcs_global_grad_matrix
    global tdcs_global_post_pro
    del tdcs_global_post_pro

#### Finished functionr to tun tDCS leadfields in parallel ####

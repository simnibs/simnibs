import copy
import os
import sys
from unittest.mock import patch
import tempfile
import pytest
import h5py
import numpy as np
import scipy.sparse.linalg as spalg
import scipy.sparse as sparse

from ... import SIMNIBSDIR
from .. import fem
from .. import analytical_solutions
from ...mesh_tools import mesh_io

@pytest.fixture
def sphere3_msh():
    fn = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn)

@pytest.fixture
def sphere_el_msh():
    fn = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere_w_electrodes.msh')
    return mesh_io.read_msh(fn)


@pytest.fixture
def cube_msh():
    fn = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'cube_w_electrodes.msh')
    return mesh_io.read_msh(fn)


@pytest.fixture
def cube_lr():
    fn = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'cube.msh')
    return mesh_io.read_msh(fn)


@pytest.fixture
def tms_sphere(sphere3_msh):
    m = sphere3_msh.crop_mesh(elm_type=4)
    dipole_pos = np.array([0., 0., 300])
    dipole_moment = np.array([1., 0., 0.])
    didt = 1e6
    r = (m.nodes.node_coord - dipole_pos) * 1e-3
    dAdt = 1e-7 * didt * np.cross(dipole_moment, r) / (np.linalg.norm(r, axis=1)[:, None] ** 3)
    dAdt = mesh_io.NodeData(dAdt, mesh=m)
    dAdt.field_name = 'dAdt'
    dAdt.mesh = m
    pos = m.elements_baricenters().value
    E_analytical = analytical_solutions.tms_E_field(dipole_pos * 1e-3,
                                                    dipole_moment, didt,
                                                    pos * 1e-3)
    cond = mesh_io.ElementData(np.ones(m.elm.nr))
    cond.mesh = m
    return m, cond, dAdt, E_analytical


def rdm(a, b):
    return np.linalg.norm(a / np.linalg.norm(a) -
                          b / np.linalg.norm(b))


def mag(a, b):
    return np.abs(np.log(np.linalg.norm(a) / np.linalg.norm(b)))


def cube_msh_all_nodes_at_tag(m, tag):
    return np.unique(m.elm[m.elm.tag1 == tag, :3])

def cube_msh_center_node_at_tag(m, tag):
    index = np.unique(m.elm[m.elm.tag1 == tag,:3])
    a = np.linalg.norm(m.nodes[index]-m.nodes.node_coord.mean(0), axis=1).argmin()
    return index[a]

def cube_msh_center_tri_nodes_at_tag(m, tag):
    index = m.elm[m.elm.tag1 == tag,:3]
    mesh = m.nodes[index]
    center_tri = np.linalg.norm(mesh.mean(1)-mesh.mean((0,1)), axis=1).argmin()
    return index[center_tri]


class TestCalcJ:
    @pytest.mark.parametrize('d_type', [np.array, mesh_io.ElementData])
    def test_calc_j_scalar(self, d_type):
        E = np.ones((10, 3))
        s = np.arange(10)
        J = fem.calc_J(d_type(E), d_type(s))
        assert np.allclose(J, np.tile(np.arange(10), (3, 1)).T)

    @pytest.mark.parametrize('d_type', [np.array, mesh_io.ElementData])
    def test_calc_j_tensor(self, d_type):
        E = np.ones((10, 3))
        s = np.tile((1, 0, 0, 0, 0.1, 0, 0, 0, 0.01), (10, 1))
        J = fem.calc_J(d_type(E), d_type(s))
        assert np.allclose(J[:, 0], 1)
        assert np.allclose(J[:, 1], 0.1)
        assert np.allclose(J[:, 2], 0.01)


class Testcalc_fields:
    def test_calc_vEJgs(self, sphere3_msh):
        phi = sphere3_msh.nodes.node_coord[:, 0] + \
            2 * sphere3_msh.nodes.node_coord[:, 1] + \
            -3 * sphere3_msh.nodes.node_coord[:, 2]
        potential = mesh_io.NodeData(phi, mesh=sphere3_msh)

        E = np.zeros((sphere3_msh.elm.nr, 3))
        E[:] = [-1., -2., 3.]
        E = mesh_io.ElementData(E, mesh=sphere3_msh)

        cond = sphere3_msh.elm.tag1
        cond = mesh_io.ElementData(cond, mesh=sphere3_msh)

        m = fem.calc_fields(potential, 'vJEgsej', cond)
        assert np.allclose(m.field['v'].value, potential.value)
        assert np.allclose(m.field['E'].value, E.value * 1e3)
        assert np.allclose(m.field['J'].value, cond.value[:, None] * E.value * 1e3)
        assert np.allclose(m.field['g'].value, -E.value * 1e3)
        assert np.allclose(m.field['conductivity'].value, cond.value)
        assert np.allclose(m.field['magnE'].value, np.linalg.norm(E.value, axis=1) * 1e3)
        assert np.allclose(m.field['magnJ'].value,
                           np.linalg.norm(cond.value[:, None] * E.value, axis=1) * 1e3)

    def test_calc_dadt(self, sphere3_msh):
        phi = sphere3_msh.nodes.node_coord[:, 0] + \
            2 * sphere3_msh.nodes.node_coord[:, 1] + \
            -3 * sphere3_msh.nodes.node_coord[:, 2]
        potential = mesh_io.NodeData(phi, mesh=sphere3_msh)

        dadt = .2 * sphere3_msh.nodes.node_coord
        dadt = mesh_io.NodeData(dadt, mesh=sphere3_msh)

        E = np.zeros((sphere3_msh.elm.nr, 3))
        E = [-1, -2, 3] - dadt.node_data2elm_data().value
        E = mesh_io.ElementData(E, mesh=sphere3_msh)
        E.assign_triangle_values()

        cond = sphere3_msh.elm.tag1
        cond = mesh_io.ElementData(cond, mesh=sphere3_msh)

        m = fem.calc_fields(potential, 'vDJEgsej', cond, dadt=dadt,
                            units='m')
        assert np.allclose(m.field['v'].value, potential.value)
        assert np.allclose(m.field['D'].value, dadt.value)
        assert np.allclose(m.field['g'].value, [1, 2, -3])
        assert np.allclose(m.field['E'].value, E.value)
        assert np.allclose(m.field['J'].value, cond.value[:, None] * E.value)
        assert np.allclose(m.field['conductivity'].value, cond.value)
        assert np.allclose(m.field['magnE'].value, np.linalg.norm(E.value, axis=1))
        assert np.allclose(m.field['magnJ'].value,
                           np.linalg.norm(cond.value[:, None] * E.value, axis=1))

    def test_calc_tensor(self, sphere3_msh):
        phi = sphere3_msh.nodes.node_coord[:, 0] + \
            2 * sphere3_msh.nodes.node_coord[:, 1] + \
            -3 * sphere3_msh.nodes.node_coord[:, 2]
        potential = mesh_io.NodeData(phi, mesh=sphere3_msh)

        o = np.ones(sphere3_msh.elm.nr)
        z = np.zeros(sphere3_msh.elm.nr)
        cond = np.reshape(np.eye(3) * np.array([1, 2, 3]), -1)
        cond = np.tile(cond, [sphere3_msh.elm.nr, 1])
        cond = mesh_io.ElementData(cond, mesh=sphere3_msh)
        m = fem.calc_fields(potential, 'vJEgsej', cond, units='m')

        assert np.allclose(m.field['v'].value, potential.value)
        assert np.allclose(m.field['g'].value, [1, 2, -3])
        assert np.allclose(m.field['E'].value, [-1, -2, 3])
        assert np.allclose(m.field['J'].value, [-1, -4, 9])

class TestDOFMap:
    def test_define_dof_map(self):
        dof_map = fem.dofMap(inverse=[1, 2, 5, 3])
        assert dof_map[1] == 0
        assert dof_map[2] == 1
        assert dof_map[5] == 2
        assert dof_map[3] == 3


class TestDirichlet:
    def test_dirichlet_apply(self):
        bc = fem.DirichletBC([1, 3, 5], [2, 4, 6])
        dof_map = fem.dofMap(np.arange(1, 7))
        A = sparse.diags(1 + np.arange(6)).tocsc()
        b = -np.arange(6)
        A, b, dof_map = bc.apply(A, b, dof_map)
        A = A.toarray()
        assert A.shape == (3, 3)
        assert b.shape == (3, 1)
        assert A[0, 0] == 2
        assert A[1, 1] == 4
        assert A[2, 2] == 6
        assert np.all(dof_map.inverse == [2, 4, 6])
        assert dof_map[2] == 0

    def test_add_to_solution(self):
        bc = fem.DirichletBC([1, 3, 5], [1, 3, 5])
        dof_map = fem.dofMap([2, 6, 4])
        x = -np.arange(3)
        x, dof_map = bc.apply_to_solution(x, dof_map)
        assert np.all(dof_map.inverse == [2, 6, 4, 1, 3, 5])
        assert np.allclose(x.T, [0, -1, -2, 1, 3, 5])

    def test_appy_to_vector(self):
        bc = fem.DirichletBC([1, 3, 5], [2, 4, 6])
        dof_map = fem.dofMap(np.arange(0, 6))
        b = np.arange(6)
        b, dof_map = bc.apply_to_vector(b, dof_map)
        assert np.all(dof_map.inverse == [0, 2, 4])
        assert np.allclose(b, [0, 2, 4])

    def test_join(self):
        bc1 = fem.DirichletBC([1], [2])
        bc2 = fem.DirichletBC([3, 5], [4, 6])
        bc = fem.DirichletBC.join([bc1, bc2])
        assert np.all(bc.nodes == [1, 3, 5])
        assert np.all(bc.values == [2, 4, 6])




class TestSolve:
    @pytest.mark.parametrize("ksp_type", ["cg"])
    @pytest.mark.parametrize("pc_type", ["ilu"])
    def test_solve_petsc_diagonal(self, ksp_type, pc_type):
        n = 5
        A = sparse.diags(2 * np.ones(n)).tocsr()
        b = np.ones(n)
        ksp = fem.KSPSolver(A, ksp_type, pc_type)
        x = ksp.solve(b).squeeze()
        assert np.allclose(A.dot(x), b)

    @pytest.mark.parametrize("ksp_type", ["gmres"])
    @pytest.mark.parametrize("pc_type", ["ilu"])
    def test_solve_petsc_random(self, ksp_type, pc_type):
        np.random.seed(0)
        n = 5
        A = np.random.random((5, 5))
        A += A.T
        A = sparse.csr_matrix(A)
        b = np.ones(n)
        ksp = fem.KSPSolver(A, ksp_type, pc_type)
        x = ksp.solve(b).squeeze()
        assert np.allclose(A.dot(x), b)

    @pytest.mark.parametrize("ksp_type", ["gmres"])
    @pytest.mark.parametrize("pc_type", ["ilu"])
    def test_multiple_rhs(self, ksp_type, pc_type):
        np.random.seed(0)
        n = 5
        A = np.random.random((n, n))
        A += A.T
        A = sparse.csr_matrix(A)
        b = np.random.random((n, 3))
        #b = np.ones((n, 3))
        ksp = fem.KSPSolver(A, ksp_type, pc_type)
        x = ksp.solve(b).squeeze()
        assert np.allclose(A.dot(x), b)


class TestAssemble:
    def test_gradient_operator(self, cube_msh):
        G = fem._gradient_operator(cube_msh)
        z = cube_msh.nodes[
            cube_msh.elm.node_number_list[
                cube_msh.elm.elm_type == 4]][:, :, 2]
        grad = np.einsum('aij,ai->aj', G, z)
        assert np.allclose(grad, np.array([0, 0, 1]))

    @pytest.mark.parametrize('split', [False, True])
    def test_grad_matrix(self, split, cube_msh):
        cube_msh = cube_msh.crop_mesh(elm_type=4)
        D = fem.grad_matrix(cube_msh, split=split)
        if not split:
            for i in range(3):
                z = cube_msh.nodes.node_coord[:, i]
                t = [0, 0, 0]
                t[i] = 1
                assert np.allclose(D.dot(z).reshape(-1, 3), t, atol=1e-2)
        if split:
            for i in range(3):
                z = cube_msh.nodes.node_coord[:, i]
                assert np.allclose(D[i].dot(z), 1, atol=1e-2)



    def test_vol(self, sphere3_msh):
        v = fem._vol(sphere3_msh)
        assert np.isclose(np.sum(v), 4./3.*np.pi*(95**3), rtol=1e-2)

    def test_mass_matrix(self, sphere3_msh):
        v = np.ones(sphere3_msh.nodes.nr)
        M = fem.assemble_diagonal_mass_matrix(sphere3_msh)
        assert np.isclose(v.dot(M.dot(v)),
                          1e-9 * 4./3. * np.pi * 95 ** 3,
                          rtol=1e-2)


class TestFEMSystem:
    def test_assemble_fem_matrix(self, sphere3_msh):
        msh = sphere3_msh
        cond = np.ones(msh.elm.nr)
        s = fem.FEMSystem(msh, cond)
        assert np.allclose(s.A.dot(np.ones(s.A.shape[0])), 0)
        assert np.allclose(s.A.dot(np.pi*np.ones(s.A.shape[0])), 0)
        assert np.allclose(s.A.T.toarray(), s.A.toarray())
        cond = np.tile(np.eye(3), (msh.elm.nr, 1, 1))
        s = fem.FEMSystem(msh, cond)
        assert np.allclose(s.A.dot(np.pi*np.ones(s.A.shape[0])), 0)
        assert np.allclose(s.A.T.toarray(), s.A.toarray())

    def test_set_up_tms(self, tms_sphere):
        m, cond, dAdt, E_analytical = tms_sphere
        S = fem.TMSFEM(m, cond)
        b = S.assemble_rhs(dAdt)
        x = spalg.spsolve(S.A, b)
        v = mesh_io.NodeData(x, 'FEM', mesh=m)
        E = -v.gradient().value * 1e3 - dAdt.node_data2elm_data().value
        m.elmdata = [mesh_io.ElementData(E_analytical, 'analytical'),
                     mesh_io.ElementData(E, 'E_FEM'),
                     mesh_io.ElementData(E_analytical + dAdt.node_data2elm_data().value, 'grad_analytical'),
                     dAdt]

        m.nodedata = [mesh_io.NodeData(x, 'FEM')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(E, E_analytical) < .2
        assert np.abs(mag(E, E_analytical)) < np.log(1.1)


    def test_dirichlet_problem_cube(self, cube_lr):
        m = cube_lr.crop_mesh(5)
        cond = np.ones(len(m.elm.tetrahedra))
        top = m.nodes.node_number[m.nodes.node_coord[:, 2] > 49]
        bottom = m.nodes.node_number[m.nodes.node_coord[:, 2] < -49]
        bcs = [fem.DirichletBC(top, np.ones(len(top))),
               fem.DirichletBC(bottom, -np.ones(len(bottom)))]
        S = fem.FEMSystem(m, cond)
        A = S.A
        b = np.zeros(m.nodes.nr)
        dof_map = S.dof_map
        for bc in bcs:
            A, b, dof_map = bc.apply(A, b, dof_map)
        x = spalg.spsolve(A, b)
        for bc in bcs:
            x, dof_map = bc.apply_to_solution(x, dof_map)
        order = dof_map.inverse.argsort()
        x = x[order]
        sol = m.nodes.node_coord[:, 2]/50.
        assert np.allclose(sol, x.T)

    def test_dirichlet_problem_sphere(self, sphere3_msh):
        m = sphere3_msh.crop_mesh(elm_type=4)
        cond = np.ones(len(m.elm.tetrahedra))
        cond[m.elm.tag1 == 4] = .01
        anode = m.nodes.node_number[m.nodes.node_coord[:, 2].argmax()]
        cathode = m.nodes.node_number[m.nodes.node_coord[:, 2].argmin()]
        bcs = [fem.DirichletBC([anode], [1]),
               fem.DirichletBC([cathode], [-1])]
        S = fem.FEMSystem(m, cond)
        A = S.A
        b = np.zeros(m.nodes.nr)
        dof_map = S.dof_map
        for bc in bcs:
            A, b, dof_map = bc.apply(A, b, dof_map)
        x = spalg.spsolve(A, b)
        for bc in bcs:
            x, dof_map = bc.apply_to_solution(x, dof_map)
        order = dof_map.inverse.argsort()
        x = x[order].squeeze()
        v_analytical = analytical_solutions.potential_3layers_surface_electrodes(
            [85, 90, 95], [1., .01, 1.], [0, 0, -95], [0, 0, 95], m.nodes.node_coord)
        v_analytical /= v_analytical[anode - 1]
        v_analytical -= v_analytical[0]
        x -= x[0]
        m.nodedata = [mesh_io.NodeData(v_analytical, 'Analytical'),
                      mesh_io.NodeData(x, 'FEM')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        m = m.crop_mesh(3)
        assert rdm(m.nodedata[0].value, m.nodedata[1].value) < .1

    def test_solve_petsc(self, tms_sphere):
        m, cond, dAdt, E_analytical = tms_sphere
        S = fem.TMSFEM(m, cond)
        b = S.assemble_rhs(dAdt)
        x = S.solve(b)
        v = mesh_io.NodeData(x, 'FEM', mesh=m)
        E = -v.gradient().value * 1e3 - dAdt.node_data2elm_data().value
        m.elmdata = [mesh_io.ElementData(E_analytical, 'analytical'),
                     mesh_io.ElementData(E, 'E_FEM'),
                     mesh_io.ElementData(E_analytical + dAdt.node_data2elm_data().value, 'grad_analytical'),
                     dAdt]

        m.nodedata = [mesh_io.NodeData(x, 'FEM')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(E, E_analytical) < .2
        assert np.abs(mag(E, E_analytical)) < np.log(1.1)

    def test_solve_dirichlet_petsc(self, cube_msh):
        m = cube_msh
        cond = np.ones(m.elm.nr)
        cond[m.elm.tag1 > 5] = 1e3
        cond = mesh_io.ElementData(cond)
        el_tags = [1100, 1101]
        potentials = [1, -1]
        S = fem.TDCSFEMDirichlet(m, cond, el_tags, potentials)
        x = S.solve()
        sol = m.nodes.node_coord[:, 1]/50.
        m.nodedata = [mesh_io.NodeData(x, 'FEM'), mesh_io.NodeData(sol, 'Analytical')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(sol, x.T) < .1
        assert np.abs(mag(x, sol)) < np.log(1.1)

    def test_solve_assemble_neumann_tag(self, cube_msh):
        m = cube_msh
        cond = np.ones(m.elm.nr)
        cond[m.elm.tag1 > 5] = 25
        cond = mesh_io.ElementData(cond)
        el_tags = [1100, 1101]
        currents = [1, -1]
        S = fem.TDCSFEMNeumann(m, cond, el_tags[0])
        b = S.assemble_rhs(el_tags[1:], currents[1:])
        x = S.solve(b)
        sol = (m.nodes.node_coord[:, 1] - 50) / 10
        m.nodedata = [mesh_io.NodeData(x, 'FEM'), mesh_io.NodeData(sol, 'Analytical')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(sol, x.T) < .1
        assert np.abs(mag(x, sol)) < np.log(1.5)

    @pytest.mark.parametrize('current_type', ["float", "list"])
    def test_solve_assemble_neumann_nodes(self, current_type, cube_msh):
        m = cube_msh
        cond = np.ones(m.elm.nr)
        cond[m.elm.tag1 > 5] = 25
        cond = mesh_io.ElementData(cond)

        nodes_top = cube_msh_all_nodes_at_tag(m, 1100)
        nodes_bottom = cube_msh_all_nodes_at_tag(m, 1101)
        electrodes = [nodes_top, nodes_bottom]

        if current_type == "float":
            currents = [1, -1]
            weigh_by_area=True
        elif current_type == "list":
            # this is similar to "float" except that we weigh each node equally
            # whereas they are weighted by area in the other condition
            currents = [np.ones_like(nodes_top)/nodes_top.size, -np.ones_like(nodes_bottom)/nodes_bottom.size]
            weigh_by_area=False
        else:
            raise ValueError
        S = fem.TDCSFEMNeumann(m, cond, electrodes[0], input_type="nodes", weigh_by_area=weigh_by_area)
        b = S.assemble_rhs(electrodes[1:], currents[1:])
        x = S.solve(b)
        sol = (m.nodes.node_coord[:, 1] - 50) / 10
        m.nodedata = [mesh_io.NodeData(x, 'FEM'), mesh_io.NodeData(sol, 'Analytical')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(sol, x.T) < .1
        assert np.abs(mag(x, sol)) < np.log(1.5)


    def test_solve_assemble_aniso(self, cube_msh):
        m = cube_msh
        cond = np.tile(np.eye(3), (m.elm.nr, 1, 1))
        cond[m.elm.tag1 > 5] *= 100
        cond[m.elm.tag1 < 5, 0, :] = 0.01
        cond[m.elm.tag1 < 5, 1, :] = 0.1
        cond = mesh_io.ElementData(cond.reshape(-1, 9))
        el_tags = [1100, 1101]
        currents = [1, -1]
        S = fem.TDCSFEMNeumann(m, cond, el_tags[0])
        b = S.assemble_rhs(el_tags[1:], currents[1:])
        x = S.solve(b)
        sol = (m.nodes.node_coord[:, 1] - 50) / 10
        m.nodedata = [mesh_io.NodeData(x, 'FEM'), mesh_io.NodeData(sol, 'Analytical')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(sol, x.T) < .1
        assert np.abs(mag(x, sol)) < np.log(1.5)

    def test_calc_gradient(self, cube_msh):
        m = cube_msh
        cond = np.ones(m.elm.nr)
        cond = mesh_io.ElementData(cond)
        S = fem.TDCSFEMNeumann(m, cond, 1100)

        z = m.nodes.node_coord[:, 0]
        assert np.allclose(S.calc_gradient(z), [1, 0, 0])
        z = m.nodes.node_coord[:, 1]
        assert np.allclose(S.calc_gradient(z), [0, 1, 0])
        z = m.nodes.node_coord[:, 2]
        assert np.allclose(S.calc_gradient(z), [0, 0, 1])

        z = m.nodes.node_coord
        grad = S.calc_gradient(z) * [2, 3, 1]
        assert np.allclose(grad[:, :, 0], [2, 0, 0])
        assert np.allclose(grad[:, :, 1], [0, 3, 0])
        assert np.allclose(grad[:, :, 2], [0, 0, 1])


class TestTDCS:
    def test_tdcs_petsc(self, cube_msh):
        m = cube_msh
        cond = np.ones(m.elm.nr)
        cond[m.elm.tag1 > 5] = 1e3
        cond = mesh_io.ElementData(cond)
        el_tags = [1100, 1101]
        currents = [1, -1]
        x = fem.tdcs(m, cond, currents, el_tags)
        sol = (m.nodes.node_coord[:, 1] - 50) / 10
        m.nodedata = [x, mesh_io.NodeData(sol, 'Analytical')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(sol, x.value) < .1
        assert np.abs(mag(x.value, sol)) < np.log(1.1)

        currents = [-1, 1]
        x = fem.tdcs(m, cond, currents, el_tags)
        sol = -(m.nodes.node_coord[:, 1] - 50) / 10
        assert rdm(sol, x.value) < .1
        assert np.abs(mag(x.value, sol)) < np.log(1.1)


    def test_tdcs_petsc_3_el(self, cube_msh):
        m = cube_msh
        areas = m.nodes_areas()
        cond = np.ones(m.elm.nr)
        cond[m.elm.tag1 > 5] = 1e3
        cond = mesh_io.ElementData(cond)
        el_tags = [1100, 1101, 1101]
        currents = [.5, -1.5, 1.]
        x = fem.tdcs(m, cond, currents, el_tags)
        sol = (m.nodes.node_coord[:, 1] - 50) / 20
        m.nodedata = [x, mesh_io.NodeData(sol, 'Analytical')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(sol, x.value) < .1
        assert np.abs(mag(x.value, sol)) < np.log(1.1)

    def test_tdcs_petsc_3_el_mp(self, cube_msh):
        m = cube_msh
        cond = np.ones(m.elm.nr)
        cond[m.elm.tag1 > 5] = 1e3
        cond = mesh_io.ElementData(cond)
        el_tags = [1100, 1101, 1101]
        currents = [.5, -1.5, 1.]
        x = fem.tdcs(m, cond, currents, el_tags, n_workers=2)
        sol = (m.nodes.node_coord[:, 1] - 50) / 20
        m.nodedata = [x, mesh_io.NodeData(sol, 'Analytical')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(sol, x.value) < .1
        assert np.abs(mag(x.value, sol)) < np.log(1.1)



class TestTDCSNeumann:
    def test_tdcs_neumann_petsc(self, cube_msh):
        m = cube_msh
        cond = np.ones(m.elm.nr)
        cond[m.elm.tag1 > 5] = 1e3
        cond = mesh_io.ElementData(cond)
        el_tags = [1100, 1101]
        currents = [1, -1]
        x = fem.tdcs_neumann(m, cond, currents, el_tags)
        sol = (m.nodes.node_coord[:, 1] - 50) / 10
        m.nodedata = [x, mesh_io.NodeData(sol, 'Analytical')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(sol, x.value) < .1
        assert np.abs(mag(x.value, sol)) < np.log(1.1)

        currents = [-1, 1]
        x = fem.tdcs_neumann(m, cond, currents, el_tags)
        sol = -(m.nodes.node_coord[:, 1] - 50) / 10
        assert rdm(sol, x.value) < .1
        assert np.abs(mag(x.value, sol)) < np.log(1.1)

    def test_tdcs_neumann_3_el(self, cube_msh):
        m = cube_msh
        cond = np.ones(m.elm.nr)
        cond[m.elm.tag1 > 5] = 1e3
        cond = mesh_io.ElementData(cond)
        el_tags = [1100, 1101, 1101]
        currents = [.5, -1.5, 1.]
        x = fem.tdcs_neumann(m, cond, currents, el_tags)
        sol = (m.nodes.node_coord[:, 1] - 50) / 20
        m.nodedata = [x, mesh_io.NodeData(sol, 'Analytical')]
        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(sol, x.value) < .1
        assert np.abs(mag(x.value, sol)) < np.log(1.1)


class TestTMS:
    def test_tms_dadt(self, tms_sphere):
        m, cond, dAdt, E_analytical = tms_sphere
        v = fem.tms_dadt(m, cond, dAdt)
        E = -v.gradient().value * 1e3 - dAdt.node_data2elm_data().value
        m.elmdata = [mesh_io.ElementData(E_analytical, 'analytical'),
                     mesh_io.ElementData(E, 'E_FEM'),
                     mesh_io.ElementData(E_analytical + dAdt.node_data2elm_data().value, 'grad_analytical'),
                     dAdt]

        #mesh_io.write_msh(m, '~/Tests/fem.msh')
        assert rdm(E, E_analytical) < .2
        assert np.abs(mag(E, E_analytical)) < np.log(1.1)

    @patch.object(fem, '_get_da_dt_from_coil')
    def test_tms_coil(self, mock_set_up, tms_sphere):
        m, cond, dAdt, E_analytical = tms_sphere
        mock_set_up.return_value = dAdt.node_data2elm_data()
        matsimnibs = 'MATSIMNIBS'
        didt = 6
        fn_out = tempfile.NamedTemporaryFile(delete=False).name
        fem.tms_coil(m, cond, None, 'coil.ccd',
                     'EJ', [matsimnibs],
                     [didt], [fn_out])
        E = mesh_io.read_msh(fn_out).field['E'].value
        os.remove(fn_out)
        mock_set_up.assert_called_once_with('coil.ccd', m,
                                            didt, matsimnibs)
        assert rdm(E, E_analytical) < .2
        assert np.abs(mag(E, E_analytical)) < np.log(1.1)

    @patch.object(fem, '_get_da_dt_from_coil')
    def test_tms_coil_parallel(self, mock_set_up, tms_sphere):
        if sys.platform in ['win32', 'darwin']:
            '''Won't run on windows or MacOS because Mock does not work through multiprocessing '''
            assert True
            return
        m, cond, dAdt, E_analytical = tms_sphere
        mock_set_up.return_value = dAdt.node_data2elm_data()
        matsimnibs = 'MATSIMNIBS'
        didt = 6
        fn_out = [tempfile.NamedTemporaryFile(delete=False).name for i in range(4)]
        fem.tms_coil(m, cond, None, 'coil.ccd',
                     'EJ', 4*[matsimnibs],
                     4*[didt], fn_out, n_workers=2)
        for f in fn_out:
            E = mesh_io.read_msh(f).field['E'].value
            os.remove(f)
            assert rdm(E, E_analytical) < .2
            assert np.abs(mag(E, E_analytical)) < np.log(1.1)


class TestLeadfield:
    @pytest.mark.parametrize('post_pro', [False, True])
    @pytest.mark.parametrize('field', ['E', 'J'])
    @pytest.mark.parametrize('n_workers', [1, 2])
    @pytest.mark.parametrize('input_type', ['tag', 'nodes'])
    def test_leadfield(self, input_type, n_workers, field, post_pro, cube_msh):
        if sys.platform in ['win32', 'darwin'] and n_workers > 1:
            ''' Same as above, does not work on windows or MacOS'''
            return

        m = cube_msh
        cond = np.ones(m.elm.nr)
        cond[m.elm.tag1 > 5] = 1e3
        cond = mesh_io.ElementData(cond, mesh=m)
        if input_type == 'tag':
            el = [1100, 1101, 1101]
            current = 1.0
            weigh_by_area = True
        elif input_type == 'nodes':
            # Multiple (three) nodes per electrode, explicit current input, no
            # area weighting (nodes of center triangle)
            el = [cube_msh_center_tri_nodes_at_tag(m, i) for i in [1100, 1101]]
            el.append(el[-1])
            current = 2*[[0.3, 0.3, 0.4]]
            weigh_by_area = False

            # Single node per electrode (center node)
            # el = [cube_msh_center_node_at_tag(m, i) for i in [1100, 1101]]
            # el.append(el[-1])
            # current = 1.
            # weigh_by_area = False

            # Multiple (all) nodes per electrode, single current
            # el = [cube_msh_all_nodes_at_tag(m, i) for i in [1100, 1101]]
            # el.append(el[-1])
            # current = 1.
            # weigh_by_area = True

        fn_hdf5 = tempfile.NamedTemporaryFile(delete=False).name
        dataset = 'leadfield'
        if post_pro:
            def post(E):
                return E[:10]*2
        else:
            post = None

        fem.tdcs_leadfield(
            m, cond, el, fn_hdf5, dataset,
            current=current,
            roi=[5],
            field=field,
            post_pro=post,
            n_workers=n_workers,
            input_type=input_type,
            weigh_by_area = weigh_by_area,
        )

        if not post_pro:
            n_roi = np.sum(m.elm.tag1 == 5)
            with h5py.File(fn_hdf5, 'r') as f:
                assert f[dataset].shape == (2, n_roi, 3)
                assert rdm(f[dataset][0, ...],
                           np.tile([0., 100, 0.], (n_roi, 1))) < .2
                assert mag(f[dataset][0, ...],
                           np.tile([0., 100, 0.], (n_roi, 1))) < np.log(1.1)
        if post_pro:
            with h5py.File(fn_hdf5, 'r') as f:
                assert f[dataset].shape == (2, 10, 3)
                assert rdm(f[dataset][0, ...],
                           np.tile([0., 200, 0.], (10, 1))) < .2
                assert mag(f[dataset][0, ...],
                           np.tile([0., 200, 0.], (10, 1))) < np.log(1.1)

        os.remove(fn_hdf5)


class TestTMSMany:
    @pytest.mark.parametrize('post_pro', [False, True])
    @pytest.mark.parametrize('n_workers', [1, 2])
    @patch.object(fem, '_get_da_dt_from_coil')
    def test_many_simulations(self, mock_set_up, n_workers, post_pro, tms_sphere):
        if sys.platform in ['win32', 'darwin'] and n_workers > 1:
            ''' Same as above, does not work on windows '''
            return
        m, cond, dAdt, E_analytical = tms_sphere
        mock_set_up.return_value = dAdt.node_data2elm_data()
        matsimnibs = np.eye(6)
        didt = 6
        fn_hdf5 = tempfile.NamedTemporaryFile(delete=False).name
        dataset = 'leadfield'
        if post_pro:
            def post(E):
                return E[:10]*2
        else:
            post = None
        fem.tms_many_simulations(
            m, cond, 'coil.ccd',
            2*[matsimnibs], 2*[didt],
            fn_hdf5, dataset, roi=[3],
            post_pro=post,
            n_workers=n_workers
        )
        roi_select = m.elm.tag1 == 3
        with h5py.File(fn_hdf5, 'r') as f:
            for E in f[dataset]:
                if post_pro:
                    assert rdm(E, post(E_analytical[roi_select])) < .3
                    assert mag(E, post(E_analytical[roi_select])) < np.log(1.1)
                else:
                    assert rdm(E, E_analytical[roi_select]) < .3
                    assert mag(E, E_analytical[roi_select]) < np.log(1.1)
        os.remove(fn_hdf5)

class TestDipole:
    # st. venant fails with dipole [80,0,0], [1,0,0]!
    @pytest.mark.parametrize('source_model', ["partial integration"])#, "st. venant"])
    @pytest.mark.parametrize('pos_target', [[50, 0, 0], [80, 0, 0]])
    @pytest.mark.parametrize('dipole_direction', [[0, 1, 0], [1, 0, 0]])
    def test_single_dipole(self, source_model, dipole_direction, pos_target, sphere3_msh):
        bar = sphere3_msh.elements_baricenters()[sphere3_msh.elm.tetrahedra]
        dipole_th = np.argmin(np.linalg.norm(pos_target - bar, axis=1))
        dipole_pos = bar[dipole_th]
        dipole_th = sphere3_msh.elm.tetrahedra[dipole_th]

        surface_nodes = np.unique(sphere3_msh.elm[sphere3_msh.elm.tag1 == 1005, :3])
        surface_nodes_pos = sphere3_msh.nodes[surface_nodes]

        analytical_v = analytical_solutions.potential_dipole_3layers(
            [85, 90, 95], 2, 0.1, dipole_pos, dipole_direction, surface_nodes_pos
        )

        cond = 2 * np.ones(sphere3_msh.elm.nr)
        cond[sphere3_msh.elm.tag1 == 4] = 0.1
        numerical_v = fem.electric_dipole(
            sphere3_msh, cond, pos_target, dipole_direction, source_model,
        )
        numerical_v = numerical_v[0, surface_nodes - 1]

        analytical_v -= np.average(analytical_v)
        numerical_v -= np.average(numerical_v)
        assert rdm(analytical_v, numerical_v) < 0.2
        assert mag(analytical_v, numerical_v) < 0.15

    def test_single_dipole_units(self, sphere3_msh):
        """Test unit specification."""
        source_model = "partial integration"
        dip_pos = np.array([50, 0, 0], dtype=float)
        dip_mom = np.array([0, 1, 0], dtype=float) # always interpreted as m!

        surface_nodes = np.unique(sphere3_msh.elm[sphere3_msh.elm.tag1 == 1005, :3])

        cond = 2 * np.ones(sphere3_msh.elm.nr)
        cond[sphere3_msh.elm.tag1 == 4] = 0.1

        # mesh and dip_pos in mm
        v_mm = fem.electric_dipole(
            sphere3_msh, cond, dip_pos, dip_mom, source_model, units="mm"
        )
        v_mm = v_mm[0, surface_nodes - 1]

        # mesh and dip_pos in m
        dip_pos *= 1e-3
        sphere3_msh_m = copy.deepcopy(sphere3_msh)
        sphere3_msh_m.nodes.node_coord *= 1e-3

        v_m = fem.electric_dipole(
            sphere3_msh_m, cond, dip_pos, dip_mom, source_model, units="m"
        )
        v_m = v_m[0, surface_nodes - 1]

        np.testing.assert_allclose(v_mm, v_m)

    def test_many_dipoles(self, sphere3_msh):
        positions = [[60, 0, 0], [0, 80, 0]]
        moments = [[0, 0, 1], [0, 0.5, 0.5]]
        bar = sphere3_msh.elements_baricenters()[sphere3_msh.elm.tetrahedra]

        surface_nodes = np.unique(sphere3_msh.elm[sphere3_msh.elm.tag1 == 1005, :3])
        surface_nodes_pos = sphere3_msh.nodes[surface_nodes]
        analytical_v = []
        dipole_pos = []

        for p, m in zip(positions, moments):
            dipole_th = np.argmin(np.linalg.norm(p - bar, axis=1))
            dipole_pos.append(bar[dipole_th])
            dipole_th = sphere3_msh.elm.tetrahedra[dipole_th]
            analytical_v.append(
                analytical_solutions.potential_dipole_3layers(
                    [85, 90, 95], 2, 0.1, p, m, surface_nodes_pos
                )
            )

        cond = 2 * np.ones(sphere3_msh.elm.nr)
        cond[sphere3_msh.elm.tag1 == 4] = 0.1
        numerical_v = fem.electric_dipole(
            sphere3_msh, cond, dipole_pos, moments, "partial integration"
        )[:, surface_nodes - 1]

        for a_v, n_v in zip(analytical_v, numerical_v):
            a_v -= np.average(a_v)
            n_v -= np.average(n_v)
            assert rdm(a_v, n_v) < 0.2
            assert mag(a_v, n_v) < 0.15

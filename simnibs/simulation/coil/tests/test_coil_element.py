import os

import numpy as np
import pytest

from simnibs.mesh_tools.mesh_io import Msh
from simnibs.simulation.coil.coil_constants import CoilElementTag
from simnibs.simulation.coil.coil_element import (CoilDipoles, CoilLinePoints,
                                                  CoilSampledGridElements, CoilLineElements)

class TestCalcdAdt:
    def test_calc_dAdt_dipoles(self):   
        dipoles = CoilDipoles(None, None, None, np.array([[0,0,0], [0,0,0]]), np.array([[0,0,1], [0,1,0]]), None)

        target_pos = np.array([[1e3, 0., 0.], [0., 1e3, 0.], [0, 0., 1e3]])
        coil_matrix = np.eye(4)
        da_dt = dipoles.get_da_dt(target_pos, coil_matrix, 1e6)
        np.testing.assert_allclose(da_dt, np.array([[0., 0.1, -0.1], [-0.1, 0, 0], [0.1, 0., 0.]]))

    def test_calc_dAdt_single_line_element(self):   
        dipoles = CoilLineElements(None, None, None, np.array([[0,0,0]]), np.array([[1,0,0]]), None)

        target_pos = np.array([[0, 100, 0]])
        coil_matrix = np.eye(4)
        da_dt = dipoles.get_da_dt(target_pos, coil_matrix, 1)
        print(da_dt)
        np.testing.assert_allclose(da_dt, np.array([[0., 0.1, -0.1], [-0.1, 0, 0], [0.1, 0., 0.]]))

    def test_calc_dAdt_line_points_circle_center(self):   
        angles = np.linspace(0, 2 * np.pi, 500, endpoint=False)
        x = 0 + 10 * np.cos(angles)
        y = 0 + 10 * np.sin(angles)
        z = np.zeros_like(x)
        coordinates = np.column_stack((x, y, z))
        print(coordinates)
        dipoles = CoilLinePoints(None, None, None, coordinates, None)

        target_pos = np.array([[0, 0, 0]])
        coil_matrix = np.eye(4)
        da_dt = dipoles.get_da_dt(target_pos, coil_matrix, 0.1)
        print(da_dt, np.linalg.norm(da_dt[0]))
        np.testing.assert_allclose(da_dt, np.array([[0., 0.1, -0.1], [-0.1, 0, 0], [0.1, 0., 0.]]))

    def test_calc_dAdt_sampled_elements(self, sphere3_msh: Msh):
        affine = np.array([[5., 0., 0., -300],
                           [0., 5., 0., -200],
                           [0., 0., 5., 0.],
                           [0., 0., 0., 1]])

        field = np.ones((121, 81, 41, 3))
        field[..., 1] = 2
        field[..., 2] = 3
        coil_matrix = np.array([[0., 1., 0., 0],
                                [1., 0., 0., 0],
                                [0., 0., 1., -100.],
                                [0., 0., 0., 1]])
        
        sampled_elements = CoilSampledGridElements(None, None, None, field, affine, None)
        da_dt = sampled_elements.get_da_dt(sphere3_msh.nodes.node_coord, coil_matrix, 1e6)
        np.testing.assert_allclose(da_dt[:, 0], 2e6, atol=1e-6)
        np.testing.assert_allclose(da_dt[:, 1], 1e6, atol=1e-6)
        np.testing.assert_allclose(da_dt[:, 2], 3e6, atol=1e-6)

class TestGetMesh:
    def test_dipole_element_mesh(self):   
        dipoles = CoilDipoles(None, None, None, np.array([[0,0,0], [0,0,0]]), np.array([[0,0,1], [0,1,0]]), None)
        mesh = dipoles.get_mesh(np.eye(4), False, False, False, True)

        assert mesh.elm.node_number_list.shape[0] == 2
        np.testing.assert_allclose(np.unique(mesh.elm.elm_type), [1])
        assert np.all(mesh.elm.tag1 < 100)
        assert np.all(mesh.elm.tag1 == CoilElementTag.DIPOLES)
        np.testing.assert_allclose(mesh.nodes.node_coord, [[0,0,0], [0,0,0], [0,0,1], [0,1,0]]) 

    def test_line_element_element_mesh(self):   
        dipoles = CoilLineElements(None, None, None, np.array([[0,0,0], [1, 0, 0], [1, 1, 0]]), np.array([[1,0,0], [0, 1, 0], [-1, -1, 0]]), None)
        mesh = dipoles.get_mesh(np.eye(4), False, False, False, True, 100)

        assert mesh.elm.node_number_list.shape[0] == 3
        np.testing.assert_allclose(np.unique(mesh.elm.elm_type), [1])
        assert np.all(mesh.elm.tag1 < 200)
        assert np.all(mesh.elm.tag1 > 99)
        assert np.all(mesh.elm.tag1 == 100 + CoilElementTag.LINE_ELEMENTS)
        np.testing.assert_allclose(mesh.nodes.node_coord, [[0,0,0], [1, 0, 0], [1, 1, 0], [1, 0, 0], [1, 1, 0], [0,0,0]]) 

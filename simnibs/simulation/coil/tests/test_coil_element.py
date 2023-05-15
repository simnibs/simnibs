import os

import numpy as np
import pytest

from simnibs.mesh_tools.mesh_io import Msh, Nodes
from simnibs.simulation.coil.coil_element import (CoilDipoles,
                                                  CoilSampledElements)

from .... import SIMNIBSDIR


@pytest.fixture(scope='module')
def sphere3_msh():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
    return Msh(fn=fn)

class TestCalcdAdt:
    def test_calc_dAdt_dipoles(self):   
        dipoles = CoilDipoles('', [], np.array([[0,0,0], [0,0,0]]), np.array([[0,0,1], [0,1,0]]), None)

        target_pos = np.array([[1e3, 0., 0.], [0., 1e3, 0.], [0, 0., 1e3]])
        coil_matrix = np.eye(4)
        da_dt = dipoles.get_da_dt(target_pos, coil_matrix, 1e6)
        print(da_dt)
        assert np.allclose(da_dt, np.array([[0., 0.1, -0.1],
                                                      [-0.1, 0, 0],
                                                      [0.1, 0., 0.]]))

    @pytest.mark.skip(reason="no way of currently testing this")
    def test_calc_dAdt_sampled_elements(self, sphere3_msh: Msh):
        affine = np.array([[5., 0., 0., -300],
                           [0., 5., 0., -200],
                           [0., 0., 5., 0.],
                           [0., 0., 0., 1]])
        limits = np.array([[-300, 300],[-200, 200],[0, 200]])
        resolution = np.array([5,5,5])

        dims = [int((max_ - min_) // res) for [min_, max_], res in zip(limits, resolution)]

        dx = np.spacing(1e4)
        x = np.linspace(limits[0][0], limits[0][1] + dx, dims[0]+1)
        y = np.linspace(limits[1][0], limits[1][1] + dx, dims[1]+1)
        z = np.linspace(limits[2][0], limits[2][1] + dx, dims[2]+1)
        points = np.array(np.meshgrid(x, y, z, indexing="ij"))
        points = points.reshape((3, -1)).T

        field = np.ones((121, 81, 41, 3))
        field[..., 1] = 2
        field[..., 2] = 3
        field = field.reshape((3, -1)).T
        coil_matrix = np.array([[0., 1., 0., 0],
                                [1., 0., 0., 0],
                                [0., 0., 1., -100.],
                                [0., 0., 0., 1]])
        
        print(field.shape)
        print(points.shape)
        
        sampled_elements = CoilSampledElements('', [], points, field, None)
        da_dt = sampled_elements.get_da_dt(sphere3_msh.nodes.node_coord, coil_matrix, 1e6)
        assert np.allclose(da_dt[:, 0], 2e6, atol=1e-6)
        assert np.allclose(da_dt[:, 1], 1e6, atol=1e-6)
        assert np.allclose(da_dt[:, 2], 3e6, atol=1e-6)

    
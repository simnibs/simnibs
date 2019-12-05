import os

import numpy as np
import pytest

import simnibs.msh.mesh_io as mesh_io
import simnibs.transformations.surface as t_surface


@pytest.fixture
def sphere3_msh():
    fn = os.path.join(os.path.dirname(os.path.realpath(
        __file__)), '..', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn)


def test_surf2surf(sphere3_msh):
    in_surf = sphere3_msh.crop_mesh(1005)
    in_nodes = in_surf.nodes.node_coord
    in_nodes /= np.average(np.linalg.norm(in_nodes, axis=1))
    field = in_nodes[:, 0]
    out_surf = sphere3_msh.crop_mesh(1004)
    out_nodes = out_surf.nodes.node_coord
    out_nodes /= np.average(np.linalg.norm(out_nodes, axis=1))
    out_field, _ = t_surface._surf2surf(field, in_surf, out_surf)
    assert np.allclose(out_field, out_nodes[:, 0], atol=1e-1)

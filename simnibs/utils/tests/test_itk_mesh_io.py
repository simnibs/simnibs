import os

import numpy as np
import pytest

from ... import SIMNIBSDIR
from .. import itk_mesh_io

@pytest.fixture
def atlas_itk_msh_fn():
    fn = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'cube_atlas', 'atlas.txt.gz')
    return fn

class TestITKReader:
    def test_reader(self, atlas_itk_msh_fn):
        v,t,d = itk_mesh_io.read_itk_tetrahedron(atlas_itk_msh_fn)
        assert(len(v)==111)
        assert(len(t)==539)
        assert(len(d)==len(v))
        assert(np.allclose(np.sum(d,axis=0), (79, 32)))

    def test_itk_to_msh(self, atlas_itk_msh_fn):
        msh = itk_mesh_io.itk_to_msh(atlas_itk_msh_fn)
        assert(np.allclose(msh.nodes.node_coord.max(0),(39,39,39)))
        assert(np.allclose(msh.nodes.node_coord.min(0),(0,0,0)))
        assert(np.all(msh.nodedata[0].value.sum(1)==1))
        
        l, p = msh.nodedata[0].interpolate_to_grid_max(
            np.max(msh.nodes.node_coord,axis=0).astype(int),
            np.identity(4))
        assert np.allclose(np.bincount(l.ravel()),(49984, 9335))
        assert np.allclose(np.bincount(p.ravel().astype(int)),(5652, 53667))
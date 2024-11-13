# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 16:08:27 2022

@author: axthi
"""

import os
import numpy as np
import pytest
import h5py 
from .. import TI_utils as TI
from ...mesh_tools import mesh_io
from ... import SIMNIBSDIR

@pytest.fixture()
def sphere_surf():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn).crop_mesh([1003, 1004])

@pytest.fixture()
def leadfield_surf(sphere_surf):
    _, lf, _ = np.meshgrid(np.arange(sphere_surf.nodes.nr),np.arange(4),np.arange(3))
    return lf.astype(float)

@pytest.fixture()
def fn_surf(sphere_surf, leadfield_surf):
    fn_leadfield = 'tmp_surf_leadfied.hdf5'
    if os.path.isfile(fn_leadfield):
        os.remove(fn_leadfield)
    sphere_surf.write_hdf5(fn_leadfield, 'mesh_leadfield')
    dset = '/mesh_leadfield/leadfields/tdcs_leadfield'
    with h5py.File(fn_leadfield, 'a') as f:
        f.create_dataset(dset, data=leadfield_surf)
        f[dset].attrs['electrode_names'] = ['a','b','c','d']
        f[dset].attrs['reference_electrode'] = 'e'
        
    yield fn_leadfield
    os.remove(fn_leadfield)


def test_load_leadfield(fn_surf):
    leadfield, mesh, idx_lf = TI.load_leadfield(fn_surf)
    assert idx_lf == {'a': 0, 'b': 1, 'c': 2, 'd': 3, 'e': None}

def test_get_field(fn_surf):
    leadfield, mesh, idx_lf = TI.load_leadfield(fn_surf)
    
    ef = TI.get_field(['a','b',1000.],leadfield,idx_lf)
    assert np.all(ef == -1000.)
    
    ef = TI.get_field(['c','d',1],leadfield,idx_lf)
    assert np.all(ef == -1.)
    
    ef = TI.get_field(['c','e',1],leadfield,idx_lf)
    assert np.all(ef == 2.)
    
    ef = TI.get_field(['e','c',1],leadfield,idx_lf)
    assert np.all(ef == -2.)


@pytest.mark.filterwarnings('ignore::RuntimeWarning')
def test_get_maxTI():
    ef1 = np.array([[1,0,0],[1,0,0],[1,0,0],[1,0,0]])
    ef2 = np.array([[1,0,0],[-2,0,0],[0,-1000,0],[0,0,3000]])
    
    TImax = TI.get_maxTI(ef1,ef2)
    assert np.all(np.isclose(TImax,2.))

def test_get_dirTI():
    ef1 = np.array([[1,0,0],[1,0,0],[1,0,0],[1,0,0]])
    ef2 = np.array([[1,0,0],[-2,0,0],[0,-1000,0],[0,0,3000]])
    
    TIamp = TI.get_dirTI(ef1,ef2,[1,0,0])
    assert np.all(np.isclose(TIamp[0:2],2.))
    assert np.all(np.isclose(TIamp[2:4],0.))
    
    TIamp = TI.get_dirTI(ef1,ef2,np.array([1,0,0]))
    assert np.all(np.isclose(TIamp[0:2],2.))
    assert np.all(np.isclose(TIamp[2:4],0.))
    
    TIamp = TI.get_dirTI(ef1,ef2,ef1)
    assert np.all(np.isclose(TIamp[0:2],2.))
    assert np.all(np.isclose(TIamp[2:4],0.))


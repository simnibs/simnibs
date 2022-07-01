# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 17:41:21 2022

@author: axthi
"""

import h5py
import numpy as np

from ..mesh_tools import mesh_io


def load_leadfield(leadfield_hdf, 
                   leadfield_path = '/mesh_leadfield/leadfields/tdcs_leadfield',
                   mesh_path = '/mesh_leadfield/'):
    """
    load leadfield, mesh on which leadfield was calculated and mapping from 
    electrode names to index in the leadfield

    Parameters
    ----------
    leadfield_hdf : string
        path to the leadfield hdf5 file
    leadfield_path : string, optional
        path inside the hdf5 file to the leadfield. 
        The default is '/mesh_leadfield/leadfields/tdcs_leadfield'.
    mesh_path : string, optional
        path inside the hdf5 file to the mesh.
        The default is '/mesh_leadfield/'.

    Returns
    -------
    leadfield : np.ndarray
        Leadfield matrix (N_elec -1 x M x 3) where M is either the number of 
        nodes (for surface-based leadfields) or the number of elements 
        (tet-based leadfields) in the mesh.
    mesh : simnibs.msh.mesh_io.Msh
        Mesh on which the leadfield was calculated
    idx_lf : dict
        Mapping from electrode name to index in the leadfield matrix. The 
        reference electrode has None as index
    """
    with h5py.File(leadfield_hdf, 'r') as f:
        lf_struct = f[leadfield_path]
        leadfield = lf_struct[:] # elecs x mesh nodes x 3
        
        # make a dict: elec name --> index in leadfield matrix
        name_elecs = lf_struct.attrs.get('electrode_names')
        name_ref = lf_struct.attrs.get('reference_electrode')
        name_elecs = name_elecs[name_elecs != name_ref]
        assert leadfield.shape[0] == len(name_elecs) 
        idx_lf = dict(zip(name_elecs,range(len(name_elecs))))
        idx_lf[name_ref] = None
        
    mesh = mesh_io.Msh.read_hdf5(leadfield_hdf, mesh_path) # usually the two hemispheres and eyes
    
    return leadfield, mesh, idx_lf
 

def get_field(elec_pair,leadfield,idx_lf):
    """
    return field for a specific electrode pair and current intensity

    Parameters
    ----------
    elec_pair : list [elec_1,elec_2,current_intensity]
        specifies an electrode pair and the current intensity into first electrode
        (current in second electrode: -current_intensity)
    leadfield : np.ndarray
         Leadfield matrix (N_elec -1 x M x 3) where M is either the number of 
         nodes (for surface-based leadfields) or the number of elements 
         (tet-based leadfields) in the mesh.
    idx_lf : dict
        Mapping from electrode name to index in the leadfield matrix. The 
        reference electrode has None as index

    Returns
    -------
    np.ndarray
       field matrix (M x 3) where M is either the number of 
       nodes (for surface-based leadfields) or the number of elements 
       (tet-based leadfields) in the mesh.

    """
    assert elec_pair[0] != elec_pair[1]
    if idx_lf[elec_pair[0]] is None:
        return -elec_pair[2]*leadfield[ idx_lf[elec_pair[1]] ]
    if idx_lf[elec_pair[1]] is None:
        return  elec_pair[2]*leadfield[ idx_lf[elec_pair[0]] ]
    return elec_pair[2]*(leadfield[ idx_lf[elec_pair[0]] ]
                         - leadfield[ idx_lf[elec_pair[1]] ])


def get_maxTI(E1_org,E2_org):
    """
    calculates the maximal modulation amplitude of the TI envelope using 
    the equation given in Grossman et al, Cell 169, 1029–1041.e6, 2017

    Parameters
    ----------
    E1 : np.ndarray
           field of electrode pair 1 (N x 3) where N is the number of 
           positions at which the field was calculated
    E2 : np.ndarray
        field of electrode pair 2 (N x 3) 

    Returns
    -------
    TImax : np.ndarray (N,)
        maximal modulation amplitude

    """
    assert E1_org.shape == E2_org.shape
    assert E1_org.shape[1] == 3
    E1 = E1_org.copy()
    E2 = E2_org.copy()
    
    # ensure E1>E2
    idx = np.linalg.norm(E2, axis=1) > np.linalg.norm(E1, axis=1)
    E1[idx] = E2[idx]
    E2[idx] = E1_org[idx]

    # ensure alpha < pi/2
    idx = np.sum(E1*E2, axis=1) < 0
    E2[idx] = -E2[idx]
    
    # get maximal amplitude of envelope
    normE1 = np.linalg.norm(E1, axis=1)
    normE2 = np.linalg.norm(E2, axis=1)
    cosalpha = np.sum(E1*E2, axis=1)/(normE1*normE2)
    
    TImax = 2*np.linalg.norm(np.cross(E2,E1-E2), axis=1) \
            /np.linalg.norm(E1-E2, axis=1)
    idx = normE2<=normE1*cosalpha
    TImax[idx] = 2*normE2[idx]
    return TImax
    

def get_dirTI(E1,E2,dirvec_org):
    """
    calculates the TI envelope amplitude along the direction 
    specified by vector n
        
    TIamp = | |(E1+E2)*n| - |(E1-E1)*n| |

    Parameters
    ----------
    E1 : np.ndarray
        field of electrode pair 1 (N x 3) where N is the number of 
        positions at which the field was calculated
    E2 : np.ndarray
        field of electrode pair 2 (N x 3) 
    n : np.ndarray or list
        can be either a single vector (1 x 3) that is applied to all positions
        or one vector per position (N x 3) 

    Returns
    -------
    TIamp : np.ndarray
        modulation amplitude along the direction specified by n
    """
    assert E1.shape == E2.shape
    assert E1.shape[1] == 3
    
    dirvec = np.array(dirvec_org)
    if dirvec.ndim == 1:
        assert len(dirvec) == 3
        dirvec = dirvec.reshape((1,3))
    if dirvec.shape[0] > 1:
        assert dirvec.shape == E1.shape
    dirvec = dirvec/np.linalg.norm(dirvec,axis=1)[:,None]

    TIamp = np.abs( np.abs(np.sum((E1+E2)*dirvec,axis=1)) 
                  - np.abs(np.sum((E1-E2)*dirvec,axis=1)) )
    return TIamp

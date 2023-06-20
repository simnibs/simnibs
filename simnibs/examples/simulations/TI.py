# -*- coding: utf-8 -*-
"""
 example script that runs two simnibs tDCS simulations
 and calculates maximal amplitude of the TI envelope from the E-fields
 
 Created on Thu Jun 23 17:41:21 2022

@author: axthi
"""

from copy import deepcopy
import os
import numpy as np

from simnibs import sim_struct, run_simnibs, mesh_io, ElementTags
from simnibs.utils import TI_utils as TI


"""
     set up and run simulations for the two electrode pairs
"""

# specify general parameters
S = sim_struct.SESSION()
S.subpath = 'm2m_ernie'  # m2m-folder of the subject
S.pathfem = 'TI'  # Directory for the simulation

# specify first electrode pair
tdcs = S.add_tdcslist()
tdcs.currents = [0.001, -0.001]  # Current flow though each channel (A)

electrode = tdcs.add_electrode()
electrode.channelnr = 1
electrode.centre = 'F5'  
electrode.shape = 'ellipse' 
electrode.dimensions = [40, 40]  # diameter in [mm]
electrode.thickness = 2  # 2 mm thickness

electrode = tdcs.add_electrode()
electrode.channelnr = 2
electrode.centre = 'P5'
electrode.shape = 'ellipse'
electrode.dimensions = [40, 40]
electrode.thickness = 2

# specify second electrode pair
tdcs = S.add_tdcslist(deepcopy(tdcs))
tdcs.electrode[0].centre = 'F6' 
tdcs.electrode[1].centre = 'P6' 

run_simnibs(S)


"""
    generate the TI field from the simulation results
"""
m1 = mesh_io.read_msh(os.path.join(S.pathfem, 'ernie_TDCS_1_scalar.msh'))
m2 = mesh_io.read_msh(os.path.join(S.pathfem, 'ernie_TDCS_2_scalar.msh'))

# remove all tetrahedra and triangles belonging to the electrodes so that
# the two meshes have same number of elements
tags_keep = np.hstack((np.arange(ElementTags.TH_START, ElementTags.SALINE_START - 1), np.arange(ElementTags.TH_SURFACE_START, ElementTags.SALINE_TH_SURFACE_START - 1)))
m1=m1.crop_mesh(tags = tags_keep)
m2=m2.crop_mesh(tags = tags_keep)

# calculate the maximal amplitude of the TI envelope
ef1=m1.field['E']
ef2=m2.field['E']
TImax = TI.get_maxTI(ef1.value, ef2.value)

# make a new mesh for visualization of the field strengths
# and the amplitude of the TI envelope
mout = deepcopy(m1)
mout.elmdata = []
mout.add_element_field(ef1.norm(), 'magnE - pair 1')
mout.add_element_field(ef2.norm(), 'magnE - pair 2')                    
mout.add_element_field(TImax,'TImax')
mesh_io.write_msh(mout,os.path.join(S.pathfem, 'TI.msh'))
v = mout.view(
    visible_tags=[1002, 1006],
    visible_fields='TImax',    
    )
v.write_opt(os.path.join(S.pathfem, 'TI.msh'))
mesh_io.open_in_gmsh(os.path.join(S.pathfem, 'TI.msh'), True)
            
            


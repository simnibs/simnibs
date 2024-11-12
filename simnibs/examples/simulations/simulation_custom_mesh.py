''' Example on running SimNIBS simulations 
    with a custom mesh

    Run with:
    simnibs_python simulation_custom_mesh.py

    NOTE: This example requires the mesh "myspheres.msh"
    Please see "How to create and use a custom mesh"
    in the SimNIBS tutorials for instructions to create
    the mesh
    
    Copyright (c) 2021 SimNIBS developers. Licensed under the GPL v3.
'''
import os
from simnibs import sim_struct, run_simnibs

S = sim_struct.SESSION()
S.fnamehead = 'myspheres.msh' # name of custom mesh
S.pathfem = 'simu_custom_mesh'
# Note: As there is no m2m_{subID} folder, postprocessing
#       options are not available.


# add a TDCS simulation
tdcs = S.add_tdcslist()
tdcs.currents = [0.001, -0.001]  # Current flow though each channel (A)

# 'myspheres.msh' contains a custom tissue with label number 17.
# We need to assign a conductivity to this tissue label.
# Note: Python indexing starts with 0, thus the conductivity has
#       to be assigned to index 16 of the conductivity list
tdcs.cond[16].value = 2 # [S/m]
tdcs.cond[16].name = 'custom_tissue'

electrode1 = tdcs.add_electrode()
electrode1.channelnr = 1  # Connect the electrode to the first channel
electrode1.centre = [10, 50, 50]  # position determined from the nifti file
electrode1.shape = 'rect'  # Rectangular shape
electrode1.dimensions = [50, 50]  # 50x50 mm
electrode1.thickness = 4  # 4 mm thickness

electrode2 = tdcs.add_electrode()
electrode2.channelnr = 2
electrode2.centre = [90, 50, 50]  
electrode2.shape = 'ellipse'  # Circiular shape
electrode2.dimensions = [50, 50]  # 50 mm diameter
electrode2.thickness = 4  # 4 mm thickness


# add a TMS simulation
tms = S.add_tmslist()
tms.fnamecoil = os.path.join('legacy_and_other','Magstim_70mm_Fig8.ccd')  # Choose a coil model

tms.cond[16].value = 2 # [S/m]
tms.cond[16].name = 'custom_tissue'

# Define the coil position
pos = tms.add_position()
pos.centre = [50, 50, 90]  
pos.pos_ydir = [50, 40, 90] # Polongation of coil handle (see documentation)
pos.distance = 4  #  4 mm distance from coil surface to head surface


# Run simulation
run_simnibs(S)


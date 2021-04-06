''' Example on how to run a SimNIBS TMS simulation in Python
    Run with:

    simnibs_python TMS.py

    Copyright (C) 2018 Guilherme B Saturnino
'''
from simnibs import sim_struct, run_simnibs

### General Information
S = sim_struct.SESSION()
S.subpath = 'm2m_ernie'  # m2m-folder of the subject
S.pathfem = 'tms_simu'  # Directory for the simulation

## Define the TMS simulation
tms = S.add_tmslist()
tms.fnamecoil = 'Magstim_70mm_Fig8.nii.gz'  # Choose a coil model

# Define the coil position
pos = tms.add_position()
pos.centre = 'C3'  # Place the coil over C3
pos.pos_ydir = 'CP3'  # Polongation of coil handle (see documentation)
pos.distance = 4  #  4 mm distance from coil surface to head surface

# Run Simulation
run_simnibs(S)
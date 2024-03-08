''' Example of how to run a SimNIBS tDCS simulation with ring electrodes in Python
    Run with:

    simnibs_python example_tDCS_ernie_ring.py

    Copyright (C) 2018 Guilherme B Saturnino
'''
from simnibs import sim_struct, run_simnibs

### General Infoarmation
S = sim_struct.SESSION()
S.subpath = 'm2m_ernie'  # subject folder
S.pathfem = 'tdcs_ring'  # Directory for the simulation


### Define tDCS simulation
tdcs = S.add_tdcslist()
tdcs.currents = [0.001, -0.001]  # Current flow though each channel (mA)

# First Electrode (center, circular)
central_electrode = tdcs.add_electrode()
central_electrode.channelnr = 1  # Connect the electrode to the first channel
central_electrode.centre = 'C3'  # Place it over C3
central_electrode.shape = 'ellipse'  # Elliptical / circular shape
central_electrode.dimensions = [34, 34]  # 34 mm diameter
central_electrode.thickness = 4  # 4 mm thickness

# Second electrode (surround, ring)
ring_electrode = tdcs.add_electrode()
ring_electrode.channelnr = 2
ring_electrode.centre = 'C3'  # Place it over C3
ring_electrode.shape = 'ellipse'  # Elliptical / circular shape
ring_electrode.dimensions = [100, 100]  # Ring outer diameter
ring_electrode.thickness = 4  # 4 mm thickness

hole = ring_electrode.add_hole()
hole.centre = 'C3'
hole.shape = 'ellipse'
hole.dimensions = [75, 75]  # Ring inner diameter
hole.thickness = 4

run_simnibs(S)

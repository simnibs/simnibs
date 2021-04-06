''' Example on running a SimNIBS tDCS simulation in Python
    Run with:

    simnibs_python tDCS.py

    Copyright (C) 2018 Guilherme B Saturnino
'''
from simnibs import sim_struct, run_simnibs

S = sim_struct.SESSION()
S.subpath = 'm2m_ernie'  # m2m-folder of the subject
S.pathfem = 'tdcs_simu'  # Directory for the simulation
S.map_to_surf = True # map to subject's middle gray matter surface (optional)

tdcs = S.add_tdcslist()
tdcs.currents = [0.001, -0.001]  # Current flow though each channel (A)

mc_electrode = tdcs.add_electrode()
mc_electrode.channelnr = 1  # Connect the electrode to the first channel
mc_electrode.centre = 'C3'  # Place it over C3
mc_electrode.shape = 'rect'  # Rectangular shape
mc_electrode.dimensions = [50, 50]  # 50x50 mm
mc_electrode.thickness = 4  # 4 mm thickness

so_electrode = tdcs.add_electrode()
so_electrode.channelnr = 2
so_electrode.centre = 'AF4'
so_electrode.shape = 'rect'
so_electrode.dimensions = [50, 70]
so_electrode.thickness = 4

run_simnibs(S)

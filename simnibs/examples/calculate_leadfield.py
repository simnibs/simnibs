''' Example for running a SimNIBS tDCS leadfield in Python
    Run with:

    simnibs_python calculate_leadfield.py

    Copyright (C) 2019 Guilherme B Saturnino
'''
from simnibs import sim_struct, run_simnibs
tdcs_lf = sim_struct.TDCSLEADFIELD()
# head mesh
tdcs_lf.fnamehead = 'ernie.msh'
# output directory
tdcs_lf.pathfem = 'tdcs_leadfield'
# Interpolate to the middle gray matter surface
tdcs_lf.map_to_surf = True

run_simnibs(tdcs_lf)

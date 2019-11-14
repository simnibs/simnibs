''' Example for running a SimNIBS tDCS leadfield in Python
    Run with:

    simnibs_python leadfield.py

    Copyright (C) 2019 Guilherme B Saturnino
'''
from simnibs import sim_struct, run_simnibs
tdcs_lf = sim_struct.TDCSLEADFIELD()
# head mesh
tdcs_lf.fnamehead = 'ernie.msh'
# output directory
tdcs_lf.pathfem = 'leadfield'

run_simnibs(tdcs_lf)

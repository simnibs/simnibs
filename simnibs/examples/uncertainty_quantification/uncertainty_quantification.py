''' Example for running a SimNIBS tDCS simulation in Python
    Run with:

    simnibs_python uncertainty_quantification.py

    Copyright (C) 2019 Guilherme B Saturnino
'''
from simnibs import sim_struct, run_simnibs

S = sim_struct.SESSION()
S.subpath = 'm2m_ernie'
S.pathfem = 'tdcs_uq'

tdcs = S.add_tdcslist()
tdcs.currents = [0.001, -0.001]

# Set-up the electrodes
mc_electrode = tdcs.add_electrode()
mc_electrode.channelnr = 1
mc_electrode.centre = 'C3'
mc_electrode.shape = 'rect'
mc_electrode.dimensions = [50, 50]
mc_electrode.thickness = 4

so_electrode = tdcs.add_electrode()
so_electrode.channelnr = 2
so_electrode.centre = 'AF4'
so_electrode.shape = 'rect'
so_electrode.dimensions = [50, 70]
so_electrode.thickness = 4

# Set-up the uncertain conductivities
# White Matter
tdcs.cond[0].distribution_type = 'beta'
tdcs.cond[0].distribution_parameters = [3, 3, .1, .4]
# Gray Matter
tdcs.cond[1].distribution_type = 'beta'
tdcs.cond[1].distribution_parameters = [3, 3, .1, .6]
# Compact Bone
tdcs.cond[6].distribution_type = 'beta'
tdcs.cond[6].distribution_parameters = [3, 3, 0.001, 0.012]

run_simnibs(S)

''' Example of a SimNIBS tDCS leadfield in Python
    Run with:

    simnibs_python leadfield.py

    Copyright (c) 2019 SimNIBS developers. Licensed under the GPL v3.
    
    place script in the main folder of the example dataset
'''
from simnibs import sim_struct, run_simnibs

tdcs_lf = sim_struct.TDCSLEADFIELD()
# subject folder
tdcs_lf.subpath = 'm2m_ernie'
# output directory
tdcs_lf.pathfem = 'leadfield'


# Uncoment to use the pardiso solver
#tdcs_lf.solver_options = 'pardiso'
# This solver is faster than the default. However, it requires much more
# memory (~12 GB)


run_simnibs(tdcs_lf)

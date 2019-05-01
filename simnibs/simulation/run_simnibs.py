import sys
import logging
import warnings
import argparse
import numpy as np
from ..utils.matlab_read import read_mat
from .. import __version__

def run_simnibs(simnibs_struct, cpus=1):
    """Runs a simnnibs problem.

    Parameters:
    --------------
    simnibs_struct: sim_struct.mat.SESSION of str
        SESSION of name of '.mat' file defining the simulation
    cpus: int
        Number of processes to run in parallel (if avaliable)
    """
    np.set_printoptions(precision=4)

    if isinstance(simnibs_struct, str):
        p = read_mat(simnibs_struct)
    else:
        p = simnibs_struct

    out = p.run(cpus=cpus)
    logging.shutdown()
    return out


def run_simulation(simnibs_struct, cpus=1):
    warnings.warn(
        'run_simulation deprecated, please use run_simnibs instead',
        DeprecationWarning)
    run_simnibs(simnibs_struct, cpus=cpus)

"""
TMS coil position/orientation optimizations functions.


Authors: Ole Numssen, Konstantin Weise
"""
import logging

import pyfempp

from simnibs.msh import mesh_io
from simnibs import sim_struct, run_simnibs
import numpy as np
import pickle
from simnibs.simulation.opt import get_opt_grid

n_cpu = 10
handle_direction_ref = [22.5, 45, 67.5]
radius = 20
resolution_pos = 3
resolution_angle = 10
angle_limits = [-30, 30]
target = np.mean(np.array([[-27.17, -17.94, 69.94],
                           [-28.15, -18.02, 69.24],
                           [-27.72, -18.95, 69.82]]), axis=0).tolist()


def read_msh_from_pckl(fn, m=None):
    ''' Reads a gmsh '.msh' file

    Parameters
    ------------
    fn: str
        File name
    m: simnibs.msh.Msh (optional)
        Mesh structure to be overwritten. If unset, will create a new structure

    Returns
    --------
    msh: simnibs.msh.Msh
        Mesh structure
    '''
    print("This is the monkey-patched version.")
    assert fn.endswith('.msh')
    fn = fn[:-3] + "pckl"
    return pickle.load(open(fn, 'rb'))


# replace mesh_io.read_msh with own version to load pickled mesh
mesh_io.read_msh = read_msh_from_pckl

# General Information
S = sim_struct.SESSION()
mesh = '/data/pt_01756/tmp/optim/15484.08_fixed.msh'
S.fnamehead = mesh  # head mesh
S.pathfem = '/data/pt_01756/tmp/optim'  # Directory for the simulation

# Define the TMS simulation
tms = S.add_tmslist()
tms.fnamecoil = '/data/pt_01756/coils/ccd/MagVenture_MC_B60_1233.nii.gz'  # Choose a coil from the ccd-files folder

print("prepare")
# mesh = mesh_io.read_msh(mesh)
# mesh.fix_surface_labels()
mesh = pickle.load(open("/data/pt_01756/tmp/optim/15484.08_fixed.pckl", 'rb'))

print("determine coil positions and orientations for optimization")
tms = get_opt_grid(tms=tms,
                   msh=mesh,
                   target=target,
                   handle_direction_ref=handle_direction_ref,
                   radius=radius,
                   resolution_pos=resolution_pos)
tms.pos
for i, pos in enumerate(tms.pos):  # pos = tms.pos[0]
    print(i)
    try:
        pos.pos_ydir = pos.pos_ydir.tolist()
    except Exception:
        pass
    try:
        pos.centre = pos.centre.tolist()

    except Exception:
        pass

    pos.matsimnibs = None
    pos.distance = .1
    pos.matsimnibs = pos.calc_matsimnibs(mesh, log=False)

pyfempp.create_stimsite_from_poslist("/data/pt_01756/tmp/optim/coil.hdf5", tms)
# Run Simulation
out = S.run(allow_multiple_runs=True, cpus=10)
logging.shutdown()

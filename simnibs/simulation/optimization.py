"""
TMS coil position/orientation optimizations functions.


Authors: Ole Numssen, Konstantin Weise
"""
import glob
import logging
import pandas as pd
import pyfempp
from importlib import reload  # Python 3.4+ only.
from simnibs import sim_struct, run_simnibs, msh
from simnibs.msh import mesh_io
import numpy as np
import pickle
from simnibs.simulation.opt import get_opt_grid

n_cpu = 60
handle_direction_ref = [22.5, 45, 67.5]
radius = 20
resolution_pos = 1.5
resolution_angle = 15
angle_limits = [-60, 60]
distance = 1.
target = np.mean(np.array([[-27.17, -17.94, 69.94],
                           [-28.15, -18.02, 69.24],
                           [-27.72, -18.95, 69.82]]), axis=0).tolist()

# replace mesh_io.read_msh with own version to load pickled mesh
from simnibs.simulation.opt import read_msh_from_pckl
mesh_io.read_msh = read_msh_from_pckl

# General Information
S = sim_struct.SESSION()
mesh = '/data/pt_01756/tmp/optim/15484.08_fixed.msh'
S.fnamehead = mesh  # head mesh
sim_folder = '/data/pt_01756/tmp/optim/simulations'
S.pathfem = sim_folder  # Directory for the simulation

# Define the TMS simulation
tms = S.add_tmslist()
tms.fnamecoil = '/data/pt_01756/coils/ccd/MagVenture_MC_B60_1248.nii.gz'  # Choose a coil from the ccd-files folder

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
                   angle_limits=angle_limits,
                   resolution_pos=resolution_pos)
for i, pos in enumerate(tms.pos):  # pos = tms.pos[0]
    print("Creating matsimnibs {0:0>4} ({1})".format(i, len(tms.pos)))

    pos.pos_ydir = pos.pos_ydir
    pos.centre = pos.centre
    pos.matsimnibs = None
    pos.distance = distance
    pos.matsimnibs = pos.calc_matsimnibs(mesh, log=False)

# pyfempp.create_stimsite_from_poslist("/data/pt_01756/tmp/optim/coil.hdf5", tms)
# # Run Simulation
out = S.run(allow_multiple_runs=True, cpus=60)
logging.shutdown()

all_simulations = glob.glob(sim_folder + "/*.msh")
all_simulations.sort()

# reverse monkey patch
mesh_io = reload(msh.mesh_io)
d = dict()
for i, simulation in enumerate(all_simulations):
    print("Load mesh {} {}".format(simulation, i))
    sim_msh = mesh_io.read_msh(simulation)
    elm, idx = sim_msh.find_closest_element(target, return_index=True)
    # target_e.append(sim_msh.field['normE'][idx-1])
    d[i] = [sim_msh.field['normE'][idx]]+ tms.pos[i].centre + tms.pos[i].pos_ydir

data = pd.DataFrame.from_dict(d)
data = data.transpose()
data.columns = ['normE', 'x', 'y', 'z', 'handle_x', 'handle_y', 'handle_z']
data.to_csv(sim_folder + '/results.csv')

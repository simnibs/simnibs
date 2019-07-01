"""
This is the run script for the validation part of the "How to localize..." paper.
"""

import os
import pickle
from simnibs.simulation.optim_tms import get_opt_grid, read_msh_from_pckl, optimize_tms_coil_pos
from simnibs.msh import mesh_io
import numpy as np

from simnibs.simulation import sim_struct, TMSLIST, TMSOPTIMIZATION

subject = "29965.48"
subject = "14102.d1"
subject = "15484.08"

n_cpu = 50
handle_direction_ref = [2.8, 7, 8.1]
radius = 20
resolution_pos = 1.5
resolution_angle = 15
angle_limits = [-60, 60]
distance = 1.
coil_fn = '/data/pt_01756/coils/ccd/MagVenture_MC_B60_1248.nii.gz'

targets = {'15484.08': np.mean(np.array([[-27.17, -17.94, 69.94],
                                         [-28.15, -18.02, 69.24],
                                         [-27.72, -18.95, 69.82]]), axis=0).tolist(),
           '14102.d1': np.mean(np.array([[-30.34, -0.02, 49.86],
                                         [-30.78, .084, 49.51],
                                         [-31.45, -0.01, 49.81]]), axis=0).tolist(),
           '29965.48': np.mean(np.array([[-30.34, -0.02, 49.86],
                                         [-30.78, .084, 49.51],
                                         [-31.45, -0.01, 49.81]]), axis=0).tolist()
           }

tensor_fn = "/data/pt_01756/probands/{0}/mesh/0/d2c_{0}/CTI_vn_tensor.nii.gz"

# replace mesh_io.read_msh with own version to load pickled mesh
mesh_io.read_msh = read_msh_from_pckl
# mesh_io.Msh.calc_matsimnibs = calc_matsimnibs

# setup session with global options
session = sim_struct.SESSION()
session.fnamehead = '/data/pt_01756/tmp/optim/{0}/{0}_fixed.msh'.format(subject)
session.pathfem = '/data/pt_01756/tmp/optim/{}/'.format(subject)
session.fname_tensor = tensor_fn.format(subject)

poslist = TMSLIST()
poslist.fnamecoil = coil_fn  # Choose a coil from the ccd-files folder
poslist.write_hdf5 = True
poslist.create_visuals = False
poslist.remove_msh = True
poslist.anisotropy_type = "vn"  # for isotropic: 'scalar'

session.add_poslist(poslist)

# session.read_mat_struct("/data/pt_01756/tmp/optim/15484.08/optim_session.mat")
# session.poslists[0].pos = session.poslists[0].pos[3059:]

results = optimize_tms_coil_pos(session=session,
                                target=targets[subject],
                                n_cpu=n_cpu,
                                resolution_pos=resolution_pos,
                                resolution_angle=resolution_angle,
                                handle_direction_ref=handle_direction_ref)

print(results['best_conds'])

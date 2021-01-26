import os

import simnibs
from simnibs.simulation import gpc

## Define the TMS simulation
tms = simnibs.sim_struct.TMSLIST()
tms.fnamecoil = 'Magstim_70mm_Fig8.nii.gz'
tms.mesh = simnibs.read_msh('ernie.msh')

# We need to set the EEG cap as we are giving the coil
# Coordinates as EEG positions
sub_files = simnibs.SubjectFiles('ernie.msh')
tms.eeg_cap = sub_files.eeg_cap_1010

# Define the coil position
pos = tms.add_position()
pos.centre = 'C3'
pos.pos_ydir = 'CP3'
pos.distance = 4

# Define the uncertain conductivities
tms.cond[0].distribution_type = 'beta'
tms.cond[0].distribution_parameters = [3, 3, .1, .4]
tms.cond[1].distribution_type = 'beta'
tms.cond[1].distribution_parameters = [3, 3, .1, .6]
tms.cond[3].distribution_type = 'beta'
tms.cond[3].distribution_parameters = [3, 3, 0.003, 0.012]

# Run the UQ calling with White and Gray matter as an ROI and tolerance of 1e-2
gpc.run_tms_gpc(tms, 'tms_gpc/TMS_UQ', tissues=[1, 2], eps=1e-2)

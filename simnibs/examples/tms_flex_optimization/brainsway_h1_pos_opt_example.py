import os
import numpy as np
from simnibs import opt_struct, mni2subject_coilpos

# Define coil
coil_path = os.path.join('flexible_coils', 'Brainsway_H1.tcd')

# Define coil-scalp distance
coil_skin_distance= 0

#Define initial position from matsimnibs
matsimnibs_MNI = np.array([[ -0.8660254 ,   0.        ,   0.5       , -44],
       [  0.        ,   1.        ,   0.        ,  40],
       [ -0.5       ,   0.        ,  -0.8660254 ,  59],
       [  0.        ,   0.        ,   0.        ,   1.        ]])

initial_matsimnibs = mni2subject_coilpos(matsimnibs_MNI, 'm2m_ernie', coil_skin_distance)

# Define deformation ranges
h1_translation_ranges = np.array([[-20, 20], [-20, 20], [-30, 30]])
h1_rotation_ranges = np.array([[-30, 30], [-10,10], [-5, 5]])
                         
                         
# Initialize structure for optimization
tms_opt = opt_struct.TmsFlexOptimization()

# Subject folder
tms_opt.subpath = 'm2m_ernie'

# Select output folder
tms_opt.path_optimization = 'tms_optimization/'

# Select the coil model
tms_opt.fnamecoil = os.path.join('flexible_coils', 'Brainsway_H1.tcd')

# Select a target for the optimization
pos1 = tms_opt.add_position()
pos1.matsimnibs = initial_matsimnibs

tms_opt.method = 'distance'
tms_opt.global_translation_ranges = h1_translation_ranges
tms_opt.global_rotation_ranges = h1_rotation_ranges

tms_opt.run()
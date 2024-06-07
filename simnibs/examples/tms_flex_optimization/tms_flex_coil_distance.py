'''
    Basic example to show the approach for a distance optimization 
    for a curved round coil
    
    The coil center will be placed as close as possible to position C1
    while avoiding skin intersections
'''

import os
from simnibs import opt_struct

# Initialize structure
tms_opt = opt_struct.TmsFlexOptimization()
# Subject folder
tms_opt.subpath = 'm2m_ernie'
# Select output folder
tms_opt.path_optimization = 'tms_optimization/'
# Select the coil model
tms_opt.fnamecoil = os.path.join('Drakaki_BrainStim_2022', 'MagVenture_MMC-140-II.ccd')
# Select a target for the optimization
pos = tms_opt.add_position()
pos.centre = 'C1'
# Pointing towards Cz
pos.pos_ydir = 'Cz'

tms_opt.method = 'distance'
# Note: translations and rotations are defined in the 
#       "coil coordinate system" having its origin in C1, 
#       y-axis pointing to Cz and z-axis pointing orthogonally away from 
#       the head surface at C1
#
# translations relative to C1 in [mm]
tms_opt.global_translation_ranges = [[0, 0],[0, 0],[-10, 30]]
# rotations of +/- degrees around all three axis
tms_opt.global_rotation_ranges = [[-20, 20],[-20, 20],[-20, 20]]

tms_opt.run()
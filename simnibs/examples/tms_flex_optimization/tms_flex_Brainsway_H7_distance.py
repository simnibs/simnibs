'''
    Distance optimization for a Brainsway H7 coil
    
    The coil center will be placed as close as possible (both in terms of 
    position and orientation) to the defined position while avoiding skin 
    intersections
'''
import os
from simnibs import opt_struct, mni2subject_coilpos

# Initialize structure for optimization
tms_opt = opt_struct.TmsFlexOptimization()
# Subject folder
tms_opt.subpath = 'm2m_ernie'
# Select output folder
tms_opt.path_optimization = 'tms_optimization_H7/'
# Select the coil model
tms_opt.fnamecoil = os.path.join('flexible_coils', 'Brainsway_H7.tcd')
# Desired distance from the coil to the head in [mm] 
# (standard: 4 mm, as rough estimate of the hair thickness)
tms_opt.distance = 0

# Select target position (here: via matsimnibs in MNI space)
matsimnibs_MNI = [[ -1.,   0.,   0.,   0.],
                  [  0.,   1.,   0.,  25.],
                  [  0.,   0.,  -1.,  87.],
                  [  0.,   0.,   0.,   1.]]
pos = tms_opt.add_position()
pos.matsimnibs = mni2subject_coilpos(matsimnibs_MNI, tms_opt.subpath, tms_opt.distance)

# Set optimization method and parameters: 'distance' minimizes distance to the skin
tms_opt.method = 'distance'
tms_opt.global_translation_ranges = [[-5, 5], [-5, 5], [-30, 30]]
# rotations of +/- degrees around all three axis
tms_opt.global_rotation_ranges = [[-30, 30], [-10,10], [-5, 5]]

opt_pos=tms_opt.run()

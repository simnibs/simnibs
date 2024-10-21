'''
    Optimization of the electric field strength
    in region of interest for the MagVenture MST coil
    
    The coil center position and arrangement of the two coil halves will be 
    placed to maximize the field strength in the ROI while avoiding 
    skin and self intersections
'''
import os
from simnibs import opt_struct, mni2subject_coilpos, ElementTags

# Initialize structure
tms_opt = opt_struct.TmsFlexOptimization()
# Subject folder
tms_opt.subpath = 'm2m_ernie'
# Select output folder
tms_opt.path_optimization = 'tms_optimization_MSTemag/'
# Select the coil model
tms_opt.fnamecoil = os.path.join('flexible_coils', 'MagVenture_MST-Twin.tcd')
# Desired distance from the coil to the head in [mm] 
# (standard: 4 mm, as rough estimate of the hair thickness)
tms_opt.distance = 0

# Select an initial coil position
matsimnibs_MNI = [[ 0.995 , -0.    ,  0.1005,  0. ],
                  [ 0.    , -0.7071, -0.7035, 61.7],
                  [ 0.0995,  0.7071, -0.7035, 67.4],
                  [ 0.    ,  0.    ,  0.    ,  1. ]]
pos = tms_opt.add_position()
pos.matsimnibs = mni2subject_coilpos(matsimnibs_MNI, tms_opt.subpath, tms_opt.distance)

# Select ROI in which electric field will be evaluated
roi = tms_opt.add_region_of_interest()
# define GM tissue ROI from MNI coordinate
roi.method = 'volume'
# domains to include in mask (here: GM volume)
roi.tissues = [ElementTags.GM]
# center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_space = ["mni"]
roi.roi_sphere_center = [0., 40., 64.8]
# radius of spherical ROI (in mm)
roi.roi_sphere_radius = [50]
# for visual control of the ROI before running the optimization:
# roi.subpath = tms_opt.subpath
# roi.write_visualization('','ROI_MSToptimization')

# Set optimization method and parameters: 'emag' maximizes electric field strength in ROI
tms_opt.method = 'emag'
# Note: translations and rotations are defined in the "coil coordinate system":
#       origin in the initial coil position,
#       z-axis pointing orthogonally into the head surface,
#       y-axis defined by pos.pos_ydir (set arbitrarily when using auto init)
#
# translations relative to initial position in [mm]
tms_opt.global_translation_ranges = [[-30, 30], [-30, 30], [-30, 30]]
tms_opt.global_rotation_ranges = [[-30, 30], [-30, 30], [-30, 30]]

opt_pos=tms_opt.run()

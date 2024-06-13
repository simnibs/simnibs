# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:37:30 2024

@author: bianru
"""

import os
import numpy as np
from simnibs import opt_struct, mni2subject_coilpos, RegionOfInterest, ElementTags

#Define coil
coil_path = os.path.join('flexible_coils', 'MagVenture_MST-Twin.tcd')

# Define coil-scalp distance
coil_skin_distance= 0

#Define initial position from matsimnibs
matsimnibs_MNI = np.array([[ 1. , -0. ,  0.1,  0],
       [ 0. , -0.7, -0.7, 61.7],
       [ 0.1,  0.7, -0.7, 67.4],
       [ 0. ,  0. ,  0. ,  1. ]])

initial_matsimnibs = mni2subject_coilpos(matsimnibs_MNI, "m2m_ernie", coil_skin_distance)

# Define ROI from coordinate
roi = RegionOfInterest()
roi.subpath = "m2m_ernie"
roi.method = "volume"

# domains to include in mask
roi.tissues = [ElementTags.GM]

# center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_space = ["mni"]
roi.roi_sphere_center = [0,0, 64.8]

# radius of spherical ROI (in mm)
roi.roi_sphere_radius = [50]

#uncomment to write / inspect ROI
#roi.write_visualization('C:\\Users\\bianru\\Downloads\\simnibs4_examples\\', 'roi_surf_vis')


# Define deformation ranges
mst_translation_ranges = np.array([[-30, 30], [-30, 30], [-30, 30]])
mst_rotation_ranges = np.array([[-30, 30], [-30, 30], [-95, 95]])

                         
# Initialize structure for optimization
tms_opt = opt_struct.TmsFlexOptimization()

# Subject folder
tms_opt.subpath = 'm2m_ernie'

# Select output folder
tms_opt.path_optimization = 'tms_optimization/'

# Select the coil model
tms_opt.fnamecoil = coil_path

# Select a target for the optimization
pos1 = tms_opt.add_position()
pos1.matsimnibs = initial_matsimnibs

tms_opt.method = 'emag'
tms_opt.global_translation_ranges = mst_translation_ranges
tms_opt.global_rotation_ranges = mst_rotation_ranges
tms_opt.add_region_of_interest(roi)

tms_opt.run()
import os
from simnibs import opt_struct
import simnibs

# Initialize structure
tms_opt = opt_struct.TmsFlexOptimization()
# Subject folder
tms_opt.subpath = 'm2m_ernie'
# Select output folder
tms_opt.path_optimization = 'tms_optimization/'
# Select the coil model
tms_opt.fnamecoil = os.path.join('Drakaki_BrainStim_2022', 'Deymed_70BF.ccd')
# Select a target for the optimization
pos = tms_opt.add_position()
pos.centre = 'C1'
# Pointing towards Cz
pos.pos_ydir = 'Cz'

tms_opt.method = 'emag'
tms_opt.global_translation_ranges = [[-10, 10],[-10, 10],[-10, 10]]
tms_opt.global_rotation_ranges = [[-10, 10],[-10, 10],[-10, 10]]

roi = tms_opt.add_region_of_interest()

roi.type = "GMmidlayer"

# center of spherical ROI in subject space (in mm)
roi.roi_sphere_center_subject = [-36.34, 14.53,  94.708]

# radius of spherical ROI (in mm)
roi.roi_sphere_radius = 100

tms_opt._prepare()
print(tms_opt._roi.n_center)
exit()

opt_pos=tms_opt.run()
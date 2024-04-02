"""
Example to run TESoptimize with a standard TES montage to optimize the field focality in the ROI vs non-ROI

Written by: Konstantin Weise (2023)
"""

import simnibs
from simnibs import ElementTags

# Initialize structure
opt = simnibs.opt_struct.TesFlexOptimization()

# path of m2m folder containing the headmodel
opt.subpath = 'm2m_ernie'

# output folder
opt.output_folder = f"tes_optimize_tes_focality"

# type of goal function
opt.goal = "focality"

# define threshold(s)
opt.threshold = [0.1, 0.2]

# postprocessing of e-fields ("magn": magnitude, "normal": normal component, "tangential": tangential component)
opt.e_postproc = "magn"

# define electrode
electrode = opt.add_electrode()
electrode.type = "ElectrodeArrayPair"                   # Pair of TES electrodes
electrode.center = [[0, 0]]                             # electrode center in reference electrode space (x-y plane)
electrode.length_x = [70]                               # x-dimension of electrodes
electrode.length_y = [50]                               # y-dimension of electrodes
electrode.dirichlet_correction_detailed = False         # node wise dirichlet correction
electrode.current = [0.002, -0.002]                     # electrode currents

# define ROI
roi = opt.add_roi()
roi.type = "GMmidlayer"

# center of spherical ROI in subject space (in mm)
roi.roi_sphere_center_subject = [-41, -13,  66]

# radius of spherical ROI (in mm)
roi.roi_sphere_radius = 20

# define non-ROI
roi = opt.add_roi()
roi.type = "volume"
roi.domains = [ElementTags.WM, ElementTags.GM]

# Run optimization
simnibs.run_simnibs(opt)

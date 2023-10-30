"""
Example to run TESoptimize for Temporal Interference (TI) to optimize the field focality in the ROI vs non-ROI

Written by: Konstantin Weise (2023)
"""

import simnibs
from simnibs import ElementTags

# Initialize structure
opt = simnibs.opt_struct.TESoptimize()

# path of m2m folder containing the headmodel
opt.subpath = '/data/pt_01756/probands/ernie/mesh/charm_4.0.1/m2m_ernie/'

# output folder
opt.output_folder = f"/data/pt_02381/studies/ttf/test/tes_optimize_ti_focality"

# type of goal function
opt.goal = "focality"

# define threshold(s)
opt.threshold = [0.1, 0.2]

# postprocessing function of e-fields
# "max_TI": maximal envelope of TI field magnitude
# "dir_TI_normal": maximize envelope of normal component
# "dir_TI_tangential": maximize envelope of tangential component
opt.e_postproc = "max_TI"

# define first pair of electrodes
electrode = opt.add_electrode()
electrode.type = "ElectrodeArrayPair"                   # Pair of TES electrodes
electrode.center = [[0, 0]]                             # electrode center in reference electrode space (x-y plane)
electrode.radius = [10]                                 # radius of electrodes
electrode.dirichlet_correction_detailed = False         # node wise dirichlet correction
electrode.current = [0.002, -0.002]                     # electrode currents

# define second pair of electrodes
electrode = opt.add_electrode()
electrode.type = "ElectrodeArrayPair"                   # Pair of TES electrodes
electrode.center = [[0, 0]]                             # electrode center in reference electrode space (x-y plane)
electrode.radius = [10]                                 # radius of electrodes
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
roi.type = "custom"
roi.domains = [ElementTags.WM, ElementTags.GM]

# Run optimization
simnibs.run_simnibs(opt)

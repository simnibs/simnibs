"""
Example to run TESoptimize for Temporal Interference (TI) to optimize the field intensity in the ROI

Written by: Konstantin Weise (2023)
"""

import simnibs

# Initialize structure
opt = simnibs.opt_struct.TesFlexOptimization()

# path of m2m folder containing the headmodel
opt.subpath = 'm2m_ernie'

# output folder
opt.output_folder = f"tes_optimize_ti_intensity"

# type of goal function
opt.goal = "mean"

# postprocessing function of e-fields
# "max_TI": maximize envelope of e-field magnitude
# "dir_TI_normal": maximize envelope of e-field normal component
# "dir_TI_tangential": maximize envelope of e-field tangential component
opt.e_postproc = "max_TI"

# define first pair of electrodes
electrode = opt.add_electrode_layout("ElectrodeArrayPair")
electrode.center = [[0, 0]]                             # electrode center in reference electrode space (x-y plane)
electrode.radius = [10]                                 # radius of electrodes
electrode.dirichlet_correction_detailed = False         # node wise dirichlet correction
electrode.current = [0.002, -0.002]                     # electrode currents

# define second pair of electrodes
electrode = opt.add_electrode_layout("ElectrodeArrayPair")
electrode.center = [[0, 0]]                             # electrode center in reference electrode space (x-y plane)
electrode.radius = [10]                                 # radius of electrodes
electrode.dirichlet_correction_detailed = False         # node wise dirichlet correction
electrode.current = [0.002, -0.002]                     # electrode currents

# define ROI
roi = opt.add_roi()
roi.method = "surface"
roi.surface_type = "central"

# center of spherical ROI in subject space (in mm)
roi.roi_sphere_center_space = "subject"
roi.roi_sphere_center = [-41.0, -13.0,  66.0]

# radius of spherical ROI (in mm)
roi.roi_sphere_radius = 20

# Run optimization
simnibs.run_simnibs(opt)

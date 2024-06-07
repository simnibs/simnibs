"""
Example to run TESoptimize with an HDTES montage to optimize the field intensity in the ROI

Written by: Konstantin Weise (2023)
"""

import simnibs

# Initialize structure
opt = simnibs.opt_struct.TesFlexOptimization()

# path of m2m folder containing the headmodel
opt.subpath = "m2m_ernie"

# output folder
opt.output_folder = "tes_optimze_hdtes_intensity"

# type of goal function
opt.goal = "mean"

# postprocessing of e-fields ("magn": magnitude, "normal": normal component, "tangential": tangential component)
opt.e_postproc = "magn"

# define electrode
electrode = opt.add_electrode()
electrode.type = "CircularArray"                                     # HDTES center surround montage
electrode.radius_inner = 10                                          # radius of inner electrode
electrode.radius_outer = 10                                          # radius of outer electrodes
electrode.distance_bounds = [25, 100]                                # distance bounds between inner and outer electrodes
electrode.n_outer = 4                                                # number of outer electrodes
electrode.dirichlet_correction = False                               # electrode wise dirichlet correction
electrode.dirichlet_correction_detailed = False                      # node wise dirichlet correction
electrode.current = [0.002, -0.002/4, -0.002/4, -0.002/4, -0.002/4]  # initial currents

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

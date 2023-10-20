"""
Example to run TESoptimize with an HDTES montage to optimize the field intensity in the ROI

Written by: Konstantin Weise (2023)
"""

import simnibs

# Initialize structure
opt = simnibs.opt_struct.TESoptimize()

# filename of headmodel
opt.mesh = "/data/pt_01756/probands/ernie/mesh/charm_4.0.1/m2m_ernie/ernie_mmg_no_sizing_field_3_bin.msh"

# output folder
opt.output_folder = f"/data/pt_02381/studies/ttf/test/tes_optimze_hdtes_intensity"

# type of goal function
opt.goal = "mean"

# postprocessing of e-fields
opt.e_postproc = "norm"

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
roi.center = [[-13.754, 13.104, 71.152],                             # list of ROI points
              [-16.096, 12.434, 71.120],
              [-16.141, 12.429, 70.581]]

# Run optimization
simnibs.run_simnibs(opt)

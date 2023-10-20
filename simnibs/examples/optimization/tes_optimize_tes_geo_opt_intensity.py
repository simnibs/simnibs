"""
Example to run TESoptimize with a standard TES montage including geometry optimization
to optimize the field intensity in the ROI

Written by: Konstantin Weise (2023)
"""

import simnibs

# Initialize structure
opt = simnibs.opt_struct.TESoptimize()

# filename of headmodel
opt.mesh = "/data/pt_01756/probands/ernie/mesh/charm_4.0.1/m2m_ernie/ernie_mmg_no_sizing_field_3_bin.msh"

# output folder
opt.output_folder = f"/data/pt_02381/studies/ttf/test/tes_optimize_tes_geo_opt_intensity"

# type of goal function
opt.goal = "mean"

# postprocessing function of e-fields ("norm": magnitude)
opt.e_postproc = "norm"

# define electrode
electrode = opt.add_electrode()
electrode.type = "ElectrodeArrayPair"                   # pair of TES electrodes
electrode.center = [[0, 0]]                             # electrode center in reference electrode space (x-y plane)
electrode.length_x_bounds = [50, 70]                    # x-dimension of electrodes [min, max]
electrode.length_y_bounds = [50, 70]                    # y-dimension of electrodes [min, max]
electrode.dirichlet_correction_detailed = False         # node wise dirichlet correction
electrode.current = [0.002, -0.002]                     # electrode currents

# define ROI
roi = opt.add_roi()
roi.center = [[-13.754, 13.104, 71.152],                # list of ROI points
              [-16.096, 12.434, 71.120],
              [-16.141, 12.429, 70.581]]

# Run optimization
simnibs.run_simnibs(opt)

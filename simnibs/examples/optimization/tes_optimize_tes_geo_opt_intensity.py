"""
Example to run TESoptimize with a standard TES montage including geometry optimization
to optimize the field intensity in the ROI

Written by: Konstantin Weise (2023)
"""

import simnibs

# Initialize structure
opt = simnibs.opt_struct.TESoptimize()

# path of m2m folder containing the headmodel
opt.subpath = 'm2m_ernie'

# output folder
opt.output_folder = f"tes_optimize_tes_geo_opt_intensity"

# type of goal function
opt.goal = "mean"

# postprocessing of e-fields ("magn": magnitude, "normal": normal component, "tangential": tangential component)
opt.e_postproc = "magn"

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

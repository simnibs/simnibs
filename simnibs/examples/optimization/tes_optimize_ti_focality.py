"""
Example to run TESoptimize for Temporal Interference (TI) to optimize the field focality in the ROI vs non-ROI

Written by: Konstantin Weise (2023)
"""

import simnibs

# Initialize structure
opt = simnibs.opt_struct.TESoptimize()

# filename of headmodel
opt.mesh = "/data/pt_01756/probands/ernie/mesh/charm_4.0.1/m2m_ernie/ernie_mmg_no_sizing_field_3_bin.msh"

# output folder
opt.output_folder = f"/data/pt_02381/studies/ttf/test/tes_optimize_ti_focality"

# type of goal function
opt.goal = "focality"

# define threshold(s)
opt.threshold = [0.1, 0.2]

# postprocessing function of e-fields ("max_TI": maximal envelope of TI field)
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
roi.center = [[-13.754, 13.104, 71.152],                # list of ROI points
              [-16.096, 12.434, 71.120],
              [-16.141, 12.429, 70.581]]

# define non-ROI
roi = opt.add_roi()
roi.domains = [1, 2]                                    # domain indices (1: WM, 2: GM)

# Run optimization
simnibs.run_simnibs(opt)

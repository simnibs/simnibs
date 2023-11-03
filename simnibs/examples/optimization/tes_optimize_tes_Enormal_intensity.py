"""
Example to run TESoptimize with a standard TES montage to optimize the field intensity of the normal component of the
electric field in the ROI.

Written by: Konstantin Weise (2023)
"""

import simnibs

# Initialize structure
opt = simnibs.opt_struct.TESoptimize()

# path of m2m folder containing the headmodel
opt.subpath = '/data/pt_01756/probands/ernie/mesh/charm_4.0.1/m2m_ernie/'

# output folder
opt.output_folder = f"/data/pt_02381/studies/ttf/test/tes_optimize_tes_Enormal_intensity"

# type of goal function
opt.goal = "mean"

# postprocessing of e-fields ("magn": magnitude, "normal": normal component, "tangential": tangential component)
opt.e_postproc = "normal"

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

# Run optimization
simnibs.run_simnibs(opt)

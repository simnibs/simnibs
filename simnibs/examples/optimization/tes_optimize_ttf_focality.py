"""
Example to run TESoptimize for Tumor Treating Fields (TTF) to optimize the field focality in the ROI vs non-ROI

Written by: Konstantin Weise (2023)
"""

import simnibs
from simnibs import ElementTags

# Initialize structure
opt = simnibs.opt_struct.TESoptimize()

# path of m2m folder containing the headmodel
opt.subpath = 'm2m_ernie'

# output folder
opt.output_folder = f"tes_optimize_ttf_focality"

# type of goal function
opt.goal = "focality"

# define threshold(s)
opt.threshold = [100, 100]

# postprocessing of e-fields ("magn": magnitude, "normal": normal component, "tangential": tangential component)
opt.e_postproc = "magn"

# define first pair of electrodes
opt.constrain_electrode_locations = True
electrode = opt.add_electrode()
electrode.type = "ElectrodeArrayPair"                    # Pair of TES electrodes
electrode.center = [[-33,  22],                          # electrode center in reference electrode space (x-y plane)
                    [  0,  22],
                    [ 33,  22],
                    [-33,   0],
                    [  0,   0],
                    [ 33,   0],
                    [-33, -22],
                    [  0, -22],
                    [ 33, -22]]
electrode.radius = [10, 10, 10, 10, 10, 10, 10, 10, 10]  # radius of electrodes
electrode.dirichlet_correction_detailed = False          # node wise dirichlet correction
electrode.current = [1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  # electrode currents
                    -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9]

# define second pair of electrodes
electrode = opt.add_electrode()
electrode.type = "ElectrodeArrayPair"                    # Pair of TES electrodes
electrode.center = [[-33,  22],                          # electrode center in reference electrode space (x-y plane)
                    [  0,  22],
                    [ 33,  22],
                    [-33,   0],
                    [  0,   0],
                    [ 33,   0],
                    [-33, -22],
                    [  0, -22],
                    [ 33, -22]]
electrode.radius = [10, 10, 10, 10, 10, 10, 10, 10, 10]  # radius of electrodes
electrode.dirichlet_correction_detailed = False          # node wise dirichlet correction
electrode.current = [1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  # electrode currents
                    -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9]

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

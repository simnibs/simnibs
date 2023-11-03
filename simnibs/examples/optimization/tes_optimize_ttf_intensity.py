"""
Example to run TESoptimize for Tumor Treating Fields (TTF) to optimize the field intensity in the ROI

Written by: Konstantin Weise (2023)
"""

import simnibs

# Initialize structure
opt = simnibs.opt_struct.TESoptimize()

# path of m2m folder containing the headmodel
opt.subpath = '/data/pt_01756/probands/ernie/mesh/charm_4.0.1/m2m_ernie/'

# output folder
opt.output_folder = f"/data/pt_02381/studies/ttf/test/tes_optimize_ttf_intensity"

# type of goal function
opt.goal = "mean"

# postprocessing of e-fields ("magn": magnitude, "normal": normal component, "tangential": tangential component)
opt.e_postproc = "magn"

# define first pair of electrodes
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

# Run optimization
simnibs.run_simnibs(opt)

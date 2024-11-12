"""
Example to run TESoptimize with an 4x1 center-surround TES montage to optimize 
the intensity-focality tradeoff between the field strengths in ROI vs non-ROI

Copyright (c) 2024 SimNIBS developers. Licensed under the GPL v3.
"""
from simnibs import opt_struct


''' Initialize structure '''
opt = opt_struct.TesFlexOptimization()
opt.subpath = "m2m_ernie"                                            # path of m2m folder containing the headmodel
opt.output_folder = "tes_optimze_4x1tes_focality"

''' Set up goal function '''
opt.goal = "focality"                                                # optimize intensity-focality tradeoff of "magn" ("magn" defined by e_postproc)
opt.threshold = [0.1, 0.2]                                           # define threshold(s) of the electric field in V/m in the non-ROI and the ROI:
                                                                     # if one threshold is defined, it is the goal that the e-field in the non-ROI is lower than this value and higher than this value in the ROI
                                                                     # if two thresholds are defined, the first one is the threshold of the non-ROI and the second one is for the ROI
opt.e_postproc = "magn"                                              # postprocessing of e-fields ("magn": magnitude, 
                                                                     # "normal": normal component, "tangential": tangential component)
''' Define electrodes and array layout '''
electrode_layout = opt.add_electrode_layout("CircularArray")                # Nx1 center surround montage
electrode_layout.radius_inner = 10                                          # radius of inner electrode
electrode_layout.radius_outer = 10                                          # radius of outer electrodes
electrode_layout.distance_bounds = [25, 100]                                # distance bounds between inner and outer electrodes
electrode_layout.n_outer = 4                                                # number of outer electrodes
electrode_layout.dirichlet_correction = False                               # set to True when all outer electrodes are connected to the same channel (slower)
electrode_layout.current = [0.002, -0.002/4, -0.002/4, -0.002/4, -0.002/4]  # currents for each electrode

''' Define ROI '''
roi = opt.add_roi()
roi.method = "surface"
roi.surface_type = "central"                                         # define ROI on central GM surfaces
roi.roi_sphere_center_space = "subject"
roi.roi_sphere_center = [-41.0, -13.0,  66.0]                        # center of spherical ROI in subject space (in mm)
roi.roi_sphere_radius = 20                                           # radius of spherical ROI (in mm)
# uncomment for visual control of ROI:
#roi.subpath = opt.subpath
#roi.write_visualization('','roi.msh')

''' Define non-ROI '''
# all of GM surface except a spherical region with 25 mm around roi center
non_roi = opt.add_roi()
non_roi.method = "surface"
non_roi.surface_type = "central"
non_roi.roi_sphere_center_space = "subject"
non_roi.roi_sphere_center = [-41.0, -13.0,  66.0]
non_roi.roi_sphere_radius = 25
non_roi.roi_sphere_operator = ["difference"]                         # take difference between GM surface and the sphere region
# uncomment for visual control of non-ROI:
#non_roi.subpath = opt.subpath
#non_roi.write_visualization('','non-roi.msh')

''' Run optimization '''
opt.run()
"""
Example to run TESoptimize for Temporal Interference (TI) to optimize the 
focality in the ROI vs non-ROI

© SimNIBS developers 2024 under the GPL v3 license
"""
from simnibs import opt_struct


''' Initialize structure '''
opt = opt_struct.TesFlexOptimization()
opt.subpath = 'm2m_ernie'                               # path of m2m folder containing the headmodel
opt.output_folder = "tes_optimize_ti_focality"

''' Set up goal function '''
opt.goal = "focality"                                   # optimize the focality of "max_TI" in the ROI ("max_TI" defined by e_postproc)
opt.threshold = [0.1, 0.2]                              # define threshold(s) of the electric field in V/m in the non-ROI and the ROI:
                                                        # if one threshold is defined, it is the goal that the e-field in the non-ROI is lower than this value and higher than this value in the ROI
                                                        # if two thresholds are defined, the first one is the threshold of the non-ROI and the second one is for the ROI
opt.e_postproc = "max_TI"                               # postprocessing of e-fields
                                                        # "max_TI": maximal envelope of TI field magnitude
                                                        # "dir_TI_normal": envelope of normal component
                                                        # "dir_TI_tangential": envelope of tangential component
''' Define first electrode pair '''
electrode_layout = opt.add_electrode_layout("ElectrodeArrayPair")   # Pair of TES electrode arrays (here: 1 electrode per array)
electrode_layout.radius = [10]                                      # radii of electrodes
electrode_layout.current = [0.002, -0.002]                          # electrode currents

''' Define second electrode pair '''
electrode_layout = opt.add_electrode_layout("ElectrodeArrayPair")
electrode_layout.radius = [10]
electrode_layout.current = [0.002, -0.002]

''' Define ROI '''
roi = opt.add_roi()
roi.method = "surface"
roi.surface_type = "central"                            # define ROI on central GM surfaces
roi.roi_sphere_center_space = "subject"
roi.roi_sphere_center = [-41.0, -13.0,  66.0]           # center of spherical ROI in subject space (in mm)
roi.roi_sphere_radius = 20                              # radius of spherical ROI (in mm)
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
non_roi.roi_sphere_operator = ["difference"]                             # take difference between GM surface and the sphere region
# uncomment for visual control of non-ROI:
#non_roi.subpath = opt.subpath
#non_roi.write_visualization('','non-roi.msh')

''' Run optimization '''
opt.run()

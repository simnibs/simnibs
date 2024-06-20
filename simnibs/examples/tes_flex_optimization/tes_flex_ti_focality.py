"""
Example to run TESoptimize for Temporal Interference (TI) to optimize the 
focality in the ROI vs non-ROI

© SimNIBS developers 2024 under the GPL v3 license
"""
from simnibs import opt_struct, ElementTags


''' Initialize structure '''
opt = opt_struct.TesFlexOptimization()
opt.subpath = 'm2m_ernie' # path of m2m folder containing the headmodel
opt.output_folder = "tes_optimize_ti_focality"

''' Set up goal function '''
opt.goal = "focality"                                   # optimize the focality of "max_TI" in the ROI ("max_TI" defined by e_postproc)
opt.threshold = [0.1, 0.2]                              # define threshold(s)
opt.e_postproc = "max_TI"                               # postprocessing of e-fields
                                                        # "max_TI": maximal envelope of TI field magnitude
                                                        # "dir_TI_normal": envelope of normal component
                                                        # "dir_TI_tangential": envelope of tangential component
''' Define first electrode pair '''
electrode = opt.add_electrode_layout("ElectrodeArrayPair")                   # Pair of TES electrode arrays (here: 1 electrode per array)
electrode.center = [[0, 0]]                             # electrode center in reference electrode space (x-y plane)
electrode.radius = [10]                                 # radius of electrodes
electrode.dirichlet_correction_detailed = False         # node wise dirichlet correction
electrode.current = [0.002, -0.002]                     # electrode currents

''' Define second electrode pair '''
electrode = opt.add_electrode_layout("ElectrodeArrayPair")                 # Pair of TES electrodes
electrode.center = [[0, 0]]                             # electrode center in reference electrode space (x-y plane)
electrode.radius = [10]                                 # radius of electrodes
electrode.dirichlet_correction_detailed = False         # node wise dirichlet correction
electrode.current = [0.002, -0.002]                     # electrode currents

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
roi = opt.add_roi()
roi.method = "volume"                                   # define non-ROI as volume
roi.tissues = [ElementTags.WM, ElementTags.GM]          # include WM and GM in non-ROI
# uncomment for visual control of non-ROI:
#roi.subpath = opt.subpath
#roi.write_visualization('','non-roi.msh')

''' Run optimization '''
opt.run(cpus=1)

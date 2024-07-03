"""
Example to run TESoptimize for Tumor Treating Fields (TTF) to 
maximize the field in the ROI while at the same time making
the field as unfocal as possible to cover most of the brain

© SimNIBS developers 2024 under the GPL v3 license
"""
from simnibs import opt_struct, ElementTags


''' Initialize structure '''
opt = opt_struct.TesFlexOptimization() 
opt.subpath = 'm2m_ernie'                                          # path of m2m folder containing the headmodel
opt.output_folder = "tes_optimize_ttf_focality"

''' Set up goal function '''
opt.goal = "focality_inv"                                          # optimize intensity - non-focality tradeoff of "magn" ("magn" defined by e_postproc)
opt.threshold = [100, 100]                                         # define threshold(s)
opt.e_postproc = "magn"                                            # postprocessing of e-fields ("magn": magnitude, "normal": normal component, "tangential": tangential component)
opt.constrain_electrode_locations = True                           # WHAT IS THIS?

''' Define first pair of electrode arrays '''
electrode_layout = opt.add_electrode_layout("ElectrodeArrayPair")  # Pair of TES electrode arrays
electrode_layout.center = [[-33,  22], [  0,  22], [ 33,  22],     # electrode center(s) in reference electrode space (x-y plane)
                           [-33,   0], [  0,   0], [ 33,   0],
                           [-33, -22], [  0, -22], [ 33, -22]]
electrode_layout.radius = [10, 10, 10, 10, 10, 10, 10, 10, 10]     # radii of electrodes
electrode_layout.dirichlet_correction = True                       # set to True when all electrodes of an array are connected to the same channel (slower)
electrode_layout.current = [1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,   # electrode currents: 1/9 for each electrode of the first array
                            -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9]  # -1/9 for each electrode of the second array

''' Define second pair of electrode arrays '''
electrode_layout = opt.add_electrode_layout("ElectrodeArrayPair")
electrode_layout.center = [[-33,  22], [  0,  22], [ 33,  22],
                           [-33,   0], [  0,   0], [ 33,   0],
                           [-33, -22], [  0, -22], [ 33, -22]]
electrode_layout.radius = [10, 10, 10, 10, 10, 10, 10, 10, 10] 
electrode_layout.dirichlet_correction = True
electrode_layout.current = [1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9, 
                            -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9] 

''' Define ROI '''
roi = opt.add_roi()
roi.method = "volume"
roi.tissues = [ElementTags.WM, ElementTags.GM]
roi.roi_sphere_center_space = "subject"                         # center of spherical ROI in subject space (in mm)
roi.roi_sphere_center = [34.7, -9.4, 49.8]                      # right parietal region
roi.roi_sphere_radius = 20                                      # radius of spherical ROI (in mm)
# uncomment for visual control of ROI:
#roi.subpath = opt.subpath
#roi.write_visualization('','roi.msh')

# define non-ROI
non_roi = opt.add_roi()
non_roi.method = "volume"
non_roi.tissues = [ElementTags.WM, ElementTags.GM]
non_roi.roi_sphere_center_space = "subject"
non_roi.roi_sphere_center = [34.7, -9.4, 49.8]
non_roi.roi_sphere_radius = 20
non_roi.roi_sphere_operator = ["difference"]                # take difference between brain volume and the sphere region
# uncomment for visual control of non-ROI:
#non_roi.subpath = opt.subpath
#non_roi.write_visualization('','non-roi.msh')

''' Run optimization '''
opt.run()

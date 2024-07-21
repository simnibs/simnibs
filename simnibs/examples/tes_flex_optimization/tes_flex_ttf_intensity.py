"""
Example to run TESoptimize for Tumor Treating Fields (TTF) to optimize 
the field intensity in the ROI

© SimNIBS developers 2024 under the GPL v3 license
"""
from simnibs import opt_struct, ElementTags


''' Initialize structure '''
opt = opt_struct.TesFlexOptimization()
opt.subpath = "m2m_ernie"                                     # path of m2m folder containing the headmodel
opt.output_folder = "tes_optimize_ttf_intensity"

''' Set up goal function '''
opt.goal = "mean"                                             # maximize the mean of field magnitude in the ROI
opt.e_postproc = "magn"                                       # postprocessing of e-fields ("magn": magnitude, 
                                                              # "normal": normal component, "tangential": tangential component)
opt.constrain_electrode_locations = True                      # electrode array locations are restricted to be frontal, parietal and occipital
                                                              # to reduce possibility of overlapping configurations, which will be sorted out anyway

''' Define first pair of electrode arrays '''
electrode_layout = opt.add_electrode_layout("ElectrodeArrayPair")    # Pair of TES electrode arrays (each with 9 electrodes)
electrode_layout.center = [[-33,  22], [  0,  22], [ 33,  22],       # electrode center(s) in reference electrode space (x-y plane)
                           [-33,   0], [  0,   0], [ 33,   0],
                           [-33, -22], [  0, -22], [ 33, -22]]
electrode_layout.radius = [10, 10, 10, 10, 10, 10, 10, 10, 10]       # radii of electrodes
electrode_layout.dirichlet_correction = True                         # set to True when all electrodes of an array are connected to the same channel (slower)
electrode_layout.current = [1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  # electrode currents: 1/9 for each electrode of the first array
                            -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9] # -1/9 for each electrode of the second array

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
roi.roi_sphere_center_space = "subject"                       # center of spherical ROI in subject space (in mm)
roi.roi_sphere_center = [34.7, -9.4, 49.8]                    # right parietal region
roi.roi_sphere_radius = 20                                    # radius of spherical ROI (in mm)
# uncomment for visual control of ROI:
#roi.subpath = opt.subpath
#roi.write_visualization('','roi.msh')

''' Run optimization '''
opt.run()

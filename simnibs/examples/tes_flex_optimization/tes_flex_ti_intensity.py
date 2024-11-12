"""
Example to run TESoptimize for Temporal Interference (TI) to optimize 
the field intensity in the ROI

Copyright (c) 2024 SimNIBS developers. Licensed under the GPL v3.
"""
from simnibs import opt_struct


''' Initialize structure '''
opt = opt_struct.TesFlexOptimization()
opt.subpath = 'm2m_ernie'                           # path of m2m folder containing the headmodel
opt.output_folder = "tes_optimize_ti_intensity"

''' Set up goal function '''
opt.goal = "mean"                                   # maximize the mean of "max_TI" in the ROI
opt.e_postproc = "max_TI"                           # postprocessing of e-fields:
                                                    # "max_TI": maximal envelope of e-field magnitude
                                                    # "dir_TI_normal": envelope of e-field normal component
                                                    # "dir_TI_tangential": envelope of e-field tangential component
''' Define first electrode pair '''
electrode_layout = opt.add_electrode_layout("ElectrodeArrayPair") # Pair of TES electrode arrays (here: 1 electrode per array)
electrode_layout.radius = [10]                                    # radii of electrodes
electrode_layout.current = [0.002, -0.002]                        # electrode currents

''' Define second electrode pair '''
electrode_layout = opt.add_electrode_layout("ElectrodeArrayPair")
electrode_layout.radius = [10]
electrode_layout.current = [0.002, -0.002]

''' Define ROI '''
roi = opt.add_roi()
roi.method = "surface"
roi.surface_type = "central"
roi.roi_sphere_center_space = "subject"
roi.roi_sphere_center = [-41.0, -13.0,  66.0]       # center of spherical ROI in subject space (in mm)
roi.roi_sphere_radius = 20                          # radius of spherical ROI (in mm)
# uncomment for visual control of ROI:
#roi.subpath = opt.subpath
#roi.write_visualization('','roi.msh')

''' Run optimization '''
opt.run()

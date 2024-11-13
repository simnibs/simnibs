"""
Example to run TESoptimize with a standard TES montage to optimize the 
field intensity in the ROI using a multistart approach

Copyright (c) 2024 SimNIBS developers. Licensed under the GPL v3.
"""
import os
import shutil
import numpy as np
from simnibs import opt_struct

''' Run optimization multiple times and keep best solution '''
n_multistart = 2
optim_funvalue_list = np.zeros(n_multistart)
output_folder_list = [f"tes_optimize_tes_intensity/{i_opt:02}" for i_opt in range(n_multistart)]

for i_opt in range(n_multistart):
    ''' Initialize structure '''
    opt = opt_struct.TesFlexOptimization()
    opt.subpath = 'm2m_ernie'                           # path of m2m folder containing the headmodel
    opt.output_folder = output_folder_list[i_opt]

    ''' Set up goal function '''
    opt.goal = "mean"                                   # maximize the mean of field magnitude in the ROI
    opt.e_postproc = "magn"                             # postprocessing of e-fields ("magn": magnitude,
                                                        # "normal": normal component, "tangential": tangential component)

    ''' Define electrodes and array layout '''
    electrode_layout = opt.add_electrode_layout("ElectrodeArrayPair") # Pair of TES electrode arrays (here: 1 electrode per array)s
    electrode_layout.length_x = [70]                                  # x-dimension of electrodes
    electrode_layout.length_y = [50]                                  # y-dimension of electrodes
    electrode_layout.dirichlet_correction_detailed = False            # account for inhomogenous current distribution at electrode-skin interface (slow)
    electrode_layout.current = [0.002, -0.002]                        # electrode currents

    ''' Define ROI '''
    roi = opt.add_roi()
    roi.method = "surface"
    roi.surface_type = "central"                        # define ROI on central GM surfaces
    roi.roi_sphere_center_space = "subject"
    roi.roi_sphere_center = [-41.0, -13.0,  66.0]       # center of spherical ROI in subject space (in mm)
    roi.roi_sphere_radius = 20                          # radius of spherical ROI (in mm)
    # uncomment for visual control of ROI:
    #roi.subpath = opt.subpath
    #roi.write_visualization('','roi.msh')

    ''' Run optimization '''
    opt.run()
    optim_funvalue_list[i_opt] = opt.optim_funvalue

print(f"{optim_funvalue_list=}")

''' Identify best solution '''
best_opt_idx = np.argmin(optim_funvalue_list)
print(f"{best_opt_idx=}")

''' Keep best solution and remove the others '''
for i_opt in range(n_multistart):
    if i_opt == best_opt_idx:
        shutil.copytree(output_folder_list[i_opt], os.path.split(output_folder_list[i_opt])[0], dirs_exist_ok=True)
    shutil.rmtree(output_folder_list[i_opt])

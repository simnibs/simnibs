"""
Example to run TESoptimize with a standard TES montage including  optimization
of the electrode geometries to optimize the field intensity in the ROI

© SimNIBS developers 2024 under the GPL v3 license
"""
from simnibs import opt_struct

''' Initialize structure '''
opt = opt_struct.TesFlexOptimization()
opt.subpath = 'm2m_ernie'                               # path of m2m folder containing the headmodel
opt.output_folder = "tes_optimize_tes_geo_opt_intensity"

''' Set up goal function '''
opt.goal = "mean"                                       # maximize the mean of e-field magnitude in the ROI
opt.e_postproc = "magn"                                 # postprocessing of e-fields ("magn": magnitude, 
                                                        # "normal": normal component, "tangential": tangential component)
''' Define electrodes and array layout '''
electrode = opt.add_electrode_layout("ElectrodeArrayPair")                  # pair of TES electrodes
electrode.center = [[0, 0]]                             # electrode center in reference electrode space (x-y plane)
electrode.length_x_bounds = [50, 70]                    # x-dimension of electrodes [min, max]
electrode.length_y_bounds = [50, 70]                    # y-dimension of electrodes [min, max]
electrode.dirichlet_correction_detailed = False         # node wise dirichlet correction
electrode.current = [0.002, -0.002]                     # electrode currents

''' Define ROI '''
roi = opt.add_roi()
roi.metShod = 'custom'
roi.nodes = [[-13.754, 13.104, 71.152],                # list of ROI points
             [-16.096, 12.434, 71.120],
             [-16.141, 12.429, 70.581]]

''' Run optimization '''
opt.run(cpus=1)

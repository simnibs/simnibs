''' Selecting a particular region to be heavily avoided  
'''

import simnibs

opt = simnibs.opt_struct.TDCSoptimize()
opt.leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'
opt.name = 'optimization/avoid_region'

opt.max_total_current = 2e-3
opt.max_individual_current = 1e-3
opt.max_active_electrodes = 8

target = opt.add_target()
target.positions = [-55.4, -20.7, 73.4]
target.intensity = 0.2

avoid = opt.add_avoid()
# Center of the region
avoid.positions = [-35, -19, 85]
# Radius of the region, in mm
avoid.radius = 10

simnibs.run_simnibs(opt)

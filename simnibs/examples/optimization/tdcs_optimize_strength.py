''' Optimize controlling electric field strength

Copyright (c) 2019 SimNIBS developers. Licensed under the GPL v3.
'''

import simnibs

opt = simnibs.opt_struct.TDCSoptimize()
opt.leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'
opt.name = 'optimization/strength'

opt.max_total_current = 4e-3
opt.max_individual_current = 2e-3
opt.max_active_electrodes = 8

# Target in the left motor cortex
target_left = opt.add_target()
target_left.positions = [-30.3, 5.4, 71.6]
target_left.intensity = 0.2
target_left.directions = None
# Target in the right motor cortex
target_right = opt.add_target()
target_right.positions = [36.0, 2.5, 72.6]
target_right.intensity = 0.2
target_right.directions = None

simnibs.run_simnibs(opt)

''' Example of an optimization with two targets
'''

import simnibs

opt = simnibs.opt_struct.TDCSoptimize()
opt.leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'
opt.name = 'optimization/two_targets'

opt.max_total_current = 4e-3
opt.max_individual_current = 2e-3
opt.max_active_electrodes = 8

# Target in the left motor cortex
target_left = opt.add_target()
target_left.positions = [-34.0, -21.4, 88.5]
target_left.intensity = 0.2
# Target in the right motor cortex
target_right = opt.add_target()
target_right.positions = [32.4, -25.5, 90.4]
target_right.intensity = -0.2 # negative value reverses direction

simnibs.run_simnibs(opt)

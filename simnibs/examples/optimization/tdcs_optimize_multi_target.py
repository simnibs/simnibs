''' Example of an optimization with two targets
'''

import simnibs

opt = simnibs.opt_struct.TDCSoptimize()
opt.leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'
opt.name = 'optimization/two_targets'

opt.max_total_current = 2e-3
opt.max_individual_current = 1e-3
opt.max_active_electrodes = 8

# Target in the left motor cortex
target_left = opt.add_target()
target_left.positions = [-55.4, -20.7, 73.4]
target_left.intensity = 0.2
# Target in the right motor cortex
target_right = opt.add_target()
target_right.positions = [46.2, -35.8, 80.1]
target_right.intensity = -0.2 # negative value revert the direction


simnibs.run_simnibs(opt)

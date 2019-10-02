''' Example of differet optimization set-ups


'''

from simnibs import msh, run_simnibs
from simnibs.internal import optimization

'''
Single target
'''
# Initialize structure
opt = optimization.TDCSoptimize()
# Select the leadfield file
opt.leadfield_hdf = 'tdcs_leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'
# Select a name for the optimization
opt.name = 'tdcs_leadfield/optimization_example'
# Select a maximum total current (in A)
opt.max_total_current = 2e-3
# Select a maximum current at each electrodes (in A)
opt.max_individual_current = 1e-3
# Select a maximum number of active electrodes (optional)
opt.max_active_electrodes = 8

target = opt.add_target()
# Position of target
target.positions = [-55.4, -20.7, 73.4]
# Intensity of the electric field (in V/m)
target.intensity = 0.2
run_simnibs(opt)

'''
Two targets
'''

# in python, ALWAYS initialize a new TDCSoptimize structure
opt = optimization.TDCSoptimize(
    leadfield_hdf='tdcs_leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5',
    max_total_current=2e-3,
    max_individual_current=1e-3,
    max_active_electrodes=8,
    name='tdcs_leadfield/multiple_targets')

# Target in the left motor cortex
target_left = opt.add_target()
target_left.positions = [-55.4, -20.7, 73.4]
target_left.intensity = 0.2
# Target in the right motor cortex
target_right = opt.add_target()
target_right.positions = [46.2, -35.8, 80.1]
target_right.intensity = -0.2
run_simnibs(opt)

'''
Avoidance region
'''

opt = optimization.TDCSoptimize(
    leadfield_hdf='tdcs_leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5',
    max_total_current=2e-3,
    max_individual_current=1e-3,
    max_active_electrodes=8,
    name='tdcs_leadfield/avoidance_region')

target = opt.add_target()
# Position of target
target.positions = [-55.4, -20.7, 73.4]
# Intensity of the electric field (in V/m)
target.intensity = 0.2

avoid = opt.add_avoid()
# Center of the region
avoid.positions = [-35, -19, 85]
# Radius of the region, in mm
avoid.radius = 10
run_simnibs(opt)

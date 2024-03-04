''' Example of an optimization punishing more the field in the eyes

    Copyright (C) 2019 Guilherme B Saturnino
'''
import simnibs

opt = simnibs.opt_struct.TDCSoptimize()
opt.leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'
opt.name = 'optimization/avoid'

opt.max_total_current = 2e-3
opt.max_individual_current = 1e-3
opt.max_active_electrodes = 8

target = opt.add_target()
target.positions = simnibs.mni2subject_coords([-37, -21, 58], 'm2m_ernie')
target.intensity = 0.2

avoid = opt.add_avoid()
avoid.tissues = 1006  # 1006 corresponds to the eye surface

# Run optimization
simnibs.run_simnibs(opt)

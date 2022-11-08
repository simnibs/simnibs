import simnibs

opt = simnibs.opt_struct.TESLFoptimize()
opt.leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5'
opt.name = 'optimization/MNI_target'

opt.max_total_current = 2e-3
opt.max_individual_current = 1e-3
opt.max_active_electrodes = 8

target = opt.add_target()
# Transfrorm a set of coordinates from MNI space to subject space.
# The second argument of the mni2subject_coords function
# is the path to the "m2m_subID" folder.
target.positions = simnibs.mni2subject_coords([-37, -21, 58], 'm2m_ernie')
target.intensity = 0.2

# Run optimization
simnibs.run_simnibs(opt)

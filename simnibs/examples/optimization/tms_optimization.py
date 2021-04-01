from simnibs import opt_struct

# Initialize structure
tms_opt = opt_struct.TMSoptimize()
# Subject folder
tms_opt.subpath = 'm2m_ernie'
# Select output folder
tms_opt.pathfem = 'tms_optimization/'
# Select the coil model
tms_opt.fnamecoil = 'Magstim_70mm_Fig8.nii.gz'
# Select a target for the optimization
tms_opt.target = [-39.7, 7.5, 65.6]
# Optional: Use the MKL PARDISO solver
# Will make the simulations much faster
# but has large (approx 12GB) memory usage
tms_opt.solver_options = 'pardiso'
# Run optimization to get optimal coil position
opt_pos=tms_opt.run()

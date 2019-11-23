% Example script to run a TMS optimization
% Guilherme Saturnino, 2019

tms_opt = opt_struct('TMSoptimize');
% Select the head mesh
tms_opt.fnamehead = 'ernie.msh';
% Select output folder
tms_opt.pathfem = 'tms_optimization/';
% Select the coil model
tms_opt.fnamecoil = 'Magstim_70mm_Fig8.nii.gz';
% Select a target for the optimization
tms_opt.target = [-43.4, -20.7, 83.4];
% Optional: Use the MKL PARDISO solver
% Will make the simulations much faster
% but has large (approx 12GB) memory usage
tms_opt.solver_options = 'pardiso';
% Run optimization
run_simnibs(tms_opt);
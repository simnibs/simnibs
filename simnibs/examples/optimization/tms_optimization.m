% Example script to run a TMS optimization using grid search
%
% Copyright (c) 2019 SimNIBS developers. Licensed under the GPL v3.

tms_opt = opt_struct('TMSoptimize');
% Subject folder
tms_opt.subpath = 'm2m_ernie';
% Select output folder
tms_opt.pathfem = 'tms_optimization/';
% Select the coil model
tms_opt.fnamecoil = fullfile('legacy_and_other','Magstim_70mm_Fig8.ccd');
% Select a target for the optimization
tms_opt.target = [-39.7, 7.5, 65.6];
% Optional: Use the MKL PARDISO solver
% Will make the simulations much faster
% but has large (approx 12GB) memory usage
tms_opt.solver_options = 'pardiso';
% Run optimization
run_simnibs(tms_opt);

tms_opt = opt_struct('TMSoptimize');
% Subject folder
tms_opt.subpath = 'm2m_ernie';
% Select output folder
tms_opt.pathfem = 'tms_optimization_adm';
% Select the coil model
% The ADM method requires a '.ccd' coil model
tms_opt.fnamecoil = fullfile('legacy','Magstim_70mm_Fig8.ccd');
% Select a target for the optimization
tms_opt.target = mni2subject_coords([-37, -21, 58], 'm2m_ernie');
% Use the ADM method
tms_opt.method = 'ADM';
% Run optimization
run_simnibs(tms_opt);

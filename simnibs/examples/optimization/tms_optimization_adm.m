tms_opt = opt_struct('TMSoptimize');
% Select the head mesh
tms_opt.fnamehead = 'ernie.msh';
% Select output folder
tms_opt.pathfem = 'tms_optimization/';
% Select the coil model
% The ADM method requires a '.ccd' coil model
tms_opt.fnamecoil = 'Magstim_70mm_Fig8.ccd';
% Select a target for the optimization
tms_opt.target = [-43.4, -20.7, 83.4];
% Use the ADM method
tms_opt.method = 'ADM';
% Run optimization
run_simnibs(tms_opt);

% Refining the TMS Optimization

% Initialize structure
tms_opt = opt_struct('TMSoptimize');
% subject folder
tms_opt.subpath = 'm2m_ernie';
% Select output folder
tms_opt.pathfem = 'tms_optimization_refined/';
% Select the coil model
tms_opt.fnamecoil = 'Magstim_70mm_Fig8.nii.gz';
% Select a target for the optimization
tms_opt.target = [-39.7, 7.5, 65.6];
% Optional: Use the MKL PARDISO solver
% Will make the simulations much faster
% but has large (approx 12GB) memory usage
tms_opt.solver_options = 'pardiso';

% Set the center of the search area
tms_opt.centre = [-52.5, 8.8, 78.7];
% Change the search radius
tms_opt.search_radius = 10;
% Change the search resolution
tms_opt.spatial_resolution = 2.5;
% Set the coil direction reference
tms_opt.pos_ydir = [-52.5, 2.9, 80.0];
% Change the angles to include in the search
tms_opt.search_angle = 90;
% Change the angular resolution
tms_opt.angle_resolution = 15;

run_simnibs(tms_opt)

%
% Example to run TESoptimize with an 4x1 center-surround TES montage to 
% optimize the field intensity in the ROI
%
% © SimNIBS developers 2024 under the GPL v3 license
%

% Initialize structure
opt = opt_struct('TesFlexOptimization');
opt.subpath = 'm2m_ernie';                                                  % path of m2m folder containing the headmodel
opt.output_folder = 'tes_optimze_4x1tes_intensity';

% Set up goal function
opt.goal = 'mean';                                                          % maximize the mean of 'magn' in the ROI ('magn' defined by e_postproc)
opt.e_postproc = 'magn';                                                    % postprocessing of e-fields ('magn': magnitude, 
                                                                            % 'normal': normal component, 'tangential': tangential component)
% Define electrodes and array layout
electrode_layout = opt_struct('CircularArray');                             % Nx1 center surround montage
electrode_layout.radius_inner = 10;                                         % radius of inner electrode
electrode_layout.radius_outer = 10;                                         % radius of outer electrodes
electrode_layout.distance_bounds = [25, 100];                               % distance bounds between inner and outer electrodes
electrode_layout.n_outer = 4;                                               % number of outer electrodes
electrode_layout.dirichlet_correction = false;                              % set to True when all outer electrodes are connected to the same channel (slower)
electrode_layout.current = [0.002, -0.002/4, -0.002/4, -0.002/4, -0.002/4]; % initial currents
opt.electrode = electrode_layout;

% Define ROI
roi = opt_struct('RegionOfInterest');
roi.method = 'surface';
roi.surface_type = 'central';                                               % define ROI on central GM surfaces
roi.roi_sphere_center_space = 'subject';
roi.roi_sphere_center = [-41.0, -13.0,  66.0];                              % center of spherical ROI in subject space (in mm)
roi.roi_sphere_radius = 20;                                                 % radius of spherical ROI (in mm)
opt.roi{1} = roi;
% uncomment for visual control of the ROI before running the optimization:
% roi.subpath = opt.subpath;
% roi.fname_visu = 'roi';
% run_simnibs(roi)

% Run optimization
run_simnibs(opt)

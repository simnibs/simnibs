%
% Example to run TESoptimize with a standard TES montage to optimize the 
% field intensity in the ROI
%
% Copyright (c) 2024 SimNIBS developers. Licensed under the GPL v3.
%

% Initialize structure
opt = opt_struct('TesFlexOptimization');
opt.subpath = 'm2m_ernie';                                % path of m2m folder containing the headmodel
opt.output_folder = 'tes_optimize_tes_intensity';

% Set up goal function
opt.goal = 'mean';                                        % maximize the mean of field magnitude in the ROI
opt.e_postproc = 'magn';                                  % postprocessing of e-fields ('magn': magnitude,
                                                          % 'normal': normal component, 'tangential': tangential component)

% Define electrodes and array layout
electrode_layout = opt_struct('ElectrodeArrayPair');      % Pair of TES electrode arrays (here: 1 electrode per array)
electrode_layout.length_x = [70];                         % x-dimension of electrodes
electrode_layout.length_y = [50];                         % y-dimension of electrodes
electrode_layout.dirichlet_correction_detailed = false;   % account for inhomogenous current distribution at electrode-skin interface (slow)
electrode_layout.current = [0.002, -0.002];               % electrode currents
opt.electrode = electrode_layout;

% Define ROI
roi = opt_struct('RegionOfInterest');
roi.method = 'surface';
roi.surface_type = 'central';                             % define ROI on central GM surfaces
roi.roi_sphere_center_space = 'subject';
roi.roi_sphere_center = [-41.0, -13.0,  66.0];            % center of spherical ROI in subject space (in mm)
roi.roi_sphere_radius = 20;                               % radius of spherical ROI (in mm)
opt.roi{1} = roi;
% uncomment for visual control of ROI:
% roi.subpath = opt.subpath;
% roi.fname_visu = 'roi';
% run_simnibs(roi)

% Run optimization
run_simnibs(opt)

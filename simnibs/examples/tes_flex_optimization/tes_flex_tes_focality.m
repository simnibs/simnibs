%
% Example to run TESoptimize with a standard TES montage to optimize the field focality 
% in the ROI vs non-ROI
%
% Copyright (c) 2024 SimNIBS developers. Licensed under the GPL v3.
%

% Initialize structure
opt = opt_struct('TesFlexOptimization');
opt.subpath = 'm2m_ernie';                               % path of m2m folder containing the headmodel
opt.output_folder = 'tes_optimize_tes_focality';

% Set up goal function
opt.goal = 'focality';                                   % optimize intensity-focality tradeoff of 'magn' ('magn' defined by e_postproc)
opt.threshold = [0.1, 0.2];                              % define threshold(s) of the electric field in V/m in the non-ROI and the ROI:
                                                         % if one threshold is defined, it is the goal that the e-field in the non-ROI is lower than this value and higher than this value in the ROI
                                                         % if two thresholds are defined, the first one is the threshold of the non-ROI and the second one is for the ROI
opt.e_postproc = 'magn';                                 % postprocessing of e-fields ('magn': magnitude, 
                                                         % 'normal': normal component, 'tangential': tangential component)
% Define electrodes and array layout
electrode_layout = opt_struct('ElectrodeArrayPair');     % Pair of TES electrode arrays (here: 1 electrode per array)
electrode_layout.length_x = [70];                        % x-dimension of electrodes
electrode_layout.length_y = [50];                        % y-dimension of electrodes
electrode_layout.dirichlet_correction_detailed = false;  % account for inhomogenous current distribution at electrode-skin interface (slow, but "true" is recommended for focality optimization for large electrodes)
electrode_layout.current = [0.002, -0.002];              % electrode currents
opt.electrode = electrode_layout;

% Define ROI
roi = opt_struct('RegionOfInterest');
roi.method = 'surface';
roi.surface_type = 'central';                            % define ROI on central GM surfaces
roi.roi_sphere_center_space = 'subject';
roi.roi_sphere_center = [-41.0, -13.0,  66.0];           % center of spherical ROI in subject space (in mm)
roi.roi_sphere_radius = 20;                              % radius of spherical ROI (in mm)
opt.roi{1} = roi;
% uncomment for visual control of ROI:
% roi.subpath = opt.subpath;
% roi.fname_visu = 'roi';
% run_simnibs(roi)

% Define non-ROI
% all of GM surface except a spherical region with 25 mm around roi center
non_roi = opt_struct('RegionOfInterest');
non_roi.method = 'surface';
non_roi.surface_type = 'central';
non_roi.roi_sphere_center_space = 'subject';
non_roi.roi_sphere_center = [-41.0, -13.0,  66.0];
non_roi.roi_sphere_radius = 25;
non_roi.roi_sphere_operator = ['difference'];            % take difference between GM surface and the sphere region
opt.roi{2} = non_roi;
% uncomment for visual control of the ROI before running the optimization:
% non_roi.subpath = opt.subpath;
% non_roi.fname_visu = 'non-roi';
% run_simnibs(non_roi)

% Run optimization
run_simnibs(opt)

%
% Example to run TESoptimize for Temporal Interference (TI) to optimize the 
% focality in the ROI vs non-ROI
%
% © SimNIBS developers 2024 under the GPL v3 license
%

% Initialize structure
opt = opt_struct('TesFlexOptimization');
opt.subpath = 'm2m_ernie';                              % path of m2m folder containing the headmodel
opt.output_folder = 'tes_optimize_ti_focality';

% Set up goal function
opt.goal = 'focality';                                  % optimize the focality of 'max_TI' in the ROI ('max_TI' defined by e_postproc)
opt.threshold = [0.1, 0.2];                             % define threshold(s)
opt.e_postproc = 'max_TI';                              % postprocessing of e-fields
                                                        % 'max_TI': maximal envelope of TI field magnitude
                                                        % 'dir_TI_normal': envelope of normal component
                                                        % 'dir_TI_tangential': envelope of tangential component
% Define first electrode pair
electrode_layout = opt_struct('ElectrodeArrayPair');    % Pair of TES electrode arrays (here: 1 electrode per array)
electrode_layout.radius = [10];                         % radii of electrodes
electrode_layout.current = [0.002, -0.002];             % electrode currents
opt.electrode{1} = electrode_layout;

% Define second electrode pair
electrode_layout = opt_struct('ElectrodeArrayPair');
electrode_layout.radius = [10];
electrode_layout.current = [0.002, -0.002];
opt.electrode{2} = electrode_layout;

% Define ROI
roi = opt_struct('RegionOfInterest');
roi.method = 'surface';
roi.surface_type = 'central';                           % define ROI on central GM surfaces
roi.roi_sphere_center_space = 'subject';
roi.roi_sphere_center = [-41.0, -13.0,  66.0];          % center of spherical ROI in subject space (in mm)
roi.roi_sphere_radius = 20;                             % radius of spherical ROI (in mm)
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
non_roi.roi_sphere_operator = ['difference'];           % take difference between GM surface and the sphere region
opt.roi{2} = non_roi;
% uncomment for visual control of the ROI before running the optimization:
% non_roi.subpath = opt.subpath;
% non_roi.fname_visu = 'non-roi';
% run_simnibs(non_roi)

% Run optimization
run_simnibs(opt)

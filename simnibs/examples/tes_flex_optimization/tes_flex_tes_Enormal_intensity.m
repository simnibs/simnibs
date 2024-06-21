%
% Example to run TESoptimize with a standard TES montage to optimize the strength 
% of the normal component of the electric field in the ROI.
%
% © SimNIBS developers 2024 under the GPL v3 license
%

% Initialize structure
opt = opt_struct('TesFlexOptimization');
opt.subpath = 'm2m_ernie';                                 % path of m2m folder containing the headmodel
opt.output_folder = 'tes_optimize_tes_Enormal_intensity';

% Set up goal function
opt.goal = 'mean';                                         % maximize the mean of 'normal' in the ROI ('normal' defined by e_postproc)
opt.e_postproc = 'normal';                                 % postprocessing of e-fields ('magn': magnitude,
                                                           % 'normal': normal component, 'tangential': tangential component)
% Define electrodes and array layout
electrode_layout = opt_struct('ElectrodeArrayPair');       % Pair of TES electrode arrays (here: 1 electrode per array)
electrode_layout.length_x = [70];                          % x-dimension of electrodes
electrode_layout.length_y = [50];                          % y-dimension of electrodes
electrode_layout.dirichlet_correction_detailed = false;    % account for inhomogenous current distribution at electrode-skin interface (slow)
electrode_layout.current = [0.002, -0.002];                % electrode currents
opt.electrode = electrode_layout;

% Define ROI
roi = opt_struct('RegionOfInterest');
roi.method = 'surface';
roi.surface_type = 'central';                              % define ROI on central GM surfaces
roi.roi_sphere_center_space = 'subject';                 
roi.roi_sphere_center = [-41, -13,  66];                   % center of spherical ROI in subject space (in mm)
roi.roi_sphere_radius = 20;                                % radius of spherical ROI (in mm)
opt.roi{1} = roi;
% uncomment for visual control of ROI:
% roi.subpath = opt.subpath;
% roi.fname_visu = 'roi';
% run_simnibs(roi)

% Run optimization
run_simnibs(opt)

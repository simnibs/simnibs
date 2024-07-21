%
% Example to run TESoptimize for Tumor Treating Fields (TTF) to 
% maximize the field in the ROI while at the same time making
% the field as unfocal as possible to cover most of the brain
%
% � SimNIBS developers 2024 under the GPL v3 license
%

% Initialize structure %
opt = opt_struct('TesFlexOptimization');
opt.subpath = 'm2m_ernie';                                            % path of m2m folder containing the headmodel
opt.output_folder = 'tes_optimize_ttf_focality';

% Set up goal function %
opt.goal = 'focality_inv';                                            % optimize intensity - non-focality tradeoff of 'magn' ('magn' defined by e_postproc)
opt.threshold = [100, 100];                                           % define threshold(s) of the electric field in V/m in the non-ROI and the ROI:
                                                                      % if one threshold is defined, it is the goal that the e-field in the non-ROI is lower than this value and higher than this value in the ROI
                                                                      % if two thresholds are defined, the first one is the threshold of the non-ROI and the second one is for the ROI
opt.e_postproc = 'magn';                                              % postprocessing of e-fields ('magn': magnitude, 'normal': normal component, 'tangential': tangential component)
opt.constrain_electrode_locations = true;                             % electrode array locations are restricted to be frontal, parietal and occipital
                                                                      % to reduce possibility of overlapping configurations, which will be sorted out anyway

% Define first pair of electrode arrays %
electrode_layout = opt_struct('ElectrodeArrayPair');                  % Pair of TES electrode arrays
electrode_layout.center = {[-33,  22], [  0,  22], [ 33,  22], ...    % electrode center(s) in reference electrode space (x-y plane)
                           [-33,   0], [  0,   0], [ 33,   0], ...
                           [-33, -22], [  0, -22], [ 33, -22]};
electrode_layout.radius = [10, 10, 10, 10, 10, 10, 10, 10, 10];       % radii of electrodes
electrode_layout.dirichlet_correction = true;                         % set to true when all electrodes of an array are connected to the same channel (slower)
electrode_layout.current = [1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9, ... % electrode currents: 1/9 for each electrode of the first array
                           -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9];    % -1/9 for each electrode of the second array
opt.electrode{1} = electrode_layout;

% Define second pair of electrode arrays %
electrode_layout = opt_struct('ElectrodeArrayPair');
electrode_layout.center = {[-33,  22], [  0,  22], [ 33,  22], ...
                    [-33,   0], [  0,   0], [ 33,   0], ...
                    [-33, -22], [  0, -22], [ 33, -22]};
electrode_layout.radius = [10, 10, 10, 10, 10, 10, 10, 10, 10]; 
electrode_layout.dirichlet_correction = true;
electrode_layout.current = [1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9,  1./9, ...
                           -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9, -1./9]; 
opt.electrode{2} = electrode_layout;

% Define ROI %
roi = opt_struct('RegionOfInterest');
roi.method = 'volume';                                                % create volume mask
roi.tissues = [1, 2];                                                 % select white matter (index 1) and gray matter (index 2)
roi.roi_sphere_center_space = 'subject';                              % center of spherical ROI in subject space (in mm)
roi.roi_sphere_center = [34.7, -9.4, 49.8];                           % right parietal region
roi.roi_sphere_radius = 20;                                           % radius of spherical ROI (in mm)
opt.roi{1} = roi;
% uncomment for visual control of ROI:
% roi.subpath = opt.subpath;
% roi.fname_visu = 'roi';
% run_simnibs(roi)

% Define non-ROI
non_roi = opt_struct('RegionOfInterest');
non_roi.method = 'volume';
non_roi.tissues = [1, 2];   
non_roi.roi_sphere_center_space = 'subject';
non_roi.roi_sphere_center = [34.7, -9.4, 49.8];
non_roi.roi_sphere_radius = 20;
non_roi.roi_sphere_operator = ['difference'];                         % take difference between brain volume and the sphere region
opt.roi{2} = non_roi;
% uncomment for visual control of the ROI before running the optimization:
% non_roi.subpath = opt.subpath;
% non_roi.fname_visu = 'non-roi';
% run_simnibs(non_roi)

% Run optimization
run_simnibs(opt)

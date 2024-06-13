%
% Optimization of the electric field strength
% in region of interest for the MagVenture MST coil
% 
% The coil center position and arrangement of the two coil halves will be 
% placed to maximize the field strength in the ROI while avoiding 
% skin and self intersections
%

% Initialize structure
tms_opt = opt_struct('TmsFlexOptimization');
% Subject folder
tms_opt.subpath = 'm2m_ernie';
% Select output folder
tms_opt.path_optimization = 'tms_optimization_MSTemag/';
% Select the coil model
tms_opt.fnamecoil = fullfile('flexible_coils', 'MagVenture_MST-Twin.tcd');
% Desired distance from the coil to the head in [mm] 
% (standard: 4 mm, as rough estimate of the hair thickness)
tms_opt.distance = 0;

% Select an initial coil position
% (here: via coil center and coil axis orientations in MNI space)
center_MNI = [0., 61.7, 67.4];
ydir_MNI   = [0., -0.7071, 0.7071];
zdir_MNI   = [0.1005, -0.7035, -0.7035];
tms_opt.pos = sim_struct('POSITION');
tms_opt.pos.matsimnibs = mni2subject_coilpos(center_MNI, ydir_MNI, zdir_MNI, ...
                                             tms_opt.subpath, tms_opt.distance);

% Select ROI in which electric field will be evaluated
tms_opt.roi = opt_struct('RegionOfInterest');
% define GM tissue ROI from MNI coordinate
tms_opt.roi.method = 'volume';
% domains to include in mask (here: GM volume)
tms_opt.roi.tissues = 2;
% center of spherical ROI in MNI space (in mm)
tms_opt.roi.roi_sphere_center_space = {'mni'};
tms_opt.roi.roi_sphere_center = [0., 0., 64.8];
% radius of spherical ROI (in mm)
tms_opt.roi.roi_sphere_radius = 50;
% for visual control of the ROI before running the optimization:
% tms_opt.roi.subpath = tms_opt.subpath;
% tms_opt.roi.fname_visu = 'ROI_MSToptimization';
% run_simnibs(tms_opt.roi)

% Set optimization method and parameters: 'emag' maximizes electric field strength in ROI
tms_opt.method = 'emag';
% Note: translations and rotations are defined in the "coil coordinate system":
%       origin in the initial coil position,
%       z-axis pointing orthogonally into the head surface,
%       y-axis defined by pos.pos_ydir (set arbitrarily when using auto init)
%
% translations relative to initial position in [mm]
tms_opt.global_translation_ranges = {[-30, 30], [-30, 30], [-30, 30]};
tms_opt.global_rotation_ranges = {[-30, 30], [-30, 30], [-95, 95]};

run_simnibs(tms_opt)

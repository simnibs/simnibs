%
% Basic example demonstrating optimization of the electric field strength
% in region of interest for a flat figure-8 coil
%
% The coil center will be placed to maximize the field strength in the ROI
% while avoiding skin intersections
%

% Initialize structure
tms_opt = opt_struct('TmsFlexOptimization');
% Subject folder
tms_opt.subpath = 'm2m_ernie';
% Select output folder
tms_opt.path_optimization = 'tms_optimization/';
% Select the coil model
tms_opt.fnamecoil = fullfile('Drakaki_BrainStim_2022', 'MagVenture_MCF-B65.ccd');
% Desired distance from the coil to the head in [mm] 
% (standard: 4 mm, as rough estimate of the hair thickness)
tms_opt.distance = 0;

% Select an initial coil position (optional, otherwise the starting position is automatically determined)
%tms_opt.pos = sim_struct('POSITION');
%tms_opt.pos.centre = 'C1';
%tms_opt.pos.pos_ydir = 'Cz'; % Pointing towards Cz

% Select ROI in which electric field will be evaluated
tms_opt.roi = opt_struct('RegionOfInterest');
% Define a ROI on the central gray matter surface(s)
tms_opt.roi.method = 'surface';
tms_opt.roi.surface_type = 'central';
% Center of spherical ROI in subject space (in mm)
tms_opt.roi.roi_sphere_center_space = 'subject';
tms_opt.roi.roi_sphere_center = [-29.90, 1.29, 72.47];
% Radius of spherical ROI (in mm)
tms_opt.roi.roi_sphere_radius = 30;

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
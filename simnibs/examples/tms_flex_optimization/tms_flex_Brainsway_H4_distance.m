%
%  Distance optimization for a Brainsway H4 coil
% 
%  The coil center will be placed as close as possible (both in terms of 
%  position and orientation) to the defined position while avoiding skin 
%  intersections
%

% Initialize structure for optimization
tms_opt = opt_struct('TmsFlexOptimization');
% Subject folder
tms_opt.subpath = 'm2m_ernie';
% Select output folder
tms_opt.path_optimization = 'tms_optimization_H4/';
% Select the coil model
tms_opt.fnamecoil = fullfile('flexible_coils', 'Brainsway_H4.tcd');
% Desired distance from the coil to the head in [mm] 
% (standard: 4 mm, as rough estimate of the hair thickness)
tms_opt.distance = 0;

% Select target position
% (here: via coil center and coil axis orientations in MNI space)
center_MNI = [0., 53., 63.];
ydir_MNI   = [0., 1., 0.];
zdir_MNI   = [0., 0., -1.];
tms_opt.pos = sim_struct('POSITION');
tms_opt.pos.matsimnibs = mni2subject_coilpos(center_MNI, ydir_MNI, zdir_MNI, ...
                                             tms_opt.subpath, tms_opt.distance);

% Set optimization method and parameters: 'distance' minimizes distance to the skin
tms_opt.method = 'distance';
tms_opt.global_translation_ranges = {[-20, 20], [-20, 20], [-30, 30]};
% rotations of +/- degrees around all three axis
tms_opt.global_rotation_ranges = {[-30, 30], [-10,10], [-5, 5]};

run_simnibs(tms_opt)

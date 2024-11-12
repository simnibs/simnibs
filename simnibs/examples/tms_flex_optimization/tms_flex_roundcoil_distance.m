%
% Basic example demonstrating a distance optimization for a curved round coil
%   
% The coil center will be placed as close as possible to position C1
% while avoiding skin intersections
%
% Copyright (c) 2024 SimNIBS developers. Licensed under the GPL v3.

% Initialize structure
tms_opt = opt_struct('TmsFlexOptimization');
% Subject folder
tms_opt.subpath = 'm2m_ernie';
% Select output folder
tms_opt.path_optimization = 'tms_optimization/';
% Select the coil model
tms_opt.fnamecoil = fullfile('Drakaki_BrainStim_2022', 'MagVenture_MMC-140-II.ccd');
% Desired distance from the coil to the head in [mm] 
% (standard: 4 mm, as rough estimate of the hair thickness)
tms_opt.distance = 0;

% Select a target position
tms_opt.pos = sim_struct('POSITION');
tms_opt.pos.centre = 'C1';
% Pointing towards Cz
tms_opt.pos.pos_ydir = 'Cz';

% Set optimization method and parameters: 'distance' minimizes distance to the skin
tms_opt.method = 'distance';
% Note: translations and rotations are defined in the "coil coordinate system":
%       origin in the target coil position (here: C1),
%       z-axis pointing orthogonally into the head surface,
%       y-axis defined by pos.pos_ydir (here: pointing to Cz)
%
% translations relative to C1 in [mm]
tms_opt.global_translation_ranges = {[0, 0], [0, 0], [-30, 10]};
% rotations of +/- degrees around all three axis
tms_opt.global_rotation_ranges = {[-20, 20], [-20, 20], [-20, 20]};

run_simnibs(tms_opt)
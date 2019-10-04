%% Single Target

% Initialize structure
opt = opt_struct('TDCSoptimize');
% Select the leadfield file
opt.leadfield_hdf = 'tdcs_leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5';
% Select a name for the optimization
opt.name = 'tdcs_leadfield/optimization_example';
% Select a maximum total current
opt.max_total_current = 2e-3;
% Select a maximum current at each electrodes
opt.max_individual_current = 1e-3;
% Select a maximum number of active electrodes (optional)
opt.max_active_electrodes = 8;

% Position of target
opt.target(1).positions = [-55.4, -20.7, 73.4];
% Intensity of the electric field (in V/m)
opt.target(1).intensity = 0.2;

run_simnibs(opt)
%% Two targets

opt.name = 'tdcs_leadfield/multiple_targets';
% Target in the left motor cortex
opt.target(1).positions = [-55.4, -20.7, 73.4];
opt.target(1).intensity = 0.2;
% Target in the right motor cortex
opt.target(2).positions = [46.2, -35.8, 80.1];
opt.target(2).intensity = -0.2;
run_simnibs(opt)

%% Avoidance region

opt.name = 'tdcs_leadfield/avoidance_region';
% Position of target
opt.target(1).positions = [-55.4, -20.7, 73.4];
% Intensity of the electric field (in V/m)
opt.target(1).intensity = 0.2;

% Center of the region
opt.avoid(1).positions = [-35, -19, 85];
% Radius of the region, in mm
opt.avoid(1).radius = 10;
run_simnibs(opt)
% Selecting a particular region to be heavily avoided  
%

opt = opt_struct('TDCSoptimize');
opt.leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5';
opt.name = 'optimization/avoid_region';

opt.max_total_current = 2e-3;
opt.max_individual_current = 1e-3;
opt.max_active_electrodes = 8;

% Set target
opt.target.positions = [-55.4, -20.7, 73.4];
opt.target.intensity = 0.2;

% Center of the region
opt.avoid.positions = [-35, -19, 85];
% Radius of the region, in mm
opt.avoid.radius = 10;

run_simnibs(opt);

% Example of an optimization with two targets
%

opt = opt_struct('TDCSoptimize');
opt.leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5';
opt.name = 'optimization/two_targets';

opt.max_total_current = 4e-3;
opt.max_individual_current = 2e-3;
opt.max_active_electrodes = 8;

% Target in the left motor cortex
opt.target(1).positions = [-34.0 -21.4 88.5];
opt.target(1).intensity = 0.2;
% Target in the right motor cortex
opt.target(2).positions  = [32.4 -25.5 90.4];
opt.target(2).intensity = -0.2;


run_simnibs(opt);

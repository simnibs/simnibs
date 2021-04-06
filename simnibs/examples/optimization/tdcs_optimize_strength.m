% Optimize controlling electric field strength
%

opt = opt_struct('TDCSoptimize');
opt.leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5';
opt.name = 'optimization/strength';

opt.max_total_current = 4e-3;
opt.max_individual_current = 2e-3;
opt.max_active_electrodes = 8;

% Target in the left motor cortex
opt.target(1).positions = [-30.3, 5.4, 71.6];
opt.target(1).intensity = 0.2;
opt.target(1).directions='none';
% Target in the right motor cortex
opt.target(2).positions  = [36.0, 2.5, 72.6];
opt.target(2).intensity = 0.2;
opt.target(2).directions='none';

run_simnibs(opt);

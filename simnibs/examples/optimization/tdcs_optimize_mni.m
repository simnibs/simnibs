%  Example on how to set the target position using a MNI coordinate 
%
% Copyright (c) 2019 SimNIBS developers. Licensed under the GPL v3.

opt = opt_struct('TDCSoptimize');
opt.leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5';
opt.name = 'optimization/MNI_target';

opt.max_total_current = 2e-3;
opt.max_individual_current = 1e-3;
opt.max_active_electrodes = 8;

% Transfrorm a set of coordinates from MNI space to subject space.
% The second argument of the mni2subject_coords function
% is the path to the "m2m_subID" folder.
opt.target.positions = mni2subject_coords([-37, -21, 58], 'm2m_ernie');
opt.target.intensity = 0.2;

% Run optimization
run_simnibs(opt);

% Example of an optimization punishing more the field in the eyes
%
% Copyright (c) 2019 SimNIBS developers. Licensed under the GPL v3.

opt = opt_struct('TDCSoptimize');
opt.leadfield_hdf = 'leadfield/ernie_leadfield_EEG10-10_UI_Jurak_2007.hdf5';
opt.name = 'optimization/avoid';

opt.max_total_current = 2e-3;
opt.max_individual_current = 1e-3;
opt.max_active_electrodes = 8;

opt.target.positions = mni2subject_coords([-37, -21, 58], 'm2m_ernie');
opt.target.intensity = 0.2;

opt.avoid.tissues = 1006; % 1006 corresponds to the eye surface

% Run optimization
run_simnibs(opt)

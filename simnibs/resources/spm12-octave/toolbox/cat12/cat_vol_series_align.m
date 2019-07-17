function out = cat_vol_series_align(job)
% Longitudinal rigid registration of image series
% FORMAT out = cat_vol_series_align(job)
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
%
% modified version of
% John Ashburner
% spm_series_align.m 5044 2012-11-09 13:40:35Z john
%
% $Id cat_vol_series_align.m $

N = numel(job.data);

noise = zeros(N,1);

for i=1:N,
    % Make an estimate of the scanner noise
    noise(i,1) = spm_noise_estimate(job.data{i});
    fprintf('Estimated noise sd for "%s" = %g\n', job.data{i}, noise(i,1));
end
prec   = noise.^(-2);

bparam    = [0 0 job.bparam];

Nii    = nifti(strvcat(job.data));

% always write realigned images and average to disk
output = [{'wimg'}, {'wavg'}];

dat    = cat_vol_groupwise_ls(Nii, output, prec, bparam);
out.avg{1} = dat.avg;
out.rimg   = dat.rimg;

return


function cat_vol_slice_overlay_ui
%__________________________________________________________________________
% Christian Gaser
% $Id: cat_vol_slice_overlay_ui.m 1234 2017-12-04 10:52:18Z gaser $

OV.reference_image = fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm','Template_T1_IXI555_MNI152_GS.nii');
OV.reference_range = [0.2 1.0];                         % intensity range for reference image
OV.opacity = Inf;                                      % transparence value for overlay (<1)
OV.cmap    = jet;                                      % colormap for overlay

% name of files
OV.name = char(fullfile(spm('dir'),'tpm','TPM.nii,1'),...
               fullfile(spm('dir'),'tpm','labels_Neuromorphometrics.nii'));
                
% range for each file
% Use range 0..0 if you want to autoscale range.
% If you are using log. scaling, check the highest p-value in the table
% and approximate the range; e.g. for a max. p-value of p=1e-7 and a
% threshold of p<0.001 use a range of [3 7]. Check cat_stat_spmT2x.m for details.
% If you are unsure, simply use the autorange option by using a range of [0 0].
% The log-scaled values are calculated by -log10(1-p):
% p-value       -log10(1-P)
%  0.1           1
%  0.05          1.3
%  0.01          2
%  0.001         3
%  0.0001        4

% Number of fields in range should be the same as number of files (see above)
% or give one field, which is valid for all.
% Be carefule: intensities below the lower range are not shown!
OV.range   =[[0.5 1]; [0.5 1]];

% OV.func can be used to set the image to defined values (e.g. NaN) for the given range
%OV.func = 'i1(i1>-1.3 & i1<1.3)=NaN;';

% selection of slices and orientations
% if OV.slices_str is an empty string then slices with local maxima are estimated automatically
OV.slices_str = char('','-30:2:30','-20:5:45');
OV.transform = char('axial','sagittal','coronal');

% define output format of slices
OV.labels.format = '%3.1f';

% define number of columns and rows
% comment this out for interactive selection
OV.xy = [3 5];

% save result as png/jpg/pdf
% comment this out for interactive selection or use 'none' for not 
% saving any file
OV.save = 'result.png';

% Comment this out if you don't wish slice labels
%OV.labels = [];

% Comment this out if you don't wish colorbar
%OV.cbar = [];

cat_vol_slice_overlay(OV)

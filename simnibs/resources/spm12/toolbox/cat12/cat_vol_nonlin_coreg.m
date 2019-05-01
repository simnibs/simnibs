%-----------------------------------------------------------------------
% Job for non-linear coregistration batch
% Christian Gaser
% $Id: cat_vol_nonlin_coreg.m 1220 2017-11-20 13:34:33Z gaser $
%-----------------------------------------------------------------------

global vox reg bb
warning('off','MATLAB:DELETE:FileNotFound');

matlabbatch{1}.spm.spatial.preproc.channel.vols = '<UNDEFINED>';

matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,1')};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,2')};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,3')};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,4')};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,5')};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {fullfile(spm('dir'),'tpm','TPM.nii,6')};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{2}.spm.spatial.smooth.data(2) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{2}.spm.spatial.smooth.data(3) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
matlabbatch{2}.spm.spatial.smooth.data(4) = cfg_dep('Segment: c4 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{4}, '.','c', '()',{':'}));
matlabbatch{2}.spm.spatial.smooth.data(5) = cfg_dep('Segment: c5 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{5}, '.','c', '()',{':'}));
matlabbatch{2}.spm.spatial.smooth.data(6) = cfg_dep('Segment: c6 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{6}, '.','c', '()',{':'}));
matlabbatch{2}.spm.spatial.smooth.fwhm = [3 3 3];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';
matlabbatch{3}.spm.util.cat.vols(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{3}.spm.util.cat.name = 'TPM_tmp.nii';
matlabbatch{3}.spm.util.cat.dtype = 16;
matlabbatch{4}.spm.spatial.normalise.estwrite.subj.vol = '<UNDEFINED>';
matlabbatch{4}.spm.spatial.normalise.estwrite.subj.resample = '<UNDEFINED>';
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.tpm(1) = cfg_dep('3D to 4D File Conversion: Concatenated 4D Volume', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','mergedfile'));
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.affreg = 'subj';
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.reg = reg*[0 1e-05 0.005 0.0005 0.002];
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.fwhm = 3;
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.bb = bb;
matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.vox = vox;
matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Segment: Seg Params', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','param', '()',{':'}));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Segment: c1 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{1}, '.','c', '()',{':'}));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('Segment: c2 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{2}, '.','c', '()',{':'}));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(4) = cfg_dep('Segment: c3 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{3}, '.','c', '()',{':'}));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(5) = cfg_dep('Segment: c4 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{4}, '.','c', '()',{':'}));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(6) = cfg_dep('Segment: c5 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{5}, '.','c', '()',{':'}));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(7) = cfg_dep('Segment: c6 Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','tiss', '()',{6}, '.','c', '()',{':'}));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(8) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.files(9) = cfg_dep('3D to 4D File Conversion: Concatenated 4D Volume', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','mergedfile'));
matlabbatch{5}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;

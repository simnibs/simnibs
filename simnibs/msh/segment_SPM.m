function segment_SPM(TPM_dir, T1, T2, T2biasreg, dsf,skip_coreg)
% 
% PARAMETERS
% ----------
% T1 : T1 weighted image.
% T2 : T2 weighted image.
% TPM_dir : 
% T2biasreg : bias regularization for the (optional) T2 image.
% dsf       : downsampling factor for segmentation (0: no downsampling).
% 
% RETURNS
% ----------

% CHECK INPUTS
if ~exist('T2biasreg','var') || isempty(T2biasreg)
    biasreg = 0.01;
end
if ~exist('dsf','var') || isempty(dsf)
    dsf = 3;
end
if exist('T2','var') && ~isempty(T2)
    fprintf('Bias regularization for T2 image     : %d\n', T2biasreg)
end
fprintf('Downsampling factor for segmentation : %d\n\n', dsf)

if nargin < 2
   T1 = spm_select(1,'image','Please select T1 image.');
end
if nargin < 1
   TPM_dir = spm_select(1,'dir','Please select the directory containing the tissue probability maps to use.');
   T1 = spm_select(1,'image','Please select T1 image.');
end
%templatesdir = fullfile(TPM_dir,'templates');
templatesdir = TPM_dir;

% GET PATH TO TISSUE PROBABILITY MAPS
TPM = cell(1,6);
for i = 1:6
    TPM(i)=cellstr(fullfile(TPM_dir,sprintf('eTPM.nii,%d',i)));
end
% get voxel size of TPM
mni = spm_vol(TPM{1});
vx_size = sqrt(sum(mni.mat(1:3,1:3)^2));

% initialize SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline',true);

% START SEGMENTATION AND (REVERSE) COREGISTRATION OF THE NORMALIZED IMAGE
% TO THE ORIGINAL SPACE
i=1;
clear matlabbatch
if exist('T2','var') && ~isempty(T2) && ~skip_coreg % coregister anatomical images and segment
    matlabbatch{i}.spm.spatial.coreg.estwrite.ref = {T1};
    matlabbatch{i}.spm.spatial.coreg.estwrite.source = {T2};
    matlabbatch{i}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{i}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{i}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{i}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{i}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{i}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{i}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{i}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{i}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    i=i+1;
end
matlabbatch{i}.spm.spatial.preproc.channel(1).vols = {T1};
matlabbatch{i}.spm.spatial.preproc.channel(1).biasreg = 0.001;
matlabbatch{i}.spm.spatial.preproc.channel(1).biasfwhm = 60;
matlabbatch{i}.spm.spatial.preproc.channel(1).write = [1 1]; % save bias field as corrected input image
if exist('T2','var') && ~isempty(T2)  % add channel for coregistered (and resliced) T2
    if ~skip_coreg
        matlabbatch{i}.spm.spatial.preproc.channel(2).vols(1) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
        matlabbatch{i}.spm.spatial.preproc.channel(2).biasreg = T2biasreg; % 0.001 (default)
        matlabbatch{i}.spm.spatial.preproc.channel(2).biasfwhm = 60;
        matlabbatch{i}.spm.spatial.preproc.channel(2).write = [1 1];
    else
        matlabbatch{i}.spm.spatial.preproc.channel(2).vols(1) = {T2};
        matlabbatch{i}.spm.spatial.preproc.channel(2).biasreg = T2biasreg; % 0.001 (default)
        matlabbatch{i}.spm.spatial.preproc.channel(2).biasfwhm = 60;
        matlabbatch{i}.spm.spatial.preproc.channel(2).write = [1 1];
    end
end
matlabbatch{i}.spm.spatial.preproc.tissue(1).tpm = TPM(1);
matlabbatch{i}.spm.spatial.preproc.tissue(1).ngaus = 2;
matlabbatch{i}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{i}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{i}.spm.spatial.preproc.tissue(2).tpm = TPM(2);
matlabbatch{i}.spm.spatial.preproc.tissue(2).ngaus = 2;
matlabbatch{i}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{i}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{i}.spm.spatial.preproc.tissue(3).tpm = TPM(3);
matlabbatch{i}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{i}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{i}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{i}.spm.spatial.preproc.tissue(4).tpm = TPM(4);
matlabbatch{i}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{i}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{i}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{i}.spm.spatial.preproc.tissue(5).tpm = TPM(5);
matlabbatch{i}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{i}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{i}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{i}.spm.spatial.preproc.tissue(6).tpm = TPM(6);
matlabbatch{i}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{i}.spm.spatial.preproc.tissue(6).native = [1 0];
matlabbatch{i}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{i}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{i}.spm.spatial.preproc.warp.cleanup = 0;
matlabbatch{i}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{i}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{i}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{i}.spm.spatial.preproc.warp.samp = dsf;
% Save forward deformation field which maps mm coordinates from subject
% space to MNI space. It contains three images stored in a 4D nifti file
% each mapping a coordinate in subject space (x, y, and z, respectively)
% to a particular voxel in MNI space
matlabbatch{i}.spm.spatial.preproc.warp.write = [1 1];

matlabbatch{i+1}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{i}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{i+1}.spm.spatial.normalise.write.subj.resample = {T1};
matlabbatch{i+1}.spm.spatial.normalise.write.woptions.bb = [-90 -126 -190; 90 90 107];
matlabbatch{i+1}.spm.spatial.normalise.write.woptions.vox = vx_size;
matlabbatch{i+1}.spm.spatial.normalise.write.woptions.interp = 4;

% pushforward uses the inverse deformation field of pullback, hence
% here we use the forward deformation to go from MNI to subject space
matlabbatch{i+2}.spm.util.defs.comp{1}.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{i}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{i+2}.spm.util.defs.out{1}.push.fnames = {fullfile(templatesdir,'spmprior_tissue.nii')
						     fullfile(templatesdir,'spmprior_eye1.nii')
						     fullfile(templatesdir,'spmprior_eye2.nii')
						     fullfile(templatesdir,'spmprior_air.nii')
                                                     fullfile(templatesdir,'spmprior_ventricles_lateral.nii')
						     fullfile(templatesdir,'spmprior_spinal.nii')};
matlabbatch{i+2}.spm.util.defs.out{1}.push.weight = {''};
matlabbatch{i+2}.spm.util.defs.out{1}.push.savedir.saveusr = {fileparts(T1)};
matlabbatch{i+2}.spm.util.defs.out{1}.push.fov.file = {T1};
matlabbatch{i+2}.spm.util.defs.out{1}.push.preserve = 0;
matlabbatch{i+2}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
matlabbatch{i+2}.spm.util.defs.out{1}.push.prefix = '';

spm_jobman('run',matlabbatch)

disp('Getting transformation matrices from MNI to subject space...')

% Get and save 6 parameter (rigid body) transform by coregistering the 
% original T1 and the normalized T1
[path,name,ext] = fileparts(T1);
wT1 = spm_vol(fullfile(path,sprintf('w%s%s',name,ext)));

options.sep = [4 2];
options.params = zeros(1,6);
options.tol = [0.02 0.02 0.02 0.001 0.001 0.001];
options.cost_fun = 'nmi';
options.fwhm = [7 7];
options.graphics = spm('CmdLine');

mni2conform6 = spm_matrix(spm_coreg(wT1, T1, options));

% Get and save 12 parameter (affine) transform by taking only the affine
% part of the transformation from the unified segmentation procedure and
% invert it

s = load(fullfile(path,sprintf('%s_seg8.mat',name)));
mni2conform12 = inv(s.Affine);

% Using these matrices it is possible to go from voxel indices in one image
% to another. Specifically, to convert point p from voxel indices in image
% 1 to voxel indices in image 2, the following steps are necessary:
% (1) convert voxel indices of image 1 to world coordinates (in mm) by
%     multiplying with M1,
% (2) multiply by R to go from coordinates in image 1 to coordinates in
%     image 2,
% (3) and finally, multiply by inv(M2) to go from world coordinates to
%     voxel indices in image 2.
%
% Thus, p2 = inv(M2)*R*M1 * p1 = M2\R*M1 * p1 with the matrix M2\R*M1
% transforming from image 1 to image 2. Here R corresponds to either
% mni2conform6 or mni2conform12.

save(fullfile(path,'MNI2conform_6DOF.txt'), 'mni2conform6', '-ascii');
save(fullfile(path,'MNI2conform_12DOF.txt'),'mni2conform12','-ascii');

end

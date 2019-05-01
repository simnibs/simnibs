function segment_CAT(T1, template_dir,no_print)
% 
% T1 : T1 weighted image for segmentation with CAT12
% template_dir : directory with templates to warp
% 

% hide all figures

if nargin > 1
    deform_templates = true;
else
    deform_templates = false;
end

% try to get CAT defaults
try
    cat = cat_get_defaults();
catch
    try
        addpath(fullfile(spm('dir'),'toolbox','cat12'))
        cat = cat_get_defaults();
    catch ME
        disp('Could not get CAT12 defaults. Is CAT12 in toolbox/cat12?')
        %throw(ME)
        exit(2);
    end
end
cat_get_defaults('output.CSF.native', true);
cat_get_defaults('output.bias.native', true);
cat_get_defaults('extopts.subfolders', true);
if(no_print)
	cat_get_defaults('extopts.print',0);
end
%cat_get_defaults('output.ROI', 0);
%cat_get_defaults('extopts.verb',1);
%cat_get_defaults('extopts.debug',0);

[path, name, ext] = fileparts(T1);

% initialize SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');
spm_get_defaults('cmdline', true);

% segment using CAT12
matlabbatch = [];
matlabbatch{1}.spm.tools.cat.estwrite.data = {T1};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 0; % block
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = cat.opts.tpm;
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.cleanupstr = 0.5;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.darteltpm = cat.extopts.darteltpm;
matlabbatch{1}.spm.tools.cat.extopts.print = 0; 
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.ignoreErrors = 1;
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.debug = 0;
%matlabbatch{1}.spm.tools.cat.estwrite.extopts.verb = 0;

matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROI = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobian.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 0]; % [forward inverse]
spm_jobman('run', matlabbatch);
%close all

if deform_templates
    forward_def = fullfile(path, 'mri', strcat('y_', name, ext));
    
    % the templates to transform
    brainstem  = fullfile(template_dir, 'cattemplate_CC_hammers.nii');
    CC         = fullfile(template_dir, 'cattemplate_brainstem_hammers.nii');
    cerebellum = fullfile(template_dir, 'cattemplate_cerebellum_hammers.nii');
    fornix     = fullfile(template_dir, 'cattemplate_fornix_JHU.nii');
    vlateral   = fullfile(template_dir, 'cattemplate_ventricles_lateral_hammers.nii');
    vthird     = fullfile(template_dir, 'cattemplate_ventricles_third_hammers.nii');
    thalamus   = fullfile(template_dir, 'cattemplate_thalamus_harvard.nii');
    
    % Deformation of templates
    matlabbatch = [];
    matlabbatch{1}.spm.util.defs.comp{1}.def = {forward_def};
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames = {brainstem
                                                       CC
                                                       cerebellum
                                                       fornix
                                                       vlateral
                                                       vthird
                                                       thalamus
                                                       cat.extopts.darteltpm{1}
                                                       };
    matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {path};
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {T1};
    matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
    matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix = '';
    spm_jobman('run', matlabbatch);
end
end

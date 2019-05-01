function nonlin_coreg = cat_conf_nonlin_coreg
% Configuration file for non-linear co-registration
%
% Christian Gaser
% $Id: cat_conf_nonlin_coreg.m 1220 2017-11-20 13:34:33Z gaser $

%--------------------------------------------------------------------------
% Reference Image
%--------------------------------------------------------------------------
ref         = cfg_files;
ref.tag     = 'ref';
ref.name    = 'Reference Image';
ref.help    = {'This is the image that is assumed to remain stationary (e.g. T1 image), while the source image is moved to match it.'};
ref.filter  = 'image';
ref.ufilter = '.*';
ref.num     = [1 1];
ref.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% Source Image
%--------------------------------------------------------------------------
source         = cfg_files;
source.tag     = 'source';
source.name    = 'Source Image';
source.help    = {'This is the image that is jiggled about to best match the reference (e.g. mean EPI, B0 image).'};
source.filter  = 'image';
source.ufilter = '.*';
source.num     = [1 1];
source.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% Other Images
%--------------------------------------------------------------------------
other = cfg_files;
other.name = 'Other images to write';
other.tag  = 'other';
other.filter = 'image';
other.num  = [1 Inf];
other.help    = {'These are any images that need to remain in alignment with the source image.'};
other.preview = @(f) spm_image('Display',char(f));

%--------------------------------------------------------------------------
% Warping regularisation
%--------------------------------------------------------------------------
reg         = cfg_entry;
reg.tag     = 'reg';
reg.name    = 'Warping Regularisation';
reg.help    = {'The objective function for registering the tissue probability maps to the image to process, involves minimising the sum of two terms. One term gives a function of how probable the data is given the warping parameters. The other is a function of how probable the parameters are, and provides a penalty for unlikely deformations. Smoother deformations are deemed to be more probable. The amount of regularisation determines the tradeoff between the terms. Start with a value around one.  However, if your normalised images appear distorted, then it may be an idea to increase the amount of regularisation (by an order of magnitude). More regularisation gives smoother deformations, where the smoothness measure is determined by the bending energy of the deformations. '};
reg.strtype = 'r';
reg.num     = [1 1];
reg.val     = {1};

%--------------------------------------------------------------------------
% Bounding box
%--------------------------------------------------------------------------
bb         = cfg_entry;
bb.tag     = 'bb';
bb.name    = 'Bounding box';
bb.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure). NaN is using the bounding box of the reference image.'};
bb.strtype = 'r';
bb.num     = [2 3];
bb.val     = {[NaN NaN NaN; NaN NaN NaN]};

%--------------------------------------------------------------------------
% Voxel sizes
%--------------------------------------------------------------------------
vox         = cfg_entry;
vox.tag     = 'vox';
vox.name    = 'Voxel sizes';
vox.help    = {'The voxel sizes (x, y & z, in mm) of the written normalised images.'};
vox.strtype = 'r';
vox.num     = [1 3];
vox.def     = @(val)spm_get_defaults('normalise.write.vox', val{:});

%------------------------------------------------------------------------

nonlin_coreg = cfg_exbranch;
nonlin_coreg.name = 'Non-linear co-registration';
nonlin_coreg.tag  = 'nonlin_coreg';
nonlin_coreg.val  = {ref,source,other,reg,bb,vox};
nonlin_coreg.prog = @cat_vol_nonlin_coreg_multi_run;
nonlin_coreg.vout = @vout_nonlin_coreg;
nonlin_coreg.help = {
    'Within-subject non-linear co-registration performed via the segmentation routine.'
    ''
    'The non-linear co-registration method used here is based on the SPM12 non-linear normalisation method that uses the segmented images to estimate deformations to match the source to the target image.'
    ''
    'The resliced images are named the same as the originals except that they are prefixed by ''w''.'
''
};

%------------------------------------------------------------------------

return;
%------------------------------------------------------------------------
 
%------------------------------------------------------------------------
function dep = vout_nonlin_coreg(job)

dep(1)            = cfg_dep;
dep(1).sname      = 'Non-linear coregistered data';
dep(1).src_output = substruct('.','ofiles');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%------------------------------------------------------------------------

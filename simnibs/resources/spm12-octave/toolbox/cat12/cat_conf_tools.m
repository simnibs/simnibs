function tools = cat_conf_tools(expert)
% wrapper for calling CAT utilities
%
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_conf_tools.m 1275 2018-02-12 22:00:05Z gaser $

%_______________________________________________________________________


data_T2x         = cfg_files;
data_T2x.tag     = 'data_T2x';
data_T2x.name    = 'Data';
data_T2x.filter  = {'image'};
data_T2x.ufilter = '^spmT.*';
data_T2x.num     = [1 Inf];
data_T2x.help    = {'Select spmT-data to transform or convert.'};

sel        = cfg_menu;
sel.name   = 'Convert t value to';
sel.tag    = 'sel';
sel.labels = {'p','-log(p)','correlation coefficient cc','effect size d','apply thresholds without conversion'};
sel.values = {1,2,3,4,5};
sel.val    = {2};
sel.help   = {'Select conversion of t-value'};

thresh05         = cfg_entry;
thresh05.tag     = 'thresh05';
thresh05.name    = 'Threshold';
thresh05.help    = {''};
thresh05.strtype = 'r';
thresh05.num     = [1 1];
thresh05.val     = {0.05};

thresh001         = cfg_entry;
thresh001.tag     = 'thresh001';
thresh001.name    = 'Threshold';
thresh001.help    = {''};
thresh001.strtype = 'r';
thresh001.num     = [1 1];
thresh001.val     = {0.001};

kthresh         = cfg_entry;
kthresh.tag     = 'kthresh';
kthresh.name    = 'Extent (voxels)';
kthresh.help    = {'Enter the extent threshold in voxels'};
kthresh.strtype = 'r';
kthresh.val     = {0};
kthresh.num     = [1 1];

noniso        = cfg_menu;
noniso.name   = 'Correct for non-isotropic smoothness';
noniso.tag    = 'noniso';
noniso.labels = {'Yes','No'};
noniso.values = {1,0};
noniso.val    = {1};
noniso.help  = {'Correct for non-isotropic smoothness for cluster extent thresholds.'};

none         = cfg_const;
none.tag     = 'none';
none.name    = 'None';
none.val     = {1};
none.help    = {'No threshold'};

k         = cfg_branch;
k.tag     = 'k';
k.name    = 'k-value';
k.val     = {kthresh, noniso };
k.help    = {''};

fwe         = cfg_branch;
fwe.tag     = 'fwe';
fwe.name    = 'FWE';
fwe.val     = {thresh05 };
fwe.help    = {''};

fdr         = cfg_branch;
fdr.tag     = 'fdr';
fdr.name    = 'FDR';
fdr.val     = {thresh05 };
fdr.help    = {''};

fwe2         = cfg_branch;
fwe2.tag     = 'fwe2';
fwe2.name    = 'FWE';
fwe2.val     = {thresh05, noniso };
fwe2.help    = {''};

uncorr         = cfg_branch;
uncorr.tag     = 'uncorr';
uncorr.name    = 'uncorrected';
uncorr.val     = {thresh001 };
uncorr.help    = {''};

kuncorr         = cfg_branch;
kuncorr.tag     = 'kuncorr';
kuncorr.name    = 'uncorrected';
kuncorr.val     = {thresh05, noniso };
kuncorr.help    = {''};

En         = cfg_branch;
En.tag     = 'En';
En.name    = 'Expected voxels per cluster';
En.val     = {noniso };
En.help    = {''};

inverse        = cfg_menu;
inverse.name   = 'Show also inverse effects (e.g. neg. values)';
inverse.tag    = 'inverse';
inverse.labels = {'Yes','No'};
inverse.values = {1,0};
inverse.val    = {0};
inverse.help   = {'Show also inverse effects (e.g. neg. values). This is not valid if you convert to (log) p-values.'};

threshdesc        = cfg_choice;
threshdesc.name   = 'Threshold type peak-level';
threshdesc.tag    = 'threshdesc';
threshdesc.values = {none uncorr fdr fwe};
threshdesc.val    = {uncorr};
threshdesc.help   = {'Select method for voxel threshold'};

cluster        = cfg_choice;
cluster.name   = 'Cluster extent threshold';
cluster.tag    = 'cluster';
cluster.values = {none k En kuncorr fwe2};
cluster.val    = {none};
cluster.help   = {'Select method for extent threshold'};

conversion         = cfg_branch;
conversion.tag     = 'conversion';
conversion.name    = 'Conversion';
conversion.val     = {sel threshdesc inverse cluster};
conversion.help    = {''};

atlas        = cfg_menu;
atlas.name   = 'Atlas Labeling';
atlas.tag    = 'atlas';
list = spm_atlas('List','installed');
atlas.labels{1} = 'None';
atlas.values{1} = 'None';
for i=1:numel(list)
  atlas.labels{i+1} = list(i).name;
  atlas.values{i+1} = list(i).name;
end
atlas.val    = {'None'};
atlas.help   = {'Select atlas for labeling. The prepending atlas name ''dartel_'' indicates that this atlas was created using Dartel spatial normalization with the Dartel IXI template as default.'
''
'Please note, that you can install additional atlases for CAT12 using the command ''cat_install_atlases''. '};

T2x      = cfg_exbranch;
T2x.tag  = 'T2x';
T2x.name = 'Threshold and transform spmT images';
T2x.val  = {data_T2x,conversion,atlas};
T2x.prog = @cat_stat_spm2x;
T2x.vout = @vout_stat_spm2x;
T2x.help = {
          'This function transforms t-maps to P, -log(P), r or d-maps.'
          'The following formulas are used:'
          '--------------------------------'
          'correlation coefficient:'
          '          sign(t)'
          'r = ------------------'
          '           df'
          '    sqrt(------ + 1)'
          '          t*t'
          'effect-size'
          '           2r'
          'd = ----------------'
          '    sqrt(1-sqr(r))'
          'p-value'
          'p = 1-spm_Tcdf'
          'log p-value'
          '-log10(1-P) = -log(1-spm_Tcdf)'
          'For the last case of log transformation this means that a p-value of p=0.99 (0.01) is transformed to a value of 2.'
          'Examples:'
          'p-value  -log10(1-P)'
          '0.1      1'
          '0.05     1.3'
          '0.01     2'
          '0.001    3'
          '0.0001   4'
          'All maps can be thresholded using height and extent thresholds and you can also apply corrections for multiple comparisons based on family-wise error (FWE) or false discovery rate (FDR). You can easily threshold and/or transform a large number of spmT-maps using the same thresholds.'
          'Naming convention of the transformed files:'
          '   Type_Contrast_Pheight_Pextent_K_Neg'
          '   Type:      P    - p-value'
          '              logP - log p-value'
          '              R    - correlation coefficient'
          '              D    - effect size'
          '              T    - t-value'
          '   Contrast:  name used in the contrast manager with replaced none valid'
          '              strings'
          '   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")'
          '              pFWE - p-value with FWE correction in %'
          '              pFDR - p-value with FDR correction in %'
          '   Pextent:   pk    - uncorr. extent p-value in % (p<0.05 coded with "p5")'
          '              pkFWE - extent p-value with FWE correction in %'
          '   K:         extent threshold in voxels'
          '   Neg:       image also shows thresholded inverse effects (e.g. neg. '
          '              values) '
}';

%------------------------------------------------------------------------
% Do not use 3D atlases for surfaces
data_T2x.filter  = {'gifti'};
T2x_surf      = T2x;
T2x_surf.val  = {data_T2x,conversion};
T2x_surf.tag  = 'T2x_surf';
T2x_surf.name = 'Threshold and transform spmT surfaces';
T2x_surf.vout = @vout_stat_spm2x_surf;

%------------------------------------------------------------------------

data_F2x         = cfg_files;
data_F2x.tag     = 'data_F2x';
data_F2x.name    = 'Data';
data_F2x.filter  = {'image'};
data_F2x.ufilter = '^spmF.*';
data_F2x.num     = [1 Inf];
data_F2x.help    = {'Select spmF-data to select.'};

sel        = cfg_menu;
sel.name   = 'Convert F value to';
sel.tag    = 'sel';
sel.labels = {'p','-log(p)','coefficient of determination R^2','apply thresholds without conversion'};
sel.values = {1,2,3,4};
sel.val    = {2};
sel.help   = {'Select conversion of F-value'};

none         = cfg_const;
none.tag     = 'none';
none.name    = 'None';
none.val     = {1};
none.help    = {'No threshold'};

cluster        = cfg_choice;
cluster.name   = 'Cluster extent threshold';
cluster.tag    = 'cluster';
cluster.values = {none k En kuncorr fwe2};
cluster.val    = {none};
cluster.help  = {'Select method for extent threshold'};

conversion         = cfg_branch;
conversion.tag     = 'conversion';
conversion.name    = 'Conversion';
conversion.val     = {sel threshdesc cluster};
conversion.help    = {''};

F2x      = cfg_exbranch;
F2x.tag  = 'F2x';
F2x.name = 'Threshold and transform spmF images';
F2x.val  = {data_F2x,conversion,atlas};
F2x.prog = @cat_stat_spm2x;
F2x.vout = @vout_stat_spm2x;
F2x.help = {
          'This function transforms F-maps to P, -log(P), or R2-maps.'
          'The following formulas are used:'
          '--------------------------------'
          'coefficient of determination R2:'
          '          F*(n-1)'
          'R2 = ------------------'
          '        n-p + F*(n-1)'
          'p-value:'
          'p = 1-spm_Fcdf'
          'log p-value:'
          '-log10(1-P) = -log(1-spm_Fcdf)'
          'For the last case of log transformation this means that a p-value of p=0.99 (0.01) is transformed to a value of 2.'
          'Examples:'
          'p-value  -log10(1-P)'
          '0.1      1'
          '0.05     1.3'
          '0.01     2'
          '0.001    3'
          '0.0001   4'
          'All maps can be thresholded using height and extent thresholds and you can also apply corrections for multiple comparisons based on family-wise error (FWE) or false discovery rate (FDR). You can easily threshold and/or transform a large number of spmT-maps using the same thresholds.'
          'Naming convention of the transformed files:'
          '   Type_Contrast_Pheight_K'
          '   Type:      P    - p-value'
          '              logP - log p-value'
          '              R2   - coefficient of determination'
          '   Contrast:  name used in the contrast manager with replaced none valid'
          '              strings'
          '   Pheight:   p    - uncorrected p-value in % (p<0.05 will coded with "p5")'
          '              pFWE - p-value with FWE correction in %'
          '              pFDR - p-value with FDR correction in %'
          '   K:         extent threshold in voxels'
}';

%------------------------------------------------------------------------
% Do not use 3D atlases for surfaces
data_F2x.filter  = {'gifti'};
F2x_surf      = F2x;
F2x_surf.val  = {data_F2x,conversion};
F2x_surf.tag  = 'F2x_surf';
F2x_surf.name = 'Threshold and transform spmF surfaces';
F2x_surf.vout = @vout_stat_spm2x_surf;

%------------------------------------------------------------------------

c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {'Vector of nuisance values'};
c.strtype = 'r';
c.num     = [Inf 1];

slice         = cfg_entry;
slice.tag     = 'slice';
slice.name    = 'Selected slice (in mm)?';
slice.strtype = 'r';
slice.num     = [1 1];
slice.val     = {0};
slice.help    = {'Choose slice in mm.'};

gap         = cfg_entry;
gap.tag     = 'gap';
gap.name    = 'Separation';
gap.strtype = 'n';
gap.num     = [1 1];
gap.val     = {3};
gap.help    = {
  'To speed up calculations you can define that covariance is estimated only every x voxel. Smaller values give slightly more accurate covariance, but will be much slower.'};

scale        = cfg_menu;
scale.tag    = 'scale';
scale.name   = 'Proportional scaling?';
scale.labels = {'No','Yes'};
scale.values = {0 1};
scale.val    = {0};
scale.help   = {'This option should be only used if image intensity is not scaled (e.g. T1 images) or if images have to be scaled during statistical analysis (e.g. modulated images).'};

nuisance         = cfg_repeat;
nuisance.tag     = 'nuisance';
nuisance.name    = 'Nuisance variable';
nuisance.values  = {c};
nuisance.num     = [0 Inf];
nuisance.help    = {'This option allows for the specification of nuisance effects to be removed from the data. A potential nuisance parameter can be TIV if you check segmented data with the default modulation. In this case the variance explained by TIV will be removed prior to the calculation of the correlation. Another meaningful nuisance effect is age.'};

data_xml = cfg_files;
data_xml.name = 'Quality measures (optional)';
data_xml.tag  = 'data_xml';
data_xml.filter = 'xml';
data_xml.ufilter = '^cat_.*';
data_xml.val  = {{''}};
data_xml.num  = [0 Inf];
data_xml.help   = {...
'Select optional the quality measures that are saved during segmentation as xml-files in the report folder. This additionally allows to analyze image quality parameters such as noise, and bias. Please note, that the order of the xml-files should be the same as the other data files.'};

data_vol = cfg_files;
data_vol.name = 'Sample data';
data_vol.tag  = 'data_vol';
data_vol.filter = 'image';
data_vol.num  = [1 Inf];
data_vol.help   = {...
'These are the (spatially registered) data. They must all have the same image dimensions, orientation, voxel size etc. Furthermore, it is recommended to use unsmoothed files.'};

sample         = cfg_repeat;
sample.tag     = 'sample';
sample.name    = 'Data';
sample.values  = {data_vol };
sample.num     = [1 Inf];
sample.help = {...
'Specify data for each sample. If you specify different samples the mean correlation is displayed in separate boxplots (or violin plots) for each sample.'};

 
check_cov      = cfg_exbranch;
check_cov.tag  = 'check_cov';
check_cov.name = 'Check sample homogeneity of 3D data';
check_cov.val  = {sample,data_xml,gap,nuisance};
check_cov.prog = @cat_stat_check_cov;
check_cov.help  = {
'In order to identify images with poor image quality or even artefacts you can use this function. Images have to be in the same orientation with same voxel size and dimension (e.g. normalized images without smoothing). The idea of this tool is to check the correlation of all files across the sample.'
''
'The correlation is calculated between all images and the mean for each image is plotted using a spmmat and the indicated filenames. The smaller the mean correlation the more deviant is this image from the sample mean. In the plot outliers from the sample are usually isolated from the majority of images which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the image order.'};

%------------------------------------------------------------------------

spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.filter  = {'mat'};
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];
spmmat.help    = {'Select the SPM.mat file that contains the design specification.'};

use_unsmoothed_data        = cfg_menu;
use_unsmoothed_data.name   = 'Use unsmoothed data if found';
use_unsmoothed_data.tag    = 'use_unsmoothed_data';
use_unsmoothed_data.labels = {'Yes','No'};
use_unsmoothed_data.values = {1,0};
use_unsmoothed_data.val    = {1};
use_unsmoothed_data.help  = {'Check for sample homogeneity results in more reliable values if unsmoothed data are used. Unsmoothed data contain more detailed information about differences and similarities between the data.'};

adjust_data        = cfg_menu;
adjust_data.name   = 'Adjust data using design matrix';
adjust_data.tag    = 'adjust_data';
adjust_data.labels = {'Yes','No'};
adjust_data.values = {1,0};
adjust_data.val    = {1};
adjust_data.help  = {'This option allows to use nuisance and group parameters from the design matrix to obtain adjusted data. In this case the variance explained by these parameters will be removed prior to the calculation of the correlation. Furthermore, global scaling (if defined) is also applied to the data.'};

none         = cfg_const;
none.tag     = 'none';
none.name    = 'No';
none.val     = {1};
none.help    = {''};

do_check_cov         = cfg_branch;
do_check_cov.tag     = 'do_check_cov';
do_check_cov.name    = 'Yes';
do_check_cov.val     = {use_unsmoothed_data adjust_data};
do_check_cov.help    = {''};

check_SPM_cov        = cfg_choice;
check_SPM_cov.name   = 'Check for sample homogeneity';
check_SPM_cov.tag    = 'check_SPM_cov';
check_SPM_cov.values = {none do_check_cov};
check_SPM_cov.val    = {do_check_cov};
check_SPM_cov.help   = {'In order to identify images with poor image quality or even artefacts you can use this function. The idea of this tool is to check the correlation of all files across the sample using the files that are already defined in SPM.mat.'
''
'The correlation is calculated between all images and the mean for each image is plotted using a boxplot (or violin plot) and the indicated filenames. The smaller the mean correlation the more deviant is this image from the sample mean. In the plot outliers from the sample are usually isolated from the majority of images which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the image order'};

check_SPM_ortho        = cfg_menu;
check_SPM_ortho.name   = 'Check for design orthogonality';
check_SPM_ortho.tag    = 'check_SPM_ortho';
check_SPM_ortho.labels = {'Yes','No'};
check_SPM_ortho.values = {1,0};
check_SPM_ortho.val    = {1};
check_SPM_ortho.help   = {'Review Design Orthogonality.'};
 
check_SPM      = cfg_exbranch;
check_SPM.tag  = 'check_SPM';
check_SPM.name = 'Check design orthogonality and homogeneity';
check_SPM.val  = {spmmat,check_SPM_cov,check_SPM_ortho};
check_SPM.prog = @cat_stat_check_SPM;
check_cov.help  = {
'Use design matrix saved in SPM.mat to check for sample homogeneity of the used data and for orthogonality of parameters.'
};

%------------------------------------------------------------------------

outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];
outdir.help    = {'Select a directory where files are written.'};


data         = cfg_files; 
data.tag     = 'data';
data.name    = 'Volumes';
data.filter  = 'image';
data.ufilter = '.*';
data.num     = [1 Inf];
data.help = {
'Select images for quality assurance.'};

qa        = cfg_exbranch;
qa.tag    = 'qa';
qa.name   = 'CAT quality assurance';
qa.val    = {data};
qa.prog   = @cat_vol_qa;
qa.vfiles = @vfiles_qa;
qa.help   = {'CAT Quality Assurance of T1 images. '};

%------------------------------------------------------------------------

data.help = {
'Select all images. Images have to be in the same orientation with same voxel size and dimension (e.g. normalized images)'};

orient        = cfg_menu;
orient.tag    = 'orient';
orient.name   = 'Spatial orientation';
orient.labels = {'axial','coronal','sagittal'};
orient.values = {3 2 1};
orient.val    = {3};
orient.help   = {'Spatial orientation of slice.'};

showslice      = cfg_exbranch;
showslice.tag  = 'showslice';
showslice.name = 'Display one slice for all images';
showslice.val  = {data_vol,scale,orient,slice};
showslice.prog = @cat_stat_showslice_all;
showslice.help = {'This function displays a selected slice for all images and indicates the respective filenames which is useful to check image quality for a large number of files in a circumscribed region (slice).'};

%------------------------------------------------------------------------

data.help = {
'Select images for filtering'};

rician         = cfg_menu;
rician.tag     = 'rician';
rician.name    = 'Rician noise';
rician.labels  = {'Yes' 'No'};
rician.values  = {1 0};
rician.val     = {0};
rician.help    = {'MRIs can have Gaussian or Rician distributed noise with uniform or nonuniform variance across the image. If SNR is high enough (>3) noise can be well approximated by Gaussian noise in the foreground. However, for SENSE reconstruction or DTI data a Rician distribution is expected.'
''
'Please note that the Rician noise estimation is sensitive for large signals in the neighbourhood and can lead to artefacts (e.g. cortex can be affected by very high values in the scalp or in blood vessels.'
''
};

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'sanlm_'};
prefix.help    = {'Specify the string to be prepended to the filenames of the smoothed image file(s). Default prefix is ''samlm_''.'};

NCstr         = cfg_entry;
NCstr.tag     = 'NCstr';
NCstr.name    = 'Strength of Noise Corrections';
NCstr.strtype = 'r';
NCstr.num     = [1 1];
NCstr.def     = @(val)cat_get_defaults('extopts.NCstr', val{:});
NCstr.help    = {
  'Strength of the SANLM noise correction. The default "-inf" uses an adaptive noise correction and was successfully tested on a variety of scans. Use smaller values (>0) for less denoising  and higher values (<=1) for stronger denoising. The value 0 will turn off any noise correction.'
''
};

sanlm        = cfg_exbranch;
sanlm.tag    = 'sanlm';
sanlm.name   = 'Spatially adaptive non-local means denoising filter';
sanlm.val    = {data prefix NCstr rician};
sanlm.prog   = @cat_vol_sanlm;
sanlm.vfiles = @vfiles_sanlm;
sanlm.help   = {
'This function applies an spatial adaptive non-local means denoising filter to the data. This filter will remove noise while preserving edges. The filter strength is automatically estimated based on the standard deviation of the noise. '
''
'This filter is internally used in the segmentation procedure anyway. Thus, it is not neccessary (and not recommended) to apply the filter before segmentation.'
''
};

%------------------------------------------------------------------------
roi_xml = cfg_files;
roi_xml.name = 'XML files';
roi_xml.tag  = 'roi_xml';
roi_xml.filter = 'xml';
roi_xml.ufilter = '^catROI.*';
roi_xml.num  = [1 Inf];
roi_xml.help   = {...
'These are the xml-files that are saved in the label folder.'};

usefolder         = cfg_menu;
usefolder.tag     = 'folder';
usefolder.name    = 'Use foldername';
usefolder.labels  = {'Yes' 'No'};
usefolder.values  = {1 0};
usefolder.val     = {0};
usefolder.help    = {
'Use foldername to describe the subject.'};

point         = cfg_menu;
point.tag     = 'point';
point.name    = 'Decimal point';
point.labels  = {',','.'};
point.values  = {',','.'};
point.val     = {'.'};
point.help    = {
'Decimal point.'};  % that has to be unequal to the column delimiter.'};

% tab "\t" does not work and so we automatically switch in case of decimal 
% point "," to delimiter ";".
%{
delimiter         = cfg_menu;
delimiter.tag     = 'delimiter';
delimiter.name    = 'column delimiter';
delimiter.labels  = {',',';',' '};
delimiter.values  = {',',';',' '};
delimiter.val     = {','};
delimiter.help    = {
'Delimiter between columns.'};
%}

calcroi_name         = cfg_entry;
calcroi_name.tag     = 'calcroi_name';
calcroi_name.name    = 'Output file';
calcroi_name.strtype = 's';
calcroi_name.num     = [1 Inf];
calcroi_name.val     = {'ROI'};
calcroi_name.help    = {
'The output file is written to the current working directory unless a valid full pathname is given. The output file will also include the name of the atlas and the measure (e.g. Vgm). The file is using tabstops to separate values in order to easily import the file into Excel or SPSS or any other software for subsequent analysis.'};

calcroi       = cfg_exbranch;
calcroi.tag   = 'calcroi';
calcroi.name  = 'Estimate mean values inside ROI';
calcroi.val   = {roi_xml,point,outdir,calcroi_name}; 
%calcroi.val   = {roi_xml,usefolder,point,outdir,calcroi_name}; % usefolder is never used
calcroi.prog  = @(job)cat_roi_fun('exportSample',job);
calcroi.help  = {
'This function reads mean values inside a ROIs from different atlases and saves values for all data in a csv-file. '
'Missed values were replaced by NaN.'
};

%------------------------------------------------------------------------
calcvol_name         = cfg_entry;
calcvol_name.tag     = 'calcvol_name';
calcvol_name.name    = 'Output file';
calcvol_name.strtype = 's';
calcvol_name.num     = [1 Inf];
calcvol_name.val     = {'TIV.txt'};
calcvol_name.help    = {
'The output file is written to current working directory unless a valid full pathname is given.'};

calcvol_TIV         = cfg_menu;
calcvol_TIV.tag     = 'calcvol_TIV';
calcvol_TIV.name    = 'Save values';
calcvol_TIV.labels  = {'TIV only' 'TIV & GM/WM/CSF/WMH'};
calcvol_TIV.values  = {1 0};
calcvol_TIV.val     = {1};
calcvol_TIV.help    = {'You can save either the total intracranial volume (TIV) only or additionally also save the global volumes for GM, WM, CSF, and WM hyperintensities.'
''
};

clear data_xml
data_xml = cfg_files;
data_xml.name = 'XML files';
data_xml.tag  = 'data_xml';
data_xml.filter = 'xml';
data_xml.ufilter = '^cat_.*';
data_xml.num  = [1 Inf];
data_xml.help   = {...
'Select xml-files that are saved during segmentation in the report folder.'};

calcvol       = cfg_exbranch;
calcvol.tag   = 'calcvol';
calcvol.name  = 'Estimate TIV and global tissue volumes';
calcvol.val   = {data_xml,calcvol_TIV,calcvol_name};
calcvol.prog  = @cat_stat_TIV;
calcvol.vout  = @vout_stat_TIV;
calcvol.help  = {
'This function reads raw volumes for TIV/GM/WM/CSF/WM hyperintensities (WMH) and saves values in a txt-file. These values can be read with the matlab command: vol = spm_load. If you choode to save all values the entries for TIV/GM/WM/CSF/WMH are now saved in vol(:,1) vol(:,2) vol(:,3), vol(:,4), and vol(:,5) respectively.'
''
'You can use TIV either as nuisance in an AnCova model or as user-specified globals with the "global calculation" option depending on your hypothesis. The use of TIV as nuisance or globals is recommended for modulated data where both the affine transformation and the non-linear warping of the registration are corrected for. '
''
};

%------------------------------------------------------------------------
iqr_name         = cfg_entry;
iqr_name.tag     = 'iqr_name';
iqr_name.name    = 'Output file';
iqr_name.strtype = 's';
iqr_name.num     = [1 Inf];
iqr_name.val     = {'IQR.txt'};
iqr_name.help    = {
'The output file is written to current working directory unless a valid full pathname is given'};

iqr       = cfg_exbranch;
iqr.tag   = 'iqr';
iqr.name  = 'Get Weighted Overall Image Quality';
iqr.val   = {data_xml,iqr_name};
iqr.prog  = @cat_stat_IQR;
iqr.help  = {
'This function reads weighted overall image quality from saved xml-files. '
''
};

%------------------------------------------------------------------------

field         = cfg_files;
field.tag     = 'field';
field.name    = 'Deformation Field';
field.filter  = 'image';
field.ufilter = '^(i)?y_.*\.nii$'; % '.*y_.*\.nii$';
field.num     = [1 Inf];
field.help    = {[
'Deformations can be thought of as vector fields. These can be represented by three-volume images.' ...
'Use the "y_*.nii" to project data from subject to template space, and the "iy_*.nii" to map data from template to individual space.' ...
'Both deformation maps can be created in the CAT preprocessing by setting the "Deformation Field" flag (no written by default).' ... 
]};

field1         = cfg_files;
field1.tag     = 'field1';
field1.name    = 'Deformation Field';
field1.filter  = 'image';
field1.ufilter = '^(i)?y_.*\.nii$'; % '.*y_.*\.nii$';
field1.num     = [1 1];
field1.help    = {[
'Deformations can be thought of as vector fields. These can be represented by three-volume images.' ...
'Use the "y_*.nii" to project data from subject to template space, and the "iy_*.nii" to map data from template to individual space.' ...
'Both deformation maps can be created in the CAT preprocessing by setting the "Deformation Field" flag (no written by default).' ... 
]};

images1         = cfg_files;
images1.tag     = 'images';
images1.name    = 'Images';
images1.help    = {'Select images to be warped. Note that there should be the same number of images as there are deformation fields, such that each flow field warps one image.'};
images1.filter = 'image';
images1.ufilter = '.*';
images1.num     = [1 Inf];

images         = cfg_repeat;
images.tag     = 'images';
images.name    = 'Images';
images.help    = {'The flow field deformations can be applied to multiple images. At this point, you are choosing how many images each flow field should be applied to.'};
images.values  = {images1 };
images.num     = [1 Inf];

interp        = cfg_menu;
interp.name   = 'Interpolation';
interp.tag    = 'interp';
interp.labels = {'Nearest neighbour','Trilinear','2nd Degree B-spline',...
'3rd Degree B-Spline ','4th Degree B-Spline ','5th Degree B-Spline',...
'6th Degree B-Spline','7th Degree B-Spline'};
interp.values = {0,1,2,3,4,5,6,7};
interp.val    = {5};
interp.help   = {
'The method by which the images are sampled when being written in a different space.'
'    Nearest Neighbour:     - Fastest, but not normally recommended.'
'    Bilinear Interpolation:     - OK for PET, or realigned fMRI.'
'    B-spline Interpolation:     - Better quality (but slower) interpolation/* \cite{thevenaz00a}*/, especially with higher degree splines.  Do not use B-splines when there is any region of NaN or Inf in the images. '
}';

modulate        = cfg_menu;
modulate.tag    = 'modulate';
modulate.name   = 'Modulate image (preserve volume)';
modulate.labels = {'No','Affine + non-linear (SPM12 default)','Non-linear only'};
modulate.values = {0 1 2};
modulate.val    = {0};
modulate.help = {
'"Modulation" is to compensate for the effect of spatial normalisation. Spatial normalisation causes volume changes due to affine transformation (global scaling) and non-linear warping (local volume change). The SPM default is to adjust spatially normalised grey matter (or other tissue class) by using both terms and the resulting modulated images are preserved for the total amount of grey matter. Thus, modulated images reflect the grey matter volumes before spatial normalisation. However, the user is often interested in removing the confound of different brain sizes and there are many ways to apply this correction. We can use the total amount of GM, GM+WM, GM+WM+CSF, or manual estimated total intracranial volume (TIV). Theses parameters can be modeled as nuisance parameters (additive effects) in an AnCova model or used to globally scale the data (multiplicative effects): '
''
'% Correction   Interpretation'
'% ----------   --------------'
'% nothing      absolute volume'
'% globals 	    relative volume after correcting for total GM or TIV (multiplicative effects)'
'% AnCova 	    relative volume that can not be explained by total GM or TIV (additive effects)'
''
'Modulated images can be optionally saved by correcting for non-linear warping only. Volume changes due to affine normalisation will be not considered and this equals the use of default modulation and globally scaling data according to the inverse scaling factor due to affine normalisation. I recommend this option if your hypothesis is about effects of relative volumes which are corrected for different brain sizes. This is a widely used hypothesis and should fit to most data. The idea behind this option is that scaling of affine normalisation is indeed a multiplicative (gain) effect and we rather apply this correction to our data and not to our statistical model. These modulated images are indicated by "m0" instead of "m". '
''
};

defs        = cfg_exbranch;
defs.tag    = 'defs';
defs.name   = 'Apply deformations (many images)';
defs.val    = {field1,images1,interp,modulate};
defs.prog   = @cat_vol_defs;
defs.vfiles = @vfiles_defs;
defs.help   = {'This is a utility for applying a deformation field of one subject to many images.'};

defs2        = cfg_exbranch;
defs2.tag    = 'defs2';
defs2.name   = 'Apply deformations (many subjects)';
defs2.val    = {field,images,interp,modulate};
defs2.prog   = @cat_vol_defs;
defs2.vfiles = @vfiles_defs2;
defs2.help   = {'This is a utility for applying deformation fields of many subjects to images.'};

data.help = {
'Select all images for this subject'};

bparam         = cfg_entry;
bparam.tag     = 'bparam';
bparam.name    = 'Bias Regularisation';
bparam.help    = {'MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images.'
                   ''
                   'An important issue relates to the distinction between variations in the difference between the images that arise because of the differential bias artifact due to the physics of MR scanning, and those that arise due to shape differences.  The objective is to model the latter by deformations, while modelling the former with a bias field. We know a priori that intensity variations due to MR physics tend to be spatially smooth. A more accurate estimate of a bias field can be obtained by including prior knowledge about the distribution of the fields likely to be encountered by the correction algorithm. For example, if it is known that there is little or no intensity non-uniformity, then it would be wise to penalise large estimates of the intensity non-uniformity.'
                   'Knowing what works best should be a matter of empirical exploration, as it depends on the scans themselves.  For example, if your data has very little of the artifact, then the bias regularisation should be increased.  This effectively tells the algorithm that there is very little bias in your data, so it does not try to model it.'
                   }';
bparam.strtype = 'e';
bparam.num      = [1 1];
bparam.val      = {1e6};

realign         = cfg_exbranch;
realign.tag     = 'series';
realign.name    = 'Longitudinal Rigid Registration';
realign.val     = {data bparam};
realign.help    = {'Longitudinal registration of series of anatomical MRI scans for a single subject.  It is based on inverse-consistent alignment among each of the subject''s scans, and incorporates a bias field correction.  Prior to running the registration, the scans should already be in very rough alignment, although because the model incorporates a rigid-body transform, this need not be extremely precise.  Note that there are a bunch of hyper-parameters to be specified.  If you are unsure what values to take, then the defaults should be a reasonable guess of what works.  Note that changes to these hyper-parameters will impact the results obtained.'
''
'The alignment assumes that all scans have similar resolutions and dimensions, and were collected on the same (or very similar) MR scanner using the same pulse sequence.  If these assumption are not correct, then the approach will not work as well.'
''
'The resliced images are named the same as the originals, except that they are prefixed by ''r''.'};
realign.prog = @cat_vol_series_align;
realign.vout = @vout_reslice;


if expert
  %------------------------------------------------------------------------
  % Ultra-High Resolution Quantitative Image Optimization
  %------------------------------------------------------------------------

  % -- Data ---
  
  r1         = cfg_files; 
  r1.tag     = 'r1';
  r1.name    = 'R1-Volumes';
  r1.filter  = 'image';
  r1.ufilter = '.*';
  r1.num     = [1 Inf];
  r1.help    = {'Select R1 weighted images.'};

  pd         = cfg_files; 
  pd.tag     = 'pd';
  pd.name    = 'PD-Volumes';
  pd.filter  = 'image';
  pd.ufilter = '.*';
  pd.num     = [1 Inf];
  pd.help    = {'Select PD weighted images.'};

  r2s         = cfg_files; 
  r2s.tag     = 'r2s';
  r2s.name    = 'R2s-Volumes';
  r2s.filter  = 'image';
  r2s.ufilter = '.*';
  r2s.num     = [1 Inf];
  r2s.help    = {'Select R2s weighted images.'};

  data        = cfg_branch;
  data.tag    = 'data';
  data.name   = 'Input data';
  data.val    = {r1 pd r2s}; 
  data.help   = {
    'Input Images.'
  };

  
  % --- Parameter ---
  
  spm        = cfg_menu;
  spm.tag    = 'spm';
  spm.name   = 'Use SPM Preprocessing';
  spm.labels = {'No','Yes'};
  spm.values = {0 1};
  spm.val    = {1};
  spm.help   = {
    'Use SPM preprocessing if the data is not skull-stripped.'
  };
  
  bc        = cfg_menu;
  bc.tag    = 'bc';
  bc.name   = 'Bias Correction';
  bc.labels = {'No','light','medium','strong'};
  bc.values = {0 0.5 1 2};
  bc.val    = {1};
  bc.help   = {
    'Additional bias correction that is important for detection and correction of blood vessels.'
    ''
    'The correction uses a simple tissue classification and local filter approaches to estimate the local signal intensity in the WM and GM segment, e.g. a minimum/maximum filter in the WM for PD/T1 images.  Next, unclassified voxels were approximated and smoothed depending on the defined strength.  '
    ''
  };

  in        = cfg_menu;
  in.tag    = 'in';
  in.name   = 'Intensity Normalization';
  in.labels = {'No','Yes'};
  in.values = {0 1};
  in.val    = {1};
  in.help   = {
    'Additional global intensity normalization that is also important for detection and correction of blood vessels.'
    ''
  };

  bvc        = cfg_menu;
  bvc.tag    = 'bvc';
  bvc.name   = 'Blood Vessel Correction';
  bvc.labels = {'No','Yes'};
  bvc.values = {0 1};
  bvc.val    = {1};
  bvc.help   = {
    'Correction of blood vessels with high intensity in T1/R1/R2s and low intensity in PD images by CSF-like intensities. '
    ''
  };

  ss        = cfg_menu;
  ss.tag    = 'ss';
  ss.name   = 'Apply Skull-Stripping';
  ss.labels = {'No','Yes'};
  ss.values = {0 1};
  ss.val    = {1};
  ss.help   = {
    'Write skull-stripped images. '
    ''
  };

  nc        = cfg_menu;
  nc.tag    = 'nc';
  nc.name   = 'Noise Correction';
  nc.labels = {'No','Yes'};
  nc.values = {0 1};
  nc.val    = {1};
  nc.help   = {
    'Noise corrections of the final images.'
    ''
  };

  prefix         = cfg_entry;
  prefix.tag     = 'prefix';
  prefix.name    = 'Filename prefix';
  prefix.strtype = 's';
  prefix.num     = [0 Inf];
  prefix.val     = {'catsyn_'};
  prefix.help    = {
    'Prefix of output files.'};


  opts        = cfg_branch;
  opts.tag    = 'opts';
  opts.name   = 'Parameter';
  opts.val    = {spm bc in bvc ss nc prefix}; 
  opts.help   = {
    'Parameter settings for image correction.'
  };

  % --- Output
  
  pdo        = cfg_menu;
  pdo.tag    = 'pd';
  pdo.name   = 'PD Output';
  pdo.labels = {'No','Yes'};
  pdo.values = {0 1};
  pdo.val    = {1}; 
  pdo.help   = {
    'Write PD output images.'
  };

  t1o        = cfg_menu;
  t1o.tag    = 't1';
  t1o.name   = 'T1 Output';
  t1o.labels = {'No','Yes'};
  t1o.values = {0 1};
  t1o.val    = {1}; 
  t1o.help   = {
    'Write synthesized T1 output images based on the PD image.'
  };

  r1o        = cfg_menu;
  r1o.tag    = 'r1';
  r1o.name   = 'R1 Output';
  r1o.labels = {'No','Yes'};
  r1o.values = {0 1};
  r1o.val    = {1}; 
  r1o.help   = {
    'Write R1 output images.'
  };

  r2so        = cfg_menu;
  r2so.tag    = 'r2s';
  r2so.name   = 'R2s Output';
  r2so.labels = {'No','Yes'};
  r2so.values = {0 1};
  r2so.val    = {1}; 
  r2so.help   = {
    'Write R2s output images.'
  };

  bvco        = cfg_menu;
  bvco.tag    = 'bv';
  bvco.name   = 'Blood Vessel Output';
  bvco.labels = {'No','Yes'};
  bvco.values = {0 1};
  bvco.val    = {0}; 
  bvco.help   = {
    'Write map of blood vessels.'
  };
    
  output        = cfg_branch;
  output.tag    = 'output';
  output.name   = 'Output';
  output.val    = {r1o r2so pdo t1o bvco}; 
  output.help   = {
    'Output images.'
  };

  
  % ---
  
  urqio         = cfg_exbranch;
  urqio.tag     = 'urqio';
  urqio.name    = 'Ultra-High Resolution Quantitative Image Optimization';
  urqio.val     = {data opts output};
  urqio.prog    = @cat_vol_urqio;
  %urqio.vout    = @vfiles_urqio;
  urqio.help   = {
    'Additional correction of high resolution PD, R1, and R2s weighted images that includes another bias correction, intensity normalization, and blood vessel correction step. '
    ''
    'WARNING: This tool is in development and was just tested on a small set of subjects!'
  };
end

%------------------------------------------------------------------------
long          = cat_conf_long;
nonlin_coreg  = cat_conf_nonlin_coreg;
%------------------------------------------------------------------------

tools = cfg_choice;
tools.name   = 'Tools';
tools.tag    = 'tools';
tools.values = {showslice,check_cov,check_SPM,calcvol,calcroi,iqr,T2x,F2x,T2x_surf,F2x_surf,sanlm,realign,long,nonlin_coreg,defs,defs2}; %,qa
if expert 
  tools.values = [tools.values,{urqio}]; 
end
return

%_______________________________________________________________________

function vf = vfiles_defs(job)

PU = job.field1;
PI = job.images;

vf = cell(numel(PI),1);
for i=1:numel(PU),
    for m=1:numel(PI),
        [pth,nam,ext,num] = spm_fileparts(PI{m});
        
        switch job.modulate
        case 2
            fname = fullfile(pth,['m0w' nam ext num]);
        case 1
            fname = fullfile(pth,['mw' nam ext num]);
        case 0
            fname = fullfile(pth,['w' nam ext num]);
        end;
        vf{m} = fname;
    end
end

return;
%_______________________________________________________________________

function vf = vfiles_defs2(job)

PU = job.field;
PI = job.images;

vf = cell(numel(PU),numel(PI));
for i=1:numel(PU),
    for m=1:numel(PI),
        [pth,nam,ext,num] = spm_fileparts(PI{m}{i});
        
        switch job.modulate
        case 2
            fname = fullfile(pth,['m0w' nam ext num]);
        case 1
            fname = fullfile(pth,['mw' nam ext num]);
        case 0
            fname = fullfile(pth,['w' nam ext num]);
        end;
        vf{i,m} = fname;
    end
end

return;
%_______________________________________________________________________
function cdep = vfiles_urqio(job)
%%
cdep = cfg_dep;
if job.output.r1
  cdep(end+1)          = cfg_dep;
  cdep(end).sname      = 'R1 Images';
  cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'r1_'],'()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
if job.output.pd
  cdep(end+1)          = cfg_dep;
  cdep(end).sname      = 'PD Images';
  cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'pd_'],'()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
if job.output.t1
  cdep(end+1)          = cfg_dep;
  cdep(end).sname      = 'T1 Images';
  cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 't1_'],'()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
if job.output.r2s==1 || job.output.r2s==3
  cdep(end+1)           = cfg_dep;
  cdep(end).sname      = 'R2s nobc Images';
  cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'nobc_r2s_'],'()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end 
if job.output.r2s==2 || job.output.r2s==3
  cdep(end+1)          = cfg_dep;
  cdep(end).sname      = 'R2s bc Images';
  cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'bc_r2s_'],'()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end 
if job.output.bv
  cdep(end+1)          = cfg_dep;
  cdep(end).sname      = 'Blood Vessel Images';
  cdep(end).src_output = substruct('.','data','()',{1},'.',[job.opts.prefix 'bv_'],'()',{':'});
  cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end
if numel(cdep)>1
  cdep(1)=[];
end
%%
return;
%_______________________________________________________________________
function vf = vfiles_sanlm(job)
s  = cellstr(char(job.data)); vf = s; 
for i=1:numel(s),
    [pth,nam,ext,num] = spm_fileparts(s{i});
    vf{i} = fullfile(pth,[job.prefix,nam,ext,num]);
end;
return;
%_______________________________________________________________________
function vf = vfiles_qa(job)
s  = cellstr(char(job.data)); vf = s; 
for i=1:numel(s),
    [pth,nam,ext,num] = spm_fileparts(s{i});
    vf{i} = fullfile(pth,[job.prefix,nam,ext,num]);
end;
return;

%------------------------------------------------------------------------
function dep = vout_stat_TIV(job)

dep            = cfg_dep;
dep.sname      = 'TIV';
dep.src_output = substruct('.','calcvol','()',{':'});
dep.tgt_spec   = cfg_findspec({{'strtype','e','strtype','r'}});

%------------------------------------------------------------------------
function cdep = vout_reslice(job)

cdep(1)            = cfg_dep;
cdep(1).sname      = 'Midpoint Average';
cdep(1).src_output = substruct('.','avg','()',{':'});
cdep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
cdep(2)            = cfg_dep;
cdep(2).sname      = 'Realigned images';
cdep(2).src_output = substruct('.','rimg','()',{':'});
cdep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

return
 
%------------------------------------------------------------------------
function dep = vout_stat_spm2x(job)

dep            = cfg_dep;
dep.sname      = 'Transform & Threshold spm volumes';
dep.src_output = substruct('.','Pname');
dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});

%------------------------------------------------------------------------
function dep = vout_stat_spm2x_surf(job)

dep            = cfg_dep;
dep.sname      = 'Transform & Threshold spm surfaces';
dep.src_output = substruct('.','Pname');
dep.tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});

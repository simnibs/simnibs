function cat = tbx_cfg_cat
% Configuration file for segment jobs
%_______________________________________________________________________
%
% Christian Gaser
% $Id: tbx_cfg_cat.m 1245 2017-12-12 10:35:00Z dahnke $
%
%#ok<*AGROW>
 
addpath(fileparts(which(mfilename)));

%% ------------------------------------------------------------------------
try
  expert = cat_get_defaults('extopts.expertgui');
catch %#ok<CTCH>
  expert = 0; 
end
if isempty(expert) 
  expert = 0;
end  

% try to estimate number of processor cores
try
  numcores = feature('numcores');
  % because of poor memory management use only half of the cores for windows
  if ispc
    numcores = round(numcores/2);
  end
  numcores = max(numcores,1);
catch
  numcores = 0;
end

% force running in the foreground if only one processor was found or for compiled version
if numcores == 1 || isdeployed, numcores = 0; end

%_______________________________________________________________________
nproc         = cfg_entry;
nproc.tag     = 'nproc';
nproc.name    = 'Split job into separate processes';
nproc.strtype = 'w';
nproc.val     = {numcores};
nproc.num     = [1 1];
nproc.help    = {
    'In order to use multi-threading the CAT12 segmentation job with multiple subjects can be split into separate processes that run in the background. You can even close Matlab, which will not affect the processes that will run in the background without GUI. If you do not want to run processes in the background then set this value to 0.'
    ''
    'Keep in mind that each process needs about 1.5..2GB of RAM, which should be considered to choose the appropriate  number of processes.'
    ''
    'Please further note that no additional modules in the batch can be run except CAT12 segmentation. Any dependencies will be broken for subsequent modules.'
  };

%_______________________________________________________________________

data          = cfg_files;
data.tag      = 'data';
data.name     = 'Volumes';
data.filter   = 'image';
data.ufilter  = '.*';
data.num      = [0 Inf];
data.help     = {
  'Select highres raw data (e.g. T1 images) for segmentation. This assumes that there is one scan for each subject. Note that multi-spectral (when there are two or more registered images of different contrasts) processing is not implemented for this method.'};
data.preview  = @(f) spm_check_registration(char(f));

data_spm          = cfg_files;
data_spm.tag      = 'data';
data_spm.name     = 'Segmentations in native space';
data_spm.filter   = 'image';
data_spm.ufilter  = '^c1.*';
data_spm.num      = [0 Inf];
data_spm.help     = {
  'Select SPM segmentations for class 1 for all subjects. Names for all other remaining classes 2 and 3 are automatically estimated.'};
data_spm.preview  = @(f) spm_check_registration(char(f));


%------------------------------------------------------------------------
% writing options for data
%------------------------------------------------------------------------

%------------------------------------------------------------------------

  surface        = cfg_menu;
  surface.tag    = 'surface';
  surface.name   = 'Surface and thickness estimation';
  surface.def    = @(val)cat_get_defaults('output.surface', val{:});
  surface.help   = {
    'Use projection-based thickness (PBT) (Dahnke et al. 2012) to estimate cortical thickness and to create the central cortical surface for left and right hemisphere. Surface reconstruction includes topology correction (Yotter et al. 2011), spherical inflation (Yotter et al.) and spherical registration.'
    ''
    'Please note, that surface reconstruction additionally requires about 20-60 min of computation time.'
    ''
  };

if expert==0
  surface.labels = {'No','Yes'};
  surface.values = {0 1};
elseif expert==1
  surface.labels = {'No','lh + rh','lh + rh + cerebellum','lh + rh (fast, no registration, only for quick quality check and not for analysis)','Full'};
  surface.values = {0 1 2 5 12};  
  surface.help   = [surface.help; {
    'Cerebellar reconstruction is still in development and is strongly limited due to the high frequency of folding and image properties! '
    ''
    'The fast reconstruction allows a _visual_ check of the processing quality by providing a reconstructed cortical surface and thickness values with slightly lower accuracy. It takes only a few minutes for both hemisheres by internally using 0.8 mm (instead of 0.5mm) spatial resolution and by skipping topology correction and surface registration. Please note that neither the surface nor the measures on the surface (e.g. thickness) can be used for any further statistical analysis! '
    ''
  }];    
else
  surface.labels = {'No','lh + rh','lh + rh + cerebellum',...
    'lh + rh (fast, no registration)','lh + rh + cerebellum (fast, no registration)', ...
    'lh + rh (fast registration)','lh + rh + cerebellum (fast registration)'};
  surface.values = {0 1 2 5 6 7 8};
end


if expert==0
  %------------------------------------------------------------------------
  % Parameter to choose between different ways to extract data by ROIs.
  % Inactive due to missing implementation and evaluation of the optimized 
  % space ROI analysis. RD 20160112
  % Use this now only for default user mode.
  %------------------------------------------------------------------------
  ROI        = cfg_menu;
  ROI.tag    = 'ROI';
  ROI.name   = 'ROI analyis';
  ROI.labels = {'No ROI analyis','ROI analysis'};
  ROI.values = {0 1};
  ROI.def    = @(val)cat_get_defaults('output.ROI', val{:});
  ROI.help   = {
    'Export of ROI data of volume, intensity, and thickness to a xml-files. '
    'For further information see atlas specific text files in "templates_1.50mm" CAT12 subdir. '
    ''
    'For thickness estimation the projection-based thickness (PBT) [Dahnke:2012] is used that average cortical thickness for each GM voxel. '
    ''
    'There are different atlas maps available: '
    '(1) Hammers (68 CSF/GM/[WM] ROIs of 20 subjects, 2003):'
    '    Alexander Hammers brain atlas from the Euripides project (www.brain-development.org).'
    '    Hammers et al. Three-dimensional maximum probability atlas of the human brain, with particular reference to the temporal lobe. Hum Brain Mapp 2003, 19: 224-247.'
    ''
    '(2) Neuromorphometrics (142 GM ROIs of 15 subjects, 2012):'
    '    Maximum probability tissue labels derived from the MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling'
    '    https://masi.vuse.vanderbilt.edu/workshop2012/index.php/Challenge_Details'
    ''
    '(3) LPBA40 (56 GM ROIs of 40 subjects, 2008):'
    '    The LONI Probabilistic Brain Atlas (LPBA40) is a series of maps of brain anatomic regions. These maps were produced from a set of whole-head MRI of 40 human volunteers. Each MRI was manually delineated to identify a set of 56 structures in the brain, most of which are within the cortex. These delineations were then transformed into a common atlas space to produce a set of coregistered anatomical labels. The original MRI data were also transformed into the atlas space. '
    '    Shattuck et al. 2008. Construction of a 3D Probabilistic Atlas of Human Cortical Structures, NeuroImage 39 (3): 1064-1070. DOI:	10.1016/j.neuroimage.2007.09.031'
    ''
    ...'(4) IBSR (32 CSF/GM ROIs of 18 subjects, 2004):'
    ...'    See IBSR terms "http://www.nitrc.org/projects/ibsr"'
    ...''
    ...'(5) AAL (122 GM ROIs of 1 subject, 2002):'
    ...'    Combining probabilistic cytoarchitectonic maps and functional imaging data of a single Brain.'
    ...'    Tzourio-Mazoyer et al., Automated anatomical labelling of activations in spm using a macroscopic anatomical parcellation of the MNI MRI single subject brain. Neuroimage 2002, 15: 273-289.'
    ...'    http://www.fz-juelich.de/inm/inm-1/DE/Forschung/_docs/SPMAnatomyToolbox/SPMAnatomyToolbox_node.html'
    ...''
    ...'(6) MORI (128 GM/WM ROIs of 1 subject, 2009):'
    ...'    Oishi et al. Atlas-based whole brain white matter analysis using large deformation diffeomorphic metric mapping: application to normal elderly and Alzheimer''s disease participants. 2009'
    ... ''
    ... '() Anatomy ():
    ... '() COBRA ():'
    ''
    ''
  };

else
  %% ------------------------------------------------------------------------
  % Parameter to choose between different ways to extract data by ROIs.
  % Inactive due to missing implementation and evaluation of the optimized 
  % space ROI analysis. 
  % RD 20160112
  %------------------------------------------------------------------------

  noROI        = cfg_branch;
  noROI.tag    = 'noROI';
  noROI.name   = 'No ROI processing';
  noROI.help   = {'No ROI processing'};

  exatlas  = cat_get_defaults('extopts.atlas'); 
  matlas = {}; mai = 1; 
  for ai = 1:size(exatlas,1)
    if exatlas{ai,2}<=expert && exist(exatlas{ai,1},'file')
      [pp,ff,ee]  = spm_fileparts(exatlas{ai,1}); 

      % if output.atlases.ff does not exist then set it by the default file value
      try
        cat_get_defaults(['output.atlases.' ff]);
      catch
        cat_get_defaults(['output.atlases.' ff], exatlas{ai,4})
      end

      matlas{mai}        = cfg_menu;
      matlas{mai}.tag    = ff;
      matlas{mai}.name   = ff; 
      matlas{mai}.labels = {'No','Yes'};
      matlas{mai}.values = {0 1};
      matlas{mai}.def    = @(val)cat_get_defaults(['output.atlases.' ff], val{:});
      txtfile = fullfile(pp,[ff '.txt']);
      if exist(txtfile,'file')
        fid = fopen(txtfile,'r');
        txt = textscan(fid,'%s','delimiter','\n');
        fclose(fid);
        matlas{mai}.help   = [{ 
          'Processing flag of this atlas map.'
          ''
          }
          txt{1}];
      else
        matlas{mai}.help   = {
          'Processing flag of this atlas map.'
          ''
          ['No atlas readme text file "' txtfile '"!']
        };
      end
      mai = mai+1; 
    end
  end
  
  atlases          = cfg_branch;
  atlases.tag      = 'atlases';
  atlases.name     = 'Atlases';
  atlases.val      = matlas;
  atlases.help     = {'Writing options of ROI atlas maps.'
  ''
  };


  ROI        = cfg_choice;
  ROI.tag    = 'ROImenu';
  ROI.name   = 'ROI analysis';
  if cat_get_defaults('output.ROI')>0
    ROI.val  = {atlases};
  else
    ROI.val  = {noROI};
  end
  ROI.values = {noROI atlases};
  ROI.help   = {
    'Export of ROI data of volume, intensity, and thickness to a xml-files. '
    'For further information see atlas specific text files in "templates_1.50mm" CAT12 subdir. '
    ''
    'For thickness estimation the projection-based thickness (PBT) [Dahnke:2012] is used that average cortical thickness for each GM voxel. '
    ''
    'There are different atlas maps available: '
    '(1) Hammers (68 CSF/GM/[WM] ROIs of 20 subjects, 2003):'
    '    Alexander Hammers brain atlas from the Euripides project (www.brain-development.org).'
    '    Hammers et al. Three-dimensional maximum probability atlas of the human brain, with particular reference to the temporal lobe. Hum Brain Mapp 2003, 19: 224-247.'
    ''
    '(2) Neuromorphometrics (142 GM ROIs of 15 subjects, 2012):'
    '    Maximum probability tissue labels derived from the MICCAI 2012 Grand Challenge and Workshop on Multi-Atlas Labeling'
    '    https://masi.vuse.vanderbilt.edu/workshop2012/index.php/Challenge_Details'
    ''
    '(3) LPBA40 (56 GM ROIs of 40 subjects, 2008):'
    '    The LONI Probabilistic Brain Atlas (LPBA40) is a series of maps of brain anatomic regions. These maps were produced from a set of whole-head MRI of 40 human volunteers. Each MRI was manually delineated to identify a set of 56 structures in the brain, most of which are within the cortex. These delineations were then transformed into a common atlas space to produce a set of coregistered anatomical labels. The original MRI data were also transformed into the atlas space. '
    '    Shattuck et al. 2008. Construction of a 3D Probabilistic Atlas of Human Cortical Structures, NeuroImage 39 (3): 1064-1070. DOI:	10.1016/j.neuroimage.2007.09.031'
    ''
    '(4) COBRA (1 GM/WM ROI in amgdala, 2 combined GM/WM ROIs in hippocampus and 13 GM/WM ROIs in cerebellum of 5 subjects):'
    '    The Cobra atlas is build from 3 atlases that are provided by the Computational Brain Anatomy Laboratory at the Douglas Institute (CoBra Lab). The 3 atlases are based on high-resolution (0.3mm isotropic voxel size) images of the amygdala, hippocampus and the cerebellum. Some of the hippocampus subfields were merged because of their small size (CA1/CA2/CA3/stratum radiatum/subiculum/stratum lacunosum/stratum moleculare). Please note that the original labels were changed in order to allow a combined atlas. '
    '    Entis JJ, Doerga P, Barrett LF, Dickerson BC. A reliable protocol for the manual segmentation of the human amygdala and its subregions using ultra-high resolution MRI. Neuroimage. 2012;60(2):1226-35.'
    '    Winterburn JL, Pruessner JC, Chavez S, et al. A novel in vivo atlas of human hippocampal subfields using high-resolution 3 T magnetic resonance imaging.  Neuroimage. 2013;74:254-65.'
    '    Park, M.T., Pipitone, J., Baer, L., Winterburn, J.L., Shah, Y., Chavez, S., Schira, M.M., Lobaugh, N.J., Lerch, J.P., Voineskos, A.N., Chakravarty, M.M. Derivation of high-resolution MRI atlases of the human cerebellum at 3T and segmentation using multiple automatically generated templates. Neurimage. 2014; 95: 217-31.'
    ''
    '(5) IBSR (32 CSF/GM ROIs of 18 subjects, 2004):'
    '    See IBSR terms "http://www.nitrc.org/projects/ibsr"'
    ''
    '(6) AAL (122 GM ROIs of 1 subject, 2002):'
    '    Combining probabilistic cytoarchitectonic maps and functional imaging data of a single Brain.'
    '    Tzourio-Mazoyer et al., Automated anatomical labelling of activations in spm using a macroscopic anatomical parcellation of the MNI MRI single subject brain. Neuroimage 2002, 15: 273-289.'
    '    http://www.fz-juelich.de/inm/inm-1/DE/Forschung/_docs/SPMAnatomyToolbox/SPMAnatomyToolbox_node.html'
    ''
    '(7) MORI (128 GM/WM ROIs of 1 subject, 2009):'
    '    Oishi et al. Atlas-based whole brain white matter analysis using large deformation diffeomorphic metric mapping: application to normal elderly and Alzheimer''s disease participants. 2009'
     ''
    '(8) Anatomy (44 GM/WM ROIs in 10 post-mortem subjects, 2014):'
    '    Eickhoff SB, Stephan KE, Mohlberg H, Grefkes C, Fink GR, Amunts K, Zilles K. A new SPM toolbox for combining probabilistic cytoarchitectonic maps and functional imaging data. NeuroImage 25(4), 1325-1335, 2005'
    ''
  }; 

end
  
 
native        = cfg_menu;
native.tag    = 'native';
native.name   = 'Native space';
native.labels = {'No','Yes'};
native.values = {0 1};
native.help   = {
  'The native space option allows you to produce a tissue class image (p*) that is in alignment with the original/* (see Figure \ref{seg1})*/. It can also be used for ''''importing'''' into a form that can be used with the DARTEL toolbox (rp*).'
''
};

warped        = cfg_menu;
warped.tag    = 'warped';
warped.name   = 'Normalized';
warped.labels = {'No','Yes'};
warped.values = {0 1};
warped.help   = {'Write image in normalized space without any modulation.'
''
};

dartel        = cfg_menu;
dartel.tag    = 'dartel';
dartel.name   = 'DARTEL export';
if expert
  dartel.labels = {'No','Rigid (SPM12 default)','Affine','Both'};
  dartel.values = {0 1 2 3};
else
  dartel.labels = {'No','Rigid (SPM12 default)','Affine'};
  dartel.values = {0 1 2};
end
dartel.help   = {
'This option is to export data into a form that can be used with DARTEL. The SPM default is to only apply rigid body transformation. However, a more appropriate option is to apply affine transformation, because the additional scaling of the images requires less deformations to non-linearly register brains to the template.'
''
};

native.def  = @(val)cat_get_defaults('output.bias.native', val{:});
warped.def  = @(val)cat_get_defaults('output.bias.warped', val{:});
dartel.def  = @(val)cat_get_defaults('output.bias.dartel', val{:});
bias        = cfg_branch;
bias.tag    = 'bias';
bias.name   = 'Bias, noise and global intensity corrected T1 image';
if expert
  bias.val    = {native warped dartel};
else
  bias.val    = {warped};
end
bias.help   = {
  'This is the option to save a bias, noise, and global intensity corrected version of the original T1 image. MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images. The bias corrected version should have more uniform intensities within the different types of tissues and can be saved in native space and/or normalised. Noise is corrected by an adaptive non-local mean (NLM) filter (Manjon 2008, Medical Image Analysis 12).'
''
};

native.def  = @(val)cat_get_defaults('output.las.native', val{:});
warped.def  = @(val)cat_get_defaults('output.las.warped', val{:});
dartel.def  = @(val)cat_get_defaults('output.las.dartel', val{:});
las        = cfg_branch;
las.tag    = 'las';
las.name   = 'Bias, noise and local intensity corrected T1 image';
las.val    = {native warped dartel};
las.help   = {
  'This is the option to save a bias, noise, and local intensity corrected version of the original T1 image. MR images are usually corrupted by a smooth, spatially varying artifact that modulates the intensity of the image (bias). These artifacts, although not usually a problem for visual inspection, can impede automated processing of the images. The bias corrected version should have more uniform intensities within the different types of tissues and can be saved in native space and/or normalised. Noise is corrected by an adaptive non-local mean (NLM) filter (Manjon 2008, Medical Image Analysis 12).'
''
};


%------------------------------------------------------------------------

warped.def    = @(val)cat_get_defaults('output.jacobian.warped', val{:});
jacobian      = cfg_branch;
jacobian.tag  = 'jacobian';
jacobian.name = 'Jacobian determinant';
jacobian.val  = {warped};
jacobian.help = {
  'This is the option to save the Jacobian determinant, which expresses local volume changes. This image can be used in a pure deformation based morphometry (DBM) design. Please note that the affine part of the deformation field is ignored. Thus, there is no need for any additional correction for different brain sizes using ICV.'
''
};

%------------------------------------------------------------------------

native.def  = @(val)cat_get_defaults('output.label.native', val{:});
warped.def  = @(val)cat_get_defaults('output.label.warped', val{:});
dartel.def  = @(val)cat_get_defaults('output.label.dartel', val{:});

label       = cfg_branch;
label.tag   = 'label';
label.name  = 'PVE label image';
label.val   = {native warped dartel};
label.help  = {
'This is the option to save a labeled version of your segmentations. Labels are saved as Partial Volume Estimation (PVE) values with different mix classes for GM-WM (2.5) and GM-CSF (0.5). BG=0, CSF=1, GM=2, WM=3, WMH=4 (if WMHC=3)'
''
};

%------------------------------------------------------------------------

modulated        = cfg_menu;
modulated.tag    = 'mod';
modulated.name   = 'Modulated normalized';
if expert
  modulated.labels = {'No','Affine + non-linear (SPM12 default)','Non-linear only'};
  modulated.values = {0 1 2};
else
  modulated.labels = {'No','Yes'};
  modulated.values = {0 1};
end
modulated.help = {
'"Modulation" is to compensate for the effect of spatial normalisation. Spatial normalisation causes volume changes due to affine transformation (global scaling) and non-linear warping (local volume change). After modulation the resulting modulated images are preserved for the total amount of grey matter signal in the normalised partitions. Thus, modulated images reflect the tissue volumes before spatial normalisation. However, the user is almost always interested in removing the confound of different brain sizes and there are many ways to apply this correction. In contrast to previous VBM versions I now recommend to use total intracranial volume (TIV) as nuisance parameter in an AnCova model. '
''
'Please note that I do not use the SPM modulation where the original voxels are projected into their new location in the warped images because this method introduces aliasing artifacts. Here, I use the scaling by the Jacobian determinants to generate "modulated" data. '
''
};

native.def    = @(val)cat_get_defaults('output.GM.native', val{:});
warped.def    = @(val)cat_get_defaults('output.GM.warped', val{:});
modulated.def = @(val)cat_get_defaults('output.GM.mod', val{:});
dartel.def    = @(val)cat_get_defaults('output.GM.dartel', val{:});
grey          = cfg_branch;
grey.tag      = 'GM';
grey.name     = 'Grey matter';
grey.help     = {'Options to produce grey matter images.'
''
};
grey_spm      = grey;
if expert
  grey.val      = {native warped modulated dartel};
  grey_spm.val  = {warped modulated dartel};
else
  grey.val      = {native modulated dartel};
  grey_spm.val  = {modulated dartel};
end

native.def    = @(val)cat_get_defaults('output.WM.native', val{:});
warped.def    = @(val)cat_get_defaults('output.WM.warped', val{:});
modulated.def = @(val)cat_get_defaults('output.WM.mod',    val{:});
dartel.def    = @(val)cat_get_defaults('output.WM.dartel', val{:});
white         = cfg_branch;
white.tag     = 'WM';
white.name    = 'White matter';
white.help    = {'Options to produce white matter images.'
''
};
white_spm     = white;
if expert
  white.val      = {native warped modulated dartel};
  white_spm.val  = {warped modulated dartel};
else
  white.val      = {native modulated dartel};
  white_spm.val  = {modulated dartel};
end


native.def    = @(val)cat_get_defaults('output.CSF.native', val{:});
warped.def    = @(val)cat_get_defaults('output.CSF.warped', val{:});
modulated.def = @(val)cat_get_defaults('output.CSF.mod',    val{:});
dartel.def    = @(val)cat_get_defaults('output.CSF.dartel', val{:});
csf           = cfg_branch;
csf.tag       = 'CSF';
csf.name      = 'Cerebro-Spinal Fluid (CSF)';
csf.help      = {'Options to produce CSF images.'
''
};
csf_spm       = csf;
csf.val       = {native warped modulated dartel};
csf_spm.val   = {warped modulated dartel};

% head/background tissue classes
native.def    = @(val)cat_get_defaults('output.TPMC.native', val{:});
warped.def    = @(val)cat_get_defaults('output.TPMC.warped', val{:});
modulated.def = @(val)cat_get_defaults('output.TPMC.mod',    val{:});
dartel.def    = @(val)cat_get_defaults('output.TPMC.dartel', val{:});
tpmc          = cfg_branch;
tpmc.tag      = 'TPMC';
tpmc.name     = 'Tissue Probability Map Classes';
tpmc.val      = {native warped modulated dartel};
tpmc.help     = {'Option to produce the SPM tissue class 4 to 6: p#*.img, wp#*.img and m[0]wp#*.img.'
''
};

% WMH
native.def    = @(val)cat_get_defaults('output.WMH.native', val{:});
warped.def    = @(val)cat_get_defaults('output.WMH.warped', val{:});
modulated.def = @(val)cat_get_defaults('output.WMH.mod',    val{:});
dartel.def    = @(val)cat_get_defaults('output.WMH.dartel', val{:});
wmh           = cfg_branch;
wmh.tag       = 'WMH';
wmh.name      = 'White matter hyperintensities (WMH)';
wmh.val       = {native warped modulated dartel};
wmh.help      = {'Options to produce WMH images, if WMHC==3: p7*.img, wp7*.img and m[0]wp7*.img.'
''
};

% main structure atlas
native.def    = @(val)cat_get_defaults('output.atlas.native', val{:});
warped.def    = @(val)cat_get_defaults('output.atlas.warped', val{:});
dartel.def    = @(val)cat_get_defaults('output.atlas.dartel', val{:});
atlas         = cfg_branch;
atlas.tag     = 'atlas';
atlas.name    = 'Atlas label maps';
atlas.val     = {native warped dartel};
atlas.help    = {
  'WARNING: The functions that create this maps are still under development! This is the option to save an atlas map with major structures (a1*). Odd numbers code the left, even numbers the right hemisphere. Furthermore, AAL and Broadman atlas maps were created based on maps from MRIcron that where adapted to the other VBM maps. Other maps are used from the IBASPM toolbox.  http://www.thomaskoenig.ch/Lester/ibaspm.htmAnatomy toolbox:Alexander Hammers brain atlas from the Euripides project:   www.brain-development.org  Hammers A, Allom R, Koepp MJ, Free SL, Myers R, Lemieux L, Mitchell   TN, Brooks DJ, Duncan JS. Three-dimensional maximum probability atlas   of the human brain, with particular reference to the temporal lobe.   Hum Brain Mapp 2003, 19: 224-247.'
''
};


warps = cfg_menu;
warps.tag    = 'warps';
warps.name   = 'Deformation Fields';
warps.labels = {
    'No'
    'Image->Template (forward)'
    'Template->Image (inverse)'
    'inverse + forward'};
warps.values = {[0 0],[1 0],[0 1],[1 1]};
warps.def    = @(val)cat_get_defaults('output.warps', val{:});
warps.help   = {
  'Deformation fields can be saved to disk, and used by the Deformations Utility and/or applied to coregistered data from other modalities (e.g. fMRI). For spatially normalising images to MNI space, you will need the forward deformation, whereas for spatially normalising (eg) GIFTI surface files, you''ll need the inverse. It is also possible to transform data in MNI space on to the individual subject, which also requires the inverse transform. Deformations are saved as .nii files, which contain three volumes to encode the x, y and z coordinates.'
''
};

%------------------------------------------------------------------------

output      = cfg_branch;
output.tag  = 'output';
output.name = 'Writing options';
if expert==2
  output.val  = {surface ROI grey white csf wmh tpmc atlas label bias las jacobian warps}; 
elseif expert==1
  output.val  = {surface ROI grey white csf wmh label bias las jacobian warps};
else
  output.val  = {surface ROI grey white bias jacobian warps};
end
output.help = {
'There are a number of options about what data you would like the routine to produce. The routine can be used for producing images of tissue classes, as well as bias corrected images. The native space option will produce a tissue class image (p*) that is in alignment with the original image. You can also produce spatially normalised versions - both with (m[0]wp*) and without (wp*) modulation. In the cat toolbox, the voxel size of the spatially normalised versions is 1.5 x 1.5 x 1.5mm as default. The produced images of the tissue classes can directly be used for doing voxel-based morphometry (both un-modulated and modulated). All you need to do is smooth them and do the stats (which means no more questions on the mailing list about how to do "optimized VBM").'
''
'Modulation is to compensate for the effect of spatial normalisation. When warping a series of images to match a template, it is inevitable that volumetric differences will be introduced into the warped images. For example, if one subject''s temporal lobe has half the volume of that of the template, then its volume will be doubled during spatial normalisation. This will also result in a doubling of the voxels labeled grey matter. In order to remove this confound, the spatially normalised grey matter (or other tissue class) is adjusted by multiplying by its relative volume before and after warping. If warping results in a region doubling its volume, then the correction will halve the intensity of the tissue label. This whole procedure has the effect of preserving the total amount of grey matter signal in the normalised partitions.'
''
'A deformation field is a vector field, where three values are associated with each location in the field. The field maps from co-ordinates in the normalised image back to co-ordinates in the original image. The value of the field at co-ordinate [x y z] in the normalised space will be the co-ordinate [x'' y'' z''] in the original volume. The gradient of the deformation field at a co-ordinate is its Jacobian matrix, and it consists of a 3x3 matrix:'
'%   /                      \%   | dx''/dx  dx''/dy dx''/dz |%   |                       |%   | dy''/dx  dy''/dy dy''/dz |%   |                       |%   | dz''/dx  dz''/dy dz''/dz |%   \                      /'
''
'The value of dx''/dy is a measure of how much x'' changes if y is changed by a tiny amount. The determinant of the Jacobian is the measure of relative volumes of warped and unwarped structures.  The modulation step simply involves multiplying by the relative volumes.'};


%% ------------------------------------------------------------------------
tools      = cat_conf_tools(expert);     % volume tools
stools     = cat_conf_stools(expert);    % surface tools
if expert 
  stoolsexp  = cat_conf_stoolsexp;       % surface expert tools
end
extopts    = cat_conf_extopts(expert);   
opts       = cat_conf_opts(expert); 
%------------------------------------------------------------------------

estwrite        = cfg_exbranch;
estwrite.tag    = 'estwrite';
estwrite.name   = 'CAT12: Segmentation';
% use multithreading only if availabe
if feature('numcores') > 1
  estwrite.val    = {data nproc opts extopts output};
else
  estwrite.val    = {data opts extopts output};
end
estwrite.prog   = @cat_run;
estwrite.vout   = @vout;
estwrite.help   = {
'This toolbox is an extension to the default segmentation in SPM12, but uses a completely different segmentation approach.'
''
'The segmentation approach is based on an Adaptive Maximum A Posterior (MAP) technique without the need for a priori information about tissue probabilities. That is, the Tissue Probability Maps (TPM) are not used constantly in the sense of the classical Unified Segmentation approach (Ashburner et. al. 2005), but just for spatial normalization. The following AMAP estimation is adaptive in the sense that local variations of the parameters (i.e., means and variance) are modeled as slowly varying spatial functions (Rajapakse et al. 1997). This not only accounts for intensity inhomogeneities but also for other local variations of intensity.'
''
'Additionally, the segmentation approach uses a Partial Volume Estimation (PVE) with a simplified mixed model of at most two tissue types (Tohka et al. 2004). We start with an initial segmentation into three pure classes: gray matter (GM), white matter (WM), and cerebrospinal fluid (CSF) based on the above described AMAP estimation. The initial segmentation is followed by a PVE of two additional mixed classes: GM-WM and GM-CSF. This results in an estimation of the amount (or fraction) of each pure tissue type present in every voxel (as single voxels - given by their size - probably contain more than one tissue type) and thus provides a more accurate segmentation.'
''
'Another important extension to the SPM12 segmentation is the integration of the Dartel normalisation (Ashburner 2007) into the toolbox by an already existing Dartel template in MNI space. This template was derived from 555 healthy control subjects of the IXI-database (http://www.brain-development.org) and provides the six Dartel iteration. Thus, for the majority of studies the creation of sample-specific Dartel templates is not necessary anymore.'};

%------------------------------------------------------------------------
% CAT surface processing with existing SPM segmentation 

extopts_spm = cat_conf_extopts(expert,1);   
output_spm  = output; 
if expert==2
  output_spm.val  = {ROI surface grey_spm white_spm csf_spm label jacobian warps}; 
elseif expert==1
  output_spm.val  = {ROI surface grey_spm white_spm csf_spm jacobian warps};
else % also CSF output because it is requiered as input ...
  output_spm.val  = {ROI surface grey_spm white_spm csf_spm jacobian warps};
end

estwrite_spm        =  cfg_exbranch;
estwrite_spm.tag    = 'estwrite_spm';
estwrite_spm.name   = 'CAT12: SPM Segmentation';
% use multithreading only if availabe
if feature('numcores') > 1
  estwrite_spm.val  = {data_spm nproc extopts_spm output_spm};
else
  estwrite_spm.val  = {data_spm extopts_spm output_spm};
end
estwrite_spm.prog   = @cat_run;
estwrite_spm.vout   = @vout;
estwrite_spm.help   = {
'CAT processing with thickness estimation and surface creation for SPM segmentation which is using the input of CSF, GM, and WM and also integrates Dartel normalisation (Ashburner 2007) into the toolbox by an already existing Dartel template in MNI space. This template was derived from 555 healthy control subjects of the IXI-database (http://www.brain-development.org) and provides the six Dartel iteration. Thus, for the majority of studies the creation of sample-specific Dartel templates is not necessary anymore.'};


%------------------------------------------------------------------------
cat        = cfg_choice;
cat.name   = 'CAT12';
cat.tag    = 'cat';
if expert==2
  cat.values = {estwrite estwrite_spm tools stools stoolsexp};
elseif expert==1
  cat.values = {estwrite estwrite_spm tools stools};
else
  cat.values = {estwrite tools stools}; 
end
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function dep = vout(job)

opts  = job.output;

if isfield(opts.GM,'warped') && isfield(opts.GM,'native')
  tissue(1).warped = [opts.GM.warped  (opts.GM.mod==1)        (opts.GM.mod==2)       ];
  tissue(1).native = [opts.GM.native  (opts.GM.dartel==1)     (opts.GM.dartel==2)    ];
  tissue(2).warped = [opts.WM.warped  (opts.WM.mod==1)        (opts.WM.mod==2)       ];
  tissue(2).native = [opts.WM.native  (opts.WM.dartel==1)     (opts.WM.dartel==2)    ];
elseif ~isfield(opts.GM,'native')
  if isfield(opts.GM,'warped')
    tissue(1).warped = [opts.GM.warped  (opts.GM.mod==1)        (opts.GM.mod==2)       ];
    tissue(2).warped = [opts.WM.warped  (opts.WM.mod==1)        (opts.WM.mod==2)       ];
  else
    tissue(1).warped = [0               (opts.GM.mod==1)        (opts.GM.mod==2)       ];
    tissue(2).warped = [0               (opts.WM.mod==1)        (opts.WM.mod==2)       ];
   end
else
  tissue(1).warped = [0               (opts.GM.mod==1)        (opts.GM.mod==2)       ];
  tissue(1).native = [opts.GM.native  (opts.GM.dartel==1)     (opts.GM.dartel==2)    ];
  tissue(2).warped = [0               (opts.WM.mod==1)        (opts.WM.mod==2)       ];
  tissue(2).native = [opts.WM.native  (opts.WM.dartel==1)     (opts.WM.dartel==2)    ];
end

if isfield(opts,'CSF')
  tissue(3).warped = [opts.CSF.warped (opts.CSF.mod==1)       (opts.CSF.mod==2)      ];
  if isfield(opts.CSF,'native')
    tissue(3).native = [opts.CSF.native (opts.CSF.dartel==1)    (opts.CSF.dartel==2)   ];
  end
end

% This depends on job contents, which may not be present when virtual
% outputs are calculated.

% CAT report XML file
cdep = cfg_dep;
cdep(end).sname      = 'CAT Report';
cdep(end).src_output = substruct('.','catreport','()',{':'});
cdep(end).tgt_spec   = cfg_findspec({{'filter','xml','strtype','e'}});

% lh/rh central surface and thickness
if isfield(opts,'surface')
  if opts.surface,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Left Central Surface';
    cdep(end).src_output = substruct('()',{1}, '.','lhcentral','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Right Central Surface';
    cdep(end).src_output = substruct('()',{1}, '.','rhcentral','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Left Thickness';
    cdep(end).src_output = substruct('()',{1}, '.','lhthickness','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Right Thickness';
    cdep(end).src_output = substruct('()',{1}, '.','rhthickness','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
  end
end;

% XML label
if isfield(opts,'ROI')
  if opts.ROI,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'ROI XML File';
    cdep(end).src_output = substruct('()',{1}, '.','roi','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','xml','strtype','e'}});
  end
end;

% bias corrected
if isfield(opts,'bias')
  if isfield(opts.bias,'native')
    if opts.bias.native,
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Native Bias Corr. Image';
      cdep(end).src_output = substruct('()',{1}, '.','biascorr','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
  end;
  if opts.bias.warped,
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Warped Bias Corr. Image';
      cdep(end).src_output = substruct('()',{1}, '.','wbiascorr','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
end

% LAS bias corrected
if isfield(opts,'las')
  if opts.las.native,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Native LAS Bias Corr. Image';
    cdep(end).src_output = substruct('()',{1}, '.','ibiascorr','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end
  if opts.las.warped,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Warped LAS Bias Corr. Image';
    cdep(end).src_output = substruct('()',{1}, '.','wibiascorr','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
  if opts.las.dartel==1,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Rigidly Registered LAS Bias Corr. Image';
    cdep(end).src_output = substruct('()',{1}, '.','ribiascorr','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
  if opts.las.dartel==2,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Affine Registered LAS Bias Corr. Image';
    cdep(end).src_output = substruct('()',{1}, '.','aibiascorr','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
end;

% label
if isfield(opts,'label')
  if opts.label.native,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Native Label Image';
    cdep(end).src_output = substruct('()',{1}, '.','label','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
  if opts.label.warped,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Warped Label Image';
    cdep(end).src_output = substruct('()',{1}, '.','wlabel','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
  if opts.label.dartel==1,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Rigidly Registered Label Image';
    cdep(end).src_output = substruct('()',{1}, '.','rlabel','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
  if opts.label.dartel==2,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Affine Registered Label Image';
    cdep(end).src_output = substruct('()',{1}, '.','alabel','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
end

% atlas
if isfield(opts,'atlas')
  if opts.atlas.native,
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Native Atlas Image';
      cdep(end).src_output = substruct('()',{1}, '.','atlas','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
  if opts.atlas.warped,
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Warped Atlas Image';
      cdep(end).src_output = substruct('()',{1}, '.','watlas','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
  if opts.atlas.dartel==1,
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Rigidly Registered Atlas Image';
      cdep(end).src_output = substruct('()',{1}, '.','ratlas','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
  if opts.atlas.dartel==2,
      cdep(end+1)          = cfg_dep;
      cdep(end).sname      = 'Affine Registered Atlas Image';
      cdep(end).src_output = substruct('()',{1}, '.','aatlas','()',{':'});
      cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
  end;
end

% jacobian
if opts.jacobian.warped,
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Jacobian Determinant Image';
    cdep(end).src_output = substruct('()',{1}, '.','jacobian','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end;

% warps
if opts.warps(1),
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Deformation Field';
    cdep(end).src_output = substruct('()',{1}, '.','fordef','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end;
if opts.warps(2),
    cdep(end+1)          = cfg_dep;
    cdep(end).sname      = 'Inverse Deformation Field';
    cdep(end).src_output = substruct('()',{1}, '.','invdef','()',{':'});
    cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end;

% tissues
for i=1:numel(tissue),
    if isfield(tissue(i),'native') && tissue(i).native(1),
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('p%d Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','p','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if isfield(tissue(i),'native') && tissue(i).native(2),
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('rp%d rigid Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','rp','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if isfield(tissue(i),'native') && tissue(i).native(3),
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('rp%d affine Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','rpa','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if tissue(i).warped(1),
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('wp%d Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','wp','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if tissue(i).warped(2),
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('mwp%d Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','mwp','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
    if tissue(i).warped(3),
        cdep(end+1)          = cfg_dep;
        cdep(end).sname      = sprintf('m0wp%d Image',i);
        cdep(end).src_output = substruct('.','tiss','()',{i},'.','m0wp','()',{':'});
        cdep(end).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    end
end



dep = cdep;




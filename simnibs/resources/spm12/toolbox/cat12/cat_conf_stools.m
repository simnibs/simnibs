function stools = cat_conf_stools(expert)
%_______________________________________________________________________
% wrapper for calling CAT surface utilities
%_______________________________________________________________________
% Robert Dahnke and Christian Gaser
% $Id: cat_conf_stools.m 1271 2018-02-08 14:28:51Z gaser $
%_______________________________________________________________________

  % try to estimate number of processor cores
  try
    numcores = max(feature('numcores'),1);
  catch
    numcores = 0;
  end
  
  if isdeployed, numcores = 0; end
  
  if ~exist('expert','var')
    expert = 0; % switch to de/activate further GUI options
  end

  % parallelize
  % ____________________________________________________________________
  nproc         = cfg_entry;
  nproc.tag     = 'nproc';
  nproc.name    = 'Split job into separate processes';
  nproc.strtype = 'w';
  nproc.val     = {numcores};
  nproc.num     = [1 1];
  nproc.help    = {
    'In order to use multi-threading the CAT12 segmentation job with multiple subjects can be split into separate processes that run in the background. You can even close Matlab, which will not affect the processes that will run in the background without GUI. If you do not want to run processes in the background then set this value to 0.'
      ''
      'Please note that no additional modules in the batch can be run except CAT12 segmentation. Any dependencies will be broken for subsequent modules.'
    };
 
  
  % do not process, if result already exist
  % ____________________________________________________________________
  lazy         = cfg_menu;
  lazy.tag     = 'lazy';
  lazy.name    = 'Lazy processing';
  lazy.labels  = {'Yes','No'};
  lazy.values  = {1,0};
  lazy.val     = {0};
  lazy.help    = {
    'Do not process data if the result exist. '
  };


%% Surface correlation and quality check
%-----------------------------------------------------------------------
  data_surf_cov         = cfg_files;
  data_surf_cov.tag     = 'data_surf';
  data_surf_cov.name    = 'Sample';
  data_surf_cov.filter  = 'gifti';
  data_surf_cov.ufilter = 'resampled';
  data_surf_cov.num     = [3 Inf];
  data_surf_cov.help    = {'Select resampled surfaces parameter files.'};
  
  data_xml = cfg_files;
  data_xml.name = 'Quality measures (optional)';
  data_xml.tag  = 'data_xml';
  data_xml.filter = 'xml';
  data_xml.ufilter = '^cat_.*';
  data_xml.val  = {{''}};
  data_xml.num  = [0 Inf];
  data_xml.help   = {...
'Select optional the quality measures that are saved during segmentation as xml-files in the report folder. This additionally allows to analyze image quality parameters such as noise, and bias. Please note, that the order of the xml-files should be the same as the other data files.'};

  sample_cov         = cfg_repeat;
  sample_cov.tag     = 'sample';
  sample_cov.name    = 'Data';
  sample_cov.values  = {data_surf_cov};
  sample_cov.num     = [1 Inf];
  sample_cov.help = {...
  'Specify data for each sample. If you specify different samples the mean correlation is displayed in separate boxplots for each sample.'};

  c         = cfg_entry;
  c.tag     = 'c';
  c.name    = 'Vector';
  c.help    = {'Vector of nuisance values'};
  c.strtype = 'r';
  c.num     = [Inf Inf];

  nuisance         = cfg_repeat;
  nuisance.tag     = 'nuisance';
  nuisance.name    = 'Nuisance variable';
  nuisance.values  = {c};
  nuisance.num     = [0 Inf];
  nuisance.help    = {...
  'This option allows for the specification of nuisance effects to be removed from the data. A potential nuisance parameter can be age. In this case the variance explained by age will be removed prior to the calculation of the correlation.'};

  check_mesh_cov      = cfg_exbranch;
  check_mesh_cov.tag  = 'check_mesh_cov';
  check_mesh_cov.name = 'Check sample homogeneity of surfaces';
  check_mesh_cov.val  = {sample_cov,data_xml,nuisance};
  check_mesh_cov.prog = @cat_stat_check_cov;
  check_mesh_cov.help = {
  'In order to identify surfaces with poor image quality or even artefacts you can use this function. Surfaces measures have to be first resampled to the template space (e.g. resampled and smoothed data) and can be then check for each hemisphere separately. Please do not mix data from both hemispehre.'
  ''
  'The idea of this tool is to check the correlation of all files across the sample. The correlation is calculated between all surfaces measures and the mean for each surface measure is plotted using a boxplot and optionally the indicated filenames. The smaller the mean correlation the more deviant is this surface measures from the sample mean. In the plot outliers from the sample are usually isolated from the majority of surface measures which are clustered around the sample mean. The mean correlation is plotted at the y-axis and the x-axis reflects the order of your measures.'};

    
  
  
%% surface measures
%-----------------------------------------------------------------------  
  data_surf_extract         = cfg_files;
  data_surf_extract.tag     = 'data_surf';
  data_surf_extract.name    = 'Central Surfaces';
  data_surf_extract.filter  = 'gifti';
  data_surf_extract.ufilter = '^lh.central';
  data_surf_extract.num     = [1 Inf];
  data_surf_extract.help    = {'Select left surfaces to extract values.'};
  
  GI        = cfg_menu;
  GI.name   = 'Gyrification index';
  GI.tag    = 'GI';
  GI.labels = {'No','Yes'};
  GI.values = {0,1};
  GI.val    = {1};
  GI.help   = {
    'Extract gyrification index (GI) based on absolute mean curvature. The method is described in Luders et al. NeuroImage, 29: 1224-1230, 2006.'
  };


  if expert
    area        = cfg_menu;
    area.name   = 'surface area';
    area.tag    = 'area';
    area.labels = {'No','Yes'};
    area.values = {0,1};
    area.val    = {1};
    area.help   = {
      'Estimation of surface area.'
    };
  
    % Two cases ... one with the stardard average surface and one were
    % another average surface metric file "area" of a group can be choosen.
    GIA        = cfg_menu;
    GIA.name   = 'Average gyrification index (only resampled output)';
    GIA.tag    = 'GIA';
    GIA.labels = {'No','Yes'};
    GIA.values = {0,1};
    GIA.val    = {1};
    GIA.help   = {
      'WARNING: This GI measures is still in development and not varified yet!\n\n It extract the average gyrification index (AGI) as local area relation between the individual central and group average surface (FS average).'
    };
  
    GII        = cfg_menu;
    GII.name   = 'Inflating gyrification index ';
    GII.tag    = 'GII';
    GII.labels = {'No','Yes'};
    GII.values = {0,1};
    GII.val    = {1};
    GII.help   = {
      'WARNING: This GI measures is still in development and not varified yet!\n\n It extract a inflating gyrification index (IGI) as local area relation between the individual central and an inflated version of it.'
    };
  
    GIS        = cfg_menu;
    GIS.name   = 'Spherical gyrification index (only resampled output)';
    GIS.tag    = 'GIS';
    GIS.labels = {'No','Yes'};
    GIS.values = {0,1};
    GIS.val    = {1};
    GIS.help   = {
      'WARNING: This GI measures is still in development and not varified yet!\n\n It extract a spherical mapping based gyrification index (SGI) as local area relation between the individual central and the hull surface.'
    };
    
    GIL        = cfg_menu;
    GIL.name   = 'Laplacian gyrification index';
    GIL.tag    = 'GIL';
    GIL.labels = {'No','Yes'};
    GIL.values = {0,1};
    GIL.val    = {1};
    GIL.help   = {
      'WARNING: This GI measures is still in development and not varified yet!\n\n Extract Laplacian gyrification index (LGI) as local area relation between the individual central and the hull surface [Dahnke:2010].'
    };
  
    OS        = cfg_menu;
    OS.name   = 'Outer surface';
    OS.tag    = 'OS';
    OS.labels = {'No','Yes'};
    OS.values = {0,1};
    OS.val    = {1};
    OS.help   = {
      'Creates outer surface by moving each CS vertex by its half thickness along the surface normal.'
    };
  
    IS        = cfg_menu;
    IS.name   = 'Inner surface';
    IS.tag    = 'IS';
    IS.labels = {'No','Yes'};
    IS.values = {0,1};
    IS.val    = {1};
    IS.help   = {
      'Creates outer surface by moving each CS vertex by its half thickness along the surface normal.'
    };
  end

  
  FD        = cfg_menu;
  FD.name   = 'Cortical complexity (fractal dimension)';
  FD.tag    = 'FD';
  FD.labels = {'No','Yes'};
  FD.values = {0,1};
  FD.val    = {0};
  FD.help   = {
    'Extract cortical complexity (fractal dimension) which is described in Yotter et al. Neuroimage, 56(3): 961-973, 2011.'
    ''
    'Warning: Estimation of cortical complexity is very slow!'
    ''
  };

  SD        = cfg_menu;
  SD.name   = 'Sulcus depth';
  SD.tag    = 'SD';
  SD.labels = {'No','Yes'};
  SD.values = {0,1};
  SD.val    = {1};
  SD.help   = {
    'Extract sqrt-transformed sulcus depth based on the euclidean distance between the central surface and its convex hull.'
    ''
    'Transformation with sqrt-function is used to render the data more normally distributed.'
    ''
  };

  SA        = cfg_menu;
  SA.name   = 'Surface area';
  SA.tag    = 'SA';
  SA.labels = {'No','Yes'};
  SA.values = {0,1};
  SA.val    = {1};
  SA.help   = {
    'Extract log10-transformed local surface area using re-parameterized tetrahedral surface. The method is described in Winkler et al. NeuroImage, 61: 1428-1443, 2012.'
    ''
    'Log-transformation is used to render the data more normally distributed.'
    ''
  };


  surfextract      = cfg_exbranch;
  surfextract.tag  = 'surfextract';
  surfextract.name = 'Extract additional surface parameters';
  if expert > 1
    surfextract.val  = {data_surf_extract,area,GI,GIA,GII,GIL,GIS,FD,SD,IS,OS,nproc,lazy};
  else
    surfextract.val  = {data_surf_extract,GI,FD,SD,nproc};
  end
  surfextract.prog = @cat_surf_parameters;
  surfextract.vout = @vout_surfextract;
  surfextract.help = {'Additional surface parameters can be extracted that can be used for statistical analysis.'};




%% map volumetric data
%-----------------------------------------------------------------------  
  datafieldname         = cfg_entry;
  datafieldname.tag     = 'datafieldname';
  datafieldname.name    = 'Output Name';
  datafieldname.strtype = 's';
  datafieldname.num     = [1 Inf];
  datafieldname.val     = {'intensity'};
  datafieldname.help    = {
    'Name that is prepended to the filename of the mapped volume.'
    ''
    };
 
  interp         = cfg_menu;
  interp.tag     = 'interp';
  interp.name    = 'Interpolation Type';
  interp.labels  = {'Nearest neighbour','Linear','Cubic'};
  interp.values  = {{'nearest_neighbour'},{'linear'},{'cubic'}};
  interp.val     = {{'linear'}};
  interp.help    = {
    'Volume extraction interpolation type. '
    ' -linear:            Use linear interpolation (default).'
    ' -nearest_neighbour: Use nearest neighbour interpolation.'
    ' -cubic:             Use cubic interpolation.'
    ''
  };
  
  % sample function 
  sample         = cfg_menu;
  sample.tag     = 'sample';
  sample.name    = 'Sample Function';
  sample.labels  = {'Mean','Weighted mean','Maximum','Minimum','Absolute maximum','Multi-values'};
  sample.values  = {{'avg'},{'weighted_avg'},{'max'},{'min'},{'maxabs'},{'multi'}};
  sample.val     = {{'maxabs'}};
  sample.help    = {
    'Sample function to combine the values of the grid along the surface normals.'
    ' Mean:          Use average for mapping along normals.'
    ' Weighted mean: Use weighted average with gaussian kernel for mapping along normals. The kernel is so defined that values at the boundary are weighted with 50% while center is weighted with 100% (useful for (r)fMRI data.'
    ' Maximum:       Use maximum value for mapping along normals.'
    ' Minimum:       Use minimum value for mapping along normals.'
    ' Absolute maximum: Use absolute maximum value for mapping along normals (useful for mapping contrast images from 1st-level fMRI analysis).'
    ' Multi-values:  Map data for each grid step separately and save files with indicated grid value. Please note that this option is intended for high-resolution (f)MRI data only (e.g. 0.5mm voxel size).'
    ''
  };

  %% -- sampling points and average function
 
  % startpoint
  abs_startpoint         = cfg_entry;
  abs_startpoint.tag     = 'startpoint';
  abs_startpoint.name    = 'Startpoint';
  abs_startpoint.strtype = 'r';
  abs_startpoint.val     = {-0.5};
  abs_startpoint.num     = [1 1];
  abs_startpoint.help    = {
    'Absolute position of the start point of the grid along the surface normals in mm according to the surface. Give negative value for a start point outside of the surface (CSF direction, outwards). '
  };
  rel_startpoint = abs_startpoint;
  rel_startpoint.val     = {-0.5};
  rel_startpoint.help    = {
    'Relative position of the start point of the grid along the surface normals according to a tissue class. A value of "-0.5" begins at the GM/CSF border and even lower values define a start point outside of the tissue class (CSF direction, outwards). A value of "0" means that the central surface is used as starting point and "0.5" is related to the GM/WM border.'
  };
  
  % steps
  abs_steps         = cfg_entry;
  abs_steps.tag     = 'steps';
  abs_steps.name    = 'Steps';
  abs_steps.strtype = 'w';
  abs_steps.val     = {7};
  abs_steps.num     = [1 1];
  abs_steps.help    = {
    'Number of grid steps. '
  };
  rel_steps = abs_steps; 

  % endpoint
  abs_endpoint         = cfg_entry;
  abs_endpoint.tag     = 'endpoint';
  abs_endpoint.name    = 'Endpoint';
  abs_endpoint.strtype = 'r';
  abs_endpoint.val     = {0.5};
  abs_endpoint.num     = [1 1];
  abs_endpoint.help    = {
    'Absolute position of the end point of the grid along the surface normals (pointing inwards) in mm according to the surface. '
  };
  rel_endpoint = abs_endpoint;
  rel_endpoint.val     = {0.5};
  rel_endpoint.help    = {
    'Relative position of the end point of the grid along the surface normals (pointing inwards) according to a tissue class. A value of "0.5" ends at the GM/WM border and values > 0.5 define an end point outside of the tissue class (WM direction, inwards). A value of "0" ends at the central surface.'
  };

  % tissue class
  rel_class         = cfg_menu;
  rel_class.tag     = 'class';
  rel_class.name    = 'Tissue Class';
  rel_class.labels  = {'GM'};
  rel_class.values  = {'GM'};
  rel_class.val     = {'GM'};
  rel_class.help    = {
    'Tissue class for which the relative positions are estimated.'
  };

  % tissue class
  abs_class         = cfg_menu;
  abs_class.tag     = 'surface';
  abs_class.name    = 'Surface';
  abs_class.labels  = {'WM Surface','Central Surface','Pial Surface'};
  abs_class.values  = {'WM','Central','Pial'};
  abs_class.val     = {'Central'};
  abs_class.help    = {
    'Surface (or tissue boundary) for which the absolute positions are estimated.'
  };

  % absolute position
  abs_mapping         = cfg_branch;
  abs_mapping.tag     = 'abs_mapping';
  abs_mapping.name    = 'Absolute Grid Position From a Surface';
  abs_mapping.val   = {
    abs_class ...
    abs_startpoint ...
    abs_steps ...
    abs_endpoint ...
  }; 
  tissue.help    = {
    'Map volumetric data from abolute grid positions from a surface (or tissue boundary).'
  };
  
  %% relative mapping with equi-distance approach
  rel_mapping         = cfg_branch;
  rel_mapping.tag     = 'rel_mapping';
  rel_mapping.name    = 'Relative Grid Position Within a Tissue Class (Equi-distance Model)';
  rel_mapping.val   = {
    rel_class ...
    rel_startpoint ...
    rel_steps ...
    rel_endpoint ...
  };
  rel_mapping.help    = {
    'Map volumetric data from relative grid positions within a tissue class using equi-distance approach. Here, the grid lines have equal distances in between the tissue.'
  };

  %% relative mapping with equi-volume approach
  rel_equivol_mapping         = cfg_branch;
  rel_equivol_mapping.tag     = 'rel_equivol_mapping';
  rel_equivol_mapping.name    = 'Relative Grid Position Within a Tissue Class (Equi-volume Model)';
  rel_equivol_mapping.val   = {
    rel_class ...
    rel_startpoint ...
    rel_steps ...
    rel_endpoint ...
  };
  rel_equivol_mapping.help    = {
    'Map volumetric data from relative positions within a tissue class using equi-volume approach. '
    'This option is using the approach by Bok (Z. Gesamte Neurol. Psychiatr. 12, 682–750, 1929). '
    'Here, the volume between the grids is constant. The correction is based on Waehnert et al. (NeuroImage, 93: 210-220, 2014).'
    'Please note that this option is intended for high-resolution (f)MRI data only'
    '' 
  };

  %% -- Mapping function

  mapping         = cfg_choice;
  mapping.tag     = 'mapping';
  mapping.name    = 'Mapping Function';
  mapping.values  = {
    abs_mapping ...
    rel_mapping ...
  }; 
  mapping.val = {rel_mapping};
  mapping.help    = {
    'Volume extration type. '
    '  Absolute Grid Position From a Surface (or Tissue Boundary):'
    '    Extract values around a surface or tissue boundary with a specified absolute sample '
    '    distance and either combine these values or save values separetely.'
    '  Relative Grid Position Within a Tissue Class (Equi-distance approach):' 
    '    Extract values within a tissue class with a specified relative sample distance'
    '    with equally distributed distances and either combine these values or save values separetely.'
    '' 
  };

  mapping_native = mapping;
  mapping_native.values{3} = rel_equivol_mapping;
  mapping_native.help    = {
    'Volume extration type. '
    '  Absolute Grid Position From a Surface (or Tissue Boundary):'
    '    Extract values around a surface or tissue boundary with a specified absolute sample '
    '    distance and either combine these values or save values separetely.'
    '  Relative Grid Position Within a Tissue Class (Equi-distance approach):' 
    '    Extract values within a tissue class with a specified relative sample distance'
    '    with equally distributed distances and either combine these values or save values separetely.'
    '  Relative Grid Position Within a Tissue Class (Equi-volume approach):' 
    '    Extract values within a tissue class with a specified relative sample distance'
    '    that is corrected for constant volume between the grids and either combine these values or save values separetely.'
    '' 
  };



% extract volumetric data in individual space 
%-----------------------------------------------------------------------  

  data_surf_sub_lh         = cfg_files;
  data_surf_sub_lh.tag     = 'data_mesh_lh';
  data_surf_sub_lh.name    = '(Left) Individual Surfaces';
  data_surf_sub_lh.filter  = 'gifti';
  data_surf_sub_lh.ufilter = '^lh.central.(?!nofix).*';
  data_surf_sub_lh.num     = [1 Inf];
  data_surf_sub_lh.help    = {
    'Select left subject surface files (do not select the *.nofix.* surface).'
    'Right side will be automatically processed.'
    };
   
  data_sub         = cfg_files; 
  data_sub.tag     = 'data_vol';
  data_sub.name    = '(Co-registered) Volumes in Native Space';
  data_sub.filter  = 'image';
  data_sub.ufilter = '^(?!wm|wp|m0wp|mwp|wc).*'; % no normalized images
  data_sub.num     = [1 Inf];
  data_sub.help    = {
    'Select volumes in native (subject) space.'
    'Please note that these images have to be in the same space as the T1-image that was used to extract the cortical surface. An optional co-registration might be necessary if you have functional or structual data that are not yet aligned to the T1-image.'
  };

  vol2surf      = cfg_exbranch;
  vol2surf.tag  = 'vol2surf';
  vol2surf.name = 'Map Volume (Native Space) to Individual Surface';
  if expert
    vol2surf.val = {
      data_sub ...
      data_surf_sub_lh ...
      sample ...
      interp ...
      datafieldname ...
      mapping_native ...
      };
  else
    vol2surf.val = {
      data_sub ...
      data_surf_sub_lh ...
      sample ...
      datafieldname ...
      mapping_native ...
    };
  end
  vol2surf.prog = @cat_surf_vol2surf;
  vol2surf.vout = @vout_vol2surf;
  vol2surf.help = {
    'Map volume (native space) to individual surface. These mapped volumes have to be finally resampled and smoothed before any statistical analysis.'
    ''
    'The ouput will be named:' 
    '  [rh|lh].OutputName_VolumeName' 
    ''
  };



%% extract volumetric data in template space
%-----------------------------------------------------------------------  
  merge_hemi         = cfg_menu;
  merge_hemi.tag     = 'merge_hemi';
  merge_hemi.name    = 'Merge hemispheres';
  merge_hemi.labels  = {
    'No -   save resampled data for each hemisphere',...
    'Yes -  merge hemispheres'
  };
  merge_hemi.values  = {0,1};
  merge_hemi.val     = {1};
  merge_hemi.help    = {
    'Meshes for left and right hemisphere can be merged to one single mesh. This simplifies the analysis because only one analysis has to be made for both hemispheres.'
    'Hoever, this also means that data size is double for one single analysis which might be too memory demanding for studies with several hundreds or even more files. If your model cannot be estimated due to memory issues you should not merge the resampled data.'
  };

  mesh32k         = cfg_menu;
  mesh32k.tag     = 'mesh32k';
  mesh32k.name    = 'Resample Size';
  mesh32k.labels  = {
    '32k  mesh (HCP)',...
    '164k mesh (Freesurfer)'
  };
  mesh32k.values  = {1,0};
  mesh32k.val     = {1};
  mesh32k.help    = {
    'Resampling can be done either to a higher resoluted 164k mesh that is compatible to Freesurfer data or to a lower resoluted 32k mesh (average vertex spacing of ~2 mm) that is compatible to the Human Connectome Project (HCP).'
    'The HCP mesh has the advantage of being processed and handled much faster and with less memory demands. Another advantage is that left and right hemispheres are aligned to optionally allow a direct comparison between hemispheres.'
  };


  data_surf_avg_lh         = cfg_files; 
  data_surf_avg_lh.tag     = 'data_mesh_lh';
  data_surf_avg_lh.name    = '(Left) Template Hemisphere';
  data_surf_avg_lh.filter  = 'gifti';
  data_surf_avg_lh.ufilter = '^lh.*';
  data_surf_avg_lh.num     = [1 1];
  data_surf_avg_lh.val{1}  = {fullfile(spm('dir'),'toolbox','cat12','templates_surfaces','lh.central.Template_T1_IXI555_MNI152_GS.gii')};
  data_surf_avg_lh.dir     = fullfile(spm('dir'),'toolbox','cat12');
  data_surf_avg_lh.help    = {
    'Select left template surface file. '
    'Right hemisphere will be automatically processed.'
    };
    
  data_norm         = cfg_files; 
  data_norm.tag     = 'data_vol';
  data_norm.name    = 'Spatially Normalized Volumes';
  data_norm.filter  = 'image';
  data_norm.ufilter = '.*';
  data_norm.num     = [1 Inf];
  data_norm.help    = {
    'Select spatially normalized volumes (in template space).'
    ''
  };

  vol2tempsurf      = cfg_exbranch;
  vol2tempsurf.tag  = 'vol2surftemp';
  vol2tempsurf.name = 'Map Normalized Volume to Template Surface';
  vol2tempsurf.val  = {
    data_norm ...
    merge_hemi ...
    mesh32k ...
    sample ...
    interp ...
    datafieldname ...
    mapping ...
  };
  vol2tempsurf.prog = @cat_surf_vol2surf;
  vol2tempsurf.vout = @vout_vol2surf;
  vol2tempsurf.help = {
    'Map spatially normalized data (in template space) to template surface.'
    'The template surface was generated by CAT12 surface processing of the average of 555 Dartel-normalized images of the IXI database that were also used to create the IXI Dartel template.   '
    ''
    'The ouput will be named:' 
    '  [rh|lh|mesh].OutputName_VolumeName.Template_T1_IXI555_MNI152_GS.gii' 
    ''
    '  WARNING: This function is primarily thought in order to project statistical results '
    '           to the template surface. Do not use the output of this function'
    '           for any statistical analysis. Statistical analysis should only use'
    '           resampled data of individual volume data (in native space) that were mapped'
    '           to individual surfaces.'
    ''
  };






%  region-based measures in subject space
%  ---------------------------------------------------------------------
if expert>1
  %%
  % surface files
    ROI.sdata         = cfg_files;
    ROI.sdata.tag     = 'sdata';
    ROI.sdata.name    = '(Left) Surface Data Files';
    ROI.sdata.filter  = 'gifti';
    ROI.sdata.ufilter = 'lh.central.*';
    ROI.sdata.num     = [1 Inf];
    ROI.sdata.help    = {'Surface data sample. Both sides will be processed'};

  % ROI files
    ROI.atlas         = cfg_files;
    ROI.atlas.tag     = 'rdata';
    ROI.atlas.name    = '(Left) ROI atlas files';
    ROI.atlas.filter  = 'any';
    if expert
      ROI.atlas.ufilter = 'lh.aparc.*';
    else
      ROI.atlas.ufilter = 'lh.aparc(_DKT40JT|_a2009s)\..*';
    end
    ROI.atlas.dir     = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces'); 
    ROI.atlas.num     = [1 Inf];
    ROI.atlas.help    = {'These are the ROI atlas files. Both sides will be processed.'};

  % check fields of different measures
    ROI.rLGI        = cfg_menu;
    ROI.rLGI.name   = 'region-wise Laplacian-based gyrification index (rLGI)';
    ROI.rLGI.tag    = 'rLGI';
    ROI.rLGI.labels = {'No','Yes'};
    ROI.rLGI.values = {0,1};
    ROI.rLGI.val    = {1};
    ROI.rLGI.help   = {[
        'The region-wise Laplacian-based gyrification index describes the surface area relation' ...
        'of a atlas region on the central surface vs. the hull surface that was generated by mapping' ...
        'the central surface by the Laplacian approach to the position of the hull what takes arount 5 minutes per subject. '
      ]};
    
    ROI.area        = cfg_menu;
    ROI.area.name   = 'surface area';
    ROI.area.tag    = 'area';
    ROI.area.labels = {'No','Yes'};
    ROI.area.values = {0,1};
    ROI.area.val    = {1};
    ROI.area.help   = {
        'Region-wise surface area.'
      };
    
    ROI.rarea        = cfg_menu;
    ROI.rarea.name   = 'relative surface area';
    ROI.rarea.tag    = 'rarea';
    ROI.rarea.labels = {'No','Yes'};
    ROI.rarea.values = {0,1};
    ROI.rarea.val    = {1};
    ROI.rarea.help   = {
        'Region-wise surface area, normalized by total area.'
      };
    
    % volume ??? comparison to VBM
    % * map all GM/WM/CSF voxel to there closest CS vertex
    % * create inner and outer surface to used delaunay for subregions ...
    
    
  % main call
    estroi      = cfg_exbranch;
    estroi.tag  = 'estroi';
    estroi.name = 'Surface-based ROI measures in subject space';
    estroi.val  = {
      ROI.sdata ...
      ROI.atlas ...
      ROI.area ...
      ROI.rarea ...
      }; 
    if expert>1
      estroi.val{end+1} = ROI.rLGI; 
      estroi.val{end+1} = nproc; 
      estroi.val{end+1} = lazy;
    end
    %estroi.prog = @cat_roi_parameters;
    estroi.help = {
      'Surface-based ROI measures that required estimation on the original rather than the template surface mesh, such as the area or the gyrification index.'
    };
end

%% surface to ROI (in template space)
%  ---------------------------------------------------------------------
%  * CxN  cdata files [ thickness , curvature , ... any other file ] in template space [ resampled ]
%  * M    atlas files [ choose files ]
%  * average measures [ mean , std , median , max , min , mean95p ]
% 
%  - csv export (multiple files) > catROIs_ATLAS_SUBJECT
%  - xml export (one file)       > catROIs_ATLAS_SUBJECT
%  ---------------------------------------------------------------------
%  surface ROIs have sidewise index?!
%  ---------------------------------------------------------------------

% set of cdata files
  if expert && 0 % this is not ready now
    s2r.cdata         = cfg_files;
    s2r.cdata.tag     = 'cdata';
    s2r.cdata.name    = '(Left) Surface Data Files';
    s2r.cdata.filter  = 'any';
    s2r.cdata.ufilter = 'lh.(?!cent|sphe|defe|hull).*';
    s2r.cdata.num     = [1 Inf];
    s2r.cdata.help    = {'Surface data sample. Both sides will be processed'};
  else % only smoothed/resampled
    s2r.cdata         = cfg_files;
    s2r.cdata.tag     = 'cdata';
    s2r.cdata.name    = '(Left) Surface Data Files';
    s2r.cdata.filter  = 'any';
    s2r.cdata.ufilter = '^lh.(?!cent|sphe|defe|hull).*';
    s2r.cdata.num     = [1 Inf];
    s2r.cdata.help    = {'Surface data sample. Both sides will be processed'};
  end
  
  s2r.cdata_sample         = cfg_repeat;
  s2r.cdata_sample.tag     = 'cdata_sub.';
  s2r.cdata_sample.name    = 'Surface Data Sample';
  s2r.cdata_sample.values  = {s2r.cdata};
  s2r.cdata_sample.num     = [1 Inf];
  s2r.cdata_sample.help = {[...
    'Specify data for each sample (i.e. thickness, gyrification, ...). ' ...
    'All samples must have the same size and same order. ' ...
    ''
  ]};

% ROI files
  s2r.ROIs         = cfg_files;
  s2r.ROIs.tag     = 'rdata';
  s2r.ROIs.name    = '(Left) ROI atlas files';
  s2r.ROIs.filter  = 'any';
  if expert 
    s2r.ROIs.ufilter = 'lh.aparc.*';
  else
    s2r.ROIs.ufilter = 'lh.aparc_(a2009s|DK40|HCP_MMP1).*'; % not yet working for all atlases
  end
  s2r.ROIs.dir     = fullfile(spm('dir'),'toolbox','cat12','atlases_surfaces'); 
  s2r.ROIs.num     = [1 Inf];
  s2r.ROIs.help    = {'These are the ROI atlas files. Both sides will be processed.'};

% ROI area
  s2r.area         = cfg_menu;
  s2r.area.tag     = 'area';
  s2r.area.name    = 'Estimate ROI Area';
  s2r.area.labels  = {'No','Yes'};
  s2r.area.values  = {0,1};
  s2r.area.val     = {1}; 
  s2r.area.help    = {'Estimate area of each ROI.'};
 
% ROI area
  s2r.vernum         = cfg_menu;
  s2r.vernum.tag     = 'vernum';
  s2r.vernum.name    = 'Count ROI Vertices';
  s2r.vernum.labels  = {'No','Yes'};
  s2r.vernum.values  = {0,1};
  s2r.vernum.val     = {1}; 
  s2r.vernum.help    = {'Count vertices of each ROI.'};
  
% average mode within a ROI  
  % mean
  s2r.avg.mean         = cfg_menu;
  s2r.avg.mean.tag     = 'mean';
  s2r.avg.mean.name    = 'Mean Estimation';
  s2r.avg.mean.labels  = {'No','Yes'};
  s2r.avg.mean.values  = {0,1};
  s2r.avg.mean.val     = {1}; 
  s2r.avg.mean.help    = {'Set mean value estimation per ROI.'}; 
  % std
  s2r.avg.std         = cfg_menu;
  s2r.avg.std.tag     = 'std';
  s2r.avg.std.name    = 'STD Estimation';
  s2r.avg.std.labels  = {'No','Yes'};
  s2r.avg.std.values  = {0,1};
  s2r.avg.std.val     = {1}; 
  s2r.avg.std.help    = {'Set standard deviation estimation per ROI.'}; 
  % min
  s2r.avg.min         = cfg_menu;
  s2r.avg.min.tag     = 'min';
  s2r.avg.min.name    = 'Minimum Estimation';
  s2r.avg.min.labels  = {'No','Yes'};
  s2r.avg.min.values  = {0,1};
  s2r.avg.min.val     = {0}; 
  s2r.avg.min.help    = {'Set minimum estimation per ROI.'};   
  % max
  s2r.avg.max         = cfg_menu;
  s2r.avg.max.tag     = 'max';
  s2r.avg.max.name    = 'Maximum Estimation';
  s2r.avg.max.labels  = {'No','Yes'};
  s2r.avg.max.values  = {0,1};
  s2r.avg.max.val     = {0}; 
  s2r.avg.max.help    = {'Set maximum estimation per ROI.'};   
  % median
  s2r.avg.median         = cfg_menu;
  s2r.avg.median.tag     = 'median';
  s2r.avg.median.name    = 'Median Estimation';
  s2r.avg.median.labels  = {'No','Yes'};
  s2r.avg.median.values  = {0,1};
  s2r.avg.median.val     = {0}; 
  s2r.avg.median.help    = {'Set median estimation per ROI.'};   
  % all functions
  s2r.avg.main         = cfg_branch;
  s2r.avg.main.tag     = 'avg';
  s2r.avg.main.name    = 'ROI Average Functions';
  s2r.avg.main.val     = {
    s2r.avg.mean ...
    s2r.avg.std ...
    s2r.avg.min ...
    s2r.avg.max ...
    s2r.avg.median ...
  };

%% main function
  surf2roi      = cfg_exbranch;
  surf2roi.tag  = 'surf2roi';
  surf2roi.name = 'Extract ROI-based surface values';
  switch expert
  case 2
    surf2roi.val  = {
      s2r.cdata_sample ...
      s2r.ROIs ...
      nproc ... 
      s2r.avg.main};
  case {0, 1}
    surf2roi.val  = {s2r.cdata_sample};
  end
  surf2roi.prog = @cat_surf_surf2roi;
  surf2roi.vout = @vout_surf_surf2roi;
  surf2roi.help = {
    'While ROI-based values for VBM (volume) data are automatically saved in the label folder as XML file it is necessary to additionally extract these values for surface data. This has to be done after preprocessing the data and creating cortical surfaces. '
    ''
    'You can extract ROI-based values for cortical thickness but also for any other surface parameter that was extracted using the ''Extract Additional Surface Parameters'' function.'
    ''
    'Please note that these values are extracted from data in native space without any smoothing. As default the mean inside a ROI is calculated and saved as XML file in the label folder.'
  };

%% roi to surface  
%  ---------------------------------------------------------------------

% ROI files
  r2s.ROIs         = cfg_files;
  r2s.ROIs.tag     = 'rdata';
  r2s.ROIs.name    = 'ROI atlas files';
  r2s.ROIs.filter  = 'xml';
  r2s.ROIs.ufilter = '.*';
  r2s.ROIs.dir     = fullfile(spm('dir'),'toolbox','cat12','templates_1.50mm'); 
  r2s.ROIs.num     = [1 Inf];
  r2s.ROIs.help    = {'These are the indivudal ROI atlas files from the label directory. Choose XML files.'};

 % atlas used for extraction .. if xml 
  r2s.atlas         = cfg_entry;
  r2s.atlas.tag     = 'atlas';
  r2s.atlas.name    = 'Atlas maps';
  r2s.atlas.help    = {'Enter the name of the atlas maps, e.g. "LPBA40", or "all".'};
  r2s.atlas.strtype = 's';
  r2s.atlas.val     = {'all'};
  r2s.atlas.num     = [1 Inf];
  
% datafield in the CSV or XML file eg the GM thickness
  r2s.data         = cfg_entry;
  r2s.data.tag     = 'fields';
  r2s.data.name    = 'Datafields';
  r2s.data.help    = {'Enter the name of the data fields, e.g. "Vgm", "Igm", "Tgm", "mySurfData" or "all"'};
  r2s.data.strtype = 's';
  r2s.data.val     = {'all'};
  r2s.data.num     = [1 Inf];
  
% surface 
  r2s.surf        = cfg_menu;
  r2s.surf.name   = 'Surface type for mapping';
  r2s.surf.tag    = 'surf';
  r2s.surf.labels = {'FS average','Dartel average','subject'};
  r2s.surf.values = {'freesurfer','dartel','subject'};
  r2s.surf.val    = {'freesurfer'};
  r2s.surf.help   = {'Surface type for value projection. '};
 
% average function? - not required 
 
%% main function
  roi2surf      = cfg_exbranch;
  roi2surf.tag  = 'roi2surf';
  roi2surf.name = 'Map ROI data to the surface';
  roi2surf.val  = {
    r2s.ROIs ...
    r2s.atlas ...
    r2s.data ...
    r2s.surf ...
    };
  roi2surf.prog = @cat_roi_roi2surf;
  roi2surf.help = {
    ''
  };
%% surface calculations 
%  ---------------------------------------------------------------------
%  estimation per subject (individual and group sampling):
%  g groups with i datafiles and i result datafile
 
  sc.cdata         = cfg_files;
  sc.cdata.tag     = 'cdata';
  sc.cdata.name    = 'Surface Data Files';
  sc.cdata.filter  = 'any';
  sc.cdata.ufilter = '(lh|rh|mesh).*';
  sc.cdata.num     = [1 Inf];
  sc.cdata.help    = {'These are the surface data files that are used by the calculator.  They are referred to as s1, s2, s3, etc in the order they are specified.'};
  
  sc.cdata_sub         = cfg_files;
  sc.cdata_sub.tag     = 'cdata';
  sc.cdata_sub.name    = 'Surface Data Files';
  sc.cdata_sub.filter  = 'gifti';
  sc.cdata_sub.ufilter = '(lh|rh|mesh).(?!cent|sphe|defe|hull).*gii';
  sc.cdata_sub.num     = [1 Inf];
  sc.cdata_sub.help    = {'These are the surface data files that are used by the calculator.  They are referred to as s1, s2, s3, etc in the order they are specified.'};
   
  sc.cdata_sample         = cfg_repeat;
  sc.cdata_sample.tag     = 'cdata_sub.';
  sc.cdata_sample.name    = 'Surface Data Sample';
  sc.cdata_sample.values  = {sc.cdata_sub};
  sc.cdata_sample.num     = [1 Inf];
  sc.cdata_sample.help = {...
    'Specify data for each sample. All samples must have the same size and same order.'};
 
  sc.outdir         = cfg_files;
  sc.outdir.tag     = 'outdir';
  sc.outdir.name    = 'Output Directory';
  sc.outdir.filter  = 'dir';
  sc.outdir.ufilter = '.*';
  sc.outdir.num     = [0 1];
  sc.outdir.val{1}  = {''};
  sc.outdir.help    = {
    'Files produced by this function will be written into this output directory.  If no directory is given, images will be written  to current working directory.  If both output filename and output directory contain a directory, then output filename takes precedence.'
  };
  
  sc.surfname         = cfg_entry;
  sc.surfname.tag     = 'dataname';
  sc.surfname.name    = 'Output Filename';
  sc.surfname.strtype = 's';
  sc.surfname.num     = [1 Inf];
  sc.surfname.val     = {'output'};
  sc.surfname.help    = {'The output surface data file is written to current working directory unless a valid full pathname is given.  If a path name is given here, the output directory setting will be ignored.'};
 
  sc.dataname         = cfg_entry;
  sc.dataname.tag     = 'dataname';
  sc.dataname.name    = 'Texture Name';
  sc.dataname.strtype = 's';
  sc.dataname.num     = [1 Inf];
  sc.dataname.val     = {'output'};
  sc.dataname.help    = {
    'Name of the texture as part of the filename.'
    ''
    '  [rh|lh].TEXTURENAME[.resampled].subjectname[.gii]' 
  };

  sc.expression         = cfg_entry;
  sc.expression.tag     = 'expression';
  sc.expression.name    = 'Expression';
  sc.expression.strtype = 's';
  sc.expression.num     = [1 Inf];
  sc.expression.val     = {'s1'};
  sc.expression.help    = {
    'Example expressions (f):'
    '  * Mean of six surface textures (select six texture files)'
    '    f = ''(s1+s2+s3+s4+s5+s6)/6'''
    '  * Make a binaray mask texture at threshold of 100'
    '    f = ''(s1>100)'''
    '  * Make a mask from one texture and apply to another'
    '    f = ''s2.*(s1>100)'''
    '        - here the first texture is used to make the mask, which is applied to the second texture'
    '  * Sum of n texures'
    '    f = ''s1 + s2 + s3 + s4 + s5 + ...'''
    '  * Sum of n texures (when reading data into a data-matrix - use dmtx arg)'
    '    f = mean(S)'
    ''
  };

  sc.dmtx         = cfg_menu;
  sc.dmtx.tag     = 'dmtx';
  sc.dmtx.name    = 'Data Matrix';
  sc.dmtx.labels  = {
    'No - don''t read images into data matrix',...
    'Yes -  read images into data matrix'
  };
  sc.dmtx.values  = {0,1};
  sc.dmtx.val     = {0};
  sc.dmtx.help    = {
    'If the dmtx flag is set, then textures are read into a data matrix S (rather than into separate variables s1, s2, s3,...). The data matrix should be referred to as S, and contains textures in rows. Computation is vertex by vertex, S is a NxK matrix, where N is the number of input textures, and K is the number of vertices per plane.'
  };


% ----------------------------------------------------------------------

  surfcalc      = cfg_exbranch;
  surfcalc.tag  = 'surfcalc';
  surfcalc.name = 'Surface Calculator';
  surfcalc.val  = {
    sc.cdata ...
    sc.surfname ...
    sc.outdir ...
    sc.expression ...
    sc.dmtx ...
  };
  surfcalc.prog = @cat_surf_calc;
  surfcalc.help = {
    'Mathematical operations for surface data (textures).'
    'It works similar to ''spm_imcalc''.  The input surface data must have the same number of entries (e.g. data of the same hemisphere of a subject or resampled data).'
  };


  surfcalcsub      = cfg_exbranch;
  surfcalcsub.tag  = 'surfcalcsub';
  surfcalcsub.name = 'Surface Calculator (subject-wise)';
  surfcalcsub.val  = {
    sc.cdata_sample ...
    sc.dataname ...
    sc.outdir ...
    sc.expression ...
    sc.dmtx ...
  };
  surfcalcsub.prog = @cat_surf_calc;
  surfcalcsub.help = {
    'Mathematical operations for surface data sets (textures).'
    'In contrast to the ''Surface Calculator'' it allows to apply the same expression to multiple subjects. Please note that a fixed name structure is expected: [rh|lh].TEXTURENAME[.resampled].subjectname[.gii]. Here TEXTURENAME will be replaced by the output name.'
  };


%% Resample and smooth surfaces 
%-----------------------------------------------------------------------
  data_surf         = cfg_files;
  data_surf.tag     = 'data_surf';
  data_surf.name    = '(Left) Surfaces Data';
  data_surf.filter  = 'any';
  if expert > 1
    data_surf.ufilter = '^lh.';
  else
    data_surf.ufilter = '^lh.(?!cent|sphe|defe|hull).*';
  end
  data_surf.num     = [1 Inf];
  data_surf.help    = {'Select surfaces data files for left hemisphere for resampling to template space.'};

  fwhm_surf         = cfg_entry;
  fwhm_surf.tag     = 'fwhm_surf';
  fwhm_surf.name    = 'Smoothing Filter Size in FWHM';
  fwhm_surf.strtype = 'r';
  fwhm_surf.num     = [1 1];
  fwhm_surf.val     = {15};
  fwhm_surf.help    = {
    'Select filter size for smoothing. For cortical thickness a good starting value is 15mm, while other surface parameters based on cortex folding (e.g. gyrification, cortical complexity) need a larger filter size of about 20-25mm. For no filtering use a value of 0.'};

  surfresamp      = cfg_exbranch;
  surfresamp.tag  = 'surfresamp';
  surfresamp.name = 'Resample and Smooth Surface Data';
  surfresamp.val  = {data_surf,merge_hemi,mesh32k,fwhm_surf,nproc};
  surfresamp.prog = @cat_surf_resamp;
  surfresamp.vout = @vout_surf_resamp;
  surfresamp.help = {
  'In order to analyze surface parameters all data have to be resampled into template space and the resampled data have to be finally smoothed. Resampling is done using the warped coordinates of the resp. sphere.'};




%% Resample and smooth FreeSurfer surfaces
%-----------------------------------------------------------------------
  data_fs         = cfg_files;
  data_fs.tag     = 'data_fs';
  data_fs.name    = 'Freesurfer Subject Directories';
  data_fs.filter  = 'dir';
  data_fs.ufilter = '.*';
  data_fs.num     = [1 Inf];
  data_fs.help    = {'Select subject folders of freesurfer data to resample data (e.g. thickness).'};

  measure_fs         = cfg_entry;
  measure_fs.tag     = 'measure_fs';
  measure_fs.name    = 'Freesurfer Measure';
  measure_fs.strtype = 's';
  measure_fs.num     = [1 Inf];
  measure_fs.val     = {'thickness'};
  measure_fs.help    = {
    'Name of surface measure that should be resampled and smoothed.'
    ''
    };

  outdir         = cfg_files;
  outdir.tag     = 'outdir';
  outdir.name    = 'Output Directory';
  outdir.filter  = 'dir';
  outdir.ufilter = '.*';
  outdir.num     = [0 1];
  outdir.help    = {'Select a directory where files are written.'};

  surfresamp_fs      = cfg_exbranch;
  surfresamp_fs.tag  = 'surfresamp_fs';
  surfresamp_fs.name = 'Resample and Smooth Existing FreeSurfer Surface Data';
  surfresamp_fs.val  = {data_fs,measure_fs,merge_hemi,mesh32k,fwhm_surf,outdir};
  surfresamp_fs.prog = @cat_surf_resamp_freesurfer;
  surfresamp_fs.help = {
  'If you have existing freesurfer data (e.g. thickness) this function can be used to resample these data, smooth the resampled data, and convert freesurfer data to gifti format.'};

%% Flipsides
%-----------------------------------------------------------------------
  flip.cdata         = cfg_files;
  flip.cdata.tag     = 'cdata';
  flip.cdata.name    = 'Surface Data Files';
  flip.cdata.filter  = 'gifti';
  flip.cdata.ufilter = '^s.*mm\.lh.*';
  flip.cdata.num     = [1 Inf];
  flip.cdata.help    = {'Texture maps that should be flipped/mirrored from right to left.'};
  

  flipsides      = cfg_exbranch;
  flipsides.tag  = 'flipsides';
  flipsides.name = 'Flip right to left hemisphere';
  flipsides.val  = {flip.cdata};
  flipsides.prog = @cat_surf_flipsides;
  flipsides.help = {
  'This function flip the right hemisphere to the left side, to allow side'
  ''
  ''};


%% Toolset
%-----------------------------------------------------------------------
  
  stools = cfg_choice;
  stools.name   = 'Surface Tools';
  stools.tag    = 'stools';
  if expert==2
    stools.values = { ...
      check_mesh_cov, ...
      surfextract, ...
      surfresamp, ...
      surfresamp_fs,...
      vol2surf, ...
      vol2tempsurf, ...
      surfcalc, ...
      surfcalcsub, ...
      surf2roi, ...
      roi2surf, ...
      ... estroi, ...
      flipsides, ...
      ... roicalc, ...
      };    
  elseif expert==1
    stools.values = { ...
      check_mesh_cov, ...
      surfextract, ...
      surfresamp, ...
      surfresamp_fs,...
      vol2surf, ...
      vol2tempsurf, ...
      surfcalc, ...
      surfcalcsub, ...
      surf2roi, ...
      ... roi2surf, ...
      ... estroi, ...
      flipsides, ...
      ... roicalc, ...
      };
  else
    stools.values = { ...
      check_mesh_cov, ...
      surfextract, ...
      surfresamp, ...
      surfresamp_fs,...
      surf2roi, ...
      vol2surf, ...
      vol2tempsurf, ...
      surfcalc, ...
      surfcalcsub, ...
      };
  end

%==========================================================================
function dep = vout_surf_surf2roi(job)

dep(1)            = cfg_dep;
dep(1).sname      = 'Extracted Surface ROIs';
dep(1).src_output = substruct('()',{1}, '.','xmlname','()',{':'});
dep(1).tgt_spec   = cfg_findspec({{'filter','xml','strtype','e'}});

%==========================================================================
function dep = vout_surf_resamp(job)

if job.merge_hemi
  dep(1)            = cfg_dep;
  dep(1).sname      = 'Merged Resample & Smooth';
  dep(1).src_output = substruct('()',{1}, '.','Psdata','()',{':'});
  dep(1).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
else
  dep(1)            = cfg_dep;
  dep(1).sname      = 'Left Resample & Smooth';
  dep(1).src_output = substruct('()',{1}, '.','lPsdata','()',{':'});
  dep(1).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
  dep(2)            = cfg_dep;
  dep(2).sname      = 'Right Resample & Smooth';
  dep(2).src_output = substruct('()',{1}, '.','rPsdata','()',{':'});
  dep(2).tgt_spec   = cfg_findspec({{'filter','gifti','strtype','e'}});
end

%==========================================================================
function dep = vout_surfextract(job)

if job.GI
  if ~exist('dep','var'), dep = cfg_dep;
  else, dep(end+1) = cfg_dep; end
  dep(end).sname      = 'Left gyrification';
  dep(end).src_output = substruct('()',{1}, '.','lPGI','()',{':'});
  dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
  dep(end+1)            = cfg_dep;
  dep(end).sname      = 'Right gyrification';
  dep(end).src_output = substruct('()',{1}, '.','rPGI','()',{':'});
  dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
end

if job.SD
  if ~exist('dep','var'), dep = cfg_dep;
  else, dep(end+1) = cfg_dep; end
  dep(end).sname      = 'Left sulcal depth';
  dep(end).src_output = substruct('()',{1}, '.','lPSD','()',{':'});
  dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
  dep(end+1)            = cfg_dep;
  dep(end).sname      = 'Right sulcal depth';
  dep(end).src_output = substruct('()',{1}, '.','rPSD','()',{':'});
  dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
end

if job.FD
  if ~exist('dep','var'), dep = cfg_dep;
  else, dep(end+1) = cfg_dep; end
  dep(end).sname      = 'Left fractal dimension';
  dep(end).src_output = substruct('()',{1}, '.','lPFD','()',{':'});
  dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
  dep(end+1)            = cfg_dep;
  dep(end).sname      = 'Right fractal dimension';
  dep(end).src_output = substruct('()',{1}, '.','rPFD','()',{':'});
  dep(end).tgt_spec   = cfg_findspec({{'filter','any','strtype','e'}});
end

%==========================================================================
function dep = vout_vol2surf(job)

if isfield(job,'merge_hemi') && job.merge_hemi
  dep(1)            = cfg_dep;
  dep(1).sname      = 'Mapped values';
  dep(1).src_output = substruct('.','mesh');
  dep(1).tgt_spec   = cfg_findspec({{'filter','mesh','strtype','e'}});
else
  dep(1)            = cfg_dep;
  dep(1).sname      = 'Left mapped values';
  dep(1).src_output = substruct('.','lh');
  dep(1).tgt_spec   = cfg_findspec({{'filter','mesh','strtype','e'}});
  dep(2)            = cfg_dep;
  dep(2).sname      = 'Right mapped values';
  dep(2).src_output = substruct('.','rh');
  dep(2).tgt_spec   = cfg_findspec({{'filter','mesh','strtype','e'}});
end
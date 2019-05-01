function PSgi = cat_surf_gyrification(type,PS,opt)
% Collection of gyrification measures
% ______________________________________________________________________
% This function includes surface relation measures to describe the 
% cortical folding by the folded original central surface and a second
% unfolded version of it.
%
% In the ideal case, the smoothing has to be done for each area and not
% the final measure, to be compatible to global or regional definitions.
% This is not standard yet! Therefore it is maybe better to combine this
% estimation with the resampling function and create only resampled 
% output data.
%
% WARNING:All measures required further test, validation and evaluation.
% 
%   cat_surf_SGI1(type,PS,opt)
% 
%   PS   .. individual central surface file (subject space)
%   opt  .. structure with different option (see code)
%   type .. 
%     'sphere':        central vs. sphere 
%                      > thin is anatomical usefull measure
%     'inflate':       central vs. inflate
%                      > this works, but depend on the kind and strength
%                        of the smooting 
%     'hullmapping':   central vs. hull  
%                      > only normalized output yet
%     'average':       central vs. fsavg 
%                      > only normalized output yet 
%     'laplacian':     central vs. hull by Laplacian approach 
%                      [Dahnke:2010,Jones:2001,Zilles:1989]
%
%  ---------------------------------------------------------------------
%
%  Sphere:
%   GI as the area-relation between the central surface and the sphere. 
%   This shows were regions have to be compressed in the registration 
%   process, but give no information about the anatomical relation. 
%
%
%  Inflate:
%   GI as the area-relation between the central surface and a smooth 
%   (areanormalized) version of it. High GI regions are the occupital 
%   and temporal lobe because these structures are overall much thinner 
%   than other parts of the brain. Subcortical structures on the other 
%   side show folding below 1. 
%   Different degree of inflateing are possible (opt.inflate = 0..10).
%   Furthermore a normalization by the global hull area is maybe usefull
%   to get a typical GI.
%
%
%  Hullmapping: 
%   GI as the arealelation between the central surface and hull. The 
%   hull is generated as separate surface and has to be maped to the
%   individual surface mesh. 
%   There are many options do to the spherical mapping. 
%   This is just for comparison to the Laplacian GI.
%
%
%  Average:
%   GI as the arealelation between the central surface and an average
%   surface. This will code the individual local volume increasement.
%   Different average surface are possible ...
%   This is my personal favorite because it describe the individual 
%   area decreasing/increasing compared to the group or healty control. 
%
%  Laplacian:
%   The function use the Laplacian mapping between two boundaries to 
%   deform the central surface to the position of the hull. 
%   Although, the laplacian GI is the best local representation of the 
%   classical GI [Zilles:1989], it has a lot of limitations: 
%   * It required strong smoothing 
%   * Due to the hull focus it is a sulcation measures, that is inverse 
%     to the development of gyris. 
%
% ______________________________________________________________________
% Robert Dahnke 
% Structural Brain Mapping Group
% University Jena
% ______________________________________________________________________
% $Id: cat_surf_gyrification.m 1070 2016-11-08 14:41:09Z dahnke $ 

  sinfo = cat_surf_info(PS);
  
  % set default options
  if ~exist('opt','var'), opt = struct(); end 
  def.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  opt = cat_io_checkinopt(opt,def); clear def; 
  
  def.debug     = 0;  % debuging information
  def.trerr     = 0;  % through error
  def.cleanup   = 1;  % delete temporary surfaces
  def.smooth    = 0;  % smoothing of original data (only for tests)
  
  def.presmooth = 5;  % inflate only: smooth area before GI estimation
  def.inflate   = 5;  % inflate only: only cat_surf_SGI_inflate with value from 1 to 10  
  def.normalize = 2;  % inflate only: normalization: 0 - none, 1 - by smooth, 2 - by smoothed and hull

  def.laplaceerr     = 0.0001; 
  def.GIpresmooth    = 30;    % laplacian only: shoold to be greater than 20 .. vary by brain size? - no, because of surf resamp 
  def.GIstreamopt(1) = 0.01; % stepsize of streamlines in mm... fast=0.05 - 0.01  
  def.GIstreamopt(2) = 30 ./ def.GIstreamopt(1); % in mm
  def.GIwritehull    = 1;  % write laplace laplacian hull surface
  
  def.avgsurf   = fullfile(opt.fsavgDir,sprintf('%s.central.freesurfer.gii',sinfo.side));      % fsaverage central
  def.Pfsavg    = fullfile(opt.fsavgDir,sprintf('%s.central.freesurfer.gii',sinfo.side));      % fsaverage central
  def.Pfsavgsph = fullfile(opt.fsavgDir,sprintf('%s.sphere.freesurfer.gii',sinfo.side));       % fsaverage sphere    
  
  opt = cat_io_checkinopt(opt,def);
   
  % call different subroutines
  switch lower(type)
    case 'sphere',      PSgi = cat_surf_SGI_sphere(sinfo,opt);
    case 'inflate',     PSgi = cat_surf_SGI_inflate(sinfo,opt);
    case 'hullmapping', PSgi = cat_surf_SGI_hullmapping(sinfo,opt);
    case 'average',     PSgi = cat_surf_SGI_average(sinfo,opt);
    case 'laplacian',   PSgi = cat_surf_SGI_laplacian(sinfo,opt);
    otherwise
      error('cat_surf_gyrification:unknown_type','Unknown type "%s".\n"',type);
  end
end

function Psgi = cat_surf_SGI_sphere(sinfo,opt) 
%% central vs. sphere
%  ---------------------------------------------------------------------
%  GI as the area-relation between the central surface and the sphere. 
%  This shows were regions have to be compressed in the registration 
%  process, but gives no information about the anatomical relation. 
%  ---------------------------------------------------------------------
  Psgi  = char(cat_surf_rename(sinfo,'resampled',0,'dataname','SGI','ee',''));

  Scs = gifti(sinfo.Pmesh);
  Ssp = gifti(sinfo.Psphere);

  % estiamte and smooth areas
  if opt.presmooth>0
    [ASsc,ASsp] = cat_surf_smoothGIarea(sinfo,Scs,Ssp,opt.presmooth);
  end
  
  % estimate GI
  GI = cat_surf_estimateGI(Scs,ASsc,ASsp,opt.normalize);
  
  % write data
  cat_io_FreeSurfer('write_surf_data',Psgi,GI);

end                      % central vs. sphere (work, but useless)

function Psgi = cat_surf_SGI_inflate(sinfo,opt)
%% central vs. inflate
%  ---------------------------------------------------------------------
%  GI as the area-relation between the central surface and a smooth 
%  version of it. High GI regions are the occipital and temporal lobe 
%  because these structures are overall much thinner than other parts 
%  of the brain. Subcortical structures on the other side show folding
%  below 1. 
%  ---------------------------------------------------------------------

  %Psgi      = char(cat_surf_rename(sinfo,'dataname',sprintf('IGI%d',opt.inflate),'ee',''));
  Psgi      = char(cat_surf_rename(sinfo,'dataname','IGI','ee',''));
  Pinflate  = char(cat_surf_rename(sinfo,'dataname','inflate'));
  
  % spherical mapping of the hull
  cmd = sprintf('CAT_Surf2Sphere "%s" "%s" %d',sinfo.Pmesh,Pinflate,opt.inflate);
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);

  % load surfaces
  Scs = gifti(sinfo.Pmesh);
  Ssp = gifti(Pinflate);

  % area smoothing and GI estimation
  ASsc = cat_surf_smootharea(Scs,Scs,opt.presmooth);
  ASsp = cat_surf_smootharea(Scs,Ssp,opt.presmooth);
  
  % estimate GI - normalization by area is not required (done by Surf2Sphere) 
  GI = cat_surf_estimateGI(Scs,ASsc,ASsp,opt.normalize);
  
  % write data
  cat_io_FreeSurfer('write_surf_data',Psgi,GI);
  
  % smoothing 
  if opt.smooth>0
    cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',sinfo.Pmesh,Psgi,opt.smooth,Psgi);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
    cat_surf_display(Psgi); 
  end
  
  % delete temporary files
  if opt.cleanup
    delete(Pinflate);
  end
end                      % central vs. inflate (work)

function Psgigii = cat_surf_SGI_average(sinfo,opt)
%% central vs. fs average
%  ---------------------------------------------------------------------
%  GI as the area-relation between the central surface and an average
%  surface. This will code the individual local volume increasement.
%  ---------------------------------------------------------------------
%  does not work yet ... only template space

  Psgi     = char(cat_surf_rename(sinfo,'resampled',0,'dataname','AGI','ee',''));
  Ptmp     = char(cat_surf_rename(sinfo,'resampled',0,'dataname','tmp','ee',''));
  Psgigii  = char(cat_surf_rename(sinfo,'resampled',0,'dataname','AGI'));
  Pcentral = char(cat_surf_rename(sinfo,'resampled',1));

  % resample hull surface in subject space
  cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s"',sinfo.Pmesh,sinfo.Psphere,opt.Pfsavgsph,Pcentral);
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug,opt.trerr); 
  
  % load surfaces
  Scs = gifti(Pcentral);
  Ssp = gifti(opt.avgsurf);

  % area smoothing and GI estimation
  ASsc = cat_surf_smootharea(Scs,Scs,opt.presmooth);
  ASsp = cat_surf_smootharea(Scs,Ssp,opt.presmooth);
  
  % estimate GI
  GI = cat_surf_estimateGI(Scs,ASsc,ASsp,opt.normalize);
  
  %% write data
  SR  = gifti(Pcentral); 
  SR.cdata = GI;
  
  save(gifti(SR),Psgigii);
  %cat_io_FreeSurfer('write_surf_data',Psgi,GI);
  
  %% resample hull surface in subject space ... 
  %cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',opt.Pfsavg,opt.Pfsavgsph,sinfo.Pmesh,sinfo.Psphere,Psgi,Psgi);
  %[ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug,opt.trerr); 
  
  %% smoothing
  if opt.smooth>0
    %cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pcentral,Ptmp,opt.smooth,Psgi);
    %[ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
    %cat_surf_display(Psgi)
    
    %%
    cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pcentral,Ptmp,opt.smooth,Psgi);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
    cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Pcentral,Ptmp,[Ptmp '.gii']);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
    cat_surf_display([Ptmp '.gii'])
  end
  
 

end                      % central vs. fs average (work, but only resampled output)

function Psgi = cat_surf_SGI_hullmapping(sinfo,opt)
%% central vs. hull
%  ---------------------------------------------------------------------
%  GI as the area-relation between the central surface and hull. 
%  The hull is generated as separate surface and has to be mapped to the
%  individual surface mesh. 
%  ---------------------------------------------------------------------

  Phull          = char(cat_surf_rename(sinfo,'dataname','hull'));
  Phullsphere    = char(cat_surf_rename(sinfo,'dataname','hullsphere'));
  Phullspherereg = char(cat_surf_rename(sinfo,'dataname','hullspherereg'));
  PhullR         = char(cat_surf_rename(sinfo,'dataname','centralhull'));
  Psgi           = char(cat_surf_rename(sinfo,'dataname','SGI','ee',''));
  
  % load surface and create hull
  Scs = gifti(sinfo.Pmesh);
  
  % create hull 
  Sh = cat_surf_fun('hull',Scs);
  
  % correct and optimze hull surface
  save(gifti(Sh),Phull);
  
  % remove some unconnected meshes
  cmd = sprintf('CAT_SeparatePolygon "%s" "%s" -1',Phull,Phull); 
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);

  % surface refinement by simple smoothing
  cmd = sprintf('CAT_BlurSurfHK "%s" "%s" %0.2f',Phull,Phull,5);
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);

  % spherical mapping of the hull
  cmd = sprintf('CAT_Surf2Sphere "%s" "%s" 5',Phull,Phullsphere);
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
  Phullsphere  = Phull;
  
  % spherical registration to central surface 
  cmd = sprintf('CAT_WarpSurf -norot -type 0 -i "%s" -is "%s" -t "%s" -ts "%s" -ws "%s"',...
    Phull,Phullsphere,sinfo.Pmesh,sinfo.Psphere,Phullspherereg);
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);

  % resample hull surface in subject space
  cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s"',Phull,Phullspherereg,sinfo.Psphere,PhullR);
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug,opt.trerr); 
  
  %% load surfaces
  Scs = gifti(sinfo.Pmesh);
  Ssp = gifti(PhullR);

  % area smoothing and GI estimation
  ASsc = cat_surf_smootharea(Scs,Scs,opt.presmooth);
  ASsp = cat_surf_smootharea(Scs,Ssp,opt.presmooth);
  
  % estimate GI
  GI = ASsc ./ ASsp;
  
  % write data
  cat_io_FreeSurfer('write_surf_data',Psgi,GI);
  
  % smoothing 
  if opt.smooth>0
    cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',sinfo.Pmesh,Psgi,opt.smooth,Psgi);
    [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
    cat_surf_display(Psgi)
  end

  if opt.cleanup
    delete(Phull);
    %delete(Phullsphere);
    delete(Phullspherereg);
    delete(PhullR); 
  end
end                  % central vs. hull (work, but only resampled output)

function Psgi = cat_surf_SGI_laplacian(sinfo,opt)
%% central vs. hull
%  ---------------------------------------------------------------------
%  The function use the Laplacian mapping between two boundaries to 
%  deform the central surface to the position of the hull. 
%  Although, the laplacian GI is the best local representation of the 
%  classical GI [Zilles:1989], it has a lot of limitations: 
%   * It required strong smoothing 
%   * Due to the hull focus it is a sulcation measures, that is inverse 
%     to the development of gyris. 
%  ---------------------------------------------------------------------

  Phull      = char(cat_surf_rename(sinfo,'dataname','hull','ee',''));
  Psgi       = char(cat_surf_rename(sinfo,'dataname','LGI','ee',''));
  
  % load surface and create hull
  Scs = gifti(sinfo.Pmesh);
  
  % create hull 
  [Sh,Yh] = cat_surf_fun('hull',Scs);
  [Yc,mat1] = cat_surf_fun('surf2vol',Scs);
  
  % estimate mapping
  opt.streamopt = opt.GIstreamopt; 
  S.vertices = Scs.vertices + repmat(mat1,size(Scs.vertices,1),1); 
  S.faces    = Scs.faces; 
  [S,SH]     = cat_surf_GI3D(S,1+Yh+Yc,Yh-Yc,opt);
  
  % write surface
  if opt.GIwritehull
    cat_io_FreeSurfer('write_surf',Phull,SH);
  end
 
  % area smoothing and GI estimation
  Ac = cat_surf_smootharea(S,S,opt.GIpresmooth);
  Ah = cat_surf_smootharea(S,SH,opt.GIpresmooth,cat_surf_smootharea(S,S,3)*2);
  
  GI  = Ac ./ max(eps,Ah); 
  % GI  = Ah ./ Ac; % inverse GI 
  % hist(GI(GI<4 & GI>0.01),0:0.1:4)
  cat_io_FreeSurfer('write_surf_data',Psgi,GI);
  
  % no smoothing here!
  if opt.smooth>0
    cat_surf_display(Psgi)
  end
  
end                    % central vs. hull (work)

function A = cat_surf_smootharea(S,SA,smooth,Amax)
%  create smooth area texture files
%  ---------------------------------------------------------------------
  debug = 0;

  % temporary file names
  Pname  = tempname; 
  Pmesh  = [Pname 'mesh'];
  Parea  = [Pname 'area'];

  % estimate areas 
  SAX.vertices = SA.vertices; SAX.faces = SA.faces; 
  A = cat_surf_fun('area',SAX);
  if exist('Amax','var');  A = min(A,Amax); end 
  
  % write surface and textures
  cat_io_FreeSurfer('write_surf',Pmesh,S);
  cat_io_FreeSurfer('write_surf_data',Parea,A);
  
  % smooth textures
  cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s"',Pmesh,Parea,smooth,Parea);
  [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,debug);
  
  % load smoothed textures
  A  = cat_io_FreeSurfer('read_surf_data',Parea);
  
  % delete temporary file
  delete(Parea);
end

function GI = cat_surf_estimateGI(Scs,ASsc,ASsp,normalize)
  % estimate GI - normalization by area is not required (done by Surf2Sphere) 
  if normalize==0 % no normalization
    GI = ASsc ./ ASsp;
  elseif normalize==1 % normalization to 1
    GI = (ASsc/sum(ASsc)) ./ (ASsp/sum(ASsp));
    %GI = (ASsc ./ ASsp) ./ (sum(ASsc)/sum(ASsp));
  elseif normalize==2 % normalization to hullarea
    Sh  = cat_surf_fun('hull',Scs);
    ASh = cat_surf_fun('area',Sh);
    GI  = ASsc ./ ASsp .* ((sum(ASsp) ./ sum(ASsc)) .* (sum(ASsc) ./ sum(ASh)));
  end
end
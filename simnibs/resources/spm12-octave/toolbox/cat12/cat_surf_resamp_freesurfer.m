function cat_surf_resamp_freesurfer(vargin)
%cat_surf_resamp_freesurfer to resample parameters to template
% space and smooth it.
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_surf_resamp_freesurfer.m 1264 2018-01-30 11:02:00Z gaser $
  
  rev = '$Rev: 1264 $';
  
  if nargin == 1
    Psubj = char(vargin.data_fs);
    fwhm_surf = vargin.fwhm_surf;
    outdir = vargin.outdir{1};
    mesh32k = vargin.mesh32k;
    pname = vargin.measure_fs;
  else
    error('Not enough parameters.');
  end
  
  opt.debug     = cat_get_defaults('extopts.verb') > 2;
      
  if mesh32k
    opt.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k');
    str_resamp = '.resampled_32k';
  else
    opt.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
    str_resamp = '.resampled';
  end

  hemi_str = char('lh','rh');
  
  for i=1:size(Psubj,1)
  
    stime = clock; 
    exist_hemi = [];
    [pp,name]   = spm_fileparts(deblank(Psubj(i,:)));
    
    % subject directory
    dname = fullfile(pp,name,'surf');
  
    % check that surf subfolder exists
    if ~exist(dname,'dir')
      fprintf('Could not find ''surf'' subfolder in %s.\n\n',Psubj(i,:));
      continue
    end
    
    for j=1:2
    
      hemi = hemi_str(j,:);
      exist_hemi = [exist_hemi j];
      
      Psmoothwm  = fullfile(dname,[hemi '.smoothwm']);
      Psphere    = fullfile(dname,[hemi '.sphere']);
      Pspherereg = fullfile(dname,[hemi '.sphere.reg']);
      Pmeasure   = fullfile(dname,[hemi '.' pname]);
      Presamp    = fullfile(dname,[hemi '.smoothwm' str_resamp]);
      Pvalue     = fullfile(dname,[hemi '.' pname str_resamp]);
      
      if fwhm_surf > 0
          Pfwhm      = fullfile(outdir,[sprintf('s%g.',fwhm_surf) hemi '.' pname str_resamp '.'  name]);
      else
          Pfwhm      = fullfile(outdir,[hemi '.' pname str_resamp '.'  name]);
      end
  
      % save fwhm name to merge meshes
      Pfwhm_all{j} = [Pfwhm '.gii'];
  
      Pfsavg     = fullfile(opt.fsavgDir,[hemi '.sphere.freesurfer.gii']);
      Pmask      = fullfile(opt.fsavgDir,[hemi '.mask']);
    
      fprintf('Resample %s in %s\n',hemi,deblank(Psubj(i,:)));
  
      % resample values using warped sphere 
      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Psmoothwm,Pspherereg,Pfsavg,Presamp,Pmeasure,Pvalue);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
    
      % smooth resampled values
      cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s" "%s"',Presamp,Pfwhm,fwhm_surf,Pvalue,Pmask);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
  
      % add values to resampled surf and save as gifti
      cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Pfwhm,[Pfwhm '.gii']);
      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,opt.debug);
    
      % remove path from metadata to allow that files can be moved (pathname is fixed in metadata) 
      [pp2,ff2,ex2]   = spm_fileparts([Pfwhm '.gii']);
      g = gifti([Pfwhm '.gii']);
      g.private.metadata = struct('name','SurfaceID','value',[ff2 ex2]);
      save(g, [Pfwhm '.gii'], 'Base64Binary');
  
      delete(Presamp);
      delete(Pfwhm);
      if fwhm_surf > 0, delete(Pvalue); end
   end
   
    % merge hemispheres
    if vargin.merge_hemi
      % replace hemi info with "mesh"    
      Pfwhm   = strrep(Pfwhm_all{1},['lh.' pname],['mesh.' pname]);
      [pp,ff,ex]   = spm_fileparts(Pfwhm);
  
      % combine left and right and optionally cerebellar meshes
      if numel(exist_hemi) > 1
        M0 = gifti({Pfwhm_all{1}, Pfwhm_all{2}});
        delete(Pfwhm_all{1}); delete(Pfwhm_all{2})
        M.faces = [M0(1).faces; M0(2).faces+size(M0(1).vertices,1)];
        M.vertices = [M0(1).vertices; M0(2).vertices];
        M.cdata = [M0(1).cdata; M0(2).cdata];
        M.mat = M0(1).mat;
        M.private.metadata = struct('name','SurfaceID','value',[ff ex]);
        save(gifti(M), Pfwhm, 'Base64Binary');
        Psdata{i} = Pfwhm;
      else
        disp('No data for opposite hemisphere found!');
      end
              
      fprintf('(%3.0f s) Display resampled %s\n',etime(clock,stime),spm_file(Psdata{i},'link','cat_surf_display(''%s'')'));
    end
  
  end

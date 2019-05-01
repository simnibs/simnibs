function varargout = cat_surf_resamp(varargin)
% ______________________________________________________________________
% Function to resample parameters to template space and smooth it.
%
% [Psdata] = cat_surf_resamp(job)
% 
% job.data_surf .. cellstr of files
% job.fwhm_surf .. filter size in mm
% job.verb      .. display command line progress
% ______________________________________________________________________
% Christian Gaser
% $Id: cat_surf_resamp.m 1251 2017-12-29 09:22:44Z gaser $


  SVNid = '$Rev: 1251 $';

  if nargin == 1
    P    = char(varargin{1}.data_surf);
    job  = varargin{1}; 
  else
    spm_clf('Interactive'); 
    P = cellstr(spm_select([1 inf],'any','Select surface data'));
    job = struct();
  end

  if ~isfield(job,'fwhm_surf')
    spm('alert!', ['Surface smoothing method has changed with release r1248 and the '...
    'recommended FWHM is now slightly smaller. For cortical thickness a good starting value '...
    'is 15mm, while other surface parameters based on cortex folding (e.g. gyrification, '...
    'cortical complexity) need a larger filter size of about 20-25mm. Please update your scripts '...
    'and replace the old field "fwhm" by "fwhm_surf" and adapt the values.'], 1);  
  end
  
  def.trerr      = 0; 
  def.fwhm_surf  = 0; 
  def.nproc      = 0; 
  def.mesh32k    = 1; 
  def.merge_hemi = 1;
  def.verb       = cat_get_defaults('extopts.verb'); 
  def.debug      = cat_get_defaults('extopts.verb')>2;
  def.fsavgDir   = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
  
  job = cat_io_checkinopt(job,def);

  if job.mesh32k
    job.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k');
    str_resamp = '.resampled_32k';
  else
    job.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
    str_resamp = '.resampled';
  end

  % split job and data into separate processes to save computation time
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
    if nargout==1
      varargout{1} = cat_parallelize(job,mfilename,'data_surf');
    else
      cat_parallelize(job,mfilename,'data_surf');
    end
    return
  end  
  
  % normal processing
  % ____________________________________________________________________
  
  % new banner
  if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end
  
  % display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',size(P,1),'Smoothed Resampled','Surfaces Completed');

  Psdata  = cell(size(P,1),1);
  lPsdata = cell(size(P,1),1);
  rPsdata = cell(size(P,1),1);
  
  for i=1:size(P,1)
    
    stime = clock; 
    [pp,ff,ex]   = spm_fileparts(deblank(P(i,:)));
    if any([strfind(ff,'.sphere.'),strfind(ff,'.central.')])
      if job.verb
        fprintf('Cannot process "%s"!\n',deblank(P(i,:)));
      end
      continue; 
    end
    
    name0 = [ff(3:end) ex];          % remove leading hemisphere information
    name0 = strrep(name0,'.gii',''); % remove .gii extension
    hemi = ff(1:2);
    hemistr = char('lh','rh','lc','rc');
    exist_hemi = [];
    
    % go through left and right and potential cerebellar hemispheres
    for j=1:length(hemistr)
    
      % add hemisphere name
      hemi = hemistr(j,:);
      name = [hemi name0];
      
      Pvalue0 = fullfile(pp,name);
      
      % check that file exists
      if ~exist(Pvalue0,'file'), continue; end
      
      exist_hemi = [exist_hemi j];

      k = strfind(name,'.');
      pname = ff(k(1)+1:k(2)-1);
      Pcentral   = [strrep(name,pname,'central') '.gii'];
      Pspherereg = fullfile(pp,strrep(Pcentral,'central','sphere.reg'));
      Pvalue     = fullfile(pp,strrep(Pcentral,'central',[pname str_resamp]));
      Pvalue     = strrep(Pvalue,'.gii',''); % remove .gii extension
      
      if job.fwhm_surf > 0
        Pfwhm    = fullfile(pp,[sprintf('s%g.',job.fwhm_surf) strrep(Pcentral,'central',[pname str_resamp])]);
        Presamp  = fullfile(pp,[sprintf('s%g.',job.fwhm_surf) strrep(Pcentral,'central',[pname '.tmp.resampled'])]);
      else
        Pfwhm    = fullfile(pp,strrep(Pcentral,'central',[pname str_resamp]));
        Presamp  = fullfile(pp,strrep(Pcentral,'central',[pname 'tmp.resampled']));
      end
      
      Pfwhm      = strrep(Pfwhm,'.gii',''); % remove .gii extension
      Pcentral   = fullfile(pp,Pcentral);
      Pfsavg     = fullfile(job.fsavgDir,[hemi '.sphere.freesurfer.gii']);
      Pmask      = fullfile(job.fsavgDir,[hemi '.mask']);
      
      % we have to rename the final files for each hemisphere if we want to merge the hemispheres 
      % to not interfere with existing files
      if job.merge_hemi
        Pfwhm_gii = [Pfwhm '_tmp.gii'];
      else
        Pfwhm_gii = [Pfwhm '.gii'];
      end
      
      % save fwhm name to merge meshes
      Pfwhm_all{j} = Pfwhm_gii;
      
      % resample values using warped sphere 
      cmd = sprintf('CAT_ResampleSurf "%s" "%s" "%s" "%s" "%s" "%s"',Pcentral,Pspherereg,Pfsavg,Presamp,Pvalue0,Pvalue);
      [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end

      % resample surface using warped sphere with better surface quality (using Spherical harmonics)
      cmd = sprintf('CAT_ResampleSphericalSurfSPH -n 327680 "%s" "%s" "%s"',Pcentral,Pspherereg,Presamp);
      [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end

      % resample surface according to freesurfer sphere
      cmd = sprintf('CAT_ResampleSurf "%s" NULL "%s" "%s"',Presamp,Pfsavg,Presamp);
      [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end

      % smooth resampled values
      cmd = sprintf('CAT_BlurSurfHK "%s" "%s" "%g" "%s" "%s"',Presamp,Pfwhm,job.fwhm_surf,Pvalue,Pmask);
      [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end

      % add values to resampled surf and save as gifti
      cmd = sprintf('CAT_AddValuesToSurf "%s" "%s" "%s"',Presamp,Pfwhm,Pfwhm_gii);
      [ST, RS] = cat_system(cmd); err = cat_check_system_output(ST,RS,job.debug,def.trerr); if err, continue; end

      if exist(Pfwhm_gii,'file'), Psname = Pfwhm_gii; end
      
      % remove path from metadata to allow that files can be moved (pathname is fixed in metadata) 
      [pp2,ff2,ex2]   = spm_fileparts(Psname);

      g = gifti(Psname);
      g.private.metadata = struct('name','SurfaceID','value',[ff2 ex2]);
      save(g, Psname, 'Base64Binary');

      delete(Presamp);
      delete(Pfwhm);
      if job.fwhm_surf > 0, delete(Pvalue); end

      if job.verb
        fprintf('Resampling %s\n',Psname);
      end
      
      if j==1, lPsdata{i} = Psname; end
      if j==2, rPsdata{i} = Psname; end
    end

    % merge hemispheres
    if job.merge_hemi
      % name for combined hemispheres
      k = strfind(name,'.');
      pname = ff(k(1)+1:k(2)-1);
      Pcentral   = strrep(['mesh' name0 '.gii'],pname,'central');
      
      if job.fwhm_surf > 0
        Pfwhm     = [sprintf('s%g.',job.fwhm_surf) strrep(Pcentral,'central',[pname str_resamp])];
      else
        Pfwhm     = strrep(Pcentral,'central',[pname str_resamp]);
      end
  
      % combine left and right and optionally cerebellar meshes
      switch numel(exist_hemi)
      case{2,4}
        M0 = gifti({Pfwhm_all{1}, Pfwhm_all{2}});
        delete(Pfwhm_all{1}); delete(Pfwhm_all{2})
        M.faces = [M0(1).faces; M0(2).faces+size(M0(1).vertices,1)];
        M.vertices = [M0(1).vertices; M0(2).vertices];
        M.cdata = [M0(1).cdata; M0(2).cdata];
      case 4
        M0 = gifti({Pfwhm_all{3}, Pfwhm_all{4}});
        delete(Pfwhm_all{3}); delete(Pfwhm_all{4})
        M.faces = [M.faces; M0(1).faces+2*size(M0(1).vertices,1); M0(2).faces+3*size(M0(1).vertices,1)];
        M.vertices = [M.vertices; M0(1).vertices; M0(2).vertices];
        M.cdata = [M.cdata; M0(1).cdata; M0(2).cdata];
      case 1
        disp('No data for opposite hemisphere found!');
      case 3
        disp('No data for opposite cerebellar hemisphere found!');
      end
      
      if numel(exist_hemi) > 1
        M.mat = M0(1).mat;
        M.private.metadata = struct('name','SurfaceID','value',Pfwhm);
        save(gifti(M), fullfile(pp,Pfwhm), 'Base64Binary');
        Psdata{i} = fullfile(pp,Pfwhm);
      end
          
      if job.verb
        fprintf('(%3.0f s) Display resampled %s\n',etime(clock,stime),spm_file(Psdata{i},'link','cat_surf_display(''%s'')'));
      end
    end
    
    spm_progress_bar('Set',i);
  end
  
  if isfield(job,'process_index')
    fprintf('Done\n'); 
  end
      
  if nargout==1
    if job.merge_hemi
      varargout{1}.Psdata = Psdata; 
    else
      varargout{1}.lPsdata = lPsdata; 
      varargout{1}.rPsdata = rPsdata; 
    end
  end
  
  spm_progress_bar('Clear');
end

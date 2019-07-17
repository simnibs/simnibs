function varargout = cat_surf_parameters(job)
% ______________________________________________________________________
%
% cat_surf_parameters to extract surface parameters such as
% gyrification and cortical complexity.
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_surf_parameters.m 1233 2017-12-03 23:26:25Z gaser $

  SVNid = '$Rev: 1233 $';
 
  if nargin == 1
    P  = char(job.data_surf);
  else
    error('Not enough parameters.');
  end
  
  def.trerr       = 0; 
  def.nproc       = 0; 
  def.verb        = cat_get_defaults('extopts.verb'); 
  def.lazy        = 0; % reprocess exist results
  def.debug       = cat_get_defaults('extopts.verb')>2;
  
  job = cat_io_checkinopt(job,def);

  % split job and data into separate processes to save computation time
  if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
    if nargout==1
      varargout{1} = cat_parallelize(job,mfilename,'data_surf');
    else
      cat_parallelize(job,mfilename,'data_surf');
    end
    return
  end

  % new banner
  if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end
  
  % display something
  spm_clf('Interactive'); 
  spm_progress_bar('Init',size(P,1),'Processed surfaces','Surfaces Completed');
  
  for i=1:size(P,1)
  
    % go through left and right hemisphere
    for j=1:2
    
      if j == 1 % lh
        Pname = deblank(P(i,:));
      else % rh
        Pname = cat_surf_rename(deblank(P(i,:)),'side','rh');
        Pname = Pname{1};
      end
      
      [pp,ff,ex]   = spm_fileparts(Pname);
  
      name = [ff ex];
  
      PGI     = fullfile(pp,strrep(ff,'central','gyrification'));         
      PGII    = fullfile(pp,strrep(ff,'central','IGI'));         
      PGIA    = fullfile(pp,[strrep(ff,'central','AGI'),'.gii']);       
      PGIS    = fullfile(pp,strrep(ff,'central','SGI'));         
      PGIL    = fullfile(pp,strrep(ff,'central','LGI'));         
      PFD     = fullfile(pp,strrep(ff,'central','fractaldimension'));
      PSD     = fullfile(pp,strrep(ff,'central','sqrtsulc'));
      PSA     = fullfile(pp,strrep(ff,'central','logarea'));
      PIS     = fullfile(pp,strrep(ff,'central','inner'));         
      POS     = fullfile(pp,strrep(ff,'central','outer'));         
      Psphere = fullfile(pp,strrep(name,'central','sphere'));
   
      fprintf('Extract parameters for %s\n',Pname);
      if job.GI 
        %% gyrification index based on absolute mean curvature
        if exist(PGI,'file') && job.lazy  
          if job.verb>=1, fprintf('  Display already processed %s\n',spm_file(PGI,'link','cat_surf_display(''%s'')')); end
        else
          cmd = sprintf('CAT_DumpCurv "%s" "%s" 0 0 1',Pname,PGI);
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
          if nargout==1 && j==1, varargout{1}.lPGI{i} = PGI; end  
          if nargout==1 && j==2, varargout{1}.rPGI{i} = PGI; end  
          if job.verb>=1, fprintf('  Display %s\n',spm_file(PGI,'link','cat_surf_display(''%s'')')); end
        end
      end
      
      
      roijob.cdata  = P; 
  %    roijob.rdata  = job.rdata; % works only for developer mode
      roijob.verb   = job.verb;
      roijob.avg    = struct('mean',1,'std',0,'min',0,'max',0,'median',0);  
      roijob.area   = 0; % separate function 
      roijob.vernum = 0; 
      
      
      % expert folding measures
      % ---------------------------------------------------------------------
      % These approaches are still in development. 
      % See cat_surf_gyrification for further information.
      % ---------------------------------------------------------------------
      if cat_get_defaults('extopts.expertgui') > 1
        if job.GIA
          %% gyrification index based on average surface
          %  basic idea is to find local increase/decreasment of surface area
          %  (and compare it to thickness) or estimate the local volume etc.
          %  > I like this idea, but it required more work ... 
          if exist(PGIA,'file') && job.lazy  
            if job.verb>=1, fprintf('  Display already processed %s\n',spm_file(PGIA,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            PGIA  = cat_surf_gyrification('average',Pname);
            if nargout==1 && j==1, varargout{1}.lPGIA{i} = PGIA; end  
            if nargout==1 && j==2, varargout{1}.rPGIA{i} = PGIA; end  
            if job.verb>=1, fprintf('  %4.0fs. Display %s\n',etime(clock,stime),spm_file(PGIA,'link','cat_surf_display(''%s'')')); end
          end
          
        end
  
        if job.GII
          %% gyrification index based on inflating
          %  simple GI approach with differnt smoothing options
          if exist(PGII,'file') && job.lazy  
            if job.verb>=1, fprintf('  Display already processed %s\n',spm_file(PGII,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            PGII = cat_surf_gyrification('inflate',Pname,struct('inflate',5));
            if nargout==1 && j==1, varargout{1}.lPGII{i} = PGII; end
            if nargout==1 && j==2, varargout{1}.rPGII{i} = PGII; end
            if job.verb>=1, fprintf('  %4.0fs. Display %s\n',etime(clock,stime),spm_file(PGII,'link','cat_surf_display(''%s'')')); end
          end
         end
  
         if job.GIL
          %% gyrification index based on laplacian GI
          %  nice model, but very close to sulcal depth
          if exist(PGIL,'file') && job.lazy  
            if job.verb>=1, fprintf('  Display already processed %s\n',spm_file(PGIL,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            PGIL = cat_surf_gyrification('laplacian',Pname,struct('verb',job.verb));
            if nargout==1 && j==1, varargout{1}.lPGIL{i} = PGIL; end
            if nargout==1 && j==2, varargout{1}.rGIL{i} = PGIL; end
            if job.verb>=1, fprintf('  %4.0fs. Display %s\n',etime(clock,stime),spm_file(PGIL,'link','cat_surf_display(''%s'')')); end
          end
        end
  
        if job.GIS
          %% gyrification index based on spericial mapping with hull
          %  similar to GII, but with stronger normalization by the sphereical
          %  mapping
          if exist(PGIS,'file') && job.lazy  
            if job.verb>=1, fprintf('  Display already processed %s\n',spm_file(PGIS,'link','cat_surf_display(''%s'')')); end
          else
            stime = clock; 
            PGIS = cat_surf_gyrification('hullmapping',Pname);
            if nargout==1 && j==1, varargout{1}.lPGIS{i} = PGIS; end  
            if nargout==1 && j==2, varargout{1}.rPGIS{i} = PGIS; end  
            if job.verb>=1, fprintf('  %4.0fs. Display %s\n',etime(clock,stime),spm_file(PGIL,'link','cat_surf_display(''%s'')')); end
          end
        end
  
      %  if job.SA
      %    %% local surface area
      %    fprintf('Not yet working.');
      %    if exist(PSA,'file') && job.lazy  
      %      if job.verb>=1, fprintf('  Display already processed %s\n',spm_file(PSA,'link','cat_surf_display(''%s'')')); end
      %    else 
      %      cmd = sprintf('CAT_DumpSurfArea -log -sphere "%s" "%s" "%s"',Psphere,Pname,PSA);
      %      [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug);
      %      if nargout==1, varargout{1}.PSA{i} = PSA; end  
      %      if job.verb>=1, fprintf('  Display %s\n',spm_file(PSA,'link','cat_surf_display(''%s'')')); end
      %    end
      %  end
  
      end
    % ----------------------------------------------------------------------
    % ----------------------------------------------------------------------
      
      
        
      
      if job.SD
        %% sulcus depth
        if exist(PSD,'file') && job.lazy  
          if job.verb>=1, fprintf('  Display already processed %s\n',spm_file(PSD,'link','cat_surf_display(''%s'')')); end
        else
          cmd = sprintf('CAT_SulcusDepth -sqrt "%s" "%s" "%s"',Pname,Psphere,PSD); %-sqrt
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
          if nargout==1 && j==1, varargout{1}.lPSD{i} = PSD; end
          if nargout==1 && j==2, varargout{1}.rPSD{i} = PSD; end
          if job.verb>=1, fprintf('  Display %s\n',spm_file(PSD,'link','cat_surf_display(''%s'')')); end
        end
      end
  
      if job.FD
        %% fractal dimension using spherical harmonics
        if exist(PFD,'file') && job.lazy  
          if job.verb>=1, fprintf('  Display already processed %s\n',spm_file(PFD,'link','cat_surf_display(''%s'')')); end
        else
          cmd = sprintf('CAT_FractalDimension -sphere "%s" -nosmooth "%s" "%s" "%s"',Psphere,Pname,Psphere,PFD);
          [ST, RS] = cat_system(cmd); cat_check_system_output(ST,RS,job.debug,job.trerr);
          if nargout==1 && j==1, varargout{1}.lPFD{i} = PFD; end  
          if nargout==1 && j==2, varargout{1}.rPFD{i} = PFD; end  
          if job.verb>=1, fprintf('  Display %s\n',spm_file(PFD,'link','cat_surf_display(''%s'')')); end
        end
      end
  
      
      
      
      %% ----------------------------------------------------------------------
      %  No measures, but I do not want another script
      %  ----------------------------------------------------------------------
      if isfield('IS',job) && job.IS
        if exist(PIS,'file') && job.lazy  
          if job.verb>=1, fprintf('  Display already processed %s\n',spm_file(PIS,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          PIS = cat_surf_fun('inner',Pname);
          if nargout==1 && j==1, varargout{1}.lPIS{i} = PIS; end  
          if nargout==1 && j==2, varargout{1}.rPIS{i} = PIS; end  
          if job.verb>=1, fprintf('  %4.0fs. Display %s\n',etime(clock,stime),spm_file(PIS,'link','cat_surf_display(''%s'')')); end
        end
      end
      
      if isfield('OS',job) && job.OS
        if exist(POS,'file') && job.lazy  
          if job.verb>=1, fprintf('  Display already processed %s\n',spm_file(POS,'link','cat_surf_display(''%s'')')); end
        else
          stime = clock; 
          POS = cat_surf_fun('outer',Pname);
          if nargout==1 && j==1, varargout{1}.lPOS{i} = POS; end  
          if nargout==1 && j==2, varargout{1}.rPOS{i} = POS; end  
          if job.verb>=1, fprintf('  %4.0fs. Display %s\n',etime(clock,stime),spm_file(POS,'link','cat_surf_display(''%s'')')); end
        end
      end
        
      spm_progress_bar('Set',i);
    end
    
    if isfield(job,'process_index')
      fprintf('Done\n');
    end  
    
  end
  spm_progress_bar('Clear');  
  
end

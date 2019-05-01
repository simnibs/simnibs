function varargout = cat_parallelize(job,func,datafield)
% ______________________________________________________________________
% Function to parallelize other functions with job structure, by the 
% following call:
% 
%   SVNid = '$Rev: 1093 $';
% 
%   ... further initialization code
%  
%   % split job and data into separate processes to save computation time
%   if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))
%     if nargout==1
%       varargout{1} = cat_parallelize(job,mfilename,'data_surf');
%     else
%       cat_parallelize(job,mfilename,'data_surf');
%     end
%     return
%   elseif isfield(job,'printPID') && job.printPID 
%     cat_display_matlab_PID
%   end 
%  
%   % new banner
%   if isfield(job,'process_index'), spm('FnBanner',mfilename,SVNid); end
%   
%   % add system dependent extension to CAT folder
%   if ispc
%     job.CATDir = [job.CATDir '.w32'];
%   elseif ismac
%     job.CATDir = [job.CATDir '.maci64'];
%   elseif isunix
%     job.CATDir = [job.CATDir '.glnx86'];
%   end  
%
%   ... main code
% ______________________________________________________________________
% Christian Gaser, Robert Dahnke
% $Id: cat_parallelize.m 1093 2016-12-15 14:37:27Z dahnke $

  def.verb      = cat_get_defaults('extopts.verb'); 
  def.lazy      = 0; % reprocess exist results
  def.debug     = cat_get_defaults('extopts.verb')>2;
  job.CATDir    = fullfile(spm('dir'),'toolbox','cat12','CAT');   
  job.fsavgDir  = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces'); 
 
  job = cat_io_checkinopt(job,def);

  if ~exist('datafield','var'), datafield = 'data'; end

  cat_io_cprintf('warn',...
    ['\nWARNING: Please note that no additional modules in the batch can be run \n' ...
     '         except CAT12 segmentation. Any dependencies will be broken for \n' ...
     '         subsequent modules if you split the job into separate processes.\n\n']);

  % rescue original subjects
  job_data = job.(datafield);
  if iscell(job.(datafield){1})
    n_subjects = numel(job.(datafield){1});
  else
    n_subjects = numel(job.(datafield));
  end
  if job.nproc > n_subjects
    job.nproc = n_subjects;
  end
  job.process_index = cell(job.nproc,1);
  job.verb = 1; 

  % initial splitting of data
  for i=1:job.nproc
    job.process_index{i} = (1:job.nproc:(n_subjects-job.nproc+1))+(i-1);
  end

  % check if all data are covered
  for i=1:rem(n_subjects,job.nproc)
    job.process_index{i} = [job.process_index{i} n_subjects-i+1];
  end

  tmp_array = cell(job.nproc,1);
  logdate   = datestr(now,'YYYYmmdd_HHMMSS');
  for i=1:job.nproc
    fprintf('Running job %d (with datafield 1):\n',i);
    if iscell(job.(datafield){1})
      for fi=1:numel(job_data{1}(job.process_index{i}))
        fprintf('  %s\n',spm_str_manip(char(job_data{1}(job.process_index{i}(fi))),'a78')); 
      end
      for ix=1:numel(job_data)
        job.(datafield){ix} = job_data{ix}(job.process_index{i});
      end
    else
      for fi=1:numel(job_data(job.process_index{i}))
        fprintf('  %s\n',spm_str_manip(char(job_data(job.process_index{i}(fi))),'a78')); 
      end
      job.(datafield) = job_data(job.process_index{i});
    end
    job.verb        = 1; 
    job.printPID    = 1; 
    % temporary name for saving job information
    tmp_name = [tempname '.mat'];
    tmp_array{i} = tmp_name; 
    def = cat_get_defaults; job = cat_io_checkinopt(job,def); % further job update required here to get the latest cat defaults
    global defaults cat12; %#ok<NUSED,TLEV>
    save(tmp_name,'job','defaults','cat12');
    clear defaults cat12;

    % matlab command, cprintferror=1 for simple printing        
    matlab_cmd = sprintf('"global cprintferror; cprintferror=1; addpath %s %s; load %s; %s(job); "',...
      spm('dir'),fullfile(spm('dir'),'toolbox','cat12'),tmp_name,func);

    % log-file for output
    log_name = ['log_' func '_' logdate '_' sprintf('%02d',i) '.txt'];

    % call matlab with command in the background
    if ispc
      % check for spaces in filenames that will not work with windows systems and background jobs
      if strfind(spm('dir'),' ')
        cat_io_cprintf('warn',...
            ['\nWARNING: No background processes possible because your SPM installation is located in \n' ...
             '         a folder that contains spaces. Please set the number of processes in the GUI \n'...
             '         to ''0''. In order to split your job into different processes,\n' ...
             '         please do not use any spaces in folder names!\n\n']);
         job.nproc = 0;
         job = update_job(job);
         varargout{1} = run_job(job);
         return; 
      end
      % prepare system specific path for matlab
      export_cmd = ['set PATH=' fullfile(matlabroot,'bin')];
      [status,result] = system(export_cmd);
      system_cmd = ['start matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name];
    else
      % -nodisplay .. nodisplay is without figure output > problem with CAT report ... was there a server problem with -nodesktop?
      system_cmd = [fullfile(matlabroot,'bin') '/matlab -nodesktop -nosplash -r ' matlab_cmd ' -logfile ' log_name ' 2>&1 & '];
    end

    [status,result] = system(system_cmd);

    test = 0; lim = 10; ptime = 0.5;
    while test<lim
      if ~exist(log_name,'file')
        pause(ptime); 
        test = test + ptime; 
        if test>=lim
          cat_io_cprintf('warn',sprintf('"%s" did not exist after %d seconds! Proceed! \n',log_name,lim))
        end
      else 
        test = inf; 
        edit(log_name);
      end
    end

    fprintf('\nCheck %s for logging information.\n',spm_file(log_name,'link','edit(''%s'')'));
    fprintf('_______________________________________________________________\n');


  end

  
  %job = update_job(job);
  varargout{1} = job; 
  %vout_job(job);

end
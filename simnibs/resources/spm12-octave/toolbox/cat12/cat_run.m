function varargout = cat_run(job)
% Segment a bunch of images
% ______________________________________________________________________
%
%   FORMAT cat_run(job)
%
%   job.channel(n).vols{m}
%   job.channel(n).biasreg
%   job.channel(n).biasfwhm
%   job.channel(n).write
%   job.tissue(k).tpm
%   job.tissue(k).ngaus
%   job.tissue(k).native
%   job.tissue(k).warped
%
% See the user interface for a description of the fields.
%
% based on John Ashburners version of
% spm_preproc8_run.m 2281 2008-10-01 12:52:50Z john $
% ______________________________________________________________________
% Christian Gaser
% $Id: cat_run.m 1245 2017-12-12 10:35:00Z dahnke $

%#ok<*AGROW>

%rev = '$Rev: 1245 $';

%  -----------------------------------------------------------------
%  Lazy processing (expert feature)
%  -----------------------------------------------------------------
%  If N>10000 files were processed the crash of one of J jobs by 
%  small errors makes it hard to find the unprocess files. 
%  The lazy processing will only process files, if one of the output
%  is missed and if the same preprocessing options were used before.
%  -----------------------------------------------------------------
if isfield(job.extopts,'admin') && isfield(job.extopts.admin,'lazy') && job.extopts.admin.lazy && ...
  ~isfield(job,'process_index') && isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))  
  jobl      = update_job(job);
  jobl.vout = vout_job(jobl);
  job.data  = remove_already_processed(jobl); 
end



% split job and data into separate processes to save computation time
if isfield(job,'nproc') && job.nproc>0 && (~isfield(job,'process_index'))  
  cat_io_cprintf('warn',...
    ['\nWARNING: Please note that no additional modules in the batch can be run \n' ...
     '         except CAT12 segmentation. Any dependencies will be broken for \n' ...
     '         subsequent modules if you split the job into separate processes.\n\n']);
    
  % rescue original subjects
  job_data = job.data;
  n_subjects = numel(job.data);
  if job.nproc > n_subjects
    job.nproc = n_subjects;
  end
  job.process_index = cell(job.nproc,1);

  % initial splitting of data
  for i=1:job.nproc
    job.process_index{i} = (1:job.nproc:(n_subjects-job.nproc+1))+(i-1);
  end

  % check if all data are covered
  for i=1:rem(n_subjects,job.nproc)
    job.process_index{i} = [job.process_index{i} n_subjects-i+1];
  end

  tmp_array = cell(job.nproc,1); job.printPID = 1; 
    
  logdate   = datestr(now,'YYYYmmdd_HHMMSS');
  for i=1:job.nproc
    fprintf('Running job %d:\n',i);
    for fi=1:numel(job_data(job.process_index{i}))
      fprintf('  %s\n',spm_str_manip(char(job_data(job.process_index{i}(fi))),'a78')); 
    end
    job.data = job_data(job.process_index{i});
         
    % temporary name for saving job information
    tmp_name = [tempname '.mat'];
    tmp_array{i} = tmp_name; 
    %def = cat_get_defaults; job = cat_io_checkinopt(job,def); % further job update required here to get the latest cat defaults
    spm12def = spm_get_defaults;    
    cat12def = cat_get_defaults; 
    save(tmp_name,'job','spm12def','cat12def');
    clear spm12def cat12;
    
    % matlab command, cprintferror=1 for simple printing         
    matlab_cmd = sprintf(...
        ['"global cprintferror; cprintferror=1; addpath %s %s %s %s;load %s; ' ...
         'global defaults; defaults=spm12def; clear defaults; '...
         'global cat; cat=cat12def; clear cat; cat_run(job); "'],...
      spm('dir'),fullfile(spm('dir'),'toolbox','cat12'),...
        fullfile(spm('dir'),'toolbox','OldNorm'),fullfile(spm('dir'),'toolbox','DARTEL'), tmp_name);

    % log-file for output
    log_name = ['catlog_main_' logdate '_log' sprintf('%02d',i) '.txt'];

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
    cat_check_system_output(status,result);
    
    
  
    test = 0; lim = 20; ptime = 0.5;
    while test<lim
      if ~exist(log_name,'file')
        pause(ptime); 
        test = test + ptime; 
        if test>=lim
          cat_io_cprintf('warn','"%s" not exist after %d seconds! Proceed! \n',log_name,lim)
        end
      else 
        test = inf; 
        edit(log_name);
      end
    end

    edit(log_name);
    fprintf('\nCheck %s for logging information.\n',spm_file(log_name,'link','edit(''%s'')'));
    fprintf('_______________________________________________________________\n');

    % starting many large jobs can cause servere MATLAB errors
    pause(1 + rand(1) + job.nproc + numel(job.data)/100);
  end

  job = update_job(job);
  varargout{1} = vout_job(job);
  return
end

if isfield(job,'printPID') && job.printPID 
  cat_display_matlab_PID
end

job = update_job(job);

varargout{1} = run_job(job);

return
%_______________________________________________________________________
function job = update_job(job)

  % set GUI specific parameter if available
  FN = {}; GUIfields = {'registration','segmentation','admin','surface'}; 
  for fnj=1:numel(GUIfields)
    if isfield(job.extopts,GUIfields{fnj})
       FN = [FN;{GUIfields{fnj} fieldnames(job.extopts.(GUIfields{fnj}) )'}];
    end
  end
  for fnj=1:size(FN,1)  
    if isfield(job.extopts,FN{fnj,1})
      for fni=1:numel(FN{fnj,2})
        if isfield(job.extopts.(FN{fnj,1}),FN{fnj,2}{fni})
          job.extopts.(FN{fnj,2}{fni}) = job.extopts.(FN{fnj,1}).(FN{fnj,2}{fni});
        %$else
        %  fprintf('err1: %s\n', FN{fnj,2}{fni});
        end
      end
      job.extopts = rmfield(job.extopts,FN{fnj,1}); % this is just a GUI field! 
    end 
  end
  
  % get defaults
  def = cat_get_defaults;
  
  if isfield(job.extopts,'restypes')
    def.extopts.restype = (char(fieldnames(job.extopts.restypes))); 
    def.extopts.resval  = job.extopts.restypes.(def.extopts.restype);
  end
  
  def.extopts.lazy   = 0;
  def.opts.fwhm      = 1;
  def.nproc          = 0; 
  
   
  % ROI atlas maps
  if isfield(job.output,'ROImenu') % expert/developer GUI that allows control each atlas map 
    if isfield(job.output.ROImenu,'atlases')
      % image output
      def.output.atlases = job.output.ROImenu.atlases;
      def.output.ROI     = any(cell2mat(struct2cell(job.output.ROImenu.atlases))); 
    else
      def.output.atlases = struct();
      def.output.ROI     = 0; 
    end
    job = cat_io_checkinopt(job,def);
  end
  
  if ~isfield(job.output,'atlases') 
    % default GUI that only allow to switch on the settings defined in the default file 
    if ~isfield(job.extopts,'atlas')
      job.extopts.atlas  = def.extopts.atlas;
    end
    
    job.output.atlases   = struct();
    if job.output.ROI 
      % if output, than use the parameter of the default file
      job.output.atlases = cell2struct(job.extopts.atlas(:,4)',spm_str_manip(job.extopts.atlas(:,1),'tr')',2);
      job.output.ROI     = any(cell2mat(struct2cell(job.output.atlases))); 
    end
  end
 
  
  % ROI export 
  for ai = 1:size(job.extopts.atlas,1)
    [pp,ff,ee]  = spm_fileparts(job.extopts.atlas{ai,1}); 
    job.extopts.atlas{ai,4} = job.extopts.atlas{ai,2}<=cat_get_defaults('extopts.expertgui') && ...
      exist(job.extopts.atlas{ai,1},'file') && isfield(def.output,'atlases') && isfield(def.output.atlases,ff) && def.output.atlases.(ff);
  end
  job = cat_io_checkinopt(job,def);
  if ~isfield(job.extopts,'restypes')
    job.extopts.restypes.(def.extopts.restype) = job.extopts.resval;  
  end

  % handling of SPM biasoptions for specific GUI entry
  if isfield(job.opts,'bias')
    if isfield(job.opts.bias,'biasfwhm')
      job.opts.biasstr  = 0; 
      job.opts.biasfwhm = job.opts.bias.biasfwhm; 
      job.opts.biasreg  = job.opts.bias.biasreg; 
    elseif isfield(job.opts.bias,'biasstr')
      job.opts.biasstr  = job.opts.bias.biasstr; 
    end
    job.opts = rmfield(job.opts,'bias'); 
  end
  % the extopts.biasstr controls and overwrites (biasstr>0) the SPM biasreg and biasfwhm parameter
  %   biasstr  = [0.01  0.25  0.50  0.75  1.00] ... result in ?
  %   biasreg  = [0.01  0.0032  0.0010  0.0003  0.0001] ? and ?
  %   biasfwhm = [30 45 60 75 90] for "30 + 60*biasstr? 
  %   biasfwhm = [30.32  42.65  60  84.39 118.71)] for "10^(5/6 + biasstr/3)?  .. allows lower fields 
  if job.opts.biasstr>0 % update biasreg and biasfwhm only if biasreg>0
    % limits only describe the SPM standard range
    job.opts.biasreg	= min(  10 , max(  0 , 10^-(job.opts.biasstr*2 + 2) ));
    job.opts.biasfwhm	= min( inf , max( 30 , 30 + 60*(1-job.opts.biasstr) ));  
  end
 
  
  %% find and check the Dartel templates
  [tpp,tff,tee] = spm_fileparts(job.extopts.darteltpm{1});
  job.extopts.darteltpm{1} = fullfile(tpp,[tff,tee]); 
  numpos = min([strfind(tff,'Template_1')]) + 8;
  if isempty(numpos)
    error('CAT:cat_main:TemplateNameError', ...
    ['Could not find the string "Template_1" in Dartel template that \n'...
     'indicates the first file of the Dartel template. \n' ...
     'The given filename is "%s.%s" \n'],tff,tee);
  end
  job.extopts.darteltpms = cat_vol_findfiles(tpp,[tff(1:numpos) '*' tff(numpos+2:end) tee],struct('depth',1));
  
  % if we also have found Template_0 we have to remove it from the list
  if numel(job.extopts.darteltpms)==7 
    if ~isempty(strfind(job.extopts.darteltpms{1},'Template_0'))
      for i=1:6, job.extopts.darteltpms{i} = job.extopts.darteltpms{i+1}; end
      job.extopts.darteltpms(7) = [];
    end
  end
  
  job.extopts.darteltpms(cellfun('length',job.extopts.darteltpms)~=length(job.extopts.darteltpm{1}))=[]; % remove to short/long files
  if numel(job.extopts.darteltpms)~=6 && any(job.extopts.regstr==4)
    %%
    files = ''; for di=1:numel(job.extopts.darteltpms), files=sprintf('%s\n  %s',files,job.extopts.darteltpms{di}); end
    error('CAT:cat_main:TemplateFileError', ...
     ['Could not find the expected 6 Dartel template files (Template_1 to Template_6). \n' ...
      'Found %d templates: %s'],numel(job.extopts.darteltpms),files);
  end

  % find and check the Shooting templates
  [tpp,tff,tee] = spm_fileparts(job.extopts.shootingtpm{1});
  job.extopts.shootingtpm{1} = fullfile(tpp,[tff,tee]); 
  numpos = min([strfind(tff,'Template_0')]) + 8;
  if isempty(numpos)
    error('CAT:cat_main:TemplateNameError', ...
    ['Could not find the string "Template_0" in Shooting template that \n'...
     'indicates the first file of the Shooting template. \n' ...
     'The given filename is "%s.%s" \n'],tff,tee);
  end
  job.extopts.shootingtpms = cat_vol_findfiles(tpp,[tff(1:numpos) '*' tff(numpos+2:end) tee],struct('depth',1));
  job.extopts.shootingtpms(cellfun('length',job.extopts.shootingtpms)~=length(job.extopts.shootingtpm{1}))=[]; % remove to short/long files
  if numel(job.extopts.shootingtpms)~=5 && any(job.extopts.regstr~=4)
    %%
    files = ''; for di=1:numel(job.extopts.shootingtpms), files=sprintf('%s\n  %s',files,job.extopts.shootingtpms{di}); end
    error('CAT:cat_main:TemplateFileError', ...
     ['Could not find the expected 5 Shooting template files (Template_0 to Template_4).\n' ...
      'Found %d templates: %s'],numel(job.extopts.shootingtpms),files);
  end
  
  
  % check range of str variables
  FN = {'WMHCstr','LASstr','BVCstr','gcutstr','cleanupstr','mrf'};
  for fni=1:numel(FN)
    if ~isfield(job.extopts,FN{fni})  
      job.extopts.(FN{fni}) = max(0,min(1,job.extopts.(FN{fni})));
    end
  end


  
  % deselect ROI output and print warning if ROI output is true and dartel template was changed
  [pth,nam] = spm_fileparts(job.extopts.darteltpm{1});
  if isempty(strfind(nam,'MNI152')) && strcmp(job.extopts.species,'human') && cat_get_defaults('output.ROI')  %~strcmp(nam,'Template_1_IXI555_MNI152')
    warning('DARTEL:template:change',...
      ['Dartel template was changed: Please be aware that ROI analysis \n' ...
       'and other template-specific options cannot be used and ROI \n ' ...
       'output has been deselected.']);
    job.output.ROI = 0;
  end
  
  
  % set boundary box by Template properties if box inf
  if job.extopts.regstr(1)==4
    Vd       = spm_vol([job.extopts.darteltpm{1} ',1']);
  else
    Vd       = spm_vol([job.extopts.shootingtpm{1} ',1']);
  end
  [bb,vox] = spm_get_bbox(Vd, 'old');  
  if bb(1)>bb(2), bbt=bb(1); bb(1)=bb(2); bb(2)=bbt; clear bbt; end
  % Removed BB defintion in GUI and default file in november 2011, because
  % it did not work (need changes in Dartel/Shooting processing) and is not required yet.
  %if job.extopts.bb(1)>job.extopts.bb(2), bbt=job.extopts.bb(1); job.extopts.bb(1)=job.extopts.bb(2); job.extopts.bb(2)=bbt; clear bbt; end
  %job.extopts.bb(isinf(job.extopts.bb))=nan; 
  %job.extopts.bb  = [ min(bb(1,1:3) , job.extopts.bb(1,1:3) ) ; max(bb(2,1:3) , job.extopts.bb(2,1:3) ) ];
  job.extopts.bb = bb; 
  
  job.extopts.vox( isinf(job.extopts.vox) | isnan(job.extopts.vox) ) = []; 
  if isempty( job.extopts.vox ),  job.extopts.vox = cat_get_defaults('extopts.vox'); end 
  job.extopts.vox = abs( job.extopts.vox );
  
  % prepare tissue priors and number of gaussians for all 6 classes
  [pth,nam,ext] = spm_fileparts(job.opts.tpm{1});
  clsn = min(6,numel(spm_vol(fullfile(pth,[nam ext])))); 
  tissue = struct();
  for i=1:clsn;
    tissue(i).ngaus = job.opts.ngaus(i);
    tissue(i).tpm = [fullfile(pth,[nam ext]) ',' num2str(i)];
  end
  
  tissue(1).warped = [job.output.GM.warped  (job.output.GM.mod==1)        (job.output.GM.mod==2)       ];
  tissue(1).native = [job.output.GM.native  (job.output.GM.dartel==1)     (job.output.GM.dartel==2)    ];
  tissue(2).warped = [job.output.WM.warped  (job.output.WM.mod==1)        (job.output.WM.mod==2)       ];
  tissue(2).native = [job.output.WM.native  (job.output.WM.dartel==1)     (job.output.WM.dartel==2)    ];
  tissue(3).warped = [job.output.CSF.warped (job.output.CSF.mod==1)       (job.output.CSF.mod==2)      ];
  tissue(3).native = [job.output.CSF.native (job.output.CSF.dartel==1)    (job.output.CSF.dartel==2)   ];

  % never write class 4-6
  for i=4:6;
    tissue(i).warped = [0 0 0];
    tissue(i).native = [0 0 0];
  end

  job.channel  = struct('vols',{job.data});
  job.tissue   = tissue;
return;

%_______________________________________________________________________
function vout = run_job(job)
  vout   = vout_job(job);

  % load tpm priors 
  tpm = char(cat(1,job.tissue(:).tpm));
  tpm = spm_load_priors8(tpm);

  for subj=1:numel(job.channel(1).vols),
    % __________________________________________________________________
    % Separation for old and new try-catch blocks of matlab. The new
    % try-catch block has to be in a separate file to avoid an error.
    % Both functions finally call cat_run_job.
    % See also cat_run_newcatch and cat_run_newcatch.
    % __________________________________________________________________
    %if job.extopts.ignoreErrors
      if cat_io_matlabversion>20072 
        cat_run_newcatch(job,tpm,subj); 
      else
        % inactive because of unclear error messages
        %cat_run_oldcatch(job,tpm,subj);
        if job.extopts.APP == 1070
          cat_run_job1070(job,tpm,subj); 
        else
          cat_run_job(job,tpm,subj); 
        end
      end
    %else
    %  cat_run_job(job,tpm,subj);
    %end
  end

  colormap(gray)
  
  if isfield(job,'nproc') && job.nproc>0 
    fprintf('\n%s',repmat('_',1,72));
    fprintf('\nCAT12 Segmentation job finished.\n');
  end
return
%_______________________________________________________________________

function vout = vout_job(job)
% ----------------------------------------------------------------------
% create output structure for SPM batch mode
% ----------------------------------------------------------------------

n     = numel(job.channel(1).vols);

parts = cell(n,4); % fileparts

biascorr    = {};
wbiascorr   = {};
ibiascorr   = {};
wibiascorr  = {};
ribiascorr  = {};
aibiascorr  = {};
label       = {};
wlabel      = {};
rlabel      = {};
alabel      = {};
catreport   = {};
lhcentral   = {};
rhcentral   = {};
lhthickness = {};
rhthickness = {};
roi         = {};
fordef      = {};
invdef      = {};
jacobian    = {};

if job.extopts.subfolders
  roifolder    = 'label';
  surffolder   = 'surf';
  mrifolder    = 'mri';
  reportfolder = 'report';
else
  roifolder    = '';
  surffolder   = '';
  mrifolder    = '';
  reportfolder = '';
end

for j=1:n,
    [parts{j,:}] = spm_fileparts(job.channel(1).vols{j});
end

% CAT report XML file
% ----------------------------------------------------------------------
roi = cell(n,1);
for j=1:n
    catreport{j} = fullfile(parts{j,1},reportfolder,['cat_',parts{j,2},'.xml']);
end


% lh/rh central surface and thickness
% ----------------------------------------------------------------------
if job.output.surface,
    lhcentral = cell(n,1);
    rhcentral = cell(n,1);
    lhthickness = cell(n,1);
    rhthickness = cell(n,1);
    if job.output.surface == 2
        lccentral = cell(n,1);
        rccentral = cell(n,1);
        lcthickness = cell(n,1);
        rcthickness = cell(n,1);
    end
    for j=1:n
        lhcentral{j} = fullfile(parts{j,1},surffolder,['lh.central.',parts{j,2},'.gii']);
        rhcentral{j} = fullfile(parts{j,1},surffolder,['rh.central.',parts{j,2},'.gii']);
        lhthickness{j} = fullfile(parts{j,1},surffolder,['lh.thickness.',parts{j,2}]);
        rhthickness{j} = fullfile(parts{j,1},surffolder,['rh.thickness.',parts{j,2}]);
        if job.output.surface == 2
            lhcentral{j} = fullfile(parts{j,1},surffolder,['lc.central.',parts{j,2},'.gii']);
            rhcentral{j} = fullfile(parts{j,1},surffolder,['rc.central.',parts{j,2},'.gii']);
            lhthickness{j} = fullfile(parts{j,1},surffolder,['lc.thickness.',parts{j,2}]);
            rhthickness{j} = fullfile(parts{j,1},surffolder,['rc.thickness.',parts{j,2}]);
        end
    end
end


% XML label
% ----------------------------------------------------------------------
if job.output.ROI,
    roi = cell(n,1);
    for j=1:n
        roi{j} = fullfile(parts{j,1},roifolder,['catROI_',parts{j,2},'.xml']);
    end
end

% bias
% ----------------------------------------------------------------------
if job.output.bias.native,
    biascorr = cell(n,1);
    for j=1:n
        biascorr{j} = fullfile(parts{j,1},mrifolder,['m',parts{j,2},'.nii']);
    end
end

if job.output.bias.warped,
    wbiascorr = cell(n,1);
    for j=1:n
        wbiascorr{j} = fullfile(parts{j,1},mrifolder,['wm',parts{j,2},'.nii']);
    end
end

if job.output.bias.dartel==1,
    rbiascorr = cell(n,1);
    for j=1:n
        rbiascorr{j} = fullfile(parts{j,1},mrifolder,['rm',parts{j,2},'.nii']);
    end
end

if job.output.bias.dartel==2,
    abiascorr = cell(n,1);
    for j=1:n
        abiascorr{j} = fullfile(parts{j,1},mrifolder,['rm',parts{j,2},'_affine.nii']);
    end
end

% intensity corrected bias
% ----------------------------------------------------------------------
if job.output.las.native,
    ibiascorr = cell(n,1);
    for j=1:n
        ibiascorr{j} = fullfile(parts{j,1},mrifolder,['mi',parts{j,2},'.nii']);
    end
end

if job.output.las.warped,
    wibiascorr = cell(n,1);
    for j=1:n
        wibiascorr{j} = fullfile(parts{j,1},mrifolder,['wmi',parts{j,2},'.nii']);
    end
end

if job.output.las.dartel==1,
    ribiascorr = cell(n,1);
    for j=1:n
        ribiascorr{j} = fullfile(parts{j,1},mrifolder,['rmi',parts{j,2},'.nii']);
    end
end

if job.output.las.dartel==2,
    aibiascorr = cell(n,1);
    for j=1:n
        aibiascorr{j} = fullfile(parts{j,1},mrifolder,['rmi',parts{j,2},'_affine.nii']);
    end
end


% label
% ----------------------------------------------------------------------
if job.output.label.native,
    label = cell(n,1);
    for j=1:n
        label{j} = fullfile(parts{j,1},mrifolder,['p0',parts{j,2},'.nii']);
    end
end

if job.output.label.warped,
    wlabel = cell(n,1);
    for j=1:n
        wlabel{j} = fullfile(parts{j,1},mrifolder,['wp0',parts{j,2},'.nii']);
    end
end

if job.output.label.dartel==1,
    rlabel = cell(n,1);
    for j=1:n
        rlabel{j} = fullfile(parts{j,1},mrifolder,['rp0',parts{j,2},'.nii']);
    end
end

if job.output.label.dartel==2,
    alabel = cell(n,1);
    for j=1:n
        alabel{j} = fullfile(parts{j,1},mrifolder,['rp0',parts{j,2},'_affine.nii']);
    end
end


% tissues
% ----------------------------------------------------------------------
tiss = struct('p',{},'rp',{},'rpa',{},'wp',{},'mwp',{},'m0wp',{});
for i=1:numel(job.tissue),
    if job.tissue(i).native(1),
        tiss(i).p = cell(n,1);
        for j=1:n
            tiss(i).p{j} = fullfile(parts{j,1},mrifolder,['p',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(2),
        tiss(i).rp = cell(n,1);
        for j=1:n
            tiss(i).rp{j} = fullfile(parts{j,1},mrifolder,['rp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).native(3),
        tiss(i).rpa = cell(n,1);
        for j=1:n
            tiss(i).rpa{j} = fullfile(parts{j,1},mrifolder,['rp',num2str(i),parts{j,2},'_affine.nii']);
        end
    end
    if job.tissue(i).warped(1),
        tiss(i).wp = cell(n,1);
        for j=1:n
            tiss(i).wp{j} = fullfile(parts{j,1},mrifolder,['wp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(2),
        tiss(i).mwp = cell(n,1);
        for j=1:n
            tiss(i).mwp{j} = fullfile(parts{j,1},mrifolder,['mwp',num2str(i),parts{j,2},'.nii']);
        end
    end
    if job.tissue(i).warped(3),
        tiss(i).m0wp = cell(n,1);
        for j=1:n
            tiss(i).m0wp{j} = fullfile(parts{j,1},mrifolder,['m0wp',num2str(i),parts{j,2},'.nii']);
        end
    end
end


% deformation fields
% ----------------------------------------------------------------------
if job.output.warps(1),
    fordef = cell(n,1);
    for j=1:n
        fordef{j} = fullfile(parts{j,1},mrifolder,['y_',parts{j,2},'.nii']);
    end
else
    fordef = {};
end

if job.output.warps(2),
    invdef = cell(n,1);
    for j=1:n
        invdef{j} = fullfile(parts{j,1},mrifolder,['iy_',parts{j,2},'.nii']);
    end
else
    invdef = {};
end


% jacobian
% ----------------------------------------------------------------------
if job.output.jacobian.warped,
    jacobian = cell(n,1);
    for j=1:n
        jacobian{j} = fullfile(parts{j,1},mrifolder,['wj_',parts{j,2},'.nii']);
    end
end


% ----------------------------------------------------------------------
vout  = struct('tiss',tiss,'label',{label},'wlabel',{wlabel},'rlabel',{rlabel},'alabel',{alabel},...
               'biascorr',{biascorr},'wbiascorr',{wbiascorr},'roi',{roi},'ibiascorr',{ibiascorr},...
               'wibiascorr',{wibiascorr},'ribiascorr',{ribiascorr},'aibiascorr',{aibiascorr},...
               'invdef',{invdef},'fordef',{fordef},'jacobian',{jacobian},'catreport',{catreport},...
               'lhcentral',{lhcentral},'rhcentral',{rhcentral},'lhthickness',{lhthickness},...
               'rhthickness',{rhthickness});
%_______________________________________________________________________
return

%=======================================================================
function [data,err] = remove_already_processed(job,verb)
  if ~exist('verb','var'), verb=0; end
  remove = []; err = zeros(size(job));
  cat_io_cprintf('warn','Lazy processing: \n');
  for subj = 1:numel(job.data)
    [lazy,err(subj)] = checklazy(job,subj,verb); 
    if lazy
      remove = [remove subj];
    end
  end
  cat_io_cprintf('warn','  Skip %d subjects!\n',numel(remove));
  data = job.data(setxor(1:numel(job.data),remove)); 
  cat_io_cprintf([0 0.4 0.6],'\n\nProcess:\n');
  for subj = 1:numel(data)
    cat_io_cprintf([0 0.4 0.6],sprintf(' Code%3d: "%s"\n',err(subj),data{subj}));
  end
  cat_io_cprintf('warn',sprintf('  Process %d subjects!\n',numel(data)));
return
%=======================================================================
function [lazy,FNok] = checklazy(job,subj,verb)
  if job.extopts.subfolders
    roifolder    = 'label';
    surffolder   = 'surf';
    mrifolder    = 'mri';
    reportfolder = 'report';
  else
    roifolder    = '';
    surffolder   = '';
    mrifolder    = '';
    reportfolder = '';
  end

  lazy = 0;
  
  [pp,ff] = spm_fileparts(job.data{subj}); 
  catxml  = fullfile(pp,reportfolder,['cat_' ff '.xml']);
  
  FNok = 0;
  if exist(catxml,'file')

    xml         = cat_io_xml(catxml);
    
    FNopts      = fieldnames(job.opts); 
    FNextopts   = fieldnames(job.extopts);
    FNok        = 1; 
    FNextopts   = setxor(FNextopts,{'LAB','lazy','mrf','atlas','NCstr','resval'});
   
   
    %% check opts
    if isempty(FNopts) || isempty(FNextopts) || ...
       ~isfield(xml.parameter,'opts') || ~isfield(xml.parameter,'extopts')
      return
    end
    for fni=1:numel(FNopts)
      if ~isfield(xml.parameter.opts,FNopts{fni})
        FNok = 2; break
      end
      if ischar(xml.parameter.opts.(FNopts{fni}))
        if ischar(job.opts.(FNopts{fni}))
          if ~strcmp(xml.parameter.opts.(FNopts{fni}),job.opts.(FNopts{fni}))
            FNok = 3; break
          end
        else
          if ~strcmp(xml.parameter.opts.(FNopts{fni}),job.opts.(FNopts{fni}){1})
            FNok = 4; break
          end
        end
      else
        if isnumeric(job.opts.(FNopts{fni}))
          if strcmp(FNopts{fni},'ngaus') && numel(xml.parameter.opts.(FNopts{fni}))==4
            % nothing to do (skull-stripped case)
          else
            if xml.parameter.opts.(FNopts{fni}) ~= job.opts.(FNopts{fni})
              FNok = 5; break
            end
          end
        elseif ischar(job.opts.(FNopts{fni}))
          if ~strcmp(xml.parameter.opts.(FNopts{fni}),job.opts.(FNopts{fni})); 
            FNok = 5; break
          end
        end
      end
    end
    if FNok~=1 % different opts
      return
    end

    %% check extopts
    for fni=1:numel(FNextopts)
      if ~isfield(xml.parameter.extopts,FNextopts{fni})
        FNok = 6; break
      end
      if ischar(xml.parameter.extopts.(FNextopts{fni}))
        if ischar(job.extopts.(FNextopts{fni}))
          if ~strcmp(xml.parameter.extopts.(FNextopts{fni}),job.extopts.(FNextopts{fni}))
            FNok = 7; break
          end
        else
          if ~strcmp(xml.parameter.extopts.(FNextopts{fni}),job.extopts.(FNextopts{fni}){1})
            FNok = 8; break
          end
        end
      elseif iscell(xml.parameter.extopts.(FNextopts{fni}))
        if numel(xml.parameter.extopts.(FNextopts{fni}))~=numel(job.extopts.(FNextopts{fni}))
          FNok = 9; break
        end
        for fnic = 1:numel(xml.parameter.extopts.(FNextopts{fni}))
          if iscell(xml.parameter.extopts.(FNextopts{fni}){fnic})
            for fnicc = 1:numel(xml.parameter.extopts.(FNextopts{fni}){fnic})
              if xml.parameter.extopts.(FNextopts{fni}){fnic}{fnicc} ~= job.extopts.(FNextopts{fni}){fnic}{fnicc}
                FNok = 10; break
              end
            end
            if FNok==10; break; end
          else
            try
              if any(xml.parameter.extopts.(FNextopts{fni}){fnic} ~= job.extopts.(FNextopts{fni}){fnic})
                FNok = 11; break
              end
            catch
                FNok = 11;
            end
            if FNok==11; break; end
          end
          if FNok==11 || FNok==10; break; end
        end
      elseif isstruct(xml.parameter.extopts.(FNextopts{fni}))
        FNX = fieldnames(xml.parameter.extopts.(FNextopts{fni}));
        for fnic = 1:numel(FNX)
          if any(xml.parameter.extopts.(FNextopts{fni}).(FNX{fnic}) ~= job.extopts.(FNextopts{fni}).(FNX{fnic}))
            FNok = 12; break
          end
          if FNok==12; break; end
        end
      else
        % this did not work anymore due to the GUI subfields :/
        if 0 %any(xml.parameter.extopts.(FNextopts{fni}) ~= job.extopts.(FNextopts{fni}))
          FNok = 13; break
        end
      end
    end
    if FNok~=1 % different extopts
      return
    end
    

    % check output
    
    % surface
    if job.output.surface && exist(fullfile(pp,surffolder),'dir')
      Pcentral = cat_vol_findfiles(fullfile(pp,surffolder),['*h.central.' ff '.gii']);
      if  numel(Pcentral)==1 % miss ROI xml
        return
      end
    end
    
    % rois
    if job.output.ROI && ~exist(fullfile(pp,roifolder,['catROI_' ff '.xml']),'file')  % miss ROI xml
      return
    end
      
    %% volumes 
    FNO = fieldnames(job.vout);
    for fnoi = 1:numel(FNO)
      if isempty(job.vout.(FNO{fnoi}))
        continue
      elseif iscell(job.vout.(FNO{fnoi}))
         if ~isempty(job.vout.(FNO{fnoi}){subj}) && ~exist(job.vout.(FNO{fnoi}){subj},'file')
           FNok = 14; break
         end
      elseif isstruct(job.vout.(FNO{fnoi}))
        for si = numel(job.vout.(FNO{fnoi}))
          FNOS = fieldnames(job.vout.(FNO{fnoi})); 
          for fnosi = 1:numel(FNOS)
            if isempty([job.vout.(FNO{fnoi})(si).(FNOS{fnosi})])
              continue
            elseif ~exist(job.vout.(FNO{fnoi})(si).(FNOS{fnosi}){subj},'file')
              FNok = 14; break
            end
          end
        end
      end
      if FNok==14 % miss volumes
        return
      end
    end
    %%
    
    lazy = FNok==1; 
      
  end
 
  if lazy 
    cat_io_cprintf('warn','  "%s" \n',job.data{subj});
  end
return
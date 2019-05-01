function varargout = cat_update(update)
% check for new CAT updates
%
% FORMAT [sts, msg] = cat_update(update)
% sts    - status code:
%        NaN - SPM server not accessible
%        Inf - no updates available
%        0   - SPM installation up to date
%        n   - new revision <n> is available for download
% msg    - string describing outcome, that would otherwise be displayed.
% update - allow installation of update
% 
% This function will connect itself to the SBM server, compare the
% version number of the updates with the one of the CAT12 installation 
% currently in the MATLAB path and will display the result.
%_______________________________________________________________________
% Christian Gaser
% $Id: cat_update.m 1219 2017-11-19 23:28:59Z gaser $

rev = '$Rev: 1219 $';

url = 'http://www.neuro.uni-jena.de/cat12/';

if ~nargin
    update = false;
else
    update = true;
end

r = 0;

% get current release number
[n, r] = cat_version;
r = str2double(r);

% get new release numbers
[s,sts] = urlread(url);
if ~sts
  sts = NaN;
  msg = sprintf('Cannot access %s. Please check your proxy and/or firewall to allow access.\nYou can download your update at %s\n',url,url); 
  if ~nargout, error(msg); else varargout = {sts, msg}; end
  return
end

n = regexp(s,'cat12_r(\d.*?)\.zip','tokens');
if isempty(n)
  sts= Inf;
  msg = 'There are no updates available yet.';
  if ~nargout, fprintf([blanks(9) msg '\n']);
  else varargout = {sts, msg}; end
  return;
else
  % get largest release number
  rnew = [];
  for i=1:length(n)
    rnew = [rnew str2double(n{i})];
  end
  rnew = max(rnew);
end

if rnew > r
    sts = n;
    msg = sprintf('         A new version of CAT12 is available on:\n');
    msg = [msg sprintf('   %s\n',url)];
    msg = [msg sprintf('        (Your version: %d - New version: %d)\n',r,rnew)];
    if ~nargout, fprintf(msg); else varargout = {sts, msg}; end
else
    sts = 0;
    msg = sprintf('Your version of CAT12 is up to date.');
    if ~nargout, fprintf([blanks(9) msg '\n']);
    else varargout = {sts, msg}; end
    return
end

if update
    overwrite = spm_input('Update',1,'yes|no',[1 0],1);
    d0 = spm('Dir');
    d = fullfile(spm('Dir'),'toolbox'); 
    if overwrite
      try
        % list mex-files and delete these files to prevent that old
        % compiled files are used
        mexfiles = dir(fullfile(d,'cat12','*.mex*'));
        for i=1:length(mexfiles)
          name = fullfile(d,'cat12',mexfiles(i).name);
          spm_unlink(name);
        end
        
        % delete old atlas files
        atlasfiles = dir(fullfile(d,'cat12','atlases_surfaces','*.*'));
        for i=1:length(atlasfiles)
          name = fullfile(d,'cat12','atlases_surfaces',atlasfiles(i).name);
          spm_unlink(name);
        end

        % delete old atlas files with 32k meshes
        atlasfiles = dir(fullfile(d,'cat12','atlases_surfaces_32k','*.*'));
        for i=1:length(atlasfiles)
          name = fullfile(d,'cat12','atlases_surfaces_32k',atlasfiles(i).name);
          spm_unlink(name);
        end

        % delete old surface template files
        templatefiles = dir(fullfile(d,'cat12','templates_surfaces','*.*'));
        for i=1:length(templatefiles)
          name = fullfile(d,'cat12','templates_surfaces',templatefiles(i).name);
          spm_unlink(name);
        end

        % delete old surface template files with 32k meshes
        templatefiles = dir(fullfile(d,'cat12','templates_surfaces_32k','*.*'));
        for i=1:length(templatefiles)
          name = fullfile(d,'cat12','templates_surfaces_32k',templatefiles(i).name);
          spm_unlink(name);
        end

        % delete old volume template files 
        templatefiles = dir(fullfile(d,'cat12','templates_1.50mm','*.*'));
        for i=1:length(templatefiles)
          name = fullfile(d,'cat12','templates_1.50mm',templatefiles(i).name);
          spm_unlink(name);
        end

        lastwarn('');
        delete(get(0,'Children')); spm('clean'); evalc('spm_rmpath'); drawnow
        m = '          Download and install CAT12...\n';
        if ~nargout, fprintf(m); else varargout = {sts, [msg m]}; end
        s = unzip([url sprintf('cat12_r%d.zip',rnew)], d);
        m = sprintf('         Success: %d files have been updated.\n',numel(s));
        if ~nargout, fprintf(m); else varargout = {sts, [msg m]}; end
        addpath(d0);
        rehash
        rehash toolboxcache;
        toolbox_path_cache
        eval(['spm fmri;clear cat_version;spm_cat12']);
      catch
        le = lasterror;
        switch le.identifier
            case 'MATLAB:checkfilename:urlwriteError'
                fprintf('          Update failed: cannot download update file.\n');
            otherwise
                fprintf('\n%s\n',le.message);
        end
      end
      [warnmsg, msgid] = lastwarn;
      switch msgid
        case ''
        case 'MATLAB:extractArchive:unableToCreate'
            fprintf('          Update failed: check folder permission.\n');
        case 'MATLAB:extractArchive:unableToOverwrite'
            fprintf('          Update failed: check file permissions.\n');
        otherwise
            fprintf('          Update failed: %s.\n',warnmsg);
      end
    end
end

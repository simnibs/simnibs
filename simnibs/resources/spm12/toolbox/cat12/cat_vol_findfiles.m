function [filesfound,numberfound] = cat_vol_findfiles(varargin)
% cat_vol_findfiles  - Linux/UNIX-like find command/function
% ______________________________________________________________________
%
% FORMAT:       [files, number] = cat_vol_findfiles(startfolder, patterns [, opts])
%
% Input fields:
%
%       startfolder 1xN char, folder where to start search
%                   may alternatively contain the pattern at the end
%       patterns    1xP cell of string or 1xN char file pattern(s)
%       opts        struct, optional parameters in the form of
%        .cellstr   1x1 double, if set and ~= 0 return as cellstr (default)
%        .chararr   1x1 double, if set and ~= 0 return as char array
%        .depth     1x1 double, sets both minimum and maximum depth
%        .dirs      1x1 double, if set and ~= 0 find dirs instead of files
%        .filesize  1x1 double, if set and > 0, only matching files
%        .maxdepth  1x1 double, maximum depth to find files in
%        .mindepth  1x1 double, minimum depth to find files in
%        .maxage    1x1 double, seconds file must have changed in
%        .minage    1x1 double, seconds file must not have changed in
%        .oneperdir 1x1 double, if set and ~= 0 only first match (per dir)
%        .relative  1xN char, prepend this instead of startfolder
%
% Output fields:
%
%       files       Fx1 cell array (or FxL char array, if requested)
%       number      1x1 number of files found
%
% when used as a command, multiple opts can be given as multiple arguments,
% separated by spaces (' '):
%
% cat_vol_findfiles /search/folder *mypattern*.txt depth=3 oneperdir=1 relative=./
%
% when used in functional context, a second return value, the number
% of matching files, can be obtained:
%
% [files, number] = cat_vol_findfiles('/where','*.txt');
%
% NOTE: the minage/maxage feature only fully works when the system
% returns English-style month in calls to dir. i.e. under Linux, set the
% LANG environmental setting to 'en_US' before starting up MatLab
%
% Version:  v0.9d
% Build:    14071015
% Date:     Jul-10 2014, 3:28 PM EST
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://neuroelf.net/
%
% ______________________________________________________________________
%
% Copyright (c) 2010 - 2014, Jochen Weber
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Columbia University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% ______________________________________________________________________
% $Id: cat_vol_findfiles.m 1117 2017-03-16 15:31:20Z dahnke $ 

%#ok<*EFIND>
%#ok<*AGROW>

% enough arguments ?
if nargin < 1 || ...
   ((~ischar(varargin{1}) || ...
     isempty(varargin{1})) && ...
    (~iscell(varargin{1}) || ...
     isempty(varargin{1}) || ...
     ~ischar(varargin{1}{1}) || ...
     isempty(varargin{1}{1})))
    help cat_vol_findfiles; return
    %{
    error( ...
        'neuroelf:TooFewArguments',...
        'Too few arguments. Try ''help %s''.',...
        mfilename ...
    );
    %}
end

% for single argument
if nargin == 1

    % check for full path
    [p, f, x] = fileparts(varargin{1}(:)');

    % local file
    if isempty(p) || ...
        strcmp(p, '.')

        % look in current directory
        varargin{1} = pwd;

    % full path
    else

        % use the path given
        varargin{1} = p;
    end

    % extend varargin
    varargin{2} = [f, x];

% allow local path with options
elseif ischar(varargin{1}) && ...
    any(varargin{1} == '*') && ...
    ischar(varargin{2}) && ...
    any(varargin{2} == '=')

    % pass on
    filesfound = cat_vol_findfiles(pwd, varargin{:});

    % and return
    if iscell(filesfound)
        filesfound = filesfound(:);
        numberfound = numel(filesfound);
    else
        numberfound = size(filesfound, 1);
    end
    return;
end

% shortcut variable
fsep = filesep;

% sanity checks: startfolder
startfolder = varargin{1};
if ischar(startfolder) && ...
   ~isempty(startfolder)

    % does startfolder contain an asterisk?
    if any(startfolder == '?' | startfolder == '*')

        % then put startfolder into a cell array to treat this case
        [filesfound, numberfound] = cat_vol_findfiles({startfolder}, varargin{2:end});
        if iscell(filesfound)
            filesfound = filesfound(:);
        end
        return;
    end

% startfolder is a list of folders
elseif iscell(startfolder) && ...
   ~isempty(startfolder)

    % generate expanded list
    nstartfolder = cell(0);

    % iterate over items
    for nelem = 1:numel(startfolder)

        % must be non-empty char
        if ~ischar(startfolder{nelem}) || ...
            isempty(startfolder{nelem})
            error( ...
                'neuroelf:BadArgument',...
                'Bad startfolder argument.' ...
            );
        end

        % expand pattern?
        if ~any(startfolder{nelem} == '?' | startfolder{nelem} == '*')

            % put into final list
            nstartfolder{end+1} = startfolder{nelem};

        % or look up matching folders
        else

            % split along possible path separators
            [pparts, cparts] = splittocell(startfolder{nelem}, '/\', 1, 1);
            if any('/\:' == startfolder{nelem}(1))
                pparts = [{''}, pparts];
                cparts = cparts + 1;
            end

            % look for first pattern
            for cpart = 1:cparts
                if any(pparts{cpart} == '?' | pparts{cpart} == '*')
                    break;
                end
            end

            % if the pattern occurs in first part, look in current dir
            if cpart==1
                [pfolders, npfolders] = cat_vol_findfiles('.', ...
                    pparts{1}, ...
                    struct('dirs', 1, 'depth', 1));

            % otherwise glue first parts together and look up matches
            else
                spart = gluetostring({pparts{1:(cpart-1)}}, fsep);
                if exist(spart, 'dir') > 0
                    [pfolders, npfolders] = cat_vol_findfiles(spart, ...
                        pparts{cpart}, ...
                        struct('dirs', 1, 'depth', 1));
                else
                    pfolders = cell(0, 1);
                    npfolders = 0;
                end
            end

            % the pattern was not in last part
            if cpart < cparts

                % put remaining parts back on
                for ppart = 1:npfolders
                    pfolders{ppart} = [pfolders{ppart} fsep ...
                        gluetostring({pparts{(cpart+1):end}},fsep)];
                end
            end

            % put results at end of list
            nstartfolder = [nstartfolder, pfolders];
        end
    end

    % we start with no files found
    filesfound = cell(0);

    % for each folder in list
    for nelem = 1:numel(nstartfolder)

        % if there (iteratively) remains a pattern char, redo this
        if any(nstartfolder{nelem} == '?' | nstartfolder{nelem} == '*')
            filesfound = [filesfound(1:end); ...
                cat_vol_findfiles({nstartfolder{nelem}}, varargin{2:end})];

        % otherwise get files in this folder and put at end of array
        elseif exist(nstartfolder{nelem}, 'dir') == 7
            filesfound = [filesfound(1:end); ...
                cat_vol_findfiles(nstartfolder{nelem}, varargin{2:end})];
        end
    end

    % report the total number found
    if iscell(filesfound)
        filesfound = filesfound(:);
        numberfound = numel(filesfound);
    else
        numberfound = size(filesfound, 1);
    end
    return;

% illegal first argument
else
    error( ...
        'neuroelf:BadArgument',...
        'Bad startfolder argument.'...
    );
end

% we're now going for the single folder case and startfolder is OK

% append missing fsep if needed
if startfolder(end) ~= fsep
    startfolder = [startfolder fsep];
end

% startfolder exists?
if exist(startfolder,'dir') ~= 7
    error( ...
        'neuroelf:FolderNotFound',...
        'Startfolder (%s) not found or no folder.',...
        strrep(startfolder,'\','\\') ...
    );
end

% - sanity checks patterns
patterns = varargin{2};

% default is cell, otherwise
if ~iscell(patterns)

    % if is character, put into cell
    if ischar(patterns)
        patterns = { patterns };

    % otherwise bail out
    else
        error( ...
            'neuroelf:BadArgument', ...
            'Patterns must either be a single string or a cell array!' ...
        );
    end
end

% check each pattern
for count = 1:numel(patterns)

    % only accept chars
    if ~ischar(patterns{count}) || ...
        isempty(patterns{count})
        patterns{count} = '*';
    end
end
patterns = unique(patterns(:));

% option argument parsing, default options
if nargin < 3
    opt.dirs = 0;
    opt.filesize = 0;
    opt.maxdepth = 0;
    opt.mindepth = 0;
    opt.maxage = -1;
    opt.minage = -1;
    opt.oneperdir = 0;
    opt.relative = 0;
    opt.return = 'cellarr';
    opt.rfolder = startfolder;

% parse options
else
    opt = varargin{3};

    % non-struct options
    if ~isstruct(opt)

        % yet start with default struct
        clear opt;
        opt.dirs = 0;
        opt.filesize = 0;
        opt.maxdepth = 0;
        opt.mindepth = 0;
        opt.maxage = -1;
        opt.minage = -1;
        opt.oneperdir = 0;
        opt.relative = 0;
        opt.return = 'cellarr';
        opt.rfolder = startfolder;

        % parse all arguments
        for acount=3:nargin

            % char option
            if ischar(varargin{acount})

                % special case: -a#A#d#Dors#
                % (minage, maxage, depth, dirs, oneperdir, relative=, size)
                if ~isempty(varargin{acount}) && ...
                    varargin{acount}(1) == '-'
                    optarg = varargin{acount}(2:end);
                    while ~isempty(optarg)
                        switch (optarg(1))
                            case {'a'}
                                if numel(optarg) > 1 && ...
                                    optarg(2) >= '0' && ...
                                    optarg(2) <= '9'
                                    optpos = findfirst(optarg(2:end) < '0' | optarg(2:end) > '9');
                                    if isempty(optpos)
                                        oval = optarg(2:end);
                                        optarg = '';
                                    else
                                        oval = optarg(2:optpos);
                                        optarg(1:optpos) = [];
                                    end
                                    opt.minage = str2double(oval);
                                else
                                    warning( ...
                                        'neuroelf:BadOption', ...
                                        'Option -a (minage) requires numeric input.' ...
                                    );
                                    optarg(1) = [];
                                end
                            case {'A'}
                                if numel(optarg) > 1 && ...
                                    optarg(2) >= '0' && ...
                                    optarg(2) <= '9'
                                    optpos = findfirst(optarg(2:end) < '0' | optarg(2:end) > '9');
                                    if isempty(optpos)
                                        oval = optarg(2:end);
                                        optarg = '';
                                    else
                                        oval = optarg(2:optpos);
                                        optarg(1:optpos) = [];
                                    end
                                    opt.maxage = str2double(oval);
                                else
                                    warning( ...
                                        'neuroelf:BadOption', ...
                                        'Option -A (maxage) requires numeric input.' ...
                                    );
                                    optarg(1) = [];
                                end
                            case {'d'}
                                if numel(optarg) > 1 && ...
                                    optarg(2) >= '0' && ...
                                    optarg(2) <= '9'
                                    optpos = findfirst(optarg(2:end) < '0' | optarg(2:end) > '9');
                                    if isempty(optpos)
                                        oval = optarg(2:end);
                                        optarg = '';
                                    else
                                        oval = optarg(2:optpos);
                                        optarg(1:optpos) = [];
                                    end
                                    opt.maxdepth = str2double(oval);
                                    opt.mindepth = opt.maxdepth;
                                else
                                    warning( ...
                                        'neuroelf:BadOption', ...
                                        'Option -d (depth) requires numeric input.' ...
                                    );
                                    optarg(1) = [];
                                end
                            case {'D'}
                                opt.dirs = 1;
                                optarg(1) = [];
                            case {'o'}
                                opt.oneperdir = 1;
                                optarg(1) = [];
                            case {'r'}
                                opt.relative = 1;
                                opt.rfolder = '';
                                optarg(1) = [];
                            case {'s'}
                                if numel(optarg) > 1 && ...
                                    optarg(2) >= '0' && ...
                                    optarg(2) <= '9'
                                    optpos = findfirst(optarg(2:end) < '0' | optarg(2:end) > '9');
                                    if isempty(optpos)
                                        oval = optarg(2:end);
                                        optarg = '';
                                    else
                                        oval = optarg(2:optpos);
                                        optarg(1:optpos) = [];
                                    end
                                    opt.filesize = str2double(oval);
                                else
                                    warning( ...
                                        'neuroelf:BadOption', ...
                                        'Option -s (filesize) requires numeric input.' ...
                                    );
                                    optarg(1) = [];
                                end
                            otherwise
                                warning( ...
                                    'neuroelf:BadOption', ...
                                    'Unknown cat_vol_findfiles option: %s.', ...
                                    optarg(1) ...
                                );
                                optarg(1) =[];
                        end
                    end
                    continue;
                end

                % get argument name
                argnv = splittocell(varargin{acount}, '=', 1);
                oname = argnv{1};

                % and possible option value
                if numel(argnv) > 1
                    oval = argnv{2};
                else
                    oval = '';
                end

                % only accept known arguments
                switch lower(oname)

                    % option: cellstr, set return type
                    case {'cellstr'}
                        oval = str2double(oval);
                        if oval ~= 0
                            opt.return = 'cellstr';
                        end

                    % option: chararr, set return type
                    case {'chararr'}
                        oval = str2double(oval);
                        if oval ~= 0
                            opt.return = 'chararr';
                        end

                    % option: depth (min and max)
                    case {'depth'}
                        if str2double(oval) >= 0
                            opt.maxdepth = str2double(oval);
                            opt.mindepth = str2double(oval);
                        else
                            opt.maxdepth = 0;
                            opt.mindepth = 0;
                        end

                    % option: dirs, set lookup type
                    case {'dirs'}
                        oval = str2double(oval);
                        if oval == 0
                            opt.dirs = 0;
                        else
                            opt.dirs = 1;
                        end

                    % option: filesize
                    case {'filesize'}
                        if str2double(oval) >= 0
                            opt.filesize = fix(str2double(oval));
                        else
                            opt.filesize = 0;
                        end

                    % option: maxdepth
                    case {'maxdepth'}
                        if str2double(oval) >= 0
                            opt.maxdepth = fix(str2double(oval));
                        else
                            opt.maxdepth = 0;
                        end

                    % option: mindepth
                    case {'mindepth'}
                        if str2double(oval) >= 0
                            opt.mindepth = fix(str2double(oval));
                        else
                            opt.mindepth = 0;
                        end

                    % option: maxage
                    case {'maxage'}
                        if str2double(oval) >= 0
                            opt.maxage = fix(str2double(oval));
                        else
                            opt.maxage = -1;
                        end

                    % option: minage
                    case {'minage'}
                        if str2double(oval) >= 0
                            opt.minage = fix(str2double(oval));
                        else
                            opt.minage = -1;
                        end

                    % option: oneperdir
                    case {'oneperdir'}
                        oval = str2double(oval);
                        if oval == 0
                            opt.oneperdir = 0;
                        else
                            opt.oneperdir = 1;
                        end

                    % option: relative
                    case 'relative'
                        noval = str2double(oval);
                        if ~isnan(noval)
                            opt.relative = noval;
                            if noval < 1
                                opt.rfolder = startfolder;
                            else
                                opt.rfolder = ['.' fsep];
                            end
                        else
                            opt.relative = 1;
                            opt.rfolder = oval;
                        end;
                end
            end
        end

    % struct option argument
    else

        % make sure options are present
        if ~isfield(opt, 'dirs')
            opt.dirs = 0;
        end
        if ~isfield(opt, 'filesize')
            opt.filesize = 0;
        end
        if ~isfield(opt, 'maxdepth')
            if isfield(opt, 'depth')
                opt.maxdepth = opt.depth;
            else
                opt.maxdepth = 0;
            end
        end
        if ~isfield(opt, 'mindepth')
            if isfield(opt, 'depth')
                opt.mindepth = opt.depth;
            else
                opt.mindepth = 0;
            end
        end
        if ~isfield(opt, 'maxage')
            opt.maxage = -1;
        end
        if ~isfield(opt, 'minage')
            opt.minage = -1;
        end
        if ~isfield(opt, 'oneperdir')
            opt.oneperdir = 0;
        end
        if  isfield(opt, 'rfolder')
            opt = rmfield(opt, 'rfolder');
        end
        if ~isfield(opt, 'relative')
            opt.relative = 0;
            opt.rfolder = startfolder;
        else
            if ischar(opt.relative)
                opt.rfolder=opt.relative;
                opt.relative=1;
            else
                if double(opt.relative) >= 1
                    opt.rfolder = ['.' fsep];
                    opt.relative = 1;
                else
                    opt.rfolder = startfolder;
                    opt.relative = 0;
                end
            end
        end
        if ~isfield(opt, 'return')
            opt.return = 'cellarr';
        end
    end
end

% more interdependent checks now
if isfield(opt,'cellstr') && ...
    opt.cellstr > 0
    opt.return = 'cellstr';
end
if isfield(opt,'chararr') && ...
    opt.chararr > 0
    opt.return = 'chararr';
end
if opt.dirs ~= 0
    opt.dirs = 1;
end

% check option types
if ~isa(opt.filesize, 'double')
    opt.filesize = 0;
end
if ~isa(opt.maxdepth, 'double')
    opt.maxdepth = 0;
end
if ~isa(opt.mindepth, 'double')
    opt.mindepth = 0;
end
if ~isa(opt.maxage, 'double')
    opt.maxage = -1;
end
if ~isa(opt.minage, 'double')
    opt.minage = -1;
end
if opt.oneperdir ~= 0
    opt.oneperdir = 1;
end
if opt.relative ~= 0
    opt.relative = 1;
else opt.rfolder = startfolder;
end

% calculate age here
opt.maxage=opt.maxage / 86400;
if opt.maxage < 0
    opt.maxage = -1;
end
opt.minage = opt.minage / 86400;
if opt.minage < 0
    opt.minage = -1;
end

% make call for files
if opt.dirs == 0
    filesfound = findsubfiles( ...
        startfolder, patterns, 1, ...
        opt.filesize, opt.mindepth, opt.maxdepth, ...
        opt.minage, opt.maxage, ...
        opt.oneperdir, opt.rfolder);

% make call for dirs
else
    filesfound = findsubdirs( ...
        startfolder, patterns, 1, ...
        opt.mindepth, opt.maxdepth, ...
        opt.minage, opt.maxage, ...
        opt.oneperdir, opt.rfolder);
end

% return the correct number of values
if nargout > 1
    numberfound = size(filesfound, 1);
end

% return correct type
if strcmpi(opt.return(:)', 'chararr')
	filesfound = char(filesfound);
end
end
% - end of cat_vol_findfiles(...)


% %%%%internal functions%%%%


% findsubfiles
function found = findsubfiles(path, patterns, adepth, fsize, sdepth, mdepth, mnage, mxage, operdir, relative)

% start with zero files found
nfound = 0;
found = cell(0, 1);
mfilesep = filesep;

% first, recursively handle all subfolders, if depth is still valid
if mdepth == 0 || ...
    adepth < mdepth

    % get list of files and folders, and size of list
    ilist = dir(path);
    slist = numel(ilist);

    % get isdir flag into array and find indices of dirs
    [ilistd(1:slist)] = [ilist(:).isdir];
    ilistd = find(ilistd > 0);

    % check items
    for count = ilistd

        % don't heed . and ..
        if strcmp(ilist(count).name, '.') || ...
            strcmp(ilist(count).name, '..')
            continue;
        end;

        % find files in subdirs
        filestoadd = findsubfiles([path ilist(count).name mfilesep], ...
            patterns, adepth + 1, fsize, sdepth, mdepth, ...
            mnage, mxage, operdir, [relative ilist(count).name mfilesep]);
        sfound = numel(filestoadd);

        % if files found
        if sfound > 0
            nfoundfrm = nfound + 1;
            nfoundnew = nfound + sfound;
            found(nfoundfrm:nfoundnew, 1) = filestoadd(:);
            nfound = nfoundnew;
        end
    end
end

% then, if depth is valid, add files to the output
if sdepth == 0 || ...
    sdepth <= adepth

    % only get time if needed
    if any([mnage, mxage] >= 0)
        rnow = now;
    end;

    % number of patterns
    spatt = numel(patterns);
    for pcount = 1:spatt

        % no "*" pattern
        if ~any(patterns{pcount} == '*') && ...
           ~any(patterns{pcount} == '?')
            ilist = dir([path patterns{pcount} '*']);
            if isempty(ilist)
                continue;
            end
            ilistn = {ilist(:).name};
            if any(strcmp(ilistn, patterns{pcount}))
                nfound = nfound + 1;
                found{nfound, 1} = [relative patterns{pcount}];
            end
            continue;

        % find matching entries with ?
        elseif any(patterns{pcount} == '?')
            ilist = dir([path strrep(strrep(patterns{pcount}, '?', '*'), '**', '*')]);
            ilistn = {ilist(:).name};
            ilist(cellfun('isempty', regexp(ilistn, ...
                [strrep(strrep(strrep(patterns{pcount}, '.', '\.'), ...
                '?', '.'), '*', '.*') '$']))) = [];

        % and without ?
        else
            ilist = dir([path patterns{pcount}]);
        end
        slist = numel(ilist);

        % get isdir flag into array and remove dirs from list
        ilistd = [];
        [ilistd(1:slist)] = [ilist(:).isdir];
        ilist(ilistd > 0) = [];
        slist = numel(ilist);

        % if only one per dir
        if operdir == 1
            count = 1;

            % reject all non-matching
            while count <= slist && ...
                 ((mnage >= 0 && (rnow - datenum(ilist(count).date)) < mnage) || ...
                  (mxage >= 0 && (rnow - datenum(ilist(count).date)) > mxage) || ...
                  (fsize ~= 0 && ilist(count).bytes ~= fsize))
                count = count + 1;
            end

            % choose first if still valid
            if count <= slist
                nfound = nfound + 1;
                found{nfound, 1} = [relative ilist(count).name];
            end

        % otherwise check all
        else

            % iterate over all
            for count = 1:slist

                % reject non-matching
                if ((mnage >= 0 && (rnow - datenum(ilist(count).date)) < mnage) || ...
                    (mxage >= 0 && (rnow - datenum(ilist(count).date)) > mxage) || ...
                    (fsize ~= 0 && ilist(count).bytes ~= fsize))
                    continue;
                end

                % accept rest
                nfound = nfound + 1;
                found{nfound, 1} = [relative ilist(count).name];
            end
        end
    end

    % linearize found
    found = found(:);
end
% end of function findsubfiles
end
% findsubdirs
function found = findsubdirs(path, patterns, adepth, sdepth, mdepth, mnage, mxage, operdir, relative)

  % start with zero dirs found
  nfound = 0;
  found = cell(0, 1);
  mfilesep = filesep;

  % first, recursively handle all subfolders, if depth is still valid
  if mdepth == 0 || ...
      adepth < mdepth

      % get list of files and folders, and size of list
      ilist = dir(path);
      slist = numel(ilist);

      % get isdir flag into array
      [ilistd(1:slist)] = [ilist(:).isdir];

      % find indices of dirs
      ilistd = find(ilistd > 0);

      % check items
      for count = ilistd

          % don't heed . and ..
          if strcmp(ilist(count).name, '.') || ...
              strcmp(ilist(count).name, '..')
              continue;
          end

          % iterate over subdirs
          filestoadd = findsubdirs([path ilist(count).name mfilesep], ...
              patterns, adepth + 1, sdepth, mdepth, ...
              mnage, mxage, operdir, [relative ilist(count).name mfilesep]);
          sfound = numel(filestoadd);

          % if dirs founds
          if sfound > 0
              nfoundfrm = nfound + 1;
              nfoundnew = nfound + sfound;
              found(nfoundfrm:nfoundnew, 1) = filestoadd(:);
              nfound = nfoundnew;
          end
      end
  end

  % then, if depth is valid, add folders to the output
  if sdepth == 0 || ...
      sdepth <= adepth

      % only get time if needed
      if any([mnage, mxage]>=0)
          rnow = now;
      end;

      % number of patterns
      spatt = numel(patterns);
      for pcount = 1:spatt

          % no "*" or "?" pattern
          if ~any(patterns{pcount} == '*') && ...
             ~any(patterns{pcount} == '?')
              if exist([path patterns{pcount}], 'dir') == 7
                  nfound = nfound + 1;
                  found{nfound} = [relative patterns{pcount}];
              end
              continue;

          % "?" pattern/s
          elseif any(patterns{pcount} == '?')
              ilist = dir([path strrep(strrep(patterns{pcount}, '?', '*'), '**', '*')]);
              ilistn = {ilist(:).name};
              ilist(cellfun('isempty', regexp(ilistn, ...
                  [strrep(strrep(strrep(patterns{pcount}, '.', '\.'), ...
                  '?', '.'), '*', '.*') '$']))) = [];

          % "*" pattern/s
          else
              ilist = dir([path patterns{pcount}]);
          end

          % get matching entries
          slist = numel(ilist);

          % get isdir flag into array and remove files from list
          ilistd = [];
          [ilistd(1:slist)] = [ilist(:).isdir];
          ilist(~ilistd) = [];
          slist = numel(ilist);

          % if only one per dir
          if operdir == 1
              count = 1;

              % reject all non-matching entries
              while count <= slist && ...
                   ((mnage >= 0 && (rnow - datenum(ilist(count).date)) < mnage) || ...
                    (mxage >= 0 && (rnow - datenum(ilist(count).date)) > mxage))
                  count = count + 1;
              end

              % find next entry
              while count <= slist

                  % still reject . and ..
                  if  strcmp(ilist(count).name, '.') || ...
                      strcmp(ilist(count).name, '..')
                      count = count + 1;
                      continue;
                  end

                  % get next entry
                  nfound = nfound + 1;
                  found{nfound, 1} = [relative ilist(count).name];
                  break;
              end

          % otherwise check all
          else

              % iterate over all
              for count = 1:slist

                  % reject non-matching
                  if ((mnage >= 0 && (rnow - datenum(ilist(count).date)) < mnage) || ...
                      (mxage >= 0 && (rnow - datenum(ilist(count).date)) > mxage))
                      continue;
                  end

                  % reject . and ..
                  if strcmp(ilist(count).name, '.') || ...
                      strcmp(ilist(count).name, '..')
                      continue;
                  end

                  % accept others
                  nfound = nfound + 1;
                  found{nfound, 1} = [relative ilist(count).name];
              end
          end
      end

      % linearize found
      found = found(:);
  end
  % end of function findsubdirs
end
function [linetocell,cellcount] = splittocell(varargin)
% splittocell  - split a delimited string into a cell array
%
% usage is straight forward:
%
% FORMAT:         [outcell,count] = splittocell(string[,delimiters,multi])
%
% Input fields:
%    string       string to split
%    delimiters   char array containing one or more delimiters
%                 if left empty -> char(9) == <TAB>
%    multi        must be '1' (numeric) to be effective, if set
%                 multiple delimiters will be treated as one
%
% Output fields:
%    outcell      cell array containing the tokens after split
%    count        number of tokens in result

  % no arguments -> help me!
  if nargin == 0, help(mfilename); return; end

  % initialize return values and varargin{3}
  linetocell=cell(0);
  cellcount =0;
  multidelim=0;



  % do we have useful input ?
  if ~ischar(varargin{1}) | length(varargin{1})==0, return; end
  line=varargin{1};
  if size(line,2) ~= prod(size(line))
      dispdebug('splittocell: input must be a 1xN shaped char array!',4);
      return;
  end

  % are any other arguments specified
  if nargin < 2 | ~ischar(varargin{2})
      delimiter = char(9);
  else
      delimiter = reshape(varargin{2},1,prod(size(varargin{2})));
      if nargin > 2 & isnumeric(varargin{3}) & varargin{3} ~= 0, multidelim = 1; end
  end

  % multi-delimitting requested ?
  if multidelim == 0

      % set initial parameters
      ldelim=size(delimiter,2);
      lline =size(line,2);

      % find occurences of delimiter
      if ldelim==1
          cpos=[(1-ldelim),find(line==delimiter)];
      else
          cpos=[(1-ldelim),findstr(line,delimiter)];
      end
      lcpos =size(cpos,2);

      % any delimiter found at all ?
      if lcpos==1, cellcount=1; linetocell={line}; return; end

      % large array?
      if lcpos < 4096

          % line doesn't end with delimiter ?
          if cpos(lcpos) <= (lline-ldelim)
              % then make it look like it was...
              cpos =[cpos lline+1];
              lcpos=lcpos+1;
          end

          % extract substrings
          for dpos=1:(lcpos-1)
              linetocell{end+1} = line(cpos(dpos)+ldelim:cpos(dpos+1)-1);
          end

      else

          % get good rate
          crate = min(384,floor(lcpos^0.666));

          % iterate over parts
          linetocell={};
          for cmpos = 1:crate:(lcpos-crate)
              linetocell = [linetocell,splittocell(line(cpos(cmpos)+ldelim:cpos(cmpos+crate)-1),delimiter,multidelim)];
          end
          linetocell = [linetocell,splittocell(line(cpos(cmpos+crate)+ldelim:cpos(end)-1),delimiter,multidelim)];

      end

  else

      % set initial parameters
      ldelim=size(delimiter,2);
      lline =size(line,2);

      % find occurences of delimiter
      pdelim = [0];
      for cdelim=1:ldelim
          pdelim = union(pdelim,find(line==delimiter(cdelim)));
      end
      if pdelim(end) ~= lline, pdelim(end+1)=lline+1; end
      lpdel = size(pdelim,2);

      % extract substrings
      if pdelim(2)==1, linetocell{end+1} = ''; end
      for ppdel=1:(lpdel-1)
          if (pdelim(ppdel+1)-1) ~= pdelim(ppdel)
              linetocell{end+1} = line(pdelim(ppdel)+1:pdelim(ppdel+1)-1);
          end
      end

  end

  cellcount=length(linetocell);
end




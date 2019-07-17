function [newjobs, uind] = cfg_load_jobs(job)

% function newjobs = cfg_load_jobs(job)
%
% Load a list of possible job files, return a cell list of jobs.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_load_jobs.m 7024 2017-02-21 17:53:55Z guillaume $

rev = '$Rev: 7024 $'; %#ok

if ischar(job)
    filenames = cellstr(job);
else
    filenames = job;
end
[ufilenames, unused, uind] = unique(filenames);
ujobs = cell(size(ufilenames));
usts  = false(size(ufilenames));
for cf = 1:numel(ufilenames)
    [ujobs{cf}, usts(cf)] = load_single_job(ufilenames{cf});
end
sts   = usts(uind);
uind  = uind(sts);
newjobs = ujobs(uind);

function [matlabbatch, sts] = load_single_job(filename)
[p,nam,ext] = fileparts(filename);
switch ext
    case '.xml'
        try
            loadxml(filename,'matlabbatch');
        catch
            cfg_message('matlabbatch:initialise:xml','LoadXML failed: ''%s''',filename);
        end
    case '.mat'
        try
            S=load(filename);
            matlabbatch = S.matlabbatch;
        catch
            cfg_message('matlabbatch:initialise:mat','Load failed: ''%s''',filename);
        end
    case '.m'
        try
            str = fileread(filename);
            if isdeployed && strncmp(str,'V1MCC',5) % mcc compiled script
                eval(nam);
            else
                eval(str);
            end
        catch
            cfg_message('matlabbatch:initialise:m','Eval failed: ''%s''',filename);
        end
        if ~exist('matlabbatch','var')
            cfg_message('matlabbatch:initialise:m','No matlabbatch job found in ''%s''', filename);
        end
    otherwise
        cfg_message('matlabbatch:initialise:unknown','Unknown extension: ''%s''', filename);
end
if exist('matlabbatch','var')
    sts = true;
else
    sts = false;
    matlabbatch = [];
end

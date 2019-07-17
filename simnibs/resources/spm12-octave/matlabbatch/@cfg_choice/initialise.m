function item = initialise(item, val, dflag)

% function item = initialise(item, val, dflag)
% Initialise a configuration tree with values. If val is a job
% struct/cell, only the parts of the configuration that are present in
% this job will be initialised. If dflag is true, then matching items
% from item.values will be initialised. If dflag is false, the matching
% item from item.values will be added to item.val and initialised after
% copying.
% If val has the special value '<DEFAULTS>', the entire configuration
% will be updated with values from .def fields. If a .def field is
% present in a cfg_leaf item, the current default value will be inserted,
% possibly replacing a previously entered (default) value. If dflag is
% true, defaults will only be set in item.values. If dflag is false,
% defaults will be set for both item.val and item.values.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: initialise.m 4898 2012-09-05 13:40:16Z volkmar $

rev = '$Rev: 4898 $'; %#ok

if strcmp(val,'<UNDEFINED>')
    val = struct([]);
end
if strcmp(val,'<DEFAULTS>')
    item = initialise_def(item, val, dflag);
elseif isstruct(val)
    item = initialise_job(item, val, dflag);
elseif iscell(val) && numel(val) == 1 && isstruct(val{1})
    item = initialise_job(item, val{1}, dflag);
else
    cfg_message('matlabbatch:initialise', ...
                'Can not initialise %s: job is not a struct.', ...
                gettag(item));
end;

function item = initialise_def(item, val, dflag)
if ~dflag
    % initialise defaults both in current job and in defaults
    citem = subsref(item, substruct('.','val'));
    for k = 1:numel(citem)
        citem{k} = initialise(citem{k}, val, dflag);
    end;
    item = subsasgn(item, substruct('.','val'), citem);
end;
for k = 1:numel(item.values)
    item.values{k} = initialise(item.values{k}, val, dflag);
end;

function item = initialise_job(item, val, dflag)
% Modify job before initialisation
if ~dflag && ~isempty(item.cfg_item.rewrite_job)
    val = feval(item.cfg_item.rewrite_job, val);
end
% Determine possible tags
vtags = fieldnames(val);
utags = false(size(vtags));

if dflag % set defaults
    for k = 1:numel(item.values)
        % find field in values that corresponds to one of the branch vals
        vi = strcmp(gettag(item.values{k}), vtags);
        if any(vi) % field names are unique, so there will be at most one match
            item.values{k} = initialise(item.values{k}, ...
                                        val.(vtags{vi}), dflag);
        end;
    end;
else
    % select matching values struct, initialise and assign it to val
    % field
    item.cfg_item.val = {};
    if ~isempty(vtags)
        for k = 1:numel(item.values)
            if strcmp(gettag(item.values{k}), vtags{1})
                item.cfg_item.val{1} = initialise(item.values{k}, ...
                                                  val.(vtags{1}), dflag);
                utags(1) = true;
                break;
            end;
        end;
    end;
end;

% Check whether some fields were not found in child tags
if ~dflag && any(~utags)
    cfg_message('matlabbatch:initialise', ...
                'Item %s: No field(s) named%s', ...
                gettag(item), sprintf('\n%s', vtags{~utags}));
end

function chanind = selectchannels(this, channels)
% Method for getting channel indices based on labels and/or types
% FORMAT  res = selectchannels(this, label)
% this       - MEEG object
% channels   - string or cell array of labels that may also include 
%              'all', or types ('EEG', 'MEG' etc.)
%
% res        - vector of channel indices matching labels
%__________________________________________________________________________
% Copyright (C) 2010-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: selectchannels.m 7253 2018-02-04 17:20:57Z vladimir $

if ischar(channels)
    channels = {channels};
end

chanind = [];

for i = 1:numel(channels)    
    if strncmpi('regexp_', channels{i}, 7)
        re        = channels{i}(8:end);
        match     = regexp(chanlabels(this), re);
        chanind   = [chanind find(~cellfun('isempty', match))];
    else
        cind    = indchannel(this, channels{i});
        if ~isempty(cind)
            chanind = [chanind cind];
        elseif ismember(upper(channels{i}), ...
                {'ALL','MEG', 'MEGPLANAR', 'MEGMAG', 'MEGGRAD', 'MEGCOMB','EEG',...
                'EOG', 'ECG', 'EMG', 'LFP', 'SRC', 'PHYS', 'ILAM', 'OTHER', 'REF', 'REFMAG', 'REFGRAD'})
            chanind = [chanind indchantype(this, upper(channels{i}))];
        end
    end
    
    if any(size(chanind) == 0)
        chanind = [];
    end
end

chanind = unique(chanind);
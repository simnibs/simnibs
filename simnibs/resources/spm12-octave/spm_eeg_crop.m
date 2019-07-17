function D = spm_eeg_crop(S)
% Reduce the data size by cutting in time and frequency
% FORMAT D = spm_eeg_crop(S)
%
% S        - optional input struct
%  fields of S:
%   D        - MEEG object or filename of M/EEG mat-file with epoched data
%   timewin  - time window to retain {in PST ms}
%   freqwin  - frequency window to retain
%   channels - cell array of channel labels or 'all' [default]
%   prefix   - prefix for the output file [default: 'p']
%
% Output:
% D        - MEEG object (also written on disk)
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_crop.m 7132 2017-07-10 16:22:58Z guillaume $

SVNrev = '$Rev: 7132 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Crop M/EEG data'); spm('Pointer','Watch');

if ~isfield(S, 'prefix'),       S.prefix   = 'p';           end
if ~isfield(S, 'timewin'),      S.timewin  = [-Inf Inf];    end
if ~isfield(S, 'freqwin'),      S.freqwin  = [-Inf Inf];    end
if ~isfield(S, 'channels'),     S.channels = 'all';         end

D = spm_eeg_load(S.D);

isTF = strncmpi(D.transformtype,'TF',2);

timeind = D.indsample(1e-3*(min(S.timewin))):D.indsample(1e-3*(max(S.timewin)));
if isempty(timeind) || any(isnan(timeind))
    error('Selected time window is invalid.');
end

if isTF
    freqind = D.indfrequency(min(S.freqwin)):D.indfrequency(max(S.freqwin));
    if isempty(freqind) || any(isnan(freqind))
        error('Selected frequency window is invalid.');
    end
end

chanind = D.selectchannels(S.channels);

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
if isTF
    Dnew = clone(D, [S.prefix fname(D)], [length(chanind) length(freqind) length(timeind) D.ntrials]);
    Dnew = frequencies(Dnew, ':', D.frequencies(freqind));
else
    Dnew = clone(D, [S.prefix fname(D)], [length(chanind) length(timeind) D.ntrials]);
end

Dnew = timeonset(Dnew, D.time(timeind(1)));

Dnew = chanlabels(Dnew, ':', D.chanlabels(chanind));
Dnew = badchannels(Dnew, ':', badchannels(D, chanind));
Dnew = chantype(Dnew, ':', chantype(D, chanind));
Dnew = units(Dnew, ':', units(D, chanind));
Dnew = coor2D(Dnew, ':', coor2D(D, chanind));

if isequal(Dnew.type, 'continuous')
    ev = Dnew.events;
    if ~isempty(ev)
        ev = ev([ev.time]>=Dnew.time(1) & [ev.time]<=Dnew.time(end));
        Dnew = events(Dnew, 1, ev);
    end
end

%-Copy data
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Channels',sprintf('%d/%d',numel(chanind),D.nchannels));%-#
if isTF, fprintf('%-40s: %30s\n','Frequencies',sprintf('%d/%d',numel(freqind),D.nfrequencies)); end %-#
fprintf('%-40s: %30s\n','Samples',sprintf('%d/%d',numel(timeind),D.nsamples)); %-#
spm_progress_bar('Init', D.ntrials, 'Trials copied');
if D.ntrials > 100, Ibar = floor(linspace(1, D.ntrials, 100));
else Ibar = 1:D.ntrials; end

onsets = zeros(1,D.ntrials);

for i = 1:D.ntrials
    
    if isTF
        Dnew(:, :, :, i) = D(chanind, freqind, timeind, i);
    else
        Dnew(:, :, i) = D(chanind, timeind, i);
    end
    
    if D.trialonset(i) ~= 0
        onsets(i) = D.trialonset(i)+ D.time(timeind(1))-D.time(1);
    end
    
    if any(Ibar == i), spm_progress_bar('Set', i); end
end

Dnew = trialonset(Dnew, ':', onsets);

spm_progress_bar('Clear');

%-Save the new M/EEG dataset
%--------------------------------------------------------------------------
Dnew = Dnew.history(mfilename, S);
save(Dnew);

D = Dnew;

%-Cleanup
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','Crop M/EEG data: done'); spm('Pointer','Arrow');

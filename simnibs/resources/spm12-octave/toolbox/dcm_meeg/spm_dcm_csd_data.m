function DCM = spm_dcm_csd_data(DCM)
% gets cross-spectral density data-features using a VAR model
% FORMAT DCM = spm_dcm_csd_data(DCM)
% DCM    -  DCM structure
% requires
%
%    DCM.xY.Dfile        - name of data file
%    DCM.M.U             - channel subspace
%    DCM.options.trials  - trial to evaluate
%    DCM.options.Tdcm    - time limits
%    DCM.options.Fdcm    - frequency limits
%    DCM.options.D       - Down-sampling
%
% sets
%
%    DCM.xY.pst     - Peristimulus Time [ms] sampled
%    DCM.xY.dt      - sampling in seconds [s] (down-sampled)
%    DCM.xY.U       - channel subspace
%    DCM.xY.y       - cross spectral density over sources
%    DCM.xY.csd     - cross spectral density over sources
%    DCM.xY.It      - Indices of time bins
%    DCM.xY.Ic      - Indices of good channels
%    DCM.xY.Hz      - Frequency bins
%    DCM.xY.code    - trial codes evaluated
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_csd_data.m 6481 2015-06-16 17:01:47Z karl $
 
% Set defaults and Get D filename
%-------------------------------------------------------------------------
try
    Dfile = DCM.xY.Dfile;
catch
    errordlg('Please specify data and trials');
    error('')
end

% ensure spatial modes have been computed (see spm_dcm_csd)
%-------------------------------------------------------------------------
try
    DCM.M.U;
catch
    Nm      = DCM.options.Nmodes;
    DCM.M.U = spm_dcm_eeg_channelmodes(DCM.M.dipfit,Nm);
end

% load D
%--------------------------------------------------------------------------
try
    D = spm_eeg_load(Dfile);
catch
    try
        [p,f]        = fileparts(Dfile);
        D            = spm_eeg_load(f);
        DCM.xY.Dfile = fullfile(pwd,f);
    catch
        try
            [f,p]        = uigetfile('*.mat','please select data file');
            name         = fullfile(p,f);
            D            = spm_eeg_load(name);
            DCM.xY.Dfile = fullfile(name);
        catch
            warndlg([Dfile ' could not be found'])
            return
        end
    end
end

 
% indices of EEG channel (excluding bad channels)
%--------------------------------------------------------------------------
if ~isfield(DCM.xY, 'modality')
    [mod, list] = modality(D, 0, 1);

    if isequal(mod, 'Multimodal')
        qstr = 'Only one modality can be modelled at a time. Please select.';
        if numel(list) < 4
            
            % Nice looking dialog
            %--------------------------------------------------------------
            options = [];
            options.Default = list{1};
            options.Interpreter = 'none';
            DCM.xY.modality = questdlg(qstr, 'Select modality', list{:}, options);
        else
            
            % accomodate more buttons
            %--------------------------------------------------------------
            ind = menu(qstr, list);
            DCM.xY.modality = list{ind};
        end
    else
        DCM.xY.modality = mod;
    end
end



if ~isfield(DCM.xY, 'Ic')
    DCM.xY.Ic  = D.indchantype(DCM.xY.modality,'GOOD');
end

Ic        = DCM.xY.Ic;
Nm        = size(DCM.M.U,2);
DCM.xY.Ic = Ic;

% options
%--------------------------------------------------------------------------
try
    DT    = DCM.options.D;
catch
    DT    = 1;
end
try
    trial = DCM.options.trials;
catch
    trial = 1:D.nconditions;
end
 
% check data are not oversampled (< 4ms)
%--------------------------------------------------------------------------
if DT/D.fsample < 0.004
    DT            = ceil(0.004*D.fsample);
    DCM.options.D = DT;
end
 
 
% get peristimulus times
%--------------------------------------------------------------------------
try
    
    % time window and bins for modelling
    %----------------------------------------------------------------------
    DCM.xY.Time = time(D, [], 'ms'); 
    T1          = DCM.options.Tdcm(1);
    T2          = DCM.options.Tdcm(2);
    [i, T1]     = min(abs(DCM.xY.Time - T1));
    [i, T2]     = min(abs(DCM.xY.Time - T2));
    
    % Time [ms] of down-sampled data
    %----------------------------------------------------------------------
    It          = [T1:DT:T2]';               % indices - bins
    DCM.xY.pst  = DCM.xY.Time(It);           % PST
    DCM.xY.It   = It;                        % Indices of time bins
    DCM.xY.dt   = DT/D.fsample;              % sampling in seconds
    Nb          = length(It);                % number of bins
    
catch
    errordlg('Please specify time window');
    error('')
end
 
% get frequency range
%--------------------------------------------------------------------------
try
    Hz1     = DCM.options.Fdcm(1);          % lower frequency
    Hz2     = DCM.options.Fdcm(2);          % upper frequency
catch
    pst     = DCM.xY.pst(end) - DCM.xY.pst(1);
    Hz1     = max(ceil(2*1000/pst),4);
    if Hz1 < 8;
        Hz2 = 48;
    else
        Hz2 = 128;
    end
end

 
% Frequencies
%--------------------------------------------------------------------------
DCM.xY.Hz  = fix(Hz1:Hz2);             % Frequencies
Nf         = length(DCM.xY.Hz);        % number of frequencies
Ne         = length(trial);            % number of trial types


% Cross spectral density for each trial type
%==========================================================================
condlabels = D.condlist;               % condition or trial type labels
DCM.xY.csd = cell(1,Ne);               % CSD for each condition

w     = min(fix(2/DCM.xY.dt),Nb);      % window length (bins)
m     = 1;                             % retain principal mode
for i = 1:Ne;
   
    % trial indices
    %----------------------------------------------------------------------
    c = D.indtrial(condlabels(trial(i)), 'GOOD');
    fprintf('\nevaluating CSD for condition %i\n',i)
    
    
    % use only the first 512 trial
    %----------------------------------------------------------------------
    try c = c(1:512); end
    Nt    = length(c);
    
    % Get data
    %----------------------------------------------------------------------
    Nw    = max(8*(fix(Nb/w) - 1),1);
    K     = zeros(Nf*Nm*Nm,Nw);
    for k = 1:Nw
        P     = zeros(Nf,Nm,Nm);
        for j = 1:Nt
            Iw  = It((1:w) + fix((k - 1)*w/8));
            Y   = full(double(D(Ic,Iw,c(j))'*DCM.M.U));
            mar = spm_mar(Y,8);
            mar = spm_mar_spectra(mar,DCM.xY.Hz,1/DCM.xY.dt);
            P   = P + mar.P;
        end
        
        % store
        %------------------------------------------------------------------
        K(:,k) = spm_vec(P/Nt);
    end
    
    % retain principal eigenmode
    %----------------------------------------------------------------------
    [u s v]       = spm_svd(K,1);
    P             = u(:,m)*s(m,m)*mean(v(:,m));
    DCM.xY.csd{i} = spm_unvec(P,mar.P);
   
end

 
% place cross-spectral density in xY.y
%==========================================================================
DCM.xY.y    = DCM.xY.csd; 
DCM.xY.U    = DCM.M.U;
DCM.xY.code = condlabels(trial);

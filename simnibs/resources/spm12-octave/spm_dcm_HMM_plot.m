function spm_dcm_HMM_plot(HMM,s)
% FORMAT spm_dcm_HMM_plot(HMM,s)
% subroutine for potting the results of a hidden Markov model of state
% transitions in the DCM
% 
% HMM(s)
%     HMM(s).X  - posterior expectation of hidden states
%     HMM(s).qB - posterior expectation of HMM parameters
%     HMM(s).qb - and Dirichlet concentration parameters
%     HMM(s).qP - posterior expectation of PEB parameters
%     HMM(s).qC - posterior covariances of PEB parameters
%     HMM(s).iP - indices of DCM parameters
%     HMM(s).Ep - posterior expectation of DCM parameters
%     HMM(s).Cp - posterior covariances of DCM parameters
%     HMM(s).L  - free energy components
%     HMM(s).F  - total free energy (model evidence)
%
% s  -  index of HMM structure (number of hidden states)
%       [default: HMM(end)]
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_HMM_plot.m 7272 2018-03-04 13:12:31Z karl $

%  preliminaries: get the hidden Markov model to report
%--------------------------------------------------------------------------
spm_figure('Getwin','HMM'); clf
if nargin < 2, s = length(HMM); end
hmm = HMM(s);

% extract posterior estimates of DCM parameters and likelihoods
%--------------------------------------------------------------------------
N     = numel(hmm.Ep);                          % number of epochs
for t = 1:N
    qp(:,t) = spm_vec(hmm.Ep{t}.A);             % expected connectivity
end
for i = 1:size(hmm.L,1)
    L(i,:)  = hmm.L(i,:) - min(hmm.L(i,:));     % free energy terms
end

% plot sequence of expected states
%--------------------------------------------------------------------------
subplot(4,1,1),imagesc(hmm.X);
title('Hidden states','FontSize',16), ylabel('State')

% epoch-specific fluctuations in DCM parameters
%--------------------------------------------------------------------------
subplot(4,1,2),imagesc(qp);
title('Parameter fluctuations','FontSize',16), ylabel('Connection')

% graphical format, with empirical priors from the HMM
%--------------------------------------------------------------------------
subplot(4,1,3),
plot((1:N),qp(hmm.iP,:)),     hold on
plot((1:N),hmm.qP*hmm.X,':'), hold off
title('Parameter fluctuations','FontSize',16)
xlabel('Time (epochs)'), ylabel('Effective connectivity'), spm_axis tight

% connection strengths associated with each state (PEB Parameters)
%--------------------------------------------------------------------------
subplot(4,2,7),imagesc(hmm.qP); axis square
title('State-dependent connectivity','FontSize',16)
xlabel('State'), ylabel('Connection')

% probability transition matrix (HMM Parameters)
%--------------------------------------------------------------------------
subplot(4,2,8),imagesc(hmm.qB);axis square
title('State transition matrix','FontSize',16)
xlabel('Hidden state'), ylabel('State')

% lower bound on model evidence (DCM, PEB and HMM components)
%--------------------------------------------------------------------------
spm_figure('Getwin','HMM-F'); clf

subplot(2,2,1)
plot([L' sum(L)']);axis square
title('Log-likelihood','FontSize',16)
xlabel('Iteration'), ylabel('(negative) free energy')
legend({'Hidden states (HMM)',...
        'State parameters (HMM)',...
        'Connectivity (DCM)',...
        'Total'},'Location','Best')
    
% and log evidence, if there is more than one model
%----------------------------------------------------------------------    
if length(HMM) > 1
    
    subplot(2,2,2)
    F = spm_cat({HMM.F});
    bar(spm_softmax(F(:)),'c')
    title('Log-evidence','FontSize',16), axis square
    xlabel('Number of states'), ylabel('Log-likehood')
  
end

drawnow

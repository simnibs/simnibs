function [L] = spm_lx_dem(P,M)
% observer matrix for a neural mass model of erps: y = G*x
% FORMAT [G] = spm_lx_dem(P,M)
% x      - state vector
% G      - where y = G*x; G = L*J
%          L = dy/dsource
%          J = dsource/dstate
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_lx_dem.m 4814 2012-07-30 19:56:05Z karl $

% parameterised lead field times source contribution to ECD
%--------------------------------------------------------------------------
L       = spm_erp_L(P,M.dipfit)*P.J;               % lead field per state


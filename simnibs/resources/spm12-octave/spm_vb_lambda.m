function [block] = spm_vb_lambda(Y,block)
% Variational Bayes for GLM-AR models - Update lambda
% FORMAT [block] = spm_vb_lambda(Y,block)
%
% Y      - [T x N] time series 
% block  - data structure (see spm_vb_glmar)
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_lambda.m 6079 2014-06-30 18:25:37Z spm $

if block.verbose
    disp('Updating lambda');
end

p = block.p;
k = block.k;
N = block.N;

for n=1:N
    if p > 0
        % Equation 77 in paper VB1
        Gn          = spm_vb_get_Gn (Y,block,n);
    else
        subblock_n  = [(n-1)*k+1:n*k];
        en          = Y(:,n) - block.X*block.w_mean(subblock_n,1);
        Gn          = trace(block.w_cov{n}*block.XTX) + en'*en;
    end
    % Equation 75 in paper VB1
    block.b_lambda(n,1)    = 1./(Gn./2 + 1./block.b_lambda_prior(n));
    block.mean_lambda(n,1) = block.c_lambda(n)*block.b_lambda(n);
end

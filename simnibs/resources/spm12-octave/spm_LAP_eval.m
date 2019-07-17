function [p,dp] = spm_LAP_eval(M,qu,qh)
% Evaluate precisions for a LAP model
% FORMAT [p dp] = spm_LAP_eval(M,qu,qh)
%
% p.h     - vector of precisions for causal states (v)
% p.g     - vector of precisions for hidden states (x)
%
% dp.h.dx - dp.h/dx
% dp.h.dv - dp.h/dv
% dp.h.dh - dp.h/dh
%
% dp.g.dx - dp.g/dx
% dp.g.dv - dp.g/dv
% dp.g.dg - dp.g/dg
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_LAP_eval.m 6290 2014-12-20 22:11:50Z karl $


% Get states {qu.v{1},qu.x{1}} in hierarchical form (v{i},x{i})
%--------------------------------------------------------------------------
N          = length(M);
v          = cell(N,1);
x          = cell(N,1);
v(1:N - 1) = spm_unvec(qu.v{1},{M(1 + 1:N).v});
x(1:N - 1) = spm_unvec(qu.x{1},{M(1:N - 1).x});


% precisions
%==========================================================================
for i = 1:N

    % precision of causal and hidden states
    %----------------------------------------------------------------------
    try
        h{i,1} = spm_vec(feval(M(i).ph,x{i},v{i},qh.h{i},M(i)));
    catch
        h{i,1} = sparse(M(i).l,1);
    end
    try
        g{i,1} = spm_vec(feval(M(i).pg,x{i},v{i},qh.g{i},M(i)));
    catch
        g{i,1} = sparse(M(i).n,1);
    end

end

% Concatenate over hierarchical levels
%--------------------------------------------------------------------------
p.h  = spm_cat(h);
p.g  = spm_cat(g);

if nargout < 2, return, end

% gradients
%==========================================================================

% assume precisions can be functions of hyper-parameters and states
%--------------------------------------------------------------------------
try method.h = M(1).E.method.h; catch, method.h = 1; end
try method.g = M(1).E.method.g; catch, method.g = 1; end
try method.x = M(1).E.method.x; catch, method.x = 1; end
try method.v = M(1).E.method.v; catch, method.v = 1; end

% number of variables
%--------------------------------------------------------------------------
nx      = numel(spm_vec(x));
nv      = numel(spm_vec(v));
hn      = numel(spm_vec(qh.h));
gn      = numel(spm_vec(qh.g));
nh      = size(p.h,1);
ng      = size(p.g,1);

dp.h.dh = sparse(nh,hn);
dp.g.dg = sparse(ng,gn);
dp.h.dx = sparse(nh,nx);
dp.h.dv = sparse(nh,nv);
dp.g.dx = sparse(ng,nx);
dp.g.dv = sparse(ng,nv);


% gradients w.r.t. h only (no state-dependent noise)
%----------------------------------------------------------------------
if method.h || method.g

    for i = 1:N

        % precision of causal and hidden states
        %--------------------------------------------------------------
        dhdh{i,i} = spm_diff(M(i).ph,x{i},v{i},qh.h{i},M(i),3);
        dgdg{i,i} = spm_diff(M(i).pg,x{i},v{i},qh.g{i},M(i),3);

    end

    % Concatenate over hierarchical levels
    %------------------------------------------------------------------
    dp.h.dh = spm_cat(dhdh);
    dp.g.dg = spm_cat(dgdg);

end


% gradients w.r.t. causal states
%----------------------------------------------------------------------
if method.v

    for i = 1:N

        % precision of causal states
        %--------------------------------------------------------------
        dhdv{i,i} = spm_diff(M(i).ph,x{i},v{i},qh.h{i},M(i),2);

        % precision of hidden states
        %--------------------------------------------------------------
        dgdv{i,i} = spm_diff(M(i).pg,x{i},v{i},qh.g{i},M(i),2);

    end

    % Concatenate over hierarchical levels
    %------------------------------------------------------------------
    dp.h.dv = spm_cat(dhdv);
    dp.g.dv = spm_cat(dgdv);

end

% gradients w.r.t. hidden states
%----------------------------------------------------------------------
if method.x

    for i = 1:N

        % precision of causal states
        %--------------------------------------------------------------
        dhdx{i,i} = spm_diff(M(i).ph,x{i},v{i},qh.h{i},M(i),1);

        % precision of hidden states
        %--------------------------------------------------------------
        dgdx{i,i} = spm_diff(M(i).pg,x{i},v{i},qh.g{i},M(i),1);

    end

    % Concatenate over hierarchical levels
    %------------------------------------------------------------------
    dp.h.dx = spm_cat(dhdx);
    dp.g.dx = spm_cat(dgdx);

end



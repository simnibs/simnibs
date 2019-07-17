function [y] = spm_int_B(P,M,U)
% Integrate a MIMO nonlinear system using a bilinear Jacobian
% FORMAT [y] = spm_int_B(P,M,U)
% P   - model parameters
% M   - model structure
% U   - input structure or matrix
%
% y   - (v x l)  response y = g(x,u,P)
%__________________________________________________________________________
%
% Integrates the MIMO system described by
%
%    dx/dt = f(x,u,P,M)
%    y     = g(x,u,P,M)
%
% using the update scheme:
%
%    x(t + dt) = x(t) + U*dx(t)/dt
%
%            U = (expm(dt*J) - I)*inv(J)
%            J = df/dx
%
% at input times.  This integration scheme evaluates the update matrix (Q)
% at each time point
%
%--------------------------------------------------------------------------
%
% SPM solvers or integrators
%
% spm_int_ode:  uses ode45 (or ode113) which are one and multi-step solvers
% respectively.  They can be used for any ODEs, where the Jacobian is
% unknown or difficult to compute; however, they may be slow.
%
% spm_int_J: uses an explicit Jacobian-based update scheme that preserves
% nonlinearities in the ODE: dx = (expm(dt*J) - I)*inv(J)*f.  If the
% equations of motion return J = df/dx, it will be used; otherwise it is
% evaluated numerically, using spm_diff at each time point.  This scheme is
% infallible but potentially slow, if the Jacobian is not available (calls
% spm_dx).
%
% spm_int_E: As for spm_int_J but uses the eigensystem of J(x(0)) to eschew 
% matrix exponentials and inversion during the integration. It is probably
% the best compromise, if the Jacobian is not available explicitly.
%
% spm_int_B: As for spm_int_J but uses a first-order approximation to J
% based on J(x(t)) = J(x(0)) + dJdx*x(t).
%
% spm_int_L: As for spm_int_B but uses J(x(0)).
%
% spm_int_U: like spm_int_J but only evaluates J when the input changes.
% This can be useful if input changes are sparse (e.g., boxcar functions).
% It is used primarily for integrating EEG models
%
% spm_int:   Fast integrator that uses a bilinear approximation to the
% Jacobian evaluated using spm_bireduce. This routine will also allow for
% sparse sampling of the solution and delays in observing outputs. It is
% used primarily for integrating fMRI models
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_int_B.m 5219 2013-01-29 17:07:07Z spm $
 
 
% convert U to U.u if necessary
%--------------------------------------------------------------------------
if ~isstruct(U), U.u = U; end
try, dt = U.dt; catch, dt = 1; end
 
% state equation; add [0] states if not specified
%--------------------------------------------------------------------------
try
    f   = fcnchk(M.f,'x','u','P','M');
catch
    f   = inline('sparse(0,1)','x','u','P','M');
    M.n = 0;
    M.x = sparse(0,0);
end
 
% and output nonlinearity
%--------------------------------------------------------------------------
try
    g   = fcnchk(M.g,'x','u','P','M');
end
 
% Initial states and inputs
%--------------------------------------------------------------------------
u = U.u(1,:);
try
    x   = M.x;
catch
    x   = sparse(0,1);
    M.x = x;
end

% check for delay operator
%--------------------------------------------------------------------------
try
    [fx,dfdx,D] = f(x,u,P,M);
catch
    D = 1;
end
 
% get Jacobian and its derivatives
%--------------------------------------------------------------------------
[dJdx,J]  = spm_diff(f,x,u,P,M,[1 1]);
[dJdu,J]  = spm_diff(f,x,u,P,M,[1 2]);
 
% integrate
%==========================================================================
x0    = spm_vec(M.x);
for i = 1:size(U.u,1)
 
    % input
    %----------------------------------------------------------------------
    u     = U.u(i,:);
 
    % motion
    %----------------------------------------------------------------------
    fx    = f(x,u,P,M);
 
    % dx(t)/dt and Jacobian df/dx
    %----------------------------------------------------------------------
    dx    = spm_vec(x) - x0;
    dfdx  = J;
    for j = 1:length(dJdx)
        dfdx = dfdx + dJdx{j}*dx(j);
    end
    for j = 1:length(dJdu)
        dfdx = dfdx + dJdu{j}*u(j);
    end
 
    % update dx = (expm(dt*J) - I)*inv(J)*fx
    %----------------------------------------------------------------------
    x  = spm_unvec(spm_vec(x) + spm_dx(D*dfdx,D*fx,dt),x);
 
    % output - implement g(x)
    %----------------------------------------------------------------------
    try
        y(:,i) = spm_vec(g(x,u,P,M));
    catch
        y(:,i) = spm_vec(x);
    end
 
end
 
% transpose
%--------------------------------------------------------------------------
y      = real(y');

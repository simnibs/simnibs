function y = cat_stat_polynomial(x,p)
% Polynomial expansion and orthogonalization of function x
% FORMAT y = cat_stat_polynomial(x,p)
% x   - data matrix
% p   - order of polynomial [default: 1]
% 
% y   - orthogonalized data matrix
%__________________________________________________________________________
%
% cat_stat_polynomial orthogonalizes a polynomial function of order p
%__________________________________________________________________________
% Christian Gaser 
% $Id: cat_stat_polynomial.m 1187 2017-09-15 14:06:21Z gaser $ 

if nargin < 2, p = 1; end

if size(x,1) < size(x,2)
    x = x';
end

y = spm_detrend(x(:));
v = zeros(size(y,1),p + 1);

for j = 0:p
    v(:,(j + 1)) = (y.^j) - v*(pinv(v)*(y.^j));
end

for j = 2:p
    y = [y v(:,(j + 1))];
end

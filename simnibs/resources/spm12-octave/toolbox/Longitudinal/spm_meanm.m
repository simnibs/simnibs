function M = spm_meanm(A)
% Compute barycentre of matrix exponentials
% FORMAT M = spm_meanm(A)
% A - A 3D array, where each slice is a matrix
% M - the resulting mean
%
% Note that matrices should not be too dissimilar to each other or the
% procedure fails.
% See http://hal.archives-ouvertes.fr/hal-00699361/
%_______________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_meanm.m 7408 2018-08-24 14:54:57Z john $

N = size(A,3);
M = eye(size(A,1),size(A,2));

for iter = 1:1024
    S = zeros(size(M));
    for i=1:N
        L = real(logm(M\A(:,:,i)));
        S = S + L;
    end
    S = S/N;
    M = M*expm(S);
    %imagesc(M); drawnow
    %fprintf('%d\t%g\n', iter,sum(S(:).^2));
    if sum(S(:).^2)<1e-20
        break;
    end
end


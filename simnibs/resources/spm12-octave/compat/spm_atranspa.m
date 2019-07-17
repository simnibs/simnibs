function C = spm_atranspa(A)
% Multiplies the transpose of a matrix by itself
% FORMAT C = spm_atranspa(A)
% A - real matrix
% C - real symmetric matrix resulting from A'A
%_______________________________________________________________________
%
% This compiled routine was written to save both memory and CPU time but
% is now deprecated. Use A'*A directly instead.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_atranspa.m 4418 2011-08-03 12:00:13Z guillaume $

persistent runonce
if isempty(runonce)
    warning('spm_atranspa is deprecated. Use A''*A instead.');
    runonce = 1;
end

C = A'*A;

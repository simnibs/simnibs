function varagout = cat_display_matlab_PID
% cat_display_matlab_PID (in development!)
% ______________________________________________________________________
%
% Display exor return PID of this (Linux/Mac) or the last started 
% (Windows) MATLAB instace.
% ______________________________________________________________________
%
%   Robert Dahnke - robert.dahnke@uni-jena.de
%   Center of Neuroimaging 
%   Department of Psychiatry and Psychotherapy 
%   University Hostpital Jena
% ______________________________________________________________________
% $Id: cat_display_matlab_PID.m 921 2016-04-08 14:57:33Z dahnke $

  % get PID
  pid = feature('getpid'); 
  
  % display PID
  if nargout==0
    if isnumeric(pid) && ~isempty(pid)
      fprintf('CAT parallel processing with MATLAB PID: %d\n',pid);
    else
      fprintf('CAT parallel processing with MATLAB PID: unknown %s PID %d\n',computer,pids);
    end
  else 
    varagout{1} = pid; 
  end
  
end

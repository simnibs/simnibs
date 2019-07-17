function [ST, RS] = cat_system(varargin)
% ______________________________________________________________________
% CAT12 wrapper for system calls
% This is necessary because windows does not allow spaces in system
% calls. Thus, we have to cd into that folder and call the command
% from this folder.
% ______________________________________________________________________
% Christian Gaser
% $Id: cat_system.m 1186 2017-09-14 14:58:48Z gaser $

rev = '$Rev: 1186 $';

if nargin == 0
  error('Argument is missing');
end

CATDir      = fullfile(spm('dir'),'toolbox','cat12','CAT');
if ispc
  CATDir = [CATDir '.w32'];
elseif ismac
  CATDir = [CATDir '.maci64'];
elseif isunix
  CATDir = [CATDir '.glnx86'];
end  

if ispc
  olddir = pwd;
  cd(CATDir);
  [ST, RS] = system(varargin{1});
  cd(olddir);
else
  cmd = fullfile(CATDir,varargin{1});
  warning off % this is to prevent warnings by calling cat12 from the shell script
  [ST, RS] = system(cmd);
  warning on
end

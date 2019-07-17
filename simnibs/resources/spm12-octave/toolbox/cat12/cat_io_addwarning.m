function warn = cat_io_addwarning(warn,id,mess,nline)
% Fuction to add warnings to a structure in cat_main. 
% ______________________________________________________________________
%
%   Robert Dahnke (robert.dahnke@uni-jena.de)
%   Structural Brain Mapping Group (http://dbm.neuro.uni-jena.de/)
%   Department of Neurology
%   University Jena
% ______________________________________________________________________
% $Id: cat_io_addwarning.m 861 2016-02-12 16:49:37Z dahnke $

  warn(end+1) = struct('identifier',id,'message',mess);
  warnstr = strrep(mess,'\\n','\n'); 
  warnstr = strrep(warnstr,'\n','\n         '); 
  if exist('nline','var') && nline, fprintf('\n'); end
  cat_io_cmd(sprintf(['WARNING: ' warnstr]),'warning');
return
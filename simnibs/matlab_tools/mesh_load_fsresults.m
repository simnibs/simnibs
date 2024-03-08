function m=mesh_load_fsresults(fnameIn, varargin)
% reads simulation results that were mapped on the GM surface
% (i.e, results in 'subject_overlays' or 'fsavg_overlays' folders)
% and stores them as node_data field in a mesh
%
% The function will automatically locate the surfaces belonging
% to the results unless specified otherwise (see options).
% 
% USAGE:
%   m=mesh_load_fsresults(fname_data [,'OptionName',OptionValue,...])
% 
%   fname_data: file name of the results (can either be lh. or rh.)
%
%   Options are set using the option name followed by the value.
%   Supported options:
%   addtomesh: add data to given mesh-structure, instead of loading surface
%   hemi: 'lh', 'rh' or 'both' (standard: 'both')
%
%   Output:
%   m: mesh structure with data added to node_data field
%
%
%   Note: wrapper function around read_curv.m from FreeSurfer
%
% A. Thielscher, 06-Nov-2017; updated 21-Sep-2018

%    This program is part of the SimNIBS package.
%    Please check on www.simnibs.org how to cite our work in publications.
%
%    Copyright (C) 2017-2018 Axel Thielscher
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% standard settings
s.hemi = 'both';
s.addtomesh = [];

% parse input
if nargin<1; error('file name needed as input'); end
s=parse_input(s,varargin{:});

load_lh = false; load_rh = false;
if strcmpi(s.hemi,'lh')||strcmpi(s.hemi,'both'); load_lh=true; end
if strcmpi(s.hemi,'rh')||strcmpi(s.hemi,'both'); load_rh=true; end

% get full path, filename without 'rh.' or 'lh.', and name of data field
if ~exist(fnameIn,'file'); error(['could not find ' fnameIn]); end

[~,fnameFull,~]=fileattrib(fnameIn);
fnameFull=fnameFull.Name;
idx=find(fnameFull==filesep,2,'last');
pnameFull=fnameFull(1:idx(2)-1);
pnamePart=fnameFull(idx(1)+1:idx(2)-1);

[~,fnameStripped,extHlp]=fileparts(fnameIn);
fnameStripped=[fnameStripped extHlp];
idx=find(fnameStripped=='.',1,'first');
if ~strcmpi(fnameStripped(1:idx-1),'lh')&&~strcmpi(fnameStripped(1:idx-1),'rh')
    error('filename has to start with lh. or rh.');
end
fnameStripped=fnameStripped(idx+1:end);

idx=find(fnameStripped=='.',2,'last');
fieldname=fnameStripped(idx(1)+1:end);
if strcmpi(fieldname,'E.magn')
    fieldname='magnE';
elseif strcmpi(fieldname,'J.magn')
    fieldname='magnJ';
end

% load surface(s) if needed
if isempty(s.addtomesh)

   if strcmpi(pnamePart,'fsavg_overlays')
       m=mesh_load_fssurf('fsaverage', 'hemi', s.hemi);
       
   elseif strcmpi(pnamePart,'subject_overlays')
       m=mesh_empty;
       if load_lh; m=mesh_load_fssurf([pnameFull filesep 'lh.central']); end
       if load_rh
           m2=mesh_load_fssurf([pnameFull filesep 'rh.central']);
           idx_firstrh=size(m.nodes,1)+1;
           m.nodes = [m.nodes; m2.nodes];
           m.triangles= [m.triangles; m2.triangles+idx_firstrh-1];
           m.triangle_regions= [m.triangle_regions; 2*m2.triangle_regions];
       end
   else
       error('folder structure unclear');
   end
   
else
    m=s.addtomesh;
end

% load data
dataidx=length(m.node_data)+1;
m.node_data{dataidx}.name=fieldname;
m.node_data{dataidx}.data=[];
disp(['adding data as node_data nr. ' num2str(dataidx)]);
if load_lh
    m.node_data{dataidx}.data=read_curv([pnameFull filesep 'lh.' fnameStripped]);
end
if load_rh
     m.node_data{dataidx}.data =[m.node_data{dataidx}.data; ...
                                 read_curv([pnameFull filesep 'rh.' fnameStripped])];
end

end



function [curv, fnum] = read_curv(fname)
%
% [curv, fnum] = read_curv(fname)
% reads a binary curvature file into a vector
%


%
% read_curv.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.4 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


%fid = fopen(fname, 'r') ;
%nvertices = fscanf(fid, '%d', 1);
%all = fscanf(fid, '%d %f %f %f %f\n', [5, nvertices]) ;
%curv = all(5, :)' ;

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
	 str = sprintf('could not open curvature file %s', fname) ;
	 error(str) ;
end
vnum = fread3(fid) ;
NEW_VERSION_MAGIC_NUMBER = 16777215;
if (vnum == NEW_VERSION_MAGIC_NUMBER)
	 vnum = fread(fid, 1, 'int32') ;
	 fnum = fread(fid, 1, 'int32') ;
	 vals_per_vertex = fread(fid, 1, 'int32') ;
   curv = fread(fid, vnum, 'float') ; 
	   	
  fclose(fid) ;
else

	fnum = fread3(fid) ;
  curv = fread(fid, vnum, 'int16') ./ 100 ; 
  fclose(fid) ;
end
%nvertices = fscanf(fid, '%d', 1);
%all = fscanf(fid, '%d %f %f %f %f\n', [5, nvertices]) ;
%curv = all(5, :)' ;
end


function [retval] = fread3(fid)
% [retval] = fd3(fid)
% read a 3 byte integer out of a file


%
% fread3.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
retval = bitshift(b1, 16) + bitshift(b2,8) + b3 ;
end


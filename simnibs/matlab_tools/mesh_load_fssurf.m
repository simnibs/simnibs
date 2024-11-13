function [m, varargout]=mesh_load_fssurf(fnameIn, varargin)
% loads triangle surfaces, and stores them in a mesh structure
%
% USAGE:
%   m=mesh_load_fssurf(fname [,'OptionName',OptionValue,...])
%
% fname: Can be one of the three:
%        1) 'fsaverage': will load the gray matter of the fsaverage
%              template
%        2) m2m_{subID} path: will load the central gray matter surface
%              of the subject specified by the m2m_... path
%        3) surface filename, including path
%
%   Options are set using the option name followed by the value.
%   Supported options:
%   hemi: 'lh', 'rh' or 'both' (standard: 'both')
%         (will be ignored when giving a surface filename)
%   label: 'a2009s', 'DK40' or 'HCP_MMP1'
%          automatically adds the labels of the stated atlas as node_data
%          to the mesh; works ONLY when 'fsaverage' is given as fname;
%          The label names are given as second output, e.g.
%           [m names]=mesh_load_fssurf('fsaverage','label','a2009s')
%
%           'a2009s': Destrieux atlas (FreeSurfer v4.5, aparc.a2009s)
%           Cite: Destrieux, C. Fischl, B. Dale, A., Halgren, E. A sulcal 
%           depth-based anatomical parcellation of the cerebral cortex. 
%           Human Brain Mapping (HBM) Congress 2009, Poster #541
%
%           'DK40': Desikan-Killiany atlas (FreeSurfer, aparc.a2005s)
%           Cite: Desikan RS, S�gonne F, Fischl B, Quinn BT, Dickerson BC,
%           Blacker D, Buckner RL, Dale AM, Maguire RP, Hyman BT, Albert MS,
%           Killiany RJ. An automated labeling system for subdividing the
%           human cerebral cortex on MRI scans into gyral based regions of
%           interest. Neuroimage. 2006 Jul 1;31(3):968-80. 
%
%           'HCP_MMP1': Human Connectome Project (HCP) Multi-Modal Parcellation
%           Cite: Glasser MF, Coalson TS, Robinson EC, et al. A multi-modal 
%           parcellation of human cerebral cortex. Nature. 2016;536(7615):171-178.
% 
%   Output:
%   m: mesh structure
%
%  Notes: This functions wraps read_surf.m from FreeSurfer and gifti from spm12
%  to read FS and CAT12 surfaces and store them in a mesh structure
%  automatically adds 1 to face indices of FS surfaces, so that they start with 1
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
s.label = '';

% parse input
if nargin<1; error('file or path name needed as input'); end
s=parse_input(s,varargin{:});

% determine what to load
path_to_avg_surf = fullfile(SIMNIBSDIR, 'resources', 'templates', 'fsaverage_surf');
path_to_labels = fullfile(SIMNIBSDIR, 'resources', 'templates', 'fsaverage_atlases');

load_lh = false; load_rh = false;
if strcmpi(s.hemi,'lh')||strcmpi(s.hemi,'both'); load_lh=true; end
if strcmpi(s.hemi,'rh')||strcmpi(s.hemi,'both'); load_rh=true; end
lh={}; rh={};

if strcmpi(fnameIn,'fsaverage')
    % look up fsaverage template in cat12 folder
    if load_lh; lh={ fullfile(path_to_avg_surf,'lh.central.freesurfer.gii') }; end
    if load_rh; rh={ fullfile(path_to_avg_surf,'rh.central.freesurfer.gii') }; end
    
elseif exist(fnameIn,'dir')
    % load subject-specific surfaces
    if exist(fullfile(fnameIn,'charm_log.html'),'file')
        if load_lh; lh={ fullfile(fnameIn,'surfaces','lh.central.gii') }; end
        if load_rh; rh={ fullfile(fnameIn,'surfaces','rh.central.gii') }; end
    else
        error(['No .._log.html found in ' fnameIn '. Unclear whether it was created by charm']);
    end
    
elseif exist(fnameIn,'file')
    % load given surface
    lh={ fnameIn };
else
    error([fnameIn ' not found']);
end

for i=1:length(lh)
    if ~exist(lh{i},'file'); error(['could not find ' lh{i}]); end
end

for i=1:length(rh)
    if ~exist(rh{i},'file'); error(['could not find ' rh{i}]); end
end

% load surfaces
m=mesh_empty;

if ~isempty(lh)
    [nodes, triangles] = load_surface(lh{1});
    m.nodes = [m.nodes; nodes];
    m.triangles= [m.triangles; triangles];
    m.triangle_regions= [m.triangle_regions; ones(size(triangles,1),1)];
end
if length(lh)>1
   [nodes, ~] = load_surface(lh{2});
   m.nodes=(m.nodes+nodes)/2; 
end
    
if ~isempty(rh) 
    idx_firstrh=size(m.nodes,1)+1;
    
    [nodes, triangles] = load_surface(rh{1});
    m.nodes = [m.nodes; nodes];
    m.triangles= [m.triangles; triangles+idx_firstrh-1];
    m.triangle_regions= [m.triangle_regions; 2*ones(size(triangles,1),1)];
end
if length(rh)>1
   [nodes, ~] = load_surface(rh{2});
   m.nodes(idx_firstrh:end,:)=(m.nodes(idx_firstrh:end,:)+nodes)/2; 
end

m.nodes=double(m.nodes);
m.triangles=double(m.triangles);

if ~isempty(s.label)
    if strcmpi(s.label,'a2009s')
        fname_label=fullfile(path_to_labels,'lh.aparc_a2009s.annot');
    elseif strcmpi(s.label,'DK40')
        fname_label=fullfile(path_to_labels,'lh.aparc_DK40.annot');
    elseif strcmpi(s.label,'HCP_MMP1')
        fname_label=fullfile(path_to_labels,'lh.aparc_HCP_MMP1.annot');
    else
        error(['unknown label file: ' s.label]);
    end
    
    [m, varargout{1}]=mesh_load_fsannot(m,fname_label);
end

end

function [nodes, triangles] = load_surface(fname)
    [~,~,extHlp] = fileparts(fname);
    if strcmpi(extHlp,'.gii') % load gifti
        s=gifti(fname);
        nodes=s.vertices;
        triangles=s.faces;
    else
        [nodes, triangles] = read_surf(fname);
        triangles=triangles+1; % FS indexing starts at 0
    end
end


function [vertex_coords, faces, magic] = read_surf(fname)
%
% [vertex_coords, faces] = read_surf(fname)
% reads a the vertex coordinates and face lists from a surface file
% note that reading the faces from a quad file can take a very long
% time due to the goofy format that they are stored in. If the faces
% output variable is not specified, they will not be read so it 
% should execute pretty quickly.
%


%
% read_surf.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: fischl $
%    $Date: 2014/04/30 12:59:03 $
%    $Revision: 1.7 $
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


%QUAD_FILE_MAGIC_NUMBER =  (-1 & 0x00ffffff) ;
%NEW_QUAD_FILE_MAGIC_NUMBER =  (-3 & 0x00ffffff) ;

TRIANGLE_FILE_MAGIC_NUMBER =  16777214 ;
QUAD_FILE_MAGIC_NUMBER =  16777215 ;
NEW_QUAD_FILE_MAGIC_NUMBER =  16777213 ;

fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
  str = sprintf('could not open surface file %s.', fname) ;
  error(str) ;
end
magic = fread3(fid) ;

if((magic == QUAD_FILE_MAGIC_NUMBER) | (magic == NEW_QUAD_FILE_MAGIC_NUMBER))
  vnum = fread3(fid) ;
  fnum = fread3(fid) ;
  vertex_coords = fread(fid, vnum*3, 'int16') ./ 100 ; 
  if (nargout > 1)
    for i=1:fnum
      for n=1:4
	faces(i,n) = fread3(fid) ;
      end
    end
  end
elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER)
  fgets(fid) ;
  fgets(fid) ;
  vnum = fread(fid, 1, 'int32') ;
  fnum = fread(fid, 1, 'int32') ;
  vertex_coords = fread(fid, vnum*3, 'float32') ; 
  if (nargout > 1)
    faces = fread(fid, fnum*3, 'int32') ;
    faces = reshape(faces, 3, fnum)' ;
  end
else
  fprintf('ERROR: magic number %d unknown\n',magic);
  vertex_coords = [];
  faces = [];
  return;
end

vertex_coords = reshape(vertex_coords, 3, vnum)' ;
fclose(fid) ;
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

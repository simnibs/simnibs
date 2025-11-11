function [m, struct_names]=mesh_load_fsannot(m,fnameIn)
%
% [m, struct_names]=mesh_load_fsannot(m, fname_annot)
%
% read annotation or label data and store it as node_data field in the mesh
% 
% Input:
%   m: mesh to which the data will be added to
%   fname_annot: filename including path. Either lh. or rh. can be given,
%                the other will be located automatically
%
% Output:
%   m: mesh, with annotation data added as node_data field
%   struct_names: names of the anatomical structures
%
% Note: wrapper function around read_annotation.m and readlabel.m from
% FreeSurfer
%
% A. Thielscher, 06-Nov-2017; updated 05-Oct-2018

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

% parse input
if ~nargin; error('mesh has to be given as first input'); end

% split into path- and filename, remove initial lh. or rh. from filename
if nargin<2||isempty(fnameIn)
    [fnameIn,pnameIn] = uigetfile({'*.annot;*.label'},'Select an annotation or label file');
    if isequal(fnameIn,0) || isequal(pnameIn,0); return; end
    fnameIn = fullfile(pnameIn,fnameIn);
end

[pnameIn,fnameHlp,extHlp] = fileparts(fnameIn);
fnameHlp=[fnameHlp extHlp];

idx=find(fnameHlp=='.',1,'first');
if isempty(idx); error('filename has to start with lh. or rh.'); end
fnamePrefix=fnameHlp(1:idx-1);
if ~strcmpi(fnamePrefix,'lh')&&~strcmpi(fnamePrefix,'rh')
    error('filename has to start with lh. or rh.');
end
fnameIn=fnameHlp(idx+1:end);

% determine hemispheres from mesh triangle regions and get corresponding
% file names
region_idx=unique(m.triangle_regions);
load_lh = false; load_rh = false; start_idx=[]; fnames={}; prefix={};
if any(region_idx == 1)
    surf_idx = 1;
    both_hemispheres = true;
elseif any(region_idx == 1002)
    surf_idx = 1002;
    both_hemispheres = true;
end
if both_hemispheres
    load_lh=true;
    load_rh=true;
     % get first lh node
    hlpIdx=m.triangles(m.triangle_regions==surf_idx,:);
    start_idx=[start_idx min(hlpIdx(:))]; 
    fnames={fnames{:} [pnameIn filesep 'lh.' fnameIn]};
    prefix={prefix{:} 'lh.'};
end
if any(region_idx == 2)
    load_rh=true;
    hlpIdx=m.triangles(m.triangle_regions==2,:); % get first rh node
    start_idx=[start_idx min(hlpIdx(:))];
    fnames={fnames{:} [pnameIn filesep 'rh.' fnameIn]};
    prefix={prefix{:} 'rh.'};
% This here handles when both hemispheres are together as a .msh file
elseif (surf_idx == 1002)
    fnames={fnames{:} [pnameIn filesep 'rh.' fnameIn]};
    prefix={prefix{:} 'rh.'};
end
if ~load_lh&&~load_rh; error('mesh has to contain triangle regions 1 or 2'); end
for i=1:length(fnames)
    if ~exist(fnames{i},'file'); error(['could not find ' fnames{i}]); end
end

% get index of node_data field
dataidx=length(m.node_data)+1;
m.node_data{dataidx}.name=fnameIn;
m.node_data{dataidx}.data=zeros(size(m.nodes,1),1);
disp(['adding annotations as node_data nr. ' num2str(dataidx) ' (name: ' fnameIn ')']);

struct_names={};
for i=1:length(fnames)
    N_struct_name=length(struct_names);
    
    if strcmpi(extHlp, '.label')
        % add vertices from label file to node_data field
        l = read_label('',fnames{i});
        vertices=l(:,1)+start_idx(i); % FS indexing starts at 0
        m.node_data{dataidx}.data(vertices)=1+N_struct_name;

        % add structure name
        struct_names{1+N_struct_name}=[prefix{i} fnameIn];
    else
        % add labels from annotation file to node_data field
        [vertices, label, colortable] = read_annotation(fnames{i});
        vertices=int32(vertices)+start_idx(i); % FS indexing starts at 0
        label_remap=zeros(size(label)); % remap labeling to be 1,2,3,4,...
        for j=1:length(colortable.table(:,5))
            label_remap(label==colortable.table(j,5))=j+N_struct_name;
        end
        m.node_data{dataidx}.data(vertices)=label_remap;
        
        % get structure names
        for j=1:length(colortable.struct_names)
            struct_names{N_struct_name+j}=[prefix{i} colortable.struct_names{j}];
        end
    end
    % This is an ugly hack dot the loading to work with the CS meshes
    if surf_idx == 1002
        start_idx=[start_idx max(vertices) + 1];
    end
end
end


function [l] = read_label(sname, lname)
% l = read_label(<sname>, lname)
%
% reads the label file 'lname' from the subject 'sname' 
% in the subject's label directory into the vector l
% l will be nvertices-by-5, where each column means:
% (1) vertex number, (2-4) xyz at each vertex, (5) stat
%
% IMPORTANT: the vertex number is 0-based.
% 


%
% read_label.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
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


l = [];

if(nargin ~= 2)
  fprintf('l = read_label(<sname>, lname)\n');
  return;
end

if(~isempty(sname))
  sdir = getenv('SUBJECTS_DIR') ;
  fname = sprintf('%s/%s/label/%s.label', sdir, sname, lname) ;
else
	ind = findstr(lname, '.label') ;
	if (length(ind) > 0)
        fname = lname ;
    else
		fname = sprintf('%s.label', lname);
	end
end

% open it as an ascii file
fid = fopen(fname, 'r') ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

fgets(fid) ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

line = fgets(fid) ;
nv = sscanf(line, '%d') ;
l = fscanf(fid, '%d %f %f %f %f\n') ;
l = reshape(l, 5, nv) ;
l = l' ;

fclose(fid) ;
end


function [vertices, label, colortable] = read_annotation(filename, varargin)
%
% NAME
%
%       function [vertices, label, colortable] = ...
%                                       read_annotation(filename [, verbosity])
%
% ARGUMENTS
% INPUT
%       filename        string          name of annotation file to read
%
% OPTIONAL
%       verbosity       int             if true (>0), disp running output
%                                       + if false (==0), be quiet and do not
%                                       + display any running output
%
% OUTPUT
%       vertices        vector          vector with values running from 0 to
%                                       + size(vertices)-1
%       label           vector          lookup of annotation values for 
%                                       + corresponding vertex index.
%       colortable      struct          structure of annotation data
%                                       + see below
%       
% DESCRIPTION
%
%       This function essentially reads in a FreeSurfer annotation file
%       <filename> and returns structures and vectors that together 
%       assign each index in the surface vector to one of several 
%       structure names.
%       
% COLORTABLE STRUCTURE
% 
%       Consists of the following fields:
%       o numEntries:   number of entries
%       o orig_tab:     filename of original colortable file
%       o struct_names: cell array of structure names
%       o table:        n x 5 matrix
%                       Columns 1,2,3 are RGB values for struct color
%                       Column 4 is a flag (usually 0)
%                       Column 5 is the structure ID, calculated from
%                       R + G*2^8 + B*2^16 + flag*2^24
%                       
% LABEL VECTOR
% 
%       Each component of the <label> vector has a structureID value. To
%       match the structureID value with a structure name, lookup the row
%       index of the structureID in the 5th column of the colortable.table
%       matrix. Use this index as an offset into the struct_names field
%       to match the structureID with a string name.      
%
% PRECONDITIONS
%
%       o <filename> must be a valid FreeSurfer annotation file.
%       
% POSTCONDITIONS
%
%       o <colortable> will be an empty struct if not embedded in a
%         FreeSurfer annotation file. 
%       
% EXAMPLE
% [vertices label ctab] = read_annotation(fname);
% stgctab = strmatch('superiortemporal',char(ctab.struct_names));
% stgcode = ctab.table(stgctab,5);
% indstg = find(label==stgcode);
% nstg = length(indstg);

%
% read_annotation.m
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2014/02/25 19:54:10 $
%    $Revision: 1.10 $
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

fp = fopen(filename, 'r', 'b');

verbosity = 0;
if length(varargin)
    verbosity       = varargin{1};  
end;

if(fp < 0)
   if verbosity, disp('Annotation file cannot be opened'); end;
   return;
end

A = fread(fp, 1, 'int');

tmp = fread(fp, 2*A, 'int');
vertices = tmp(1:2:end);
label = tmp(2:2:end);

bool = fread(fp, 1, 'int');
if(isempty(bool)) %means no colortable
   if verbosity, disp('No Colortable found.'); end;
   colortable = struct([]);
   fclose(fp);
   return; 
end

if(bool)
    
    %Read colortable
    numEntries = fread(fp, 1, 'int');

    if(numEntries > 0)
        
        if verbosity, disp(['Reading from Original Version']); end;
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        for i = 1:numEntries
            len = fread(fp, 1, 'int');
            colortable.struct_names{i} = fread(fp, len, '*char')';
            colortable.struct_names{i} = colortable.struct_names{i}(1:end-1);
            colortable.table(i,1) = fread(fp, 1, 'int');
            colortable.table(i,2) = fread(fp, 1, 'int');
            colortable.table(i,3) = fread(fp, 1, 'int');
            colortable.table(i,4) = fread(fp, 1, 'int');
            colortable.table(i,5) = colortable.table(i,1) + colortable.table(i,2)*2^8 + colortable.table(i,3)*2^16 + colortable.table(i,4)*2^24;
        end
        if verbosity
            disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
        end
    else
        version = -numEntries;
        if verbosity
          if(version~=2)    
            disp(['Error! Does not handle version ' num2str(version)]);
          else
            disp(['Reading from version ' num2str(version)]);
          end
        end
        numEntries = fread(fp, 1, 'int');
        colortable.numEntries = numEntries;
        len = fread(fp, 1, 'int');
        colortable.orig_tab = fread(fp, len, '*char')';
        colortable.orig_tab = colortable.orig_tab(1:end-1);
        
        colortable.struct_names = cell(numEntries,1);
        colortable.table = zeros(numEntries,5);
        numEntriesToRead = fread(fp, 1, 'int');
        for i = 1:numEntriesToRead
            structure = fread(fp, 1, 'int')+1;
            if (structure < 0)
              if verbosity, disp(['Error! Read entry, index ' num2str(structure)]); end;
            end
            if(~isempty(colortable.struct_names{structure}))
              if verbosity, disp(['Error! Duplicate Structure ' num2str(structure)]); end;
            end
            len = fread(fp, 1, 'int');
            colortable.struct_names{structure} = fread(fp, len, '*char')';
            colortable.struct_names{structure} = colortable.struct_names{structure}(1:end-1);
            colortable.table(structure,1) = fread(fp, 1, 'int');
            colortable.table(structure,2) = fread(fp, 1, 'int');
            colortable.table(structure,3) = fread(fp, 1, 'int');
            colortable.table(structure,4) = fread(fp, 1, 'int');
            colortable.table(structure,5) = colortable.table(structure,1) + colortable.table(structure,2)*2^8 + colortable.table(structure,3)*2^16 + colortable.table(structure,4)*2^24;       
	end
        if verbosity 
          disp(['colortable with ' num2str(colortable.numEntries) ' entries read (originally ' colortable.orig_tab ')']);
        end
    end    
else
    if verbosity
        disp('Error! Should not be expecting bool = 0');    
    end;
end

fclose(fp);

% This makes it so that each empty entry at least has a string, even
% if it is an empty string. This can happen with average subjects.
for i = 1:numEntries
  if(isempty(colortable.struct_names{i}))
    colortable.struct_names{i}='';
  end
end

return;
end

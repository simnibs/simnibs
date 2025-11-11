function [m]=mesh_extract_regions(m, varargin)
% extracts the triangles and/or tetrahedra with the given region number
% deletes unused nodes, elment_data entries and updates the node numbers
%
% USAGE:
%   m=mesh_extract_regions(m, [,'OptionName',OptionValue,...]) 
%
%   m: input mesh
%  
%   Options are set using the option name followed by the value.
%   elemtype: determines the element type to extract
%             ('tri', 'tet' or 'both'; standard: 'both')
%   region_idx: extracts the elements with the given region numbers
%   keepAllNodes: when set to true, do not remove unused nodes and 
%                 do not update node numbers (standard: false)
%   node_idx: indices of the nodes to keep (standard: keep all)
%   tri_idx: indices of the triangles to keep
%   tet_idx: indices of the tetrahedra to keep
%
% Examples:
%  m=mesh_extract_regions(m, 'elemtype','tri'); % delete all tetrahedra
%  m=mesh_extract_regions(m, 'elemtype','tet'); % delete all triangles
%  m=mesh_extract_regions(m, 'region_idx', [5 1005]);
%             % keep tetrahedra and triangles with region numbers 5 or 1005
%         
% A. Thielscher 11-Apr-2018, based on prior code from M. Windhoff
% A.Thielscher: updated 03-Oct-2018, added 'node_idx', 'tri_idx' and
%               'tet_idx' options

%    This program is part of the SimNIBS package.
%    Please check on www.simnibs.org how to cite our work in publications.
%
%    Copyright (C) 2018 Axel Thielscher, M. Windhodff
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
s.elemtype='both';
s.region_idx=[];
s.keepAllNodes=false;
s.node_idx=[];
s.tri_idx=[];
s.tet_idx=[];

% parse input
if nargin<1
    error('at least one argument is needed');
end
s=parse_input(s,varargin{:});

keepTri=false;
keepTet=false;
if strcmpi(s.elemtype,'tri')
    keepTri=true;
elseif strcmpi(s.elemtype,'tet')
    keepTet=true;
elseif strcmpi(s.elemtype,'both')
    keepTri=true;
    keepTet=true;
else
    error(['unclear elemtype ' s.elemtype]);
end

% determine what to keep
idx_kept_nodes=false(size(m.nodes,1),1);
idx_kept_tri=false(size(m.triangles,1),1);
idx_kept_tet=false(size(m.tetrahedra,1),1);

if keepTri && size(m.triangles,1) > 0
    if isempty(s.region_idx)
        idx_kept_tri(:)=true;
    else
        for i=1:length(s.region_idx)
             idx_kept_tri = idx_kept_tri | m.triangle_regions==s.region_idx(i);
        end
    end
    
    if ~isempty(s.tri_idx)
        if islogical(s.tri_idx)&&length(s.tri_idx)~=size(m.triangles,1)
            error('logical index must have same length as number of triangles')
        end
        if ~islogical(s.tri_idx)
            idx_hlp=false(size(m.triangles,1),1);
            idx_hlp(s.tri_idx)=true;
            s.tri_idx=idx_hlp;
        end
        idx_kept_tri = idx_kept_tri&s.tri_idx;
    end
end

if keepTet && size(m.tetrahedra,1) > 0
    if isempty(s.region_idx)
        idx_kept_tet(:)=true;
    else
        for i=1:length(s.region_idx)
             idx_kept_tet = idx_kept_tet | m.tetrahedron_regions==s.region_idx(i);
        end
    end
    
    if ~isempty(s.tet_idx)
        if islogical(s.tet_idx)&&length(s.tet_idx)~=size(m.tetrahedra,1)
            error('logical index must have same length as number of tets')
        end
        if ~islogical(s.tet_idx)
            idx_hlp=false(size(m.tetrahedra,1),1);
            idx_hlp(s.tet_idx)=true;
            s.tet_idx=idx_hlp;
        end
        idx_kept_tet = idx_kept_tet&s.tet_idx;
    end
end

if s.keepAllNodes
    if ~isempty(s.node_idx); error('keepAllNodes and node_idx cannot be used together'); end
    idx_kept_nodes(:)=true;
    disp('keeping all nodes');
else
    idx_kept_nodes(unique(m.triangles(idx_kept_tri,:)))=true;
    idx_kept_nodes(unique(m.tetrahedra(idx_kept_tet,:)))=true;
    
    if ~isempty(s.node_idx)
        if islogical(s.node_idx)&&length(s.node_idx)~=size(m.nodes,1)
            error('logical index must have same length as number of nodes')
        end
        if ~islogical(s.node_idx)
            idx_hlp=false(size(m.nodes,1),1);
            idx_hlp(s.node_idx)=true;
            s.node_idx=idx_hlp;
        end
        idx_kept_nodes = idx_kept_nodes&s.node_idx;
        
        % delete affected triangles and tets
        idx_kept_tri=idx_kept_tri&(sum(idx_kept_nodes(m.triangles),2)==3);
        idx_kept_tet=idx_kept_tet&(sum(idx_kept_nodes(m.tetrahedra),2)==4);
    end
end

% update node indices, delete unused nodes, triangles and tetrahedra
ab_map=zeros(size(m.nodes,1),1);
ab_map(idx_kept_nodes)=1:sum(idx_kept_nodes);

m.nodes=m.nodes(idx_kept_nodes,:);

m.triangles=ab_map(m.triangles(idx_kept_tri,:));
m.triangle_regions=m.triangle_regions(idx_kept_tri);

m.tetrahedra=ab_map(m.tetrahedra(idx_kept_tet,:));
m.tetrahedron_regions=m.tetrahedron_regions(idx_kept_tet);


% update node and element data
for i=1:length(m.node_data)
    m.node_data{i}.data=m.node_data{i}.data(idx_kept_nodes,:);
end

for i=1:length(m.element_data)
    if ~isempty(m.element_data{i}.tridata)
        m.element_data{i}.tridata=m.element_data{i}.tridata(idx_kept_tri,:);
    end
    
    if ~isempty(m.element_data{i}.tetdata)
        m.element_data{i}.tetdata=m.element_data{i}.tetdata(idx_kept_tet,:);
    end
end

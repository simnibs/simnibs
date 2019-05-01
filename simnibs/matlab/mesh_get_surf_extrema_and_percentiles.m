function [results, varargout] = mesh_get_surf_max_and_percentiles(m,data_index,varargin)
%
% results = mesh_get_surf_max_and_percentiles(m,data_index {, percentiles })
% [results results_annot] = mesh_get_surf_max_and_percentiles(m,data_index, percentiles, annot_index, struct_names)
% 
% get maximum and some percentiles of field data that was 
% read in by mesh_load_fsresults.m
%
% m: mesh with data
% data_index: index of node_data containing the field data
% percentiles: vectors listing the percentiles that will be determined
%              (standard: [50 75 90 95 99 99.5])
% annot_index: index of node_data containing annotation labels
% struct_names: list of the names of the anatomical structures
% 
% A. Thielscher, 06-Nov-2017

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



disp(['extracting data from ' m.node_data{data_index}.name ' (index ' num2str(data_index) ')']);

percentiles=[50 75 90 95 99 99.5];
if nargin>2&&~isempty(varargin{1})
    percentiles=varargin{1};
end

annot_index=[];
if nargin>3&&~isempty(varargin{2})
    annot_index=varargin{2};
    
    if nargin<4||isempty(varargin{3}); error('struct_names has to be given as input when annotation labels are used'); end
    struct_names=varargin{3};
    
    disp(['using annotations from ' m.node_data{annot_index}.name ' (index ' num2str(annot_index) ')']);
end

% get node areas (for each node, this is the average of the areas of the
% surrounding triangles)
triareas=mesh_get_triangle_sizes(m);
nodeareas=zeros(size(m.nodes,1),1);
for i=1:size(m.triangles,1)
    nodeareas(m.triangles(i,1))=nodeareas(m.triangles(i,1))+triareas(i)/3;
    nodeareas(m.triangles(i,2))=nodeareas(m.triangles(i,2))+triareas(i)/3;
    nodeareas(m.triangles(i,3))=nodeareas(m.triangles(i,3))+triareas(i)/3;
end

nodedata=m.node_data{data_index}.data;


% get results for overall surface
results(1).name='overall';
results(1).fieldname=m.node_data{data_index}.name;
results(1).percentiles=percentiles;
[results(1).max, results(1).min, results(1).perc_values, results(1).perc_area] = get_max_and_perc(nodedata, nodeareas, percentiles);


% get results for lh and rh separately
if any(m.triangle_regions==2) % when lh and rh are stored in the mesh, then rh has region number 2
    % get index of first rh node
    idx_startrh=m.triangles(m.triangle_regions==2,:); 
    idx_startrh=min(idx_startrh(:));
    
    % results for lh
    results(2).name='lh';
    results(2).fieldname=m.node_data{data_index}.name;
    results(2).percentiles=percentiles;
    [results(2).max, results(2).min, results(2).perc_values, results(2).perc_area] = ...
            get_max_and_perc(nodedata(1:idx_startrh-1), nodeareas(1:idx_startrh-1), percentiles);
    
    % results for rh
    results(3).name='rh';
    results(3).fieldname=m.node_data{data_index}.name;
    results(3).percentiles=percentiles;
    [results(3).max, results(3).min, results(3).perc_values, results(3).perc_area] = ...
            get_max_and_perc(nodedata(idx_startrh:end), nodeareas(idx_startrh:end), percentiles);
end
results = orderfields(results, [1 2 5 4 3 6 7]);


% get results for annotated areas
if ~isempty(annot_index)
    for i=1:length(struct_names)
        res_annot(i).name=struct_names{i};
        res_annot(i).fieldname=m.node_data{data_index}.name;
        res_annot(i).percentiles=percentiles; 
        
        idx=m.node_data{annot_index}.data==i;
        if any(idx)
            [res_annot(i).max, res_annot(i).min, res_annot(i).perc_values, res_annot(i).perc_area] = ...
                get_max_and_perc(nodedata(idx), nodeareas(idx), percentiles);
        end
    end
    
    varargout{1} = res_annot;
end

end


function [ maxval, minval, percvals, percareas] = get_max_and_perc( nodedata, nodeareas, percentiles )

idx_notnan=~isnan(nodedata);
nodedata=nodedata(idx_notnan);
nodeareas=nodeareas(idx_notnan);

maxval=max(nodedata);
minval=min(nodedata);

[nodedata,idx] = sort(nodedata);
nodeareas=nodeareas(idx);
nodeareas=cumsum(nodeareas);
nodeareasNorm=nodeareas/nodeareas(end); % normalize to one

for i=1:length(percentiles)
    idx=find(nodeareasNorm>percentiles(i)/100,1,'first');
    percvals(i)=nodedata(idx);
    percareas(i)=nodeareas(idx);
end

percareas=nodeareas(end)-percareas;

end

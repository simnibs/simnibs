function [areas,mean_edge_sizes]=mesh_get_triangle_sizes(m,varargin)
  % Calculates the areas of the triangles for the regions given in
  % varargin. Optional output: Mean edge sizes for each triangle.
  % USAGE:
  % [areas,mean_edge_sizes]=MESH_GET_TRIANGLE_SIZES(m[, regions])
  % where varargin=regions e.g. [1,2,3]
  % Mirko Windhoff, 2009
  
  % triangle area = 1/2*norm((cross(a,b)), where a,b edge vectors

%    This program is part of the SimNIBS package.
%    Please check on www.simnibs.org how to cite our work in publications.
%
%    Copyright (C) 2009-2018 Mirko Windhoff
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


  if nargin>1
    regions_values=unique(m.triangle_regions)';
    regions=varargin{1};
    regions=regions_values(regions); % map numbers
    % get affected tetrahedra, all regions tetrahedra in ascending region
    % order
    disp(['Affected regions' sprintf(' %d,',regions)]);
    triangles=m.triangles(~any(abs(repmat(m.triangle_regions,1,size(regions,2))-repmat(regions,size(m.triangles,1),1)),2),:);
  else
    triangles=m.triangles;
  end
  nt=size(triangles,1);
  % dim(edges)=[3 cords,2 edges, nt]
  edges=reshape(m.nodes(reshape(triangles(:,2:3)',nt*2,1),:)',3,2,nt)-repmat(reshape(m.nodes(triangles(:,1),:)',3,1,nt),[1,2,1]);
  areas=squeeze(1/2*sqrt(sum(cross(edges(:,1,:),edges(:,2,:),1).^2,1)));
  if nargout>1
    % dim(edges)=[3 coords, 3 edges, nt]
    edges=cat(2,edges,edges(:,1,:)-edges(:,2,:));
    mean_edge_sizes=squeeze(mean(sqrt(sum(edges.^2,1)),2));
  end
end

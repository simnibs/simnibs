function [volumes,mean_edge_sizes]=mesh_get_tetrahedron_sizes(m)
% Calculates the volumes of the tetrahedra for the regions given in
% varargin. Optional output: mean edge sizes for each tetrahedron.
%
% USAGE:
% [volumes,mean_edge_sizes]=mesh_get_tetrahedron_sizes(m)
%
% Mirko Windhoff, 2009
% AT 11-Apr-2018: simplified i/o code
 
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

 

nt=size(m.tetrahedra,1);
% dim(edges)=[3 coords,3 edges, nt]
edges=reshape(m.nodes(reshape(m.tetrahedra(:,2:4)',nt*3,1),:)',3,3,nt)-repmat(reshape(m.nodes(m.tetrahedra(:,1),:)',3,1,nt),[1,3,1]);
volumes=squeeze(1/6*sqrt(sum(dot(cross(edges(:,1,:),edges(:,2,:),1),edges(:,3,:),1).^2,1)));

if nargout>1
    % dim(edges)=[3 coords, 6 edges, nt]
    edges=cat(2,edges,cat(2,edges(:,1,:)-edges(:,2,:),edges(:,1,:)-edges(:,3,:),edges(:,2,:)-edges(:,3,:)));
    mean_edge_sizes=squeeze(mean(sqrt(sum(edges.^2,1)),2));
end;

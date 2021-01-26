function [centers]=mesh_get_triangle_centers(m)
%
% returns the center of mass of the triangles
%
% M. Windhoff, 2009
centers=(m.nodes(m.triangles(:,1),:) + ...
    m.nodes(m.triangles(:,2),:) + ...
    m.nodes(m.triangles(:,3),:))./3;
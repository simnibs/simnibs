function areas = mesh_get_node_areas(m)
% Calculates the areas of the nodes, defined as 1/3 of the area of all
% triangles containing said node
  % USAGE:
  %  areas =MESH_GET_NODE_AREAS(m[)
  % m: mesh structure

  % Guilherme Saturnino 2019

tr_areas = mesh_get_triangle_sizes(m);
tr_flat = reshape(m.triangles', length(m.triangles)*3, 1);
areas = accumarray(...
    tr_flat, repelem(tr_areas, 3), ...
    [length(m.nodes), 1]) * (1./3.);
end

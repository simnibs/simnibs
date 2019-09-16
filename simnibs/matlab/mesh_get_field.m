function field = mesh_get_field(mesh, field_name)
% Returns only the selected field from the mesh
% Usage:
%   mesh: mesh structure
%   field name: name of field to be selected

% Returns:
%   field: structure with the field named field_name

occurences = 0;
field.name = field_name;

for i=1:length(mesh.node_data)
    if strcmp(mesh.node_data{i}.name, field_name)
        occurences = occurences + 1;
        field = mesh.node_data{i};
    end
end

for i=1:length(mesh.element_data)
    if strcmp(mesh.element_data{i}.name, field_name)
        occurences = occurences + 1;
        field = mesh.element_data{i};
    end
end

if occurences == 0
    error(['Could not find the field ' field_name ' in the mesh']);
elseif occurences > 1
    warning(['Found more than one occurance of the field ' field_name]);
end
end
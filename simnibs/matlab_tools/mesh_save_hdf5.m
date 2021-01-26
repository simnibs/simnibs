function mesh_save_hdf5(m, fn, group_name)
% Saves a mesh as an hdf5 file

% mesh_save_hdf5(m, fn, group_name);
% m: mesh structure
% fn: name of HDF5 file
% group_name: Name of the group where to save the mesh. Example: '/mesh'

if isfile(fn)
    error(['File ' fn ' exists'])
end
fid = H5F.create(fn);
plist = 'H5P_DEFAULT';
% Create datasets
m_group = H5G.create(fid, group_name, plist, plist, plist);
n_group = H5G.create(m_group, 'nodes', plist, plist, plist);
H5G.close(n_group)
e_group = H5G.create(m_group, 'elm', plist, plist, plist);
H5G.close(e_group)
ed_group = H5G.create(m_group, 'elmdata', plist, plist, plist);
H5G.close(ed_group)
nd_group = H5G.create(m_group, 'nodedata', plist, plist, plist);
H5G.close(nd_group)

H5G.close(m_group)
H5F.close(fid)

% Write nodes
create_and_write(fn, [group_name '/nodes/node_coord'], m.nodes')
create_and_write(fn, [group_name '/nodes/node_number'], int32(1:length(m.nodes)))

%Write elements
node_number_list = [m.triangles, -1*ones(size(m.triangles, 1), 1); m.tetrahedra];
create_and_write(fn, [group_name '/elm/node_number_list'], node_number_list')
elm_number = int32(1:length(node_number_list));
create_and_write(fn, [group_name '/elm/elm_number'], elm_number)
elm_type = int32([2*ones(size(m.triangle_regions)); 4*ones(size(m.tetrahedron_regions))]);
create_and_write(fn, [group_name '/elm/elm_type'], elm_type')
tag = int32([m.triangle_regions; m.tetrahedron_regions]);
create_and_write(fn, [group_name '/elm/tag1'], tag')
create_and_write(fn, [group_name '/elm/tag2'], tag')

% Write node data
for i=1:length(m.node_data)
    create_and_write(fn, [group_name '/nodedata/' m.node_data{i}.name], m.node_data{i}.data') 
end


% Write element data
for i=1:length(m.element_data)
    if isfield(m.element_data{i}, 'tridata') && isfield(m.element_data{i}, 'tetdata')
        data=[m.element_data{i}.tridata; m.element_data{i}.tetdata];
    elseif isfield(m.element_data{i}, 'data')
        data= m.element_data{i}.data;
    else
        error('element_data needs a tridata/tetdata pair or a data field')
    end
    create_and_write(fn, [group_name '/elmdata/' m.element_data{i}.name], data') 
end


end

function create_and_write(fn, datasetname, value)
h5create(fn, datasetname, size(value), 'Datatype', class(value))
h5write(fn, datasetname, value)
end
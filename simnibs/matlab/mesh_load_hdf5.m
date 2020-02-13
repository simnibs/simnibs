function meshes = mesh_load_hdf5(varargin)
% Load meshes stored in hdf5 format
%
% meshes=mesh_load_hdf5(fn); Reads all the mesh structures in the root of
% an hdf5 file named 'fn' and returns a struct with all mesh structures
%
% mesh=mesh_load_hdf5(fn, group_name); Reads the mesh structure in the
% group_name and returns a single mesh structure

% Guilherme Saturnino, 2018

if nargin==2
    fn = varargin{1};
    group_name = varargin{2};
    meshes = read_group_mesh(fn, group_name);
elseif nargin==1
    fn = varargin{1};
    meshes = struct('name', {}, 'mesh', {});
    % Read the mesh in the root of the file, if any
    if is_mesh_struct(fn, '/')
        meshes(end+1).name = '/';
        meshes(end).mesh = read_group_mesh(fn, '/'); 
    end
    % Scan the groups in the root and find the meshes
    ifo = h5info(fn);
    for i = 1:length(ifo.Groups)
        g_name = ifo.Groups(i).Name;
        if is_mesh_struct(fn, g_name)
            meshes(end+1).name = g_name;
            meshes(end).mesh = read_group_mesh(fn, g_name);
        end
    end
    if isempty(meshes)
        error('Could not find any meshes in HDF5 file')
    end
else
    error('Invalid number of arguments')
end


end

function v = is_mesh_struct(fn, group_name)
    v = true;
    required_fields = {
        '/nodes/node_coord',...
        '/elm/node_number_list',...
        '/elm/tag1', '/elm/elm_type'};
    for i = 1:length(required_fields)
        try
            h5read(fn, [group_name, required_fields{i}]);
        catch
            v = false;
        return
        end
    end
end

function m = read_group_mesh(fn, group_name)
    % Initialize structure
    m.node_data = {};
    m.element_data = {};
    m.nodes = zeros(0, 3, 'double');
    m.triangles = zeros(0, 3, 'int32');
    m.triangle_regions = zeros(0, 1, 'int32');
    m.tetrahedra = zeros(0, 4, 'int32');
    m.tetrahedron_regions = zeros(0, 1, 'int32');
    % Read nodes
    m.nodes = h5read(fn, [group_name, '/nodes/node_coord'])';
    % Read elements
    elm_types = h5read(fn, [group_name, '/elm/elm_type'])';
    elm_tags = h5read(fn, [group_name, '/elm/tag1'])';
    node_number_list = h5read(fn, [group_name, '/elm/node_number_list'])';
    m.triangles = int32(node_number_list(elm_types==2, 1:3));
    m.triangle_regions = elm_tags(elm_types==2);
    m.tetrahedra =  int32(node_number_list(elm_types==4, 1:4));
    m.tetrahedron_regions = elm_tags(elm_types==4);
    % read NodeData
    ifo = h5info(fn, [group_name, '/nodedata']);
    for i = 1:length(ifo.Datasets)
        dset = ifo.Datasets(i);
        m.node_data{i}.name = dset.Name;
        m.node_data{i}.data =  double(h5read(fn, [group_name, '/nodedata/' dset.Name])');
    end
    % read ElementData
    ifo = h5info(fn, [group_name, '/elmdata']);
    for i = 1:length(ifo.Datasets)
        dset = ifo.Datasets(i);
        m.node_data{i}.name = dset.Name;
        m.element_data{i}.name = dset.Name;
        data =  h5read(fn, [group_name, '/elmdata/' dset.Name])';
        m.element_data{i}.tridata = double(data(elm_types==2, 1:end));
        m.element_data{i}.tetdata = double(data(elm_types==4, 1:end));
    end
end
function field_values = get_fields_at_coordinates(mesh, coords, out_fill)
% Evaluates the fields of the mesh in a given set of coordinates
% This function calls the command line tool "get_fields_at_coordinates". It has
% a large overhead.

% mesh: mesh structure or mesh file name. Must contain volumetric elements!
% coords: Coordinates where to evaluate the field, in Nx3 format
% out_fill: value to assign for coordinates outside the mesh, optional.
%           nan: NaN, default
%           nearest: Use the nearest value inside the mesh
% Returns: A cell with each field in the mesh evaluated at each given coordinate

% Guilherme B Saturnino, 2019
assert(size(coords, 2) == 3, 'coords must be in Nx3 format');
if nargin < 3
    out_fill = 'nan';
end

if ~any(strcmp(out_fill, {'nan', 'nearest'}))
    error('out_fill must be nan, or nearest')
end

if ischar(mesh)
    assert(exist(mesh, 'file') == 2, ['Could not find mesh file ' mesh])
    fn_mesh = mesh;
    is_temp = false; 
else
    is_temp = true;
    fn_mesh = [tempname '.msh'];
    mesh_save_gmsh4(mesh, fn_mesh)
end

fn_in = [tempname,'.csv'];
csvwrite(fn_in, coords);

% Run mni2subject_coords
[status,result] = system([simnibs_cli_call('get_fields_at_coordinates')...
                         ' -m ' fn_mesh ' -s ' fn_in, ' --out_fill ' out_fill]);
if status ~= 0
    if is_temp
        delete(fn_mesh);
    end
    delete(fn_in);
    error('There was an error running get_fields_at_coordinates:\n %s',result)
end

field_values = {};
out_names = dir([fn_in(1:end-4) '_*.csv']);
[tmp_dir, fname_main] = fileparts(fn_in);
for i = 1:length(out_names)
    match = textscan(out_names(i).name, [fname_main '_%s']);
    field.name = match{1}{1}(1:end-4); % of course matlab returns weird stuff
    field.data = csvread([tmp_dir filesep out_names(i).name]);
    field_values{end + 1} = field;
    delete([tmp_dir filesep out_names(i).name])
end

if is_temp
    delete(fn_mesh);
end
delete(fn_in);

end
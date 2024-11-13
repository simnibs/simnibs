function [data, d_type] = mesh_load_data_matrix_entry(fn, entry, data_path)
% Reads a single entry in a data matrix

% [data, d_type] = mesh_load_leadfield_entry(fn, entry, data_path);
% fn: name of hdf5 file with leadfield
% entry: entries in the leadfield to be read
% data_path: path to the the dataset within the hdf5 file

% Returns:
%  data: data in the leadfield entry
%  d_type: either 'element_data', 'node_data' or 'uknown'
%
% Guilherme Saturnino, 2018

ifo = h5info(fn, data_path);
dset_size = ifo.Dataspace.Size;
st = ones(size(dset_size));
st(end) = entry;
ct = dset_size;
ct(end) = 1;
data = h5read(fn, data_path, st, ct)';

try
    d_type = h5readatt(fn,data_path,'definition');
    d_type = d_type{1};
catch
    d_type = 'unknown';
end
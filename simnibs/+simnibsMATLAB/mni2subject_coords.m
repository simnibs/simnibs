function coords_sub = mni2subject_coords(coords_mni, subdir, transformation_type)
% Transforms a set of coordinates in MNI space to subject space
% This function calls the command line tool "mni2subject_coords", so it has
% a large overhead.
% coords_mni: MNI coordinates, in Nx3 format
% subdir: path to subject directory, example: 'm2m_ernie/'
% transformation_type (optional): type of trasfomation to use,
%                                 can be {'nonl' (non-linear),
%                                '12dof', '6dof' (6 and 12 degrees of freedom)}
% Returns: coordinates in subject space
% Guilherme B Saturnino, 2019

if nargin < 3
    transformation_type = 'nonl';
end
if ~any(strcmp(transformation_type, {'nonl', '12dof', '6dof'}))
    error('transformation_type must be nonl, 12dof, or 6dof')
end
assert(size(coords_mni, 2) == 3, 'coords_subject must be in Nx3 format');
assert(exist(subdir, 'dir') == 7, ['Could not find directory ' subdir])

% Write CSV file in SimNIBS format
fn_in = [tempname,'.csv'];
fid = fopen( fn_in, 'w' );
for i = 1 : size(coords_mni, 1)
    fprintf(...
        fid, 'Generic, %f, %f, %f\n', ...
        coords_mni(i,1), coords_mni(i, 2),  coords_mni(i, 3));
end
fclose( fid );
fn_out = [tempname,'.csv'];

% Run mni2subject_coords
[status,result] = system([simnibs_cli_call('mni2subject_coords')...
                         ' -m "' subdir '" -s ' fn_in ' -o ' fn_out ...
                         ' -t ' transformation_type]);
% Check if call was successefull
if status ~= 0
    delete(fn_in);
    delete(fn_out);
    error('There was an error running mni2subject_coords:\n %s',result)
end

% Read output
coords_sub = csvread(fn_out, 0, 1);
delete(fn_in);
delete(fn_out);

end

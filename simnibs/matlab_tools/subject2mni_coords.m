function coords_mni = subject2mni_coords(coords_subject, subdir, transformation_type)
% Transforms a set of coordinates in subject space to MNI space
% This function calls the command line tool "subject2mni_coords", so it has
% a large overhead.
% coords_subject: Subject coordinates, in Nx3 format
% subdir: path to subject directory, example: 'm2m_ernie/'
% transformation_type (optional): type of trasfomation to use,
%                                 can be {'nonl' (non-linear),
%                                '12dof', '6dof' (6 and 12 degrees of freedom)}
% Returns: coordinates in MNI space
% Guilherme B Saturnino, 2019

if nargin < 3
    transformation_type = 'nonl';
end
if ~any(strcmp(transformation_type, {'nonl', '12dof', '6dof'}))
    error('transformation_type must be nonl, 12dof, or 6dof')
end
assert(size(coords_subject, 2) == 3, 'coords_subject must be in Nx3 format');
assert(exist(subdir, 'dir')==7, ['Could not find directory ' subdir])

% Write CSV file in SimNIBS format
fn_in = [tempname,'.csv'];
fid = fopen( fn_in, 'w' );
for i = 1 : size(coords_subject, 1)
    fprintf(...
        fid, 'Generic, %f, %f, %f\n', ...
        coords_subject(i,1), coords_subject(i, 2),  coords_subject(i, 3));
end
fclose( fid );
fn_out = [tempname,'.csv'];

% Run subject2mni_coords
[status,result] = system([simnibs_cli_call('subject2mni_coords')...
                         ' -m "' subdir '" -s ' fn_in ' -o ' fn_out ...
                         ' -t ' transformation_type]);
% Check if call was successefull
if status ~= 0
    delete(fn_in);
    delete(fn_out);
    error('There was an error running subject2mni_cords:\n %s',result)
end

% Read output
coords_mni = csvread(fn_out, 0, 1);
delete(fn_in);
delete(fn_out);

end

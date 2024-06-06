function matsimnibs = mni2subject_coilpos(center_MNI, ydir_MNI, zdir_MNI, ...
                                          subdir, coil_skin_distance, ...
                                          transformation_type)
% Transforms a coil position in MNI space to subject space
%
% NOTE: Coil centers are projected on the head surface as standard. Set
% coil_skin_dist to ensure a specific distance of the coil center
% to the head surface.
%
% Coil centers away from the head are best defined by providing the
% closest positon on the surface of the MNI head and the desired
% coil-skin-distance.
%
% Non-linear transformations (used as default) are accurate for
% positions in the head and on its surface, but not outside the
% head. Thus, coil centers are automatically projected on the head
% surface during transformation.
% This function calls the command line tool "mni2subject_coords", so
% it has a large overhead.
%
% INPUTS:
%   center_MNI: coil center in MNI space (vector of size 1x3)
%   ydir_MNI: unit vector (size 1x3) indicating the y-direction
%   zdir_MNI: unit vector (size 1x3) indicating the z-direction
%   subdir: path to subject directory, example: 'm2m_ernie/'
%   coil_skin_distance (optional): skin-coil distance in [mm]. 
%                                  (standard: 0)
%   transformation_type (optional): type of transfomation to use, can be 
%                                   'nonl' (non-linear), '12dof' or '6dof' 
%                                   (6 and 12 degrees of freedom);
%                                   (standard: 'nonl')
%
% RETURNS:
%   coil position as matsimnibs (size 4x4) in subject space
%
% Axel Thielscher, 2024
if nargin < 5
    coil_skin_distance = 0;
end
if nargin < 6
    transformation_type = 'nonl';
end

assert(all(size(center_MNI) == [1, 3]) | all(size(center_MNI) == [3, 1]), ...
       'center_MNI must be a 1x3 vector');
assert(all(size(ydir_MNI) == [1, 3]) | all(size(ydir_MNI) == [3, 1]), ...
       'ydir_MNI must be a 1x3 vector');
assert(all(size(zdir_MNI) == [1, 3]) | all(size(zdir_MNI) == [3, 1]), ...
       'zdir_MNI must be a 1x3 vector');
assert(exist(subdir, 'dir') == 7, ['Could not find directory ' subdir])
assert(isscalar(coil_skin_distance), 'coil_skin_distance must be a scalar') 
if ~any(strcmp(transformation_type, {'nonl', '12dof', '6dof'}))
    error('transformation_type must be nonl, 12dof, or 6dof')
end

% Write CSV file in SimNIBS format
% Type, pos_x, pos_y, pos_z, ez_x, ez_y, ez_z, ey_x, ey_y, ey_z, dist, name
fn_in = [tempname,'.csv'];
fid = fopen( fn_in, 'w' );
fprintf( ...
        fid, 'CoilPos, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, pos1\n', ...
        center_MNI(1), center_MNI(2),  center_MNI(3), ...
        zdir_MNI(1), zdir_MNI(2), zdir_MNI(3), ...
        ydir_MNI(1), ydir_MNI(2), ydir_MNI(3), ...
        coil_skin_distance ...
        );
fclose( fid );
fn_out = tempname;
fn_geo = [fn_out,'.geo'];
fn_out = [fn_out,'.csv'];

% Run mni2subject_coords
[status,result] = system([simnibs_cli_call('mni2subject_coords')...
                         ' -m "' subdir '" -s ' fn_in ' -o ' fn_out ...
                         ' -t ' transformation_type]);
% Check if call was successful
if status ~= 0
    delete(fn_in);
    delete(fn_out);
    delete(fn_geo);
    error('There was an error running mni2subject_coords:\n %s',result)
end

% Read output
coilpos = readmatrix(fn_out,'range',[1 2 1 10]);
delete(fn_in);
delete(fn_out);
delete(fn_geo);

% make matsimnibs
matsimnibs = eye(4);
matsimnibs(1:3,4) = coilpos(1:3);
matsimnibs(1:3,3) = coilpos(4:6);
matsimnibs(1:3,2) = coilpos(7:9);
matsimnibs(1:3,1) = cross(coilpos(7:9),coilpos(4:6));
assert(det(matsimnibs) > 0, 'matsimnibs is not a right-handed coordinate system')

end

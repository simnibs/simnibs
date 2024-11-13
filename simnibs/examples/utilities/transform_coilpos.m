%
% Transform a coil position from MNI space to subject space
% 
% run within the directory containing the m2m_ernie example datatset
%
% Copyright (c) 2019 SimNIBS developers. Licensed under the GPL v3.

% coil center in MNI space (here: above left handknob)
center_MNI = [-45.8, -19.1,  84.];
% unit vectors for the y- and z-axes of the coil in MNI space
ydir_MNI = [0.379, 0.923, 0.062];
zdir_MNI = [0.546, -0.169, -0.827];

% coil-skin distance in [mm]; 
% (optional; standard: 0, i.e. the coil center is projected 
%  on the skin of the subject)
coil_skin_distance = 4;

matsimnibs = mni2subject_coilpos(center_MNI, ydir_MNI, zdir_MNI, ...
                                 'm2m_ernie', ...
                                 coil_skin_distance);

disp(matsimnibs)
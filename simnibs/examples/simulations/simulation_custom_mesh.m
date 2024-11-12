%   Example on running SimNIBS simulations with a custom mesh
% 
%   NOTE: This example requires the mesh "myspheres.msh"
%   Please see "How to create and use a custom mesh"
%   in the SimNIBS tutorials for instructions to create the mesh
%     
%   Copyright (c) 2021 SimNIBS developers. Licensed under the GPL v3.


S = sim_struct('SESSION');
S.fnamehead = 'myspheres.msh'; % name of custom mesh
S.pathfem = 'simu_custom_mesh'; % Folder for the simulation output
S.open_in_gmsh = 1;
% Note: As there is no m2m_{subID} folder, postprocessing
%       options are not available.


% add a TDCS simulation
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = [1e-3, -1e-3]; % Current going through each channel, in Ampere

% 'myspheres.msh' contains a custom tissue with label number 17.
% We need to assign a conductivity to this tissue label.
S.poslist{1}.cond(17).value = 2; % in S/m
S.poslist{1}.cond(17).name = 'custom_tissue';

% define first electrode
S.poslist{1}.electrode(1).channelnr = 1;  % Connect the first electrode to the first channel
S.poslist{1}.electrode(1).centre = [10, 50, 50];  % position determined from the nifti file
S.poslist{1}.electrode(1).shape = 'rect';  % Rectangular shape
S.poslist{1}.electrode(1).dimensions = [50, 50];  % 50x50 mm
S.poslist{1}.electrode(1).thickness = 4;  % 4 mm thickness

% define second electrode
S.poslist{1}.electrode(2).channelnr = 2;
S.poslist{1}.electrode(2).centre = [90, 50, 50];
S.poslist{1}.electrode(2).shape = 'ellipse';  % Circular shape
S.poslist{1}.electrode(2).dimensions = [50, 50];  % 50 mm diameter
S.poslist{1}.electrode(2).thickness = 4; % 4 mm thickness


% add a TMS simulation
S.poslist{2} = sim_struct('TMSLIST');
S.poslist{2}.fnamecoil = fullfile('legacy_and_other','Magstim_70mm_Fig8.ccd');  % Choose a coil model

S.poslist{2}.cond(17).value = 2; % in S/m
S.poslist{2}.cond(17).name = 'custom_tissue';

% Define the coil position
S.poslist{2}.pos(1).centre = [50, 50, 90];
S.poslist{2}.pos(1).pos_ydir = [50, 40, 90]; % Polongation of coil handle (see documentation)
S.poslist{2}.pos(1).distance = 4; % 4 mm distance from coil surface to head surface

% Run simulation
run_simnibs(S)


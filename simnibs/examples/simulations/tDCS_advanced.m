%
% example script that runs a simnibs tDCS simulation to demonstrate 
% several post-processing options and the modeling of more 
% complex electrode shapes
% 
% Copyright (c) 2021 SimNIBS developers. Licensed under the GPL v3.
%

%% General information
S = sim_struct('SESSION');
S.subpath = 'm2m_ernie'; % subject folder
S.pathfem = 'tdcs_advanced_simu'; % output folder

S.fields = 'eEjJ'; % Save the following results:
                   %  e: Electric field magnitude
                   %  E: Electric field vector
                   %  j: Current density magnitude
                   %  J: Current density vector
                   
S.map_to_surf = true;   %  Map to subject's middle gray matter surface
S.map_to_fsavg = true;  %  Map to FreeSurfer's FSAverage group template
S.map_to_vol = true;    %  Save as nifti volume
S.map_to_MNI = true;    %  Save in MNI space
S.tissues_in_niftis = [1,2,3]; % Results in the niftis will be masked 
                               % to only show WM (1), GM (2), CSF(3)
                               % (standard: only GM)
                               % To get fields everywhere: 
                               %    S.tissues_in_niftis = 'all'
                               
S.open_in_gmsh = true; % show results in gmsh (not for the the niftis)


%% Define tDCS simulation
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = [0.001, -0.001]; % Current flow through each channel (in Ampere)

% First electrode
S.poslist{1}.electrode(1).channelnr = 1; % Connect the electrode to the first channel
S.poslist{1}.electrode(1).centre = 'C3'; % Place it over C3
S.poslist{1}.electrode(1).pos_ydir = 'C4'; % Position on the scalp to define 
                                           % the electrode orientation. The electrode's 
                                           % y-axis will point from the centre to pos_ydir 
S.poslist{1}.electrode(1).shape = 'rect'; % Rectangular electrode
                                          % Other shapes: 
                                          %  'ellipse' (includes circular electrodes)
                                          %  'custom': custom shape. In this case, define
                                          %  a 2D shape using .vertices
S.poslist{1}.electrode(1).dimensions = [50, 50]; % Dimension in mm
S.poslist{1}.electrode(1).thickness = 4; % 4 mm thickness
                                         %  1 number: electrode with 1 (gel) layer.
                                         %  2 numbers: electrode with 2 layers
                                         %  (gel, conductive rubber) 
                                         %  3 numbers: electrode with 3 layers
                                         %  (sponge, conductive rubber, sponge) 
                                         %  also .dimensions_sponge must be 
                                         %  then set

% Second electrode
S.poslist{1}.electrode(2).channelnr = 1;
S.poslist{1}.electrode(2).centre = 'C4';
S.poslist{1}.electrode(2).pos_ydir = 'C3';
S.poslist{1}.electrode(2).shape = 'rect';
S.poslist{1}.electrode(2).dimensions = [40, 40];
S.poslist{1}.electrode(2).thickness = [3.5, 1, 3.5]; % sponge electrode
S.poslist{1}.electrode(2).dimensions_sponge = [50, 70];
% define where the cable is plugged in:
S.poslist{1}.electrode(2).plug = sim_struct('ELECTRODE');
S.poslist{1}.electrode(2).plug.centre = [-10, 0]; % relative to electrode center
S.poslist{1}.electrode(2).plug.shape = 'rect';
S.poslist{1}.electrode(2).plug.dimensions = [10 10];

% Third electrode
S.poslist{1}.electrode(3).channelnr = 2; % also connected to 2nd channel
S.poslist{1}.electrode(3).centre = 'Oz';
S.poslist{1}.electrode(3).pos_ydir = 'Pz';
S.poslist{1}.electrode(3).shape = 'custom';
S.poslist{1}.electrode(3).vertices = [ 0, 20; ...
                                     -40, 40; ...
                                     -20, 0; ...
                                     -40, -40; ...
                                       0, -20; ...
                                      40, -40; ...
                                      20, 0; ...
                                      40, 40];
S.poslist{1}.electrode(3).thickness = [3, 1];
% add a 15 mm hole in the center:
S.poslist{1}.electrode(3).holes = sim_struct('ELECTRODE');
S.poslist{1}.electrode(3).holes.centre = 'Oz'; 
S.poslist{1}.electrode(3).holes.shape = 'ellipse';
S.poslist{1}.electrode(3).holes.dimensions = [15, 15];

%% Run simulation
run_simnibs(S);

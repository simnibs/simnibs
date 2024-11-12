%
% example script that runs a simulation for a Nx1 montage 
% 
% Copyright (c) 2021 SimNIBS developers. Licensed under the GPL v3.
%

%% General information
S = sim_struct('SESSION');
S.subpath = 'm2m_ernie'; % m2m-folder of the subject
S.pathfem = 'tdcs_Nx1'; %Folder for the simulation output
S.map_to_surf = true; % map to subject's middle gray matter surface (optional)
S.open_in_gmsh = true; % show results in gmsh (optional; standard: false)

%% Define tDCS simulation
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = 0.001; % Current flow through center channel (A)

%Center Electrode
S.poslist{1}.electrode(1).channelnr = 1; % Connect the electrode to the first channel
S.poslist{1}.electrode(1).centre = 'C3'; % Place it over C3
S.poslist{1}.electrode(1).shape = 'ellipse'; %Rectangular electrode
S.poslist{1}.electrode(1).dimensions = [30, 30]; % Dimension in mm
S.poslist{1}.electrode(1).thickness = [2, 1]; % 2 mm rubber electrodes on top of 1 mm gel layer

% when standard parameters are OK, using the following line 
% is enough to set up the surround electrodes:
% S.poslist{1} = expand_to_center_surround(S.poslist{1},S.subpath);


% optional parameters:
radius_surround = 60; % distance (centre-to-centre) between the centre and 
                     % surround electrodes (optional; standard: 50 mm)
                     % either a single number or an array with N entries
                     % (N: number of electrodes)             
pos_dir_1stsurround = 'C4'; % a position indicating the direction in which the 
                           % first surround electrode should be placed 
                          % (optional; standard: [])                        
N = 4; % number of surround electrodes (optional; standard: 4)
multichannel = false; % when set to true: Simulation of multichannel stimulator 
                     % with each suround channel receiving 1/N-th of the
                     % center channel (optional; standard: false, i.e. all 
                     % surround electrodes connected to the same channel)
phis_surround = []; % Angles in degree at which the electrodes will be place 
                     % relative to the direction defined by pos_dir_1stsurround.
                     % (optional; standard: [], resulting in equal distances
                     % between surround electrodes)
                     
S.poslist{1} = expand_to_center_surround(S.poslist{1}, S.subpath, ...
                                         'radius_surround', radius_surround, ...
                                         'N', N, ...
                                         'pos_dir_1stsurround', pos_dir_1stsurround, ...
                                         'multichannel', multichannel, ...
                                         'phis_surround', phis_surround);
                                    
%% Run Simulation
run_simnibs(S);
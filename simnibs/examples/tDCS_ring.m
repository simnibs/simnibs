%
% example script that runs a simulation
% for a center-surround ring montage
% 
% G. Saturnino, A. Thielscher, 2018
%

%% General information

S = sim_struct('SESSION'); % Define a stimulation sessions
S.fnamehead = 'ernie.msh'; % Choose the head mesh
S.pathfem = 'tdcs_ring'; % Folder for the simulation output

%% Define the tDCS simulation
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = [1e-3, -1e-3]; % Current going through each channel, in Ampere

%First Electrode (circular, central)
S.poslist{1}.electrode(1).channelnr = 1; % Connect it to the fisrt channel (1 mA)
S.poslist{1}.electrode(1).centre = 'C3'; % Place it over C3
S.poslist{1}.electrode(1).shape = 'ellipse'; % Define it to be elliptical/circular
S.poslist{1}.electrode(1).dimensions = [34, 34]; % Electrode diameter (in mm)
S.poslist{1}.electrode(1).thickness = 4; % Electrode 

%Second Electrode (ring, surrounding)
S.poslist{1}.electrode(2).channelnr = 2; %Connect it to the second channel (-1mA)
S.poslist{1}.electrode(2).centre = 'C3'; % Place it over C3
S.poslist{1}.electrode(2).shape = 'ellipse'; % Elliptical shape
S.poslist{1}.electrode(2).dimensions = [100, 100]; %Diameter of 100mm
S.poslist{1}.electrode(2).thickness = 4;

% Hole
S.poslist{1}.electrode(2).holes = sim_struct('ELECTRODE'); % Define the hole
S.poslist{1}.electrode(2).holes.centre = 'C3'; % Hole is also centered in C3
S.poslist{1}.electrode(2).holes.shape = 'ellipse'; % Shape of the hole
S.poslist{1}.electrode(2).holes.dimensions = [75, 75]; % Diameter of 75mm

%% Run Simulation
run_simnibs(S);

%% Visualize Simulations
m = mesh_load_gmsh4(fullfile(S.pathfem, 'ernie_TDCS_1_scalar.msh'));
mesh_show_surface(m);

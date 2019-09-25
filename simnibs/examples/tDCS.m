%
% example script that runs a simple simnibs tDCS simulation
% 
% G. Saturnino, A. Thielscher, 2018
%

%% General information

S = sim_struct('SESSION');
S.fnamehead = 'ernie.msh'; % head mesh
S.pathfem = 'tdcs'; %Folder for the simulation output

%% Define tDCS simulation
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = [0.001, -0.001]; % Current flow though each channel (mA)

%First Electrode
S.poslist{1}.electrode(1).channelnr = 1; % Connect the electrode to the first channel
S.poslist{1}.electrode(1).centre = 'C3'; % Place it over C3
S.poslist{1}.electrode(1).shape = 'rect'; %Rectangular electrode
S.poslist{1}.electrode(1).dimensions = [50, 50]; % Dimension in mm
S.poslist{1}.electrode(1).thickness = 4; % 4 mm thickness

%Second Electrode
S.poslist{1}.electrode(2).channelnr = 2;
S.poslist{1}.electrode(2).centre = 'AF4';
S.poslist{1}.electrode(2).shape = 'rect';
S.poslist{1}.electrode(2).dimensions = [50, 70];
S.poslist{1}.electrode(2).thickness = 4;

%% Run Simulation
run_simnibs(S);

%% Visualize Simulations
m = mesh_load_gmsh4(fullfile(S.pathfem, 'ernie_TDCS_1_scalar.msh'));
mesh_show_surface(m);

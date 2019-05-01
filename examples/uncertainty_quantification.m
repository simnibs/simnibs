% example script that runs an Uncertainty Quantification analysis
% 
% G. Saturnino, 2019
%

%% General information

S = sim_struct('SESSION');
S.fnamehead = 'ernie.msh'; % head mesh
S.pathfem = 'tdcs_uq'; %Folder for the simulation output

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

% Set-up the uncertain conductivities
% White Matter
S.poslist{1}.cond(1).distribution_type = 'beta';
S.poslist{1}.cond(1).distribution_parameters = [3, 3, .1, .4];
% Gray Matter
S.poslist{1}.cond(2).distribution_type = 'beta';
S.poslist{1}.cond(2).distribution_parameters = [3, 3, .1, .6];
% Skull
S.poslist{1}.cond(4).distribution_type = 'beta';
S.poslist{1}.cond(4).distribution_parameters = [3, 3, 0.003, 0.012];

%% Run Simulation
run_simnibs(S)

%
% script for simulating the TACS challenge montage
% see @TACSchallenge on Twitter
% 
% A. Thielscher, 2020
%

%% General information
S = sim_struct('SESSION');
S.subpath = 'm2m_ernie'; % subject folder
S.pathfem   = 'TACSchallenge'; % folder for the simulation output
S.map_to_surf = true;
S.open_in_gmsh = true; % open results once they are ready
                       % (set to false if you are annoyed by the popup windows)

                       
%% Define Condition A

% Set current flow though each channel:
% 2mA peak-to-peak --> 1 mA baseline-to-peak
% The third entry is a "pseudochannel" to which electrodes 
% are connected that are not used in the simulated condition
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = [0.001, -0.001, 0.0]; 

% Define occipital electrodes
S.poslist{1}.electrode(1).channelnr = 1; % Connect the electrode to the first channel
S.poslist{1}.electrode(1).centre = 'O1'; % Place it over O1
S.poslist{1}.electrode(1).shape = 'ellipse'; % Define it to be elliptical/circular
S.poslist{1}.electrode(1).dimensions = [25, 25]; % Electrode diameter (in mm)
S.poslist{1}.electrode(1).thickness = [1.5, 1.5]; % 1.5 mm gel and 1.5 mm rubber 

S.poslist{1}.electrode(2) = S.poslist{1}.electrode(1);
S.poslist{1}.electrode(2).centre = 'Oz';

S.poslist{1}.electrode(3) = S.poslist{1}.electrode(1);
S.poslist{1}.electrode(3).centre = 'O2';

% Define return electrodes
S.poslist{1}.electrode(4).channelnr = 2; % Connect to the second channel
S.poslist{1}.electrode(4).centre = 'CP5';
S.poslist{1}.electrode(4).shape = 'ellipse';
S.poslist{1}.electrode(4).dimensions = [50, 50];
S.poslist{1}.electrode(4).thickness = [1.5, 1.5];

S.poslist{1}.electrode(5) = S.poslist{1}.electrode(4);
S.poslist{1}.electrode(5).centre = 'CPz';

S.poslist{1}.electrode(6) = S.poslist{1}.electrode(4);
S.poslist{1}.electrode(6).centre = 'CP6';

% Define retinal control electrodes
S.poslist{1}.electrode(7).channelnr = 3; % Connect to the "pseudochannel"
S.poslist{1}.electrode(7).centre = [35.4, 106.0, -40.6]; % coordinates determined in simnibs_gui
S.poslist{1}.electrode(7).shape = 'ellipse';
S.poslist{1}.electrode(7).dimensions = [25, 25];
S.poslist{1}.electrode(7).thickness = [1.5, 1.5];

S.poslist{1}.electrode(8) = S.poslist{1}.electrode(7);
S.poslist{1}.electrode(8).centre = [-33.2, 108, -40.6]; % coordinates determined in simnibs_gui


%% Define Condition B (retinal control)
S.poslist{2} = S.poslist{1};

S.poslist{2}.electrode(1).channelnr = 3;
S.poslist{2}.electrode(2).channelnr = 3;
S.poslist{2}.electrode(3).channelnr = 3;
S.poslist{2}.electrode(7).channelnr = 1;
S.poslist{2}.electrode(8).channelnr = 1;


%% Define Condition C (cutaneous control)
S.poslist{3} = S.poslist{1};

S.poslist{3}.electrode(1).channelnr = 2;
S.poslist{3}.electrode(2).channelnr = 1;
S.poslist{3}.electrode(3).channelnr = 2;
S.poslist{3}.electrode(4).channelnr = 3;
S.poslist{3}.electrode(5).channelnr = 3;
S.poslist{3}.electrode(6).channelnr = 3;
S.poslist{3}.electrode(7).channelnr = 3;
S.poslist{3}.electrode(8).channelnr = 3;


%% Run Simulation
run_simnibs(S);

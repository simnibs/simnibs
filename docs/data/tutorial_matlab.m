%% General Settings
% Initialize a session
s = sim_struct('SESSION');
% Name of head mesh
s.fnamehead = 'ernie.msh';
% Output folder
s.pathfem = 'tutorial/';

%% TMS simulations
% Initialize a list of TMS simulations
s.poslist{1} = sim_struct('TMSLIST');
% Select coil
s.poslist{1}.fnamecoil = fullfile('legacy_and_other','Magstim_70mm_Fig8.ccd');

%% First coil position
% Select coil centre
s.poslist{1}.pos(1).centre = 'C1';
% Select coil direction
s.poslist{1}.pos(1).pos_ydir = 'CP1';

%% Second coil position
% Centred at C1
s.poslist{1}.pos(2).centre = 'C1';
% Pointing towards Cz
s.poslist{1}.pos(2).pos_ydir = 'Cz';

%% tDCS Simulation
% Initialize a tDCS simulation
s.poslist{2} = sim_struct('TDCSLIST');
% Set currents
s.poslist{2}.currents = [-1e-3 1e-3];

%% Define cathode
% Connect electrode to first channel (-1e-3 mA, cathode)
s.poslist{2}.electrode(1).channelnr = 1;
% Electrode dimension
s.poslist{2}.electrode(1).dimensions = [50 70];
% Rectangular shape
s.poslist{2}.electrode(1).shape = 'rect';
% 5mm thickness
s.poslist{2}.electrode(1).thickness = 5;
% Electrode Position
s.poslist{2}.electrode(1).centre = 'C3';
% Electrode direction
s.poslist{2}.electrode(1).pos_ydir = 'Cz';

%% Define anode
% Assign the electrode to the second channel
s.poslist{2}.electrode(2).channelnr = 2;
% Electrode diameter
s.poslist{2}.electrode(2).dimensions = [30 30];
% Electrode shape
s.poslist{2}.electrode(2).shape = 'ellipse';
% Electrode thickness
s.poslist{2}.electrode(2).thickness = 5;
% Electrode position
s.poslist{2}.electrode(2).centre = 'C4';

%% Run Simulations
run_simnibs(s)

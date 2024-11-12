%
% example script that runs a simple simnibs TMS simulation
% 
% Copyright (c) 2018 SimNIBS developers. Licensed under the GPL v3.
%

%% General information

S = sim_struct('SESSION');
S.subpath = 'm2m_ernie'; % subject folder
S.pathfem = 'tms_simu'; % folder for the simulation output

%% Define TMS simulation
S.poslist{1} = sim_struct('TMSLIST');
S.poslist{1}.fnamecoil = fullfile('legacy_and_other','Magstim_70mm_Fig8.ccd'); % Choose a coil model

%Define Position
S.poslist{1}.pos(1).centre = 'C3'; % Place it over C3
S.poslist{1}.pos(1).pos_ydir = 'CP3'; % Polongation of coil handle (see documentation)
S.poslist{1}.pos(1).distance = 4; % 4 mm distance from coil surface to head surface


%% Run Simulation
run_simnibs(S);

%% Visualize Simulations
m = mesh_load_gmsh4(fullfile(S.pathfem, 'ernie_TMS_1-0001_Magstim_70mm_Fig8_scalar.msh'));
mesh_show_surface(m);

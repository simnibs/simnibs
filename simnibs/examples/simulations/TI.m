%
% example script that runs two simnibs tDCS simulations
% and calculates maximal amplitude of the TI envelope from the E-fields 
% 
% A. Thielscher, 2020
%

%% General information

S = sim_struct('SESSION');
S.fnamehead = 'ernie.msh'; % head mesh
S.pathfem = 'TI'; %Folder for the simulation output

%% Define tDCS simulation: First electrode pair
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = [0.001, -0.001]; % Current flow though each channel (mA)

%First Electrode
S.poslist{1}.electrode(1).channelnr = 1; % Connect the electrode to the first channel
S.poslist{1}.electrode(1).centre = 'F5'; % Place it over F5
S.poslist{1}.electrode(1).shape = 'ellipse'; %round electrodes
S.poslist{1}.electrode(1).dimensions = [40, 40]; % Dimension in mm
S.poslist{1}.electrode(1).thickness = 2; % 4 mm thickness

%Second Electrode
S.poslist{1}.electrode(2).channelnr = 2;
S.poslist{1}.electrode(2).centre = 'P5';
S.poslist{1}.electrode(2).shape = 'ellipse';
S.poslist{1}.electrode(2).dimensions = [40, 40];
S.poslist{1}.electrode(2).thickness = 2;

%% Define second electrode pair: identical, except for electrode positions
S.poslist{2} = S.poslist{1};

%First Electrode
S.poslist{2}.electrode(1).centre = 'F6'; % update Place it over F6
S.poslist{2}.electrode(2).centre = 'P6';


%% Run Simulation
run_simnibs(S);

%% analyse simulations
m1 = mesh_load_gmsh4(fullfile(S.pathfem, 'ernie_TDCS_1_scalar.msh'));
m2 = mesh_load_gmsh4(fullfile(S.pathfem, 'ernie_TDCS_2_scalar.msh'));

% remove all tetrahedra and triangles belonging to the electrodes so that
% the two meshes have same number of elements
m1=mesh_extract_regions(m1, 'region_idx', [1:99 1001:1099]);
m2=mesh_extract_regions(m2, 'region_idx', [1:99 1001:1099]);

% calculate the maximal amplitude of the TI envelope for all elements
maxTI=get_maxTI( m1.element_data{get_field_idx(m1,'E','elem')}, ...
                 m2.element_data{get_field_idx(m2,'E','elem')} );

% make a new mesh that contains the field strengths and the 
% amplitude of the TI envelope for visualization
mout = m1;
mout.element_data={};
mout.element_data{1} = m1.element_data{get_field_idx(m1,'normE','elem')};
mout.element_data{1}.name = 'normE_1';
mout.element_data{2} = m2.element_data{get_field_idx(m2,'normE','elem')};
mout.element_data{2}.name = 'normE_2';
mout.element_data{3} = maxTI;

% save for visualization in gmsh
mesh_save_gmsh4(mout, fullfile(S.pathfem, 'ernie_TIenvelope.msh'));

% show field strengths and amplitude of TI envelope on GM surface
mesh_show_surface(mout,'field_idx','normE_1');
mesh_show_surface(mout,'field_idx','normE_2');
mesh_show_surface(mout,'field_idx','maxTIamplitude');


% Setting-up a tDCS Simulation using SimNIBS
%% Input and Output Names
S = sim_struct('SESSION');
% file name oh head mesh
S.fnamehead = 'ernie.msh';
% name of "m2m_sub_id"
S.subpath = 'm2m_ernie';
% Output Folder
S.pathfem = 'tdcs_simulation';

% Output Options:
%    Fields: e: Electric field norm
%            E: Electric field vector
%            j: Current density norm
%            J: Current density vector
S.fields = 'eEjJ';
% Transformations:
%  map to subject's middle gray matter surface
%  (only possible if model is ran with headreco --cat or mri2mesh)
S.map_to_surf = true;
%  map to FreeSurfer's FSAverage middle gray matter surface
%  (only possible if model is ran with headreco --cat or mri2mesh)
S.map_to_fsavg = true;
%  map to a nifti volume
S.map_to_vol = true;
%  map to the MNI template
S.map_to_MNI = true;


%% Set-up first tDCS simulation

S.poslist{1} = sim_struct('TDCSLIST');
% Current flow through each channel (in Ampere)
S.poslist{1}.currents = [0.001, -0.001];

%% Set-up first electrode
% Centre of electrode
S.poslist{1}.electrode(1).centre = [-50.0 -20.6 100.3];
% pos_ydir is a second point in the scalp.
% The electrode's y axis is defined starting from the centre and goes to
% pos_ydir 
S.poslist{1}.electrode(1).pos_ydir = [-44.2 -7.7 103.6];
% how the electrode is defined: 'plane': Define an electrode in 2D plus the
%                                        centre and the y axis of the
%                                        electrode
%                               'conf': Define the electrode vertices
%                                       directly in the head mesh
S.poslist{1}.electrode(1).definition = 'plane';
% Shape of the electrode: 'rect': rectangle
%                         'ellipse': elliptical (a circular electrode is an
%                                    elliptical electrode)
%                         'custom': custom shape. In this case, the vertices
%                                   variable must be filled out describing
%                                   each vertice in 2D
S.poslist{1}.electrode(1).shape = 'rect';
% Dimensions of the electrode, in mm
S.poslist{1}.electrode(1).dimensions = [40, 40];
% Thickness of the electrode layers.
%   1 number: electrode with 1 (gel) layer.
%   2 numbers:electrode with 2 (gel, conductive rubber) layers
%   3 numbers:electrode with 3 (sponge, conductive rubber, sponge) layers
%             in this case, dimensions_sponge must also be filled out
S.poslist{1}.electrode(1).thickness = [3.5, 1];
% Connect the electrode to the first channel (current = 0.001 A)
S.poslist{1}.electrode(1).channelnr = 1;

%% Set-up second electrode
S.poslist{1}.electrode(2).centre = [32.1 81.1 56.9];
S.poslist{1}.electrode(2).pos_ydir = [42.2 73.5 56.4];
S.poslist{1}.electrode(2).definition = 'plane';
S.poslist{1}.electrode(2).shape = 'rect';
S.poslist{1}.electrode(2).dimensions = [40, 40];
S.poslist{1}.electrode(2).thickness = [3.5, 1, 3.5];
% Now, as we defined a sponge electrode using a thickness with 3 elements,
% we need to define dimensions_sponge
S.poslist{1}.electrode(2).dimensions_sponge = [50, 70];
% We also want to define a plug for the electrode
% the plug is defined using the same structure of the electrode,
% buth with several fields, such as 'channelnr', 'holes', 'plug' left blank
S.poslist{1}.electrode(2).plug = sim_struct('ELECTRODE');
% The plug centre can be defined using electrode coordinates
S.poslist{1}.electrode(2).plug.centre = [-10, 0];
S.poslist{1}.electrode(2).plug.shape = 'rect';
S.poslist{1}.electrode(2).plug.dimensions = [10 10];
S.poslist{1}.electrode(2).channelnr = 2;

%% Set-up second tDCS simulation
% Now, we will set-up a simulation using a ring set-up
S.poslist{2} = sim_struct('TDCSLIST');
S.poslist{2}.currents = [0.001, -0.001];

% Central circular electrode, with a diameter of 34 mm
S.poslist{2}.electrode(1).centre = [-50.0 -20.6 100.3];
S.poslist{2}.electrode(1).definition = 'plane';
S.poslist{2}.electrode(1).shape = 'ellipse';
S.poslist{2}.electrode(1).dimensions = [34, 34];
S.poslist{2}.electrode(1).channelnr = 1;
S.poslist{2}.electrode(1).thickness = [2, 2];

% Outer ring electrode, with an outer diameter of 100mm and an inner
% diameter of 75mm
S.poslist{2}.electrode(2).centre = [-50.0 -20.6 100.3];
S.poslist{2}.electrode(2).definition = 'plane';
S.poslist{2}.electrode(2).shape = 'ellipse';
S.poslist{2}.electrode(2).dimensions = [100, 100];
S.poslist{2}.electrode(2).channelnr = 2;
S.poslist{2}.electrode(2).thickness = [2, 2];
S.poslist{2}.electrode(2).holes = sim_struct('ELECTRODE');
S.poslist{2}.electrode(2).holes.centre = [0 0];
S.poslist{2}.electrode(2).definition = 'plane';
S.poslist{2}.electrode(2).holes.shape = 'ellipse';
S.poslist{2}.electrode(2).holes.dimensions = [75, 75];

%% Run Simulation
run_simnibs(S);

%% Visualize Simulations
m = mesh_load_gmsh4(fullfile('tdcs_simulation', 'ernie_TDCS_1_scalar.msh'));
mesh_show_surface(m);
m = mesh_load_gmsh4(fullfile('tdcs_simulation', 'ernie_TDCS_2_scalar.msh'));
mesh_show_surface(m);
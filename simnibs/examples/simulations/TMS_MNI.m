% General Infoarmation
S = sim_struct('SESSION');
S.subpath = 'm2m_ernie'; % subject folder
S.pathfem = 'tms_hand';  % Directory for the simulation


% Define the TMS simulation
S.poslist{1} = sim_struct('TMSLIST');
S.poslist{1}.fnamecoil = fullfile('legacy','Magstim_70mm_Fig8.ccd');  % Choose a coil from the resources/coil_models folder

% Define the coil position
% Place coil over the hand knob
% Here, the hand knob is defined in MNI coordinates (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2034289/)
% And transformed to subject coordinates
% We can use positions in the cortex. SimNIBS will automatically project them to the skin surface
% and add the specified distance
S.poslist{1}.pos(1).centre = mni2subject_coords([-37, -21, 58], 'm2m_ernie');
% Point the coil handle posteriorly, we just add 10 mm to the original M1 "y" coordinate
S.poslist{1}.pos(1).pos_ydir = mni2subject_coords([-37, -21-10, 58], 'm2m_ernie');
S.poslist{1}.pos(1).distance = 4; % 4 mm distance from coil surface to head surface


% Run Simulation
run_simnibs(S);

%% Visualize Simulations
m = mesh_load_gmsh4(fullfile(S.pathfem, 'ernie_TMS_1-0001_Magstim_70mm_Fig8_nii_scalar.msh'));
mesh_show_surface(m);

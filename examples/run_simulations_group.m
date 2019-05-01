% Set the subjects
subjects = {'sub01', 'sub09', 'sub10', 'sub12', 'sub15'};

% Start a SESSION
S = sim_struct('SESSION');
S.map_to_fsavg = true;

% Set a TDCSLIST with the simulation set-up
S.poslist{1} = sim_struct('TDCSLIST');
S.poslist{1}.currents = [0.001, -0.001];

S.poslist{1}.electrode(1).channelnr = 1;
S.poslist{1}.electrode(1).centre = 'C3';
S.poslist{1}.electrode(1).pos_ydir = 'C1';
S.poslist{1}.electrode(1).shape = 'rect';
S.poslist{1}.electrode(1).dimensions = [50, 50];
S.poslist{1}.electrode(1).thickness = 4;

S.poslist{1}.electrode(2).channelnr = 2;
S.poslist{1}.electrode(2).centre = 'AF4';
S.poslist{1}.electrode(2).pos_ydir = 'F6';
S.poslist{1}.electrode(2).shape = 'rect';
S.poslist{1}.electrode(2).dimensions = [50, 70];
S.poslist{1}.electrode(2).thickness = 4

% Run the simulation in each subject
for i = 1:length(subjects)
     sub = subjects{i};
     S.fnamehead = fullfile(sub, [sub '.msh']);  % head mesh
     S.pathfem = fullfile(sub, 'bipolar'); % Output directory
     run_simnibs(S)
end

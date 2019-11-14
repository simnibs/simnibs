% Example for running a SimNIBS tDCS leadfield 
%    Copyright (C) 2019 Guilherme B Saturnino

tdcs_lf = sim_struct('TDCSLEADFIELD');
% Head mesh
tdcs_lf.fnamehead = 'ernie.msh';
% Output directory
tdcs_lf.pathfem = 'tdcs_leadfield';
% Interpolate to the middle gray matter surface
tdcs_lf.map_to_surf = true; 
run_simnibs(tdcs_lf)

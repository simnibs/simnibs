% Example for running a SimNIBS tDCS leadfield 
%    Copyright (C) 2019 Guilherme B Saturnino

tdcs_lf = sim_struct('TDCSLEADFIELD');
% Head mesh
tdcs_lf.fnamehead = 'ernie.msh';
% Output directory
tdcs_lf.pathfem = 'leadfield';

run_simnibs(tdcs_lf)

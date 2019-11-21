% Example of a SimNIBS tDCS leadfield 
% Copyright (C) 2019 Guilherme B Saturnino

tdcs_lf = sim_struct('TDCSLEADFIELD');
% Head mesh
tdcs_lf.fnamehead = 'ernie.msh';
% Output directory
tdcs_lf.pathfem = 'leadfield';

% Uncoment to use the pardiso solver
%tdcs_lf.solver_options = 'pardiso';
% This solver is faster than the default. However, it requires much more
% memory (~12 GB)

run_simnibs(tdcs_lf)

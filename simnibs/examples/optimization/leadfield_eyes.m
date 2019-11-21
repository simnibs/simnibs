% Example of a tDCS leadfield with Eyes
%
% Copyright (C) 2019 Guilherme B Saturnino

tdcs_lf = sim_struct('TDCSLEADFIELD');
% head mesh
tdcs_lf.fnamehead = 'ernie.msh';
% output directory
tdcs_lf.pathfem = 'leadfield_eyes';

% Store field in eyes (1006) and gray matter (1002) surfaces
tdcs_lf.tissues = [1006, 1002];

% Do not map to middle gray matter surface (Overwrites TISSUE)
tdcs_lf.map_to_surf = false;

tdcs_lf.solver_options = 'pardiso';

% Uncoment to use the pardiso solver
%tdcs_lf.solver_options = 'pardiso';
% This solver is faster than the default. However, it requires much more
% memory (~12 GB)



run_simnibs(tdcs_lf)

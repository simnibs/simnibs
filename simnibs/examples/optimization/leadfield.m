% Example of a SimNIBS tDCS leadfield 
% Copyright (c) 2019 SimNIBS developers. Licensed under the GPL v3.

% place script in the main folder of the example dataset

tdcs_lf = sim_struct('TDCSLEADFIELD');
% subject folder
tdcs_lf.subpath = 'm2m_ernie';
% Output directory
tdcs_lf.pathfem = 'leadfield';

% Uncomment to use the pardiso solver
%tdcs_lf.solver_options = 'pardiso';
% This solver is much faster than the default. However, it requires much more
% memory (~12 GB)

run_simnibs(tdcs_lf)

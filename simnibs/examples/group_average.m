% This example wil go throgh simulations and calculate
% the average and the standard deviation of the normal component
% of the electric field in FsAverage space
% 
% It is a follow-up to the "run_simulations" example

subjects = {'sub01', 'sub09', 'sub10', 'sub12', 'sub15'};
results_folder = fullfile('bipolar', 'fsavg_overlays');

normals = {};
for i = 1:length(subjects)
    sub = subjects{i};
    % Load normal field data
    m = mesh_load_gmsh4(fullfile(sub, results_folder, [sub '_TDCS_1_scalar_fsavg.msh']));
    % Add to cell
    normals{i} = m.node_data{get_field_idx(m, 'E_normal', 'node')}.data;
end
% Calculate average and standard deviation of the normal at each node
normals = cell2mat(normals);
avg_normal = mean(normals, 2);
std_normal = std(normals, 0, 2);
% Place the fields in the mesh structure
m.node_data{1}.data = avg_normal;
m.node_data{1}.name = 'E_normal_avg';
m.node_data{2}.data = std_normal;
m.node_data{2}.name = 'E_normal_std';
% Plot the fields
mesh_show_surface(m, 'field_idx', 'E_normal_avg')
mesh_show_surface(m, 'field_idx', 'E_normal_std')

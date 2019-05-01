subjects = {'sub01', 'sub09', 'sub10', 'sub12', 'sub15'};
results_folder = 'bipolar/fsavg_overlays';

normals = {};
for i = 1:length(subjects)
    sub = subjects{i};
    % Load normal field data
    normal_surf = ['lh.' sub '_TDCS_1_scalar.fsavg.E.normal'];
    m = mesh_load_fsresults(fullfile(sub, results_folder, normal_surf));
    % Add to cell
    normals{i} = m.node_data{1}.data;
end
% Calculate average and standard deviation of the normal at each node
normals = cell2mat(normals);
avg_normal = mean(normals, 2);
std_normal = std(normals, 0, 2);
% Place the fields in the mesh structure
m.node_data{1}.data = avg_normal;
m.node_data{1}.name = 'E.normal.avg';
m.node_data{2}.data = std_normal;
m.node_data{2}.name = 'E.normal.std';
% Plot the fields
mesh_show_surface(m, 'field_idx', 'E.normal.avg')
mesh_show_surface(m, 'field_idx', 'E.normal.std')
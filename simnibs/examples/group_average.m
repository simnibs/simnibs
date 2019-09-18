% This example wil go throgh simulations and calculate
% the average and the standard deviation of the normal component
% of the electric field in FsAverage space
% 
% It is a follow-up to the "run_simulations" example

subjects = {'sub01', 'sub09', 'sub10', 'sub12', 'sub15'};
results_folder = fullfile('bipolar', 'fsavg_overlays');
fsavg_msh_name = '_TDCS_1_scalar_fsavg.msh';

normals = {};
for i = 1:length(subjects)
    sub = subjects{i};
    % Load normal field data
    m = mesh_load_gmsh4(fullfile(sub, results_folder, [sub fsavg_msh_name]));
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

% Calculate average in an ROI defined using the a2009s atlas
[labels, snames]=mesh_load_fssurf('fsaverage','label','a2009s');
region_name = 'lh.S_precentral-sup-part';
roi_idx=find(strcmpi(snames, region_name));
node_idx = labels.node_data{1}.data==roi_idx;
% visualize region
labels.node_data{2}.data = int8(node_idx);
labels.node_data{2}.name = region_name;
mesh_show_surface(labels, 'field_idx', region_name)
% weight fields using node areas
nodes_areas = mesh_get_node_areas(m);
avg_field_roi = ...
    sum(avg_normal(node_idx).*nodes_areas(node_idx))/sum(nodes_areas(node_idx));
sprintf('Average field in left superior precental sulcus: %f',avg_field_roi)


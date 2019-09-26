% This example wil go throgh simulations and calculate
% the average and the standard deviation of the normal component
% of the electric field in FsAverage space
% 
% It is a follow-up to the "run_simulations" example

%% Load simulation results
subjects = {'sub01', 'sub09', 'sub10', 'sub12', 'sub15'};
results_folder = fullfile('bipolar', 'fsavg_overlays');
fsavg_msh_name = '_TDCS_1_scalar_fsavg.msh';
field_name = 'E_normal';

fields = {};
for i = 1:length(subjects)
    sub = subjects{i};
    % load mesh with results transformed to fsaverage space
    m = mesh_load_gmsh4(fullfile(pwd, sub, results_folder, [sub fsavg_msh_name]));
    % Save the field of each subject
    fields{i} = m.node_data{get_field_idx(m, field_name, 'node')}.data;
end
%% Calculate and plot averages
% Calculate
fields = cell2mat(fields);
avg_field = mean(fields, 2);
std_field = std(fields, 0, 2);
% Plot
m.node_data = {}; %cleanup fields
m.node_data{1}.data = avg_field; % add average field
m.node_data{1}.name = [field_name '_avg'];
m.node_data{2}.data = std_field; % add std field
m.node_data{2}.name = [field_name '_std'];

% show surfaces with fields
mesh_show_surface(m, 'field_idx', [field_name '_avg'])
mesh_show_surface(m, 'field_idx', [field_name '_std'])

%% Calculate average in an ROI defined using an atlas
% load atlas and define a region
[m, snames]=mesh_load_fssurf('fsaverage','label','a2009s');
region_name = 'lh.S_precentral-sup-part';
roi_idx=find(strcmpi(snames, region_name));
node_idx = m.node_data{end}.data==roi_idx;
% visualize region
m.node_data{end+1}.data = int8(node_idx);
m.node_data{end}.name = region_name;
mesh_show_surface(m, 'field_idx', region_name)

% calculate a weighted mean using the node areas
nodes_areas = mesh_get_node_areas(m);
avg_field_roi = sum(avg_field(node_idx).*nodes_areas(node_idx))/sum(nodes_areas(node_idx));
fprintf('Average %s in %s: %f\n', field_name, region_name, avg_field_roi)

% ROI analysis of the electric field from a simulation using an atlas.
% We will calculate the mean electric field in a gray matter 
% ROI defined using an atlas
%
% Copyright (c) 2019 SimNIBS developers. Licensed under the GPL v3.

% Read the simulation results interpolated to the middle gray matter surface
surf = mesh_load_gmsh4(...
    fullfile('tdcs_simu', 'subject_overlays', ...
             'ernie_TDCS_1_scalar_central.msh')...
);

% Load the atlas and define the brain region of interest
[labels, snames] = subject_atlas(surf, fullfile(pwd, 'm2m_ernie'), 'HCP_MMP1');
region_name = 'lh.4';
roi_idx=find(strcmpi(snames, region_name));
node_idx = labels.node_data{end}.data==roi_idx;

% Plot the ROI
surf.node_data{end+1}.data = int8(node_idx);
surf.node_data{end}.name = region_name;
mesh_show_surface(surf, 'field_idx', region_name)

% calculate the node areas, we will use those later for averaging
nodes_areas = mesh_get_node_areas(surf);

% Get the field of interest
field_name = 'E_magn';
field_idx = get_field_idx(surf, field_name, 'node');
field = surf.node_data{field_idx}.data;

% Calculate the mean
avg_field_roi = sum(field(node_idx).*nodes_areas(node_idx))/sum(nodes_areas(node_idx));
fprintf('Average %s in %s: %f\n', field_name, region_name, avg_field_roi)


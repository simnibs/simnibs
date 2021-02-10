% Simple ROI analysis of the electric field from a simulation.
% We will calculate the mean electric field in a gray matter ROI
% The ROI is defined using an MNI coordinates

%% Load Simulation Result

% Read the simulation result
head_mesh = mesh_load_gmsh4(...
    fullfile('tdcs', 'ernie_TDCS_1_scalar.msh') ...
);

% Crop the mesh so we only have gray matter volume elements (tag 2 in the mesh)
gray_matter = mesh_extract_regions(head_mesh, 'region_idx', 2);


%% Define the ROI

% Define M1 from MNI coordinates (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2034289/)
% the first argument is the MNI coordinates
% the second argument is the subject "m2m" folder
ernie_coords = mni2subject_coords([-37, -21, 58], fullfile(pwd, 'm2m_ernie'));
% we will use a sphere of radius 10 mm
r = 10;

% Electric fields are defined in the center of the elements
% get element centers
elm_centers = mesh_get_tetrahedron_centers(gray_matter);
% determine the elements in the ROI
roi = sqrt(sum(bsxfun(@minus, elm_centers, ernie_coords).^2, 2)) < r;
% get element volumes, we will use those for averaging
elm_vols = mesh_get_tetrahedron_sizes(gray_matter);

%% Get field and calculate the mean
% Get the field of interest
field_name = 'magnE';
field_idx = get_field_idx(gray_matter, field_name, 'elements');
field = gray_matter.element_data{field_idx}.tetdata;

% Calculate the mean
avg_field_roi = sum(field(roi) .* elm_vols(roi))/sum(elm_vols(roi));
fprintf('mean %s in ROI: %f\n', field_name, avg_field_roi)

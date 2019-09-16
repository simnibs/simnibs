%Simple ROI analysis of the electric field from a simulation.
%
%We will calculate the mean electric field in a gray matter ROI defined around M1

%% Input

% Read the simulation result
head_mesh = mesh_load_gmsh4('tdcs/ernie_TDCS_1_scalar.msh');

% Crop the mesh so we only have gray matter volume elements (tag 2 in the mesh)
gray_matter = mesh_extract_regions(head_mesh, 'region_idx', 2);


%% Define the ROI
% We will now define our ROI

% Define M1 from MNI coordinates (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2034289/)
% the first argument is the MNI coordinates
% the second argument is the segmentation output folder
ernie_M1 = mni2subject_coords([-37, -21, 58], 'm2m_ernie/');
% we will use a sphere of radius 10 mm
r = 10;

% Electric fields are defined in the center of the elements
% Therefore, we will select all elements which have their centers inside the ROI
elm_centers = mesh_get_tetrahedron_centers(gray_matter);
% calculate distances to M1
roi = sqrt(sum(bsxfun(@minus, elm_centers, ernie_M1).^2, 2)) < r;

% finally, calculate the mean of the electric field norm
normE = mesh_get_field(gray_matter, 'normE');
mean_normE = mean(normE.tetdata(roi));
fprintf('mean normE in M1 ROI: %f\n', mean_normE)


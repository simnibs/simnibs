%
% Examples to define ROIs 
%

% Surface ROI
% ======================================================================
% Define a midlayer ROI from an MNI center coordinate and a sphere radius
% ====================================================
roi = opt_struct('RegionOfInterest');
roi.subpath = 'm2m_ernie';
roi.method = 'surface';
roi.surface_type = 'central';

% center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_space = {'mni', 'mni'};
roi.roi_sphere_center = {[-38.6, -18.7, 64.8], [-38.6, -18.7, 64.8]};
roi.roi_sphere_operator = {'intersection', 'difference'};
% radius of spherical ROI (in mm)
roi.roi_sphere_radius = [20, 10];

% for visual control:
% define a mesh filename for ROI visualization, and call run_simnibs
roi.fname_visu = 'roi_surf_vis';
run_simnibs(roi)
clear roi


% Surface ROI
% ======================================================================
% Define a midlayer ROI from a mask in fsaverage space
% ====================================================   
roi = opt_struct('RegionOfInterest');
roi.subpath = 'm2m_ernie';
roi.method = 'surface';
roi.surface_type = 'central';

roi.mask_path = fullfile(SIMNIBSDIR, 'examples','utilities','P1_LH_M1_control');
roi.mask_space = 'fs_avg_lh';

% for visual control:
% define a mesh filename for ROI visualization, and call run_simnibs
roi.fname_visu = 'roi_P1_LH_M1_control';
run_simnibs(roi)
clear roi


% Volume ROI
% ======================================================================
% Define a volumetric ROI from a volume_from_surface MNI center coordinate and a sphere radius
% =================================
roi = opt_struct('RegionOfInterest');
roi.subpath = 'm2m_ernie';
roi.method = 'volume';

% domains to include in mask
roi.tissues = 2; % select gray matter

% center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_space = {'mni', 'mni'};
roi.roi_sphere_center = {[-38.6, -18.7, 64.8], [-38.6, -18.7, 64.8]};
roi.roi_sphere_operator = {'intersection', 'difference'};
% radius of spherical ROI (in mm)
roi.roi_sphere_radius = [40, 30];

% for visual control:
% define a mesh filename for ROI visualization, and call run_simnibs
roi.fname_visu = 'roi_vol_vis';
run_simnibs(roi)
clear roi


% Volume ROI from Surface ROI
% ======================================================================
% Define a volumetric ROI from a surface ROI
% =================================
roi = opt_struct('RegionOfInterest');
roi.subpath = 'm2m_ernie';
roi.method = 'volume_from_surface';
roi.surface_type = 'central';

% domains to include
roi.tissues = 2; % select gray matter

% center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_space = 'mni';
roi.roi_sphere_center = [-38.6, -18.7, 64.8];
% radius of spherical ROI (in mm)
roi.roi_sphere_radius = 20;
% radius around surface in which volume elements will be included
roi.surface_inclusion_radius = 5;

% for visual control:
% define a mesh filename for ROI visualization, and call run_simnibs
roi.fname_visu = 'roi_vol_from_surf_vis';
run_simnibs(roi)
clear roi


% Custom ROI
% ======================================================================
% Define a custom ROI from coordinates
% =============================
roi = opt_struct('RegionOfInterest');
roi.subpath = 'm2m_ernie';
roi.method = 'custom';

% points where the e-field is evaluated (interpolated), e.g. centres of custom surface or volume
roi.nodes = {[-13.7, 13.1, 71.1], [-16.1, 12.4, 71.1], [-16.1, 12.4, 70.6]};

% for visual control:
% define a mesh filename for ROI visualization, and call run_simnibs
roi.fname_visu = 'roi_point_vis';
run_simnibs(roi)
clear roi

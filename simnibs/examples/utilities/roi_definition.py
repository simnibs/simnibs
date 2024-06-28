"""
Examples to define ROIs 
"""
import os
import simnibs
from simnibs import ElementTags


# Surface ROI
# ======================================================================
# Define a midlayer ROI from an MNI center coordinate and a sphere radius
# ====================================================
roi = simnibs.RegionOfInterest()
roi.subpath = "m2m_ernie"
roi.method = "surface"
roi.surface_type = "central"

# center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_space = ["mni", "mni"]
roi.roi_sphere_center = [[-38.6, -18.7, 64.8], [-38.6, -18.7, 64.8]]
roi.roi_sphere_operator = ["intersection", "difference"]
# radius of spherical ROI (in mm)
roi.roi_sphere_radius = [20, 10]

roi.write_visualization('', 'roi_surf_vis')
del roi


# Surface ROI
# ======================================================================
# Define a midlayer ROI from a mask in fsaverage space
# ====================================================   
roi = simnibs.RegionOfInterest()
roi.subpath = "m2m_ernie"
roi.method = "surface"
roi.surface_type = "central"
roi.mask_path = os.path.join(simnibs.SIMNIBSDIR, 'examples','utilities','P1_LH_M1_control')
roi.mask_space = 'fs_avg_lh'

roi.write_visualization('', 'roi_P1_LH_M1_control')
del roi


# Volume ROI
# ======================================================================
# Define a volumetric ROI from a volume_from_surface MNI center coordinate and a sphere radius
# =================================
roi = simnibs.RegionOfInterest()
roi.subpath = "m2m_ernie"
roi.method = "volume"

# domains to include in mask
roi.tissues = [ElementTags.GM]

# center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_space = ["mni", "mni"]
roi.roi_sphere_center = [[-38.6, -18.7, 64.8], [-38.6, -18.7, 64.8]]
roi.roi_sphere_operator = ["intersection", "difference"]
# radius of spherical ROI (in mm)
roi.roi_sphere_radius = [40, 30]

roi.write_visualization('', 'roi_vol_vis')
del roi


# Volume ROI from Surface ROI
# ======================================================================
# Define a volumetric ROI from a surface ROI
# =================================
roi = simnibs.RegionOfInterest()
roi.subpath = "m2m_ernie"
roi.method = "volume_from_surface"
roi.surface_type = "central"

# domains to include
roi.tissues = [ElementTags.GM]

# center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_space = "mni"
roi.roi_sphere_center = [-38.6, -18.7, 64.8]
# radius of spherical ROI (in mm)
roi.roi_sphere_radius = 20
# radius around surface in which volume elements will be included
roi.surface_inclusion_radius = 5

roi.write_visualization('', 'roi_vol_from_surf_vis')
del roi


# Volume ROI from Surface ROI mask
# ======================================================================
# Define a volumetric ROI from a mask in fsaverage space
# =================================
roi = simnibs.RegionOfInterest()
roi.subpath = "m2m_ernie"
roi.method = "volume_from_surface"
roi.surface_type = "central"
#path to mask in fsaverage space
roi.mask_path = os.path.join(simnibs.SIMNIBSDIR, 'examples','utilities','P1_LH_M1_control')
roi.mask_space = 'fs_avg_lh'
# radius around surface in which volume elements will be included
roi.surface_inclusion_radius = 2

roi.write_visualization('', 'roi_P1_LH_M1_control_volume')
del roi


# Custom ROI
# ======================================================================
# Define a custom ROI from coordinates
# =============================
roi = simnibs.RegionOfInterest()
roi.subpath = "m2m_ernie"
roi.method = "custom"

# points where the e-field is evaluated (interpolated), e.g. centres of custom surface or volume
roi.nodes = [[-13.7, 13.1, 71.1], [-16.1, 12.4, 71.1], [-16.1, 12.4, 70.6]]
roi.write_visualization('', 'roi_point_vis')
del roi

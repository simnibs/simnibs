"""
Example to define ROIs 
"""

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
roi.roi_sphere_center = [[-43.6, -35.4, 64.3],[-43.6, -35.4, 64.3]]
roi.roi_sphere_operator = ["intersection", "difference"]

# radius of spherical ROI (in mm)
roi.roi_sphere_radius = [20, 10]

roi.write_visualization('', 'roi_surf_vis_2')
del roi

# Volume ROI
# ======================================================================
# Define a volumetric ROI from a volume_from_surface MNI center coordinate and a sphere radius
# =================================
roi = simnibs.RegionOfInterest()
roi.subpath = "m2m_ernie"
roi.method = "volume"

# center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_space = ["mni", "mni"]
roi.roi_sphere_center = [[-43.6, -35.4, 64.3],[-43.6, -35.4, 64.3]]
roi.roi_sphere_operator = ["intersection", "difference"]
# radius of spherical ROI (in mm)
roi.roi_sphere_radius = [40, 30]

# domains to include
roi.tissues = [ElementTags.GM]
roi.write_visualization('', 'roi_vis_2')
del roi

# Volume ROI from Surface ROI
# ======================================================================
# Define a volumetric ROI from a surface ROI
# =================================
roi = simnibs.RegionOfInterest()
roi.subpath = "m2m_ernie"
roi.method = "volume_from_surface"
roi.surface_type = "central"

# center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_space = "mni"
roi.roi_sphere_center = [-43.6, -35.4, 64.3]

# radius of spherical ROI (in mm)
roi.roi_sphere_radius = 20

# domains to include
roi.tissues = [ElementTags.GM]

roi.surface_inclusion_radius = 5

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

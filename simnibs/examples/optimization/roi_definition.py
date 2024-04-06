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
roi.roi_sphere_center_space = "mni"
roi.roi_sphere_center = [-43.6, -35.4, 64.3]

# radius of spherical ROI (in mm)
roi.roi_sphere_radius = 20

# Volume ROI
# ======================================================================
# Define a volumetric ROI from a volume_from_surface MNI center coordinate and a sphere radius
# =================================
roi = simnibs.RegionOfInterest()
roi.subpath = "m2m_ernie"
roi.method = "volume"

# center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_space = "mni"
roi.roi_sphere_center = [-43.6, -35.4, 64.3]

# alternatively, the center can also be defined in subject space (in mm)
# roi.roi_sphere_center_space = "subject"

# radius of spherical ROI (in mm)
roi.roi_sphere_radius = 20

# domains to include
roi.tissues = [ElementTags.WM, ElementTags.GM]

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

roi.write_visualization('', 'roi_vis')

# Custom ROI
# ======================================================================
# Define a custom ROI from coordinates
# =============================
roi = simnibs.RegionOfInterest()
roi.subpath = "m2m_ernie"
roi.method = "custom"

# points where the e-field is evaluated (interpolated), e.g. centres of custom surface or volume
roi.nodes = [[-13.7, 13.1, 71.1], [-16.1, 12.4, 71.1], [-16.1, 12.4, 70.6]]

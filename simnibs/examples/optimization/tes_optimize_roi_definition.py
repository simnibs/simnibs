"""
Example to define ROIs for TESoptimize

Written by: Konstantin Weise (2023)
"""

import simnibs
from simnibs import ElementTags

# Initialize structure
opt = simnibs.opt_struct.TESoptimize()

# Option 1: Define a midlayer ROI
# ===============================
roi = opt.add_roi()
roi.type = "GMmidlayer"

# center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_mni = [-43.6, -35.4,  64.3]

# alternatively, the center can also be defined in subject space (in mm)
# roi.roi_sphere_center_subject = [-41, -13,  66]

# radius of spherical ROI (in mm)
roi.roi_sphere_radius = 20

# Option 2: Define a volumetric ROI
# =================================
roi = opt.add_roi()
roi.type = "volume"

# center of spherical ROI in MNI space (in mm)
roi.roi_sphere_center_mni = [-43.6, -35.4,  64.3]

# alternatively, the center can also be defined in subject space (in mm)
# roi.roi_sphere_center_subject = [-41, -13,  66]

# radius of spherical ROI (in mm)
roi.roi_sphere_radius = 20

# domains to include
roi.domains = [ElementTags.WM, ElementTags.GM]

# Option 3: Define a custom ROI
# =============================
roi = opt.add_roi()
roi.type = "custom"

# points where the e-field is evaluated (interpolated), e.g. centres of custom surface or volume
roi.center = [[-13.7, 13.1, 71.1],
              [-16.1, 12.4, 71.1],
              [-16.1, 12.4, 70.6]]


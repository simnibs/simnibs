'''
ROI analysis of the electric field from a simulation using an atlas.

We will calculate the mean electric field in a gray matter ROI defined using an atlas
'''
import os
import numpy as np
import simnibs

## Input ##

# Read the simulation result mapped to the gray matter surface
gm_surf = simnibs.read_msh(
    os.path.join('tdcs', 'subject_overlays',
                 'ernie_TDCS_1_scalar_central.msh')
)

# Load the atlas and define the brain region of interest
atlas = simnibs.subject_atlas('HCP_MMP1', 'm2m_ernie')
region_name = 'lh.4'
roi = atlas[region_name]

# plot the roi
gm_surf.add_node_field(roi, 'ROI')
gm_surf.view(visible_fields='ROI').show()

# calculate the node areas, we will use those later for averaging
node_areas = gm_surf.nodes_areas()

# finally, calculate the mean of the field strength
field_name = 'E_magn'
mean_magnE = np.average(gm_surf.field[field_name][roi], weights=node_areas[roi])
print('mean ', field_name, ' in ', region_name, ': ', mean_magnE)

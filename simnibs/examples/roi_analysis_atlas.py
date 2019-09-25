'''
ROI analysis of the electric field from a simulation using an atlas.

We will calculate the mean electric field in a gray matter ROI defined using an atlas
'''
import os
import numpy as np
import simnibs

## Input ##

# Read the simulation result mapped to the gray matter surface
gm_surf = simnibs.read_msh(os.path.join(
    'tdcs', 'subject_overlays', 'ernie_TDCS_1_scalar_central.msh'
))

## Define the ROI ##
# We will now define our ROI

# Load the HCP_MMP atlas transformed to subject space
atlas = simnibs.subject_atlas('HCP_MMP1', 'm2m_ernie/')
roi = atlas['lh.4']
# SOME PROBLEMS WITH THE SIZE OF THE ROI VECTOR??

# calculate the node areas
node_areas = gm_surf.nodes_areas()
breakpoint()
# finally, calculate the mean of the electric field norm
max_normE = np.average(gm_surf.field['E_normal'][roi], weights=node_areas[roi])
print('mean normE in M1 ROI: ', mean_normE)
# plot the roi
gm_surf.add_node_field(roi, 'ROI')
gm_surf.view(visible_fields='ROI').show()
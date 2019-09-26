'''
Simple ROI analysis of the electric field from a simulation.

We will calculate the mean electric field in a gray matter ROI defined around M1
'''
import os
import numpy as np
import simnibs

## Load simulation result

# Read the simulation result
head_mesh = simnibs.read_msh(
    os.path.join('tdcs_example_run', 'ernie_TDCS_1_scalar.msh')
)

# Crop the mesh so we only have gray matter volume elements (tag 2 in the mesh)
gray_matter = head_mesh.crop_mesh(2)


## Define the ROI

# Define M1 from MNI coordinates (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2034289/)
# the first argument is the MNI coordinates
# the second argument is the subject "m2m" folder
ernie_coords = simnibs.mni2subject_coords([-37, -21, 58], 'm2m_ernie/')
# we will use a sphere of radius 10 mm
r = 10.

# Electric fields are defined in the center of the elements
# get element centers
elm_centers = gray_matter.elements_baricenters()[:]
# determine the elements in the ROI
roi = np.linalg.norm(elm_centers - ernie_coords, axis=1) < r
# get the element volumes, we will use those for averaging
elm_vols = gray_matter.elements_volumes_and_areas()[:]

## Plot the ROI
gray_matter.add_element_field(roi, 'roi')
gray_matter.view(visible_fields='roi').show()

## Get field and calculate the mean
# get the field of interest
field_name = 'normE'
field = gray_matter.field[field_name][:]

# Calculate the mean
mean_normE = np.average(field[roi], weights=elm_vols[roi])
print('mean ', field_name, ' in M1 ROI: ', mean_normE)

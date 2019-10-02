'''
Simple ROI analysis of the electric field from a simulation.

We will calculate the mean electric field in a gray matter ROI defined around M1
'''
import numpy as np
import simnibs

## Input ##

# Read the simulation result
head_mesh = simnibs.read_msh('tdcs/ernie_TDCS_1_scalar.msh')

# Crop the mesh so we only have gray matter volume elements (tag 2 in the mesh)
gray_matter = head_mesh.crop_mesh(tags=2)


## Define the ROI ##
# We will now define our ROI

# Define M1 from MNI coordinates (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2034289/)
# the first argument is the MNI coordinates
# the second argument is the segmentation output folder
ernie_M1 = simnibs.mni2subject_coords([-37, -21, 58], 'm2m_ernie/')
# we will use a sphere of radius 10 mm
r = 10.

# Electric fields are defined in the center of the elements
# Therefore, we will select all elements which have their centers inside the ROI
elm_centers = gray_matter.elements_baricenters()[:]
# calculate distances to M1
roi = np.linalg.norm(elm_centers - ernie_M1, axis=1) < r

# finally, calculate the mean of the electric field norm
mean_normE = np.mean(gray_matter.field['normE'][roi])
print('mean normE in M1 ROI: ', mean_normE)


## Alternative: define ROI from EEG electrode position
eeg_pos = simnibs.eeg_positions('m2m_ernie/')

# Find the point in gray matter closest to C3
c3_gm = gray_matter.find_closest_element(eeg_pos['C3'])
# re-define the ROI
roi = np.linalg.norm(elm_centers - c3_gm, axis=1) < r
# finally, calculate the mean of the electric field norm
mean_normE = np.mean(gray_matter.field['normE'][roi])
print('mean normE in C3 ROI: ', mean_normE)


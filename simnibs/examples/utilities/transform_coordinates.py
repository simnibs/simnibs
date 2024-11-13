'''
Transform coordinates from MNI space to subject space and back

Please run with the command

    simnibs_python transform_coordinates.py

from the directory containing the m2m_ernie example dataset

Copyright (c) 2019 SimNIBS developers. Licensed under the GPL v3.
'''

import simnibs


# MNI position for M1 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2034289/)
mni_M1 = [-37, -21, 58]

# Calculate the M1 position in the "ernie" model
# we need to provide the path to the segmentation output folder
ernie_M1 = simnibs.mni2subject_coords(mni_M1, 'm2m_ernie/')
print('Ernie M1 corrdinates:', ernie_M1)
# You can open the ernie_T1fs_conform.nii.gz and check if the coordinates match

# You can also do the inverse transformation (from subject to MNI)
mni_M1 = simnibs.subject2mni_coords(ernie_M1, 'm2m_ernie/')
print('Transforming coordinates back to MNI space:', mni_M1)
# Because SimNIBS uses non-linear transformations, you will not recover the
# original coordinates exactly

'''
    This example wil go throgh simulations and calculate
    the average of the electric field magnitude in MNI space

    It is a follow-up to the "run_simulations" example
'''

import os
import numpy as np
import nibabel as nib

subjects = ['sub01', 'sub09', 'sub10', 'sub12', 'sub15']
results_folder = os.path.join('bipolar', 'mni_volumes')
field_name = 'magnE'

mni_image_suffix = f'_TDCS_1_scalar_MNI_{field_name}.nii.gz'


# Go though each subject and load the images
images = []
for sub in subjects:
    images.append(
        nib.load(os.path.join(sub, results_folder, sub + mni_image_suffix))
    )

# calculate mean
avg = np.mean([img.get_fdata() for img in images], axis=0)

# save results
fn_out = f'mean_{field_name}.nii.gz'
nib.save(nib.Nifti1Image(avg, images[0].affine), fn_out)
print('Average field file: ', fn_out)

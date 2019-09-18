'''
    This example wil go throgh simulations and calculate
    the average and the standard deviation of the normal component
    of the electric field in FsAverage space

    It is a follow-up to the "run_simulations" example
'''
import os
import numpy as np
import simnibs

subjects = ['sub01', 'sub09', 'sub10', 'sub12', 'sub15']
results_folder = os.path.join('bipolar', 'fsavg_overlays')
fsavg_msh_name = '_TDCS_1_scalar_fsavg.msh'

# Go though each subject and extract E_normal
normals = []
for sub in subjects:
    ## read the msh file
    results_fsavg = simnibs.read_msh(
        os.path.join(sub, results_folder, sub + fsavg_msh_name)
    )
    ## extract the normals
    normals.append(results_fsavg.field['E_normal'].value)

# Calculate the average and standard deviation
normals = np.vstack(normals)
avg = np.mean(normals, axis=0)
std = np.std(normals, axis=0)

# Visualize Mean and Std
## cleanup the last model by removing the fields
results_fsavg.nodedata = []
## add the average and standard deviations as nodal data
results_fsavg.add_node_field(avg, 'E_normal_avg')
results_fsavg.add_node_field(std, 'E_normal_std')
view = results_fsavg.view(visible_fields='E_normal_avg')
view.show()

## Calculate the average field in a brain region
# load atlas
atlas = simnibs.get_atlas('a2009s')
region_name = 'lh.S_precentral-sup-part'
roi = atlas[region_name]
# calculate the average in the superior part of the precental sulcus of the left
# hemisphere
node_areas = results_fsavg.nodes_areas()
avg_field_roi = np.average(avg[roi], weights=node_areas[roi])
print('Average field in left superior precental sulcus: ', avg_field_roi)
results_fsavg.add_node_field(roi, region_name)
view = results_fsavg.view(visible_fields=region_name)
view.show()

'''
    This example wil go throgh simulations and calculate
    the average and the standard deviation of the normal component
    of the electric field in FsAverage space

    It is a follow-up to the "run_simulations" example
'''
import os
import numpy as np
import simnibs

## Load simulation results
subjects = ['sub01', 'sub09', 'sub10', 'sub12', 'sub15']
results_folder = 'fsavg_overlays'
fsavg_msh_name = '_TDCS_1_scalar_fsavg.msh'
field_name = 'E_normal'

fields = []
for sub in subjects:
    # read mesh with results transformed to fsaverage space
    results_fsavg = simnibs.read_msh(
        os.path.join('bipolar', sub, results_folder, sub + fsavg_msh_name)
    )
    # save the field in each subject
    fields.append(results_fsavg.field[field_name].value)

## Calculate and plot averages
# Calculate
fields = np.vstack(fields)
avg_field = np.mean(fields, axis=0)
std_field = np.std(fields, axis=0)

# Plot
results_fsavg.nodedata = [] # cleanup fields
results_fsavg.add_node_field(avg_field, 'E_normal_avg') # add average field
results_fsavg.add_node_field(std_field, 'E_normal_std') # add std field

# show surface with the fields
results_fsavg.view(visible_fields='E_normal_avg').show()

## Calculate average in an ROI defined using an atlas
# load atlas and define a region
atlas = simnibs.get_atlas('HCP_MMP1')
region_name = 'lh.4'
roi = atlas[region_name]
# visualize region
results_fsavg.add_node_field(roi, region_name)
results_fsavg.view(visible_fields=region_name).show()

# calculate mean field using a weighted mean
node_areas = results_fsavg.nodes_areas()
avg_field_roi = np.average(avg_field[roi], weights=node_areas[roi])
print(f'Average {field_name} in {region_name}: ', avg_field_roi)
results_fsavg.add_node_field(roi, region_name)

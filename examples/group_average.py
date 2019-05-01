import os
import numpy as np
import simnibs

# Set the subjects
subjects = ['sub01', 'sub09', 'sub10', 'sub12', 'sub15']
results_folder = "bipolar_py/fsavg_overlays"

# Read the normals
normals = []
for sub in subjects:
    # Read the normal component in the left hemisphere
    normal_lh = simnibs.msh.read_curv(
        os.path.join(
            sub, results_folder,
            'lh.' + sub + '_TDCS_1_scalar.fsavg.E.normal'))
    # Read the normatl component in the right hemisphere
    normal_rh = simnibs.msh.read_curv(
        os.path.join(
            sub, results_folder,
            'rh.' + sub + '_TDCS_1_scalar.fsavg.E.normal'))
    # Join the results in a single vector
    joined = np.append(normal_lh, normal_rh)
    normals.append(joined)

# Calculate the average and standard deviation
normals = np.vstack(normals)
avg = np.mean(normals, axis=0)
std = np.std(normals, axis=0)

# Visualize Mean and Std

# Read left and right hemisphere
mesh_lh = simnibs.msh.read_gifti_surface(
    simnibs.file_finder.templates.cat_lh_cortex_ref)
mesh_rh = simnibs.msh.mesh_io.read_gifti_surface(
    simnibs.msh.templates.cat_rh_cortex_ref)
# Join both for a complete gray matter model
gm = mesh_lh.join_mesh(mesh_rh)
# Add the average and standard deviations as nodal data
gm.add_node_field(avg, 'avg_normal')
gm.add_node_field(std, 'std_normal')
view = gm.view(surfaces=False, visible_fields=['avg_normal'])
view.show()

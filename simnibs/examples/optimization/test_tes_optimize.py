''' Example of a SimNIBS tDCS optimization in Python
    Run with:

    simnibs_python teslf_optimize.py

    Written by Konstantin Weise 2022
'''

import os
import h5py
import pynibs
import simnibs
import numpy as np
import matplotlib.pyplot as plt

# location of example data
example_data_folder = os.path.join(simnibs.SIMNIBSDIR, '_internal_resources', 'testing_files')

# fn_mesh = os.path.join(example_data_folder, 'sphere3.msh')
fn_mesh = "/data/pt_01756/probands/15484.08/mesh/charm_beta_fine/m2m_15484.08/15484.08.msh"

# create a circular array with 1 center electrode and 6 outer electrodes
circular_array = simnibs.CircularArray(radius_inner=5, distance=15, n_outer=6, radius_outer=3)

# load mesh
msh = simnibs.read_msh(fn_mesh)

# load roi points
fn_roi = "/data/pt_01756/probands/15484.08/mesh/charm_beta_fine/roi/midlayer_m1s1pmd/geo.hdf5"
with h5py.File(fn_roi, "r") as f:
    points = f["mesh/nodes/node_coord"][:]

# create a region of interest
# points = np.zeros((10, 3))
# points[:, 0] = np.linspace(1, 10, 10)

# msh_cropped = msh.crop_mesh(elm_type=4)
# roi = simnibs.RegionOfInterest(points=points, msh=msh)
roi = None

# plt.spy(roi.sF)
# Initialize structure
opt = simnibs.opt_struct.TESoptimize(msh=msh,
                                     roi=roi,
                                     electrode=circular_array,
                                     init_pos=["C3"],
                                     output_folder="/data/pt_01756/probands/15484.08/opt/ttf",
                                     plot=True)

# pynibs.plot_surface(data=np.arange(len(opt.skin_surface.surf2msh_triangles)),
#                     con=opt.mesh.elm.node_number_list[opt.skin_surface.surf2msh_triangles, :3],
#                     points=opt.mesh.nodes.node_coord)

# write hdf5 _geo file
# pynibs.write_geo_hdf5_surf(out_fn="/home/kporzig/tmp/test_geo.hdf5",
#                            points=opt.mesh.nodes.node_coord,
#                            con=opt.mesh.elm.node_number_list[opt.skin_surface.surf2msh_triangles, :3]-1,
#                            replace=True,
#                            hdf5_path='/mesh')
#
# pynibs.write_data_hdf5_surf(data=[np.arange(len(opt.skin_surface.surf2msh_triangles))],
#                             data_names=["test"],
#                             data_hdf_fn_out="/home/kporzig/tmp/test_data.hdf5",
#                             geo_hdf_fn="/home/kporzig/tmp/test_geo.hdf5",
#                             replace=True)


# # Select a name for the optimization
# opt.name = 'optimization/single_target'
#
# # Select a maximum total current (in A)
# opt.max_total_current = 2e-3
# # Select a maximum current at each electrodes (in A)
# opt.max_individual_current = 1e-3
# # Select a maximum number of active electrodes (optional)
# opt.max_active_electrodes = 8
#
# # Define optimization target
# target = opt.add_target()
# # Position of target, in subject space!
# # please see teslf_optimize_mni.py for how to use MNI coordinates
# target.positions = [-50.7, 5.1, 55.5]
# # Intensity of the electric field (in V/m)
# target.intensity = 0.2
#
# # Run optimization
# simnibs.run_simnibs(opt)

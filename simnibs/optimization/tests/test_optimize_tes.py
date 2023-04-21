''' Example of a SimNIBS tDCS optimization in Python
    Run with:

    simnibs_python teslf_optimize.py

    Written by Konstantin Weise 2022
'''

import os
import h5py
import pynibs
import simnibs
import nibabel
import numpy as np
import matplotlib.pyplot as plt

from simnibs.optimization.optimize_tes import create_new_connectivity_list_point_mask

import matplotlib

# matplotlib.use("Qt5Agg")

output_folder = "/data/pt_01756/studies/ttf"
# output_folder = "/home/kporzig/tmp"

fn_mesh = "/data/pt_01756/probands/15484.08/mesh/charm_beta_coarse/m2m_15484.08/15484.08.msh"
# fn_mesh = "/home/kporzig/tmp/charm_beta_coarse/m2m_15484.08/15484.08.msh"

fn_roi_1 = "/data/pt_01756/studies/ttf/roi_0_geo.hdf5"
fn_roi_2 = "/data/pt_01756/studies/ttf/roi_1_geo.hdf5"
# fn_roi = "/data/pt_01756/probands/15484.08/mesh/charm_beta_coarse/roi/midlayer_m1s1pmd/geo.hdf5"
# fn_roi = "/home/kporzig/tmp/charm_beta_coarse/roi/midlayer_m1s1pmd/geo.hdf5"

# location of example data
example_data_folder = os.path.join(simnibs.SIMNIBSDIR, '_internal_resources', 'testing_files')

optimizer = "differential_evolution"  # "direct"  "Nelder-Mead"  "differential_evolution" "shgo"
optimize_init_vals = True
goal = ["mean", "mean"]  # "mean" "mean_max_TI"
dataType = [0, 0]
weights = [1., 0.]
constrain_electrode_locations = False
overlap_factor = 1.
polish = False
locally_biased = True
init_pos = None
current_estimator_method = "gpc"
dirichlet_correction_detailed = False
# init_pos = ["C3"]
# init_pos = ["C3", "C4"]

print("Initializing Electrode ...")
# create a circular array with 1 center electrode and 6 outer electrodes
########################################################################################################################
# electrode = simnibs.CircularArray(radius_inner=10, distance=40, n_outer=4, radius_outer=10,
#                                   current_estimator_method=current_estimator_method,
#                                   dirichlet_correction_detailed=dirichlet_correction_detailed)

# create 1 x 1 standard TES montage
########################################################################################################################
# center = np.array([[0, 0, 0]])
# radius = np.array([0])
# length_x = np.array([40])
# length_y = np.array([40])
# electrode = simnibs.ElectrodeArrayPair(center=center, radius=radius, length_x=length_x, length_y=length_y,
#                                        current_estimator_method=current_estimator_method,
#                                        dirichlet_correction_detailed=dirichlet_correction_detailed)

# create 2 channel of 1 x 1 standard TES montage (Temporal interference)
########################################################################################################################
# center = np.array([[0, 0, 0]])
# radius = np.array([0])
# length_x = np.array([40])
# length_y = np.array([40])
#
# electrode_1 = simnibs.ElectrodeArrayPair(center=center, radius=radius, length_x=length_x, length_y=length_y, current_estimator_method=current_estimator_method)
# electrode_2 = simnibs.ElectrodeArrayPair(center=center, radius=radius, length_x=length_x, length_y=length_y, current_estimator_method=current_estimator_method)
#
# electrode = [electrode_1, electrode_2]

# create 3 x 3 circular electrode array pair
########################################################################################################################
# center = np.array([[-30, 20, 0],
#                      [0, 20, 0],
#                      [30, 20, 0],
#                      [-30, 0, 0],
#                      [0, 0, 0],
#                      [30, 0, 0],
#                      [-30, -20, 0],
#                      [0, -20, 0],
#                      [30, -20, 0]])
# radius = np.array([7, 7, 7, 7, 7, 7, 7, 7, 7])
# length_x = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
# length_y = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
# electrode = simnibs.ElectrodeArrayPair(center=center, radius=radius, length_x=length_x, length_y=length_y,
#                                        current_estimator_method=current_estimator_method)

# create two 3 x 3 circular electrode array pairs (2 channels)
########################################################################################################################
center = np.array([[-30, 0, 0],
                   [0, 0, 0],
                   [30, 0, 0]])
radius = np.array([7, 7, 7])
length_x = np.array([0, 0, 0])
length_y = np.array([0, 0, 0])
center = np.array([[-30, 20, 0],
                    [0, 20, 0],
                    [30, 20, 0],
                    [-30, 0, 0],
                    [0, 0, 0],
                    [30, 0, 0],
                    [-30, -20, 0],
                    [0, -20, 0],
                    [30, -20, 0]])
radius = np.array([7, 7, 7, 7, 7, 7, 7, 7, 7])
length_x = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
length_y = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
electrode_1 = simnibs.ElectrodeArrayPair(center=center, radius=radius, length_x=length_x, length_y=length_y,
                                         current_estimator_method=current_estimator_method,
                                         dirichlet_correction_detailed=dirichlet_correction_detailed)
electrode_2 = simnibs.ElectrodeArrayPair(center=center, radius=radius, length_x=length_x, length_y=length_y,
                                         current_estimator_method=current_estimator_method,
                                         dirichlet_correction_detailed=dirichlet_correction_detailed)

electrode = [electrode_1, electrode_2]

# create 3 x 3 mixed circular and rectangular electrode array pair
#######################################################################################################################
#center = np.array([[-30, 20, 0],
#                   [0, 20, 0],
#                   [30, 20, 0],
#                   [-30, 0, 0],
#                   [0, 0, 0],
#                   [30, 0, 0],
#                   [-30, -20, 0],
#                   [0, -20, 0],
#                   [30, -20, 0]])
#radius = np.array([0, 0, 0, 7, 7, 7, 0, 0, 0])
#length_x = np.array([20, 20, 20, 0, 0, 0, 20, 20, 20])
#length_y = np.array([20, 20, 20, 0, 0, 0, 20, 20, 20])
#electrode = simnibs.ElectrodeArrayPair(center=center, radius=radius, length_x=length_x, length_y=length_y, current_estimator_method=current_estimator_method)

# load mesh
mesh = simnibs.read_msh(fn_mesh)

# # load roi points
# with h5py.File(fn_roi, "r") as f:
#     points_1 = f["mesh/nodes/node_coord"][:]
#     con_1 = f["mesh/elm/triangle_number_list"][:]
#
# lh = nibabel.load("/data/pt_01756/probands/15484.08/mesh/charm_beta_coarse/m2m_15484.08/surf/lh.central.gii")
# rh = nibabel.load("/data/pt_01756/probands/15484.08/mesh/charm_beta_coarse/m2m_15484.08/surf/rh.central.gii")
# lh_points, lh_con = lh.agg_data()
# rh_points, rh_con = rh.agg_data()
# points_2 = np.vstack((lh_points, rh_points))
# con_2 = np.vstack((lh_con, rh_con + lh_points.shape[0]))
# # tri_center = np.average(points[con, ], axis=1)
#
# # remove points from ROI #2 (non-roi) which are close to ROI #1 (roi)
# mask = np.ones(points_2.shape[0]).astype(bool)
# for i_p2, _p2 in enumerate(points_2):
#     if (np.linalg.norm(points_1 - _p2, axis=1) < 0.1).any():
#         mask[i_p2] = False
#
# points_2_masked, con_2_masked = create_new_connectivity_list_point_mask(points=points_2,
#                                                                         con=con_2,
#                                                                         point_mask=mask)

with h5py.File(fn_roi_1, "r") as f:
    points_1 = f["mesh/nodes/node_coord"][:]
    con_1 = f["mesh/elm/triangle_number_list"][:]

with h5py.File(fn_roi_2, "r") as f:
    points_2 = f["mesh/nodes/node_coord"][:]
    con_2 = f["mesh/elm/triangle_number_list"][:]

# calculate ROI volumes (areas)
p1_tri_1 = points_1[con_1[:, 0], :]
p2_tri_1 = points_1[con_1[:, 1], :]
p3_tri_1 = points_1[con_1[:, 2], :]
vol_1 = np.sum(0.5 * np.linalg.norm(np.cross(p2_tri_1 - p1_tri_1, p3_tri_1 - p1_tri_1), axis=1))

p1_tri_2 = points_2[con_2[:, 0], :]
p2_tri_2 = points_2[con_2[:, 1], :]
p3_tri_2 = points_2[con_2[:, 2], :]
vol_2 = np.sum(0.5 * np.linalg.norm(np.cross(p2_tri_2 - p1_tri_2, p3_tri_2 - p1_tri_2), axis=1))

# import pynibs
# # write hdf5 _geo file
# pynibs.write_geo_hdf5_surf(out_fn="/data/pt_01756/studies/ttf/plots/e_roi_0_geo.hdf5",
#                            points=points_1,
#                            con=con_1,
#                            replace=True,
#                            hdf5_path='/mesh')
#
# pynibs.write_data_hdf5_surf(data=[np.zeros(con_1.shape[0])],
#                             data_names=["test"],
#                             data_hdf_fn_out="/data/pt_01756/studies/ttf/plots/e_roi_0_data_test.hdf5",
#                             geo_hdf_fn="/data/pt_01756/studies/ttf/plots/e_roi_0_geo.hdf5",
#                             replace=True)
#
# pynibs.write_geo_hdf5_surf(out_fn="/data/pt_01756/studies/ttf/plots/e_roi_1_geo.hdf5",
#                            points=points_2,
#                            con=con_2,
#                            replace=True,
#                            hdf5_path='/mesh')
#
# pynibs.write_data_hdf5_surf(data=[np.zeros(con_2.shape[0])],
#                             data_names=["test"],
#                             data_hdf_fn_out="/data/pt_01756/studies/ttf/plots/e_roi_1_data_test.hdf5",
#                             geo_hdf_fn="/data/pt_01756/studies/ttf/plots/e_roi_1_geo.hdf5",
#                             replace=True)

# create a region of interest
print("Initializing ROI #1 ...")
roi_1 = simnibs.RegionOfInterest(points=points_1, con=con_1, mesh=mesh, vol=vol_1)

print("Initializing ROI #2 ...")
roi_2 = simnibs.RegionOfInterest(points=points_2, con=con_2, mesh=mesh, vol=vol_2)

roi = [roi_1, roi_2]

min_electrode_distance = 10
optimizer_options = {"maxiter": 1000,
                     "disp": True,
                     "recombination": 0.3,              # differential evolution
                     "mutation": (0.01, 0.1),           # differential evolution
                     "popsize": 15,                     # differential evolution
                     "tol": 0.5,
                     "locally_biased": locally_biased}           # DIRECT

# Initialize TESoptimize class
opt = simnibs.opt_struct.TESoptimize(mesh=mesh,
                                     roi=roi,
                                     electrode=electrode,
                                     init_pos=init_pos,
                                     output_folder=output_folder,
                                     plot=True,
                                     optimizer=optimizer,
                                     optimizer_options=optimizer_options,
                                     min_electrode_distance=min_electrode_distance,
                                     weights=weights,
                                     goal=goal,  # mean"
                                     constrain_electrode_locations=constrain_electrode_locations,
                                     overlap_factor=overlap_factor,
                                     optimize_init_vals=optimize_init_vals,
                                     dataType=dataType,
                                     polish=polish)

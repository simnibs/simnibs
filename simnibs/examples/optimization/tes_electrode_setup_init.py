import simnibs
import numpy as np

# create 3 x 3 circular electrode array pair
center = np.array([[-20,  20, 0],
                   [  0,  20, 0],
                   [ 20,  20, 0],
                   [-20,   0, 0],
                   [  0,   0, 0],
                   [ 20,   0, 0],
                   [-20, -20, 0],
                   [  0, -20, 0],
                   [ 20, -20, 0]])
radius   = np.array([ 0,  5,  0,  5,  5,  5,  0,  5,  0])
length_x = np.array([10,  0, 10,  0,  0,  0, 10,  0, 10])
length_y = np.array([10,  0, 10,  0,  0,  0, 10,  0, 10])
electrode_array_pair = simnibs.ElectrodeArrayPair(center=center, radius=radius, length_x=length_x, length_y=length_y)
electrode_array_pair.electrode_arrays[0].plot(show=True, fn_plot="/home/kporzig/tmp/electrode_array_pair.png")

# update geometry (keep centers untouched)
radius   = np.array([ 5,  0,  5,  0,  0,  0,  5,  0, 5])
length_x = np.array([ 0, 10,  0, 10, 10, 10,  0, 10, 0])
length_y = np.array([ 0, 10,  0, 10, 10, 10,  0, 10, 0])
electrode_array_pair.update_geometry(center=None, radius=radius, length_x=length_x, length_y=length_y)
electrode_array_pair.electrode_arrays[0].plot(show=True, fn_plot="/home/kporzig/tmp/electrode_array_pair_modified.png")

# create a circular array with 1 center electrode and 6 outer electrodes
circular_array = simnibs.CircularArray(radius_inner=5, distance=15, n_outer=6, radius_outer=3)
circular_array.electrode_arrays[0].plot(fn_plot="/home/kporzig/tmp/circular_array.png")

# update geometry
circular_array.update_geometry(distance=10)
circular_array.electrode_arrays[0].plot(fn_plot="/home/kporzig/tmp/circular_array_modified.png")

# create 3 x 3 mixed electrode array (circular and rectangular)
center = np.array([[-20,  25, 0],
                   [  0,  22, 0],
                   [ 20,  20, 0],
                   [-20,   0, 0],
                   [  0,   0, 0],
                   [ 20,   0, 0],
                   [-20, -30, 0],
                   [  0, -25, 0],
                   [ 20, -20, 0]])
radius   = np.array([ 0,  3,  0,  5,  7, 9, 0, 3, 0])
length_x = np.array([10,  0,  5,  0,  0, 0, 8, 0, 5])
length_y = np.array([5,  0,  10,  0,  0, 0, 8, 0, 5])
electrode_array_mixed = simnibs.ElectrodeArray(channel_id=np.arange(center.shape[0]),
                                               center=center,
                                               radius=radius,
                                               length_x=length_x,
                                               length_y=length_y)

# plot array
electrode_array_mixed.plot(show=True, fn_plot="/home/kporzig/tmp/electrode_array_mixed.png")

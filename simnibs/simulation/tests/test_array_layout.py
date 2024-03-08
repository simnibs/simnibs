import pytest
import simnibs
import numpy as np


class TestElectrode:
    def test_electrode_array_pair(self):
        # create 3 x 3 circular electrode array pair
        center = np.array([[-20,  20],
                           [  0,  20],
                           [ 20,  20],
                           [-20,   0],
                           [  0,   0],
                           [ 20,   0],
                           [-20, -20],
                           [  0, -20],
                           [ 20, -20]])
        radius   = np.array([ 0,  5,  0,  5,  5,  5,  0,  5,  0])
        length_x = np.array([10,  0, 10,  0,  0,  0, 10,  0, 10])
        length_y = np.array([10,  0, 10,  0,  0,  0, 10,  0, 10])
        electrode_array_pair = simnibs.ElectrodeArrayPair(center=center, radius=radius, length_x=length_x, length_y=length_y)
        # electrode_array_pair.electrode_arrays[0].plot(show=False, fn_plot="/home/kporzig/tmp/electrode_array_pair.png")

        # update geometry (keep centers untouched)
        radius   = np.array([ 5,  0,  5,  0,  0,  0,  5,  0, 5])
        length_x = np.array([ 0, 10,  0, 10, 10, 10,  0, 10, 0])
        length_y = np.array([ 0, 10,  0, 10, 10, 10,  0, 10, 0])
        electrode_array_pair.update_geometry(center=None, radius=radius, length_x=length_x, length_y=length_y)
        # electrode_array_pair.electrode_arrays[0].plot(show=False, fn_plot="/home/kporzig/tmp/electrode_array_pair_modified.png")

    def test_circular_array(self):
        # create a circular array with 1 center electrode and 6 outer electrodes
        circular_array = simnibs.CircularArray(radius_inner=5, distance=15, n_outer=6, radius_outer=3)
        # circular_array.electrode_arrays[0].plot(show=False, fn_plot="/home/kporzig/tmp/circular_array.png")

        # update geometry
        circular_array.update_geometry(distance=10)
        # circular_array.electrode_arrays[0].plot(show=False, fn_plot="/home/kporzig/tmp/circular_array_modified.png")

    def test_mixed_array(self):
        # create 3 x 3 mixed electrode array (circular and rectangular)
        center = np.array([[-20,  25],
                           [  0,  22],
                           [ 20,  20],
                           [-20,   0],
                           [  0,   0],
                           [ 20,   0],
                           [-20, -30],
                           [  0, -25],
                           [ 20, -20]])
        radius   = np.array([ 0,  3,  0,  5,  7, 9, 0, 3, 0])
        length_x = np.array([10,  0,  5,  0,  0, 0, 8, 0, 5])
        length_y = np.array([5,  0,  10,  0,  0, 0, 8, 0, 5])
        electrode_array_mixed = simnibs.ElectrodeArray(channel_id=np.arange(center.shape[0]),
                                                       center=center,
                                                       radius=radius,
                                                       length_x=length_x,
                                                       length_y=length_y)
        # electrode_array_mixed.plot(show=True, fn_plot="/home/kporzig/tmp/electrode_array_mixed.png")

    def test_ElectrodeInitializer_array_pair(self):
        electrode_i = simnibs.ElectrodeInitializer()
        electrode_i.type = "ElectrodeArrayPair"             # Pair of TES electrodes
        electrode_i.center = [[0, 0]]                       # electrode center in reference electrode space (x-y plane)
        electrode_i.radius = [10]                           # radius of electrodes
        electrode_i.dirichlet_correction_detailed = False   # node wise dirichlet correction
        electrode_i.current = [0.002, -0.002]               # electrode currents
        electrode = electrode_i.initialize()

    def test_ElectrodeInitializer_array_pair_geo_opt_1(self):
        electrode_i = simnibs.ElectrodeInitializer()
        electrode_i.type = "ElectrodeArrayPair"             # Pair of TES electrodes
        electrode_i.center = [[0, 0]]                       # electrode center in reference electrode space (x-y plane)
        electrode_i.radius_bounds = [10, 20]                # radius of electrodes
        electrode_i.dirichlet_correction_detailed = False   # node wise dirichlet correction
        electrode_i.current = [0.002, -0.002]               # electrode currents
        electrode = electrode_i.initialize()


    def test_ElectrodeInitializer_array_pair_geo_opt_2(self):
        electrode_i = simnibs.ElectrodeInitializer()
        electrode_i.type = "ElectrodeArrayPair"             # Pair of TES electrodes
        electrode_i.center = [[0, 0]]                       # electrode center in reference electrode space (x-y plane)
        electrode_i.length_x_bounds = [10, 20]              # x-dimension of electrodes
        electrode_i.length_y_bounds = [10, 20]              # y-dimension of electrodes
        electrode_i.dirichlet_correction_detailed = False   # node wise dirichlet correction
        electrode_i.current = [0.002, -0.002]               # electrode currents
        electrode = electrode_i.initialize()

    def test_ElectrodeInitializer_circular_array(self):
        electrode_i = simnibs.ElectrodeInitializer()
        electrode_i.type = "CircularArray"          # HDTES center surround montage
        electrode_i.radius_inner = 10               # radius of inner electrode
        electrode_i.radius_outer = 10               # radius of outer electrodes
        electrode_i.distance = 25                   # distance bounds between inner and outer electrodes
        electrode_i.n_outer = 4                     # number of outer electrodes
        electrode_i.dirichlet_correction = False    # electrode wise dirichlet correction
        electrode_i.dirichlet_correction_detailed = False  # node wise dirichlet correction
        electrode_i.current = [0.002, -0.002 / 4, -0.002 / 4, -0.002 / 4, -0.002 / 4]  # initial currents
        electrode = electrode_i.initialize()

    def test_ElectrodeInitializer_circular_array_geo_opt(self):
        electrode_i = simnibs.ElectrodeInitializer()
        electrode_i.type = "CircularArray"          # HDTES center surround montage
        electrode_i.radius_inner_bounds = [10, 20]  # radius of inner electrode
        electrode_i.radius_outer_bounds = [10, 20]  # radius of outer electrodes
        electrode_i.distance_bounds = [10, 20]      # distance bounds between inner and outer electrodes
        electrode_i.n_outer_bounds = [4, 6]         # number of outer electrodes
        electrode_i.dirichlet_correction = False    # electrode wise dirichlet correction
        electrode_i.dirichlet_correction_detailed = False  # node wise dirichlet correction
        electrode_i.current = 0.002                 # current
        electrode = electrode_i.initialize()

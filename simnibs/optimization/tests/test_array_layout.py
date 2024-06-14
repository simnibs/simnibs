import os
import pytest
import scipy
import simnibs
import numpy as np

from simnibs.optimization.tes_flex_optimization.electrode_layout import CircularArray, ElectrodeArray, ElectrodeArrayPair
from simnibs.utils.matlab_read import dict_from_matlab


class TestElectrode:
    def test_electrode_array_pair(self):
        # create 3 x 3 circular electrode array pair
        electrode_array_pair = ElectrodeArrayPair()

        electrode_array_pair.center = np.array([[-20,  20],
                           [  0,  20],
                           [ 20,  20],
                           [-20,   0],
                           [  0,   0],
                           [ 20,   0],
                           [-20, -20],
                           [  0, -20],
                           [ 20, -20]])
        electrode_array_pair.radius   = np.array([ 0,  5,  0,  5,  5,  5,  0,  5,  0])
        electrode_array_pair.length_x = np.array([10,  0, 10,  0,  0,  0, 10,  0, 10])
        electrode_array_pair.length_y = np.array([10,  0, 10,  0,  0,  0, 10,  0, 10])
        electrode_array_pair._prepare()
        # electrode_array_pair.electrode_arrays[0].plot(show=False, fn_plot="/home/kporzig/tmp/electrode_array_pair.png")

        # update geometry (keep centers untouched)
        radius   = np.array([ 5,  0,  5,  0,  0,  0,  5,  0, 5])
        length_x = np.array([ 0, 10,  0, 10, 10, 10,  0, 10, 0])
        length_y = np.array([ 0, 10,  0, 10, 10, 10,  0, 10, 0])
        electrode_array_pair.update_geometry(center=None, radius=radius, length_x=length_x, length_y=length_y)
        # electrode_array_pair.electrode_arrays[0].plot(show=False, fn_plot="/home/kporzig/tmp/electrode_array_pair_modified.png")

    def test_circular_array(self):
        # create a circular array with 1 center electrode and 6 outer electrodes
        circular_array = CircularArray()
        circular_array.radius_inner=5
        circular_array.distance=15
        circular_array.n_outer=6
        circular_array.radius_outer=3
        circular_array._prepare()
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
        electrode_array_mixed = ElectrodeArray(channel_id=np.arange(center.shape[0]),
                                                       center=center,
                                                       radius=radius,
                                                       length_x=length_x,
                                                       length_y=length_y)
        # electrode_array_mixed.plot(show=True, fn_plot="/home/kporzig/tmp/electrode_array_mixed.png")

    def test_ElectrodeInitializer_array_pair(self, tmp_path):
        electrode_i = ElectrodeArrayPair()
        electrode_i.center = [[0, 0]]                       # electrode center in reference electrode space (x-y plane)
        electrode_i.radius = [10]                           # radius of electrodes
        electrode_i.dirichlet_correction_detailed = False   # node wise dirichlet correction
        electrode_i.current = [0.002, -0.002]               # electrode currents

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, electrode_i.to_dict())
        electrode_loaded = ElectrodeArrayPair(dict_from_matlab(scipy.io.loadmat(mat_path)))
        electrode_i._prepare()
        electrode_loaded._prepare()
        np.testing.assert_equal(electrode_i.to_dict(), electrode_loaded.to_dict())


    def test_ElectrodeInitializer_array_pair_geo_opt_1(self, tmp_path):
        electrode_i = ElectrodeArrayPair()
        electrode_i.center = [[0, 0]]                       # electrode center in reference electrode space (x-y plane)
        electrode_i.radius_bounds = [10, 20]                # radius of electrodes
        electrode_i.dirichlet_correction_detailed = False   # node wise dirichlet correction
        electrode_i.current = [0.002, -0.002]               # electrode currents
        
        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, electrode_i.to_dict())
        electrode_loaded = ElectrodeArrayPair(dict_from_matlab(scipy.io.loadmat(mat_path)))
        electrode_i._prepare()
        electrode_loaded._prepare()
        np.testing.assert_equal(electrode_i.to_dict(), electrode_loaded.to_dict())


    def test_ElectrodeInitializer_array_pair_geo_opt_2(self, tmp_path):
        electrode_i = ElectrodeArrayPair()
        electrode_i.center = [[0, 0]]                       # electrode center in reference electrode space (x-y plane)
        electrode_i.length_x_bounds = [10, 20]              # x-dimension of electrodes
        electrode_i.length_y_bounds = [10, 20]              # y-dimension of electrodes
        electrode_i.dirichlet_correction_detailed = False   # node wise dirichlet correction
        electrode_i.current = [0.002, -0.002]               # electrode currents
        
        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, electrode_i.to_dict())
        electrode_loaded = ElectrodeArrayPair(dict_from_matlab(scipy.io.loadmat(mat_path)))
        electrode_i._prepare()
        electrode_loaded._prepare()
        np.testing.assert_equal(electrode_i.to_dict(), electrode_loaded.to_dict())

    def test_ElectrodeInitializer_circular_array(self, tmp_path):
        electrode_i = CircularArray()
        electrode_i.radius_inner = 10               # radius of inner electrode
        electrode_i.radius_outer = 10               # radius of outer electrodes
        electrode_i.distance = 25                   # distance bounds between inner and outer electrodes
        electrode_i.n_outer = 4                     # number of outer electrodes
        electrode_i.dirichlet_correction = False    # electrode wise dirichlet correction
        electrode_i.dirichlet_correction_detailed = False  # node wise dirichlet correction
        electrode_i.current = [0.002, -0.002 / 4, -0.002 / 4, -0.002 / 4, -0.002 / 4]  # initial currents

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, electrode_i.to_dict())
        electrode_loaded = CircularArray(dict_from_matlab(scipy.io.loadmat(mat_path)))
        electrode_i._prepare()
        electrode_loaded._prepare()
        np.testing.assert_equal(electrode_i.to_dict(), electrode_loaded.to_dict())

    def test_ElectrodeInitializer_circular_array_geo_opt(self, tmp_path):
        electrode_i = CircularArray()
        electrode_i.radius_inner_bounds = [10, 20]  # radius of inner electrode
        electrode_i.radius_outer_bounds = [10, 20]  # radius of outer electrodes
        electrode_i.distance_bounds = [10, 20]      # distance bounds between inner and outer electrodes
        electrode_i.n_outer_bounds = [4, 6]         # number of outer electrodes
        electrode_i.dirichlet_correction = False    # electrode wise dirichlet correction
        electrode_i.dirichlet_correction_detailed = False  # node wise dirichlet correction
        electrode_i.current = 0.002                 # current

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, electrode_i.to_dict())
        electrode_loaded = CircularArray(dict_from_matlab(scipy.io.loadmat(mat_path)))
        electrode_i._prepare()
        electrode_loaded._prepare()
        np.testing.assert_equal(electrode_i.to_dict(), electrode_loaded.to_dict())



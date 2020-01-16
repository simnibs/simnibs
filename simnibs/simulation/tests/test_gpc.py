import numpy as np
import h5py
import os
import tempfile
from mock import Mock, patch, call

import pytest

from simnibs import SIMNIBSDIR
import simnibs.simulation.gpc as simnibs_gpc
import simnibs.simulation.sim_struct as sim_struct
import simnibs.msh.mesh_io as mesh_io


@pytest.fixture(scope='module')
def sphere3():
    return mesh_io.read_msh(
        os.path.join(
            SIMNIBSDIR, 'resources',
            'testing_files', 'sphere3.msh'))

@pytest.fixture
def cube_msh():
    fn = os.path.join(SIMNIBSDIR, 'resources', 'testing_files', 'cube_w_electrodes.msh')
    return mesh_io.read_msh(fn)



@pytest.fixture(scope='module')
def gpc_regression_instance():
    pdftype = ['beta', 'normal']
    pdfshape = [[1, 0], [1, 2]]
    limits = [[-2, None], [2, None]]
    poly_idx = np.array([[0, 0], [1, 0], [0, 1]])
    coords_norm = np.array([[0.5, 0.2], [0.1, 0.9]])
    regobj = simnibs_gpc.gPC_regression([1, 2],
                                        pdftype, pdfshape,
                                        limits, poly_idx,
                                        coords_norm, 'TMS')
    return regobj

class TestHDF5:
    def test_write_hdf5(self):
        if os.path.isfile('test.hdf5'): os.remove('test.hdf5')
        simnibs_gpc.write_data_hdf5(np.ones((10, 1)), 'ones', 'test.hdf5', 'path/to/')
        with h5py.File('test.hdf5', 'r') as f:
            data = f['path/to/ones'][()]
        os.remove('test.hdf5')
        assert np.all(data == 1)

    def test_read_hdf5(self):
        if os.path.isfile('test.hdf5'): os.remove('test.hdf5')
        with h5py.File('test.hdf5', 'a') as f:
            f.create_dataset('path/to/ones', data=np.ones((10, 1)))
        data = simnibs_gpc.read_data_hdf5('ones', 'test.hdf5', 'path/to/')
        os.remove('test.hdf5')
        assert np.all(data == 1)

class TestgPC_Regression:
    def test_gPC_Regression_initialization(self, gpc_regression_instance):
        pdftype = ['beta', 'normal']
        pdfshape = [[1, 0], [1, 2]]
        limits = [[-2, None], [2, None]]
        poly_idx = np.array([[0, 0], [1, 0], [0, 1]])
        coords_norm = np.array([[0.5, 0.2], [0.1, 0.9]])
        assert gpc_regression_instance.pdftype == pdftype
        assert gpc_regression_instance.pdfshape == pdfshape
        assert gpc_regression_instance.limits == limits
        np.testing.assert_equal(gpc_regression_instance.poly_idx, poly_idx)
        np.testing.assert_equal(gpc_regression_instance.grid.coords_norm, coords_norm)
        np.testing.assert_almost_equal(gpc_regression_instance.grid.coords, coords_norm*2)

    def test_save_hdf5(self, gpc_regression_instance):
        if os.path.isfile('test.hdf5'):
            os.remove('test.hdf5')

        gpc_regression_instance.save_hdf5('test.hdf5')
        with h5py.File('test.hdf5', 'r') as f:
            assert np.all(f['gpc_object/random_vars'][()] == [b"1", b"2"])
            assert np.all(f['gpc_object/pdftype'][()] == [b"beta", b"normal"])
            assert np.allclose(f['gpc_object/pdfshape'][()], np.array([[1, 0],[1, 2]]))
            assert np.allclose(f['gpc_object/limits'], np.array([[-2, -1e10], [2, 1e10]]))
            assert np.allclose(f['gpc_object/poly_idx'],
                               np.array(np.array([[0, 0], [1, 0], [0, 1]])))
            assert np.allclose(f['gpc_object/grid/coords_norm'],
                               np.array([[0.5, 0.2], [0.1, 0.9]]))

        os.remove('test.hdf5')

    def test_load_hdf5(self, gpc_regression_instance):
        if os.path.isfile('test.hdf5'):
            os.remove('test.hdf5')
        gpc_regression_instance.save_hdf5('test.hdf5')
        gpc_regression = simnibs_gpc.gPC_regression.read_hdf5('test.hdf5')
        assert gpc_regression.random_vars == [1, 2]
        assert gpc_regression_instance.pdftype == gpc_regression.pdftype
        assert gpc_regression_instance.pdfshape == gpc_regression.pdfshape
        assert gpc_regression_instance.limits == gpc_regression.limits
        assert np.allclose(gpc_regression_instance.poly_idx,
                           gpc_regression.poly_idx)
        assert np.allclose(gpc_regression_instance.grid.coords_norm,
                           gpc_regression.grid.coords_norm)

        os.remove('test.hdf5')

    def test_postprocessing(self, gpc_regression_instance, sphere3):
        fn = 'gpc_testing.hdf5'
        if os.path.isfile(fn):
            os.remove(fn)
        msh = sphere3.crop_mesh(elm_type=4)
        msh.write_hdf5(fn, 'mesh_roi/')
        # "normalized" coordinates
        coords_norm = np.array([[-.5], [0], [.7]])
        # Define regression object.
        # The random variable is the conductivity of the layer 3
        # Uniform distribution, between 0 and 10
        gpc_reg = \
            simnibs_gpc.gPC_regression([3],
                                       ['beta'], [[1], [1]],
                                       [[0], [10]], [[0], [1]],
                                       coords_norm, 'TCS',
                                       data_file=fn)
        gpc_reg.regularization_factors = [0]
        # Create potentials
        pot = np.tile(msh.nodes.node_coord[None, :, 0],
                      (coords_norm.shape[0], 1))
        lst = sim_struct.SimuList()
        lst.cond[3].value = 1   # Layer 4 conductivity
        lst.cond[4].value = 10  # Layer 5 conductivity
        lst._write_conductivity_to_hdf5(fn)
        for i, c in enumerate(gpc_reg.grid.coords):
            pot[i] *= c  # Linear relationship between potential and random variable
        #  Record potentials
        with h5py.File(fn, 'a') as f:
            f.create_group('mesh_roi/data_matrices')
            f['mesh_roi/data_matrices'].create_dataset('v_samples', data=pot)
        #  Postprocess
        gpc_reg.postprocessing('eEjJ')
        with h5py.File(fn, 'r') as f:
            # Linear relationship between random variable and potential, uniform
            # distribution
            mean_E = f['mesh_roi/elmdata/E_mean'][()]
            mean_E_expected = -np.array([5, 0, 0]) * 1e3
            assert np.allclose(mean_E, mean_E_expected)
            mean_e = f['mesh_roi/elmdata/normE_mean'][()]
            assert np.allclose(mean_e, 5e3)
            # J is a bit more complicated
            mean_J = f['mesh_roi/elmdata/J_mean'][()]
            assert np.allclose(mean_J[msh.elm.tag1==4],
                               mean_E_expected)
            assert np.allclose(mean_J[msh.elm.tag1==5],
                               10 * mean_E_expected)
            coeffs = gpc_reg.expand(gpc_reg.grid.coords ** 2, return_error=False)
            mean = gpc_reg.mean(coeffs)
            assert np.allclose(mean_J[msh.elm.tag1==3],
                               mean * -np.array([1, 0, 0]) * 1e3)
            dsets = f['mesh_roi/elmdata/'].keys()

            std_E = f['mesh_roi/elmdata/E_std'][()]
            assert np.allclose(std_E, np.array([np.sqrt(1e8 / 12), 0., 0.]))
            assert 'normE_std' in dsets
            assert 'J_std' in dsets

            assert 'E_sobol_3' in dsets
            assert 'normE_sobol_3' in dsets
            assert 'J_sobol_3' in dsets

            assert 'E_sensitivity_3' in dsets
            assert 'normE_sensitivity_3' in dsets
            assert 'J_sensitivity_3' in dsets

        os.remove(fn)



@pytest.fixture
def sampler_args(sphere3):
    poslist = sim_struct.SimuList()

    poslist.cond[2].name = 'inner'
    poslist.cond[2].distribution_type = 'beta'
    poslist.cond[2].distribution_parameters = [2, 3, 0.3, 0.4]
    poslist.cond[3].name = 'middle'
    poslist.cond[3].value = 1
    poslist.cond[4].name = 'outer'
    poslist.cond[4].value = 10
    with tempfile.NamedTemporaryFile(suffix='.hdf5') as f:
        fn_hdf5 = f.name
    if os.path.isfile(fn_hdf5):
        os.remove(fn_hdf5)

    return sphere3, poslist, fn_hdf5, [3]


class TestSampler:
    def test_set_up_sampler(self, sampler_args):
        mesh, poslist, fn_hdf5, roi = sampler_args
        S = simnibs_gpc.gPCSampler(mesh, poslist, fn_hdf5, roi)

        assert np.all(S.mesh_roi.elm.tag1 == 3)

        S.create_hdf5()
        with h5py.File(fn_hdf5) as f:
            assert f['mesh']
            assert f['mesh_roi']
            assert f['cond']
            assert np.all(f['roi'][()] == 3)

        S2 = simnibs_gpc.gPCSampler.load_hdf5(fn_hdf5)
        assert S2.mesh.nodes.nr == mesh.nodes.nr
        assert S2.mesh.elm.nr == mesh.elm.nr
        assert S2.roi == roi
        for i, c in enumerate(S2.poslist.cond):
            assert c.name == poslist.cond[i].name
            assert c.value == poslist.cond[i].value
            assert c.distribution_type == poslist.cond[i].distribution_type
            assert c.distribution_parameters == poslist.cond[i].distribution_parameters

    def test_record_data_matrix(self, sampler_args):
        mesh, poslist, fn_hdf5, roi = sampler_args
        S = simnibs_gpc.gPCSampler(mesh, poslist, fn_hdf5, roi)

        rand_vars1 = [.1, .2]
        potential1 = np.random.rand(100)
        E1 = np.random.rand(100, 3)
        S.record_data_matrix(potential1, 'potential', 'data')
        S.record_data_matrix(rand_vars1, 'random_vars', 'data')
        S.record_data_matrix(E1, 'E', 'data')

        rand_vars2 = [.2, .4]
        potential2 = np.random.rand(100)
        E2 = np.random.rand(100, 3)
        S.record_data_matrix(potential2, 'potential', 'data')
        S.record_data_matrix(rand_vars2, 'random_vars', 'data')
        S.record_data_matrix(E2, 'E', 'data')

        with h5py.File(fn_hdf5) as f:
            assert np.allclose(f['data/random_vars'][0, :], rand_vars1)
            assert np.allclose(f['data/random_vars'][1, :], rand_vars2)
            assert np.allclose(f['data/potential'][0, :], potential1)
            assert np.allclose(f['data/potential'][1, :], potential2)
            assert np.allclose(f['data/E'][0, ...], E1)
            assert np.allclose(f['data/E'][1, ...], E2)

    def test_calc_E(self, sampler_args):
        mesh, poslist, fn_hdf5, roi = sampler_args
        S = simnibs_gpc.gPCSampler(mesh, poslist, fn_hdf5, roi)
        v = mesh_io.NodeData(mesh.nodes.node_coord[:, 0], mesh=mesh)
        E = S._calc_E(v, None)
        assert np.allclose(E, [-1e3, 0, 0])

    def test_update_poslist(self, sampler_args):
        mesh, poslist, fn_hdf5, roi = sampler_args
        S = simnibs_gpc.gPCSampler(mesh, poslist, fn_hdf5, roi)
        random_vars = [0.35]
        poslist = S._update_poslist(random_vars)
        assert np.isclose(poslist.cond[2].value, 0.35)
        assert np.isclose(poslist.cond[3].value, 1)
        assert np.isclose(poslist.cond[4].value, 10)

    def test_run_N_random_simulations(self, sampler_args):
        mesh, poslist, fn_hdf5, roi = sampler_args
        S = simnibs_gpc.gPCSampler(mesh, poslist, fn_hdf5, roi)
        S.run_simulation = Mock()
        S.run_N_random_simulations(1000)
        for c in S.run_simulation.call_args_list:
            assert c[0][0][0] <= .4
            assert c[0][0][0] >= .3

    def test_tdcs_set_up(self, sampler_args):
        mesh, poslist, fn_hdf5, roi = sampler_args
        S = simnibs_gpc.TDCSgPCSampler(mesh, poslist, fn_hdf5, [1101, 1102], [-1, 1],
                                       roi)
        S.create_hdf5()
        with h5py.File(fn_hdf5) as f:
            assert f['mesh']
            assert f['mesh_roi']
            assert f['cond']
            assert np.all(f['roi'][()] == 3)
            assert np.all(f['el_tags'][()] == [1101, 1102])
            assert np.allclose(f['el_currents'][()], [-1, 1])

        S2 = simnibs_gpc.TDCSgPCSampler.load_hdf5(fn_hdf5)
        assert S2.mesh.nodes.nr == mesh.nodes.nr
        assert S2.mesh.elm.nr == mesh.elm.nr
        assert S2.roi == roi
        for i, c in enumerate(S2.poslist.cond):
            assert c.name == poslist.cond[i].name
            assert c.value == poslist.cond[i].value
            assert c.distribution_type == poslist.cond[i].distribution_type
            assert c.distribution_parameters == poslist.cond[i].distribution_parameters

    @patch.object(simnibs_gpc, 'fem')
    def test_tdcs_run(self, mock_fem, sampler_args):
        mesh, poslist, fn_hdf5, roi = sampler_args
        v = mesh.nodes.node_coord[:, 0]
        v_roi = mesh.crop_mesh(roi).nodes.node_coord[:, 0]

        mock_fem.tdcs.side_effect = [
            mesh_io.NodeData(v, mesh=mesh),
            mesh_io.NodeData(-v, mesh=mesh)]

        S = simnibs_gpc.TDCSgPCSampler(
            mesh, poslist, fn_hdf5, [1101, 1102], [-1, 1], roi)

        if extra_qoi:
            S.qoi_function = extra_qoi + S.qoi_function

        E1 = S.run_simulation([1])
        assert E1.shape == (3 * np.sum(mesh.elm.tag1 == 3), )
        assert np.allclose(E1.reshape(-1, 3), [-1e3, 0, 0])

        S.run_simulation([2])
        with h5py.File(fn_hdf5) as f:
            assert np.allclose(f['random_var_samples'][()], [[1], [2]])
            assert np.allclose(f['mesh_roi/data_matrices/v_samples'][0, :], v_roi)
            assert np.allclose(f['mesh_roi/data_matrices/v_samples'][1, :],-v_roi)
            assert np.allclose(f['mesh_roi/data_matrices/E_samples'][0, :],[-1e3, 0., 0.])
            assert np.allclose(f['mesh_roi/data_matrices/E_samples'][1, :],[1e3, 0., 0.])

    @patch.object(simnibs_gpc, 'fem')
    def test_tdcs_run(self, mock_fem, sampler_args):
        mesh, poslist, fn_hdf5, roi = sampler_args
        v = mesh.nodes.node_coord[:, 0]
        v_roi = mesh.crop_mesh(roi).nodes.node_coord[:, 0]

        mock_fem.tdcs.side_effect = [
            mesh_io.NodeData(v, mesh=mesh),
            mesh_io.NodeData(-v, mesh=mesh)]

        S = simnibs_gpc.TDCSgPCSampler(
            mesh, poslist, fn_hdf5, [1101, 1102], [-1, 1], roi)

        S.qoi_function['rand'] = lambda v, rand: rand

        E1 = S.run_simulation([1])
        assert E1.shape == (3 * np.sum(mesh.elm.tag1 == 3), )
        assert np.allclose(E1.reshape(-1, 3), [-1e3, 0, 0])

        S.run_simulation([2])
        with h5py.File(fn_hdf5) as f:
            assert np.allclose(f['random_var_samples'][()], [[1], [2]])
            assert np.allclose(f['mesh_roi/data_matrices/v_samples'][0, :], v_roi)
            assert np.allclose(f['mesh_roi/data_matrices/v_samples'][1, :],-v_roi)
            assert np.allclose(f['mesh_roi/data_matrices/E_samples'][0, :],[-1e3, 0., 0.])
            assert np.allclose(f['mesh_roi/data_matrices/E_samples'][1, :],[1e3, 0., 0.])
            assert np.allclose(f['mesh_roi/data_matrices/rand_samples'][:],[[1], [2]])


    def test_tms_set_up(self, sampler_args):
        mesh, poslist, fn_hdf5, roi = sampler_args
        matsimnibs = np.eye(4)
        didt = 1e5
        coil = 'coil.nii.gz'
        S = simnibs_gpc.TMSgPCSampler(mesh, poslist, fn_hdf5, coil, matsimnibs, didt, roi)
        S.create_hdf5()
        with h5py.File(fn_hdf5) as f:
            assert f['mesh']
            assert f['mesh_roi']
            assert f['cond']
            assert np.all(f['roi'][()] == 3)
            assert np.allclose(f['matsimnibs'], np.eye(4))
            assert np.allclose(f['didt'][()], 1e5)
            assert f['fnamecoil'][()] == coil.encode()

        S2 = simnibs_gpc.TMSgPCSampler.load_hdf5(fn_hdf5)
        assert S2.mesh.nodes.nr == mesh.nodes.nr
        assert S2.mesh.elm.nr == mesh.elm.nr
        assert S2.roi == roi
        assert np.allclose(S2.matsimnibs, S.matsimnibs)
        assert np.allclose(S2.didt, S.didt)
        assert S2.fnamecoil == S.fnamecoil
        for i, c in enumerate(S2.poslist.cond):
            assert c.name == poslist.cond[i].name
            assert c.value == poslist.cond[i].value
            assert c.distribution_type == poslist.cond[i].distribution_type
            assert c.distribution_parameters == poslist.cond[i].distribution_parameters



'''
class TestRunGPC:
    def test_prep_gpc(self):
        cond = [sim_struct.COND(), sim_struct.COND(), sim_struct.COND(), sim_struct.COND(), sim_struct.COND()]
        cond[0].distribution_type = 'uniform'
        cond[0].distribution_parameters = [0.2, 0.3]
        cond[2].distribution_type = 'beta'
        cond[2].distribution_parameters = [2, 3, 0.3, 0.4]
        cond[4].distribution_type = 'normal'
        cond[4].distribution_parameters = [0.2, 0.3]
        simlist = sim_struct.SimuList()
        simlist.cond = cond
        simlist.anisotropy_distribution_type = 'normal'
        simlist.anisotropy_distribution_parameters = [0, 1]
        random_vars, pdf_type, pdfshape, limits = simnibs_gpc.prep_gpc(simlist)
        assert np.all(random_vars == ['excentricity', 1, 3, 5])
        assert np.all(pdf_type == ['normal', 'beta', 'beta', 'normal'])
        assert np.all(pdfshape == [[0, 1, 2, 0.2], [1, 1, 3, 0.3]])
        assert np.all(limits == [[None, 0.2, 0.3, None], [None, 0.3, 0.4, None]])
 

    @patch.object(simnibs_gpc, 'pygpc')
    @patch.object(simnibs_gpc, 'gPC_regression')
    def test_run_gpc(self,mock_gpc_regression,  mock_pygpc):
        poslist = sim_struct.POSLIST()
        poslist.cond[0].distribution_type = 'uniform'
        poslist.cond[0].distribution_parameters = [0.2, 0.3]
        poslist.cond[1].distribution_type = 'beta'
        poslist.cond[1].distribution_parameters = [2, 3, 0.3, 0.4]
        poslist.cond[2].distribution_type = 'normal'
        poslist.cond[2].distribution_parameters = [0.2, 0.3]
        poslist.pos = [sim_struct.POSITION()]
        
        reg_mock = Mock(simnibs_gpc.pygpc.reg)
        reg_mock.poly_idx = None
        reg_mock.grid = Mock(simnibs_gpc.pygpc.grid)
        reg_mock.grid.coords_norm = None
        mock_pygpc.run_reg_adaptive.return_value = (reg_mock, None)

        simnibs_gpc.run_tms_gpc(poslist, 'test_msh.msh', 'test_simu', tissues=[1, 3], eps=1e-4)

        mock_pygpc.run_reg_adaptive.assert_called_with(
               ['beta', 'beta', 'normal'],
               [[1, 2, 0.2], [1, 3, 0.3]],
               [[0.2, 0.3, None], [0.3, 0.4, None]],
               poslist._run_gpc_simulation,
               args=([1, 2, 3], 'test_msh.msh', 'test_simu_0001', 0, [1, 3]),
               order_start=0,
               order_end=10,
               eps=1e-4,
               print_out=True) 
'''

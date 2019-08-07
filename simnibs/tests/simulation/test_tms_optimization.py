"""Testsuite to test the TMSOptimization class, including Neuronavigation I/O and some helper functions.

All filesystem I/O is done through tempfile package, that writes files to /tmp/* and clears things afterwards.
"""

import numpy as np
import os
import tempfile
import shutil
import pytest

import simnibs.simulation as simulation
from simnibs.utils.nnav import write_tms_navigator_im, simnibs2nnav

def sphere3_msh_fn():
    """Helper functions to get files and filenames."""
    return os.path.join(os.path.dirname(os.path.realpath(
        __file__)), '..', 'testing_files', 'sphere3.msh')


def example_nii_fn():
    """Helper functions to get files and filenames."""
    return os.path.join(os.path.dirname(os.path.realpath(
        __file__)), '..','..', 'resources', 'templates', 'MNI152_T1_1mm.nii.gz')


def coil_fn():
    """Helper functions to get files and filenames."""
    return os.path.join(os.path.dirname(os.path.realpath(
        __file__)), '..', '..', 'ccd-files', 'Magstim_70mm_Fig8.ccd')


class TestTMSOptimizationUtils:
    """Tests for helper functions for the TMS coil optimization"""

    def test_simnibs2nnav(self):
        """Tests simnibs2nnav for the different input types. No FSL processing is tested here."""

        with tempfile.TemporaryDirectory() as d:
            fn_nii = d + 'example.nii.gz'
            shutil.copy(example_nii_fn(), fn_nii)
            mat = np.array([[00., 01., 02., -10.],
                            [03., 04., 05., -11.],
                            [06., 07., 08., -12.],
                            [00., 00., 00., 01.]])
            im_mat_exp = np.array([[2., -1., 0., -11.],
                                  [5., -4., 3., -28.],
                                  [8., -7., 6., 6.],
                                  [0., 0., 0., 1.]])

            im = simnibs2nnav(fn_nii, fn_nii, mat)
            assert np.all(np.squeeze(im) == im_mat_exp), "simnibs2nnav() with np.ndarray() failed"

            pos = simulation.sim_struct.POSITION()
            pos.matsimnibs = mat
            im = simnibs2nnav(fn_nii, fn_nii, pos)
            assert np.all(np.squeeze(im) == im_mat_exp), "simnibs2nnav() with Position() failed"

            tmslist = simulation.sim_struct.TMSLIST()
            tmslist.add_position(pos)
            tmslist.add_position(pos)
            im = simnibs2nnav(fn_nii, fn_nii, tmslist)
            assert np.all(im[:, :, 0] == im_mat_exp), "simnibs2nnav() with TMSLIST() failed"
            assert np.all(im[:, :, 1] == im_mat_exp), "simnibs2nnav() with TMSLIST() failed"

            session = simulation.sim_struct.SESSION()
            session.add_poslist(tmslist)
            im = simnibs2nnav(fn_nii, fn_nii, session)
            assert np.all(im[:, :, 0] == im_mat_exp), "simnibs2nnav() with Session() failed"
            assert np.all(im[:, :, 1] == im_mat_exp), "simnibs2nnav() with Session() failed"

            tms_optim = simulation.TMSOPTIMIZATION()
            tms_optim.add_poslist(tmslist)
            im = simnibs2nnav(fn_nii, fn_nii, session)
            assert np.all(im[:, :, 0] == im_mat_exp), "simnibs2nnav() with TMSOPTIMIZATION() failed"
            assert np.all(im[:, :, 1] == im_mat_exp), "simnibs2nnav() with TMSOPTIMIZATION() failed"

    def test_write_write_tms_navigator_im(self):
        """Test the Localite instrument marker file generation"""

        im_mat = np.array([[00., 01., 02., -10.],
                           [03., 04., 05., -11.],
                           [06., 07., 08., -12.],
                           [09., 10., 11., -13.]])

        with tempfile.NamedTemporaryFile('w', suffix='.xml') as f:
            write_tms_navigator_im(im_mat, f.name, overwrite=True)

            cor = ['<?xml version="1.0" encoding="UTF-8"?>',
                   '<InstrumentMarkerList coordinateSpace="RAS">',
                   '    <!--This InstrumentMarkerList was written by SimNIBS-->',
                   '    <InstrumentMarker alwaysVisible="false" index="0" selected="false">',
                   '        <Marker additionalInformation="" color="#ff0000" description="opt_0" set="true">',
                   '            <Matrix4D',
                   '                data00="+0.00000000000000000" data01="+1.00000000000000000" '
                   'data02="+2.00000000000000000" data03="-10.00000000000000000"',
                   '                data10="+3.00000000000000000" data11="+4.00000000000000000" '
                   'data12="+5.00000000000000000" data13="-11.00000000000000000"',
                   '                data20="+6.00000000000000000" data21="+7.00000000000000000" '
                   'data22="+8.00000000000000000" data23="-12.00000000000000000"',
                   '                data30="+0.00000000000000000" data31="+0.00000000000000000" '
                   'data32="+0.00000000000000000" data33="+1.00000000000000000"/>',
                   '        </Marker>',
                   '    </InstrumentMarker>',
                   '</InstrumentMarkerList>']

            with open(f.name) as fp:
                line = fp.readline().rstrip()
                cnt = 0
                while line:
                    assert line == cor[cnt]
                    line = fp.readline().rstrip()
                    cnt += 1


class TestTMSOptimizationClass:
    """Tests for the simulation.sim_struct.TMSOPTIMIZATION class"""

    def test_mat_read(self):
        """Check I/O with read_mat() """
        tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
        tms_optim.add_poslist(simulation.sim_struct.TMSLIST())
        tms_optim.optimlist.add_position()
        tms_optim.optimlist.add_position()

        tms_optim.angle_limits = [-180, 180]  # list of num
        tms_optim.resolution_angle = 60.
        tms_optim.target = np.array((-30., 0.01, 49.))  # np.ndarray
        tms_optim.handle_direction_ref = np.array([-1, 0, 1.1])
        tms_optim.save_fields = ['normE', 'E']  # list of str
        tms_optim.optimlist.postprocess = ['E', 'e', 'J', 'j']  # list of str
        tms_optim.target_coil_matsim = np.array([[0, 1, 2],
                                                 [np.inf, -np.inf, np.nan],
                                                 [-1, -1., 0],
                                                 [1.1, 1.2, 1.3]])

        with tempfile.NamedTemporaryFile('w', suffix='.mat') as f:
            simulation.save_matlab_sim_struct(tms_optim, f.name)

            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.read_mat_struct(f.name)
            assert tms_optim.angle_limits == [-180, 180]  # list of num
            assert tms_optim.resolution_angle == 60.
            assert np.allclose(tms_optim.target, np.array((-30., 0.01, 49.)))  # np.ndarray

            assert tms_optim.save_fields == ['normE', 'E']  # list of str
            assert type(tms_optim.optimlist) == simulation.TMSLIST
            assert len(tms_optim.optimlist.pos) == 2
            assert tms_optim.optimlist.postprocess == ['E', 'e', 'J', 'j']  # list of str
            assert np.allclose(tms_optim.target_coil_matsim,
                               np.array([[0, 1, 2],
                                         [np.inf, -np.inf, np.nan],
                                         [-1, -1., 0],
                                         [1.1, 1.2, 1.3]]),
                               equal_nan=True)
            assert np.all(tms_optim.handle_direction_ref == np.array([-1, 0, 1.1]))

    def test_no_tmslist(self):
        """Check for AssertionError if not TMSLIST"""
        logger = simulation.logging.getLogger("simnibs")
        logger.handlers.clear()
        with tempfile.TemporaryDirectory() as pathfem:
            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.fnamehead = sphere3_msh_fn()
            tms_optim.pathfem = pathfem
            tms_optim._set_logger()
            with pytest.raises(AssertionError):
                tms_optim._prepare()

    def test_tmslist_no_poslist(self):
        """Check for error on no poslist"""
        logger = simulation.logging.getLogger("simnibs")
        logger.handlers.clear()
        with tempfile.TemporaryDirectory() as pathfem:
            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.fnamehead = sphere3_msh_fn()
            tms_optim.pathfem = pathfem
            tms_optim.optimlist = simulation.sim_struct.TMSLIST()
            tms_optim.optimlist.fnamecoil = coil_fn()
            tms_optim._set_logger()
            with pytest.raises(AssertionError):
                tms_optim._prepare()

    def test_tmslist_prepare(self):
        """Check TMSOPTIMIZATION._prepare() """
        logger = simulation.logging.getLogger("simnibs")
        logger.handlers.clear()
        with tempfile.TemporaryDirectory() as pathfem:
            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.fnamehead = sphere3_msh_fn()
            tms_optim.pathfem = pathfem
            tms_optim.optimlist = simulation.sim_struct.TMSLIST()
            tms_optim.optimlist.fnamecoil = coil_fn()
            tms_optim.optimlist.add_position()
            tms_optim.optim_name = 'test'
            tms_optim._set_logger()
            tms_optim._prepare()

            assert tms_optim.fname_hdf5 == os.path.join(pathfem, 'test.hdf5')
            assert tms_optim.n_sim == 1

    def test_tmslist_different_distances_in_optimlist(self):
        """Check for error on different coil distances in optimlist"""
        logger = simulation.logging.getLogger("simnibs")
        logger.handlers.clear()
        with tempfile.TemporaryDirectory() as pathfem:
            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.fnamehead = sphere3_msh_fn()
            tms_optim.pathfem = pathfem
            tms_optim.optimlist = simulation.sim_struct.TMSLIST()
            tms_optim.optimlist.fnamecoil = coil_fn()

            pos = simulation.sim_struct.POSITION()
            pos.pos_ydir = [0, 0, 0]
            pos.centre = [0, 0, 0]
            pos.distance = 1
            tms_optim.optimlist.add_position(pos)
            pos = simulation.sim_struct.POSITION()
            pos.pos_ydir = [0, 0, 0]
            pos.centre = [0, 0, 0]
            pos.distance = 2
            tms_optim.optimlist.add_position(pos)
            tms_optim._set_logger()
            with pytest.raises(AssertionError):
                tms_optim._prepare()


class TestTMSOptimization:
    """Tests for the optimization."""
    def test_optimization(self,n_cpu=1):
        """Test optimization and output files."""
        with tempfile.TemporaryDirectory() as pathfem:
            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.fnamehead = sphere3_msh_fn()
            tms_optim.pathfem = pathfem
            tms_optim.fnamecoil = coil_fn()
            tms_optim.optim_name = 'optimization'

            target = [-1, 1, 1]
            radius = 20
            resolution_pos = 20
            resolution_angle = 20
            angle_limits = [0, 15]
            handle_direction_ref = [-1, 2, .5]

            results = simulation.optim_tms.optimize_tms_coil_pos(tms_optim=tms_optim,
                                                                 target=target,
                                                                 n_cpu=n_cpu,
                                                                 radius=radius,
                                                                 resolution_pos=resolution_pos,
                                                                 resolution_angle=resolution_angle,
                                                                 handle_direction_ref=handle_direction_ref,
                                                                 angle_limits=angle_limits)

            res_matsimnibs = np.array([[-0.71645712, 0.53550334, 0.44713014, -41.88353921],
                                       [0.26607118, 0.80222911, -0.53444793, 52.22674481],
                                       [-0.64489947, -0.26394058, -0.71724476, 68.80491212],
                                       [0., 0., 0., 1.]])
            assert len(results['tms_optim'].optimlist.pos) == 1
            assert np.allclose(results['tms_optim'].optimlist.pos[0].matsimnibs, res_matsimnibs)
            assert results['tms_optim'].n_sim == 5
            assert os.path.exists(os.path.join(pathfem, 'optimization_' + results['tms_optim'].time_str + '.csv'))
            assert os.path.exists(os.path.join(pathfem, 'optimization_' + results['tms_optim'].time_str + '.mat'))
            assert os.path.exists(os.path.join(pathfem, 'optimization_' + results['tms_optim'].time_str + '.log'))
            assert os.path.exists(os.path.join(pathfem, 'optimization.hdf5'))

    def test_optimization_multicore(self):
        """Test the same optimization as above with 2 cpus in parallel"""
        self.test_optimization(2)

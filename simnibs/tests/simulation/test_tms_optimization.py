import pytest
import numpy as np
import os
import tempfile
from nose.tools import *
import logging

from simnibs.simulation import sim_struct, mesh_io, logging, TMSLIST, save_matlab_sim_struct
from simnibs.simulation.optim_tms import optimize_tms_coil_pos
import simnibs.simulation as simulation

"""Some helper functions to get files and filenames."""


def sphere3_msh():
    fn = sphere3_msh_fn()
    return mesh_io.read_msh(fn)


def sphere3_msh_fn():
    return os.path.join(os.path.dirname(os.path.realpath(
        __file__)), '..', 'testing_files', 'sphere3.msh')


def coil_fn():
    return os.path.join(os.path.dirname(os.path.realpath(
        __file__)), '..', '..', 'ccd-files', 'Magstim_70mm_Fig8.ccd')


class TestTMSOptimizationClass:
    """Tests for the sim_struct.TMSOPTIMIZATION class"""

    def test_mat_read(self):
        """Check additions to read_mat() """
        tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
        tms_optim.add_poslist(sim_struct.TMSLIST())
        tms_optim.optimlist.add_position()
        tms_optim.optimlist.add_position()
        tms_optim.angle_limits = [-180, 180]  # list of num
        tms_optim.resolution_angle = 60.
        tms_optim.target = np.array((-30., 0.01, 49.))  # np.ndarray
        tms_optim.save_fields = ['normE', 'E']  # list of str
        tms_optim.optimlist.postprocess = ['E', 'e', 'J', 'j']  # list of str
        with tempfile.NamedTemporaryFile('w', suffix='.mat') as f:
            save_matlab_sim_struct(tms_optim, f.name)

            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.read_mat_struct(f.name)
            assert tms_optim.angle_limits == [-180, 180]  # list of num
            assert tms_optim.resolution_angle == 60.
            assert np.allclose(tms_optim.target, np.array((-30., 0.01, 49.)))  # np.ndarray

            assert tms_optim.save_fields == ['normE', 'E']  # list of str
            assert type(tms_optim.optimlist) == TMSLIST
            assert len(tms_optim.optimlist.pos) == 2
            assert tms_optim.optimlist.postprocess == ['E', 'e', 'J', 'j']  # list of str

    @raises(AssertionError)
    def test_no_tmslist(self):
        """Check for assertionError if not TMSLIST"""
        logger = logging.getLogger("simnibs")
        logger.handlers.clear()
        with tempfile.TemporaryDirectory() as pathfem:
            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.fnamehead = sphere3_msh_fn()
            tms_optim.pathfem = pathfem
            tms_optim._set_logger()
            tms_optim._prepare()

    @raises(AssertionError)
    def test_tmslist_no_poslist(self):
        """Check for error on no poslist"""
        logger = logging.getLogger("simnibs")
        logger.handlers.clear()
        with tempfile.TemporaryDirectory() as pathfem:
            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.fnamehead = sphere3_msh_fn()
            tms_optim.pathfem = pathfem
            tms_optim.optimlist = sim_struct.TMSLIST()
            tms_optim.optimlist.fnamecoil = coil_fn()
            tms_optim._set_logger()
            tms_optim._prepare()

    def test_tmslist_prepare(self):
        """Check TMSOPTIMIZATION._prepare() """
        logger = logging.getLogger("simnibs")
        logger.handlers.clear()
        with tempfile.TemporaryDirectory() as pathfem:
            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.fnamehead = sphere3_msh_fn()
            tms_optim.pathfem = pathfem
            tms_optim.optimlist = sim_struct.TMSLIST()
            tms_optim.optimlist.fnamecoil = coil_fn()
            tms_optim.optimlist.add_position()
            tms_optim.optim_name = 'test'
            tms_optim._set_logger()
            tms_optim._prepare()

            assert tms_optim.fname_hdf5 == os.path.join(pathfem, 'test.hdf5')
            assert tms_optim.n_sim == 1

    @raises(AssertionError)
    def test_tmslist_different_distances_in_optimlist(self):
        """Check for error on different coil distances in optimlist"""
        logger = logging.getLogger("simnibs")
        logger.handlers.clear()
        with tempfile.TemporaryDirectory() as pathfem:
            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.fnamehead = sphere3_msh_fn()
            tms_optim.pathfem = pathfem
            tms_optim.optimlist = sim_struct.TMSLIST()
            tms_optim.optimlist.fnamecoil = coil_fn()

            pos = sim_struct.POSITION()
            pos.pos_ydir = [0, 0, 0]
            pos.centre = [0, 0, 0]
            pos.distance = 1
            tms_optim.optimlist.add_position(pos)
            pos = sim_struct.POSITION()
            pos.pos_ydir = [0, 0, 0]
            pos.centre = [0, 0, 0]
            pos.distance = 2
            tms_optim.optimlist.add_position(pos)
            tms_optim._set_logger()
            tms_optim._prepare()


class TestTMSOptimization:
    """Tests for the optimization."""
    @raises(AssertionError)
    def test_optimization_empty_matsimnibs(self):
        """Test for error on empty position"""
        logger = logging.getLogger("simnibs")
        logger.handlers.clear()
        with tempfile.TemporaryDirectory() as pathfem:
            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.fnamehead = sphere3_msh_fn()
            tms_optim.pathfem = pathfem
            tms_optim.optimlist = sim_struct.TMSLIST()
            tms_optim.optimlist.fnamecoil = coil_fn()

            target = [-1, 1, 1]
            # pos = sim_struct.POSITION()
            # pos.pos_ydir = [-1, 0, 0]
            # pos.centre = [0, 0, 0]
            # pos.distance = 1
            # pos.calc_matsimnibs(tms_optim.fnamehead)
            tms_optim.optimlist.add_position()  # this should not work

            resolution_pos = 10
            resolution_angle = 15
            handle_direction_ref = [-1, 2, .5]
            optimize_tms_coil_pos(tms_optim=tms_optim,
                                  target=target,
                                  n_cpu=1,
                                  resolution_pos=resolution_pos,
                                  resolution_angle=resolution_angle,
                                  handle_direction_ref=handle_direction_ref)

    def test_optimization(self,n_cpu=1):
        """Test optimization and output files."""
        with tempfile.TemporaryDirectory() as pathfem:
            tms_optim = simulation.optim_tms.TMSOPTIMIZATION()
            tms_optim.fnamehead = sphere3_msh_fn()
            tms_optim.pathfem = pathfem
            tms_optim.optimlist = sim_struct.TMSLIST()
            tms_optim.optimlist.fnamecoil = coil_fn()
            tms_optim.optim_name = 'optimization'

            target = [-1, 1, 1]
            radius = 20
            resolution_pos = 20
            resolution_angle = 20
            angle_limits = [0, 15]
            handle_direction_ref = [-1, 2, .5]

            results = optimize_tms_coil_pos(tms_optim=tms_optim,
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
        # pass
import os
import pytest
import simnibs
import numpy as np
from ... mesh_tools import mesh_io
from ... import SIMNIBSDIR
from .. region_of_interest import RegionOfInterest


class TestRegionOfInterest:
    def test_RegionOfInterestInitializer_custom_center(self):
        nodes = np.array([[-1.2, 1.4, 7.1],
                          [-.9, 1.4, 7.2],
                          [-1.0, 1.3, 7.1],
                          [-.6, 1.3, 7.2],
                          [-.7, 1.2, 7.1]])
        con = np.array([[0, 1, 2],
                        [2, 3, 1],
                        [4, 3, 2]])
        center = np.mean(nodes[con, ], axis=1)
        fn_mesh = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
        mesh = mesh_io.read_msh(fn_mesh)
        roi = RegionOfInterest(center=center,
                               nodes=nodes,
                               con=con,
                               domains=None,
                               mesh=mesh)

    def test_RegionOfInterestInitializer_custom_domains(self):
        domains = [3, 4]
        fn_mesh = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
        mesh = mesh_io.read_msh(fn_mesh)
        roi = RegionOfInterest(center=None,
                               nodes=None,
                               con=None,
                               domains=domains,
                               mesh=mesh)

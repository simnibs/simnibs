import os
from pathlib import Path
import numpy as np
import pytest

from simnibs.simulation.coil.coil import Coil
from simnibs.simulation.coil.coil_element import CoilDipoles

from .... import SIMNIBSDIR

@pytest.fixture(scope='module')
def testcoil_ccd():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'testcoil.ccd')
    return fn

class TestReadCoil:
    def test_read_simple_ccd(self, tmp_path: Path):
        with open(tmp_path / 'test.ccd', 'w') as f:
            f.write('# number of elements\n')
            f.write('3\n')
            f.write('1e-003 3e-003 -4e-003 0e+000 0e+000 -4e-006\n')
            f.write('3e-003 1e-003  1e-003 1e+000 0e+000 0e-000\n')
            f.write('4e-003 2e-003 -5e-003 0e+000 1e+000 0e-000\n')
        coil = Coil.from_ccd(str(tmp_path / 'test.ccd'))
        dipole_elements : CoilDipoles = coil.coil_elements[0]
        assert np.allclose(dipole_elements.points, np.array([[1e-3, 3e-3, -4e-3],
                                               [3e-3, 1e-3, 1e-3],
                                               [4e-3, 2e-3, -5e-3]]))
        assert np.allclose(dipole_elements.values, np.array([[0, 0, -4e-6],
                                             [1, 0, 0],
                                             [0, 1, 0]], dtype=float))

    def test_read_ccd(self, testcoil_ccd: str):
        coil = Coil.from_ccd(testcoil_ccd)
        dipole_elements : CoilDipoles = coil.coil_elements[0]
        assert np.allclose(coil.limits,((-100,100),(-100,100),(-100,100)))
        assert np.allclose(coil.resolution,(10,10,10))
        assert coil.name =='Test coil'
        assert dipole_elements.stimulator.max_di_dt==100.0
        #assert float(info['dIdtstim'])==100.0
        assert coil.brand=='None'
        assert dipole_elements.stimulator.name=='None'
        assert np.allclose(dipole_elements.points, ((-1e-2,0,-1e-3),
                                        (1e-2,0,-1e-3),
                                        (-2e-2,2e-2,-2e-3),
                                        (2e-2,-2e-2,-2e-3)))
        assert np.allclose(dipole_elements.values, ((1e-6,0,2e-6),
                                      (-1e-6,0,-2e-6),
                                      (0,1e-6,0),
                                      (0,-1e-6,0)))
import json
import os
from pathlib import Path
from typing import Any
import numpy as np
import nibabel as nib
import pytest
from simnibs.mesh_tools.mesh_io import Elements, Msh, Nodes

from simnibs.simulation.tms_coil.tms_coil import TmsCoil
from simnibs.simulation.tms_coil.tms_coil_deformation import TmsCoilRotation
from simnibs.simulation.tms_coil.tms_coil_element import DipoleElements, LinePointElements, SampledGridPointElements
from simnibs.simulation.tms_coil.tms_coil_model import TmsCoilModel

from .... import SIMNIBSDIR

@pytest.fixture(scope='module')
def testcoil_ccd():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'testcoil.ccd')
    return fn

class TestReadCoil:
    def test_read_minimal_tcd(self, minimal_tcd_coil_dict: dict[str, Any], tmp_path: Path):
        with open(tmp_path / 'test.tcd', 'w') as f:
            json.dump(minimal_tcd_coil_dict, f)
        coil = TmsCoil.from_tcd(str(tmp_path / 'test.tcd'))
        assert coil.name is None
        assert coil.limits is None
        assert coil.casing is None
        assert len(coil.elements) == 1
        assert isinstance(coil.elements[0], LinePointElements)
        np.testing.assert_array_almost_equal(coil.elements[0].points, [[1, 2, 3]])
        assert len(coil.elements[0].deformations) == 0
        assert coil.elements[0].casing is None

    def test_read_simplest_ccd(self, tmp_path: Path):
        with open(tmp_path / 'test.ccd', 'w') as f:
            f.write('# number of elements\n')
            f.write('1\n')
            f.write('1e-003 3e-003 -4e-003 0e+000 0e+000 -4e-006\n')
        coil = TmsCoil.from_ccd(str(tmp_path / 'test.ccd'))
        dipole_elements : DipoleElements = coil.elements[0]
        assert np.allclose(dipole_elements.points, np.array([[1e-3, 3e-3, -4e-3]]) * 1e3)
        assert np.allclose(dipole_elements.values, np.array([[0, 0, -4e-6]], dtype=float))

    def test_read_simple_ccd(self, tmp_path: Path):
        with open(tmp_path / 'test.ccd', 'w') as f:
            f.write('# number of elements\n')
            f.write('3\n')
            f.write('1e-003 3e-003 -4e-003 0e+000 0e+000 -4e-006\n')
            f.write('3e-003 1e-003  1e-003 1e+000 0e+000 0e-000\n')
            f.write('4e-003 2e-003 -5e-003 0e+000 1e+000 0e-000\n')
        coil = TmsCoil.from_ccd(str(tmp_path / 'test.ccd'))
        dipole_elements : DipoleElements = coil.elements[0]
        assert np.allclose(dipole_elements.points, np.array([[1e-3, 3e-3, -4e-3],
                                               [3e-3, 1e-3, 1e-3],
                                               [4e-3, 2e-3, -5e-3]]) * 1e3)
        assert np.allclose(dipole_elements.values, np.array([[0, 0, -4e-6],
                                             [1, 0, 0],
                                             [0, 1, 0]], dtype=float))

    def test_read_ccd(self, testcoil_ccd: str):
        coil = TmsCoil.from_ccd(testcoil_ccd)
        dipole_elements : DipoleElements = coil.elements[0]
        assert coil.limits is not None
        assert coil.resolution is not None
        assert np.allclose(coil.limits,((-100,100),(-100,100),(-100,100)))
        assert np.allclose(coil.resolution,(10,10,10))
        assert coil.name =='Test coil'
        assert dipole_elements.stimulator is not None
        assert dipole_elements.stimulator.max_di_dt==100.0
        assert coil.brand=='None'
        assert dipole_elements.stimulator.name=='None'
        assert np.allclose(dipole_elements.points, np.array([[-1e-2,0,-1e-3],
                                        [1e-2,0,-1e-3],
                                        [-2e-2,2e-2,-2e-3],
                                        [2e-2,-2e-2,-2e-3]]) * 1e3)
        assert np.allclose(dipole_elements.values, ((1e-6,0,2e-6),
                                      (-1e-6,0,-2e-6),
                                      (0,1e-6,0),
                                      (0,-1e-6,0)))
        

class TestWriteCoil:
    def test_write_minimal_tcd(self, minimal_tcd_coil: TmsCoil, minimal_tcd_coil_dict: dict[str, Any], tmp_path: Path):
        minimal_tcd_coil.write(str(tmp_path / 'test.tcd'))

        with open(tmp_path / 'test.tcd', "r") as fid:
            coil = json.loads(fid.read())

        assert str(coil) == str(minimal_tcd_coil_dict)

    def test_write_medium_tcd(self, medium_tcd_coil: TmsCoil, medium_tcd_coil_dict: dict[str, Any], tmp_path: Path):
        medium_tcd_coil.write(str(tmp_path / 'test.tcd'))

        with open(tmp_path / 'test.tcd', "r") as fid:
            coil = json.loads(fid.read())
            
        assert str(coil) == str(medium_tcd_coil_dict)

    def test_write_full_tcd(self, full_tcd_coil: TmsCoil, full_tcd_coil_dict: dict[str, Any], tmp_path: Path):
        full_tcd_coil.write(str(tmp_path / 'test.tcd'))

        with open(tmp_path / 'test.tcd', "r") as fid:
            coil = json.loads(fid.read())

        assert str(coil) == str(full_tcd_coil_dict)

    def test_write_nifti(self, tmp_path: Path):
        affine = np.array(
            [
                [5.0, 0.0, 0.0, -300],
                [0.0, 5.0, 0.0, -200],
                [0.0, 0.0, 5.0, 0.0],
                [0.0, 0.0, 0.0, 1],
            ]
        )

        field = np.ones((31, 21, 11, 3))
        field[..., 1] = 2
        field[..., 2] = 3

        image = nib.Nifti1Image(field, affine)
        nib.save(image, tmp_path / 'test.nii')

        coil = TmsCoil.from_nifti(str(tmp_path / 'test.nii'))
        coil.write(str(tmp_path / 'test.tcd'))

        coil = TmsCoil.from_tcd(str(tmp_path / 'test.tcd'))

        sampled_grid : SampledGridPointElements = coil.elements[0]

        np.testing.assert_allclose(sampled_grid.data, field)
        np.testing.assert_allclose(sampled_grid.affine, affine)


class TestPositionOptimization:
    def test_simple_optimization(self, sphere3_msh: Msh):
        skin_surface = sphere3_msh.crop_mesh(tags=[1005])
        coil_affine = np.array([[1,0,0,0], [0,1,0,0], [0,0,1,100.1], [0,0,0,1]])
        casings = [
            TmsCoilModel(Msh(Nodes(np.array([[-20, -20, 0], [-20, 20, 0], [20, -20, 0], [20, 20, 0]])), Elements(triangles=np.array([[0,1,2], [3, 2, 1]]) + 1)), None, None),
            TmsCoilModel(Msh(Nodes(np.array([[-20, -60, 0], [-20, -20, 0], [20, -60, 0], [20, -20, 0]])), Elements(triangles=np.array([[0,1,2], [3, 2, 1]]) + 1)), None, None),
            TmsCoilModel(Msh(Nodes(np.array([[-20, 20, 0], [-20, 60, 0], [20, 20, 0], [20, 60, 0]])), Elements(triangles=np.array([[0,1,2], [3, 2, 1]]) + 1)), None, None),
        ]
        deformations = [
            TmsCoilRotation(0, (0, 90), np.array([0, -20, 0]), np.array([40, -20, 0])),
            TmsCoilRotation(0, (0, 90), np.array([40, 20, 0]), np.array([0, 20, 0]),) 
        ]
        coil = TmsCoil(None, None, None, None, None, None,
                    [LinePointElements(None, casings[0], None, np.array([[1, 2, 3]]), None, None),
                    LinePointElements(None, casings[1], [deformations[0]], np.array([[1, 2, 3]]), None, None),
                    LinePointElements(None, casings[2], [deformations[1]], np.array([[1, 2, 3]]), None, None),
                    ])

        assert coil.optimize_deformations(skin_surface, coil_affine)[1] < 6.2
 
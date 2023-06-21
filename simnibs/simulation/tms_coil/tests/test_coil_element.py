import numpy as np

from simnibs.mesh_tools.mesh_io import Msh
from simnibs.simulation.tms_coil.tms_coil_constants import TmsCoilElementTag
from simnibs.simulation.tms_coil.tms_coil_element import (
    DipoleElements,
    LineSegmentElements,
    SampledGridPointElements,
)
from simnibs.simulation.tms_coil.tms_stimulator import TmsStimulator

class TestCalcdAdt:
    def test_calc_dAdt_dipoles(self):
        dipoles = DipoleElements(
            None,
            None,
            None,
            np.array([[0, 0, 0], [0, 0, 0]]),
            np.array([[0, 0, 1], [0, 1, 0]]),
            TmsStimulator(None, None, None, None),
        )

        dipoles.stimulator.di_dt = 1e6

        target_pos = np.array([[1e3, 0.0, 0.0], [0.0, 1e3, 0.0], [0, 0.0, 1e3]])
        coil_matrix = np.eye(4)
        da_dt = dipoles.get_da_dt(target_pos, coil_matrix)
        np.testing.assert_allclose(
            da_dt, np.array([[0.0, 0.1, -0.1], [-0.1, 0, 0], [0.1, 0.0, 0.0]])
        )

    def test_calc_dAdt_triangle_dipoles(self):
        dipoles = DipoleElements(
            None,
            None,
            None,
            np.array([[0, 0, 0], [0, 0, 1], [0, 1 , 0.5]]),
            np.array([[0, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0.5, 0]]),
            TmsStimulator(None, None, None, None),
        )
        dipoles.stimulator.di_dt = 1

        target_pos = np.array([[0, 0.5, 0.5]])
        coil_matrix = np.eye(4)

        A = np.zeros_like(target_pos, dtype=float)
        for p, m in zip(dipoles.points, dipoles.values):
            r = (target_pos - p) * 1e-3
            A += 1e-7 * np.cross(m, r) / (np.linalg.norm(r, axis=1)[:, None] ** 3)

        np.testing.assert_allclose(
            dipoles.get_da_dt(target_pos, coil_matrix), A
        )

    def test_calc_dAdt_single_line_segment(self):
        seg_pos = np.array([[0, 0, 0]])
        seg_dir = np.array([[1, 0, 0]])
        line_segments = LineSegmentElements(
            None, None, None, seg_pos, seg_dir, TmsStimulator(None, None, None, None)
        )
        line_segments.stimulator.di_dt = 1

        target_pos = np.array([[-1, 1, 0]])
        coil_matrix = np.eye(4)
        
        r1 = np.sqrt(np.sum(((target_pos - seg_pos)**2)))
        A = seg_dir / r1[None] * 1e-7

        np.testing.assert_allclose(
            line_segments.get_da_dt(target_pos, coil_matrix), A
        )

    def test_calc_dAdt_triangle_line_segments(self):
        seg_pos = np.array([[0, 0, 0], [0, 0, 1], [0, 1 , 0.5]])
        line_segments = LineSegmentElements(
            None, None, None, seg_pos, None, TmsStimulator(None, None, None, None)
        )
        line_segments.stimulator.di_dt = 1

        target_pos = np.array([[0, 0.5, 0.5]])
        coil_matrix = np.eye(4)

        r1 = np.sqrt(np.sum((target_pos.T[:,:,None]-seg_pos.T[:,None,:])**2,axis=0))
        A = np.sum(line_segments.values.T[:,None,:] / r1[None],axis=2).T * 1e-7

        np.testing.assert_allclose(
            line_segments.get_da_dt(target_pos, coil_matrix), A
        )

    def test_calc_dAdt_sampled_elements(self, sphere3_msh: Msh):
        affine = np.array(
            [
                [5.0, 0.0, 0.0, -300],
                [0.0, 5.0, 0.0, -200],
                [0.0, 0.0, 5.0, 0.0],
                [0.0, 0.0, 0.0, 1],
            ]
        )

        field = np.ones((121, 81, 41, 3))
        field[..., 1] = 2
        field[..., 2] = 3
        coil_matrix = np.array(
            [
                [0.0, 1.0, 0.0, 0],
                [1.0, 0.0, 0.0, 0],
                [0.0, 0.0, 1.0, -100.0],
                [0.0, 0.0, 0.0, 1],
            ]
        )

        sampled_elements = SampledGridPointElements(
            None, None, None, field, affine, TmsStimulator(None, None, None, None)
        )
        sampled_elements.stimulator.di_dt = 1e6
        da_dt = sampled_elements.get_da_dt(
            sphere3_msh.nodes.node_coord, coil_matrix
        )
        np.testing.assert_allclose(da_dt[:, 0], 2e6, atol=1e-6)
        np.testing.assert_allclose(da_dt[:, 1], 1e6, atol=1e-6)
        np.testing.assert_allclose(da_dt[:, 2], 3e6, atol=1e-6)


class TestGetMesh:
    def test_dipole_element_mesh(self):
        dipoles = DipoleElements(
            None,
            None,
            None,
            np.array([[0, 0, 0], [0, 0, 0]]),
            np.array([[0, 0, 1], [0, 1, 0]]),
            None,
        )
        mesh = dipoles.get_mesh(np.eye(4), False, False, False, True)

        assert mesh.elm.node_number_list.shape[0] == 2
        np.testing.assert_allclose(np.unique(mesh.elm.elm_type), [15])
        assert np.all(mesh.elm.tag1 < 100)
        assert np.all(mesh.elm.tag1 == TmsCoilElementTag.DIPOLES)
        np.testing.assert_allclose(
            mesh.nodes.node_coord, [[0, 0, 0], [0, 0, 0], [0, 0, 1], [0, 1, 0]]
        )

    def test_line_element_element_mesh(self):
        dipoles = LineSegmentElements(
            None,
            None,
            None,
            np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]]),
            np.array([[1, 0, 0], [0, 1, 0], [-1, -1, 0]]),
            None,
        )
        mesh = dipoles.get_mesh(np.eye(4), False, False, False, True, 100)

        assert mesh.elm.node_number_list.shape[0] == 3
        np.testing.assert_allclose(np.unique(mesh.elm.elm_type), [1])
        assert np.all(mesh.elm.tag1 < 200)
        assert np.all(mesh.elm.tag1 > 99)
        assert np.all(mesh.elm.tag1 == 100 + TmsCoilElementTag.LINE_ELEMENTS)
        np.testing.assert_allclose(
            mesh.nodes.node_coord,
            [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 0, 0], [1, 1, 0], [0, 0, 0]],
        )

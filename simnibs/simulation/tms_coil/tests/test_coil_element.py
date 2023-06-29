import numpy as np
import pytest

from simnibs.mesh_tools.mesh_io import Msh
from simnibs.simulation.tms_coil.tms_coil_constants import TmsCoilElementTag
from simnibs.simulation.tms_coil.tms_coil_deformation import TmsCoilRotation, TmsCoilTranslation
from simnibs.simulation.tms_coil.tms_coil_element import (
    DipoleElements,
    LineSegmentElements,
    SampledGridPointElements,
)
from simnibs.simulation.tms_coil.tms_stimulator import TmsStimulator

class TestCalcdAdt:
    def test_calc_dAdt_dipoles(self):
        dipoles = DipoleElements(
            TmsStimulator(None),
            np.array([[0, 0, 0], [0, 0, 0]]),
            np.array([[0, 0, 1], [0, 1, 0]]),
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
            TmsStimulator(None),
            np.array([[0, 0, 0], [0, 0, 1], [0, 1 , 0.5]]),
            np.array([[0, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0.5, 0]]),
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
            TmsStimulator(None), seg_pos, seg_dir
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
            TmsStimulator(None), seg_pos
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
            TmsStimulator(None), field, affine, 
        )
        sampled_elements.stimulator.di_dt = 1e6
        da_dt = sampled_elements.get_da_dt(
            sphere3_msh.nodes.node_coord, coil_matrix
        )
        np.testing.assert_allclose(da_dt[:, 0], 2e6, atol=1e-6)
        np.testing.assert_allclose(da_dt[:, 1], 1e6, atol=1e-6)
        np.testing.assert_allclose(da_dt[:, 2], 3e6, atol=1e-6)

class TestCombinedTransformation:

    @pytest.mark.parametrize("rotation_amount", [-121.99, -24, 0, 44, 180.23, 720])
    @pytest.mark.parametrize("translation_amount", [-199.9, -34, 0, 21, 122])
    @pytest.mark.parametrize("rot_axis", [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])])
    @pytest.mark.parametrize("trans_axis", [0, 1, 2])
    def test_rotation_translation(self, rotation_amount, translation_amount, rot_axis, trans_axis):
        dipoles_rot_trans = DipoleElements(
            None, [[1,2,3]], [[1,2,3]],
            deformations=[TmsCoilRotation(rotation_amount, (-5000, 5000), np.array([0, 0, 0]), rot_axis), TmsCoilTranslation(translation_amount, (-200, 200), trans_axis)],
        )

        dipoles_trans_rot = DipoleElements(
            None,
            [[1,2,3]],
            [[1,2,3]],
            deformations=[TmsCoilTranslation(translation_amount, (-200, 200), trans_axis), TmsCoilRotation(rotation_amount, (-5000, 5000), np.array([0, 0, 0]), rot_axis)],
        )

        if rot_axis[0] == 1:
            rot = np.array([
                [1, 0, 0, 0],
                [0, np.cos(np.radians(rotation_amount)), -np.sin(np.radians(rotation_amount)), 0],
                [0, np.sin(np.radians(rotation_amount)), np.cos(np.radians(rotation_amount)), 0],
                [0, 0, 0, 1]
            ])
        elif rot_axis[1] == 1:
             rot = np.array([
                [np.cos(np.radians(rotation_amount)), 0, np.sin(np.radians(rotation_amount)), 0],
                [0, 1, 0, 0],
                [-np.sin(np.radians(rotation_amount)), 0, np.cos(np.radians(rotation_amount)), 0],
                [0, 0, 0, 1]
            ])
        else:
             rot = np.array([
                [np.cos(np.radians(rotation_amount)), -np.sin(np.radians(rotation_amount)), 0, 0],
                [np.sin(np.radians(rotation_amount)), np.cos(np.radians(rotation_amount)), 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
            ])
             
        trans = np.eye(4)
        trans[trans_axis, 3] = translation_amount

        np.testing.assert_allclose(dipoles_rot_trans.get_combined_transformation(), trans @ rot, atol=1e-5)
        np.testing.assert_allclose(dipoles_trans_rot.get_combined_transformation(), rot @ trans, atol=1e-5)

    @pytest.mark.parametrize("rotation_amount", [-121.99, 0, 720])
    @pytest.mark.parametrize("translation_amount1", [-34, 0, 122])
    @pytest.mark.parametrize("translation_amount2", [-3.1, 0, 192.64])
    @pytest.mark.parametrize("rot_axis", [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])])
    @pytest.mark.parametrize("trans_axis1", [0, 1, 2])
    @pytest.mark.parametrize("trans_axis2", [0, 1, 2])
    def test_rotation_translation_translation(self, rotation_amount, translation_amount1, translation_amount2, rot_axis, trans_axis1, trans_axis2):
        dipoles_rot_trans_trans = DipoleElements(
            None,
            [[1,2,3]],
            [[1,2,3]],
            deformations=[TmsCoilRotation(rotation_amount, (-5000, 5000), np.array([0, 0, 0]), rot_axis),
             TmsCoilTranslation(translation_amount1, (-200, 200), trans_axis1),
             TmsCoilTranslation(translation_amount2, (-200, 200), trans_axis2)],
        )

        dipoles_trans_trans_rot = DipoleElements(
            None,
            [[1,2,3]],
            [[1,2,3]],
            deformations=[TmsCoilTranslation(translation_amount1, (-200, 200), trans_axis1),
             TmsCoilTranslation(translation_amount2, (-200, 200), trans_axis2),
             TmsCoilRotation(rotation_amount, (-5000, 5000), np.array([0, 0, 0]), rot_axis)],            
        )

        dipoles_trans_rot_trans = DipoleElements(
            None,
            [[1,2,3]],
            [[1,2,3]],
            deformations=[TmsCoilTranslation(translation_amount1, (-200, 200), trans_axis1),
             TmsCoilRotation(rotation_amount, (-5000, 5000), np.array([0, 0, 0]), rot_axis),
             TmsCoilTranslation(translation_amount2, (-200, 200), trans_axis2)],
        )

        if rot_axis[0] == 1:
            rot = np.array([
                [1, 0, 0, 0],
                [0, np.cos(np.radians(rotation_amount)), -np.sin(np.radians(rotation_amount)), 0],
                [0, np.sin(np.radians(rotation_amount)), np.cos(np.radians(rotation_amount)), 0],
                [0, 0, 0, 1]
            ])
        elif rot_axis[1] == 1:
             rot = np.array([
                [np.cos(np.radians(rotation_amount)), 0, np.sin(np.radians(rotation_amount)), 0],
                [0, 1, 0, 0],
                [-np.sin(np.radians(rotation_amount)), 0, np.cos(np.radians(rotation_amount)), 0],
                [0, 0, 0, 1]
            ])
        else:
             rot = np.array([
                [np.cos(np.radians(rotation_amount)), -np.sin(np.radians(rotation_amount)), 0, 0],
                [np.sin(np.radians(rotation_amount)), np.cos(np.radians(rotation_amount)), 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
            ])
             
        trans1 = np.eye(4)
        trans1[trans_axis1, 3] = translation_amount1

        trans2 = np.eye(4)
        trans2[trans_axis2, 3] = translation_amount2

        np.testing.assert_allclose(dipoles_rot_trans_trans.get_combined_transformation(), trans2 @ trans1 @ rot, atol=1e-5)
        np.testing.assert_allclose(dipoles_trans_trans_rot.get_combined_transformation(), rot @ trans2 @ trans1, atol=1e-5)
        np.testing.assert_allclose(dipoles_trans_rot_trans.get_combined_transformation(), trans2 @ rot @ trans1 , atol=1e-5)


    @pytest.mark.parametrize("rotation_amount", [-121.99, 0, 720])
    @pytest.mark.parametrize("translation_amount1", [-34, 0, 122])
    @pytest.mark.parametrize("rotation_amount2", [-3.1, 0, 192.64])
    @pytest.mark.parametrize("rot_axis", [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])])
    @pytest.mark.parametrize("trans_axis1", [0, 1, 2])
    @pytest.mark.parametrize("rot_axis2", [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])])
    def test_rotation_translation_translation(self, rotation_amount, translation_amount1, rotation_amount2, rot_axis, trans_axis1, rot_axis2):
        dipoles_rot_trans_rot = DipoleElements(
            None,
            [[1,2,3]],
            [[1,2,3]],
            deformations=[TmsCoilRotation(rotation_amount, (-5000, 5000), np.array([0, 0, 0]), rot_axis),
             TmsCoilTranslation(translation_amount1, (-200, 200), trans_axis1),
             TmsCoilRotation(rotation_amount2, (-5000, 5000), np.array([0, 0, 0]), rot_axis2)],
        )

        dipoles_trans_rot_rot = DipoleElements(
            None,
            [[1,2,3]],
            [[1,2,3]],
            deformations=[TmsCoilTranslation(translation_amount1, (-200, 200), trans_axis1),
             TmsCoilRotation(rotation_amount, (-5000, 5000), np.array([0, 0, 0]), rot_axis),
             TmsCoilRotation(rotation_amount2, (-5000, 5000), np.array([0, 0, 0]), rot_axis2)],            
        )

        dipoles_rot_rot_trans = DipoleElements(
            None,
            [[1,2,3]],
            [[1,2,3]],
            deformations=[TmsCoilRotation(rotation_amount, (-5000, 5000), np.array([0, 0, 0]), rot_axis),
             TmsCoilRotation(rotation_amount2, (-5000, 5000), np.array([0, 0, 0]), rot_axis2),
             TmsCoilTranslation(translation_amount1, (-200, 200), trans_axis1)],
        )

        if rot_axis[0] == 1:
            rot = np.array([
                [1, 0, 0, 0],
                [0, np.cos(np.radians(rotation_amount)), -np.sin(np.radians(rotation_amount)), 0],
                [0, np.sin(np.radians(rotation_amount)), np.cos(np.radians(rotation_amount)), 0],
                [0, 0, 0, 1]
            ])
        elif rot_axis[1] == 1:
             rot = np.array([
                [np.cos(np.radians(rotation_amount)), 0, np.sin(np.radians(rotation_amount)), 0],
                [0, 1, 0, 0],
                [-np.sin(np.radians(rotation_amount)), 0, np.cos(np.radians(rotation_amount)), 0],
                [0, 0, 0, 1]
            ])
        else:
             rot = np.array([
                [np.cos(np.radians(rotation_amount)), -np.sin(np.radians(rotation_amount)), 0, 0],
                [np.sin(np.radians(rotation_amount)), np.cos(np.radians(rotation_amount)), 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
            ])
             
        if rot_axis2[0] == 1:
            rot2 = np.array([
                [1, 0, 0, 0],
                [0, np.cos(np.radians(rotation_amount2)), -np.sin(np.radians(rotation_amount2)), 0],
                [0, np.sin(np.radians(rotation_amount2)), np.cos(np.radians(rotation_amount2)), 0],
                [0, 0, 0, 1]
            ])
        elif rot_axis2[1] == 1:
             rot2 = np.array([
                [np.cos(np.radians(rotation_amount2)), 0, np.sin(np.radians(rotation_amount2)), 0],
                [0, 1, 0, 0],
                [-np.sin(np.radians(rotation_amount2)), 0, np.cos(np.radians(rotation_amount2)), 0],
                [0, 0, 0, 1]
            ])
        else:
             rot2 = np.array([
                [np.cos(np.radians(rotation_amount2)), -np.sin(np.radians(rotation_amount2)), 0, 0],
                [np.sin(np.radians(rotation_amount2)), np.cos(np.radians(rotation_amount2)), 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
            ])
             
        trans1 = np.eye(4)
        trans1[trans_axis1, 3] = translation_amount1

        np.testing.assert_allclose(dipoles_rot_trans_rot.get_combined_transformation(), rot2 @ trans1 @ rot, atol=1e-5)
        np.testing.assert_allclose(dipoles_trans_rot_rot.get_combined_transformation(), rot2 @ rot @ trans1, atol=1e-5)
        np.testing.assert_allclose(dipoles_rot_rot_trans.get_combined_transformation(), trans1 @ rot2 @ rot, atol=1e-5)




class TestGetMesh:
    def test_dipole_element_mesh(self):
        dipole_locations = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        dipole_moments = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]])
        dipoles = DipoleElements(
            None,
            dipole_locations,
            dipole_moments,
        )
        mesh = dipoles.get_mesh(np.eye(4), False, False, False, True)

        assert mesh.elm.node_number_list.shape[0] == 4
        np.testing.assert_allclose(np.unique(mesh.elm.elm_type), [15])
        assert np.all(mesh.elm.tag1 == TmsCoilElementTag.DIPOLES)
        np.testing.assert_allclose(
            mesh.nodes.node_coord, dipole_locations
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, dipole_moments
        )

    def test_dipole_element_mesh_rotation(self):
        dipole_locations = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        dipole_moments = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]])
        dipoles = DipoleElements(
            None,
            dipole_locations,
            dipole_moments,
            deformations=[TmsCoilRotation(90, (-5000, 5000), np.array([0, 0, 0]), np.array([1, 0, 0]))],
        )
        mesh = dipoles.get_mesh(np.eye(4), True, False, False, True)

        rotated_points = np.zeros_like(dipole_locations)
        rotated_points[:, 0] = dipole_locations[:, 0]
        rotated_points[:, 1] = -dipole_locations[:, 2]
        rotated_points[:, 2] = dipole_locations[:, 1]

        rotated_moments = np.zeros_like(dipole_moments)
        rotated_moments[:, 0] = dipole_moments[:, 0]
        rotated_moments[:, 1] = -dipole_moments[:, 2]
        rotated_moments[:, 2] = dipole_moments[:, 1]

        np.testing.assert_allclose(
            mesh.nodes.node_coord, rotated_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, rotated_moments, atol=1e-5
        )

    def test_dipole_element_mesh_translation(self):
        dipole_locations = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        dipole_moments = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]])
        dipoles = DipoleElements(
            None,
            dipole_locations,
            dipole_moments,
            deformations=[TmsCoilTranslation(123, (-200, 200), 0)],
        )
        mesh = dipoles.get_mesh(np.eye(4), True, False, False, True)

        translated_points = np.array(dipole_locations)
        translated_points[:, 0] = dipole_locations[:, 0] + 123

        np.testing.assert_allclose(
            mesh.nodes.node_coord, translated_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, dipole_moments, atol=1e-5
        )

    def test_dipole_element_mesh_affine(self):
        affine = np.array([
            [ 0.25  ,0.5   ,0.25  ,2 ],
            [ 0.25  ,0.5   ,0.25  ,3 ],
            [ 0.25  ,0.5   ,0.25  ,4 ],
            [  0    , 0    , 0    ,1 ]
        ])
        dipole_locations = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        dipole_moments = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]])
        dipoles = DipoleElements(
            None,
            dipole_locations,
            dipole_moments,
            deformations=[TmsCoilTranslation(123, (-200, 200), 0)],
        )
        mesh = dipoles.get_mesh(affine, False, False, False, True)

        transformed_points = dipole_locations @ affine[:3, :3].T + affine[None, :3, 3]
        transformed_moments = dipole_moments @ affine[:3, :3].T

        np.testing.assert_allclose(
            mesh.nodes.node_coord, transformed_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, transformed_moments, atol=1e-5
        )

    def test_dipole_element_mesh_deformation_and_affine(self):
        affine = np.array([
            [ 0.25  ,0.5   ,0.25  ,121 ],
            [ 0.25  ,0.5   ,0.25  ,-2 ],
            [ 0.25  ,0.5   ,0.25  ,4 ],
            [  0    , 0    , 0    ,1 ]
        ])
        dipole_locations = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        dipole_moments = np.array([[1.0, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]])
        dipoles = DipoleElements(
            None,
            dipole_locations,
            dipole_moments,
            deformations=[TmsCoilRotation(90, (-5000, 5000), np.array([0, 0, 0]), np.array([0, 1, 0])), TmsCoilTranslation(-33.33, (-200, 200), 2)],
        )
        mesh = dipoles.get_mesh(affine, True, False, False, True)

        transformed_points = np.zeros_like(dipole_locations)
        transformed_points[:, 0] = dipole_locations[:, 2]
        transformed_points[:, 1] = dipole_locations[:, 1]
        transformed_points[:, 2] = -dipole_locations[:, 0]

        transformed_moments = np.zeros_like(dipole_moments)
        transformed_moments[:, 0] = dipole_moments[:, 2]
        transformed_moments[:, 1] = dipole_moments[:, 1]
        transformed_moments[:, 2] = -dipole_moments[:, 0]

        transformed_points[:, 2] = transformed_points[:, 2] - 33.33

        transformed_points = transformed_points @ affine[:3, :3].T + affine[None, :3, 3]
        transformed_moments = transformed_moments @ affine[:3, :3].T

        np.testing.assert_allclose(
            mesh.nodes.node_coord, transformed_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, transformed_moments, atol=1e-5
        )

    def test_line_segment_element_mesh(self):
        line_segment = LineSegmentElements(
            None,
            np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        )
        mesh = line_segment.get_mesh(np.eye(4), False, False, False, True, 1)

        assert mesh.elm.node_number_list.shape[0] == 4
        np.testing.assert_allclose(np.unique(mesh.elm.elm_type), [1])
        assert np.all(mesh.elm.tag1 == 100 + TmsCoilElementTag.LINE_ELEMENTS)
        np.testing.assert_allclose(
            mesh.nodes.node_coord,
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]],
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, [[1, 0, 0], [-1, 1, 0], [0, -1, 1], [0, 0, -1], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0],]
        )

    def test_line_segment_element_mesh_translation(self):
        line_segment_locations = np.array([[0, 0.0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        line_segment_directions = np.array([[1, 0.0, 0], [-1, 1, 0], [0, -1, 1], [0, 0, -1]])
        line_segment = LineSegmentElements(
            None,
            line_segment_locations,
            deformations=[TmsCoilTranslation(45.34, (-200, 200), 1)],
        )
        mesh = line_segment.get_mesh(np.eye(4), True, False, False, True)

        translated_points = np.concatenate((np.array(line_segment_locations), line_segment_locations + line_segment_directions))
        translated_points[:, 1] = translated_points[:, 1] + 45.34

        np.testing.assert_allclose(
            mesh.nodes.node_coord, translated_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, np.concatenate((line_segment_directions, np.zeros_like(line_segment_directions))), atol=1e-5
        )

    def test_line_segment_element_mesh_rotation(self):
        line_segment_locations = np.array([[0, 0.0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        line_segment_directions = np.array([[1, 0.0, 0], [-1, 1, 0], [0, -1, 1], [0, 0, -1]])
        line_segment = LineSegmentElements(
            None,
            line_segment_locations,
            deformations=[TmsCoilRotation(-90, (-18-5000, 5000), np.array([0, 0, 0]), np.array([0, 0, 1]))],
        )
        mesh = line_segment.get_mesh(np.eye(4), True, False, False, True)

        points = np.concatenate((np.array(line_segment_locations), line_segment_locations + line_segment_directions))
        rotated_points = np.zeros_like(points)
        rotated_points[:, 0] = points[:, 1]
        rotated_points[:, 1] = -points[:, 0]
        rotated_points[:, 2] = points[:, 2]

        rotated_directions = np.zeros_like(line_segment_directions)
        rotated_directions[:, 0] = line_segment_directions[:, 1]
        rotated_directions[:, 1] = -line_segment_directions[:, 0]
        rotated_directions[:, 2] = line_segment_directions[:, 2]

        np.testing.assert_allclose(
            mesh.nodes.node_coord, rotated_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, np.concatenate((rotated_directions, np.zeros_like(rotated_directions))), atol=1e-5
        )

    def test_line_segment_element_mesh_affine(self):
        affine = np.array([
            [ 0.25  ,0.5   ,0.25  ,121 ],
            [ 0.25  ,0.5   ,0.25  ,-2 ],
            [ 0.25  ,0.5   ,0.25  ,4 ],
            [  0    , 0    , 0    ,1 ]
        ])
        line_segment_locations = np.array([[0, 0.0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        line_segment_directions = np.array([[1, 0.0, 0], [-1, 1, 0], [0, -1, 1], [0, 0, -1]])
        line_segment = LineSegmentElements(
            None,
            line_segment_locations,
            deformations=[TmsCoilRotation(-90, (-18-5000, 5000), np.array([0, 0, 0]), np.array([0, 0, 1]))],
        )
        mesh = line_segment.get_mesh(affine, False, False, False, True)

        points = np.concatenate((np.array(line_segment_locations), line_segment_locations + line_segment_directions))
        
        transformed_points = points @ affine[:3, :3].T + affine[None, :3, 3]
        transformed_direction = line_segment_directions @ affine[:3, :3].T


        np.testing.assert_allclose(
            mesh.nodes.node_coord, transformed_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, np.concatenate((transformed_direction, np.zeros_like(transformed_direction))), atol=1e-5
        )

    def test_line_segment_element_mesh_deformation_and_affine(self):
        affine = np.array([
            [ 0.25  ,0.5   ,0.25  ,121 ],
            [ 0.25  ,0.5   ,0.25  ,-2 ],
            [ 0.25  ,0.5   ,0.25  ,4 ],
            [  0    , 0    , 0    ,1 ]
        ])
        line_segment_locations = np.array([[0, 0.0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        line_segment_directions = np.array([[1, 0.0, 0], [-1, 1, 0], [0, -1, 1], [0, 0, -1]])
        line_segment = LineSegmentElements(
            None,
            line_segment_locations,
            deformations=[TmsCoilRotation(90, (-5000, 5000), np.array([0, 0, 0]), np.array([0, 1, 0])), TmsCoilTranslation(22.22, (-200, 200), 2)],
        )
        mesh = line_segment.get_mesh(affine, True, False, False, True)

        points = np.concatenate((np.array(line_segment_locations), line_segment_locations + line_segment_directions))
        
        transformed_points = np.zeros_like(points)
        transformed_points[:, 0] = points[:, 2]
        transformed_points[:, 1] = points[:, 1]
        transformed_points[:, 2] = -points[:, 0]

        transformed_direction = np.zeros_like(line_segment_directions)
        transformed_direction[:, 0] = line_segment_directions[:, 2]
        transformed_direction[:, 1] = line_segment_directions[:, 1]
        transformed_direction[:, 2] = -line_segment_directions[:, 0]

        transformed_points[:, 2] = transformed_points[:, 2] + 22.22

        transformed_points = transformed_points @ affine[:3, :3].T + affine[None, :3, 3]
        transformed_direction = transformed_direction @ affine[:3, :3].T


        np.testing.assert_allclose(
            mesh.nodes.node_coord, transformed_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, np.concatenate((transformed_direction, np.zeros_like(transformed_direction))), atol=1e-5
        )

    def test_sampled_grid_elements_element_mesh(self):
        data = np.array([[[[1,0.0,0],[0,1,0]], [[0,0,1],[1,1,0]]], [[[0,1,1],[1,0,1]], [[1,1,1],[1,0,0]]]])
        sample_affine = np.eye(4)
        sample_locations = np.array([[0,0.0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]])
        sampled_element = SampledGridPointElements(
            None,
            data,
            sample_affine,
        )
        mesh = sampled_element.get_mesh(np.eye(4), False, False, False, True)

        assert mesh.elm.node_number_list.shape[0] == 8
        np.testing.assert_allclose(np.unique(mesh.elm.elm_type), [15])
        assert np.all(mesh.elm.tag1 == TmsCoilElementTag.SAMPLED_GRID_ELEMENTS)
        np.testing.assert_allclose(
            mesh.nodes.node_coord, sample_locations
        )

        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, data.reshape(-1, 3)
        )

    def test_sampled_grid_elements_mesh_rotation(self):
        data = np.array([[[[1,0.0,0],[0,1,0]], [[0,0,1],[1,1,0]]], [[[0,1,1],[1,0,1]], [[1,1,1],[1,0,0]]]])
        sample_affine = np.eye(4)
        sample_locations = np.array([[0,0.0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]])
        sampled_element = SampledGridPointElements(
            None,
            data,
            sample_affine,  
            deformations=[TmsCoilRotation(90, (-5000, 5000), np.array([0, 0, 0]), np.array([1, 0, 0]))],
        )
        mesh = sampled_element.get_mesh(np.eye(4), True, False, False, True)

        rotated_points = np.zeros_like(sample_locations)
        rotated_points[:, 0] = sample_locations[:, 0]
        rotated_points[:, 1] = -sample_locations[:, 2]
        rotated_points[:, 2] = sample_locations[:, 1]

        location_data = data.reshape(-1, 3)
        rotated_data = np.zeros_like(location_data)
        rotated_data[:, 0] = location_data[:, 0]
        rotated_data[:, 1] = -location_data[:, 2]
        rotated_data[:, 2] = location_data[:, 1]

        np.testing.assert_allclose(
            mesh.nodes.node_coord, rotated_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, rotated_data, atol=1e-5
        )

    def test_sampled_grid_elements_mesh_translation(self):
        data = np.array([[[[1,0.0,0],[0,1,0]], [[0,0,1],[1,1,0]]], [[[0,1,1],[1,0,1]], [[1,1,1],[1,0,0]]]])
        sample_affine = np.eye(4)
        sample_locations = np.array([[0,0.0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]])
        sampled_element = SampledGridPointElements(
            None,
            data,
            sample_affine, 
            deformations=[TmsCoilTranslation(12.21, (-200, 200), 0)],
        )
        mesh = sampled_element.get_mesh(np.eye(4), True, False, False, True)

        translated_points = np.array(sample_locations)
        translated_points[:, 0] = sample_locations[:, 0] + 12.21

        np.testing.assert_allclose(
            mesh.nodes.node_coord, translated_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, data.reshape(-1, 3), atol=1e-5
        )

    def test_sampled_grid_elements_mesh_affine(self):
        affine = np.array([
            [ 0.25  ,0.5   ,0.25  ,2 ],
            [ 0.25  ,0.5   ,0.25  ,3 ],
            [ 0.25  ,0.5   ,0.25  ,4 ],
            [  0    , 0    , 0    ,1 ]
        ])
        data = np.array([[[[1,0.0,0],[0,1,0]], [[0,0,1],[1,1,0]]], [[[0,1,1],[1,0,1]], [[1,1,1],[1,0,0]]]])
        sample_affine = np.eye(4)
        sample_locations = np.array([[0,0.0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]])
        sampled_element = SampledGridPointElements(
            None,
            data,
            sample_affine,
            deformations=[TmsCoilTranslation(123, (-200, 200), 0)], 
        )
        mesh = sampled_element.get_mesh(affine, False, False, False, True)

        transformed_points = sample_locations @ affine[:3, :3].T + affine[None, :3, 3]
        transformed_moments = data.reshape(-1, 3) @ affine[:3, :3].T

        np.testing.assert_allclose(
            mesh.nodes.node_coord, transformed_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, transformed_moments, atol=1e-5
        )

    def test_sampled_grid_elements_mesh_data_affine(self):
        affine = np.array([
            [ 0.25  ,0.5   ,0.25  ,2 ],
            [ 0.25  ,0.5   ,0.25  ,3 ],
            [ 0.25  ,0.5   ,0.25  ,4 ],
            [  0    , 0    , 0    ,1 ]
        ])
        data = np.array([[[[1,0.0,0],[0,1,0]], [[0,0,1],[1,1,0]]], [[[0,1,1],[1,0,1]], [[1,1,1],[1,0,0]]]])
        sample_locations = np.array([[0,0.0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]])
        sampled_element = SampledGridPointElements(
            None,
            data,
            affine,
            deformations=[TmsCoilTranslation(123, (-200, 200), 0)],
        )
        mesh = sampled_element.get_mesh(np.eye(4), False, False, False, True)

        transformed_points = sample_locations @ affine[:3, :3].T + affine[None, :3, 3]

        np.testing.assert_allclose(
            mesh.nodes.node_coord, transformed_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, data.reshape(-1, 3), atol=1e-5
        )

    def test_sampled_grid_elements_mesh_deformation_and_affine(self):
        affine = np.array([
            [ 0.25  ,0.5   ,0.25  ,121 ],
            [ 0.25  ,0.5   ,0.25  ,-2 ],
            [ 0.25  ,0.5   ,0.25  ,4 ],
            [  0    , 0    , 0    ,1 ]
        ])
        data = np.array([[[[1,0.0,0],[0,1,0]], [[0,0,1],[1,1,0]]], [[[0,1,1],[1,0,1]], [[1,1,1],[1,0,0]]]])
        sample_affine = np.eye(4)
        sample_locations = np.array([[0,0.0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]])
        sampled_element = SampledGridPointElements(
            None,
            data,
            sample_affine,
            deformations=[TmsCoilRotation(90, (-5000, 5000), np.array([0, 0, 0]), np.array([0, 1, 0])), TmsCoilTranslation(-33.33, (-200, 200), 2)],
        )
        mesh = sampled_element.get_mesh(affine, True, False, False, True)

        transformed_points = np.zeros_like(sample_locations)
        transformed_points[:, 0] = sample_locations[:, 2]
        transformed_points[:, 1] = sample_locations[:, 1]
        transformed_points[:, 2] = -sample_locations[:, 0]

        location_data = data.reshape(-1, 3)
        transformed_data = np.zeros_like(location_data)
        transformed_data[:, 0] = location_data[:, 2]
        transformed_data[:, 1] = location_data[:, 1]
        transformed_data[:, 2] = -location_data[:, 0]

        transformed_points[:, 2] = transformed_points[:, 2] - 33.33

        transformed_points = transformed_points @ affine[:3, :3].T + affine[None, :3, 3]
        transformed_moments = transformed_data @ affine[:3, :3].T

        np.testing.assert_allclose(
            mesh.nodes.node_coord, transformed_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, transformed_moments, atol=1e-5
        )

    def test_sampled_grid_elements_mesh_deformation_and_affine_and_affine(self):
        affine = np.array([
            [ 0.25  ,0.5   ,0.25  ,121 ],
            [ 0.25  ,0.5   ,0.25  ,-2 ],
            [ 0.25  ,0.5   ,0.25  ,4 ],
            [  0    , 0    , 0    ,1 ]
        ])
        data = np.array([[[[1,0.0,0],[0,1,0]], [[0,0,1],[1,1,0]]], [[[0,1,1],[1,0,1]], [[1,1,1],[1,0,0]]]])
        sample_affine = np.array([
            [ 0.36  ,0.12   ,0.89  ,-22 ],
            [ 0.36  ,0.12   ,0.89  ,-10 ],
            [ 0.36  ,0.12   ,0.89  ,-30 ],
            [  0    , 0    , 0    ,1 ]
        ])
        sample_locations = np.array([[0,0.0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]])
        sampled_element = SampledGridPointElements(
            None,
            data,
            sample_affine,
            deformations=[TmsCoilRotation(90, (-5000, 5000), np.array([0, 0, 0]), np.array([0, 1, 0])), TmsCoilTranslation(-33.33, (-200, 200), 2)],
        )
        mesh = sampled_element.get_mesh(affine, True, False, False, True)

        transformed_locations = sample_locations @ sample_affine[:3, :3].T + sample_affine[None, :3, 3]

        transformed_points = np.zeros_like(transformed_locations)
        transformed_points[:, 0] = transformed_locations[:, 2]
        transformed_points[:, 1] = transformed_locations[:, 1]
        transformed_points[:, 2] = -transformed_locations[:, 0]

        location_data = data.reshape(-1, 3)
        transformed_data = np.zeros_like(location_data)
        transformed_data[:, 0] = location_data[:, 2]
        transformed_data[:, 1] = location_data[:, 1]
        transformed_data[:, 2] = -location_data[:, 0]

        transformed_points[:, 2] = transformed_points[:, 2] - 33.33

        transformed_points = transformed_points @ affine[:3, :3].T + affine[None, :3, 3]
        transformed_moments = transformed_data @ affine[:3, :3].T

        np.testing.assert_allclose(
            mesh.nodes.node_coord, transformed_points, atol=1e-5
        )
        assert len(mesh.field) == 1
        np.testing.assert_allclose(
            mesh.field[list(mesh.field.keys())[0]].value, transformed_moments, atol=1e-5
        )
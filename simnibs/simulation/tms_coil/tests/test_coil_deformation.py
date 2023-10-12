import numpy as np
import pytest

from simnibs.simulation.tms_coil.tms_coil_deformation import (
    TmsCoilDeformationRange,
    TmsCoilRotation,
    TmsCoilTranslation,
)


class TestInit:
    def test_wrong_shape_range(self):
        with pytest.raises(ValueError):
            TmsCoilTranslation(TmsCoilDeformationRange(0, [0, 2, 3]), 0)

        with pytest.raises(ValueError):
            TmsCoilTranslation(TmsCoilDeformationRange(0, [[0, 2, 3]]), 0)

    def test_min_greater_max_range(self):
        with pytest.raises(ValueError):
            TmsCoilTranslation(TmsCoilDeformationRange(0, [0, -3]), 0)

    def test_initial_out_of_range(self):
        with pytest.raises(ValueError):
            TmsCoilTranslation(TmsCoilDeformationRange(0, [1, 5]), 0)

    def test_axis_out_of_range(self):
        with pytest.raises(ValueError):
            TmsCoilTranslation(TmsCoilDeformationRange(3, [1, 5]), 4)

    def test_shape_point_1(self):
        with pytest.raises(ValueError):
            TmsCoilRotation(TmsCoilDeformationRange(3, [1, 5]), [], [3, 4, 5])

        with pytest.raises(ValueError):
            TmsCoilRotation(TmsCoilDeformationRange(3, [1, 5]), [1, 2], [3, 4, 5])

        with pytest.raises(ValueError):
            TmsCoilRotation(TmsCoilDeformationRange(3, [1, 5]), [1, 2, 4, 4], [3, 4, 5])

    def test_shape_point_2(self):
        with pytest.raises(ValueError):
            TmsCoilRotation(TmsCoilDeformationRange(3, [1, 5]), [3, 4, 5], [])

        with pytest.raises(ValueError):
            TmsCoilRotation(TmsCoilDeformationRange(3, [1, 5]), [1, 2, 3], [4, 5])

        with pytest.raises(ValueError):
            TmsCoilRotation(TmsCoilDeformationRange(3, [1, 5]), [1, 2, 4], [3, 4, 5, 6])


class TestCurrent:
    def test_reset(self):
        trans = TmsCoilDeformationRange(3, [1, 5])
        trans.current = 4.5
        trans.reset()
        assert trans.current == trans.initial

    def test_current_out_of_range(self):
        trans = TmsCoilDeformationRange(3, [1, 5])
        with pytest.raises(ValueError):
            trans.current = 5.1


class TestCoilTranslation:
    def test_get_translation(self):
        initial = 0.0
        range = (-1.0, 1.0)
        axis = 1
        translation = TmsCoilTranslation(TmsCoilDeformationRange(initial, range), axis)
        translation.deformation_range.current = 0.5
        expected_translation = np.array([0.0, 0.5, 0.0])
        np.testing.assert_allclose(translation.get_translation(), expected_translation)

    def test_apply_basic(self):
        translation = TmsCoilTranslation(TmsCoilDeformationRange(0, (-10.0, 10.0)), 1)
        translation.deformation_range.current = 3
        points = np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [0, 2.3, 0],
                [0, 0, 10.1],
                [4, 0.3, 0],
                [3.7, 0, 20.4],
                [0, 1.376, 1],
                [9, 21, 45],
                [-21, 0, 0],
                [0, -2.3, 0],
                [0, 0, 10.1],
                [4, -0.3, 0],
                [-3.7, 0, -20.4],
                [0, 1.376, -1],
                [-9, -21, 45],
            ]
        )
        expected_result = np.copy(points)
        expected_result[:, 1] += 3
        np.testing.assert_allclose(translation.apply(points), expected_result)

    def test_apply_negative(self):
        translation = TmsCoilTranslation(TmsCoilDeformationRange(0, (-1.0, 1.0)), 2)
        translation.deformation_range.current = -1
        points = np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [0, 2.3, 0],
                [0, 0, 10.1],
                [4, 0.3, 0],
                [3.7, 0, 20.4],
                [0, 1.376, 1],
                [9, 21, 45],
                [-21, 0, 0],
                [0, -2.3, 0],
                [0, 0, 10.1],
                [4, -0.3, 0],
                [-3.7, 0, -20.4],
                [0, 1.376, -1],
                [-9, -21, 45],
            ]
        )
        expected_result = np.copy(points)
        expected_result[:, 2] -= 1
        np.testing.assert_allclose(translation.apply(points), expected_result)

    def test_apply_zero(self):
        translation = TmsCoilTranslation(TmsCoilDeformationRange(0, (-1.0, 1.0)), 2)
        translation.deformation_range.current = 0
        points = np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [0, 2.3, 0],
                [0, 0, 10.1],
                [4, 0.3, 0],
                [3.7, 0, 20.4],
                [0, 1.376, 1],
                [9, 21, 45],
                [-21, 0, 0],
                [0, -2.3, 0],
                [0, 0, 10.1],
                [4, -0.3, 0],
                [-3.7, 0, -20.4],
                [0, 1.376, -1],
                [-9, -21, 45],
            ]
        )
        expected_result = np.copy(points)
        np.testing.assert_allclose(translation.apply(points), expected_result)

    def test_apply_fractional(self):
        translation = TmsCoilTranslation(TmsCoilDeformationRange(0, (-1.0, 1.0)), 0)
        translation.deformation_range.current = 0.3232323
        points = np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [0, 2.3, 0],
                [0, 0, 10.1],
                [4, 0.3, 0],
                [3.7, 0, 20.4],
                [0, 1.376, 1],
                [9, 21, 45],
                [-21, 0, 0],
                [0, -2.3, 0],
                [0, 0, 10.1],
                [4, -0.3, 0],
                [-3.7, 0, -20.4],
                [0, 1.376, -1],
                [-9, -21, 45],
            ]
        )
        expected_result = np.copy(points)
        expected_result[:, 0] += 0.3232323
        np.testing.assert_allclose(translation.apply(points), expected_result)

    def test_apply_large(self):
        translation = TmsCoilTranslation(
            TmsCoilDeformationRange(0, (-10000000.0, 10000000.0)), 1
        )
        translation.deformation_range.current = 3182382
        points = np.array(
            [
                [0, 0, 0],
                [1, 0, 0],
                [0, 2.3, 0],
                [0, 0, 10.1],
                [4, 0.3, 0],
                [3.7, 0, 20.4],
                [0, 1.376, 1],
                [9, 21, 45],
                [-21, 0, 0],
                [0, -2.3, 0],
                [0, 0, 10.1],
                [4, -0.3, 0],
                [-3.7, 0, -20.4],
                [0, 1.376, -1],
                [-9, -21, 45],
            ]
        )
        expected_result = np.copy(points)
        expected_result[:, 1] += 3182382
        np.testing.assert_allclose(translation.apply(points), expected_result)

    def test_as_matrix(self):
        initial = 0.0
        range = (-1.0, 1.0)
        axis = 1
        translation = TmsCoilTranslation(TmsCoilDeformationRange(initial, range), axis)
        translation.deformation_range.current = 0.5
        expected_matrix = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.5],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        np.testing.assert_allclose(translation.as_matrix(), expected_matrix)


class TestCoilRotation:
    def test_basic(self):
        rotation = TmsCoilRotation(
            TmsCoilDeformationRange(0, (-199.0, 199.0)),
            np.array([0.0, 0.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        )
        rotation.deformation_range.current = 90.0

        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        expected_result = np.array([[0, 0, 0], [0, 1, 0], [-1, 0, 0], [0, 0, 1]])
        np.testing.assert_allclose(
            rotation.apply(points), expected_result, rtol=1e-9, atol=1e-9
        )

    def test_negative(self):
        rotation = TmsCoilRotation(
            TmsCoilDeformationRange(0, (-199.0, 199.0)),
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
        )
        rotation.deformation_range.current = -45

        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        expected_result = np.array(
            [[0, 0, 0], [1, 0, 0], [0, 0.707, -0.707], [0, 0.707, 0.707]]
        )
        np.testing.assert_allclose(
            rotation.apply(points), expected_result, rtol=1e-3, atol=1e-3
        )

    def test_zero(self):
        rotation = TmsCoilRotation(
            TmsCoilDeformationRange(0, (-199.0, 199.0)),
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 1.0, 1.0]),
        )
        rotation.deformation_range.current = 0

        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        expected_result = np.copy(points)
        np.testing.assert_allclose(
            rotation.apply(points), expected_result, rtol=1e-9, atol=1e-9
        )

    def test_arbitrary_origin_axis(self):
        rotation = TmsCoilRotation(
            TmsCoilDeformationRange(0, (-199.0, 199.0)),
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 1.0, 0.0]),
        )
        rotation.deformation_range.current = 180

        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        expected_result = np.array([[0, 0, 0], [0, 1, 0], [1, 0, 0], [0, 0, -1]])
        np.testing.assert_allclose(
            rotation.apply(points), expected_result, rtol=1e-9, atol=1e-9
        )

    def test_arbitrary_axis(self):
        rotation = TmsCoilRotation(
            TmsCoilDeformationRange(0, (-199.0, 199.0)),
            np.array([1.0, 0.0, 1.0]),
            np.array([1.0, 1.0, 1.0]),
        )
        rotation.deformation_range.current = 180

        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        expected_result = np.array([[2, 0, 2], [1, 0, 2], [2, 1, 2], [2, 0, 1]])
        np.testing.assert_allclose(
            rotation.apply(points), expected_result, rtol=1e-9, atol=1e-9
        )

    def test_as_matrix(self):
        rotation = TmsCoilRotation(
            TmsCoilDeformationRange(0, (-199.0, 199.0)),
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
        )
        rotation.deformation_range.current = 90.0
        expected_rotation = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, -1.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        np.testing.assert_allclose(
            rotation.get_rotation(), expected_rotation, rtol=1e-5, atol=1e-5
        )

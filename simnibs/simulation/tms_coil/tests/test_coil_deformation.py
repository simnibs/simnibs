import numpy as np

from simnibs.simulation.tms_coil.tms_coil_deformation import (
    TmsCoilRotation,
    TmsCoilTranslation,
)


class TestCoilTranslation:
    def test_get_translation(self):
        initial = 0.0
        range = (-1.0, 1.0)
        axis = 1
        translation = TmsCoilTranslation(initial, range, axis)
        translation.current = 0.5
        expected_translation = np.array([0.0, 0.5, 0.0])
        np.testing.assert_allclose(translation.get_translation(), expected_translation)

    def test_apply_basic(self):
        translation = TmsCoilTranslation(0, (-10.0, 10.0), 1)
        translation.current = 3
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
        translation = TmsCoilTranslation(0, (-1.0, 1.0), 2)
        translation.current = -1
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
        translation = TmsCoilTranslation(0, (-1.0, 1.0), 2)
        translation.current = 0
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
        translation = TmsCoilTranslation(0, (-1.0, 1.0), 0)
        translation.current = 0.3232323
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
        translation = TmsCoilTranslation(0, (-10000000.0, 10000000.0), 1)
        translation.current = 3182382
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
        translation = TmsCoilTranslation(initial, range, axis)
        translation.current = 0.5
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
            0, (-1.0, 1.0), np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0])
        )
        rotation.current = 90.0

        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        expected_result = np.array([[0, 0, 0], [0, 1, 0], [-1, 0, 0], [0, 0, 1]])
        np.testing.assert_allclose(
            rotation.apply(points), expected_result, rtol=1e-9, atol=1e-9
        )

    def test_negative(self):
        rotation = TmsCoilRotation(
            0, (-1.0, 1.0), np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])
        )
        rotation.current = -45

        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        expected_result = np.array(
            [[0, 0, 0], [1, 0, 0], [0, 0.707, -0.707], [0, 0.707, 0.707]]
        )
        np.testing.assert_allclose(
            rotation.apply(points), expected_result, rtol=1e-3, atol=1e-3
        )

    def test_zero(self):
        rotation = TmsCoilRotation(
            0, (-1.0, 1.0), np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 1.0])
        )
        rotation.current = 0

        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        expected_result = np.copy(points)
        np.testing.assert_allclose(
            rotation.apply(points), expected_result, rtol=1e-9, atol=1e-9
        )

    def test_arbitrary_origin_axis(self):
        rotation = TmsCoilRotation(
            0, (-1.0, 1.0), np.array([0.0, 0.0, 0.0]), np.array([1.0, 1.0, 0.0])
        )
        rotation.current = 180

        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        expected_result = np.array([[0, 0, 0], [0, 1, 0], [1, 0, 0], [0, 0, -1]])
        np.testing.assert_allclose(
            rotation.apply(points), expected_result, rtol=1e-9, atol=1e-9
        )

    def test_arbitrary_axis(self):
        rotation = TmsCoilRotation(
            0, (-1.0, 1.0), np.array([1.0, 0.0, 1.0]), np.array([1.0, 1.0, 1.0])
        )
        rotation.current = 180

        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        expected_result = np.array([[2, 0, 2], [1, 0, 2], [2, 1, 2], [2, 0, 1]])
        np.testing.assert_allclose(
            rotation.apply(points), expected_result, rtol=1e-9, atol=1e-9
        )

    def test_as_matrix(self):
        rotation = TmsCoilRotation(
            0, (-1.0, 1.0), np.array([0.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])
        )
        rotation.current = 90.0
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

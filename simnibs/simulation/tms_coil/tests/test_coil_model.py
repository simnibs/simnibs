import numpy as np
import pytest

from simnibs.mesh_tools.mesh_io import Elements, Msh, Nodes
from simnibs.simulation.tms_coil.tms_coil_model import TmsCoilModel


class TestInit:
    def test_wrong_shape_min_distance_points(self):
        with pytest.raises(ValueError):
            TmsCoilModel(Msh(), min_distance_points=[1, 2, 3])

        with pytest.raises(ValueError):
            TmsCoilModel(Msh(), min_distance_points=[[1, 2, 3, 4]])

        with pytest.raises(ValueError):
            TmsCoilModel(Msh(), min_distance_points=[[1, 2], [3, 4]])


class TestGetPoints:
    def test_get_points_empty_mesh(self):
        casing = TmsCoilModel(Msh())
        assert len(casing.get_points()) == 0

    def test_get_min_distance_points_no_transformation(self):
        casing = TmsCoilModel(Msh(), min_distance_points=[[1, 2, 3], [4, 5, 6]])

        np.testing.assert_allclose(
            casing.get_min_distance_points(), [[1, 2, 3], [4, 5, 6]]
        )

    def test_apply_deformations(self):
        casing = TmsCoilModel(
            Msh(Nodes([[1, 2, 3], [4, 5, 6], [7, 8, 9]]), Elements()),
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        )

        affine = np.array([[2, 0, 0, 1], [0, 1, 0, 2], [0, 0, 3, 3], [0, 0, 0, 1]])

        casing_after = casing.apply_deformations(affine)

        np.testing.assert_allclose(
            casing_after.mesh.nodes.node_coord,
            casing.mesh.nodes.node_coord @ affine[:3, :3].T + affine[None, :3, 3],
        )
        np.testing.assert_allclose(
            casing_after.min_distance_points,
            casing.min_distance_points @ affine[:3, :3].T + affine[None, :3, 3],
        )


class TestFromPoints:
    def test_from_points(self):
        x = np.arange(-5, 5 + 0.1, 0.1)
        y = np.arange(-5, 5 + 0.1, 0.1)

        X, Y = np.meshgrid(x, y)
        Z = np.zeros_like(X)

        coords = np.stack((X, Y, Z), axis=-1).reshape(-1, 3)
        print(coords)
        casing = TmsCoilModel.from_points(coords, 1, 0.1)
        aabb_tree = casing.mesh.get_AABBTree()

        assert np.min(aabb_tree.min_sqdist(coords)) > 0.9

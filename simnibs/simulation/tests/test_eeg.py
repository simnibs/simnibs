import numpy as np

from .. import eeg


def test_get_triangle_neighbors():
    """Triangulate an array of points and test neighbors like

        .  .  .
        .  .  .
        .  .  .

    e.g.,

        x, y, z = (-1, 0, 1), (-1, 0, 1), (0, )
        points = np.dstack(np.meshgrid(x, y, z)).reshape(-1,3)

    """
    tris = np.array(
        [
            [0, 3, 1],
            [1, 3, 4],
            [3, 7, 4],
            [3, 6, 7],
            [1, 5, 2],
            [1, 4, 5],
            [4, 7, 5],
            [5, 7, 8],
        ]
    )
    pttris = eeg.get_triangle_neighbors(tris, nr=9)
    pttris_expected = [
        [0],
        [0, 1, 4, 5],
        [4],
        [0, 1, 2, 3],
        [1, 2, 5, 6],
        [4, 5, 6, 7],
        [3],
        [2, 3, 6, 7],
        [7],
    ]
    for i, j in zip(pttris, pttris_expected):
        assert all(i == j)


def test_project_points_to_surface():
    """Make a surface consisting of one triangle and project a point from each
    'region' to it.
    """
    surf = dict(
        points=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]]), tris=np.array([[0, 1, 2]])
    )

    # Test a point in each region
    points = np.array(
        [
            [0.25, 0.25, 0.67],  # Region 0
            [1, 1, 0.27],  # Region 1
            [-0.25, 2, 1.1],  # Region 2
            [-1, 0.5, -1.4],  # Region 3
            [-1, -1, 0.53],  # Region 4
            [0.5, -1, -0.77],  # Region 5
            [2, -0.25, -0.16],
        ]
    )  # Region 6
    pttris = np.atleast_2d(np.zeros(len(points), dtype=int)).T

    # Expected output
    tris = np.zeros(len(points), dtype=int)
    weights = np.array(
        [
            [0.5, 0.25, 0.25],
            [0, 0.5, 0.5],
            [0, 0, 1],
            [0.5, 0, 0.5],
            [1, 0, 0],
            [0.5, 0.5, 0],
            [0, 1, 0],
        ]
    )
    projs = np.array(
        [
            [0.25, 0.25, 0],
            [0.5, 0.5, 0],
            [0, 1, 0],
            [0, 0.5, 0],
            [0, 0, 0],
            [0.5, 0, 0],
            [1, 0, 0],
        ]
    )
    dists = np.linalg.norm(points - projs, axis=1)

    # Actual output
    t, w, p, d = eeg.project_points_to_surface(points, surf, pttris)

    np.testing.assert_array_equal(t, tris)
    np.testing.assert_allclose(w, weights)
    np.testing.assert_allclose(p, projs)
    np.testing.assert_allclose(d, dists)

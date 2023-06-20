from nibabel.affines import apply_affine
import numpy as np
import pytest

from simnibs.utils import spatial_transform


@pytest.fixture
def params():
    return dict(
        rotation={0: (np.pi / 4, 0, 0), 1: (np.pi / 4, np.pi / 4, 0)},
        scaling={0: (0.5, 0.5, 1), 1: (2, 0.5, 2)},
        translation={0: (2, 3, 0), 1: (0, 0, 5)},
    )


@pytest.fixture
def affines():
    return dict(
        rotation={
            0: np.array(
                [
                    [1.0, 0.0, 0.0, 0.0],
                    [0.0, 0.70710678, -0.70710678, 0.0],
                    [0.0, 0.70710678, 0.70710678, 0.0],
                    [0.0, 0.0, 0.0, 1.0],
                ]
            ),
            1: np.array(
                [
                    [0.70710678, 0.5, 0.5, 0.0],
                    [0.0, 0.70710678, -0.70710678, 0.0],
                    [-0.70710678, 0.5, 0.5, 0.0],
                    [0.0, 0.0, 0.0, 1.0],
                ]
            ),
        },
        scaling={
            0: np.array(
                [
                    [0.5, 0.0, 0.0, 0.0],
                    [0.0, 0.5, 0.0, 0.0],
                    [0.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0],
                ]
            ),
            1: np.array(
                [
                    [2.0, 0.0, 0.0, 0.0],
                    [0.0, 0.5, 0.0, 0.0],
                    [0.0, 0.0, 2.0, 0.0],
                    [0.0, 0.0, 0.0, 1.0],
                ]
            ),
        },
        translation={
            0: np.array([[1, 0, 0, 2], [0, 1, 0, 3], [0, 0, 1, 0], [0, 0, 0, 1]]),
            1: np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 5], [0, 0, 0, 1]]),
        },
    )


@pytest.fixture
def src_pts():
    rng = np.random.default_rng(0)
    return rng.uniform(0, 10, (5, 3))


@pytest.fixture
def match_params():
    return {
        False: (0.2, 0.3, 0.4, 1, 3, 5),
        True: (0.2, 0.3, 0.4, 1, 3, 5, 1.5, 1.5, 1.5),
    }


def test_rotation(params, affines):
    k = "rotation"
    for i in params[k]:
        np.testing.assert_allclose(
            spatial_transform.rotation(*params[k][i]), affines[k][i], atol=1e-9
        )


def test_scaling(params, affines):
    k = "scaling"
    for i in params[k]:
        np.testing.assert_allclose(
            spatial_transform.scaling(*params[k][i]), affines[k][i]
        )


def test_translation(params, affines):
    k = "translation"
    for i in params[k]:
        np.testing.assert_allclose(
            spatial_transform.translation(*params[k][i]), affines[k][i]
        )


@pytest.mark.parametrize("scale", [False, True])
def fit_matched_points_analytical(scale, src_pts, match_params):
    """Recover true transformation matrix."""

    trans = spatial_transform.matrix_from_params(*match_params[scale])
    tgt_pts = apply_affine(trans, src_pts)

    out = spatial_transform.fit_matched_points_analytical(src_pts, tgt_pts, scale=scale)

    np.testing.assert_allclose(tgt_pts, apply_affine(out, src_pts))


@pytest.mark.parametrize("scale", [False, True])
def fit_matched_points_generic(scale, src_pts, match_params):
    """Recover true transformation matrix."""

    trans = spatial_transform.matrix_from_params(*match_params[scale])
    tgt_pts = apply_affine(trans, src_pts)

    out = spatial_transform.fit_matched_points_generic(src_pts, tgt_pts, scale=scale)

    np.testing.assert_allclose(tgt_pts, apply_affine(out, src_pts))

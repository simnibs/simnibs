import numpy as np

from .. import _thickness

def create_rings(radii, img_size):
    img_size = 100

    coords = np.meshgrid(
        np.arange(img_size) - img_size/2,
        np.arange(img_size) - img_size/2,
        np.arange(img_size) - img_size/2,
        indexing='xy'
    )

    R = np.sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2)

    rings = np.zeros((img_size, img_size, img_size), dtype=np.uint8)
    for i, r_inner, r_outer in zip(range(len(radii)), radii[:-1], radii[1:]):
        rings[(R > r_inner) * (R < r_outer)] = i + 1

    return rings


def test_thickness_ring():
    outer_radius = 50
    inner_radius = 30
    img_size = 150

    coords = np.meshgrid(
        np.arange(img_size) - img_size/2,
        np.arange(img_size) - img_size/2,
        indexing='xy'
    )

    R = np.sqrt(coords[0]**2 + coords[1]**2)

    ring = (R < outer_radius) * (R > inner_radius)

    thick = _thickness._thickness_slice(ring.astype(np.uint8))
    assert np.allclose(thick[ring], 10, rtol=1e-1)

def test_calc_thickness():
    rings = create_rings([10, 20, 26], 60)
    thick = _thickness._calc_thickness(rings)
    # take a slab
    thick = thick[:, : 20:40]
    rings = rings[:, : 20:40]
    assert np.allclose(thick[rings == 1], 5, atol=1)
    assert np.allclose(thick[rings == 2], 3, atol=1)

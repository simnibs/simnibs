import numpy as np

from simnibs.segmentation import _thickness


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

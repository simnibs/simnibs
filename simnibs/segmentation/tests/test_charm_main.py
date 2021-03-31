import numpy as np

from .. import charm_main

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

class TestMeshing:
    def test_sizing_field_from_thicknes(self):
        thickness = np.arange(100)
        sf = charm_main._sizing_field_from_thickness(
            thickness, 2, (50, 150)
        )
        assert np.isclose(np.min(sf), 50)
        assert np.isclose(np.max(sf), 150)
        assert np.allclose(sf[25:75], np.arange(50, 150, 2))
        assert sf.flags['F_CONTIGUOUS']
        assert sf.dtype == np.float32

    def test_mesh(self):
        rings = create_rings([10, 15, 20, 25], 60)
        rings[rings == 2] = 4
        m = charm_main.create_mesh(rings, np.eye(4))
        # volumes
        vols = m.elements_volumes_and_areas()
        assert np.isclose(
            np.sum(vols[m.elm.tag1 == 1]),
            4/3*np.pi*(15**3-10**3), rtol=1e-1
        )
        assert np.isclose(
            np.sum(vols[m.elm.tag1 == 4]),
            4/3*np.pi*(20**3-15**3), rtol=1e-1
        )
        assert np.isclose(
            np.sum(vols[m.elm.tag1 == 3]),
            4/3*np.pi*(25**3-20**3), rtol=1e-1
        )

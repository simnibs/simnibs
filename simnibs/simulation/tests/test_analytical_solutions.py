from __future__ import division


import pytest
import numpy as np

from .. import analytical_solutions
import scipy.special


@pytest.fixture(scope='module')
def sensor_setup(nr_theta=4, nr_phi=4, r_min=0.1, r_max=.9, nr_r=5):
    theta = np.linspace(0, np.pi, nr_theta, endpoint=False)
    phi = np.linspace(0, 2 * np.pi, nr_phi, endpoint=False)
    r = np.linspace(r_min, r_max, nr_r, endpoint=False)
    theta, phi, r = np.meshgrid(theta, phi, r)
    x = (r * np.sin(theta) * np.cos(phi)).reshape(-1,)
    y = (r * np.sin(theta) * np.sin(phi)).reshape(-1,)
    z = (r * np.cos(theta)).reshape(-1,)
    return np.array((x, y, z), dtype=float).T


@pytest.fixture
def quadrant_points():
    q = np.array([[1, 1, 1],
                  [-1, 1, 1],
                  [-1, -1, 1],
                  [1, -1, 1],
                  [1, 1, -1],
                  [-1, 1, -1],
                  [-1, -1, -1],
                  [1, -1, -1]], dtype=float)
    return q / np.sqrt(3)


class Testlpmn:
    @pytest.mark.parametrize("m", [1, 10])
    @pytest.mark.parametrize("n", [10, 100])
    def test_lpmn_scipy(self, n, m):
        """Test lpmn against the one in scipy.special."""

        x = np.linspace(-1, 1, 100)
        A = np.stack([scipy.special.lpmn(m, n, j)[0] for j in x], axis=2)
        B = analytical_solutions.lpmn(m, n, x)
        np.testing.assert_allclose(A, B)

    def test_lpmn_analytical(self):
        """Test first few degree polynomials."""
        m = 1
        n = 5
        x = np.linspace(-1, 1, 100)
        P = analytical_solutions.lpmn(m, n, x)

        np.testing.assert_allclose(P[0,0], 1.0)
        np.testing.assert_allclose(P[0,1], x)
        np.testing.assert_allclose(P[0,2], 0.5 * (3*x**2 - 1))
        np.testing.assert_allclose(P[0,3], 0.5 * (5 * x**3 - 3 * x))
        np.testing.assert_allclose(P[0,4], 1/8.0 * (35 * x**4 - 30 * x**2 + 3))
        np.testing.assert_allclose(P[0,5], 1/8.0 * (63 * x**5 - 70 * x**3 + 15 * x))


class TestSphereSurfaceElectrodes:

    def test_sphere_homogeneous_cond(self, sensor_setup):
        """ The potential should be the same, independent of the radii  """
        cathode = [1, 0, 0]
        anode = [-1, 0, 0]
        cond = [1, 1, 1]
        radii = [0.7, 0.85, 1]
        pot1 = analytical_solutions.potential_3layers_surface_electrodes(
            radii, cond, anode, cathode, sensor_setup)
        radii = [0.1, 0.2, 1]
        pot2 = analytical_solutions.potential_3layers_surface_electrodes(
            radii, cond, anode, cathode, sensor_setup)

        np.testing.assert_allclose(pot1, pot2)

    def test_sphere_pot_continuity(self):
        """ V should be continuous across interfaces """
        cathode = [1, 0, 0]
        anode = [-1, 0, 0]

        cond = [1.0, 0.5, 1]
        radii = [0.4, 0.5, 1]
        sensor = np.array(((0.5000001, 0, 0),
                           (0.4999999, 0, 0)), dtype=float)
        pot = analytical_solutions.potential_3layers_surface_electrodes(
            radii, cond, anode, cathode, sensor)

        np.testing.assert_allclose(pot[0], pot[1], rtol=1e-3)

    def test_sphere_pot_symmetry_xy(self, sensor_setup):
        """ V should be symmetric to a reflection along the XZ plane """
        cathode = [1, 0, 0]
        anode = [1 / np.sqrt(2), 1 / np.sqrt(2), 0]

        cond = [1.0, 0.5, 1]
        radii = [0.4, 0.5, 1]

        pot1 = analytical_solutions.potential_3layers_surface_electrodes(
            radii, cond, anode, cathode, sensor_setup)

        anode = [1 / np.sqrt(2), -1 / np.sqrt(2), 0]
        p2 = sensor_setup
        p2[:, 1] *= -1
        pot2 = analytical_solutions.potential_3layers_surface_electrodes(
            radii, cond, anode, cathode, p2)

        np.testing.assert_array_almost_equal(pot1, pot2)

    def test_sphere_pot_symmetry_xz(self, sensor_setup):
        """ V should be symmetric to a reflection along the XY plane """
        cathode = [1, 0, 0]
        anode = [1 / np.sqrt(2), 0, 1 / np.sqrt(2)]

        cond = [1.0, 0.5, 1]
        radii = [0.4, 0.5, 1]

        pot1 = analytical_solutions.potential_3layers_surface_electrodes(
            radii, cond, anode, cathode, sensor_setup)

        anode = [1 / np.sqrt(2), 0, -1 / np.sqrt(2)]
        p2 = sensor_setup
        p2[:, 2] *= -1
        pot2 = analytical_solutions.potential_3layers_surface_electrodes(
            radii, cond, anode, cathode, p2)

        np.testing.assert_array_almost_equal(pot1, pot2)

    def test_shpere_pot_rotation_x(self, sensor_setup):
        """ V should be invariant to rotations """
        tx = np.pi / 3
        Mx = np.array(((1, 0, 0),
                       (0, np.cos(tx), -np.sin(tx)),
                       (0, np.sin(tx), np.cos(tx))), dtype=float)

        ty = np.pi / 2
        My = np.array(((np.cos(ty), 0, np.sin(ty)),
                       (0, 1, 0),
                       (-np.sin(ty), 0, np.cos(ty))), dtype=float)

        tz = np.pi / 3
        Mz = np.array(((np.cos(tz), -np.sin(tz), 0),
                       (np.sin(tz), np.cos(tz), 0),
                       (0, 0, 1)), dtype=float)

        R = Mx.dot(My.dot(Mz))

        cathode = [1, 0, 0]
        anode = [-1, 0, 0]

        cond = [1.0, 0.5, 1]
        radii = [0.4, 0.5, 1]
        pot1 = analytical_solutions.potential_3layers_surface_electrodes(
            radii, cond, anode, cathode, sensor_setup)

        cathode = R.dot(cathode)
        anode = R.dot(anode)
        sensors = np.array([R.dot(p) for p in sensor_setup])
        pot2 = analytical_solutions.potential_3layers_surface_electrodes(
            radii, cond, anode, cathode, sensors)

        np.testing.assert_allclose(pot1, pot2, rtol=1e-3, atol=1e-5)


class TestDipole:

    def test_dipole_source_at_center(self):
        points = np.array([[0, 0.7071, 0.7071]])
        v = analytical_solutions.potential_homogeneous_dipole(
            1., 1., [0, 0, 0], [0, 0, 1], points)
        v_simple = 3 / (4 * np.pi * 1e-6) * (1 / np.sqrt(2))
        assert np.abs((v[0] - v_simple) / v_simple) < 1e-3

    def test_symmetry_z_axis(self, quadrant_points):
        v = analytical_solutions.potential_homogeneous_dipole(
            1., 1., [0, 0, 0], [0, 0, 1], quadrant_points)
        np.testing.assert_almost_equal(v[quadrant_points[:, 2] > 0],
                                       -v[quadrant_points[:, 2] < 0])

    def test_tilt_dipole_x(self, quadrant_points):
        v = analytical_solutions.potential_homogeneous_dipole(
            1., 1., [0, 0, 0], [1, 0, 0], quadrant_points)
        np.testing.assert_almost_equal(v[quadrant_points[:, 0] > 0],
                                       -v[quadrant_points[:, 0] < 0])

    def test_move_dipole_xyz(self, quadrant_points):
        v = analytical_solutions.potential_homogeneous_dipole(
            1., 1., [0.5, 0.5, 0.5], [1, 0, 0], quadrant_points)
        v2 = analytical_solutions.potential_homogeneous_dipole(
            1., 1., [-0.5, -0.5, -0.5], [-1, 0, 0], -quadrant_points)
        np.testing.assert_almost_equal(v, v2)

    def test_rotate_and_move_dipole(self, quadrant_points):
        vx = analytical_solutions.potential_homogeneous_dipole(
            1., 1., [0.5, 0, 0], [1, 0, 0], quadrant_points)
        M = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        vy = analytical_solutions.potential_homogeneous_dipole(
            1., 1., [0, 0.5, 0], [0, 1, 0], M.dot(quadrant_points.T).T)
        M = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        vz = analytical_solutions.potential_homogeneous_dipole(
            1., 1., [0, 0, 0.5], [0, 0, 1], M.dot(quadrant_points.T).T)
        np.testing.assert_allclose(vx, vy)
        np.testing.assert_allclose(vx, vz)

class TestDipole3Layers:
    def test_dipole_source_at_center(self):
        points = np.array([[0, 0.7071, 0.7071]])
        v = analytical_solutions.potential_dipole_3layers(
            [0.8, 0.9, 1], 1., 1., [0, 0, 0], [0, 0, 1], points, 1)
        v_simple = 3 / (4 * np.pi * 1e-6) * (1 / np.sqrt(2))
        assert np.abs((v - v_simple) / v_simple) < 1e-3


    def test_symmetry_z_axis(self, quadrant_points):
        v = analytical_solutions.potential_dipole_3layers(
            [0.8, 0.9, 1], 1., 1., [0, 0, 0], [0, 0, 1], quadrant_points)
        np.testing.assert_almost_equal(v[quadrant_points[:, 2] > 0],
                                       -v[quadrant_points[:, 2] < 0])


    def test_tilt_dipole_x(self, quadrant_points):
        v = analytical_solutions.potential_dipole_3layers(
            [0.8, 0.9, 1], 1., 1., [0, 0, 0], [1, 0, 0], quadrant_points)
        np.testing.assert_almost_equal(v[quadrant_points[:, 0] > 0],
                                       -v[quadrant_points[:, 0] < 0])

    def test_move_dipole_xyz(self, quadrant_points):
        v = analytical_solutions.potential_dipole_3layers(
            [0.8, 0.9, 1], 1., 1., [0.3, 0.3, 0.3], [1, 0, 0], quadrant_points)
        v2 = analytical_solutions.potential_dipole_3layers(
            [0.8, 0.9, 1], 1., 1., [-0.3, -0.3, -0.3], [-1, 0, 0], -quadrant_points)
        np.testing.assert_almost_equal(v, v2)

    def test_rotate_and_move_dipole(self, quadrant_points):
        vx = analytical_solutions.potential_dipole_3layers(
            [0.8, 0.9, 1.], 1., 1., [0.5, 0, 0], [1, 0, 0], quadrant_points)
        M = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        vy = analytical_solutions.potential_dipole_3layers(
            [0.8, 0.9, 1.], 1., 1., [0, 0.5, 0], [0, 1, 0], M.dot(quadrant_points.T).T)
        M = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        vz = analytical_solutions.potential_dipole_3layers(
            [0.8, 0.9, 1.], 1., 1., [0, 0, 0.5], [0, 0, 1], M.dot(quadrant_points.T).T)
        np.testing.assert_allclose(vx, vy)
        np.testing.assert_allclose(vx, vz)

    def test_compare_homogeneous(self, quadrant_points):
        p = analytical_solutions.fibonacci_sphere(11, 2)
        v1 = analytical_solutions.potential_dipole_3layers(
                [0.8, 0.9, 2], 1., 1., [0, 0, 0.5], [1, 1, 1], p, 500)
        v2 = analytical_solutions.potential_homogeneous_dipole(
                2., 1., [0, 0, 0.5], [1, 1, 1], p)
        assert np.allclose(v1, v2)


    def test_radial_dipole(self, quadrant_points):
        ''' Tests the relationship between a 3-sphere and single sphere model, Ari et. al
        1998 '''
        b_bar = 0.5
        b_tilde = 0.317
        m_bar = 1.
        m_tilde = 1./1.536
        V_3layers = analytical_solutions.potential_dipole_3layers([.87, .92, 1.], 1., .0125,
                                                    [b_bar, 0, 0], [m_bar, 0, 0],
                                                    quadrant_points)
        V_homo = analytical_solutions.potential_homogeneous_dipole(1., 1.,
                                                     [b_tilde, 0, 0], [m_tilde, 0, 0],
                                                     quadrant_points)
        assert np.allclose(V_3layers, V_homo, rtol=5e-2)

    def test_tangential_dipole(self, quadrant_points):
        ''' Tests the relationship between a 3-sphere and single sphere model, Ari et. al
        1998 '''
        b_bar = 0.5
        b_tilde = 0.316
        m_bar = 1.
        m_tilde = 1./1.531
        V_3layers = analytical_solutions.potential_dipole_3layers([.87, .92, 1.], 1., .0125,
                                                    [b_bar, 0, 0], [0, m_bar, 0],
                                                    quadrant_points)
        V_homo = analytical_solutions.potential_homogeneous_dipole(1., 1.,
                                                     [b_tilde, 0, 0], [0, m_tilde, 0],
                                                     quadrant_points)
        assert np.allclose(V_3layers, V_homo, rtol=5e-2)




def test_reciprocity():
    cond = [1, 1, 1]
    radii = [0.7, 0.85, 1]
    electrodes = np.array(((1, 0, 0),
                           (-1, 0, 0)), dtype=float)
    dipole_pos = np.array(((0.5, 0, 0),
                           (0.5 + 1e-10, 0., 0),
                           (0.5, 1e-10, 0),
                           (0.5, 0, 1e-10),
                           (0.5 - 1e-10, 0., 0),
                           (0.5, -1e-10, 0),
                           (0.5, 0, -1e-10)), dtype=float)
    dipole_moment = np.array((1, 1, 1), dtype=float)

    pot = analytical_solutions.potential_3layers_surface_electrodes(radii, cond,
                                                      electrodes[0, :],
                                                      electrodes[1, :],
                                                      dipole_pos)
    E_field = (pot[1:4] - pot[4:7]) / 2e-13
    dipole_pot = analytical_solutions.potential_homogeneous_dipole(1., 1., dipole_pos[0],
                                                     dipole_moment, electrodes)

    np.testing.assert_allclose(
        dipole_pot[0] - dipole_pot[1], -E_field.dot(dipole_moment), rtol=1e-4)


class Test_B_field:

    def test_center_source(self, quadrant_points):
        """ A source in the center should not be visible"""
        B = analytical_solutions.B_outside_sphere(0.8, [0, 0, 0], [1, 0, 0], quadrant_points)
        np.testing.assert_allclose(B, np.zeros(B.shape))

    def test_radial_source(self, quadrant_points):
        """A radial source should not be visible """
        B = analytical_solutions.B_outside_sphere(0.8, [0.2, 0.1, 0.1], [
                                    0.2, 0.1, 0.1], quadrant_points)
        np.testing.assert_allclose(B, np.zeros(B.shape))

    def test_radial_component(self, quadrant_points):
        """ The radial component of a tangential source cant be detected """
        B1 = analytical_solutions.B_outside_sphere(
            0.8, [0, 0, 0.75], [1, 0, 0], quadrant_points)
        B2 = analytical_solutions.B_outside_sphere(
            0.8, [0, 0, 0.75], [1, 0, 1], quadrant_points)
        np.testing.assert_allclose(B1, B2)

class Test_TMS_field:
    def test_radial_source(self, quadrant_points):
        """ The field of a radial source is given by a 3.10 (Heller and Hulsteyn) """
        dipole_pos = np.array([3, 0, 0])
        dipole_moment = np.array([1, 0, 0])
        E = analytical_solutions.tms_E_field(dipole_pos, dipole_moment, 1, quadrant_points)
        E_simple = -1e-7 * np.cross(dipole_pos - quadrant_points, dipole_moment) \
                / np.linalg.norm(quadrant_points - dipole_pos, axis=1)[:, None]**3
        # Explanation for the minus: see appendix 1 in Heller and Hulsteyn
        assert np.allclose(E, E_simple)

    def test_2_dipoles(self, quadrant_points):
        dipole_pos = np.array([[3, 0, 0], [3, 0, 0]])
        dipole_moment = np.array([[1, 0, 0], [-1, 0, 0]])
        E = analytical_solutions.tms_E_field(dipole_pos, dipole_moment, 1, quadrant_points)
        assert np.allclose(E, 0)

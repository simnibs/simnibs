import math
from numbers import Integral
import warnings

import numpy as np
import numpy.typing as npt
import scipy.special as sp


def potential_3layers_surface_electrodes(radii, conductivities,
                                         anode_pos, cathode_pos,
                                         p, nbr_polynomials=50):

    """Calculates the electric potential in a 3-layered sphere caused by point-like electrodees


    Parameters
    -------------------------------
    radii: list of length = 3
        Radius of each of the 3 layers (innermost to outermost), in mm
    conductivities: list of lenght = 3
        Conductivity of the 3 layers (innermost to outermost), in S/m
    anode_pos: 3x1 ndarray
        position of the anode_pos, in mm
    cathode_pos: 3x1 ndarray
        position of cathode_pos, in mm
    p: (Nx3)ndarray
        List of positions where the poteitial should be calculated, in mm
    nbr_polynomials: int
        Number of of legendre polynomials to use

    Returns
    ------------------------------
    (Nx1) ndarray
        Values of the electric potential, in V

    Reference
    ------------------------------
    S.Rush, D.Driscol EEG electrode sensitivity--an application of reciprocity.
    """
    assert len(radii) == 3
    assert radii[0] < radii[1] and radii[1] < radii[2]
    assert len(conductivities) == 3
    assert len(anode_pos) == 3
    assert len(cathode_pos) == 3
    assert p.shape[1] == 3

    b_over_s = float(conductivities[0]) / float(conductivities[1])
    s_over_t = float(conductivities[1]) / float(conductivities[2])
    radius_brain = radii[0] * 1e-3
    radius_skull = radii[1] * 1e-3
    radius_skin = radii[2] * 1e-3

    r = np.linalg.norm(p, axis=1) * 1e-3
    theta = np.arccos(p[:, 2] * 1e-3 / r)
    phi = np.arctan2(p[:, 1], p[:, 0])

    p_r = np.vstack((r, theta, phi)).T

    cathode_pos = (np.sqrt(cathode_pos[0]**2 + cathode_pos[1]**2 + cathode_pos[2]**2) * 1e-3,
                   np.arccos(cathode_pos[2] /
                             np.sqrt(cathode_pos[0]**2 + cathode_pos[1]**2 + cathode_pos[2]**2)),
                   np.arctan2(cathode_pos[1], cathode_pos[0]))

    anode_pos = (np.sqrt(anode_pos[0]**2 + anode_pos[1]**2 + anode_pos[2]**2) * 1e-3,
                 np.arccos(anode_pos[2] /
                           np.sqrt(anode_pos[0] ** 2 + anode_pos[1]**2 + anode_pos[2]**2)),
                 np.arctan2(anode_pos[1], anode_pos[0]))

    A = lambda n: ((2 * n + 1)**3 / (2 * n)) / (((b_over_s + 1) * n + 1) * ((s_over_t + 1) * n + 1) +
                                                (b_over_s - 1) * (s_over_t - 1) * n * (n + 1) * (radius_brain / radius_skull)**(2 * n + 1) +
                                                (s_over_t - 1) * (n + 1) * ((b_over_s + 1) * n + 1) * (radius_skull / radius_skin)**(2 * n + 1) +
                                                (b_over_s - 1) * (n + 1) * ((s_over_t + 1) * (n + 1) - 1) * (radius_brain / radius_skin)**(2 * n + 1))
    # All of the bellow are modified: division by raidus_skin moved to the
    # coefficients calculations due to numerical constraints
    # THIS IS DIFFERENT FROM THE PAPER (there's a sum instead of difference)
    S = lambda n: (A(n)) * ((1 + b_over_s) * n + 1) / (2 * n + 1)
    U = lambda n: (A(n) * radius_skin) * n * (1 - b_over_s) * \
        radius_brain**(2 * n + 1) / (2 * n + 1)
    T = lambda n: (A(n) / ((2 * n + 1)**2)) *\
        (((1 + b_over_s) * n + 1) * ((1 + s_over_t) * n + 1) +
         n * (n + 1) * (1 - b_over_s) * (1 - s_over_t) * (radius_brain / radius_skull)**(2 * n + 1))
    W = lambda n: ((n * A(n) * radius_skin) / ((2 * n + 1)**2)) *\
        ((1 - s_over_t) * ((1 + b_over_s) * n + 1) * radius_skull**(2 * n + 1) +
         (1 - b_over_s) * ((1 + s_over_t) * n + s_over_t) * radius_brain**(2 * n + 1))

    brain_region = np.where(p_r[:, 0] <= radius_brain)[0]
    skull_region = np.where(
        (p_r[:, 0] > radius_brain) * (p_r[:, 0] <= radius_skull))[0]
    skin_region = np.where((p_r[:, 0] > radius_skull)
                           * (p_r[:, 0] <= radius_skin))[0]
    inside_sphere = np.where((p_r[:, 0] <= radius_skin))[0]
    outside_sphere = np.where((p_r[:, 0] > radius_skin))[0]

    cos_theta_a = np.cos(cathode_pos[1]) * np.cos(p_r[:, 1]) +\
        np.sin(cathode_pos[1]) * np.sin(p_r[:, 1]) * \
        np.cos(p_r[:, 2] - cathode_pos[2])
    cos_theta_b = np.cos(anode_pos[1]) * np.cos(p_r[:, 1]) +\
        np.sin(anode_pos[1]) * np.sin(p_r[:, 1]) * \
        np.cos(p_r[:, 2] - anode_pos[2])

    potentials = np.zeros((p.shape[0]), dtype='float64')

    coefficients = np.zeros((nbr_polynomials, p.shape[0]), dtype='float64')

    # accelerate
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        for ii in range(1, nbr_polynomials):
            n = float(ii)
            coefficients[ii, brain_region] = np.nan_to_num(
                A(n) * ((p_r[brain_region, 0] / radius_skin)**n))

            coefficients[ii, skull_region] = np.nan_to_num(S(n) * (p_r[skull_region, 0] / radius_skin)**n +
                                                           U(n) * (p_r[skull_region, 0] * radius_skin)**(-n - 1))

            coefficients[ii, skin_region] = np.nan_to_num(T(n) * (p_r[skin_region, 0] / radius_skin)**n
                                                          + W(n) * (p_r[skin_region, 0] * radius_skin)**(-n - 1))

        potentials[inside_sphere] = np.nan_to_num(
            np.polynomial.legendre.legval(cos_theta_a[inside_sphere], coefficients[:, inside_sphere], tensor=False) -
            np.polynomial.legendre.legval(cos_theta_b[inside_sphere], coefficients[:, inside_sphere], tensor=False))

    potentials *= 1.0 / (2 * np.pi * conductivities[2] * radius_skin)

    potentials[outside_sphere] = 0.0

    return potentials
    # plot_scatter(points_cart[:,0],points_cart[:,1],points_cart[:,2],potentials)


def potential_homogeneous_dipole(sphere_radius, conductivity, dipole_pos, dipole_moment,
                                 detector_positions):
    """ Calculates the surface potential generated by a dipole inside a homogeneous conducting sphere

    Calculates the surface potential generated by a dipole inside a homogeneous conducting sphere.
    Parameters
    -------------------------------------------
    sphere_radius - float
        Radius of sphere, in mm
    conductivity - float
        Conductivity of medium, in S/m
    dipole_pos - 3x1 ndarray
        Position of dipole, in mm
    dipole_moment - 3x1 ndarray
        Moment of dipole, in C.m
    detector_positions - nx3 ndarray
        Position of detectors, will be projected into the sphere surface, in mm

    Returns
    -------------------------------
    v: nx1 ndarray
       Potential at the points

    References
    ------------------------------------------------
    Dezhong Yao, Electric Potential Produced by a Dipole in a Homogeneous
    Conducting Sphere
    """
    detector_positions = np.atleast_2d(detector_positions)
    assert detector_positions.shape[1] == 3
    assert np.linalg.norm(dipole_pos) < sphere_radius

    sphere_radius = float(sphere_radius * 1e-3)
    dipole_pos = np.array(dipole_pos, dtype=float) * 1e-3
    dipole_moment = np.array(dipole_moment, dtype=float)
    detector_positions = np.array(detector_positions, dtype=float) * 1e-3

    R = sphere_radius
    r0 = np.linalg.norm(dipole_pos)
    r = np.linalg.norm(detector_positions, axis=1)
    rp = np.linalg.norm(dipole_pos - detector_positions, axis=1)

    if not np.allclose(r, R):
        warnings.warn('Some points are not in the surface!!')

    if np.isclose(r0, 0):
        cos_phi = np.zeros(len(detector_positions), dtype=float)
    else:
        cos_phi = dipole_pos.dot(detector_positions.T) / \
            (np.linalg.norm(dipole_pos) * np.linalg.norm(detector_positions, axis=1))
    second_term = 1. / (rp[:, None] * R ** 2) * \
            (detector_positions + (detector_positions * r0 * cos_phi[:, None] - R * dipole_pos) /
                                     (R + rp - r0 * cos_phi)[:, None])
    V = dipole_moment.dot((2 * (detector_positions - dipole_pos) / (rp ** 3)[:, None] +
                           second_term).T).T
    V /= 4 * np.pi * conductivity
    return V



def B_outside_sphere(sphere_radius, dipole_pos, dipole_moment, detector_positions):
    """Calculates the B field outside a sphere, does not depend on conductivity

    Parameters
    -------------------------------------------
    sphere_radius - float
        Radius of sphere
    dipole_pos - 3x1 ndarray
        Position of dipole
    dipole_moment - 3x1 ndarray
        Moment of dipole
    detector_positions - nx3 ndarray
        Position of detectors, must lie outside sphere

    Returns
    ---------------------------------------------
    nx3 ndarray
        Array with B fields in detector positions

    Notes
    ---------------------------------------------
    Dipole in SI units, positions in mm

    References
    ---------------------------------------------
        J.Savras - Basic mathematical and electromagnetic concepts of the biomagnetic inverse problem.

    """


    pos = np.array(dipole_pos, dtype=float)*1e-3
    moment = np.array(dipole_moment, dtype=float)
    detector = np.array(detector_positions, dtype=float)*1e-3

    assert np.all(np.linalg.norm(detector_positions, axis=1) >sphere_radius), "All points must be outside the sphere"

    assert np.all(np.linalg.norm(dipole_pos) <sphere_radius), "Dipole must be outside sphere"

    B = np.zeros(detector_positions.shape, dtype=float)

    for ii, r in enumerate(detector):
        norm_r = np.linalg.norm(r)

        r_0 = pos
        norm_r0 = np.linalg.norm(pos)

        a = r - r_0
        norm_a = np.linalg.norm(a)

        F = norm_a * (norm_r * norm_a + norm_r**2 - r_0.dot(r))

        grad_F = (norm_r**(-1) * norm_a**2 + norm_a**(-1) * a.dot(r) + 2 * norm_a + 2 * norm_r) * r - \
            (norm_a + 2 * norm_r + norm_a**(-1) * a.dot(r)) * r_0

        B[ii, :] = (4 * np.pi * 1e-7) /\
                   (4 * np.pi * F**2) * (F * np.cross(moment, r_0) - np.dot(np.cross(moment, r_0), r) * grad_F)

    return B


def fibonacci_sphere(nr_points, R=1):
    """ Creates N points around evenly spread through a unit sphere

    Parameters
    ------------------------------
    nr_points: int
        number of points to be spread, must be odd
    R: float
        radius of sphere

    Returns
    ------------------------------
    Nx3 ndarray
        Array with the points
    """
    assert nr_points % 2 == 1, "The number of points must be odd"
    points = []
    # The golden ratio
    phi = (1 + math.sqrt(5)) / 2.
    N = int((nr_points - 1)/2)
    for i in range(-N, N+1):
        lat = math.asin(2 * i / nr_points)
        lon = 2 * math.pi * i / phi
        x = R * math.cos(lat) * math.cos(lon)
        y = R * math.cos(lat) * math.sin(lon)
        z = R * math.sin(lat)
        points.append((x, y, z))
    return np.array(points, dtype = float)


def tms_E_field(dipole_pos, dipole_moment, didt, positions):
    """Calculates the E field in a sphere caused by external magnetic dipoles

    Everything should be in SI!!
    Independent of conductivity. See references

    Parameters
    -------------------------------------------
    dipole_pos: mx3 ndarray
        Position of dipoles, must be outside sphere
    dipole_moment: mx3 ndarray
        Moment of dipoles
    didt: float
        Variation rate of current in the coil
    positions: nx3 ndarray
        Position where fields should be calculated, must lie inside sphere

    Returns
    ---------------------------------------------
    nx3 ndarray
        Array with B fields in detector positions

    Notes
    ---------------------------------------------
    Dipole in SI units, positions in mm

    References
    ---------------------------------------------
    L. Heller and D. van Hulsteyn, Brain stimulation using electromagnetic sources:
        theoretical aspects
    """
    if dipole_pos.shape != dipole_moment.shape:
        raise ValueError('List of dipole position and moments should have the same'
                         'lengths')
    mu0_4pi = 1e-7

    E = np.zeros(positions.shape, dtype=float)
    dp = np.atleast_2d(dipole_pos)
    dm = np.atleast_2d(dipole_moment)

    r1 = positions

    for m, r2 in zip(dm, dp):
        a = r2 - r1
        norm_a = np.linalg.norm(a, axis=1)[:, None]

        norm_r1 = np.linalg.norm(r1, axis=1)[:, None]
        norm_r2 = np.linalg.norm(r2)

        r2_dot_a = np.sum(r2 * a, axis=1)[:, None]
        F = norm_a * (norm_r2 * norm_a + r2_dot_a)
        grad_F = (norm_a ** 2 / norm_r2 + 2 * norm_a + 2 * norm_r2 + r2_dot_a / norm_a)\
                * r2 - (norm_a + 2 * norm_r2 + r2_dot_a / norm_a) * r1
        E += -didt * mu0_4pi / F ** 2 * \
            (F * np.cross(r1, m) - np.cross(np.sum(m * grad_F, axis=1)[:, None] * r1, r2) )

        # Why use -didt? Take a look at the appendix 1 of the reference. It says "negative
        # time rate of change"
    return E


def lpmn(m: int, n: int, x: float | npt.ArrayLike):
    """Like scipy.special.lpmn but vectorized in x. Does not calculate
    derivative values. Does not accept negative values of m.


    PARAMETERS
    ----------
    m : int
        Order.
    n : int
        Degree.
    x : float | ndarray
        Value(s) for which to calculate function values.

    RETURNS
    -------
    pmn : (m+1, n+1[, x.shape]) array
        Function values for orders [0, 1, ..., m] and degrees [0, 1, ..., n].


    References
    ----------
    Identities taken from
        https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
    """
    assert isinstance(m, Integral) and isinstance(n, Integral)
    assert 0 <= m <= n
    x = np.array(x)

    pmn = np.zeros((m+1, n+1, *x.shape))

    # diagonal (identity 1)
    pmn[0,0] = 1
    for i in range(m):
        pmn[i+1, i+1] = -(2*i+1) * np.sqrt(1-x**2) * pmn[i, i]

    # n+1 diagonal (identity 3)
    for i in range(m if m == n else m+1):
        pmn[i, i+1] = x * (2*i+1) * pmn[i, i]

    # upper triangle (recurrence formula 1)
    for i in range(m+1):
        for j in range(i+1, n):
            pmn[i, j+1] = ((2*j+1) * x * pmn[i, j] - (j+i) * pmn[i, j-1]) / (j-i+1)

    return pmn


def potential_dipole_3layers(
    radii: npt.ArrayLike | list[float] | tuple,
    cond_brain_scalp: float,
    cond_skull: float,
    dipole_pos: npt.ArrayLike,
    dipole_moment: npt.ArrayLike,
    surface_points: npt.ArrayLike,
    nbr_polynomials: int = 100,
):
    """Calculates the electric potential in a 3-layered sphere caused by a
    dipole.

    Inputs in mm - calculations performed in SI units.

    Implements eqs. 2 and 2a from Ary et al.


    Parameters
    ----------
    radii : npt.ArrayLike | list[float] | tuple
        _description_ (in mm).
    cond_brain_scalp : float
        Conductivity of brain (inner) and scalp (outer) (in S/m).
    cond_skull : float
        Conductivity of skull (middle) (in S/m).
    dipole_pos : npt.ArrayLike
        Positions of dipoles (in mm).
    dipole_moment : npt.ArrayLike
        Dipole moment (in C x m)
    surface_points : npt.ArrayLike
        Points on the surface (in mm) to estimate the potential at (electrodes).
    nbr_polynomials : int, optional
        Number of Legendre polynomials to use (default = 100).

    Returns
    -------
    _type_
        The potentials at the surface points for each dipole.

    Notes
    -----
    Vectorized (in dipole positions and surface points) version of Guilherme's
    original implementation.


    Reference
    ---------
    Ary, J. et al. Location of sources of evoked scalp potentials: corrections
        for skull and scalp thicknesses. Biomedical Engineering 28.6 (1981).


    Raises
    ------
    ValueError
        _description_
    """
    dipole_pos = np.atleast_2d(dipole_pos).astype(float)
    dipole_moment = np.array(dipole_moment, dtype=float)
    surface_points = np.atleast_2d(surface_points).astype(float)

    assert len(radii) == 3
    assert radii[0] < radii[1] and radii[1] < radii[2]
    assert dipole_moment.ndim == 1 and dipole_moment.shape[0] == 3
    assert dipole_pos.shape[1] == 3
    assert surface_points.shape[1] == 3

    dipole_moment_norm = np.linalg.norm(dipole_moment)
    dipole_pos_norm = np.linalg.norm(dipole_pos, axis=1)

    assert np.all(dipole_pos_norm < radii[0]), "All dipoles must be inside the inner sphere"

    xi = cond_skull / cond_brain_scalp
    R = radii[2] * 1e-3
    f1 = radii[0] * 1e-3 / R
    f2 = radii[1] * 1e-3 / R
    b = dipole_pos_norm * 1e-3 / R

    if not np.allclose(np.linalg.norm(surface_points, axis=1), R * 1e3):
        raise ValueError('Some points are not on the surface!')

    at_origin = np.isclose(b, 0)
    if at_origin.any():
        r_dir = np.zeros_like(dipole_pos)
        r_dir[at_origin] = dipole_moment / dipole_moment_norm
        r_dir[~at_origin] = dipole_pos[~at_origin] / dipole_pos_norm[~at_origin, None]
    else:
        r_dir = dipole_pos / dipole_pos_norm[:, None]

    m_r = r_dir @ dipole_moment
    cos_alpha = surface_points @ r_dir.T / R * 1e-3
    t_dir = dipole_moment[None] - m_r[:, None] * r_dir

    is_radial = np.isclose(dipole_moment_norm, np.abs(r_dir @ dipole_moment))
    if is_radial.any():
        # try to set an axis in x, if the dipole is not in x
        in_x = np.isclose(np.abs(r_dir @ np.array([1, 0, 0])), 1)
        t_dir[is_radial & ~in_x] = np.array([1, 0, 0], dtype=float)

        # otherwise, set it in y
        t_dir[is_radial & in_x] = np.array([0, 1, 0], dtype=float)
        t_dir -= np.sum(r_dir * t_dir, 1)[:, None]

    t_dir /= np.linalg.norm(t_dir, axis=1, keepdims=True)
    t2_dir = np.cross(r_dir, t_dir)
    m_t = t_dir @ dipole_moment
    beta = np.arctan2(surface_points @ t2_dir.T, surface_points @ t_dir.T)
    cos_beta = np.cos(beta)

    def d(nbrs, f1, f2):
        d_n = ((nbrs + 1) * xi + nbrs) * ((nbrs * xi) / (nbrs + 1) + 1) + \
              (1 - xi) * ((nbrs + 1) * xi + nbrs) * (f1 ** (2 * nbrs + 1) - f2 ** (2 * nbrs + 1)) - \
              nbrs * (1 - xi) ** 2 * (f1 / f2) ** (2 * nbrs + 1)
        return d_n

    P = lpmn(1, nbr_polynomials, cos_alpha)
    P = P[:, 1:] # discard 0th order polynomium

    nbrs = np.arange(1, nbr_polynomials+1)

    # dimensions of the lines in the sum:
    # poly, dipole
    # poly
    # poly, elec, dipole
    potentials = np.nan_to_num(
        np.sum(
            ((2 * nbrs[:, None] + 1) / nbrs[:, None] * b[None] ** (nbrs[:, None] - 1))[:, None] * \
            ((xi * (2 * nbrs + 1) ** 2) / (d(nbrs, f1, f2) * (nbrs + 1)))[:, None, None] * \
            (nbrs[:, None, None] * m_r[None, None] * P[0] - m_t[None, None] * P[1] * cos_beta[None]),
            axis=0)
        )
        # Why should it be a minus there?

    potentials /= 4 * np.pi * cond_brain_scalp * R ** 2

    return np.squeeze(potentials)
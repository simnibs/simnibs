import numpy as np

from simnibs.optimization import ADMlib

# Direct implementation of eq S5 in
# https://www.biorxiv.org/content/10.1101/2020.05.27.120022v2
def calculate_Hprim(source_positions, source_currents, query_positions):
    H = np.zeros_like(query_positions)
    for i in range(len(query_positions)):
        r = query_positions[i] - source_positions
        dist = np.linalg.norm(r, axis=1)
        H[i] = 1e-7 * np.sum(np.cross(source_currents, r)/dist[:, None]**3, axis=0)

    return H

def test_computeHprimary():
    np.random.seed(1)
    eps = 1e-8
    current_elm_positions = np.random.rand(1000, 3)
    currents = np.random.rand(len(current_elm_positions), 3)
    observation_pos = np.random.rand(100, 3)

    fmm_Hprimary = ADMlib.computeHprimary(
        current_elm_positions.T,
        currents.T,
        observation_pos.T,
        eps
    ).T

    Hprim = calculate_Hprim(
        current_elm_positions,
        currents,
        observation_pos
    )

    assert np.allclose(fmm_Hprimary, Hprim, rtol=eps)


def test_recipcode():
    np.random.seed(2)
    current_elm_positions = np.random.rand(1000, 3)
    currents = np.random.rand(len(current_elm_positions), 3)
    observation_pos = np.random.rand(100, 3)
    observation_dir = np.random.rand(100, 3)
    observation_dir /= np.linalg.norm(observation_dir, axis=1)[:, None]

    Hprim = calculate_Hprim(
        current_elm_positions,
        currents,
        observation_pos
    )
    E_brute_force = -np.sum(Hprim*observation_dir, axis=1) # not sure about this "-"

    # create MATSIMNIBS matrices
    coil_matrices = np.zeros((4, 4, len(observation_pos)))
    for i in range(len(observation_pos)):
        A = np.random.rand(3, 3)
        A[:, 0] = observation_dir[i]
        coil_matrices[:3, :3, i] = np.linalg.qr(A)[0]
        coil_matrices[:3, :3, i] *= np.sign(
            np.dot(coil_matrices[:3, 0, i], observation_dir[i])
        )
        coil_matrices[:3, 3, i] = observation_pos[i]

    coil_dipole_pos = np.array([[0, 0, 0]])
    coil_dipole_mom = np.array([[1, 0, 0]])
    # Projem setting Nj
    E_adm = ADMlib.recipcode(
        current_elm_positions.T, currents.T,
        coil_dipole_pos.T, coil_dipole_mom.T, coil_matrices
    )
    assert np.allclose(E_brute_force, E_adm, rtol=1e-3)


def test_resample_coil_no_rotation():
    coil_dipole_pos = np.array(np.meshgrid(
        np.linspace(-1, 1, 11),
        np.linspace(-1, 1, 11),
        np.linspace(-.2, .2, 3)
    )).reshape(3, -1).T
    coil_dipole_weights = np.ones_like(coil_dipole_pos)
    #coil_dipole_weights[:, 2] = 1

    resampled_positions, resampled_weights = ADMlib.resamplecoil(
        coil_dipole_pos.T,
        coil_dipole_weights.T,
        [17, 17, 2],
        1, np.array([[0, 1, 0]]).T
    )
    resampled_positions = resampled_positions.T
    resampled_weights = resampled_weights[..., 0].T

    # Define some arbitrary function to be integrated
    def func(pos):
        return np.array([
            pos[:, 0] ** 2 + pos[:, 1],
            pos[:, 1] ** 3 + pos[:, 2],
            pos[:, 0] ** 3 + pos[:, 2],
        ]).T


    assert np.isclose(
        np.sum(coil_dipole_weights*func(coil_dipole_pos)),
        np.sum(resampled_weights*func(resampled_positions))
    )


def test_ADM():
    np.random.seed(2)
    current_elm_positions = np.random.rand(1000, 3)
    currents = np.random.rand(len(current_elm_positions), 3)
    observation_pos = np.random.rand(1, 3)
    observation_dir = np.zeros_like(observation_pos)
    observation_dir[:, 1] = 1

    #coildir = np.array([[-1, 0, 0], [0, 1, 0], [1, 0, 0]])
    coildir = np.array([[0, 1, 0]])

    # Total H field in the coil
    Htot = np.zeros_like(observation_pos)
    for sign in (-1, 1):
        for d in ([1, 0, 0], [0, 1, 0], [0, 0, 1]):
            Htot += calculate_Hprim(
                current_elm_positions,
                currents,
                observation_pos + sign * 0.1 * np.array(d)
            )

    coil_matrices = np.repeat(np.eye(4)[..., None], len(observation_pos), axis=2)
    coil_matrices[:3, 3, :] = observation_pos.T
    coil_dipole_pos = 0.1*np.array([
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1]
    ], dtype=float)
    coil_dipole_mom = np.zeros_like(coil_dipole_pos)
    coil_dipole_mom[:, 1] = 1

    E_adm = ADMlib.ADM(
        current_elm_positions.T, currents.T,
        coil_dipole_pos.T, coil_dipole_mom.T, coil_matrices,
        coildir.T
    )

    assert np.allclose(Htot[:, 1], E_adm[0])

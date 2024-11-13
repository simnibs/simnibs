import numpy as np

from .. import ADMlib

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


def test_resample_coil_rotate():
    coil_dipole_pos = np.array(np.meshgrid(
        np.linspace(-1, 1, 11),
        np.linspace(-.5, .5, 5),
        np.linspace(-.2, .2, 3)
    )).reshape(3, -1).T
    coil_dipole_weights = np.zeros_like(coil_dipole_pos)
    coil_dipole_weights[:, 1] = np.linalg.norm(coil_dipole_pos, axis=1)

    coil_dir = []
    for angle in np.linspace(np.pi/2, np.pi, 7):
        coil_dir.append([-np.sin(angle), np.cos(angle), 0])
    coil_dir = np.array(coil_dir)

    resampled_positions, resampled_weights = ADMlib.resamplecoil(
        coil_dipole_pos.T,
        coil_dipole_weights.T,
        [17, 17, 2],
        len(coil_dir), coil_dir.T
    )
    resampled_positions = resampled_positions

    # Define some arbitrary function to be integrated
    def func(pos):
        return np.array([
            pos[0] ** 2 + pos[1],
            pos[1] ** 3 + pos[2],
            pos[0] ** 3 + pos[2],
        ])

    z = np.array([0, 0, 1])
    for i, cd in enumerate(coil_dir):
        rotation_matrix = np.array([np.cross(cd, z), cd, z]).T
        positions_rotated = rotation_matrix.dot(coil_dipole_pos.T)
        weights_rotated = rotation_matrix.dot(coil_dipole_weights.T)

        assert np.isclose(
            np.sum(weights_rotated*func(positions_rotated)),
            np.sum(resampled_weights[:, :, i]*func(resampled_positions))
        )


def test_ADM_no_rotate():
    np.random.seed(2)
    coil_dipole_pos = np.array(np.meshgrid(
        np.linspace(-.5, .5, 5),
        np.linspace(-.5, .5, 5),
        np.linspace(-.2, .2, 3)
    )).reshape(3, -1).T
    coil_dipole_weights = np.zeros_like(coil_dipole_pos)
    coil_dipole_weights[:, 2] = 1

    current_elm_positions = np.random.rand(1000, 3) + 10
    currents = np.random.rand(len(current_elm_positions), 3)

    observation_pos = np.random.rand(10, 3)

    coildir = np.array([[0, 1, 0]])

    coil_matrices = np.repeat(np.eye(4)[..., None], len(observation_pos), axis=2)
    coil_matrices[:3, 3, :] = observation_pos.T

    E_adm = ADMlib.ADM(
        current_elm_positions.T, currents.T,
        coil_dipole_pos.T, coil_dipole_weights.T, coil_matrices,
        coildir.T
    )

    # Total H field in the coil
    for i, p in enumerate(observation_pos):
        H = calculate_Hprim(
            current_elm_positions,
            currents,
            coil_dipole_pos + p
        )
        E = -np.sum(H*coil_dipole_weights)
        assert np.allclose(E_adm[0, i], E)

def test_ADM_rotate():
    np.random.seed(2)
    coil_dipole_pos = np.array(np.meshgrid(
        np.linspace(-.5, .5, 5),
        np.linspace(-.5, .5, 5),
        np.linspace(-.2, .2, 3)
    )).reshape(3, -1).T
    coil_dipole_weights = np.zeros_like(coil_dipole_pos)
    coil_dipole_weights[:, 2] = 1

    current_elm_positions = np.random.rand(1000, 3) + 10
    currents = np.random.rand(len(current_elm_positions), 3)

    observation_pos = np.random.rand(10, 3)

    coil_dir = []
    for angle in np.linspace(np.pi/2, np.pi, 7):
        coil_dir.append([-np.sin(angle), np.cos(angle), 0])
    coil_dir = np.array(coil_dir)


    coil_matrices = np.repeat(np.eye(4)[..., None], len(observation_pos), axis=2)
    coil_matrices[:3, 3, :] = observation_pos.T

    E_adm = ADMlib.ADM(
        current_elm_positions.T, currents.T,
        coil_dipole_pos.T, coil_dipole_weights.T, coil_matrices,
        coil_dir.T
    )

    # Total H field in the coil
    z = np.array([0, 0, 1])
    for i, p in enumerate(observation_pos):
        for j, cd in enumerate(coil_dir):
            rotation_matrix = np.array([np.cross(cd, z), cd, z]).T
            coil_dipole_pos_moved = rotation_matrix.dot(coil_dipole_pos.T).T + p
            coil_dipole_weights_rot = rotation_matrix.dot(coil_dipole_weights.T).T
            H = calculate_Hprim(
                current_elm_positions,
                currents,
                coil_dipole_pos_moved
            )
            E = -np.sum(H*coil_dipole_weights_rot)
            assert np.allclose(E_adm[j, i], E)


def test_recipcodemag():
    np.random.seed(2)
    current_elm_positions = np.random.rand(1000, 3)

    currents = np.random.rand(3, len(current_elm_positions), 3)

    observation_pos = np.random.rand(100, 3)
    observation_dir = np.random.rand(100, 3)
    observation_dir /= np.linalg.norm(observation_dir, axis=1)[:, None]

    magnE_brute_force = np.zeros(len(observation_pos))
    for c in currents: 
        Hprim = calculate_Hprim(
            current_elm_positions, c,
            observation_pos
        )
        magnE_brute_force += np.sum(Hprim*observation_dir, axis=1)**2
    magnE_brute_force = np.sqrt(magnE_brute_force)

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
    magnE_recip = ADMlib.recipcodemag(
        current_elm_positions.T,
        currents[0, ...].T, currents[1, ...].T, currents[2, ...].T,
        coil_dipole_pos.T, coil_dipole_mom.T, coil_matrices
    )
    assert np.allclose(magnE_brute_force, magnE_recip, rtol=1e-3)


def test_ADMmag():
    np.random.seed(2)
    coil_dipole_pos = np.array(np.meshgrid(
        np.linspace(-.5, .5, 3),
        np.linspace(-.5, .5, 5),
        np.linspace(-.2, .2, 3)
    )).reshape(3, -1).T
    coil_dipole_weights = np.zeros_like(coil_dipole_pos)
    coil_dipole_weights[:, 2] = 1

    current_elm_positions = np.random.rand(1000, 3) + 10
    currents = np.random.rand(3, len(current_elm_positions), 3)

    observation_pos = np.random.rand(10, 3)

    coil_dir = []
    for angle in np.linspace(np.pi/2, np.pi, 7):
        coil_dir.append([-np.sin(angle), np.cos(angle), 0])
    coil_dir = np.array(coil_dir)


    coil_matrices = np.repeat(np.eye(4)[..., None], len(observation_pos), axis=2)
    coil_matrices[:3, 3, :] = observation_pos.T

    magnE_adm = ADMlib.ADMmag(
        current_elm_positions.T,
        currents[0, ...].T, currents[1, ...].T, currents[2, ...].T,
        coil_dipole_pos.T, coil_dipole_weights.T, coil_matrices,
        coil_dir.T
    )

    # Total H field in the coil
    z = np.array([0, 0, 1])
    for i, p in enumerate(observation_pos):
        for j, cd in enumerate(coil_dir):
            rotation_matrix = np.array([np.cross(cd, z), cd, z]).T
            coil_dipole_pos_moved = rotation_matrix.dot(coil_dipole_pos.T).T + p
            coil_dipole_weights_rot = rotation_matrix.dot(coil_dipole_weights.T).T
            magnE = 0
            for c in currents:
                H = calculate_Hprim(
                    current_elm_positions, c,
                    coil_dipole_pos_moved
                )
                magnE += np.sum(H*coil_dipole_weights_rot)**2
            magnE = np.sqrt(magnE)
            assert np.allclose(magnE_adm[j, i], magnE)



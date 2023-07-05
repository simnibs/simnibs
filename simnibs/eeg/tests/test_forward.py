from pathlib import Path
from tempfile import TemporaryDirectory
import shutil

import numpy as np
import pytest

from simnibs import SIMNIBSDIR
from simnibs.utils.cond_utils import standard_cond
from simnibs.simulation.analytical_solutions import potential_dipole_3layers
from simnibs.eeg import forward
from simnibs.utils.file_finder import SubjectFiles
from simnibs.utils.csv_reader import write_csv_positions
from simnibs.mesh_tools import mesh_io
from simnibs.segmentation.brain_surface import fibonacci_sphere, fibonacci_sphere_points

from simnibs.simulation.tests.test_fem import rdm, mag

@pytest.fixture
def sphere3_msh():
    return mesh_io.read_msh(
        Path(SIMNIBSDIR) / "_internal_resources" / "testing_files" / "sphere3.msh"
    )

# False gives large errors with this coarse mesh
@pytest.mark.parametrize("point_electrodes", [True])
@pytest.mark.parametrize("solver", ["petsc", "pardiso"])
def test_prepare_forward(point_electrodes, solver, sphere3_msh):

    s = standard_cond()
    # Special tags for sphere3_msh
    s[2].value = 0.3  # gray matter
    s[3].value = 0.006  # bone (average)
    s[4].value = 0.3  # skin

    # we cannot easily test subsampling with this faked structure
    subsampling = None
    # inner compartment is 3 in sphere3_msh
    init_kwargs = dict(interpolation_tissue=[3], cond=s)
    if solver == "pardiso":
        init_kwargs["solver_options"] = "pardiso"

    radii = [85, 90, 95]
    for tag, radius in zip([1003, 1004, 1005], radii):
        tris = np.unique(
            sphere3_msh.elm.node_number_list[sphere3_msh.elm.tag1 == tag, :3] - 1
        )
        assert np.allclose(
            np.array(3 * [-radius]), sphere3_msh.nodes.node_coord[tris].min(0)
        )
        assert np.allclose(
            np.array(3 * [radius]), sphere3_msh.nodes.node_coord[tris].max(0)
        )

    n = 50
    p_pos = fibonacci_sphere_points(n) * radii[-1]
    p_type = n * ["Electrode"]
    p_name = np.arange(n).astype(str)

    # move sources quite far away from boundaries due to low mesh resolution
    source_space = dict(
        lh=dict(zip(("points", "tris"), fibonacci_sphere(20, radius=75))),
        rh=dict(zip(("points", "tris"), fibonacci_sphere(20, radius=70))),
    )

    dip_pos = np.concatenate([v["points"] for v in source_space.values()])
    # Only calculate and compare dipoles oriented along x axis
    dip_mom = np.repeat(np.array([[1, 0, 0]]), len(dip_pos), axis=0)
    cond_brain_scalp = s[2].value
    cond_skull = s[3].value

    sol = np.zeros((len(dip_pos), n))
    for i, (dp, dm) in enumerate(zip(dip_pos, dip_mom)):
        sol[i] = potential_dipole_3layers(
            radii, cond_brain_scalp, cond_skull, dp, dm, p_pos, nbr_polynomials=100
        )
    sol = sol.T

    with TemporaryDirectory(prefix="m2m_") as d:
        m2m_dir = Path(d)
        m2m = SubjectFiles(subpath=m2m_dir)
        Path(m2m.surface_folder).mkdir()
        fem_dir = m2m_dir / "fem"
        f_montage = m2m_dir / "montage.csv"

        shutil.copyfile(sphere3_msh.fn, m2m.fnamehead)

        write_csv_positions(f_montage, p_type, p_pos, p_name)

        for k, v in source_space.items():
            m = mesh_io.Msh(mesh_io.Nodes(v["points"]), mesh_io.Elements(v["tris"] + 1))
            mesh_io.write_gifti_surface(m, m2m.get_surface(k, "central"))

        forward.compute_tdcs_leadfield(
            m2m_dir, fem_dir, str(f_montage), subsampling, point_electrodes, init_kwargs
        )
        f_leadfield = fem_dir / f"{m2m.subid}_leadfield_{f_montage.stem}.hdf5"

        for apply_average_proj in (False, True):
            lf = forward.prepare_forward(f_leadfield, apply_average_proj)

            gain = lf["data"][..., 0] # dip_mom for x axis

            if apply_average_proj:
                sol_ = sol - sol.mean(0)
            else:
                sol_ = sol[1:] - sol[0]

            assert np.all(rdm(sol_, gain) < 0.1), f"Unacceptable error using apply_average_proj={apply_average_proj}"
            assert np.all(mag(sol_, gain) < 0.1), f"Unacceptable error using apply_average_proj={apply_average_proj}"

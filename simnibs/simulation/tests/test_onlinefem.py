import os
import pytest
import numpy as np

from simnibs.optimization.tes_flex_optimization.electrode_layout import ElectrodeArrayPair
from simnibs.simulation import analytical_solutions, onlinefem
from simnibs.simulation import tms_coil
import sys
from simnibs.mesh_tools import mesh_io
from simnibs import SIMNIBSDIR
from simnibs.simulation.onlinefem import OnlineFEM, FemTargetPointCloud

@pytest.fixture
def sphere3_msh():
    fn = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn)


@pytest.fixture
def tms_sphere(sphere3_msh):
    m = sphere3_msh.crop_mesh(elm_type=4)
    dipole_pos = np.array([0., 0., 300])
    dipole_moment = np.array([1., 0., 0.])
    didt = 1e6
    r = (m.nodes.node_coord - dipole_pos) * 1e-3
    dAdt = 1e-7 * didt * np.cross(dipole_moment, r) / (np.linalg.norm(r, axis=1)[:, None] ** 3)
    dAdt = mesh_io.NodeData(dAdt, mesh=m)
    dAdt.field_name = 'dAdt'
    dAdt.mesh = m
    pos = m.elements_baricenters().value
    E_analytical = analytical_solutions.tms_E_field(dipole_pos * 1e-3,
                                                    dipole_moment, didt,
                                                    pos * 1e-3)
    cond = mesh_io.ElementData(np.ones(m.elm.nr))
    cond.mesh = m
    stimulator = tms_coil.tms_stimulator.TmsStimulator("Example Stimulator", "Example Stimulator Brand", 100e6)
    dipole_elm = tms_coil.tms_coil_element.DipoleElements(stimulator=stimulator, points=dipole_pos[None], values=dipole_moment[None])
    coil = tms_coil.tms_coil.TmsCoil([dipole_elm])
    return m, cond, dAdt, E_analytical, coil


@pytest.mark.parametrize(
    ["solver_type",],
    [
        pytest.param(
            "petsc_pardiso",
            marks=pytest.mark.skipif(
                sys.platform == "darwin",
                reason="PETSc is not built with Intel MKL on macos."
            ),
        ),
        pytest.param(
            "pardiso",
            marks=pytest.mark.skipif(
                sys.platform == "darwin",
                reason="MKL Pardiso is available on macos."
            ),
        ),
        pytest.param(
            "mumps",
            ),
        pytest.param(
            "hypre",
            ),
    ]
)

def rdm(a, b):
    return np.linalg.norm(a / np.linalg.norm(a) -
                          b / np.linalg.norm(b))


def mag(a, b):
    return np.abs(np.log(np.linalg.norm(a) / np.linalg.norm(b)))


class TestOnlineFEM:
    def test_tes(self):
        # mesh
        fn_mesh = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
        mesh = mesh_io.read_msh(fn_mesh)

        # roi
        domain = 3
        roi = FemTargetPointCloud(mesh,
                                  mesh.elements_baricenters()[mesh.elm.tag1 == domain]
                               )

        # electrode
        electrode_i = ElectrodeArrayPair()
        electrode_i.center = [[0, 0]]                       # electrode center in reference electrode space (x-y plane)
        electrode_i.radius = [10]                           # radius of electrodes
        electrode_i.dirichlet_correction_detailed = False   # node wise dirichlet correction
        electrode_i.current = [1, -1]               # electrode currents
        electrode_i._prepare()

        node_indices = mesh.elm.node_number_list[np.where(mesh.elm.tag1 == 1005)[0][0]][:3]
        electrode_i._electrode_arrays[0].electrodes[0].n_nodes = 1
        electrode_i._electrode_arrays[1].electrodes[0].n_nodes = 1
        electrode_i._electrode_arrays[0].electrodes[0].node_current = np.array([1])
        electrode_i._electrode_arrays[1].electrodes[0].node_current = np.array([-1])
        electrode_i._electrode_arrays[0].electrodes[0].node_idx = np.array([node_indices[0]])
        electrode_i._electrode_arrays[1].electrodes[0].node_idx = np.array([node_indices[1]])

        # online fem
        ofem = OnlineFEM(mesh=mesh,
                         electrode=[electrode_i],
                         method="TES",
                         roi=[roi],
                         anisotropy_type="scalar",
                         solver_options="pardiso",
                         fn_logger=None,
                         useElements=True,
                         dataType=[1],
                         dirichlet_node=1)

        # e-field
        e = ofem.update_field(electrode=[electrode_i])[0][0]

        assert (e.shape == np.array([np.sum(mesh.elm.tag1 == 3), 3])).all()
        assert (e != 0).all()

    def test_tes_dirichlet(self):
        # mesh
        fn_mesh = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
        mesh = mesh_io.read_msh(fn_mesh)

        # roi
        domain = 3
        roi = FemTargetPointCloud(mesh,
                                  mesh.elements_baricenters()[mesh.elm.tag1 == domain]
                               )

        # electrode
        electrode_i = ElectrodeArrayPair()
        electrode_i.center = [[0, 0]]                       # electrode center in reference electrode space (x-y plane)
        electrode_i.radius = [10]                           # radius of electrodes
        electrode_i.dirichlet_correction_detailed = False   # node wise dirichlet correction
        electrode_i.current = [1, -1]               # electrode currents
        electrode_i._prepare()

        node_indices = np.unique(mesh.elm.node_number_list[np.where(mesh.elm.tag1 == 1005)[0]][:, :3].flatten())[:20]
        electrode_i._electrode_arrays[0].electrodes[0].n_nodes = 10
        electrode_i._electrode_arrays[1].electrodes[0].n_nodes = 10
        electrode_i._electrode_arrays[0].electrodes[0].node_current = np.ones(10)
        electrode_i._electrode_arrays[1].electrodes[0].node_current = -np.ones(10)
        electrode_i._electrode_arrays[0].electrodes[0].node_idx = node_indices[:10]
        electrode_i._electrode_arrays[1].electrodes[0].node_idx = node_indices[10:]
        electrode_i._electrode_arrays[0].electrodes[0].node_coords = np.random.rand(10, 3)
        electrode_i._electrode_arrays[1].electrodes[0].node_coords = np.random.rand(10, 3)
        electrode_i._electrode_arrays[0].electrodes[0].node_area = np.ones(10)
        electrode_i._electrode_arrays[1].electrodes[0].node_area = np.ones(10)

        # online fem
        ofem = OnlineFEM(mesh=mesh,
                         electrode=[electrode_i],
                         method="TES",
                         roi=[roi],
                         anisotropy_type="scalar",
                         solver_options="pardiso",
                         fn_logger=None,
                         useElements=True,
                         dataType=[1],
                         dirichlet_node=1)

        # e-field
        e = ofem.update_field(electrode=[electrode_i], dirichlet_correction=True)[0][0]

        assert (e.shape == np.array([np.sum(mesh.elm.tag1 == 3), 3])).all()
        assert (e != 0).all()

    def test_tms(self):
        # mesh
        fn_mesh = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
        mesh = mesh_io.read_msh(fn_mesh)

        # roi
        domain = 3
        roi = FemTargetPointCloud(mesh,
                                  mesh.elements_baricenters()[mesh.elm.tag1 == domain]
                               )

        # coil
        fn_coil = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'testcoil.nii.gz')
        matsimnibs = np.eye(4)
        matsimnibs[2, 3] = np.max(mesh.nodes.node_coord[:, 2]) + 1e-6

        # online fem
        ofem = OnlineFEM(mesh=mesh,
                         fn_coil=fn_coil,
                         method="TMS",
                         roi=[roi],
                         anisotropy_type="scalar",
                         solver_options="pardiso",
                         fn_logger=None,
                         useElements=True,
                         dataType=[1],
                         dirichlet_node=1)

        # e-field
        e = ofem.update_field(matsimnibs=matsimnibs)[0][0]

        assert (e.shape == np.array([np.sum(mesh.elm.tag1 == 3), 3])).all()
        assert (e != 0).all()

    def test_tms_sphere(self, tms_sphere, solver_type):
        # get analytical solution and dipole coil object
        m, cond, dAdt, E_analytical, coil = tms_sphere

        # create the ROI
        center_points = m.elements_baricenters().value
        center_points = np.concatenate((center_points, np.array((0, 0, 100))), axis=0)
        point_cloud = onlinefem.FemTargetPointCloud(m, center_points, nearest_neighbor=True)
        
        #set dummy mesh filename, not used but needed
        m.fn = 'temp.msh'
        
        #prepare and setup OnlineFEM
        ofem = onlinefem.OnlineFEM(m, 'TMS', roi=[point_cloud], coil=coil, useElements=False, solver_options=solver_type, cond=cond)
        
        # v = fem.tms_dadt(m, cond, dAdt)
        # E = -v.gradient().value * 1e3 - dAdt.node_data2elm_data().value
        
        #caluculate vector E-field
        ofem.dataType = [1]
        #Solve the FEM
        E = ofem.update_field(matsimnibs=np.identity(4), didt=1e6)[0][0]        
        
        # print(f'RDM: {rdm(E, E_analytical)}')
        # print(f'MAG: {mag(E, E_analytical)}')
        
        assert rdm(E, E_analytical) < .2
        assert np.abs(mag(E, E_analytical)) < np.log(1.1)
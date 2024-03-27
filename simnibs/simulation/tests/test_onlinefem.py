import os
import pytest
import simnibs
import numpy as np
from ... mesh_tools import mesh_io
from ... import SIMNIBSDIR
from .. onlinefem import OnlineFEM
from ..onlinefem import FemTargetPointCloud


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
        electrode_i = simnibs.ElectrodeInitializer()
        electrode_i.type = "ElectrodeArrayPair"             # Pair of TES electrodes
        electrode_i.center = [[0, 0]]                       # electrode center in reference electrode space (x-y plane)
        electrode_i.radius = [10]                           # radius of electrodes
        electrode_i.dirichlet_correction_detailed = False   # node wise dirichlet correction
        electrode_i.current = [1, -1]               # electrode currents
        electrode = electrode_i.initialize()

        node_indices = mesh.elm.node_number_list[np.where(mesh.elm.tag1 == 1005)[0][0]][:3]
        electrode.electrode_arrays[0].electrodes[0].n_nodes = 1
        electrode.electrode_arrays[1].electrodes[0].n_nodes = 1
        electrode.electrode_arrays[0].electrodes[0].node_current = np.array([1])
        electrode.electrode_arrays[1].electrodes[0].node_current = np.array([-1])
        electrode.electrode_arrays[0].electrodes[0].node_idx = np.array([node_indices[0]])
        electrode.electrode_arrays[1].electrodes[0].node_idx = np.array([node_indices[1]])

        # online fem
        ofem = OnlineFEM(mesh=mesh,
                         electrode=[electrode],
                         method="TES",
                         roi=[roi],
                         anisotropy_type="scalar",
                         solver_options="pardiso",
                         fn_results=None,
                         useElements=True,
                         dataType=[1],
                         dirichlet_node=1)

        # e-field
        e = ofem.update_field(electrode=[electrode])[0][0]

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
        electrode_i = simnibs.ElectrodeInitializer()
        electrode_i.type = "ElectrodeArrayPair"             # Pair of TES electrodes
        electrode_i.center = [[0, 0]]                       # electrode center in reference electrode space (x-y plane)
        electrode_i.radius = [10]                           # radius of electrodes
        electrode_i.dirichlet_correction_detailed = False   # node wise dirichlet correction
        electrode_i.current = [1, -1]               # electrode currents
        electrode = electrode_i.initialize()

        node_indices = np.unique(mesh.elm.node_number_list[np.where(mesh.elm.tag1 == 1005)[0]][:, :3].flatten())[:20]
        electrode.electrode_arrays[0].electrodes[0].n_nodes = 10
        electrode.electrode_arrays[1].electrodes[0].n_nodes = 10
        electrode.electrode_arrays[0].electrodes[0].node_current = np.ones(10)
        electrode.electrode_arrays[1].electrodes[0].node_current = -np.ones(10)
        electrode.electrode_arrays[0].electrodes[0].node_idx = node_indices[:10]
        electrode.electrode_arrays[1].electrodes[0].node_idx = node_indices[10:]
        electrode.electrode_arrays[0].electrodes[0].node_coords = np.random.rand(10, 3)
        electrode.electrode_arrays[1].electrodes[0].node_coords = np.random.rand(10, 3)
        electrode.electrode_arrays[0].electrodes[0].node_area = np.ones(10)
        electrode.electrode_arrays[1].electrodes[0].node_area = np.ones(10)

        # online fem
        ofem = OnlineFEM(mesh=mesh,
                         electrode=[electrode],
                         method="TES",
                         roi=[roi],
                         anisotropy_type="scalar",
                         solver_options="pardiso",
                         fn_results=None,
                         useElements=True,
                         dataType=[1],
                         dirichlet_node=1)

        # e-field
        e = ofem.update_field(electrode=[electrode], dirichlet_correction=True)[0][0]

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
                         fn_results=None,
                         useElements=True,
                         dataType=[1],
                         dirichlet_node=1)

        # e-field
        e = ofem.update_field(matsimnibs=matsimnibs)[0][0]

        assert (e.shape == np.array([np.sum(mesh.elm.tag1 == 3), 3])).all()
        assert (e != 0).all()

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


@pytest.fixture
def cube_msh():
    fn = os.path.join(SIMNIBSDIR, '_internal_resources', 'testing_files', 'cube_w_electrodes.msh')
    m = mesh_io.read_msh(fn)
    m = m.crop_mesh(tags = [5, 1005])
    return m


def rdm(a, b):
    return np.linalg.norm(a / np.linalg.norm(a) -
                          b / np.linalg.norm(b))


def mag(a, b):
    return np.abs(np.log(np.linalg.norm(a) / np.linalg.norm(b)))


class TestOnlineFEM:
    
    #loop over 4 solver types if available
    @pytest.mark.parametrize(
        "solver_type", 
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
                    reason="MKL Pardiso is not available on macos."
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
    #loop over nearest neighbor or SPR interpolation
    @pytest.mark.parametrize(
        "nearest_neighbor",
        [
            pytest.param(
                False),
            pytest.param(
                True)
        ]
    )
    #loop over nearest neighbor or SPR interpolation
    @pytest.mark.parametrize(
        "dirichlet_correction",
        [
            pytest.param(
                False),
            pytest.param(
                True)
        ]
    )
    def test_tes_cube(self, cube_msh, solver_type, nearest_neighbor, dirichlet_correction):
        m = cube_msh        
        cond = mesh_io.ElementData(1.0*np.ones(m.elm.nr))
        cond.mesh = m
        
        idx_dirinode = np.argmin(np.sum(m.nodes.node_coord**2,axis=1))
        idx_node_1 = np.where(m.nodes.node_coord[:,1] >= 50.0)[0] # top plate
        idx_node_2 = np.where(m.nodes.node_coord[:,1] <= -50.0)[0] # bottom plate
        
        n_areas = m.nodes_areas().value
        n_current1 =  1.0 * n_areas[idx_node_1] / np.sum(n_areas[idx_node_1])
        n_current2 = -1.0 * n_areas[idx_node_2] / np.sum(n_areas[idx_node_2])
        
        # electrode
        electrode_i = ElectrodeArrayPair()
        electrode_i.center = [[0, 0]]                       
        electrode_i.radius = [10]                          
        electrode_i.dirichlet_correction = True
        electrode_i.dirichlet_correction_detailed = True  
        electrode_i.current = [1, -1]             
        electrode_i._prepare()
        
        electrode_i._electrode_arrays[0].electrodes[0].n_nodes = len(idx_node_1)
        electrode_i._electrode_arrays[1].electrodes[0].n_nodes = len(idx_node_2)
        electrode_i._electrode_arrays[0].electrodes[0].node_current = n_current1
        electrode_i._electrode_arrays[1].electrodes[0].node_current = n_current2
        electrode_i._electrode_arrays[0].electrodes[0].node_idx = idx_node_1
        electrode_i._electrode_arrays[1].electrodes[0].node_idx = idx_node_2
        electrode_i._electrode_arrays[0].electrodes[0].node_coords = m.nodes.node_coord[idx_node_1]
        electrode_i._electrode_arrays[1].electrodes[0].node_coords = m.nodes.node_coord[idx_node_2]
        electrode_i._electrode_arrays[0].electrodes[0].node_area = n_areas[idx_node_1]
        electrode_i._electrode_arrays[1].electrodes[0].node_area = n_areas[idx_node_2]
        
        # roi
        idx_roi_elm = np.where(m.elm.tag1 == 5)[0][0:-1:100]
        roi = FemTargetPointCloud(m, m.elements_baricenters().value[idx_roi_elm], nearest_neighbor=nearest_neighbor)
                
        # online fem
        ofem = OnlineFEM(mesh=m,
                         electrode=[electrode_i],
                         method="TES",
                         roi=[roi],
                         anisotropy_type="scalar",
                         solver_options="pardiso",
                         fn_logger=None,
                         dataType=[1],
                         dirichlet_node=idx_dirinode+1,
                         cond=cond)
        
        # e-field
        e = ofem.update_field(electrode=[electrode_i], dirichlet_correction=dirichlet_correction )[0][0]
        
        # check results
        avg_v1 =  np.sum(ofem.v[idx_node_1] * n_areas[idx_node_1]) / np.sum(n_areas[idx_node_1])
        avg_v2 =  np.sum(ofem.v[idx_node_2] * n_areas[idx_node_2]) / np.sum(n_areas[idx_node_2])
        sol = (m.nodes.node_coord[:, 1] - 25) / 10
        sol_e = np.zeros_like(e)
        sol_e[:,1] = -100
        
        if dirichlet_correction==True and nearest_neighbor == False:
            limits = [0.005, 0.05, 0.05, 1e-04, 0.03, 0.001]
        else:
            limits = [0.02, 0.1, 0.1, 1e-04, 0.06, 0.01]
        
        assert e.shape == (len(idx_roi_elm), 3)
        assert abs(avg_v1-avg_v2 - 10) < limits[0]
        assert np.std(ofem.v[idx_node_1]) < limits[1]
        assert np.std(ofem.v[idx_node_2]) < limits[2]
        assert 1-np.corrcoef(ofem.v, sol)[0,1] < limits[3]
        assert rdm(e, sol_e) < limits[4]
        assert mag(e, sol_e) < limits[5]
        
    
    #loop over 4 solver types if available
    @pytest.mark.parametrize(
        "solver_type", 
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
                    reason="MKL Pardiso is not available on macos."
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
    #loop over nearest neighbor or SPR interpolation
    @pytest.mark.parametrize(
        "dirichlet_correction_detailed",
        [
            pytest.param(
                False),
            pytest.param(
                True)
        ]
    )
    def test_tes_cube_dirichlet_channel(self, cube_msh, solver_type, dirichlet_correction_detailed):
        dirichlet_correction=True
    
        m = cube_msh
    
        cond = mesh_io.ElementData(1.0*np.ones(m.elm.nr))
        cond.mesh = m
        elm_bari = m.elements_baricenters().value
        cond.value[elm_bari[:,0]<0] = 2. # divide into 2 halves with different conductivities
    
        idx_dirinode = np.argmin(np.sum(m.nodes.node_coord**2,axis=1))
        idx_node_1_l = np.where((m.nodes.node_coord[:,1] >= 50.0) * (m.nodes.node_coord[:,0] < 0.) )[0] # top plate left
        idx_node_1_r = np.where((m.nodes.node_coord[:,1] >= 50.0) * (m.nodes.node_coord[:,0] >= 0.) )[0] # top plate right
    
        idx_node_2_l = np.where((m.nodes.node_coord[:,1] <= -50.0) * (m.nodes.node_coord[:,0] < 0.) )[0] # bottom plate left
        idx_node_2_r = np.where((m.nodes.node_coord[:,1] <= -50.0) * (m.nodes.node_coord[:,0] >= 0.))[0] # bottom plate right
    
        n_areas = m.nodes_areas().value
        n_current1_l =  0.5 * n_areas[idx_node_1_l] / np.sum(n_areas[idx_node_1_l])
        n_current1_r =  0.5 * n_areas[idx_node_1_r] / np.sum(n_areas[idx_node_1_r])
        n_current2_l = -0.5 * n_areas[idx_node_2_l] / np.sum(n_areas[idx_node_2_l])
        n_current2_r = -0.5 * n_areas[idx_node_2_r] / np.sum(n_areas[idx_node_2_r])
    
        # electrode
        electrode_i = ElectrodeArrayPair()
        electrode_i.center = [[0, 0], [0, 0]]                       
        electrode_i.radius = [10, 10]                          
        electrode_i.dirichlet_correction = True
        electrode_i.dirichlet_correction_detailed = dirichlet_correction_detailed  
        electrode_i.current = [0.5, 0.5, -0.5, -0.5]             
        electrode_i._prepare()
    
        electrode_i._electrode_arrays[0].electrodes[0].n_nodes = len(idx_node_1_l)
        electrode_i._electrode_arrays[0].electrodes[0].node_current = n_current1_l
        electrode_i._electrode_arrays[0].electrodes[0].node_idx = idx_node_1_l
        electrode_i._electrode_arrays[0].electrodes[0].node_coords = m.nodes.node_coord[idx_node_1_l]
        electrode_i._electrode_arrays[0].electrodes[0].node_area = n_areas[idx_node_1_l]
    
        electrode_i._electrode_arrays[0].electrodes[1].n_nodes = len(idx_node_1_r)
        electrode_i._electrode_arrays[0].electrodes[1].node_current = n_current1_r
        electrode_i._electrode_arrays[0].electrodes[1].node_idx = idx_node_1_r
        electrode_i._electrode_arrays[0].electrodes[1].node_coords = m.nodes.node_coord[idx_node_1_r]
        electrode_i._electrode_arrays[0].electrodes[1].node_area = n_areas[idx_node_1_r]
    
        electrode_i._electrode_arrays[1].electrodes[0].n_nodes = len(idx_node_2_l)
        electrode_i._electrode_arrays[1].electrodes[0].node_current = n_current2_l
        electrode_i._electrode_arrays[1].electrodes[0].node_idx = idx_node_2_l
        electrode_i._electrode_arrays[1].electrodes[0].node_coords = m.nodes.node_coord[idx_node_2_l]
        electrode_i._electrode_arrays[1].electrodes[0].node_area = n_areas[idx_node_2_l]
    
        electrode_i._electrode_arrays[1].electrodes[1].n_nodes = len(idx_node_2_r)
        electrode_i._electrode_arrays[1].electrodes[1].node_current = n_current2_r
        electrode_i._electrode_arrays[1].electrodes[1].node_idx = idx_node_2_r
        electrode_i._electrode_arrays[1].electrodes[1].node_coords = m.nodes.node_coord[idx_node_2_r]
        electrode_i._electrode_arrays[1].electrodes[1].node_area = n_areas[idx_node_2_r]
    
        # roi
        idx_roi_elm = np.where(m.elm.tag1 == 5)[0][0:-1:100]
        roi = FemTargetPointCloud(m, m.elements_baricenters().value[idx_roi_elm], nearest_neighbor=False)
                
        # online fem
        ofem = OnlineFEM(mesh=m,
                         electrode=[electrode_i],
                         method="TES",
                         roi=[roi],
                         anisotropy_type="scalar",
                         solver_options="pardiso",
                         fn_logger=None,
                         dataType=[1],
                         dirichlet_node=idx_dirinode+1,
                         cond=cond)
    
        # e-field
        e = ofem.update_field(electrode=[electrode_i], dirichlet_correction=dirichlet_correction )[0][0]
    
        # check results
        avg_v1_l =  np.sum(ofem.v[idx_node_1_l] * n_areas[idx_node_1_l]) / np.sum(n_areas[idx_node_1_l])
        avg_v1_r =  np.sum(ofem.v[idx_node_1_r] * n_areas[idx_node_1_r]) / np.sum(n_areas[idx_node_1_r])
        avg_v2_l =  np.sum(ofem.v[idx_node_2_l] * n_areas[idx_node_2_l]) / np.sum(n_areas[idx_node_2_l])
        avg_v2_r =  np.sum(ofem.v[idx_node_2_r] * n_areas[idx_node_2_r]) / np.sum(n_areas[idx_node_2_r])
    
        assert e.shape == (len(idx_roi_elm), 3)
        
        assert abs(avg_v1_l-avg_v2_l-6.66666) < 0.05
        assert abs(avg_v1_r-avg_v2_r-6.66666) < 0.05
        assert abs(avg_v1_l-avg_v1_r) < 0.05
        assert abs(avg_v2_l-avg_v2_r) < 0.05
    
        assert(np.std(ofem.v[idx_node_1_l])) < 0.06
        assert(np.std(ofem.v[idx_node_1_r])) < 0.06
        assert(np.std(ofem.v[idx_node_2_l])) < 0.06
        assert(np.std(ofem.v[idx_node_2_r])) < 0.06
    

    #loop over 4 solver types if available
    @pytest.mark.parametrize(
        "solver_type", 
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
                    reason="MKL Pardiso is not available on macos."
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
    #loop over filling outside values with nearest neighbor or zero
    @pytest.mark.parametrize(
        "fill",
        [
            pytest.param(
                False),
            pytest.param(
                True)
        ]
    )
    #loop over nearest neighbor or SPR interpolation
    @pytest.mark.parametrize(
        "nearest",
        [
            pytest.param(
                False),
            pytest.param(
                True)
        ]
    )
    #loop over dadt in elements or nodes
    @pytest.mark.parametrize(
        "useElements",
        [
            pytest.param(
                False),
            pytest.param(
                True)
        ]
    )
    def test_tms_sphere(self, tms_sphere, solver_type, fill, nearest, useElements):
        # get analytical solution and dipole coil object
        m, cond, dAdt, E_analytical, coil = tms_sphere

        # create the ROI
        center_points = m.elements_baricenters().value
        out_point = np.array((0, 0, 100))[None]
        center_points = np.concatenate((center_points, out_point), axis=0)
        point_cloud = onlinefem.FemTargetPointCloud(m, center_points, nearest_neighbor=nearest, fill_nearest=fill)
                
        #prepare and setup OnlineFEM
        ofem = onlinefem.OnlineFEM(m, 'TMS', roi=[point_cloud], coil=coil, useElements=useElements, solver_options=solver_type, cond=cond)
        
        #calculate vector E-field
        ofem.dataType = [1]
        #Solve the FEM
        E = ofem.update_field(matsimnibs=np.identity(4), didt=1e6)[0][0]        
        
        assert rdm(E[:-1,:], E_analytical) < .2
        assert np.abs(mag(E[:-1,:], E_analytical)) < np.log(1.1)
        
        if fill:
            #find nearest point to the extra one outside
            nearest_idx = np.argmin(np.sqrt(((
                m.elements_baricenters().value - out_point)**2).sum(axis=1)))
            if nearest:
                assert np.all(E[-1, :] == E[nearest_idx])
            else:
                assert ofem.roi[0].sF[-1, :].sum() == 1
                assert ofem.roi[0].sF[-1, nearest_idx] == 1
        else:
            assert np.all(E[-1, :] == 0)
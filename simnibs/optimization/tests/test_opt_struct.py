import os
import csv
from mock import patch, MagicMock
import tempfile

import pytest
import numpy as np
import nibabel
import h5py
import scipy.io

from ... import SIMNIBSDIR
from ...mesh_tools import mesh_io
from ...simulation import sim_struct
from ...simulation import analytical_solutions
from .. import opt_struct


@pytest.fixture()
def sphere_surf():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn).crop_mesh([1003, 1004])


@pytest.fixture()
def sphere_vol():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn).crop_mesh([4, 5])

@pytest.fixture()
def sphere_msh():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files',  'sphere3.msh')
    return mesh_io.read_msh(fn)

@pytest.fixture
def sphere_elec():
    fn = os.path.join(
        SIMNIBSDIR, '_internal_resources', 'testing_files',
        'sphere_w_electrodes.msh')
    return mesh_io.read_msh(fn).crop_mesh([1005, 1100, 1101])

@pytest.fixture()
def leadfield_surf(sphere_surf):
    np.random.seed(0)
    return np.random.random((4, sphere_surf.nodes.nr, 3))

@pytest.fixture()
def leadfield_vol(sphere_vol):
    np.random.seed(0)
    return np.random.random((4, sphere_vol.elm.nr, 3))


@pytest.fixture()
def leadfield_elec(sphere_vol):
    np.random.seed(0)
    return np.random.random((1, sphere_vol.elm.nr, 3))


@pytest.fixture()
def fn_surf(sphere_surf, leadfield_surf):
    fn_leadfield = 'tmp_surf_leadfied.hdf5'
    if os.path.isfile(fn_leadfield):
        os.remove(fn_leadfield)
    sphere_surf.write_hdf5(fn_leadfield, 'mesh_leadfield')
    dset = '/mesh_leadfield/leadfields/tdcs_leadfield'
    with h5py.File(fn_leadfield, 'a') as f:
        f.create_dataset(dset, data=leadfield_surf)
    yield fn_leadfield
    os.remove(fn_leadfield)


@pytest.fixture()
def fn_vol(sphere_vol, leadfield_vol):
    fn_leadfield = 'tmp_vol_leadfied.hdf5'
    if os.path.isfile(fn_leadfield):
        os.remove(fn_leadfield)
    sphere_vol.write_hdf5(fn_leadfield, 'mesh_leadfield')
    dset = '/mesh_leadfield/leadfields/tdcs_leadfield'
    with h5py.File(fn_leadfield, 'a') as f:
        f.create_dataset(dset, data=leadfield_vol)
    yield fn_leadfield
    os.remove(fn_leadfield)


@pytest.fixture()
def fn_elec(sphere_vol, leadfield_elec):
    fn_hdf5 = 'tmp_elec.hdf5'
    if os.path.isfile(fn_hdf5):
        os.remove(fn_hdf5)
    sphere_vol.write_hdf5(fn_hdf5, 'mesh_leadfield')
    dset = '/mesh_leadfield/leadfields/tdcs_leadfield'
    with h5py.File(fn_hdf5, 'a') as f:
        f.create_dataset(dset, data=leadfield_elec)
        f[dset].attrs['electrode_tags'] = [1100, 1101]
        f[dset].attrs['electrode_names'] = ['A'.encode(), 'B'.encode()]
    yield fn_hdf5
    os.remove(fn_hdf5)

def simple_coil():
    dipole_pos = np.array([
        [10, 0, 0],
        [-10, 0, 0],
        [0, 10, 0],
        [0, -10, 0],
        [0, 0, 0],
        [0, 0, 1]],
        dtype=float
    )
    dipole_vec = np.repeat([[0., 0., 1.]], len(dipole_pos), axis=0)
    return dipole_pos, dipole_vec

@pytest.fixture()
def simple_coil_ccd():
    dipole_pos, dipole_vec = simple_coil()
    with tempfile.NamedTemporaryFile(suffix='.ccd', delete=False, mode='w') as f:
        fn_ccd = f.name
        f.write('# simple coil number of elements\n')
        f.write(f'{len(dipole_pos)}\n')
        f.write('# centers and weighted directions of the elements (magnetic dipoles)\n')
        for l in np.hstack([dipole_pos*1e-3, dipole_vec]):
            f.write(' '.join(str(x) for x in l))
            f.write('\n')
    yield fn_ccd
    os.remove(fn_ccd)

def rdm(a, b):
    return np.linalg.norm(
        a/np.linalg.norm(a) - b/np.linalg.norm(b)
    )

class TestTMSOpt:
    def test_get_coil_positions(self, sphere_msh):
        tms_opt = opt_struct.TMSoptimize()
        tms_opt.mesh = sphere_msh
        tms_opt.centre = np.array([95, 0, 0])
        tms_opt.distance = 5.
        tms_opt.pos_ydir = None
        tms_opt.search_radius = 10
        tms_opt.angle_resolution = 30
        positions = np.array(tms_opt._get_coil_positions())
        coil_centers = positions[:, :3, 3]
        y_dirs = positions[:, :3, 1]
        assert np.allclose(
            np.linalg.norm(coil_centers, axis=1),
            100, rtol=0.1)
        assert np.allclose(
            coil_centers[0::12], coil_centers[11::12]
        )
        for i in range(11):
            assert np.allclose(
                np.rad2deg(np.arccos(np.sum(y_dirs[i::12]*y_dirs[i+1::12], axis=1))),
                30)

    def test_get_target_region(self, sphere_msh):
        tms_opt = opt_struct.TMSoptimize()
        tms_opt.mesh = sphere_msh
        tms_opt.target = np.array([85, 0, 0])
        tms_opt.target_size = 10
        tms_opt.tissues = [3]
        target_region = tms_opt._get_target_region()
        bar = sphere_msh.elements_baricenters()
        assert np.all(
            np.linalg.norm(tms_opt.target - bar[target_region], axis=1) <= 10
        )
        assert np.all(sphere_msh.elm.tag1[target_region - 1] == 3)

    def test_direct(self, sphere_msh, simple_coil_ccd):
        tms_opt = opt_struct.TMSoptimize()
        tms_opt.fnamecoil = simple_coil_ccd
        tms_opt.mesh = sphere_msh
        tms_opt.didt = 1e6
        fn_hdf5 = tempfile.mktemp(".hdf5")
        tms_opt._name_hdf5 = MagicMock(return_value=fn_hdf5)

        coil_centers = [
            [150., 0., 0.],
            [-150., 0., 0.],
            [0., 150., 0.],
            [0., -150., 0.],
            [0., 0, 150.],
            [0., 0, -150.],
        ]
        pos_matrices = []
        for cc in coil_centers:
            z_dir = -np.array(cc)/np.linalg.norm(cc)
            y_dir = np.array([0., 1., 0.])
            if np.isclose(np.abs(z_dir.dot(y_dir)), 1):
                y_dir = np.array([1., 0., 0.])
            p = np.eye(4)
            p[:3, 0] = np.cross(y_dir, z_dir)
            p[:3, 1] = y_dir
            p[:3, 2] = z_dir
            p[:3, 3] = cc
            pos_matrices.append(p)

        target_pos, target_region = sphere_msh.find_closest_element(
            [85, 0, 0],
            elements_of_interest=sphere_msh.elm.tetrahedra,
            return_index=True
        )
        cond_field = sim_struct.SimuList.cond2elmdata(tms_opt)

        E_fem = tms_opt._direct_optimize(
            cond_field,
            np.atleast_1d(target_region),
            pos_matrices, 1
        )

        dipole_pos, dipole_moment = simple_coil()
        E_analytical = []
        for p in pos_matrices:
            dp = p[:3, :3].dot(dipole_pos.T).T + p[:3, 3]
            dm = p[:3, :3].dot(dipole_moment.T).T
            E = analytical_solutions.tms_E_field(
                dp * 1e-3, dm, tms_opt.didt,
                np.atleast_2d(target_pos) * 1e-3
            )
            E_analytical.append(np.linalg.norm(E))
        assert np.allclose(E_analytical, E_fem, rtol=0.1)


    @patch('simnibs.optimization.optimize_tms.get_opt_grid_ADM')
    def test_reciprocal(self, get_opt_grid_mock, sphere_msh, simple_coil_ccd):
        tms_opt = opt_struct.TMSoptimize()
        tms_opt.fnamecoil = simple_coil_ccd
        tms_opt.mesh = sphere_msh
        tms_opt.didt = 1e6

        coil_centers = [
            [150., 0., 0.],
            [-150., 0., 0.],
            [0., 150., 0.],
            [0., -150., 0.],
            [0., 0, 150.],
            [0., 0, -150.],
        ]
        center_matrices = []
        for cc in coil_centers:
            z_dir = -np.array(cc)/np.linalg.norm(cc)
            y_dir = np.array([0., 1., 0.])
            if np.isclose(np.abs(z_dir.dot(y_dir)), 1):
                y_dir = np.array([1., 0., 0.])
            p = np.eye(4)
            p[:3, 0] = np.cross(y_dir, z_dir)
            p[:3, 1] = y_dir
            p[:3, 2] = z_dir
            p[:3, 3] = cc
            center_matrices.append(p)
        coil_dir = []
        for angle in np.linspace(-np.pi/2, np.pi/2, 7):
            coil_dir.append([-np.sin(angle), np.cos(angle), 0])
        coil_dir = np.array(coil_dir)

        get_opt_grid_mock.return_value = (
                np.array(center_matrices).transpose(1, 2, 0),
                coil_dir.T
        )

        target_pos, target_region = sphere_msh.find_closest_element(
            [85, 0, 0],
            elements_of_interest=sphere_msh.elm.tetrahedra,
            return_index=True
        )
        cond_field = sim_struct.SimuList.cond2elmdata(tms_opt)

        E_recp, pos_matrices = tms_opt._ADM_optimize(cond_field, target_region)

        dipole_pos, dipole_moment = simple_coil()
        E_analytical = []
        for p in pos_matrices:
            dp = p[:3, :3].dot(dipole_pos.T).T + p[:3, 3]
            dm = p[:3, :3].dot(dipole_moment.T).T
            E = analytical_solutions.tms_E_field(
                dp * 1e-3, dm, tms_opt.didt,
                np.atleast_2d(target_pos) * 1e-3
            )
            E_analytical.append(np.linalg.norm(E))
        assert np.allclose(E_analytical, E_recp, rtol=0.1)


class TestFindIndexes:
    @pytest.mark.parametrize('indexes', [3, [5, 2]])
    def test_find_indexes_idx_node(self, indexes, sphere_surf):
        idx, mapping = opt_struct._find_indexes(sphere_surf, 'node',
                                                 indexes=indexes)

        assert np.all(np.atleast_1d(idx) == indexes)
        assert np.all(mapping == np.arange(len(idx)))

    @pytest.mark.parametrize('indexes', [3, [5, 2]])
    def test_find_indexes_idx_element(self, indexes, sphere_surf):
        idx, mapping = opt_struct._find_indexes(sphere_surf, 'element',
                                                 indexes=indexes)

        assert np.all(np.atleast_1d(idx) == indexes)
        assert np.all(mapping == np.arange(len(idx)))


    @pytest.mark.parametrize('r', [0, 10])
    @pytest.mark.parametrize('pos', [[85., 0., 0.], [[85., 0., 0.], [0., 85., 0.]]])
    @pytest.mark.parametrize('tissues', [None, [1004]])
    def test_find_indexes_pos_node(self, tissues, pos, sphere_surf, r):
        index, mapping = opt_struct._find_indexes(
            sphere_surf, 'node',
            positions=pos, radius=r,
            tissues=tissues)

        if tissues is None:
            nodes_with_tag = sphere_surf.nodes.node_number
        else:
            nodes_with_tag = sphere_surf.elm.nodes_with_tag(tissues)
        nodes = sphere_surf.nodes[nodes_with_tag]
        roi = []
        mp = []
        pos = np.atleast_2d(pos)
        for i, p in enumerate(pos):
            dist = np.linalg.norm(p - nodes, axis=1)
            idx_nearest = nodes_with_tag[np.argmin(dist)]
            if r > 0:
                center = nodes[np.argmin(dist)]
                dist_center = np.linalg.norm(center - nodes, axis=1)
                n_roi = np.sum(dist_center <= r)
                roi.extend(nodes_with_tag[dist_center <= r])
            else:
                n_roi = 1
                roi.append(idx_nearest)

            mp.extend(n_roi*(i,))

        assert np.all(np.sort(index) == np.sort(roi))
        assert np.all(np.sort(mapping) == np.sort(mp))

    @pytest.mark.parametrize('r', [0, 10])
    @pytest.mark.parametrize('pos', [[85., 0., 0.], [[85., 0., 0.], [0., 85., 0.]]])
    @pytest.mark.parametrize('tissues', [None, [1004]])
    def test_find_indexes_pos_elm(self, tissues, pos, sphere_surf, r):
        index, mapping = opt_struct._find_indexes(
            sphere_surf, 'element',
            positions=pos, radius=r,
            tissues=tissues)

        if tissues is None:
            elm_with_tag = sphere_surf.elm.elm_number
        else:
            elm_with_tag = sphere_surf.elm.elm_number[np.isin(sphere_surf.elm.tag1, tissues)]
        elms = sphere_surf.elements_baricenters()[elm_with_tag]
        roi = []
        mp = []
        pos = np.atleast_2d(pos)
        for i, p in enumerate(pos):
            dist = np.linalg.norm(p - elms, axis=1)
            idx_nearest = elm_with_tag[np.argmin(dist)]
            if r > 0:
                center = elms[np.argmin(dist)]
                dist_center = np.linalg.norm(center - elms, axis=1)
                n_roi = np.sum(dist_center <= r)
                roi.extend(elm_with_tag[dist_center <= r])
            else:
                n_roi = 1
                roi.append(idx_nearest)

            mp.extend(n_roi*(i,))

        assert np.all(np.sort(index) == np.sort(roi))
        assert np.all(np.sort(mapping) == np.sort(mp))


class TestFindDirections:

    @pytest.mark.parametrize('idx', [np.array([1]), np.array([1, 2])])
    @pytest.mark.parametrize('lf_type', ['node', 'element'])
    def test_find_directions_normal(self, idx, lf_type, sphere_surf):
        directions = opt_struct._find_directions(
            sphere_surf, lf_type, 'normal', idx
        )
        if lf_type == 'node':
            normals = sphere_surf.nodes_normals()[idx]
        elif lf_type == 'element':
            normals = sphere_surf.triangle_normals()[idx]
        assert np.allclose(directions, -normals)


    def test_find_directions_defined_1d(self):
        directions = opt_struct._find_directions(
            None, None, [1, 0, 0], [1]
        )
        assert directions.shape == (1, 3)
        assert np.allclose(directions, [[1, 0, 0]])

    def test_find_directions_defined_2d(self):
        directions = opt_struct._find_directions(
            None, None, [[1, 0, 0], [0, 1, 0]], [1, 2]
        )
        assert directions.shape == (2, 3)
        assert np.allclose(directions, [[1, 0, 0], [0, 1, 0]])

    def test_find_directions_defined_1_to_2(self):
        directions = opt_struct._find_directions(
            None, None, [1, 0, 0], [1, 2], [1, 2]
        )
        assert directions.shape == (2, 3)
        assert np.allclose(directions, [[1, 0, 0], [1, 0, 0]])


    def test_find_directions_defined_map(self):
        directions = opt_struct._find_directions(
            None, None, [[1, 0, 0], [0, 1, 0]], [1, 2, 3, 4], [0, 0, 1, 1]
        )
        assert directions.shape == (4, 3)
        assert np.allclose(directions, [[1, 0, 0], [1, 0, 0], [0, 1, 0], [0, 1, 0]])



class TestTDCSTarget:
    def test_create_mat_struct(self):
        targets = [opt_struct.TDCStarget(indexes=1, directions=[0, 1, 0], radius=5),
                   opt_struct.TDCStarget(indexes=[1, 2], intensity=.3, max_angle=30,
                                           tissues=[3, 4])]
        m = opt_struct._save_TDCStarget_mat(targets)
        assert np.all(m[0]['indexes'] == 1)
        assert np.all(m[0]['directions'] == [0, 1, 0])
        assert np.all(m[0]['radius'] == 5)
        assert np.all(m[1]['indexes'] == [1, 2])
        assert np.all(m[1]['intensity'] == .3)
        assert np.all(m[1]['max_angle'] == 30)
        assert np.all(m[1]['tissues'] == [3, 4])

    def test_read_mat_struct(self):
        m = {'indexes': [[1]], 'positions': [''],
             'directions': ['normal'],
             'intensity': [[0.5]], 'max_angle': [[30]],
             'radius': [[4.]], 'tissues': [[3, 2]]}
        t = opt_struct.TDCStarget.read_mat_struct(m)
        assert t.indexes == [1]
        assert t.positions is None
        assert t.directions == 'normal'
        assert t.intensity == 0.5
        assert t.max_angle == 30.
        assert t.radius == 4.
        assert t.tissues == [3, 2]

    def test_read_mat_directions_array(self):
        m = {'indexes': [''],
             'positions': [[1., 2., 3.]],
             'directions': [[0., 0., 1.]]}
        t = opt_struct.TDCStarget.read_mat_struct(m)
        assert t.indexes is None
        assert np.allclose(t.directions, [[0, 0, 1]])
        assert np.allclose(t.positions, [[1, 2, 3]])
        assert np.allclose(t.positions, [[1, 2, 3]])
        assert t.tissues is None

    def test_read_mat_directions_none(self):
        m = {'directions': ['none']}
        t = opt_struct.TDCStarget.read_mat_struct(m)
        assert t.directions is None


    def test_mat_io(self):
        targets = [opt_struct.TDCStarget(indexes=1, directions=[0, 1, 0]),
                   opt_struct.TDCStarget(indexes=[1, 2], intensity=.3, max_angle=30)]
        m = opt_struct._save_TDCStarget_mat(targets)
        scipy.io.savemat('tmp.mat', {'targets': m})
        m = scipy.io.loadmat('tmp.mat', struct_as_record=True, squeeze_me=False)
        os.remove('tmp.mat')
        t = opt_struct.TDCStarget.read_mat_struct(m['targets'][0][0])
        assert t.indexes == [1]
        assert np.allclose(t.directions, [0, 1, 0])
        t = opt_struct.TDCStarget.read_mat_struct(m['targets'][0][1])
        assert np.all(t.indexes == [1, 2])
        assert t.directions == 'normal'
        assert t.intensity == .3
        assert t.max_angle == 30


    def test_get_indexes_and_directions(self, sphere_surf):
        idx = [1]
        directions = [2., 0., 0.]
        t = opt_struct.TDCStarget(
            indexes=idx, directions=directions,
            mesh=sphere_surf, lf_type='node')
        id_, dir_ = t.get_indexes_and_directions()
        assert np.all(id_ == [0])
        assert np.allclose(dir_, [1., 0., 0.])


    def test_get_indexes_and_directions_none(self, sphere_surf):
        idx = [1]
        directions = None
        t = opt_struct.TDCStarget(
            indexes=idx, directions=directions,
            mesh=sphere_surf, lf_type='node')
        id_, dir_ = t.get_indexes_and_directions()
        assert np.all(id_ == [0])
        assert dir_ is None


    def test_get_indexes_and_directions_2_targets(self, sphere_surf):
        idx = [1, 2]
        directions = [[2., 0., 0.], [0., 3., 0.]]
        t = opt_struct.TDCStarget(
            indexes=idx, directions=directions,
            mesh=sphere_surf, lf_type='node')
        id_, dir_ = t.get_indexes_and_directions()
        assert np.all(id_ == [0, 1])
        assert np.allclose(dir_, [[1., 0., 0.], [0., 1., 0.]])

    def test_get_indexes_and_directions_2_targets_1_dir(self, sphere_surf):
        idx = [1, 2]
        directions = [[2., 0., 0.]]
        t = opt_struct.TDCStarget(
            indexes=idx, directions=directions,
            mesh=sphere_surf, lf_type='node')
        id_, dir_ = t.get_indexes_and_directions()
        assert np.all(id_ == [0, 1])
        assert np.allclose(dir_, [[1., 0., 0.], [1., 0., 0.]])

    def test_get_indexes_and_directions_radius(self, sphere_vol):
        bar = sphere_vol.elements_baricenters().value
        idx_ = np.where(
            (np.linalg.norm(bar - bar[0], axis=1) < 20) *
            (sphere_vol.elm.tag1 == 4))[0]
        directions = [[2., 0., 0.]]

        t = opt_struct.TDCStarget(
            positions=bar[0], directions=directions,
            mesh=sphere_vol, lf_type='element', radius=20, tissues=4)

        id_, dir_ = t.get_indexes_and_directions()

        directions_ = [(1., 0., 0.)] * len(idx_)

        assert np.all(id_ == idx_)
        assert np.allclose(dir_, directions_)

    @pytest.mark.parametrize('lf_type', ['node', 'element'])
    @pytest.mark.parametrize('intensity', [0.2, -0.2])
    def test_as_field(self, intensity, lf_type, sphere_surf):
        t = opt_struct.TDCStarget(
            indexes=[1, 2], directions=[[1, 0, 0], [0, 2, 0]],
            mesh=sphere_surf, lf_type=lf_type, intensity=intensity)
        d = t.as_field()
        assert np.allclose(d[1], [intensity, 0, 0])
        assert np.allclose(d[2], [0, intensity, 0])
        assert np.allclose(d[3:], 0)

    def test_as_field_radius(self, sphere_vol):
        bar = sphere_vol.elements_baricenters().value
        t = opt_struct.TDCStarget(
            positions=bar[0], directions=[[1, 0, 0]],
            mesh=sphere_vol, lf_type='element', intensity=1.,
            radius=20, tissues=4)
        d = t.as_field()
        in_r = (
            (np.linalg.norm(bar - bar[0], axis=1) < 20) *
            (sphere_vol.elm.tag1 == 4))
        assert np.allclose(d[in_r], [1., 0, 0])
        assert np.allclose(d[~in_r], [0, 0, 0])

    @pytest.mark.parametrize('lf_type', ['node', 'element'])
    def test_as_field_none(self, lf_type, sphere_surf):
        t = opt_struct.TDCStarget(
            indexes=[1, 2], directions=None,
            mesh=sphere_surf, lf_type=lf_type, intensity=2
        )
        d = t.as_field()
        assert np.allclose(d[[1, 2]], 2)
        assert np.allclose(d[3:], 0)

    def test_mean_intensity(self, sphere_vol):
        t = opt_struct.TDCStarget(
            indexes=[1, 2],
            directions=[[1, 0, 0], [0, 1, 0]],
            mesh=sphere_vol, lf_type='element',
            radius=0)
        f = mesh_io.ElementData(
            np.zeros((sphere_vol.elm.nr, 3)))
        f[1] = [2, 1, 3]
        f[2] = [1, 4, 2]
        vols = sphere_vol.elements_volumes_and_areas()[:]
        m = np.average([2, 4], weights=vols[:2])
        assert np.isclose(t.mean_intensity(f), m)

    def test_mean_intensity_none(self, sphere_vol):
        t = opt_struct.TDCStarget(
            indexes=[1, 2],
            directions=None,
            mesh=sphere_vol, lf_type='element',
            radius=0)
        f = mesh_io.ElementData(
            np.zeros((sphere_vol.elm.nr, 3))
        )
        f[1] = [np.sqrt(2), np.sqrt(2), 0]
        f[2] = [np.sqrt(3), np.sqrt(3), np.sqrt(3)]
        vols = sphere_vol.elements_volumes_and_areas()[:]
        m = np.average([2, 3], weights=vols[:2])
        assert np.isclose(t.mean_intensity(f), m)


    def test_mean_angle(self, sphere_vol):
        t = opt_struct.TDCStarget(
            indexes=[1, 2],
            directions=[[1, 0, 0], [0, 1, 0]],
            mesh=sphere_vol, lf_type='element',
            radius=0)
        f = mesh_io.ElementData(
            np.zeros((sphere_vol.elm.nr, 3)))
        f[1] = [1, 1, 0]
        f[2] = [2, 2, 0]
        assert np.isclose(t.mean_angle(f), 45)



class TestTDCSAvoid:
    def test_create_mat_struct(self):
        targets = [opt_struct.TDCSavoid(indexes=1, radius=5),
                   opt_struct.TDCSavoid(positions=[1., 0., 3.], weight=1e4,
                                           tissues=[3, 4])]
        m = opt_struct._save_TDCSavoid_mat(targets)
        assert np.all(m[0]['indexes'] == 1)
        assert np.all(m[0]['radius'] == 5)
        assert np.all(m[1]['positions'] == [1, 0., 3.])
        assert np.all(m[1]['weight'] == 1e4)
        assert np.all(m[1]['tissues'] == [3, 4])

    def test_read_mat_struct(self):
        m = {'indexes': [[1]], 'positions': [''],
             'weight': [[1e4]],
             'radius': [[4.]],
             'tissues': [[3, 2]]}
        t = opt_struct.TDCSavoid.read_mat_struct(m)
        assert t.indexes == [1]
        assert t.positions is None
        assert t.weight == 1e4
        assert t.radius == 4.
        assert t.tissues == [3, 2]
        m = {'positions': [[1., 2., 3.]]}
        t = opt_struct.TDCSavoid.read_mat_struct(m)
        assert t.indexes is None
        assert np.allclose(t.positions, [[1, 2, 3]])
        assert t.tissues is None

    def test_avoid_field_node(self, sphere_surf):
        a = opt_struct.TDCSavoid(indexes=2,
                                 weight=1e4,
                                 lf_type='node',
                                 mesh=sphere_surf)
        f = a.avoid_field()
        in_r = np.zeros(sphere_surf.nodes.nr, dtype=bool)
        in_r[1] = True
        assert np.allclose(f[in_r], 1e4)
        assert np.allclose(f[~in_r], 1)


    def test_avoid_field_elm(self, sphere_surf):
        a = opt_struct.TDCSavoid(indexes=2,
                                 weight=1e4,
                                 lf_type='element',
                                 mesh=sphere_surf)
        f = a.avoid_field()
        in_r = np.zeros(sphere_surf.elm.nr, dtype=bool)
        in_r[1] = True
        assert np.allclose(f[in_r], 1e4)
        assert np.allclose(f[~in_r], 1)

    def test_avoid_field_none_elm(self, sphere_surf):
        a = opt_struct.TDCSavoid(tissues=[1003],
                                 weight=1e4,
                                 lf_type='element',
                                 mesh=sphere_surf)
        f = a.avoid_field()
        assert np.allclose(f[sphere_surf.elm.tag1 == 1003], 1e4)
        assert np.allclose(f[sphere_surf.elm.tag1 != 1003], 1)

    def test_avoid_field_elm_radius(self, sphere_surf):
        bar = sphere_surf.elements_baricenters()[:]
        a = opt_struct.TDCSavoid(positions=bar[0],
                                 radius=10,
                                 weight=1e4,
                                 lf_type='element',
                                 mesh=sphere_surf)
        f = a.avoid_field()
        in_r = np.linalg.norm(bar - bar[0], axis=1) < 10
        assert np.allclose(f[in_r], 1e4)
        assert np.allclose(f[~in_r], 1)


    def test_avoid_field_none_node(self, sphere_surf):
        a = opt_struct.TDCSavoid(tissues=[1003],
                                 weight=1e4,
                                 lf_type='node',
                                 mesh=sphere_surf)
        f = a.avoid_field()
        roi = np.linalg.norm(sphere_surf.nodes[:], axis=1) < 86
        assert np.allclose(f[roi], 1e4)
        assert np.allclose(f[~roi], 1)

    def test_avoid_mean_field_norm_elm(self, sphere_surf):
        a = opt_struct.TDCSavoid(tissues=[1003],
                                 weight=1e4,
                                 lf_type='element',
                                 mesh=sphere_surf)
        field = mesh_io.ElementData(np.ones((sphere_surf.elm.nr, 3)))
        field[sphere_surf.elm.tag1 == 1003] = [2, 0, 0]
        assert np.isclose(a.mean_field_norm_in_region(field), 2)

    def test_avoid_mean_field_norm_node(self, sphere_surf):
        a = opt_struct.TDCSavoid(tissues=[1003],
                                 weight=1e4,
                                 lf_type='node',
                                 mesh=sphere_surf)
        field = mesh_io.NodeData(np.ones((sphere_surf.nodes.nr, 3)))
        assert np.isclose(a.mean_field_norm_in_region(field), np.sqrt(3))



class TestTDCSoptimize:
    def test_prepare_read_mesh_surf(self, fn_surf, sphere_surf):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_surf)
        assert np.all(np.isclose(
            p.mesh.nodes.node_coord, sphere_surf.nodes.node_coord))
        assert np.all(
            p.mesh.elm.node_number_list==sphere_surf.elm.node_number_list)

    def test_prepare_read_mesh_vol(self, fn_vol, sphere_vol):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_vol)
        assert np.all(np.isclose(
            p.mesh.nodes.node_coord, sphere_vol.nodes.node_coord))
        assert np.all(
            p.mesh.elm.node_number_list==sphere_vol.elm.node_number_list)

    def test_prepare_set_mesh_vol(self, fn_vol, sphere_surf):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_vol)
        p.mesh = sphere_surf
        assert np.all(np.isclose(
            p.mesh.nodes.node_coord, sphere_surf.nodes.node_coord))
        assert np.all(
            p.mesh.elm.node_number_list==sphere_surf.elm.node_number_list)

    def test_read_lf(self, fn_surf, leadfield_surf):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_surf)
        assert np.all(np.isclose(p.leadfield, leadfield_surf))
        p.leadfield

    def test_set_lf(self, fn_surf, leadfield_vol):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_surf)
        p.leadfield = leadfield_vol
        assert np.all(np.isclose(p.leadfield, leadfield_vol))

    def test_read_lftype(self, fn_surf):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_surf)
        assert p.lf_type == 'node'

    def test_read_lftype_elm(self, fn_vol):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_vol)
        assert p.lf_type == 'element'

    def test_read_lftype_wrong(self, fn_vol, sphere_surf):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_vol)
        p.mesh = sphere_surf
        with pytest.raises(ValueError):
            p.lf_type


    def test_get_avoid_field_vol(self, fn_vol, sphere_vol):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_vol)
        t = p.add_avoid()
        t.tissues = 4
        t.weight = 1e5
        avoid_field = p._get_avoid_field()
        assert np.allclose(avoid_field[sphere_vol.elm.tag1 == 4], 1e5)
        assert np.allclose(avoid_field[sphere_vol.elm.tag1 != 4], 1)

    def test_add_target(self):
        p = opt_struct.TDCSoptimize()
        t = p.add_target()
        assert t is p.target[-1]

    def test_add_target_arg(self):
        t = 'taget'
        p = opt_struct.TDCSoptimize()
        p.add_target(t)
        assert t is p.target[-1]

    def test_create_mat_struct(self):
        p = opt_struct.TDCSoptimize(
            leadfield_hdf='a.hdf5',
            max_total_current=2.,
            max_individual_current=.1,
            max_active_electrodes=3,
            name='a',
            open_in_gmsh=False)
        m = p.to_mat()
        assert m['leadfield_hdf'] == 'a.hdf5'
        assert m['max_total_current'] == 2.
        assert m['max_individual_current'] == .1
        assert m['max_active_electrodes'] == 3
        assert m['name'] == 'a'
        assert m['open_in_gmsh'] == False

    def test_read_mat_struct(self):
        m = {'leadfield_hdf': ['a.hdf5'],
             'max_total_current': [2.],
             'max_individual_current': [.1],
             'max_active_electrodes': [3],
             'name': ['aaa'],
             'target': [],
             'avoid': [],
             'open_in_gmsh': [False]}
        p = opt_struct.TDCSoptimize.read_mat_struct(m)
        assert p.leadfield_hdf == 'a.hdf5'
        assert p.max_total_current == 2.
        assert p.max_individual_current == .1
        assert p.max_active_electrodes == 3
        assert p.name == 'aaa'
        assert p.open_in_gmsh == False
        m = {'leadfield_hdf': ['a.hdf5'],
             'max_total_current': [2.],
             'max_individual_current': [.1],
             'max_active_electrodes': [''],
             'target': [],
             'avoid': []}
        p = opt_struct.TDCSoptimize.read_mat_struct(m)
        assert p.max_active_electrodes is None

    @pytest.mark.parametrize('max_active_electrodes', [3, None])
    def test_mat_io(self, max_active_electrodes):
        p = opt_struct.TDCSoptimize(
            leadfield_hdf='a.hdf5',
            max_total_current=2.,
            max_individual_current=.1,
            max_active_electrodes=max_active_electrodes,
            name='aaa')
        m = p.to_mat()
        scipy.io.savemat('tmp.mat', m)
        m = scipy.io.loadmat('tmp.mat', struct_as_record=True, squeeze_me=False)
        os.remove('tmp.mat')
        p = opt_struct.TDCSoptimize.read_mat_struct(m)
        assert p.leadfield_hdf == 'a.hdf5'
        assert p.max_total_current == 2.
        assert p.max_individual_current == .1
        assert p.max_active_electrodes == max_active_electrodes
        assert p.name == 'aaa'


    @pytest.mark.parametrize('intensity', [3e-4, -3e-4])
    @pytest.mark.parametrize('max_el_c', [1e-3, None])
    @pytest.mark.parametrize('max_tot_c', [2e-3, None])
    @pytest.mark.parametrize('max_ac', [None, 3])
    @pytest.mark.parametrize('max_angle', [None, 30])
    @pytest.mark.parametrize('n_targets', [1, 3])
    def test_optimize(self, intensity, max_el_c, max_tot_c, max_ac,
                      max_angle, n_targets, sphere_surf, fn_surf, leadfield_surf):

        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_surf,
                              max_individual_current=max_el_c,
                              max_total_current=max_tot_c,
                              max_active_electrodes=max_ac)

        for i in range(n_targets):
            t = p.add_target()
            t.indexes = i + 1
            t.directions = [1., 0., 0.]
            t.intensity = intensity

        if max_el_c is None and max_tot_c is None:
            pass
        elif n_targets > 1 and max_angle is not None:
            pass
        else:
            currents = p.optimize()
            assert np.isclose(np.sum(currents), 0, atol=1e-6)
            if max_el_c is not None:
                assert np.max(np.abs(currents)) < max_el_c * 1.05
            if max_tot_c is not None:
                assert np.linalg.norm(currents, 1) < 2 * max_tot_c * 1.05
            if max_ac is not None:
                assert np.linalg.norm(currents, 0) <= max_ac
            for i in range(n_targets):
                field = currents[1:].dot(leadfield_surf[:, i, :])
                assert np.sign(field[0]) == np.sign(intensity)

    @pytest.mark.parametrize('max_el_c', [1e-3, None])
    @pytest.mark.parametrize('max_tot_c', [2e-3, None])
    @pytest.mark.parametrize('max_ac', [None, 3])
    @pytest.mark.parametrize('n_targets', [1, 3])
    def test_optimize_norm(self, max_el_c, max_tot_c, max_ac, n_targets, sphere_surf, fn_surf, leadfield_surf):

        intensity = 3e-5
        p = opt_struct.TDCSoptimize(
            leadfield_hdf=fn_surf,
            max_individual_current=max_el_c,
            max_total_current=max_tot_c,
            max_active_electrodes=max_ac
        )

        for i in range(n_targets):
            t = p.add_target()
            t.indexes = i + 1
            t.directions = None
            t.intensity = intensity

        if max_el_c is None and max_tot_c is None:
            pass
        else:
            currents = p.optimize()
            assert np.isclose(np.sum(currents), 0, atol=1e-6)
            if max_el_c is not None:
                assert np.max(np.abs(currents)) < max_el_c * 1.05
            if max_tot_c is not None:
                assert np.linalg.norm(currents, 1) < 2 * max_tot_c * 1.05
            if max_ac is not None:
                assert np.linalg.norm(currents, 0) <= max_ac

    def test_field_node(self, leadfield_surf, fn_surf):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_surf)
        c = [1., -1., 0, 0., 0.]
        E = p.field(c)
        assert isinstance(E, mesh_io.NodeData)
        assert np.allclose(E.value, -leadfield_surf[0])

    def test_field_elm(self, leadfield_vol, fn_vol):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_vol)
        c = [1., -1., 0, 0., 0.]
        E = p.field(c)
        assert isinstance(E, mesh_io.ElementData)
        assert np.allclose(E.value, -leadfield_vol[0])

    @pytest.mark.parametrize('names', [None, ['A', 'B']])
    def test_currents_csv(self, names, fn_elec):
        csv_fn = 'test.csv'
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_elec)
        currents = [-0.2, 0.2]
        p.write_currents_csv(currents, csv_fn, electrode_names=names)
        with open(csv_fn) as f:
            reader = csv.reader(f)
            csv_rows = [row for row in reader]
        os.remove(csv_fn)
        assert csv_rows[0][0] == 'A'
        assert np.isclose(float(csv_rows[0][1]), currents[0])
        assert csv_rows[1][0] == 'B'
        assert np.isclose(float(csv_rows[1][1]), currents[1])



class TestTDCSDistributedoptimize:
    def test_prepare_read_mesh_surf(self, fn_surf, sphere_surf):
        p = opt_struct.TDCSDistributedOptimize(leadfield_hdf=fn_surf)
        assert np.all(np.isclose(
            p.mesh.nodes.node_coord, sphere_surf.nodes.node_coord))
        assert np.all(
            p.mesh.elm.node_number_list==sphere_surf.elm.node_number_list)

    def test_prepare_read_mesh_vol(self, fn_vol, sphere_vol):
        p = opt_struct.TDCSDistributedOptimize(leadfield_hdf=fn_vol)
        assert np.all(np.isclose(
            p.mesh.nodes.node_coord, sphere_vol.nodes.node_coord))
        assert np.all(
            p.mesh.elm.node_number_list==sphere_vol.elm.node_number_list)

    def test_prepare_set_mesh_vol(self, fn_vol, sphere_surf):
        p = opt_struct.TDCSDistributedOptimize(leadfield_hdf=fn_vol)
        p.mesh = sphere_surf
        assert np.all(np.isclose(
            p.mesh.nodes.node_coord, sphere_surf.nodes.node_coord))
        assert np.all(
            p.mesh.elm.node_number_list==sphere_surf.elm.node_number_list)

    def test_read_lf(self, fn_surf, leadfield_surf):
        p = opt_struct.TDCSDistributedOptimize(leadfield_hdf=fn_surf)
        assert np.all(np.isclose(p.leadfield, leadfield_surf))
        p.leadfield

    def test_set_lf(self, fn_surf, leadfield_vol):
        p = opt_struct.TDCSDistributedOptimize(leadfield_hdf=fn_surf)
        p.leadfield = leadfield_vol
        assert np.all(np.isclose(p.leadfield, leadfield_vol))

    def test_read_lftype(self, fn_surf):
        p = opt_struct.TDCSDistributedOptimize(leadfield_hdf=fn_surf)
        assert p.lf_type == 'node'

    def test_read_lftype_elm(self, fn_vol):
        p = opt_struct.TDCSDistributedOptimize(leadfield_hdf=fn_vol)
        assert p.lf_type == 'element'

    def test_read_lftype_wrong(self, fn_vol, sphere_surf):
        p = opt_struct.TDCSDistributedOptimize(leadfield_hdf=fn_vol)
        p.mesh = sphere_surf
        with pytest.raises(ValueError):
            p.lf_type

    def test_target_distribution_subject(self, sphere_surf, fn_surf):
        target_image = np.meshgrid(
            np.arange(-100, 100, 2),
            np.arange(-100, 100, 2),
            np.arange(-100, 100, 2),
            indexing='ij'
        )[0].astype(np.float)
        affine = np.eye(4)
        affine[:3, 3] = -100
        affine[:3, :3] *= 2
        min_img_value = 10
        p = opt_struct.TDCSDistributedOptimize(
            leadfield_hdf=fn_surf,
            target_image=(target_image, affine),
            mni_space=False,
            min_img_value=min_img_value,
            intensity=2
        )
        field, W = p._target_distribution()
        assert np.allclose(
            field[np.abs(field) > min_img_value],
            2*sphere_surf.nodes[np.abs(field) > min_img_value, 0],
            atol=1e-2
        )
        assert np.allclose(
            field[np.abs(sphere_surf.nodes[:, 0]) < min_img_value], 0
        )
        assert np.allclose(
            W[np.abs(sphere_surf.nodes[:, 0]) > min_img_value],
            np.abs(sphere_surf.nodes[np.abs(sphere_surf.nodes[:, 0]) > min_img_value, 0]),
            atol=1e-2
        )
        assert np.allclose(
            W[np.abs(sphere_surf.nodes[:, 0]) < min_img_value], min_img_value
        )

    def test_target_distribution_subject_file(self, sphere_surf, fn_surf):
        target_image = np.meshgrid(
            np.arange(-100, 100, 2),
            np.arange(-100, 100, 2),
            np.arange(-100, 100, 2),
            indexing='ij'
        )[0].astype(np.float)
        affine = np.eye(4)
        affine[:3, 3] = -100
        affine[:3, :3] *= 2
        fn_nii = 'tmp.nii.gz'
        nibabel.save(
            nibabel.Nifti1Image(target_image, affine),
            fn_nii
        )
        p = opt_struct.TDCSDistributedOptimize(
            leadfield_hdf=fn_surf,
            target_image=fn_nii,
            mni_space=False,
            min_img_value=0,
            intensity=1
        )
        field, W = p._target_distribution()
        os.remove(fn_nii)
        assert np.allclose(field, sphere_surf.nodes[:, 0], atol=1e-3)
        assert np.allclose(W, np.abs(sphere_surf.nodes[:, 0]), atol=1e-3)

    def test_target_distribution_subject_elm(self, sphere_vol, fn_vol):
        target_image = np.meshgrid(
            np.arange(-100, 100, 2),
            np.arange(-100, 100, 2),
            np.arange(-100, 100, 2),
            indexing='ij'
        )[0].astype(np.float)
        affine = np.eye(4)
        affine[:3, 3] = -100
        affine[:3, :3] *= 2
        p = opt_struct.TDCSDistributedOptimize(
            leadfield_hdf=fn_vol,
            target_image=(target_image, affine),
            mni_space=False,
            min_img_value=0,
            intensity=1
        )
        field, _ = p._target_distribution()
        assert np.allclose(field, sphere_vol.elements_baricenters()[:, 0], atol=1e-1)

    @patch('simnibs.utils.transformations.subject2mni_coords')
    def test_terget_field_mni(self, s2mni_coords_mock, sphere_surf, fn_surf):
        target_image = np.meshgrid(
            np.arange(-100, 100, 2),
            np.arange(-100, 100, 2),
            np.arange(-100, 100, 2),
            indexing='ij'
        )[0].astype(np.float)
        affine = np.eye(4)
        affine[:3, 3] = -100
        affine[:3, :3] *= 2
        s2mni_coords_mock.return_value = -sphere_surf.nodes[:]
        p = opt_struct.TDCSDistributedOptimize(
            leadfield_hdf=fn_surf,
            target_image=(target_image, affine),
            subpath='',
            min_img_value=0,
            intensity=1
        )
        field, _ = p._target_distribution()
        assert np.allclose(field, -sphere_surf.nodes[:, 0], atol=1e-2)

    def test_normals_surf(self, fn_surf, sphere_surf):
        p = opt_struct.TDCSDistributedOptimize(leadfield_hdf=fn_surf)
        n = p.normal_directions()
        assert np.allclose(n, -sphere_surf.nodes_normals()[:])

    def test_normals_vol(self, fn_vol, sphere_vol):
        p = opt_struct.TDCSDistributedOptimize(leadfield_hdf=fn_vol)
        with pytest.raises(ValueError):
            p.normal_directions()

    def test_field_node(self, leadfield_surf, fn_surf):
        p = opt_struct.TDCSDistributedOptimize(leadfield_hdf=fn_surf)
        c = [1., -1., 0, 0., 0.]
        E = p.field(c)
        assert isinstance(E, mesh_io.NodeData)
        assert np.allclose(E.value, -leadfield_surf[0])

    def test_field_msh_node(self, leadfield_surf, fn_surf, sphere_surf):
        target_distribution = np.meshgrid(
            np.arange(-100, 100, 2),
            np.arange(-100, 100, 2),
            np.arange(-100, 100, 2),
            indexing='ij'
        )[0].astype(np.float)
        affine = np.eye(4)
        affine[:3, 3] = -100
        affine[:3, :3] *= 2
        p = opt_struct.TDCSDistributedOptimize(
            leadfield_hdf=fn_surf,
            target_image=(target_distribution, affine),
            mni_space=False,
            min_img_value=0.,
            intensity=1,
        )
        c = [1., -1., 0, 0., 0.]
        m = p.field_mesh(c)
        normals = sphere_surf.nodes_normals()[:]
        assert isinstance(m.field['Field'], mesh_io.NodeData)
        assert np.allclose(m.field['Field'][:], -leadfield_surf[0])
        assert np.allclose(m.field['magnField'][:], np.linalg.norm(leadfield_surf[0], axis=1))
        assert np.allclose(m.field['normalField'][:], np.sum(leadfield_surf[0]*normals, axis=1))
        assert np.allclose(m.field['target_map'][:], sphere_surf.nodes[:, 0], atol=1e-3)

    @pytest.mark.parametrize('intensity', [3e-5, -2e-5])
    @pytest.mark.parametrize('max_el_c', [1e-3, None])
    @pytest.mark.parametrize('max_tot_c', [2e-3, None])
    @pytest.mark.parametrize('max_ac', [None, 3])
    def test_optimize(self, intensity, max_el_c, max_tot_c, max_ac,
                      sphere_surf, fn_surf, leadfield_surf):

        target_img = np.random.rand(100, 100, 100)
        affine = np.eye(4)
        affine[:3, 3] = -100
        affine[:3, :3] *= 2

        p = opt_struct.TDCSDistributedOptimize(
            leadfield_hdf=fn_surf,
            max_individual_current=max_el_c,
            max_total_current=max_tot_c,
            max_active_electrodes=max_ac,
            target_image=(target_img, affine),
            intensity=intensity,
            min_img_value=0,
            mni_space=False
        )

        if max_el_c is None and max_tot_c is None:
            pass
        else:
            currents = p.optimize()
            assert np.isclose(np.sum(currents), 0, atol=1e-6)
            if max_el_c is not None:
                assert np.max(np.abs(currents)) < max_el_c * 1.05
            if max_tot_c is not None:
                assert np.linalg.norm(currents, 1) < 2 * max_tot_c * 1.05
            if max_ac is not None:
                assert np.linalg.norm(currents, 0) <= max_ac


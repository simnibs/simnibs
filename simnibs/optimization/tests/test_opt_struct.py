import os
import csv
from mock import patch, call
import pytest
import numpy as np
import h5py
import scipy.io

from simnibs import SIMNIBSDIR
import simnibs.mesh_tools.mesh_io as mesh_io
from simnibs.optimization import opt_struct
import simnibs.optimization.optimization_methods as methods


@pytest.fixture()
def sphere_surf():
    fn = os.path.join(
        SIMNIBSDIR, 'resources', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn).crop_mesh([1003, 1004])


@pytest.fixture()
def sphere_vol():
    fn = os.path.join(
        SIMNIBSDIR, 'resources', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn).crop_mesh([4, 5])


@pytest.fixture
def sphere_elec():
    fn = os.path.join(
        SIMNIBSDIR, 'resources', 'testing_files',
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
def fn_elec(sphere_elec, sphere_vol, leadfield_elec):
    fn_hdf5 = 'tmp_elec.hdf5'
    if os.path.isfile(fn_hdf5):
        os.remove(fn_hdf5)
    sphere_elec.write_hdf5(fn_hdf5, 'mesh_electrodes')
    sphere_vol.write_hdf5(fn_hdf5, 'mesh_leadfield')
    dset = '/mesh_leadfield/leadfields/tdcs_leadfield'
    with h5py.File(fn_hdf5, 'a') as f:
        f.create_dataset(dset, data=leadfield_elec)
        f[dset].attrs['electrode_tags'] = [1100, 1101]
        f[dset].attrs['electrode_names'] = ['A'.encode(), 'B'.encode()]
    yield fn_hdf5
    os.remove(fn_hdf5)



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
        m = {'indexes': [''],
             'positions': [[1., 2., 3.]],
             'directions': [[0., 0., 1.]]}
        t = opt_struct.TDCStarget.read_mat_struct(m)
        assert t.indexes is None
        assert np.allclose(t.directions, [[0, 0, 1]])
        assert np.allclose(t.positions, [[1, 2, 3]])
        assert np.allclose(t.positions, [[1, 2, 3]])
        assert t.tissues is None

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


    def test_calc_target_matrices(self, sphere_surf, leadfield_surf):
        idx = [1]
        directions = [2., 0., 0.]
        t = opt_struct.TDCStarget(
            indexes=idx, directions=directions,
            mesh=sphere_surf, lf_type='node')
        t_mat, Q_mat = t.calc_target_matrices(leadfield_surf)
        areas = sphere_surf.nodes_volumes_or_areas().value
        idx_ = [0]
        directions_ = [1., 0., 0.]
        t, Q = methods.target_matrices(leadfield_surf, idx_, directions_, areas)
        assert np.all(np.isclose(t_mat, t))
        assert np.all(np.isclose(Q_mat, Q))

    def test_calc_target_matrices_2_targets(self, sphere_surf, leadfield_surf):
        idx = [1, 2]
        directions = [[2., 0., 0.], [0., 3., 0.]]
        t = opt_struct.TDCStarget(
            indexes=idx, directions=directions,
            mesh=sphere_surf, lf_type='node', radius=0)
        t_mat, Q_mat = t.calc_target_matrices(leadfield_surf)
        areas = sphere_surf.nodes_volumes_or_areas().value
        idx_ = [0, 1]
        directions_ = [[1., 0., 0.], [0., 1., 0.]]
        t, Q = methods.target_matrices(leadfield_surf, idx_, directions_, areas)
        assert np.all(np.isclose(t_mat, t))
        assert np.all(np.isclose(Q_mat, Q))

    def test_calc_target_matrices_2_targets_1_dir(self, sphere_surf, leadfield_surf):
        idx = [1, 2]
        directions = [[2., 0., 0.]]
        t = opt_struct.TDCStarget(
            indexes=idx, directions=directions,
            mesh=sphere_surf, lf_type='node', radius=0)
        t_mat, Q_mat = t.calc_target_matrices(leadfield_surf)
        areas = sphere_surf.nodes_volumes_or_areas().value
        idx_ = [0, 1]
        directions_ = [[1., 0., 0.], [1., 0., 0.]]
        t, Q = methods.target_matrices(leadfield_surf, idx_, directions_, areas)
        assert np.all(np.isclose(t_mat, t))
        assert np.all(np.isclose(Q_mat, Q))

    def test_calc_target_matrices_radius(self, sphere_vol, leadfield_vol):
        vols = sphere_vol.elements_volumes_and_areas().value
        bar = sphere_vol.elements_baricenters().value
        idx_ = np.where(
            (np.linalg.norm(bar - bar[0], axis=1) < 20) *
            (sphere_vol.elm.tag1 == 4))[0]
        directions = [[2., 0., 0.]]

        t = opt_struct.TDCStarget(
            positions=bar[0], directions=directions,
            mesh=sphere_vol, lf_type='element', radius=20, tissues=4)
        t_mat, Q_mat = t.calc_target_matrices(leadfield_vol)

        directions_ = [(1., 0., 0.)] * len(idx_)
        t, Q = methods.target_matrices(leadfield_vol, idx_, directions_, vols)
        assert np.all(np.isclose(t_mat, t))
        assert np.all(np.isclose(Q_mat, Q))

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

    def test_calc_energy_matrix_surf(self, fn_surf, sphere_surf, leadfield_surf):
        p = opt_struct.TDCSoptimize(leadfield_hdf=fn_surf)
        p.mesh = sphere_surf
        energy_mat = p.calc_energy_matrix()
        areas = sphere_surf.nodes_volumes_or_areas().value
        e_matrix = methods.energy_matrix(leadfield_surf, areas)
        assert np.all(np.isclose(energy_mat, e_matrix))

    def test_get_avoid_field(self, fn_vol, sphere_vol):
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


    @pytest.mark.parametrize('intensity', [3e-5, -2e-5])
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


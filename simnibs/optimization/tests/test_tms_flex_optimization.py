from copy import deepcopy
import os
from pathlib import Path

import numpy as np
import pytest
import scipy

from simnibs.mesh_tools.mesh_io import Msh
from simnibs.simulation.onlinefem import FemTargetPointCloud
from simnibs.simulation.tms_coil.tms_coil import TmsCoil

from simnibs.optimization.tms_flex_optimization import (
    optimize_distance,
    optimize_e_mag,
    TmsFlexOptimization,
)


class TestPositionOptimization:
    def test_initial_intersection(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh
    ):
        skin_surface = sphere3_msh.crop_mesh(tags=[1005])
        coil_affine = np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 50], [0, 0, 0, 1]]
        )
        try:
            optimize_distance(
                small_functional_3_element_coil, skin_surface, coil_affine, global_optimization=False, local_optimization=True
            )
        except Exception as e:
            raise pytest.fail("DID RAISE {0}".format(e))

    def test_no_deformations(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh
    ):
        coil = deepcopy(small_functional_3_element_coil)
        skin_surface = sphere3_msh.crop_mesh(tags=[1005])
        coil_affine = np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 100], [0, 0, 0, 1]]
        )
        for element in coil.elements:
            element.deformations = []

        with pytest.raises(ValueError):
            optimize_distance(coil, skin_surface, coil_affine)

    def test_no_casings(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh
    ):
        coil = deepcopy(small_functional_3_element_coil)
        skin_surface = sphere3_msh.crop_mesh(tags=[1005])
        coil_affine = np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 100], [0, 0, 0, 1]]
        )
        for element in coil.elements:
            element.casing = None

        with pytest.raises(ValueError):
            optimize_distance(coil, skin_surface, coil_affine)

    def test_simple_optimization(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh
    ):
        coil = deepcopy(small_functional_3_element_coil)
        skin_surface = sphere3_msh.crop_mesh(tags=[1005])
        cost_surface_tree = skin_surface.get_AABBTree()
        coil_affine = np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 100], [0, 0, 0, 1]]
        )

        distance_before = np.mean(
            np.sqrt(
                cost_surface_tree.min_sqdist(
                    coil.get_casing_coordinates(coil_affine)[0]
                )
            )
        )

        (before, after, affine_after, opt_ret) = optimize_distance(
            coil, skin_surface, coil_affine
        )

        distance_after = np.mean(
            np.sqrt(
                cost_surface_tree.min_sqdist(
                    coil.get_casing_coordinates(coil_affine)[0]
                )
            )
        )

        assert distance_before > distance_after
        assert (
            len(
                cost_surface_tree.points_inside(
                    coil.get_mesh(coil_affine).nodes.node_coord
                )
            )
            == 0
        )
        np.testing.assert_allclose(coil_affine, affine_after)

    def test_dithering_optimization(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh
    ):
        coil = deepcopy(small_functional_3_element_coil)
        skin_surface = sphere3_msh.crop_mesh(tags=[1005])
        cost_surface_tree = skin_surface.get_AABBTree()
        coil_affine = np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 100], [0, 0, 0, 1]]
        )

        distance_before = np.mean(
            np.sqrt(
                cost_surface_tree.min_sqdist(
                    coil.get_casing_coordinates(coil_affine)[0]
                )
            )
        )

        (before, after, affine_after, opt_ret) = optimize_distance(
            coil, skin_surface, coil_affine, dither_skip=0
        )

        distance_after_no_dithering = np.mean(
            np.sqrt(
                cost_surface_tree.min_sqdist(
                    coil.get_casing_coordinates(coil_affine)[0]
                )
            )
        )

        coil = deepcopy(small_functional_3_element_coil)
        (before, after, affine_after, opt_ret) = optimize_distance(
            coil, skin_surface, coil_affine, dither_skip=4
        )

        distance_after_dithering = np.mean(
            np.sqrt(
                cost_surface_tree.min_sqdist(
                    coil.get_casing_coordinates(coil_affine)[0]
                )
            )
        )

        assert distance_before > distance_after_no_dithering
        assert np.allclose(distance_after_dithering, distance_after_no_dithering)
        assert (
            len(
                cost_surface_tree.points_inside(
                    coil.get_mesh(coil_affine).nodes.node_coord
                )
            )
            == 0
        )
        np.testing.assert_allclose(coil_affine, affine_after)

    def test_optimization_with_global_translation(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh
    ):
        cost_surface_tree = sphere3_msh.crop_mesh(tags=[1005]).get_AABBTree()

        coil_affine = np.array(
            [[1, 0, 0, -4], [0, 1, 0, 3], [0, 0, 1, 110], [0, 0, 0, 1]]
        )

        (before, after, affine_after, opt_ret) = optimize_distance(
            small_functional_3_element_coil,
            sphere3_msh,
            coil_affine,
            np.array([[-5, 5], [-5, 5], [-20, 20]]),
        )

        assert before > after
        assert after < before * 0.5
        assert not np.allclose(coil_affine, affine_after)
        assert (
            len(
                cost_surface_tree.points_inside(
                    small_functional_3_element_coil.get_mesh(
                        coil_affine
                    ).nodes.node_coord
                )
            )
            == 0
        )

    def test_optimization_with_global_rotation(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh
    ):
        small_functional_3_element_coil_copy = deepcopy(small_functional_3_element_coil)
        skin_surface = sphere3_msh.crop_mesh(tags=[1005])
        cost_surface_tree = skin_surface.get_AABBTree()

        coil_affine = np.array(
            [[1, 0, 0, -4], [0, 1, 0, 3], [0, 0, 1, 110], [0, 0, 0, 1]]
        )

        del small_functional_3_element_coil_copy.elements[2].deformations[0]
        del small_functional_3_element_coil_copy.elements[1].deformations[0]

        (before, after, affine_after, opt_ret) = optimize_distance(
            small_functional_3_element_coil_copy,
            skin_surface,
            coil_affine,
            coil_rotation_ranges=np.array([[-90, 90], [-90, 90], [-90, 90]]),
        )

        assert before > after
        assert not np.allclose(coil_affine, affine_after)
        assert (
            len(
                cost_surface_tree.points_inside(
                    small_functional_3_element_coil.get_mesh(
                        affine_after
                    ).nodes.node_coord
                )
            )
            == 0
        )

    def test_self_intersection_optimization(
        self, small_self_intersecting_2_element_coil: TmsCoil, sphere3_msh: Msh
    ):
        coil_affine = np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 100], [0, 0, 0, 1]]
        )
        cost_surface_tree = (
            small_self_intersecting_2_element_coil.elements[0]
            .get_mesh(coil_affine)
            .get_AABBTree()
        )

        self_intersection_before = np.sum(
            cost_surface_tree.points_inside(
                small_self_intersecting_2_element_coil.elements[1]
                .get_mesh(coil_affine)
                .nodes.node_coord
            )
        )

        (before, after, affine_after, opt_ret) = optimize_distance(
            small_self_intersecting_2_element_coil,
            sphere3_msh,
            coil_affine,
            local_optimization=True,
        )

        self_intersection_after = np.sum(
            cost_surface_tree.points_inside(
                small_self_intersecting_2_element_coil.elements[1]
                .get_mesh(coil_affine)
                .nodes.node_coord
            )
        )

        assert self_intersection_before > self_intersection_after
        assert self_intersection_after <= 0.0000000001
        np.testing.assert_allclose(coil_affine, affine_after)


class TestEMagOptimization:
    @pytest.mark.slow
    def test_initial_intersection(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh, tmp_path: Path
    ):
        sphere_path = tmp_path.joinpath("sphere.msh")
        sphere3_msh.write(str(sphere_path.absolute()))
        sphere3_msh.fn = str(sphere_path.absolute())
        coil_affine = np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 50], [0, 0, 0, 1]]
        )
        roi = FemTargetPointCloud(
            sphere3_msh,
            center=sphere3_msh.elements_baricenters()[sphere3_msh.elm.tag1 == 4],
        )
        try:
            optimize_e_mag(
                small_functional_3_element_coil, sphere3_msh, roi, coil_affine
            )
        except Exception as e:
            raise pytest.fail("DID RAISE {0}".format(e))

    def test_no_deformations(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh
    ):
        coil = deepcopy(small_functional_3_element_coil)
        coil_affine = np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 100], [0, 0, 0, 1]]
        )
        for element in coil.elements:
            element.deformations = []

        with pytest.raises(ValueError):
            optimize_e_mag(coil, sphere3_msh, sphere3_msh.elm.tag1 == 4, coil_affine)

    def test_no_casings(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh
    ):
        coil = deepcopy(small_functional_3_element_coil)
        coil_affine = np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 100], [0, 0, 0, 1]]
        )
        for element in coil.elements:
            element.casing = None

        with pytest.raises(ValueError):
            optimize_e_mag(coil, sphere3_msh, sphere3_msh.elm.tag1 == 4, coil_affine)

    @pytest.mark.slow
    def test_simple_optimization(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh, tmp_path: Path
    ):
        sphere_path = tmp_path.joinpath("sphere.msh")
        sphere3_msh.write(str(sphere_path.absolute()))
        sphere3_msh.fn = str(sphere_path.absolute())
        coil_affine = np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 100], [0, 0, 0, 1]]
        )
        roi = FemTargetPointCloud(
            sphere3_msh,
            center=sphere3_msh.elements_baricenters()[sphere3_msh.elm.tag1 == 4],
        )
        (before, after, affine_after, e_mag, opt_ret) = optimize_e_mag(
            small_functional_3_element_coil, sphere3_msh, roi, coil_affine
        )

        assert before > after
        np.testing.assert_allclose(coil_affine, affine_after)

class TestMatlab:
    def test_writ_read_mat_no_deps(self, tmp_path: Path):
        tms_flex_opt = TmsFlexOptimization()

        tms_flex_opt.fnamecoil = "path/to/coil"

        tms_flex_opt.subpath = "path/to/m2m"
        tms_flex_opt.path_optimization = "path/to/output"
        tms_flex_opt.run_simulation = True

        tms_flex_opt.method = "emag"
        tms_flex_opt.global_translation_ranges = [[10, 20, 30]]
        tms_flex_opt.global_rotation_ranges = [[40, 50, 60]]

        tms_flex_opt.dither_skip = 7
        tms_flex_opt.fem_evaluation_cutoff = 2000

        tms_flex_opt.run_global_optimization = True
        tms_flex_opt.run_local_optimization = False

        tms_flex_opt.direct_args = {
            "maxiter": 3000,
            "len": True,
            "test": 1.0,
            "ar": "gs",
        }
        tms_flex_opt.l_bfgs_b_args = {"options": {"maxls": 20}}

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, tms_flex_opt.to_mat())

        tms_flex_opt_loaded = TmsFlexOptimization(scipy.io.loadmat(mat_path))

        assert tms_flex_opt.__dict__ == tms_flex_opt_loaded.__dict__

    def test_writ_read_mat(self, tmp_path: Path):
        tms_flex_opt = TmsFlexOptimization()

        tms_flex_opt.fnamecoil = "path/to/coil"

        tms_flex_opt.subpath = "path/to/m2m"
        tms_flex_opt.path_optimization = "path/to/output"
        tms_flex_opt.run_simulation = True

        pos = tms_flex_opt.add_position()
        pos.centre = [1.0, 2.0, 3.0]
        pos.pos_ydir = [1.2, 3.4, 5.6]
        pos.matsimnibs = []

        roi = tms_flex_opt.add_region_of_interest()
        roi.method = "surface"

        roi.surface_type = "central"
        roi.subpath = "path/to/m2m"

        roi.mask_space = ["subject_lh"]
        roi.mask_path = ["path_to_file"]
        roi.mask_value = [2]
        roi.mask_operator = ["intersection"]

        roi.roi_sphere_center = [[1.0,2,3], [4,5,6]]
        roi.roi_sphere_radius = [3, 45.0]
        roi.roi_sphere_center_space = ['subject', 'mni']
        roi.roi_sphere_operator = ["union", "intersection"]

        tms_flex_opt.method = "emag"
        tms_flex_opt.global_rotation_ranges = [[40, 50, 60]]

        tms_flex_opt.dither_skip = 7
        tms_flex_opt.fem_evaluation_cutoff = 2000

        tms_flex_opt.run_global_optimization = True
        tms_flex_opt.run_local_optimization = False

        tms_flex_opt.direct_args = {
            "maxiter": 3000,
            "len": True,
            "test": 1.0,
            "ar": "gs",
        }
        tms_flex_opt.l_bfgs_b_args = {"options": {"maxls": 20}}

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, tms_flex_opt.to_mat())

        tms_flex_opt_loaded = TmsFlexOptimization(scipy.io.loadmat(mat_path))

        dict_before = tms_flex_opt.__dict__.copy()
        dict_after = tms_flex_opt_loaded.__dict__.copy()
        del dict_before['pos']
        del dict_before['roi']
        del dict_after['pos']
        del dict_after['roi']

        assert dict_before == dict_after
        print(tms_flex_opt.roi.__dict__)
        print(tms_flex_opt_loaded.roi.__dict__)
        assert tms_flex_opt.roi.__dict__ == tms_flex_opt_loaded.roi.__dict__
        assert tms_flex_opt.pos.__dict__ == tms_flex_opt_loaded.pos.__dict__


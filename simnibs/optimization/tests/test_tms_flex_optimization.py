from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from simnibs.mesh_tools.mesh_io import Msh
from simnibs.simulation.onlinefem import FemTargetPointCloud
from simnibs.simulation.tms_coil.tms_coil import TmsCoil
from simnibs.simulation.tms_coil.tms_coil_element import (
    DipoleElements,
)
from simnibs.simulation.tms_coil.tms_stimulator import TmsStimulator

from simnibs.optimization.tms_flex_optimization import optimize_distance, optimize_e_mag


class TestPositionOptimization:
    def test_initial_intersection(
        self, small_functional_3_element_coil: TmsCoil, sphere3_msh: Msh
    ):
        skin_surface = sphere3_msh.crop_mesh(tags=[1005])
        coil_affine = np.array(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 50], [0, 0, 0, 1]]
        )
        try:
            optimize_distance(small_functional_3_element_coil, skin_surface, coil_affine)
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

        (
            before,
            after,
            affine_after,
            opt_ret
        ) = optimize_distance(coil, skin_surface, coil_affine)

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

        (
            before,
            after,
            affine_after,
            opt_ret
        ) = optimize_distance(coil, skin_surface, coil_affine, dither_skip=0)

        distance_after_no_dithering = np.mean(
            np.sqrt(
                cost_surface_tree.min_sqdist(
                    coil.get_casing_coordinates(coil_affine)[0]
                )
            )
        )

        coil = deepcopy(small_functional_3_element_coil)
        (
            before,
            after,
            affine_after,
            opt_ret
        ) = optimize_distance(coil, skin_surface, coil_affine, dither_skip=4)

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

        (
            before,
            after,
            affine_after,
            opt_ret
        ) = optimize_distance(small_functional_3_element_coil, 
            sphere3_msh, coil_affine, np.array([[-5, 5], [-5, 5], [-20, 20]])
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

        (
            before,
            after,
            affine_after,
            opt_ret
        ) = optimize_distance(small_functional_3_element_coil_copy, 
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

        (
            before,
            after,
            affine_after,
            opt_ret
        ) = optimize_distance(small_self_intersecting_2_element_coil, 
            sphere3_msh, coil_affine, local_optimization=True
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
            optimize_e_mag(small_functional_3_element_coil,
                sphere3_msh, roi, coil_affine
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
        (
            before,
            after,
            affine_after,
            e_mag,
            opt_ret
        ) = optimize_e_mag(small_functional_3_element_coil, 
            sphere3_msh, roi, coil_affine
        )

        assert before > after
        np.testing.assert_allclose(coil_affine, affine_after)


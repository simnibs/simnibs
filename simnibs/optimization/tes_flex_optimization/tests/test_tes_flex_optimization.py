import os
from pathlib import Path

import numpy as np
import pytest
import scipy

from simnibs.optimization.tes_flex_optimization.tes_flex_optimization import TesFlexOptimization
from simnibs.utils.matlab_read import dict_from_matlab
from simnibs.mesh_tools.mesh_io import read_msh
from simnibs.mesh_tools import surface
from simnibs.optimization.tes_flex_optimization.tes_flex_optimization import valid_skin_region
from simnibs.optimization.tes_flex_optimization.ellipsoid import Ellipsoid
from simnibs.utils.file_finder import Templates


class TestToFromDict:
    def test_write_read_mat(self, tmp_path: Path):
        opt = TesFlexOptimization()
        opt.subpath = "m2m_ernie"
        opt.output_folder = "tes_optimze_4x1tes_focality"

        ''' Set up goal function '''
        opt.goal = "focality"
        opt.threshold = [0.1, 0.2]
        opt.e_postproc = "magn"

        electrode = opt.add_electrode_layout("CircularArray")
        electrode.radius_inner = 10
        electrode.radius_outer = 10
        electrode.distance_bounds = [25, 100]
        electrode.n_outer = 4
        electrode.dirichlet_correction = False
        electrode.dirichlet_correction_detailed = False
        electrode.current = [0.002, -0.002/4, -0.002/4, -0.002/4, -0.002/4]

        roi = opt.add_roi()
        roi.method = "surface"
        roi.surface_type = "central"
        roi.roi_sphere_center_space = "subject"
        roi.roi_sphere_center = [-41.0, -13.0,  66.0]
        roi.roi_sphere_radius = 20

        roi = opt.add_roi()
        roi.method = "surface"
        roi.surface_type = "central"
        roi.roi_sphere_center_space = "subject"
        roi.roi_sphere_center = [-41.0, -13.0,  66.0]
        roi.roi_sphere_radius = 25
        roi.roi_sphere_operator = "difference"

        mat_path = os.path.join(tmp_path, "test.mat")

        scipy.io.savemat(mat_path, opt.to_dict())
        tes_flex_opt_loaded = TesFlexOptimization(dict_from_matlab(scipy.io.loadmat(mat_path)))

        np.testing.assert_equal(opt.to_dict(), tes_flex_opt_loaded.to_dict())

        dict_before = opt.__dict__.copy()
        dict_after = tes_flex_opt_loaded.__dict__.copy()
        del dict_before["_ff_templates"]
        del dict_before["roi"]
        del dict_before["_ellipsoid"]
        del dict_before["electrode"]

        del dict_after["_ff_templates"]
        del dict_after["roi"]
        del dict_after["_ellipsoid"]
        del dict_after["electrode"]

        np.testing.assert_equal(dict_before, dict_after)
        np.testing.assert_equal(opt.roi[0].__dict__, tes_flex_opt_loaded.roi[0].__dict__)
        np.testing.assert_equal(opt.roi[1].__dict__, tes_flex_opt_loaded.roi[1].__dict__)
        np.testing.assert_equal(opt.electrode[0].__dict__, tes_flex_opt_loaded.electrode[0].__dict__)

    @pytest.mark.slow
    def test_write_read_mat_after_prepare(self, tmp_path: Path, example_dataset):
        opt = TesFlexOptimization()
        # path of m2m folder containing the headmodel
        opt.subpath = os.path.join(example_dataset, 'm2m_ernie')
        opt.seed = 42                   # seed optimizer for reproducibility

        # output folder
        opt.output_folder = f"tes_optimize_ti_intensity"

        # type of goal function
        opt.goal = "mean"

        # postprocessing function of e-fields
        # "max_TI": maximize envelope of e-field magnitude
        # "dir_TI_normal": maximize envelope of e-field normal component
        # "dir_TI_tangential": maximize envelope of e-field tangential component
        opt.e_postproc = "max_TI"

        # define first pair of electrodes
        electrode = opt.add_electrode_layout("ElectrodeArrayPair")
        electrode.center = [[0, 0],                             # electrode center in reference electrode space (x-y plane)
                            [0, 20]]
        electrode.radius = [10, 10]                             # radius of electrodes
        electrode.dirichlet_correction_detailed = False         # node wise dirichlet correction
        electrode.current = [0.002, 0.002, -0.002, -0.002]      # electrode currents

        # define second pair of electrodes
        electrode = opt.add_electrode_layout("ElectrodeArrayPair")
        electrode.center = [[0, 0],                             # electrode center in reference electrode space (x-y plane)
                            [0, 20]]
        electrode.radius = [10, 10]                             # radius of electrodes
        electrode.dirichlet_correction_detailed = False         # node wise dirichlet correction
        electrode.current = [0.002, 0.002, -0.002, -0.002]      # electrode currents

        # define ROI
        roi = opt.add_roi()
        roi.method = "surface"
        roi.surface_type = "central"

        # center of spherical ROI in subject space (in mm)
        roi.roi_sphere_center_space = "subject"
        roi.roi_sphere_center = [-41.0, -13.0,  66.0]

        # radius of spherical ROI (in mm)
        roi.roi_sphere_radius = 20

        # prepare optimization
        opt._prepare()

        # save opt instance as .mat structure
        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, opt.to_dict())

        # load .mat structure and initialize new opt instance
        tes_flex_opt_loaded = TesFlexOptimization(dict_from_matlab(scipy.io.loadmat(mat_path)))

        # prepare loaded opt instance
        tes_flex_opt_loaded._prepare()

        # check if both opt instances are equal
        np.testing.assert_equal(opt.to_dict(), tes_flex_opt_loaded.to_dict())

    @pytest.mark.slow
    def test_run_opt_from_mat(self, tmp_path: Path, example_dataset):
        opt = TesFlexOptimization()
        opt.seed = 42                   # seed optimizer for reproducibility
        opt.open_in_gmsh = False
        opt.detailed_results = True

        # path of m2m folder containing the headmodel
        opt.subpath = os.path.join(example_dataset, 'm2m_ernie')

        # output folder
        opt.output_folder = f"tes_optimize_tes_intensity_org"

        # type of goal function
        opt.goal = "mean"

        # postprocessing function of e-fields
        opt.e_postproc = "magn"

        # define first pair of electrodes
        electrode = opt.add_electrode_layout("ElectrodeArrayPair")
        electrode.center = [[0, 0]]  # electrode center in reference electrode space (x-y plane)
        electrode.radius = [10]  # radius of electrodes
        electrode.dirichlet_correction_detailed = False  # node wise dirichlet correction
        electrode.current = [0.002, -0.002]  # electrode currents

        # define ROI
        roi = opt.add_roi()
        roi.method = "surface"
        roi.surface_type = "central"
        roi.roi_sphere_center_space = "subject"
        roi.roi_sphere_center = [-41.0, -13.0, 66.0]
        roi.roi_sphere_radius = 20

        # prepare optimization
        opt._prepare()

        # save opt instance as .mat structure
        mat_path = os.path.join(tmp_path, "opt.mat")
        scipy.io.savemat(mat_path, opt.to_dict())

        # load .mat structure and initialize new opt instance
        tes_flex_opt_loaded = TesFlexOptimization(dict_from_matlab(scipy.io.loadmat(mat_path)))
        tes_flex_opt_loaded.output_folder = f"tes_optimize_tes_intensity_mat"

        # prepare loaded opt instance
        tes_flex_opt_loaded._prepare()

        # run optimization (original)
        opt.run()

        # run optimization (.mat)
        tes_flex_opt_loaded.run()

        # compare results
        msh_org = read_msh(os.path.join(opt.output_folder, "ernie_tes_flex_opt_surface_mesh.msh"))
        msh_mat = read_msh(os.path.join(tes_flex_opt_loaded.output_folder, "ernie_tes_flex_opt_surface_mesh.msh"))

        np.testing.assert_allclose(msh_org.field["magnE"].value, msh_mat.field["magnE"].value)

    @pytest.mark.slow
    def test_compute_goal(self, tmp_path: Path, example_dataset):
        opt = TesFlexOptimization()

        # path of m2m folder containing the headmodel
        opt.subpath = os.path.join(example_dataset, 'm2m_ernie')

        # output folder
        opt.output_folder = f"tes_optimize_compute_goal"

        # postprocessing function of e-fields
        opt.e_postproc = "magn"

        # type of goal function
        opt.goal = "mean"

        # define first pair of electrodes
        electrode = opt.add_electrode_layout("ElectrodeArrayPair")
        electrode.center = [[0, 0]]  # electrode center in reference electrode space (x-y plane)
        electrode.radius = [10]  # radius of electrodes
        electrode.dirichlet_correction_detailed = False  # node wise dirichlet correction
        electrode.current = [0.002, -0.002]  # electrode currents

        # define ROI
        roi = opt.add_roi()
        roi.method = "surface"
        roi.surface_type = "central"
        roi.roi_sphere_center_space = "subject"
        roi.roi_sphere_center = [-41.0, -13.0, 66.0]
        roi.roi_sphere_radius = 20

        # prepare optimization
        opt._prepare()

        # mean
        opt.goal = ["mean"]
        e_test = [[np.array([1, 2, 3])]]
        goal_mean = opt.compute_goal(e_test)
        np.testing.assert_equal(goal_mean, -2)

        # neg_mean
        opt.goal = ["neg_mean"]
        e_test = [[np.array([1, 2, 3])]]
        goal_neg_mean = opt.compute_goal(e_test)
        np.testing.assert_equal(goal_neg_mean, 2)

        # mean_abs
        opt.goal = ["mean_abs"]
        e_test = [[np.array([-1, 2, 3])]]
        goal_mean_abs = opt.compute_goal(e_test)
        np.testing.assert_equal(goal_mean_abs, -2)

        # max
        opt.goal = ["max"]
        e_test = [[np.array([1, 2, 3])]]
        goal_max = opt.compute_goal(e_test)
        np.testing.assert_allclose(goal_max, -3, rtol=1e-2)

        # neg_max
        opt.goal = ["neg_max"]
        e_test = [[np.array([1, 2, 3])]]
        goal_neg_max = opt.compute_goal(e_test)
        np.testing.assert_allclose(goal_neg_max, 3, rtol=1e-2)

        # max_abs
        opt.goal = ["max_abs"]
        e_test = [[np.array([1, 2, -3])]]
        goal_max_abs = opt.compute_goal(e_test)
        np.testing.assert_allclose(goal_max_abs, -3, rtol=1e-2)

        # focality (one threshold)
        opt.goal = ["focality"]
        opt.threshold = [2]
        e_test = [[np.array([3, 4, 5]), np.array([0, 1, 2, 3])]]
        goal_focality_threshold_1 = opt.compute_goal(e_test)
        np.testing.assert_equal(np.round(goal_focality_threshold_1), -91)

        # focality (two thresholds)
        opt.goal = ["focality"]
        opt.threshold = [2, 5]
        e_test = [[np.array([3, 4, 5]), np.array([0, 1, 2, 3])]]
        goal_focality_threshold_2 = opt.compute_goal(e_test)
        np.testing.assert_equal(np.round(goal_focality_threshold_2), -58)

        # focality_inv (one threshold)
        opt.goal = ["focality_inv"]
        opt.threshold = [2]
        e_test = [[np.array([3, 4, 5]), np.array([0, 1, 2, 3])]]
        goal_focality_inv_threshold_1 = opt.compute_goal(e_test)
        np.testing.assert_equal(np.round(goal_focality_inv_threshold_1), -112)

        # focality_inv (two thresholds)
        opt.goal = ["focality_inv"]
        opt.threshold = [2, 5]
        e_test = [[np.array([3, 4, 5]), np.array([0, 1, 2, 3])]]
        goal_focality_inv_threshold_2 = opt.compute_goal(e_test)
        np.testing.assert_equal(np.round(goal_focality_inv_threshold_2), -60)

    @pytest.mark.slow
    def test_valid_skin_region(self, tmp_path: Path, example_dataset):
        ellipsoid = Ellipsoid()
        fn_electrode_mask = Templates().mni_volume_upper_head_mask

        mesh = read_msh(os.path.join(example_dataset, 'm2m_ernie', 'ernie.msh'))

        # relabel internal air
        mesh_relabel = mesh.relabel_internal_air()

        # make final skin surface including some additional distance
        skin_surface = surface.Surface(mesh=mesh_relabel, labels=1005)
        skin_surface = valid_skin_region(
            skin_surface=skin_surface,
            fn_electrode_mask=fn_electrode_mask,
            mesh=mesh_relabel,
            additional_distance=0,
        )

        np.testing.assert_equal(skin_surface.nodes.shape[0], 5427)
        np.testing.assert_allclose(skin_surface.nodes[0, :], [9.80649432, 104.20439836,  58.96501272], rtol=1e-6)

        # fit optimal ellipsoid to valid skin points
        ellipsoid.fit(points=skin_surface.nodes)

        np.testing.assert_allclose(ellipsoid.radii, [113.60465987, 105.76692858,  86.19796459], rtol=1e-6)
        np.testing.assert_allclose(ellipsoid.center, [ 1.47022674, 16.43981549, -4.79728615], rtol=1e-6)
        np.testing.assert_allclose(ellipsoid.rotmat, np.array([[-0.02356704, -0.06572103, -0.99755969],
                                                               [-0.77416654, -0.63015041,  0.05980489],
                                                               [ 0.63254309, -0.77368676,  0.03602824]]), rtol=1e-6)

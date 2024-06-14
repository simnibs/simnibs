import os
from pathlib import Path

import numpy as np
import scipy

from simnibs.optimization.tes_flex_optimization.tes_flex_optimization import TesFlexOptimization
from simnibs.utils.matlab_read import dict_from_matlab


class TestToFromDict:
    def test_writ_read_mat(self, tmp_path: Path):
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
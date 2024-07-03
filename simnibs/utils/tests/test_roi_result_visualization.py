from copy import deepcopy
import os

import numpy as np
from simnibs.mesh_tools import gmsh_view, mesh_io
from simnibs.simulation import sim_struct
from simnibs.utils.region_of_interest import RegionOfInterest
from simnibs.utils.roi_result_visualization import RoiResultVisualization


class TestReadGeo:
    def test_read_geo(self, tmp_path):
        geo_file_name = os.path.join(tmp_path, 'test.geo')
        mesh_io.write_geo_lines([[1, 3, 4], [8,9,10]], [[5, 6, 7], [11,12,13]], geo_file_name, name='lines')
        mesh_io.write_geo_spheres([[1,2,3], [3, 4, 5]], geo_file_name, name='spheres', mode='ba')
        mesh_io.write_geo_triangles([[0,1,2]], [[1,2,3], [4,5,6], [7,8,9]], geo_file_name, name='triangles', mode='ba')

        geo_dict = RoiResultVisualization._read_geo(geo_file_name)

        assert len(geo_dict) == 3
        assert 'lines' in geo_dict
        assert 'spheres' in geo_dict
        assert 'triangles' in geo_dict

        with open(geo_file_name,'r') as f:
            geo_file_contend = f.read()

        assert geo_file_contend == ''.join(geo_dict.values())

class TestReadOpt:
        def test_read_opt_base(self, tmp_path):
            opt_path = os.path.join(tmp_path, "surf.msh")

            opt = gmsh_view.Visualization(None, None)
            opt.add_view()

            opt.write_opt(opt_path)

            opt_read = RoiResultVisualization._read_opt(f'{opt_path}.opt')

            assert str(opt_read) == str(opt)

        def test_read_opt(self, tmp_path):
            opt_path = os.path.join(tmp_path, "surf.msh")

            opt = gmsh_view.Visualization(None, None)
            opt.add_view(
                CenterGlyphs = 2,
                GlyphLocation = 3,
                VectorType = 0,
                Visible = 1,
                CustomMax = 10,
                CustomMin = 20,
                SaturateValues = 1,
                RangeType = 2,
                ShowScale = 0,
                ColormapNumber = 3,
            )

            opt.add_view(
                CenterGlyphs = 0,
                GlyphLocation = 1,
                VectorType = 0,
                Visible = 1,
                CustomMax = 0,
                CustomMin = 10,
                SaturateValues = 1,
                RangeType = 3,
                ShowScale = 0,
                ColormapNumber = 4,
            )

            opt.add_view(
                CenterGlyphs = 2,
                GlyphLocation = 3,
                VectorType = 0,
                Visible = 0,
                CustomMax = -10,
                CustomMin = 10,
                SaturateValues = 1,
                RangeType = 2,
                ShowScale = 0,
                ColormapNumber = 5,
                ColorTable=gmsh_view._gray_red_lightblue_blue_cm()
            )

            opt.write_opt(opt_path)
            opt_read = RoiResultVisualization._read_opt(f'{opt_path}.opt')
            assert str(opt_read) == str(opt)

from copy import deepcopy
import os

import numpy as np
from simnibs.mesh_tools import gmsh_view, mesh_io
from simnibs.simulation import sim_struct
from simnibs.utils.region_of_interest import RegionOfInterest
from simnibs.utils.roi_result_visualization import RoiResultVisualization


class TestReadGeo:
    def test_read_geo(self, tmp_path):
        geo_file_name = os.path.join(tmp_path, "test.geo")
        mesh_io.write_geo_lines(
            [[1, 3, 4], [8, 9, 10]],
            [[5, 6, 7], [11, 12, 13]],
            geo_file_name,
            name="lines",
        )
        mesh_io.write_geo_spheres(
            [[1, 2, 3], [3, 4, 5]], geo_file_name, name="spheres", mode="ba"
        )
        mesh_io.write_geo_triangles(
            [[0, 1, 2]],
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
            geo_file_name,
            name="triangles",
            mode="ba",
        )

        geo_dict = RoiResultVisualization._read_geo(geo_file_name)

        assert len(geo_dict) == 3
        assert "lines" in geo_dict
        assert "spheres" in geo_dict
        assert "triangles" in geo_dict

        with open(geo_file_name, "r") as f:
            geo_file_contend = f.read()

        assert geo_file_contend == "".join(geo_dict.values())


class TestReadOpt:
    def test_read_opt_base(self, tmp_path):
        opt_path = os.path.join(tmp_path, "surf.msh")

        opt = gmsh_view.Visualization(None, None)
        opt.add_view()

        opt.write_opt(opt_path)

        opt_read = RoiResultVisualization._read_opt(f"{opt_path}.opt")

        assert str(opt_read) == str(opt)

    def test_read_opt(self, tmp_path):
        opt_path = os.path.join(tmp_path, "surf.msh")

        opt = gmsh_view.Visualization(None, None)
        opt.add_view(
            CenterGlyphs=2,
            GlyphLocation=3,
            VectorType=0,
            Visible=1,
            CustomMax=10,
            CustomMin=20,
            SaturateValues=1,
            RangeType=2,
            ShowScale=0,
            ColormapNumber=3,
        )

        opt.add_view(
            CenterGlyphs=0,
            GlyphLocation=1,
            VectorType=0,
            Visible=1,
            CustomMax=0,
            CustomMin=10,
            SaturateValues=1,
            RangeType=3,
            ShowScale=0,
            ColormapNumber=4,
        )

        opt.add_view(
            CenterGlyphs=2,
            GlyphLocation=3,
            VectorType=0,
            Visible=0,
            CustomMax=-10,
            CustomMin=10,
            SaturateValues=1,
            RangeType=2,
            ShowScale=0,
            ColormapNumber=5,
            ColorTable=gmsh_view._gray_red_lightblue_blue_cm(),
        )

        opt.write_opt(opt_path)
        opt_read = RoiResultVisualization._read_opt(f"{opt_path}.opt")
        assert str(opt_read) == str(opt)


class TestRoiResultVisualization:
    def test_volume_rois(self, sphere3_msh: mesh_io.Msh, tmp_path):
        input_mesh = deepcopy(sphere3_msh)
        #Create Roi 1
        roi_1 = RegionOfInterest()
        roi_1.load_mesh(input_mesh)
        roi_1.apply_tissue_mask(1003, "intersection")

        #Create Roi 2
        roi_2 = RegionOfInterest()
        roi_2.load_mesh(input_mesh)
        roi_2.apply_tissue_mask(1004, "intersection")

        #Create result mesh 1
        result_mesh_1_file_path = os.path.join(tmp_path, "test1_TMS_scalar.msh")
        result_mesh_1 = deepcopy(sphere3_msh)
        opt_1 = gmsh_view.Visualization(result_mesh_1)
        result_mesh_1.add_node_field(-np.arange(result_mesh_1.nodes.nr), "node_test_1")
        opt_1.add_view(ColormapNumber=1)
        result_mesh_1.add_node_field(
            -np.arange(result_mesh_1.nodes.nr) * 2, "node_test_2"
        )
        opt_1.add_view(ColormapNumber=2)
        result_mesh_1.add_element_field(np.arange(result_mesh_1.elm.nr), "elm_test_1")
        opt_1.add_view(ColormapNumber=3)
        result_mesh_1.add_element_field(
            np.arange(result_mesh_1.elm.nr) * 2, "elm_test_2"
        )
        opt_1.add_view(ColormapNumber=4)
        geo_file_path = os.path.join(tmp_path, "test1_TMS_coil_pos.geo")
        opt_1.add_merge(geo_file_path)
        mesh_io.write_geo_lines(
            [[1, 3, 4], [8, 9, 10]],
            [[5, 6, 7], [11, 12, 13]],
            geo_file_path,
            name="lines",
        )
        opt_1.add_view(ColormapNumber=5)
        mesh_io.write_geo_triangles(
            [[0, 1, 2]],
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
            geo_file_path,
            name="scalp",
            mode="ba",
        )
        opt_1.add_view(ColormapNumber=6)
        result_mesh_1.write(result_mesh_1_file_path)
        opt_1.write_opt(result_mesh_1_file_path)

        #Create result mesh 2
        result_mesh_2_file_path = os.path.join(tmp_path, "test2_TMS_scalar.msh")
        result_mesh_2 = deepcopy(sphere3_msh)
        opt_2 = gmsh_view.Visualization(result_mesh_2)
        result_mesh_2.add_node_field(
            -np.arange(result_mesh_2.nodes.nr) * 3, "node_test_1"
        )
        opt_2.add_view(ColormapNumber=7)
        result_mesh_2.add_node_field(
            -np.arange(result_mesh_2.nodes.nr) * 4, "node_test_2"
        )
        opt_2.add_view(ColormapNumber=8)
        result_mesh_2.add_element_field(
            np.arange(result_mesh_2.elm.nr) * 3, "elm_test_1"
        )
        opt_2.add_view(ColormapNumber=9)
        result_mesh_2.add_element_field(
            np.arange(result_mesh_2.elm.nr) * 4, "elm_test_2"
        )
        opt_2.add_view(ColormapNumber=10)
        geo_file_path = os.path.join(tmp_path, "test2_TMS_coil_pos.geo")
        opt_2.add_merge(geo_file_path)
        mesh_io.write_geo_spheres(
            [[1, 2, 3], [3, 4, 5]], geo_file_path, name="spheres", mode="ba"
        )
        opt_2.add_view(ColormapNumber=11)
        mesh_io.write_geo_triangles(
            [[0, 1, 2]],
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
            geo_file_path,
            name="scalp",
            mode="ba",
        )
        opt_2.add_view(ColormapNumber=12)
        result_mesh_2.write(result_mesh_2_file_path)
        opt_2.write_opt(result_mesh_2_file_path)

        result_vis = RoiResultVisualization(
            [roi_1, roi_2],
            [result_mesh_1_file_path, result_mesh_2_file_path],
            tmp_path,
            "test_vis",
            ["roi", "anti_roi"],
            ["mesh_1", "mesh_2"],
        )

        result_vis.get_head_mesh().add_node_field(
            -np.arange(result_mesh_2.nodes.nr) * 10, "node_test_y"
        )

        result_vis.get_head_mesh().add_element_field(
            np.arange(result_mesh_2.elm.nr) * 10, "elm_test_x"
        )

        result_vis.create_visualization()

        result_vis.head_mesh_data_name_to_gmsh_view['mesh_1__elm_test_2'].VectorType = 6
        result_vis.remove_field_from_head_mesh('mesh_2__elm_test_2')

        assert result_vis.head_mesh is not None
        assert result_vis.head_mesh_opt is not None
        assert isinstance(result_vis.head_mesh_opt.View, list)
        assert result_vis.surface_mesh is None

        #Test for all node data fields 
        assert result_vis.head_mesh.nodedata[0].field_name == 'mesh_1__node_test_1'
        assert int(result_vis.head_mesh_opt.View[0].ColormapNumber) == 1
        assert result_vis.head_mesh.nodedata[1].field_name == 'mesh_1__node_test_2'
        assert int(result_vis.head_mesh_opt.View[1].ColormapNumber) == 2
        assert result_vis.head_mesh.nodedata[2].field_name == 'mesh_2__node_test_1'
        assert int(result_vis.head_mesh_opt.View[2].ColormapNumber) == 7
        assert result_vis.head_mesh.nodedata[3].field_name == 'mesh_2__node_test_2'
        assert int(result_vis.head_mesh_opt.View[3].ColormapNumber) == 8
        assert result_vis.head_mesh.nodedata[4].field_name == 'node_test_y'
        assert int(result_vis.head_mesh_opt.View[4].ColormapNumber) == 2 # default value
        assert len(result_vis.head_mesh.nodedata) == 5

        #Test for the Roi element data fields
        assert result_vis.head_mesh.elmdata[0].field_name == 'roi'
        assert int(result_vis.head_mesh_opt.View[5].ColormapNumber) == 2 # default value
        assert int(result_vis.head_mesh_opt.View[5].ShowScale) == 0
        assert result_vis.head_mesh.elmdata[1].field_name == 'anti_roi'
        assert int(result_vis.head_mesh_opt.View[6].ColormapNumber) == 2 # default value
        assert int(result_vis.head_mesh_opt.View[6].ShowScale) == 0

        #Test for all element data fields
        assert result_vis.head_mesh.elmdata[2].field_name == 'mesh_1__elm_test_1'
        assert int(result_vis.head_mesh_opt.View[7].ColormapNumber) == 3
        assert result_vis.head_mesh.elmdata[3].field_name == 'mesh_1__elm_test_2'
        assert int(result_vis.head_mesh_opt.View[8].VectorType) == 6
        assert int(result_vis.head_mesh_opt.View[8].ColormapNumber) == 4
        assert result_vis.head_mesh.elmdata[4].field_name == 'mesh_2__elm_test_1'
        assert int(result_vis.head_mesh_opt.View[9].ColormapNumber) == 9
        assert result_vis.head_mesh.elmdata[5].field_name == 'elm_test_x'
        assert int(result_vis.head_mesh_opt.View[10].ColormapNumber) == 2 # default value
        assert len(result_vis.head_mesh.elmdata) == 6

        #Test for geo views
        assert int(result_vis.head_mesh_opt.View[11].ColormapNumber) == 5
        assert int(result_vis.head_mesh_opt.View[12].ColormapNumber) == 11
        assert int(result_vis.head_mesh_opt.View[13].ColormapNumber) == 6
        assert len(result_vis.head_mesh_opt.View) == 14


    def test_surface_rois(self, sphere3_msh: mesh_io.Msh, tmp_path):
        #Create Roi 1
        roi_1 = RegionOfInterest()
        surface_path = os.path.join(tmp_path, "surf.msh")
        sphere3_msh.crop_mesh(tags=[1003]).write(surface_path)
        roi_1.load_surfaces("custom", surface_path=surface_path)

        #Create Roi 2
        roi_2 = RegionOfInterest()
        surface_path = os.path.join(tmp_path, "surf2.msh")
        sphere3_msh.crop_mesh(tags=[1004]).write(surface_path)
        roi_2.load_surfaces("custom", surface_path=surface_path)

        #Create result mesh 1
        result_mesh_1_file_path = os.path.join(tmp_path, "test1_TMS_scalar.msh")
        result_mesh_1 = deepcopy(sphere3_msh)
        opt_1 = gmsh_view.Visualization(result_mesh_1)
        result_mesh_1.add_node_field(-np.arange(result_mesh_1.nodes.nr), "node_test_1")
        opt_1.add_view(ColormapNumber=1)
        result_mesh_1.add_node_field(
            -np.arange(result_mesh_1.nodes.nr) * 2, "node_test_2"
        )
        opt_1.add_view(ColormapNumber=2)
        result_mesh_1.add_element_field(np.arange(result_mesh_1.elm.nr), "elm_test_1")
        opt_1.add_view(ColormapNumber=3)
        result_mesh_1.add_element_field(
            np.arange(result_mesh_1.elm.nr) * 2, "elm_test_2"
        )
        opt_1.add_view(ColormapNumber=4)
        geo_file_path = os.path.join(tmp_path, "test1_TMS_coil_pos.geo")
        opt_1.add_merge(geo_file_path)
        mesh_io.write_geo_lines(
            [[1, 3, 4], [8, 9, 10]],
            [[5, 6, 7], [11, 12, 13]],
            geo_file_path,
            name="lines",
        )
        opt_1.add_view(ColormapNumber=5)
        mesh_io.write_geo_triangles(
            [[0, 1, 2]],
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
            geo_file_path,
            name="scalp",
            mode="ba",
        )
        opt_1.add_view(ColormapNumber=6)
        result_mesh_1.write(result_mesh_1_file_path)
        opt_1.write_opt(result_mesh_1_file_path)

        #Create result mesh 2
        result_mesh_2_file_path = os.path.join(tmp_path, "test2_TMS_scalar.msh")
        result_mesh_2 = deepcopy(sphere3_msh)
        opt_2 = gmsh_view.Visualization(result_mesh_2)
        result_mesh_2.add_node_field(
            -np.arange(result_mesh_2.nodes.nr) * 3, "node_test_1"
        )
        opt_2.add_view(ColormapNumber=7)
        result_mesh_2.add_node_field(
            -np.arange(result_mesh_2.nodes.nr) * 4, "node_test_2"
        )
        opt_2.add_view(ColormapNumber=8)
        result_mesh_2.add_element_field(
            np.arange(result_mesh_2.elm.nr) * 3, "elm_test_1"
        )
        opt_2.add_view(ColormapNumber=9)
        result_mesh_2.add_element_field(
            np.arange(result_mesh_2.elm.nr) * 4, "elm_test_2"
        )
        opt_2.add_view(ColormapNumber=10)
        geo_file_path = os.path.join(tmp_path, "test2_TMS_coil_pos.geo")
        opt_2.add_merge(geo_file_path)
        mesh_io.write_geo_spheres(
            [[1, 2, 3], [3, 4, 5]], geo_file_path, name="spheres", mode="ba"
        )
        opt_2.add_view(ColormapNumber=11)
        mesh_io.write_geo_triangles(
            [[0, 1, 2]],
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
            geo_file_path,
            name="scalp",
            mode="ba",
        )
        opt_2.add_view(ColormapNumber=12)
        result_mesh_2.write(result_mesh_2_file_path)
        opt_2.write_opt(result_mesh_2_file_path)

        result_vis = RoiResultVisualization(
            [roi_1, roi_2],
            [result_mesh_1_file_path, result_mesh_2_file_path],
            tmp_path,
            "test_vis",
            ["roi", "anti_roi"],
            ["mesh_1", "mesh_2"],
        )

        result_vis.get_surface_mesh().add_node_field(
            -np.arange(result_vis.get_surface_mesh().nodes.nr) * 10, "node_test_y"
        )

        result_vis.get_surface_mesh().add_element_field(
            np.arange(result_vis.get_surface_mesh().elm.nr) * 10, "elm_test_x"
        )

        result_vis.create_visualization()

        result_vis.surface_mesh_data_name_to_gmsh_view['mesh_1__elm_test_2'].VectorType = 6
        result_vis.remove_field_from_surface_mesh('mesh_2__elm_test_2')

        assert result_vis.head_mesh is None
        assert result_vis.surface_mesh is not None
        assert result_vis.surface_mesh_opt is not None
        assert isinstance(result_vis.surface_mesh_opt.View, list)

        #Test for the Roi node data fields
        assert result_vis.surface_mesh.nodedata[0].field_name == 'roi'
        assert int(result_vis.surface_mesh_opt.View[0].ColormapNumber) == 2 # default value
        assert int(result_vis.surface_mesh_opt.View[0].ShowScale) == 0
        assert result_vis.surface_mesh.nodedata[1].field_name == 'anti_roi'
        assert int(result_vis.surface_mesh_opt.View[1].ColormapNumber) == 2 # default value
        assert int(result_vis.surface_mesh_opt.View[1].ShowScale) == 0

        #Test for all node data fields 
        assert result_vis.surface_mesh.nodedata[2].field_name == 'mesh_1__node_test_1'
        assert int(result_vis.surface_mesh_opt.View[2].ColormapNumber) == 1
        assert result_vis.surface_mesh.nodedata[3].field_name == 'mesh_1__node_test_2'
        assert int(result_vis.surface_mesh_opt.View[3].ColormapNumber) == 2
        assert result_vis.surface_mesh.nodedata[4].field_name == 'mesh_1__elm_test_1'
        assert int(result_vis.surface_mesh_opt.View[4].ColormapNumber) == 3
        assert result_vis.surface_mesh.nodedata[5].field_name == 'mesh_1__elm_test_2'
        assert int(result_vis.surface_mesh_opt.View[5].VectorType) == 6
        assert int(result_vis.surface_mesh_opt.View[5].ColormapNumber) == 4

        assert result_vis.surface_mesh.nodedata[6].field_name == 'mesh_2__node_test_1'
        assert int(result_vis.surface_mesh_opt.View[6].ColormapNumber) == 7
        assert result_vis.surface_mesh.nodedata[7].field_name == 'mesh_2__node_test_2'
        assert int(result_vis.surface_mesh_opt.View[7].ColormapNumber) == 8
        assert result_vis.surface_mesh.nodedata[8].field_name == 'mesh_2__elm_test_1'
        assert int(result_vis.surface_mesh_opt.View[8].ColormapNumber) == 9
        assert result_vis.surface_mesh.nodedata[9].field_name == 'node_test_y'
        assert int(result_vis.surface_mesh_opt.View[9].ColormapNumber) == 2 # default value
        assert len(result_vis.surface_mesh.nodedata) == 10

        #Test for all element data fields
        assert result_vis.surface_mesh.elmdata[0].field_name == 'elm_test_x'
        assert int(result_vis.surface_mesh_opt.View[10].ColormapNumber) == 2 # default value
        assert len(result_vis.surface_mesh.elmdata) == 1

        #Test for geo views
        assert int(result_vis.surface_mesh_opt.View[11].ColormapNumber) == 5
        assert int(result_vis.surface_mesh_opt.View[12].ColormapNumber) == 11
        assert int(result_vis.surface_mesh_opt.View[13].ColormapNumber) == 6
        assert len(result_vis.surface_mesh_opt.View) == 14


    def test_volume_and_surface_rois(self, sphere3_msh: mesh_io.Msh, tmp_path):
        #Create Roi 1
        roi_1 = RegionOfInterest()
        roi_1.load_mesh(deepcopy(sphere3_msh))
        roi_1.apply_tissue_mask(1003, "intersection")

        #Create Roi 2
        roi_2 = RegionOfInterest()
        surface_path = os.path.join(tmp_path, "surf2.msh")
        sphere3_msh.crop_mesh(tags=[1004]).write(surface_path)
        roi_2.load_surfaces("custom", surface_path=surface_path)

        #Create result mesh 1
        result_mesh_1_file_path = os.path.join(tmp_path, "test1_TMS_scalar.msh")
        result_mesh_1 = deepcopy(sphere3_msh)
        opt_1 = gmsh_view.Visualization(result_mesh_1)
        result_mesh_1.add_node_field(-np.arange(result_mesh_1.nodes.nr), "node_test_1")
        opt_1.add_view(ColormapNumber=1)
        result_mesh_1.add_node_field(
            -np.arange(result_mesh_1.nodes.nr) * 2, "node_test_2"
        )
        opt_1.add_view(ColormapNumber=2)
        result_mesh_1.add_element_field(np.arange(result_mesh_1.elm.nr), "elm_test_1")
        opt_1.add_view(ColormapNumber=3)
        result_mesh_1.add_element_field(
            np.arange(result_mesh_1.elm.nr) * 2, "elm_test_2"
        )
        opt_1.add_view(ColormapNumber=4)
        geo_file_path = os.path.join(tmp_path, "test1_TMS_coil_pos.geo")
        opt_1.add_merge(geo_file_path)
        mesh_io.write_geo_lines(
            [[1, 3, 4], [8, 9, 10]],
            [[5, 6, 7], [11, 12, 13]],
            geo_file_path,
            name="lines",
        )
        opt_1.add_view(ColormapNumber=5)
        mesh_io.write_geo_triangles(
            [[0, 1, 2]],
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
            geo_file_path,
            name="scalp",
            mode="ba",
        )
        opt_1.add_view(ColormapNumber=6)
        result_mesh_1.write(result_mesh_1_file_path)
        opt_1.write_opt(result_mesh_1_file_path)

        #Create result mesh 2
        result_mesh_2_file_path = os.path.join(tmp_path, "test2_TMS_scalar.msh")
        result_mesh_2 = deepcopy(sphere3_msh)
        opt_2 = gmsh_view.Visualization(result_mesh_2)
        result_mesh_2.add_node_field(
            -np.arange(result_mesh_2.nodes.nr) * 3, "node_test_1"
        )
        opt_2.add_view(ColormapNumber=7)
        result_mesh_2.add_node_field(
            -np.arange(result_mesh_2.nodes.nr) * 4, "node_test_2"
        )
        opt_2.add_view(ColormapNumber=8)
        result_mesh_2.add_element_field(
            np.arange(result_mesh_2.elm.nr) * 3, "elm_test_1"
        )
        opt_2.add_view(ColormapNumber=9)
        result_mesh_2.add_element_field(
            np.arange(result_mesh_2.elm.nr) * 4, "elm_test_2"
        )
        opt_2.add_view(ColormapNumber=10)
        geo_file_path = os.path.join(tmp_path, "test2_TMS_coil_pos.geo")
        opt_2.add_merge(geo_file_path)
        mesh_io.write_geo_spheres(
            [[1, 2, 3], [3, 4, 5]], geo_file_path, name="spheres", mode="ba"
        )
        opt_2.add_view(ColormapNumber=11)
        mesh_io.write_geo_triangles(
            [[0, 1, 2]],
            [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
            geo_file_path,
            name="scalp",
            mode="ba",
        )
        opt_2.add_view(ColormapNumber=12)
        result_mesh_2.write(result_mesh_2_file_path)
        opt_2.write_opt(result_mesh_2_file_path)

        result_vis = RoiResultVisualization(
            [roi_1, roi_2],
            [result_mesh_1_file_path, result_mesh_2_file_path],
            tmp_path,
            "test_vis",
            ["roi", "anti_roi"],
            ["mesh_1", "mesh_2"],
        )

        result_vis.get_head_mesh().add_node_field(
            -np.arange(result_mesh_2.nodes.nr) * 10, "node_test_v"
        )

        result_vis.get_head_mesh().add_element_field(
            np.arange(result_mesh_2.elm.nr) * 10, "elm_test_w"
        )

        result_vis.get_surface_mesh().add_node_field(
            -np.arange(result_vis.get_surface_mesh().nodes.nr) * 10, "node_test_y"
        )

        result_vis.get_surface_mesh().add_element_field(
            np.arange(result_vis.get_surface_mesh().elm.nr) * 10, "elm_test_x"
        )

        result_vis.create_visualization()

        result_vis.head_mesh_data_name_to_gmsh_view['mesh_1__elm_test_1'].VectorType = 5
        result_vis.remove_field_from_head_mesh('mesh_2__elm_test_1')

        result_vis.surface_mesh_data_name_to_gmsh_view['mesh_1__elm_test_2'].VectorType = 6
        result_vis.remove_field_from_surface_mesh('mesh_2__elm_test_2')

        assert result_vis.head_mesh is not None
        assert result_vis.head_mesh_opt is not None
        assert isinstance(result_vis.head_mesh_opt.View, list)
        assert result_vis.surface_mesh is not None
        assert result_vis.surface_mesh_opt is not None
        assert isinstance(result_vis.surface_mesh_opt.View, list)

        #Test surface mesh
        #Test for the Roi node data fields
        assert result_vis.surface_mesh.nodedata[0].field_name == 'anti_roi'
        assert int(result_vis.surface_mesh_opt.View[0].ColormapNumber) == 2 # default value
        assert int(result_vis.surface_mesh_opt.View[0].ShowScale) == 0

        #Test for all node data fields 
        assert result_vis.surface_mesh.nodedata[1].field_name == 'mesh_1__node_test_1'
        assert int(result_vis.surface_mesh_opt.View[1].ColormapNumber) == 1
        assert result_vis.surface_mesh.nodedata[2].field_name == 'mesh_1__node_test_2'
        assert int(result_vis.surface_mesh_opt.View[2].ColormapNumber) == 2
        assert result_vis.surface_mesh.nodedata[3].field_name == 'mesh_1__elm_test_1'
        assert int(result_vis.surface_mesh_opt.View[3].ColormapNumber) == 3
        assert result_vis.surface_mesh.nodedata[4].field_name == 'mesh_1__elm_test_2'
        assert int(result_vis.surface_mesh_opt.View[4].VectorType) == 6
        assert int(result_vis.surface_mesh_opt.View[4].ColormapNumber) == 4

        assert result_vis.surface_mesh.nodedata[5].field_name == 'mesh_2__node_test_1'
        assert int(result_vis.surface_mesh_opt.View[5].ColormapNumber) == 7
        assert result_vis.surface_mesh.nodedata[6].field_name == 'mesh_2__node_test_2'
        assert int(result_vis.surface_mesh_opt.View[6].ColormapNumber) == 8
        assert result_vis.surface_mesh.nodedata[7].field_name == 'mesh_2__elm_test_1'
        assert int(result_vis.surface_mesh_opt.View[7].ColormapNumber) == 9
        assert result_vis.surface_mesh.nodedata[8].field_name == 'node_test_y'
        assert int(result_vis.surface_mesh_opt.View[8].ColormapNumber) == 2 # default value
        assert len(result_vis.surface_mesh.nodedata) == 9

        #Test for all element data fields
        assert result_vis.surface_mesh.elmdata[0].field_name == 'elm_test_x'
        assert int(result_vis.surface_mesh_opt.View[9].ColormapNumber) == 2 # default value
        assert len(result_vis.surface_mesh.elmdata) == 1

        #Test for geo views
        assert int(result_vis.surface_mesh_opt.View[10].ColormapNumber) == 5
        assert int(result_vis.surface_mesh_opt.View[11].ColormapNumber) == 11
        assert int(result_vis.surface_mesh_opt.View[12].ColormapNumber) == 6
        assert len(result_vis.surface_mesh_opt.View) == 13

        #Test volume mesh
        #Test for all node data fields 
        assert result_vis.head_mesh.nodedata[0].field_name == 'mesh_1__node_test_1'
        assert int(result_vis.head_mesh_opt.View[0].ColormapNumber) == 1
        assert result_vis.head_mesh.nodedata[1].field_name == 'mesh_1__node_test_2'
        assert int(result_vis.head_mesh_opt.View[1].ColormapNumber) == 2
        assert result_vis.head_mesh.nodedata[2].field_name == 'mesh_2__node_test_1'
        assert int(result_vis.head_mesh_opt.View[2].ColormapNumber) == 7
        assert result_vis.head_mesh.nodedata[3].field_name == 'mesh_2__node_test_2'
        assert int(result_vis.head_mesh_opt.View[3].ColormapNumber) == 8
        assert result_vis.head_mesh.nodedata[4].field_name == 'node_test_v'
        assert int(result_vis.head_mesh_opt.View[4].ColormapNumber) == 2 # default value
        assert len(result_vis.head_mesh.nodedata) == 5

        #Test for the Roi element data fields
        assert result_vis.head_mesh.elmdata[0].field_name == 'roi'
        assert int(result_vis.head_mesh_opt.View[5].ColormapNumber) == 2 # default value
        assert int(result_vis.head_mesh_opt.View[5].ShowScale) == 0
    
        #Test for all element data fields
        assert result_vis.head_mesh.elmdata[1].field_name == 'mesh_1__elm_test_1'
        assert int(result_vis.head_mesh_opt.View[6].VectorType) == 5
        assert int(result_vis.head_mesh_opt.View[6].ColormapNumber) == 3
        assert result_vis.head_mesh.elmdata[2].field_name == 'mesh_1__elm_test_2'
        assert int(result_vis.head_mesh_opt.View[7].ColormapNumber) == 4
        assert result_vis.head_mesh.elmdata[3].field_name == 'mesh_2__elm_test_2'
        assert int(result_vis.head_mesh_opt.View[8].ColormapNumber) == 10
        assert result_vis.head_mesh.elmdata[4].field_name == 'elm_test_w'
        assert int(result_vis.head_mesh_opt.View[9].ColormapNumber) == 2 # default value
        assert len(result_vis.head_mesh.elmdata) == 5

        #Test for geo views
        assert int(result_vis.head_mesh_opt.View[10].ColormapNumber) == 5
        assert int(result_vis.head_mesh_opt.View[11].ColormapNumber) == 11
        assert int(result_vis.head_mesh_opt.View[12].ColormapNumber) == 6
        assert len(result_vis.head_mesh_opt.View) == 13
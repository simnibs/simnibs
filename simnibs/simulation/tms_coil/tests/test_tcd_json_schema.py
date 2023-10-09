from copy import deepcopy
from typing import Any

import jsonschema
import pytest


class TestPassJsonSchema:
    def test_pass_minimal_coil(
        self, minimal_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        jsonschema.validate(minimal_tcd_coil_dict, tcd_json_schema)

    def test_pass_medium_coil(
        self, medium_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        jsonschema.validate(medium_tcd_coil_dict, tcd_json_schema)

    def test_pass_full_coil(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        jsonschema.validate(full_tcd_coil_dict, tcd_json_schema)


class TestFailJsonSchema:
    def test_missing_property_coilElementList(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property.pop("coilElementList")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_element_type(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["coilElementList"][0].pop("type")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_points_and_values(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["coilElementList"][0].pop("points")
        coil_missing_property["coilElementList"][0].pop("values")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_data(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["coilElementList"][2].pop("data")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_affine(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["coilElementList"][2].pop("affine")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_from_stimulator(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["stimulatorList"][0].pop("maxdIdt")
        coil_missing_property["stimulatorList"][0].pop("name")
        coil_missing_property["stimulatorList"][0].pop("brand")
        coil_missing_property["stimulatorList"][0].pop("waveformList")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_waveform_time(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["stimulatorList"][0]["waveformList"][0].pop("time")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_waveform_signal(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["stimulatorList"][0]["waveformList"][0].pop("signal")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_initial(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["deformRangeList"][0].pop("initial")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_deform_type(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["deformList"][0].pop("type")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_point1(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["deformList"][3].pop("point1")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_point2(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["deformList"][3].pop("point2")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_model_points(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["coilModels"][0].pop("points")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_missing_property_faces(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_missing_property = deepcopy(full_tcd_coil_dict)
        coil_missing_property["coilModels"][0].pop("faces")
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_missing_property, tcd_json_schema)

    def test_wrong_list_size_coilElementList(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilElementList"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_stimulatorList(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["stimulatorList"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_waveformList(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["stimulatorList"][0]["waveformList"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_waveformList_time(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["stimulatorList"][0]["waveformList"][0]["time"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_waveformList_signal(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["stimulatorList"][0]["waveformList"][0]["signal"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_waveformList_fit(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["stimulatorList"][0]["waveformList"][0]["fit"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_coilModels(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_limits(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["limits"] = [[255], [255], [255]]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["limits"] = [[255], [255], [255]]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["limits"] = [[0, 255, 555], [0, 255, 555], [0, 255, 555]]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["limits"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["limits"] = [
            [0, 255, 555],
            [0, 255, 555],
            [0, 255, 555],
            [0, 255, 555],
        ]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_resolution(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["resolution"] = [1, 1, 1, 1]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["resolution"] = [1, 1]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_deformList(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["deformList"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_deform_range(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["deformRangeList"][0]["range"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["deformRangeList"][0]["range"] = [0, 10, 20]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_deform_points(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["deformList"][3]["point1"] = [10, 20]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["deformList"][3]["point1"] = [10, 20, 20, 3]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["deformList"][3]["point2"] = [10, 20]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["deformList"][3]["point2"] = [10, 20, 20, 4]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_model_points(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["points"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["points"] = [[1, 2]]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["points"] = [[1, 2, 3, 4]]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_model_faces(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["faces"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["faces"] = [[1, 2]]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["faces"] = [[1, 2, 3, 4]]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_model_min_distance_points(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["minDistancePoints"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["minDistancePoints"] = [[1, 2]]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["minDistancePoints"] = [[1, 2, 3, 4]]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_list_size_model_intersect_points(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["intersectPoints"] = []
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["intersectPoints"] = [[1, 2]]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

        coil_wrong_list_sizes = deepcopy(full_tcd_coil_dict)
        coil_wrong_list_sizes["coilModels"][0]["intersectPoints"] = [[1, 2, 3, 4]]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_list_sizes, tcd_json_schema)

    def test_wrong_type_coil_casing(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilCasing"] = 1.2
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_wrong_type_stimulator(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilElementList"][0]["stimulator"] = 1.2
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_wrong_type_element_casing(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilElementList"][0]["elementCasing"] = 1.2
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_wrong_type_deformations(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilElementList"][0]["deformations"][0] = 1.2
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_wrong_type_faces(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilModels"][0]["faces"][0] = 1.2
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_wrong_type_element_type(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilElementList"][0]["type"] = 1.2
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_wrong_type_deform_type(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["deformList"][0]["type"] = 3
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_out_of_bounds_coil_casing(self, full_tcd_coil_dict, tcd_json_schema):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilCasing"] = -1
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_out_of_bounds_resolution(self, full_tcd_coil_dict, tcd_json_schema):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["resolution"] = [0, 0, 0]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["resolution"] = [-1, -1, -1]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_out_of_bounds_stimulator(self, full_tcd_coil_dict, tcd_json_schema):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilElementList"][0]["stimulator"] = -1
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_out_of_bounds_element_casing(self, full_tcd_coil_dict, tcd_json_schema):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilElementList"][0]["elementCasing"] = -1
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_out_of_bounds_deformations(self, full_tcd_coil_dict, tcd_json_schema):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilElementList"][0]["deformations"][0] = -1
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_out_of_bounds_faces(self, full_tcd_coil_dict, tcd_json_schema):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilModels"][0]["faces"][0] = -1
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_out_of_bounds_maxdIdt(self, full_tcd_coil_dict, tcd_json_schema):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["stimulatorList"][0]["maxdIdt"] = -1
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["stimulatorList"][0]["maxdIdt"] = 0
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_out_of_bounds_waveform_signal(self, full_tcd_coil_dict, tcd_json_schema):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["stimulatorList"][0]["waveformList"][0]["signal"] = [-1, 1]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["stimulatorList"][0]["waveformList"][0]["signal"] = [0, 2]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_out_of_bounds_waveform_fit(self, full_tcd_coil_dict, tcd_json_schema):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["stimulatorList"][0]["waveformList"][0]["fit"] = [0, 2]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["stimulatorList"][0]["waveformList"][0]["fit"] = [-1, 1]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_out_of_bounds_element_type(self, full_tcd_coil_dict, tcd_json_schema):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["coilElementList"][0]["type"] = 0
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_out_of_bounds_deform_type(self, full_tcd_coil_dict, tcd_json_schema):
        coil_wrong_types = deepcopy(full_tcd_coil_dict)
        coil_wrong_types["deformList"][0]["type"] = "w"
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_wrong_types, tcd_json_schema)

    def test_unsupported_properties(
        self, full_tcd_coil_dict: dict[str, Any], tcd_json_schema: Any
    ):
        coil_unsupported_property = deepcopy(full_tcd_coil_dict)
        coil_unsupported_property["namé"] = "Super NAME"
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_unsupported_property, tcd_json_schema)

        coil_unsupported_property = deepcopy(full_tcd_coil_dict)
        coil_unsupported_property["coilElementList"][0]["metaInfo"] = "Super meta info"
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_unsupported_property, tcd_json_schema)

        coil_unsupported_property = deepcopy(full_tcd_coil_dict)
        coil_unsupported_property["stimulatorList"][0][" "] = "Super secret"
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_unsupported_property, tcd_json_schema)

        coil_unsupported_property = deepcopy(full_tcd_coil_dict)
        coil_unsupported_property["stimulatorList"][0]["waveformList"][0][
            "type"
        ] = "There is no type"
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_unsupported_property, tcd_json_schema)

        coil_unsupported_property = deepcopy(full_tcd_coil_dict)
        coil_unsupported_property["deformList"][0]["min"] = 42
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_unsupported_property, tcd_json_schema)

        coil_unsupported_property = deepcopy(full_tcd_coil_dict)
        coil_unsupported_property["coilModels"][0]["name"] = [1, 2]
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(coil_unsupported_property, tcd_json_schema)

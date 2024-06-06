import os
import scipy

from simnibs.utils.matlab_read import matlab_field_to_list, matlab_struct_to_dict, matlab_sub_struct_to_matlab_struct, try_to_read_matlab_field


class TestTryToReadMatlabField:
    def test_read_types(self, tmp_path):
        test_dict = {"string": "Test", "int": 3, "float": 55.43, "bool": True}

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, test_dict)

        dict_loaded = scipy.io.loadmat(mat_path)

        assert test_dict["string"] == try_to_read_matlab_field(
            dict_loaded, "string", str
        )
        assert test_dict["int"] == try_to_read_matlab_field(dict_loaded, "int", int)
        assert test_dict["float"] == try_to_read_matlab_field(
            dict_loaded, "float", float
        )
        assert test_dict["bool"] == try_to_read_matlab_field(dict_loaded, "bool", bool)


class TestMatlabStructToDict:
    def test_read_nested(self, tmp_path):
        test_dict = {
            "dict": {
                "string": "Test",
                "int": 3,
                "float": 55.43,
                "bool": True,
                "dict": {"string": "Test2", "int": 6, "float": 65.43, "bool": False},
            }
        }

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, test_dict)

        dict_loaded = matlab_struct_to_dict(scipy.io.loadmat(mat_path)['dict'])

        assert test_dict["dict"] == dict_loaded

class TestMatlabSubStructToMatlabStruct:
    def test_simple(self, tmp_path):
        test_dict = {
            "dict": {
                "string": "Test",
                "int": 3,
                "float": 55.43,
                "bool": True,
                "dict": {"string": "Test2", "int": 6, "float": 65.43, "bool": False},
            }
        }

        mat_path = os.path.join(tmp_path, "test.mat")
        mat_2_path = os.path.join(tmp_path, "test2.mat")
        scipy.io.savemat(mat_path, test_dict)
        scipy.io.savemat(mat_2_path, test_dict["dict"])

        dict_loaded = matlab_sub_struct_to_matlab_struct(scipy.io.loadmat(mat_path)['dict'])
        dict_2_loaded = scipy.io.loadmat(mat_2_path)

        del dict_2_loaded['__globals__']
        del dict_2_loaded['__header__']
        del dict_2_loaded['__version__']

        assert dict_loaded == dict_2_loaded

class TestMatlabFieldToList:
    def test_number_lists(self, tmp_path):
        test_dict = {
            'ints': [1,2,3,4,5],
            'nested_ints': [[1,2], [3,4]],
            'double_nested_ints': [[[1,2], [3,4]], [[5,6], [7,8]]],
            'floats': [1.2, 3.4, 5.6, 7.8],
            'nested_floats': [[1.2, 3.4], [5.6, 7.8]],
            'double_nested_floats': [[[1.2, 3.4], [5.6, 7.8]], [[9.1, 11.12], [13.14, 15.16]]]
        }

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, test_dict)

        dict_loaded = scipy.io.loadmat(mat_path)

        assert matlab_field_to_list(dict_loaded, 'ints', 1) == test_dict["ints"]
        assert matlab_field_to_list(dict_loaded, 'nested_ints', 2) == test_dict["nested_ints"]
        assert matlab_field_to_list(dict_loaded, 'double_nested_ints', 3) == test_dict["double_nested_ints"]

        assert matlab_field_to_list(dict_loaded, 'floats', 1) == test_dict["floats"]
        assert matlab_field_to_list(dict_loaded, 'nested_floats', 2) == test_dict["nested_floats"]        
        assert matlab_field_to_list(dict_loaded, 'double_nested_floats', 3) == test_dict["double_nested_floats"]

    def test_string_lists(self, tmp_path):
        test_dict = {
            'strings': ['a', 'bc', 'def', 'ghij'],
            'nested_strings': [['a', 'bc'], ['def', 'ghij']]
        }

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, test_dict)

        dict_loaded = scipy.io.loadmat(mat_path)

        assert matlab_field_to_list(dict_loaded, 'strings', 2) == test_dict["strings"]
        assert matlab_field_to_list(dict_loaded, 'nested_strings', 2) == test_dict["nested_strings"]
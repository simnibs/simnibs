import os
import scipy

from simnibs.utils.matlab_read import dict_from_matlab, try_to_read_matlab_field


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
        dict_loaded = scipy.io.loadmat(mat_path)

        assert test_dict["dict"] == dict_from_matlab(dict_loaded["dict"])

        del dict_loaded["__globals__"]
        del dict_loaded["__header__"]
        del dict_loaded["__version__"]
        assert test_dict == dict_from_matlab(dict_loaded)

    def test_simple_number_lists(self, tmp_path):
        test_dict = {
            "ints": [1, 2, 3, 4, 5],
            "nested_ints": [[1, 2], [3, 4]],
            "double_nested_ints": [[[1, 2], [3, 4]], [[5, 6], [7, 8]]],
            "floats": [1.2, 3.4, 5.6, 7.8],
            "nested_floats": [[1.2, 3.4], [5.6, 7.8]],
            "double_nested_floats": [
                [[1.2, 3.4], [5.6, 7.8]],
                [[9.1, 11.12], [13.14, 15.16]],
            ],
        }

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, test_dict)

        dict_loaded = scipy.io.loadmat(mat_path)
        del dict_loaded["__globals__"]
        del dict_loaded["__header__"]
        del dict_loaded["__version__"]

        assert test_dict == dict_from_matlab(dict_loaded)

    def test_reduce_lists(self, tmp_path):
        test_dict = {
            "ints": [1],
            "nested_ints": [[1], [3]],
            "double_nested_ints": [[[1], [3]], [[5], [7]]],
            "floats": [1.2],
            "nested_floats": [[1.2], [5.6]],
            "double_nested_floats": [
                [[1.2], [5.6]],
                [[9.1], [13.14]],
            ],
        }

        expected_result_dict = {
            "ints": 1,
            "nested_ints": [1, 3],
            "double_nested_ints": [[1, 3], [5, 7]],
            "floats": 1.2,
            "nested_floats": [1.2, 5.6],
            "double_nested_floats": [
                [1.2, 5.6],
                [9.1, 13.14],
            ],
        }

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, test_dict)

        dict_loaded = scipy.io.loadmat(mat_path)
        del dict_loaded["__globals__"]
        del dict_loaded["__header__"]
        del dict_loaded["__version__"]

        assert expected_result_dict == dict_from_matlab(dict_loaded)

    def test_nested_number_lists(self, tmp_path):
        test_dict = {
            "ints": [1, 2, 3, 4, 5],
            "first": {
                "nested_ints": [[1, 2], [3, 4]],
                "double_nested_ints": [[[1, 2], [3, 4]], [[5, 6], [7, 8]]],
                "second": {
                    "floats": [1.2, 3.4, 5.6, 7.8],
                    "nested_floats": [[1.2, 3.4], [5.6, 7.8]],
                    "double_nested_floats": [
                        [[1.2, 3.4], [5.6, 7.8]],
                        [[9.1, 11.12], [13.14, 15.16]],
                    ],
                },
            },
        }

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, test_dict)

        dict_loaded = scipy.io.loadmat(mat_path)
        del dict_loaded["__globals__"]
        del dict_loaded["__header__"]
        del dict_loaded["__version__"]

        assert test_dict == dict_from_matlab(dict_loaded)

    def test_simple_string_lists(self, tmp_path):
        test_dict = {
            "strings": ["a", "bc", "def", "ghij"],
            "nested_strings": [["a", "bc"], ["def", "ghij"]],
        }

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, test_dict)

        dict_loaded = scipy.io.loadmat(mat_path)
        del dict_loaded["__globals__"]
        del dict_loaded["__header__"]
        del dict_loaded["__version__"]

        assert test_dict == dict_from_matlab(dict_loaded)

    def test_nested_string_lists(self, tmp_path):
        test_dict = {
            "strings": ["a", "bc", "def", "ghij"],
            "nested_strings": [["a", "bc"], ["def", "ghij"]],
            "first":{
                "strings1": ["a1", "bc1", "def1", "ghij1"],
                "nested_strings1": [["a1", "bc1"], ["def1", "ghij1"]],
                "second":{
                    "strings2": ["a2", "bc2", "def2", "ghij2"],
                    "nested_strings1": [["a2", "bc2"], ["def2", "ghij2"]],
                }
            }
        }

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, test_dict)

        dict_loaded = scipy.io.loadmat(mat_path)
        del dict_loaded["__globals__"]
        del dict_loaded["__header__"]
        del dict_loaded["__version__"]

        assert test_dict == dict_from_matlab(dict_loaded)


    def test_list_of_dicts(self, tmp_path):
        test_dict = {
            'test': [{"key1": 1}, {"key2": 2}, {"key3": 3}],
            'test2':{
                'test3': [{"key4": 4}, {"key5": 5}, {"key6": 6}]
            }
        }

        mat_path = os.path.join(tmp_path, "test.mat")
        scipy.io.savemat(mat_path, test_dict)

        dict_loaded = scipy.io.loadmat(mat_path)
        del dict_loaded["__globals__"]
        del dict_loaded["__header__"]
        del dict_loaded["__version__"]
        print(dict_loaded)
        assert test_dict == dict_from_matlab(dict_loaded)
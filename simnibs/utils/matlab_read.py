import os
import numpy as np
import scipy.io


def try_to_read_matlab_field(matlab_structure, field_name, field_type, alternative=None):
    """
    Function for flexibilly reading a field from a matlab .mat file
    Tries to read the field with the specified name
    if successful, returns the value read if not, returns the alternative

    Parameters
    ----------
    matlab_structure: dict
        matlab structure as read by scipy.io, without squeeze
    field_name: str
        name of field in mat structure
    field_type: function
        function that transforms the field into the desired type
        'int', 'float', 'str', ...
    alternative: any
        if the field could not be read, return alternative. (Default: None)
    """

    try:
        return field_type(matlab_structure[field_name][0])
    except (TypeError, KeyError, IndexError, ValueError):
        pass
    try:
        return field_type(matlab_structure[field_name][0][0])
    except (TypeError, KeyError, IndexError, ValueError):
        pass
    return alternative

def dict_from_matlab(matlab_structure):
    """Turns a matlab structure into a normal python dict

    Limitations
    -----------
    1. All lists are reduced as much as possible
        Example: [1] -> 1, [[1], [1]] -> [1, 1]
    2. Empty elements will be removed:
        Example: [] -> del, [[None]] -> del, [[], []] -> del
    3. Leading and trailing spaces will be removed from strings
        Example: "  hello" -> "hello", "buzz    " -> "buz"
    4. Only lists with a single type are allowed
        Example: [1, "a"] -> Error
    Parameters
    ----------
    matlab_structure : dict
        matlab sub structure as read by scipy.io

    Returns
    -------
    dict
        the dictionary containing the information from the matlab structure
    """
    result = {}
    # get keys from numpy array or dict
    if isinstance(matlab_structure, dict):
        keys = matlab_structure.keys()
    else:
        keys = matlab_structure.dtype.names

    for name in keys:
        value = matlab_structure[name]
        if isinstance(value, np.ndarray):
            #unpack structure if necessary
            if value.size == 1 and isinstance(value[0], np.ndarray) and np.issubdtype(value[0][0].dtype, np.void):
                value = value[0][0]
            
            #handle list of dicts
            if value.size > 1 and isinstance(value[0], np.ndarray) and np.issubdtype(value[0][0].dtype, np.void):
                result[name] = []
                for elm in value[0]:
                    result[name].append(dict_from_matlab(elm))
                continue

            #create dict from structure or reduce list
            if np.issubdtype(value.dtype, np.void):
                result[name] = dict_from_matlab(value)
            else:
                reduced_array = _reduce_array(value)
                if reduced_array is not None:
                    result[name] = reduced_array
        else:
            result[name] = value

        #remove empty dicts
        if name in result and isinstance(result[name], dict):
            if len(result[name]) == 0:
                del result[name]
    return result

def _reduce_array(array):
    # double reduce for arrays inside of lists
    result = np.squeeze(np.array(np.squeeze(array).tolist())).tolist()

    # turn tuple into lists before reducing again
    if isinstance(result, tuple):
        result = _reduce_array(list(result))

    # if array hat object type, reduce all sub arrays
    if isinstance(result, list):
        for i, elm in enumerate(result):
            if isinstance(elm, np.ndarray):
                result[i] = _reduce_array(elm)

        if len(result) == 0:
            return None
        
    if np.issubdtype(np.array(result).dtype, np.str_):
        result = strip_strings(result)

    return result

def strip_strings(input):
    if isinstance(input, str):
        return input.strip()
    elif isinstance(input, list):
        return [strip_strings(item) for item in input]
    else:
        return input

def remove_None(src):
    '''Substitutes None by an empty char array '''
    if src is None:
        src = ''
    return src


def read_mat(fn):
    ''' Reads a ".mat" file

    Parameters
    ----------------
    fn: str
        File name

    Returns
    -----------
    struct:
        Structure defined by .mat file
    '''
    if not os.path.isfile(fn):
        raise IOError('File: {0} does not exist'.format(fn))
    if os.path.splitext(fn)[1] != '.mat':
        raise IOError('SimNIBS only accepts matlab ".mat" files')
    try:
        mat = scipy.io.loadmat(fn, struct_as_record=True, squeeze_me=False)
    except:
        raise IOError("Could not open file. It was possibly saved with -v7.3")
    try:
        structure_type = mat['type'][0]
    except:
        try:
            keys = [k for k in mat.keys() if not k.startswith('__')]
            if len(keys) > 1:
                raise IOError(
                    'Could not open .mat file. Nested structure?')
            structure_type = mat[keys[0]][0]['type'][0][0]
            mat = mat[keys[0]][0][0]
        except:
            raise IOError(
                "Could not access structure type in this .mat file")

    if structure_type.lower() == 'session':
        from ..simulation.sim_struct import SESSION
        structure = SESSION(matlab_struct=mat)
    elif structure_type.lower() == 'tdcsleadfield':
        from ..simulation.sim_struct import TDCSLEADFIELD
        structure = TDCSLEADFIELD(matlab_struct=mat)
    elif structure_type.lower() == 'tmsoptimize':
        from ..optimization.opt_struct import TMSoptimize
        structure = TMSoptimize.read_mat_struct(mat)
    elif structure_type.lower() == 'tdcsoptimize':
        from ..optimization.opt_struct import TDCSoptimize
        structure = TDCSoptimize.read_mat_struct(mat)
    elif structure_type.lower() == 'tdcsdistributedoptimize':
        from ..optimization.opt_struct import TDCSDistributedOptimize
        structure = TDCSDistributedOptimize.read_mat_struct(mat)
    elif structure_type.lower() == 'tmsflexoptimization':
        from ..optimization.opt_struct import TmsFlexOptimization
        structure = TmsFlexOptimization(dict_from_matlab(mat))
    elif structure_type.lower() == 'tesflexoptimization':
        from ..optimization.opt_struct import TesFlexOptimization
        structure = TesFlexOptimization(dict_from_matlab(mat))
    elif structure_type.lower() == 'regionofinterest':
        from .region_of_interest import RegionOfInterest
        structure = RegionOfInterest(dict_from_matlab(mat))
    else:
        raise IOError('Not a valid structure type!')

    return structure
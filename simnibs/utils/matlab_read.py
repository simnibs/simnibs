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

def matlab_struct_to_dict(matlab_structure):
    """Turns a matlab structure (keys described in dtype) into a normal python dict

    Parameters
    ----------
    matlab_structure : dict
        matlab sub structure as read by scipy.io

    Returns
    -------
    dict
        the dictionary containing the information from the matlab sub structure
    """
    result = {}
    if matlab_structure.dtype.names is None:
        return {}
    for name in matlab_structure.dtype.names:
        if isinstance(matlab_structure[0][name][0].dtype, np.dtypes.VoidDType):
            result[name] = matlab_struct_to_dict(matlab_structure[0][name][0])
        else:
            if matlab_structure[0][name][0].size == 0:
                result[name] = None
            elif isinstance(matlab_structure[0][name][0][0], np.ndarray):
                result[name] = matlab_structure[0][name][0][0][0]
            else:
                result[name] = matlab_structure[0][name][0][0]

    return result

def matlab_sub_struct_to_matlab_struct(matlab_structure):
    """Turns a matlab sub structure (dtype contains key names) into a matlab structure (dict)

    Parameters
    ----------
    matlab_structure : dict
        matlab sub structure as read by scipy.io

    Returns
    -------
    dict
        matlab structure as if it would have been directly read by scipy.io
    """
    result = {}
    if matlab_structure.dtype.names is None:
        return {}
    for name in matlab_structure.dtype.names:
        result[name] = matlab_structure[0][name][0]
    return result

def matlab_field_to_list(matlab_structure, field_name, dim) -> list | None:
    """Returns a matlab field as a list with dim as the number of dimensions

    Parameters
    ----------
    matlab_structure : dict
        matlab structure as read by scipy.io, without squeeze
    field_name: str
        name of field in mat structure
    dim : int
        Number of dimensions that the list has

    Returns
    -------
    list | None
        The list with the number of dimensions specified
    """
    if field_name not in matlab_structure.dtype.names or matlab_structure[field_name].size == 0:
        return None
    if np.issubdtype(matlab_structure[field_name].dtype, np.str_):
        nested = np.array(matlab_structure[field_name].tolist())
        return np.array([string.strip() for string in nested.flatten()]).reshape(nested.shape).tolist()
    elif np.issubdtype(matlab_structure[field_name].dtype, 'O'):
        field_as_list = matlab_structure[field_name].flatten().tolist()
        for idx, subfield in enumerate(field_as_list):
            if np.issubdtype(subfield.dtype, np.str_):
                field_as_list[idx] = subfield.flatten()[0]
            else:
                field_as_list[idx] = subfield.flatten().tolist()
        return field_as_list
    else:
        if dim > 1:
            return matlab_structure[field_name].tolist()
        else:
            return matlab_structure[field_name].tolist()[0]

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
        structure = TmsFlexOptimization()
        structure.read_mat_struct(mat)
    elif structure_type.lower() == 'regionofinterest':
        from .region_of_interest import RegionOfInterest
        structure = RegionOfInterest()
        structure.read_mat_struct(mat)

    else:
        raise IOError('Not a valid structure type!')

    return structure

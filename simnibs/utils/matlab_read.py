import os
import scipy.io
import numpy as np


def try_to_read_matlab_field(matlab_structure, field_name, field_type, alternative=None):
    """
    Function for flexibilly reading a field from the mesh file
    Tries to read the field with the specified name
    if sucesseful, returns the read
    if not, returns the alternative

    Parameters
    ----------
    matlab_structure: dict
        matlab structure as read by scipy.io, without squeeze
    field_name: str
        name of field in mat structure
    field_type: function
        function that transforms the field into the desired type
        'int', 'float', 'str',....
    alternative: any
        if the field could not be read, return alternative. (Default: None)
    """
    if field_type == list:
        try:
            res = np.squeeze(matlab_structure[field_name]).tolist()
            try:
                # lists of string are white-padded to largest element in list. Remove padding:
                res = [re.strip() for re in res]
            except AttributeError:
                pass
            return res
        except (TypeError, KeyError, IndexError, ValueError):
            pass
    elif field_type == np.array:
        try:
            return np.squeeze(np.array(matlab_structure[field_name]))
        except (TypeError, KeyError, IndexError, ValueError):
            pass
    else:
        try:
            return field_type(matlab_structure[field_name][0])
        except (TypeError, KeyError, IndexError, ValueError):
            pass
        try:
            return field_type(matlab_structure[field_name][0][0])
        except (TypeError, KeyError, IndexError, ValueError):
            pass
    return alternative


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

    if structure_type == 'SESSION':
        from ..simulation.sim_struct import SESSION
        structure = SESSION(matlab_struct=mat)
    elif structure_type == 'TDCSLEADFIELD':
        from ..simulation.sim_struct import TDCSLEADFIELD
        structure = TDCSLEADFIELD(matlab_struct=mat)
    elif structure_type == 'TDCSoptimize':
        from ..optimize.optimization import TDCSoptimize
        structure = TDCSoptimize.read_mat_struct(mat)
    else:
        raise IOError('Not a valid structure type!')

    return structure

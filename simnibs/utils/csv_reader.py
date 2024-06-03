import os
import warnings
import csv

import numpy as np

from .file_finder import SubjectFiles

def read_csv_positions(fn):
    ''' Reads positions from a .csv file

    Parameters
    ------------
    fn: str
        Name of csv file

    Returns
    --------
    type: list
        Type of position in each row ('Generic', 'Fiducial', 'Electrode', 'ReferenceElectrode' or
        'Coil')

    coordinates: list
        Coordinates in each row

    extra: list
        extra coordinates in each row (eg: electrode or coil axes)

    name: list
        Name of position

    extra_cols: list
        Any extra information stored in the columns

    header: str
        Any information in the header

    '''
    with open(os.path.expanduser(fn), 'r') as f:
        reader = csv.reader(f)
        rows = [row for row in reader]

    if len(rows[-1]) < 3:
        raise IOError('CSV file must have at least 4 rows')

    coordinates = []
    try:
        float(rows[0][1])
        header = []
        start = 0
    except:
        header = rows[0]
        start = 1

    rows = rows[start:]
    type_ = [r[0] for r in rows]

    extra = []
    name = []
    extra_cols = []
    type_filtered = []
    rows_filtered = []

    for t, r, i in zip(type_, rows, range(len(type_))):
        if t in ['Generic', 'Fiducial'] or len(r) == 4:
            name += [r[4] if len(r) >= 5 else None]
            extra_cols += [r[5:] if len(r) > 5 else None]
            extra += [None]
            type_filtered.append(t)
            rows_filtered.append(r)
        elif t in ['Electrode', 'ReferenceElectrode']:
            try:
                extra_ = np.array(r[4:7], float)
                assert len(extra_) == 3
                extra += [extra_]
                name += [r[7] if len(r) >= 8 else None]
                extra_cols += [r[8:] if len(r) > 8 else None]
            except:
                extra += [None]
                name += [r[4] if len(r) >= 5 else None]
                extra_cols += [r[5:] if len(r) > 5 else None]
            type_filtered.append(t)
            rows_filtered.append(r)
        elif t == 'CoilPos':
            extra += [np.array([float(d) for d in r[4:11]])]
            name += [r[11] if len(r) >= 12 else None]
            extra_cols += [r[12:] if len(r) > 12 else None]
            type_filtered.append(t)
            rows_filtered.append(r)
        else:
            warnings.warn('Unrecognized column type: {0}'.format(t))
    type_ = type_filtered
    rows = rows_filtered
    try:
        coordinates = np.array(
            [[float(d) for d in r[1:4]] for r in rows],
            dtype=float)
    except:
        raise IOError('Could not read coordinates from CSV file')

    return type_, coordinates, extra, name, extra_cols, header



def write_csv_positions(filename, types, coordinates, name, extra=None, extra_cols=None, header=None):
    ''' Write positions to a .csv file

    Parameters
    ------------
    fn: str
        Name of csv file
    type: list
        Type of position in each row ('Generic', 'Fiducial', 'Electrode', 'ReferenceElectrode' or
        'Coil')
    coordinates: numpy array
        Coordinates in each row
    name: list
        Name of position
    extra: list
        extra coordinates in each row (eg: electrode or coil axes)
    extra_cols: list
        Any extra information stored in the columns
    header: str
        Any information in the header

    '''
    n = len(types)

    coordinates = coordinates.tolist()

    name = [[n] if n else [] for n in name]

    if extra is None:
        extra = [None]*n
    extra = [[] if e is None else e.tolist() for e in extra]

    if extra_cols is None or len(extra_cols) == 0:
        extra_cols = [None]*n
    extra_cols = [e_c or [] for e_c in extra_cols]

    if header is None:
        header = []

    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        if header != []:
            writer.writerow(header)
        for t, c, e, n, e_c in zip(types, coordinates, extra, name, extra_cols):
            writer.writerow([t] + c + e + n + e_c)

def _get_eeg_positions(fn_csv):
    if not os.path.isfile(fn_csv):
        raise IOError('Could not find EEG cap file: {0}'.format(fn_csv))
    type_, coordinates, _, name, _, _ = read_csv_positions(fn_csv)
    eeg_pos = {}
    for i, t in enumerate(type_):
        if t in ['Electrode', 'ReferenceElectrode', 'Fiducial']:
            eeg_pos[name[i]] = coordinates[i]
    return eeg_pos


def eeg_positions(m2m_folder, cap_name='EEG10-10_UI_Jurak_2007.csv'):
    ''' Returns a directory with EEG electrode positions

    Parameters
    -----------

    m2m_folder: str
        Path to the m2m_{subject_id} folder, generated during the segmantation

    cap_name: str
        Name of EEG cap. Default: 'EEG10-10_UI_Jurak_2007.csv'


    Returns
    --------
    eeg_caps: dict
        Dictionary with cap position
    '''
    sub_files = SubjectFiles(subpath=m2m_folder)
    if not cap_name.endswith('.csv'):
        cap_name += '.csv'
    fn_cap = sub_files.get_eeg_cap(cap_name)
    return _get_eeg_positions(fn_cap)


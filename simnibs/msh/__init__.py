'''
SimNIBS mesh handling functions

'''

__all__ = [
    'read_msh',
    'write_msh',
    'read_freesurfer_surface',
    'write_freesurfer_surface',
    'read_gifti_surface',
    'read_curv',
    'write_curv',
    'read_stl',
    'write_geo_spheres',
    'write_geo_text',
    'Msh',
    'Nodes',
    'Elements',
    'ElementData',
    'subject2mni_coords',
    'mni2subject_coords',
    'eeg_positions'
]
from .mesh_io import *
from .transformations import *

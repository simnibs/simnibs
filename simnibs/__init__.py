import os
import sys
SIMNIBSDIR = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

from ._version import __version__
from .mesh_tools import mesh_io
from .mesh_tools.mesh_io import read_msh
from .utils import transformations, file_finder, region_of_interest
from .utils.transformations import mni2subject_coilpos, mni2subject_coords, subject2mni_coords, subject_atlas
from .utils.file_finder import *
from .utils.nnav import localite, softaxic, brainsight, ant
from .utils.mesh_element_properties import ElementTags
from .simulation import sim_struct
from .simulation import fem, biot_savart
from .utils.region_of_interest import *
from .simulation.run_simnibs import run_simnibs
from .optimization import opt_struct

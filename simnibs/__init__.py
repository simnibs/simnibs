import os
SIMNIBSDIR = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
from ._version import __version__
from .msh import *
from .transformations import *
from .utils import file_finder
from .utils.file_finder import *
from .simulation import sim_struct
from .simulation import cond
from .simulation.run_simnibs import run_simnibs
from .optimization import opt_struct

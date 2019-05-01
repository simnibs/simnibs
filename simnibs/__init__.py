import os
SIMNIBSDIR = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
from ._version import __version__
from .utils import file_finder
from .simulation import sim_struct
from .simulation import cond
from .simulation.run_simnibs import run_simnibs
try:
    from .internal import optimize
except ImportError:
    pass

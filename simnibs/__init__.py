import os
import sys
SIMNIBSDIR = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
# Add the external/lib/win folder to system path so that the dlls there can be found
if sys.platform == 'win32':
    os.environ['PATH'] = os.pathsep.join([
        os.path.join(SIMNIBSDIR, 'external', 'lib', 'win'),
        os.environ['PATH']
    ])
    # windows only: starting from python 3.8, the ctypes module only searches
    # for DLLs in certain (trusted) locations, hence we need to add this
    # explicitly to pick up the runtime dependencies for petsc
    os.add_dll_directory(os.path.join(SIMNIBSDIR, "external", "lib", "win"))
elif sys.platform == 'linux':
    os.environ['PATH'] = os.pathsep.join([
        os.path.join(SIMNIBSDIR, 'external', 'lib', 'linux'),
        os.environ['PATH']
    ])
elif sys.platform == 'darwin':
    os.environ['PATH'] = os.pathsep.join([
        os.path.join(SIMNIBSDIR, 'external', 'lib', 'osx'),
        os.environ['PATH']
    ])
from ._version import __version__
from .mesh_tools import mesh_io
from .utils import transformations
from .utils.transformations import *
from .utils import file_finder
from .utils.file_finder import *
from .utils.nnav import localite, softaxic, brainsight, ant
from .utils.mesh_element_properties import ElementTags
from .simulation import sim_struct
from .simulation import fem
from .utils.region_of_interest import *
from .simulation.run_simnibs import run_simnibs
from .optimization import opt_struct

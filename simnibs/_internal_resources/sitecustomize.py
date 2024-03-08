## Begining of SimNIBS sitecustomize.py ##
''' SimNIBS sitecustomize configuration
    This file is meant to be copied to the site-packages folder in Windows, as part of the postinstall process
    It allows for dynamically modifying the PATH to add the paths to the conda libraries
'''
import os

prefix = os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..'))
# Notice: this is where conda stores DLLs. If that changes in the future, this here should also change
# From conda-pack/conda_pack/scripts/windows/activate.bat
os.environ["PATH"] = os.pathsep.join([
    os.path.join(prefix, 'Library', 'mingw-w64', 'bin'),
    os.path.join(prefix, 'Library', 'usr', 'bin'),
    os.path.join(prefix, 'Library', 'bin')
]) + os.pathsep + os.environ["PATH"]
## End of SimNIBS sitecustomize.py ##
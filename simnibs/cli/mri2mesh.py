'''
    Simple wraper to call mri2mesh
'''

import sys
import os
import subprocess
from simnibs import SIMNIBSDIR


def main():
    if sys.platform == 'win32':
        raise OSError('mri2mesh does not run on windows!')

    subprocess.call(
        [os.path.join(SIMNIBSDIR, 'bash_modules', 'mri2mesh')] +
        sys.argv[1:])
    
if __name__ == '__main__':
    main()

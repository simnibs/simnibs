'''
    Simple wraper to call dwi2cond
'''

import sys
import os
import subprocess
from simnibs import SIMNIBSDIR


def main():
    if sys.platform == 'win32':
        raise OSError('dwi2cond does not run on windows!')

    subprocess.call(
        [os.path.join(SIMNIBSDIR, 'bash_modules', 'dwi2cond')] +
        sys.argv[1:])
    
if __name__ == '__main__':
    main()

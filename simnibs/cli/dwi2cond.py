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

    ret = subprocess.run(
        [os.path.join(SIMNIBSDIR, 'bash_modules', 'external', 'dwi2cond')] +
        sys.argv[1:])

    sys.exit(ret.returncode)

if __name__ == '__main__':
    main()

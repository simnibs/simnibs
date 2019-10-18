'''
    Simple wraper to call gmsh
'''

import sys
import subprocess
from simnibs.utils.file_finder import path2bin


def main():
    ret = subprocess.run(
        path2bin('gmsh') + ' '  + ' '.join(sys.argv[1:]),
        shell=True,
    )
    sys.exit(ret.returncode)

if __name__ == '__main__':
    main()

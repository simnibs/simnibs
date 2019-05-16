'''
    Simple wraper to call gmsh
'''

import sys
import subprocess
from simnibs.utils.file_finder import path2bin


def main():
    gmsh = path2bin('gmsh')
    subprocess.call(gmsh + ' ' + ' '.join(sys.argv[1:]), shell=True)


if __name__ == '__main__':
    main()

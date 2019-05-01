'''
    Simple wraper to call gmsh
'''

import sys
import subprocess
from simnibs.utils.file_finder import path2bin


def main():
    subprocess.call([path2bin('gmsh')] + sys.argv[1:])


if __name__ == '__main__':
    main()

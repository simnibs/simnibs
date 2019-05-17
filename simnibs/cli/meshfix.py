'''
    Simple wraper to call meshfix
'''

import sys
import subprocess
from simnibs.utils.file_finder import path2bin


def main():
    subprocess.call([path2bin('meshfix')] + sys.argv[1:])


if __name__ == '__main__':
    main()

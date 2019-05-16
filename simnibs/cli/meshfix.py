'''
    Simple wraper to call meshfix
'''

import sys
import subprocess
from simnibs.utils.file_finder import path2bin


def main():
    subprocess.call(path2bin('meshfix') + ' ' + ' '.join(sys.argv[1:]), shell=True)


if __name__ == '__main__':
    main()

'''
    Simple wrapper script ro call headreco
'''

import sys
from simnibs.msh.headreco import headmodel


def main():
    headmodel(sys.argv)

if __name__ == '__main__':
    main()

'''
    Simple wrapper script ro call headreco
'''

import sys
from simnibs.mesh_tools.headreco import headmodel


def main():
    headmodel(sys.argv)

if __name__ == '__main__':
    main()

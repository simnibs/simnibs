# -*- coding: utf-8 -*-
'''
    command line tool to merge Cerebellum, WM and venctricle labels to a single WM label
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2018 Andre Antunes, Guilherme B Saturnino

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
from __future__ import print_function
import simnibs.msh.mesh_io as mesh_io
from simnibs.utils.simnibs_logger import logger
import sys
import os


def usage():
    print("correct usage:")
    print("merge_labels_WM <fn.msh> <fn_out.msh>")

def main():
    if len(sys.argv) not in [3, 5]:
        print("Incorrect number of arguments")
        usage()
        sys.exit(1)


    # check extension
    if os.path.splitext(sys.argv[1])[1] != '.msh':
        print('1st argument must have a .msh extension')
        sys.exit(1)

    if os.path.splitext(sys.argv[2])[1] != '.msh':
        print('2nd argument must have a .msh extension')
        sys.exit(1)

    m = mesh_io.read_msh(sys.argv[1])

    m.fn = sys.argv[2]
    logger.info("Relabeled ventricle regions to CSF and cerebellum regions to WM")
    m.elm.tag1[m.elm.tag1 == 6] = 1
    m.elm.tag1[m.elm.tag1 == 7] = 3
    m.elm.tag2[m.elm.tag2 == 6] = 1
    m.elm.tag2[m.elm.tag2 == 7] = 3

    logger.info("Fixing surface labeling")
    m.fix_surface_labels()
    logger.info("Fixing thin tetrahedra")
    m.fix_thin_tetrahedra()
    logger.info("Fixing tetrahedra node ordering")
    m.fix_th_node_ordering()
    logger.info("Fixing triangle node ordering")
    m.fix_tr_node_ordering()

    mesh_io.write_msh(m, mode='binary')

if __name__ == '__main__':
    main()

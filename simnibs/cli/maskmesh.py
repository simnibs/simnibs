# -*- coding: utf-8 -*-
'''
    command line tool to convert NIfTI files to .msh to be read by gmsh
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2019  Guilherme B Saturnino

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

'''


import os
import sys
import argparse
import nibabel
import numpy as np

import simnibs.mesh_tools.mesh_io as mesh_io
from simnibs.utils.simnibs_logger import logger
from simnibs import __version__


def parseArguments(argv):

    parser = argparse.ArgumentParser(
        prog="maskmesh",
        description="Relabel elements using a nifti binary mask"
    )
    parser.add_argument("fn_mask", help="Mask file")
    parser.add_argument("fn_mesh", help="Mesh file")
    parser.add_argument('tag', help="Tag to use in relabled tissue", type=int)
    parser.add_argument('fn_out', help="Output file name")
    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def main():
    args = parseArguments(sys.argv[1:])
    if not os.path.isfile(args.fn_mask):
        raise IOError('Could not find file: {0}'.format(args.fn_mask))
    if not os.path.isfile(args.fn_mesh):
        raise IOError('Could not find file: {0}'.format(args.fn_mesh))
    image = nibabel.load(args.fn_mask)
    mesh = mesh_io.read_msh(args.fn_mesh)
    vol = image.dataobj
    if not np.issubdtype(vol.dtype, np.integer):
        logger.warning('Volume is not an integer type, masking may fail')

    affine = image.affine

    logger.info('Applying Mask')
    ed = mesh_io.ElementData.from_data_grid(
        mesh, vol, affine, 'from_volume'
    )
    # Re-label tetrahedra
    mesh.elm.tag1[(ed.value > 0) * mesh.elm.elm_type == 4] = args.tag
    mesh.elm.tag2[(ed.value > 0) * mesh.elm.elm_type == 4] = args.tag
    # Remove triangles
    mesh.elm.tag1[(ed.value > 0) * mesh.elm.elm_type == 2] = 99999
    mesh.elm.tag2[(ed.value > 0) * mesh.elm.elm_type == 2] = 99999
    mesh = mesh.remove_from_mesh(99999)

    logger.info(f'Writing {args.fn_out}')
    mesh_io.write_msh(mesh, args.fn_out)


if __name__ == '__main__':
    main()

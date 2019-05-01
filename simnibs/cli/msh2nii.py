# -*- coding: utf-8 -*-
'''
    command-line interface to convert ".msh" files to nifit files
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2018  Guilherme B Saturnino, Kristoffer H Madsen, Axel Thieslcher,
    Jesper D Nielsen, Andre Antunes

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

import simnibs.msh.transformations as transformations
from simnibs import __version__



def parse_arguments(argv):

    parser = argparse.ArgumentParser(prog="msh2nii",
                                     description="Interpolate fields in \".msh\" file to"
                                     " NifTI or create masks from a \".msh\" file")
    parser.add_argument("fn_mesh", help="Input mesh file")
    parser.add_argument("fn_reference",
                        help="Path to m2m folder or "
                        "reference nifti file. Used to get the dimensions "
                        "and affine transformation")
    parser.add_argument('fn_out', help="Output file name, "
                        "name of fields or mask will be appended")
    parser.add_argument("--create_masks", action="store_true",
                        help="Create Masks for each volume in the mesh instead of "
                        "interpolating fields")
    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def main():
    args = parse_arguments(sys.argv[1:])
    if not os.path.isfile(args.fn_mesh):
        raise IOError('Could not find file: {0}'.format(args.fn_mesh))
    transformations.interpolate_to_volume(
        args.fn_mesh, args.fn_reference, args.fn_out, create_masks=args.create_masks)

if __name__ == '__main__':
    main()

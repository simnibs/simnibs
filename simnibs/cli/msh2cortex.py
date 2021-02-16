# -*- coding: utf-8 -*-
'''
    command line tool to interpolates fields from a mesh into a middle grey/white matter surface
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
import argparse
import sys
import os

from simnibs import transformations
from simnibs import __version__


def parse_arguments(argv):

    # argument parsing exception handling
    q_help = "(Optional) Scalar quantities to calculate from vector fields"
    q_help += " (N: magnitude,"
    q_help += " n: normal component,"
    q_help += " t: tangent component,"
    q_help += " a: angle between vector and surface normal)"
    q_help += " Default: calculate all quantities"

    def allowed_fields(s):
        for ch in s:
            if ch not in "Nnta":
                raise argparse.ArgumentTypeError(q_help)
        return s



    parser = argparse.ArgumentParser(prog="msh2cortex",
                                     description="Interpolates fields from the gray"
                                     " matter volume to a cortical surface located between the gray"
                                     " and white matter surfaces."
                                     " Outputs freesurfer overlay files or gifti files")
    parser.add_argument("-i", '--in', dest='fn_in', required=True,
                        help="Input mesh with simulation results")
    parser.add_argument('-m', "--m2mpath", dest='m2mpath', required=True,
                        help="path to m2m_{subjectID} directory, created in the "
                        "segmentation")
    parser.add_argument('-o', '--out_folder', dest='out_folder', required=True,
                        help="Folder where output files will be saved")
    parser.add_argument('-d', '--depth', dest='depth', required=False, type=float,
                        help="(Optional) Depth where the field is to be interpolated."
                        " 0 means at the gray matter surface, 1 at the white matter"
                        " surface. Default: 0.5. This argument is only used if the"
                        " head mesh was generated with mri2mesh", default=0.5)
    parser.add_argument('-f', '--fsaverage_folder', dest='fsaverage_folder',
                        required=False,
                        help="(Optional) Folder where output files in fsaverage space will"
                        " be saved. If not set, the fields will not be transformed to "
                        " FsAverage",
                        default=None)
    parser.add_argument('--quantities', type=allowed_fields, required=False,
                        help=q_help, default='Nnta')
    parser.add_argument('--open-in-gmsh', action='store_true',
                        help="(Optional) If set, opens a gmsh window with the overlays after"
                        " performing the transformations")
    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def main():
    args = parse_arguments(sys.argv[1:])
    quantities = []
    for q in args.quantities:
        if q == 'N':
            quantities += ['magn']
        elif q == 'n':
            quantities += ['normal']
        elif q == 't':
            quantities += ['tangent']
        elif q == 'a':
            quantities += ['angle']
    transformations.middle_gm_interpolation(
        args.fn_in, args.m2mpath, args.out_folder,
        out_fsaverage=args.fsaverage_folder,
        depth=args.depth,
        quantities=quantities,
        fields=None,
        open_in_gmsh=args.open_in_gmsh)


if __name__ == '__main__':
    main()

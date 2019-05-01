# -*- coding: utf-8 -*-
'''
    command line tool to convert coordinates from MNI space to subject space
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
from __future__ import print_function
import argparse
from argparse import RawTextHelpFormatter
import sys
import os
import textwrap

import numpy as np
import simnibs.msh.transformations as transformations
from simnibs import __version__


def parse_arguments(argv):

    transf_types_str = '\nnonl: Non-linear transformation.\n' + \
                       '6dof: 6 degrees of freedom affine transformation.\n' +\
                       '12dof: 12 degrees of freedom affine transformation.\n'

    csv_descrip = textwrap.dedent(
            '''
            CSV file with coordinates. The csv file must have the format
            Generic:
                Generic, pos_x, pos_y, pos_z, name, ...
            Positions will not be changed after transformation.

            Fiducial, Electrode, ReferenceElectrode:
                 Type, pos_x, pos_y, pos_z, name, whatever
                 Type, pos_x, pos_y, pos_z, pos2_x, pos2_y, pos2_z, name, ...
            Positions will be projected on skin after transformation.
            Type must be Fiducial, Electrode, or ReferenceElectrode.

            CoilPos:
                Type, pos_x, pos_y, pos_z, ez_x, ez_y, ez_z, ey_x, ey_y, ey_z, dist, name, ...
            Position will be adjusted after transformation to have specified distance to skin
            ''')

    parser = argparse.ArgumentParser(prog="mni2subject_coords",
                                     description="Transform coordinates "
                                     "from MNI space to subject space using the "
                                     "deformation fields calculated during the "
                                     "segmentation.",
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-m', "--m2mpath", dest='m2mpath', required=True,
                        help="path to m2m_{subjectID} directory, created in the "
                        "segmentation")
    parser.add_argument("-c", '--coords', dest='coords', nargs=3, action='store',
                        help="Input coordinates. 3 values separated by spaces",
                        default=None)
    parser.add_argument('-s', '--csv', dest='csv', required=False,
                        help=csv_descrip, default=None)
    parser.add_argument('-o', '--out', dest='out',
                        help="CSV file where the transformed coordinates will be saved. "
                        "Also prints a '.geo' file that can be visualized in gmsh "
                        "If not set, the values will be printed to the screen", default=None)
    parser.add_argument('-t', "--transformation_type", dest='t',
                        help='(optional) Type of transformation to use.'
                        '{0} Default: nonl'.format(transf_types_str), default='nonl',
                        choices=['nonl', '6dof', '12dof'])
    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def main():
    args = parse_arguments(sys.argv[1:])
    m2m_dir = os.path.abspath(os.path.realpath(os.path.expanduser(args.m2mpath)))
    if not os.path.isdir(m2m_dir):
        raise IOError('Could not find directory: {0}'.format(args.m2mpath))
    if args.out is not None:
        fn_out = os.path.abspath(os.path.realpath(os.path.expanduser(args.out)))
        fn_geo = os.path.splitext(fn_out)[0] + '.geo'
    else:
        fn_out = None
        fn_geo = None

    if args.coords is not None:
        coords = [float(d) for d in args.coords]

    elif args.csv is not None:
        coords = os.path.abspath(os.path.realpath(os.path.expanduser(args.csv)))
        if not os.path.isfile(coords):
            raise IOError('Could not find CSV file: {0}'.format(args.csv))
    else:
        raise argparse.ArgumentTypeError(
            'Plase use either -c or -s')

    if args.coords is not None and args.csv is not None:
        raise argparse.ArgumentError(
            'Please use only -c or -s')

    coords = transformations.warp_coordinates(
        coords, m2m_dir,
        out_name=fn_out,
        transformation_direction='mni2subject',
        transformation_type=args.t,
        out_geo=fn_geo)[1]
    if fn_out is None:
        np.set_printoptions(precision=2,
                            formatter={'float': lambda x: '{0:.2f}'.format(x)})
        print('Transformed coodinates:\n{0}'.format(coords.squeeze()))


if __name__ == '__main__':
    main()

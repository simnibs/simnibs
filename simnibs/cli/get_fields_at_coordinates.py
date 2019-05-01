# -*- coding: utf-8 -*-\
'''
    Interpolates field values at given coordinates
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
import sys
import os
import textwrap
from argparse import RawTextHelpFormatter
import simnibs.msh.mesh_io as mesh_io
from simnibs.utils.simnibs_logger import logger
import numpy as np
from simnibs import __version__


def parse_arguments(argv):

    csv_descrip = textwrap.dedent(
            '''
            CSV file with coordinates. The csv file must have the format
            pos_x, pos_y, pos_z
            ''')

    out_descrip = textwrap.dedent(
            '''
            (Optional) what value to fill in for points outside the volume:
            nearest:    Use the nearest value inside the volume
            nan:        fill with NaNs
            number:     Use the given numerical value (example: 0.0)
            (Default = nan)
            ''')

    type_descrip = textwrap.dedent(
            '''
            (Optional) Type of interpolation to use
            linear:     For element data (eg: E, e, J, j), first interpolate nodal values
                        using SPR (Zienkiewicz and Zhu, 1992) and proceed to use linear
                        interpolation. v is already nodal
            assign:     use the raw, discontinuous values
            (Defaul: linear)
            ''')

    def check_out_fill(value):
        if value == 'nearest':
            return 'nearest'
        if value == 'nan':
            return np.nan
        try:
            return float(value)
        except:
            raise argparse.ArgumentTypeError(
                "invalid argument: {0} please use get_fields_at_coordinates --help for more "
                "information".format(value))

    parser = argparse.ArgumentParser(prog="get_fields_at_coordinates",
                                     description="Interpolate field values "
                                     "at given points. Outputs one csv file per field "
                                     "in the mesh with field values",
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-s', '--csv', dest='csv',
                        help=csv_descrip, required=True)
    parser.add_argument("-m", '--mesh', dest='mesh',
                        help='Mesh file with the fields to be interpolated',
                        required=True)
    parser.add_argument('-l', '--labels', action='append', dest='l',
                        help='(Optional) Volume labels where the interpolation is to be '
                        'performed. Should be separated by a comma.',
                        default=None)
    parser.add_argument('--out_fill',
                        help=out_descrip,
                        type=check_out_fill, default='nan')
    parser.add_argument('--method',
                        help=type_descrip,
                        default='linear', choices=['linear', 'assign'])

    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def main():
    args = parse_arguments(sys.argv[1:])
    msh = mesh_io.read_msh(os.path.abspath(os.path.realpath(os.path.expanduser(args.mesh))))
    fn_csv = os.path.expanduser(args.csv)
    if args.l is not None:
        labels = [int(n) for n in args.l[0].split(',')]
        msh = msh.crop_mesh(tags=labels)

    points = np.atleast_2d(np.loadtxt(fn_csv, delimiter=',', dtype=float))
    if points.shape[1] != 3:
        raise IOError('CSV file should have 3 columns')

    csv_basename, _ = os.path.splitext(fn_csv)
    for ed in msh.elmdata:
        c = ed.field_name in ['normJ', 'J']
        fn_out = '{0}_{1}.csv'.format(csv_basename, ed.field_name)
        logger.info('Interpolating field: {0} and writing to file: {1}'.format(
            ed.field_name, fn_out))
        f = ed.interpolate_scattered(
            points, out_fill=args.out_fill,
            method=args.method, continuous=c, squeeze=False)
        if f.ndim == 1:
            f = f[:, None]
        np.savetxt(fn_out, f, delimiter=',')

    for nd in msh.nodedata:
        fn_out = '{0}_{1}.csv'.format(csv_basename, nd.field_name)
        logger.info('Interpolating field: {0} and writing to file: {1}'.format(
            nd.field_name, fn_out))
        f = nd.interpolate_scattered(
            points, out_fill=args.out_fill, squeeze=False)
        if f.ndim == 1:
            f = f[:, None]
        np.savetxt(fn_out, f, delimiter=',')


if __name__ == '__main__':
    main()

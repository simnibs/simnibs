# -*- coding: utf-8 -*-
'''
    Command line tool to calculate magnetic fields 

    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2019 Hassan Yazdanian, Guilherme B Saturnino

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

import nibabel as nib
from simnibs import msh
from simnibs.simulation.biot_savart import calc_B
from simnibs import __version__

simnibs_version = str(__version__)


def parse_arguments(argv):
    parser = argparse.ArgumentParser(prog="calc_B", description="Calculate magnetic "
                                     "fields from current desity field")
    parser.add_argument("-i", '--in', dest='fn_in', required=True,
                        help="Input msh file. Must have a field called 'J'")
    parser.add_argument("-r", '--reference_nifti', dest='ref', required=True,
                        help="Reference volume for the B field calculation")
    parser.add_argument('-o', '-out', dest='out', required=True, help="Output nifti file name")
    parser.add_argument('--res', dest='res', default=2., type=float,
                        help="Resolution for B field calculations. Default: 2mm")
    parser.add_argument('--mask', action='store_true',
                        help='Masks the output volume using the mesh')
    parser.add_argument('--version', action='version', version=simnibs_version)
    return parser.parse_args(argv)


def main():
    args = parse_arguments(sys.argv[1:])
    fn_in = os.path.abspath(os.path.realpath(os.path.expanduser(args.fn_in)))
    if not os.path.isfile(fn_in):
        raise IOError('Could not find file: {0}'.format(args.fn_in))
    fn_out = os.path.abspath(os.path.realpath(os.path.expanduser(args.out)))
    if os.path.splitext(fn_out)[1] == '':
        fn_out += '.nii.gz'
    reference = os.path.abspath(os.path.realpath(os.path.expanduser(args.ref)))

    affine = nib.load(reference).affine
    nvox = nib.load(reference).shape
    if len(nvox) == 2:
        nvox = nvox + (1,)

    m = msh.read_msh(fn_in)
    if 'J' not in m.field.keys():
        raise IOError('Input mesh file must have a current density (J) field')

    B = calc_B(m.field['J'], nvox, affine, calc_res=args.res, mask=args.mask)

    img = nib.Nifti1Pair(B, affine)
    img.header.set_xyzt_units('mm')
    img.set_qform(affine)
    img.update_header()
    nib.save(img, fn_out)


if __name__ == '__main__':
    main()

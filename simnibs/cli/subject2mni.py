# -*- coding: utf-8 -*-
'''
    converts fields from subject space to mni space
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

    transf_types_str = '\n nonl: Non-linear transformation.\n' + \
                       '6dof: 6 degrees of freedom affine transformation.\n' +\
                       '12dof: 12 degrees of freedom affine transformation.\n'

    parser = argparse.ArgumentParser(prog="subject2mni", description="Transform fields"
                                     "from subject space to MNI space using the "
                                     "deformation fields calculated during the "
                                     "segmentation. The input can be a mesh or a "
                                     "nifiti file. The output is always a nifti file ")
    parser.add_argument("-i", '--in', dest='fn_in', required=True,
                        help="Input nifti or mesh file to be transformed to MNI space")
    parser.add_argument("-m", "--m2mpath", dest='m2mpath', required=True,
                        help="path to m2m_{subjectID} directory, created in the "
                        "segmentation")
    parser.add_argument('-o', "--out", dest='out', required=True,
                        help="Output nifti file name. If the input is a mesh, "
                        "the name of the fields will be appended")
    parser.add_argument("--transformation_type", '-t', dest='t',
                        help='(optional) Type of transformation to use.'
                        '{0} Default: nonl'.format(transf_types_str), default='nonl',
                        choices=['nonl', '6dof', '12dof'])
    parser.add_argument('--reference_nifti', '-r', dest='r',
                        help='(optional) Reference nifti file in  MNI space. '
                        'Default: $SIMNIBSDIR/resources/templates/MNI152_T1_1mm.nii.gz',
                        default=None)
    parser.add_argument('--mask', dest='mask',
                        help='(optional) File name of mask to be applied to volume before transformation to MNI '
                        'space. Only valid if the input is a nifti file',
                        default=None)
    parser.add_argument('--labels', '-l', action='append', dest='l',
                        help='(optional) Volume labels be transformed to MNI space. '
                        'Should be separeted by a comma. Only valid if '
                        'the input is a mesh file', default=None)
    parser.add_argument('--interpolation_order', type=int, dest='int_order',
                        help='(optional) Interpolation order of spline interpolation to be used. '
                        'should be between 0 and 5. Default: 1', default=1)
    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def main():
    args = parse_arguments(sys.argv[1:])
    fn_in = os.path.abspath(os.path.realpath(os.path.expanduser(args.fn_in)))
    if not os.path.isfile(fn_in):
        raise IOError('Could not find file: {0}'.format(args.fn_in))
    m2m_dir = os.path.abspath(os.path.realpath(os.path.expanduser(args.m2mpath)))
    if not os.path.isdir(m2m_dir):
        raise IOError('Could not find directory: {0}'.format(args.m2mpath))
    fn_out = os.path.abspath(os.path.realpath(os.path.expanduser(args.out)))
    if os.path.splitext(fn_out)[1] == '':
        fn_out += '.nii.gz'
    if args.l is not None:
        labels = [int(n) for n in args.l[0].split(',')]
    else:
        labels = None
    if args.r is not None:
        reference = os.path.abspath(os.path.realpath(os.path.expanduser(args.r)))
    else:
        reference = None
    if args.mask is not None:
        mask = os.path.abspath(os.path.realpath(os.path.expanduser(args.mask)))
    else:
        mask = None
    transformations.warp_volume(
        fn_in, m2m_dir, fn_out,
        transformation_direction='subject2mni',
        transformation_type=args.t,
        reference=reference,
        mask=mask,
        labels=labels,
        order=args.int_order)


if __name__ == '__main__':
    main()

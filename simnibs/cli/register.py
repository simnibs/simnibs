# -*- coding: utf-8 -*-
'''
    Command line tool registering scans either rigidly or affinely.
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2022  Oula Puonti, Kristoffer H Madsen, Axel Thieslcher,
    Jesper D Nielsen

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
from samseg import gems
from simnibs import __version__


def parseArguments(argv):

    parser = argparse.ArgumentParser(prog="register",
                                     description="Register the two input scans either rigidly or affinely."
                                     "The fixed scan is the target space and the moving scan is the one "
                                     "being registered to that space. Defaults to rigid."
                                     "Note that the rigid registration uses mutual information "
                                     "as a registration metric and is thus suitable for cross-contrast"
                                     "registration, whereas the affine uses cross-correlation and thus only"
                                     "support intra-contrast.")
    parser.add_argument("-f", "--fixed", dest="fixed_scan", required=True, help="Fixed scan.")
    parser.add_argument("-m", "--moving", dest="moving_scan", required=True, help="Moving scan.")
    parser.add_argument("-dof", "--degrees_of_freedom", dest="dof", help="Degrees of freedom: 6 (rigid (default)) or 12 (affine).")
    parser.add_argument("-o", "--output", dest="output", required=True, help="Output name, including path.")
    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def main():
    args = parseArguments(sys.argv[1:])
    if not os.path.isfile(args.fixed_scan):
        raise IOError('Could not find file: {0}'.format(args.fixed_scan))
    if not os.path.isfile(args.moving_scan):
        raise IOError('Could not find file: {0}'.format(args.moving_scan))
    if args.output==None:
        raise IOError('Error specify an output name.')
    if args.dof == None:
        args.dof = '6'
    else:
        if not (int(args.dof) == 6 or int(args.dof) == 12):
            args.dof = '6'

    # Figure out output names and path for writing out registration matrices
    path_name = os.path.split(args.output)
    filename = path_name[1]
    filename = filename.split(".")
    filename = filename[0]
    RAS2LPS = np.diag([-1, -1, 1, 1])
    if int(args.dof) == 6:
        reg = gems.KvlRigidRegistration()
        reg.read_images(args.fixed_scan, args.moving_scan)
        reg.initialize_transform()
        reg.register()
        trans_mat = RAS2LPS@reg.get_transformation_matrix()@RAS2LPS
        reg.write_out_result(args.output)
        mat_path = os.path.join(path_name[0], filename + '_dof6.dat')
        np.savetxt(mat_path, trans_mat)
    else:
        reg = gems.KvlAffineRegistration(-100,
                                         300,
                                         0,
                                         [2.0, 1.0, 0],
                                         0,
                                         [4.0, 2.0, 0],
                                         True,
                                         1.0,
                                         "b")
        reg.read_images(args.fixed_scan, args.moving_scan)
        reg.initialize_transform()
        reg.register()
        trans_mat = RAS2LPS@reg.get_transformation_matrix()@RAS2LPS
        reg.write_out_result(args.output)
        mat_path = os.path.join(path_name[0], filename + '_dof12.dat')
        np.savetxt(mat_path, trans_mat)

if __name__ == '__main__':
    main()

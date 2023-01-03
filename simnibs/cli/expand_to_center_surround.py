# -*- coding: utf-8 -*-
'''
    expand a TDCSLIST with one electrode to center-surround montage

    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2018  Guilherme B Saturnino
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
import scipy
from argparse import RawTextHelpFormatter
from simnibs import __version__, sim_struct


def parse_arguments(argv):

    parser = argparse.ArgumentParser(prog="expand_to_center_surround",
                                     description="expand a TDCSLIST with one electrode to center-surround montage",
                                     formatter_class=RawTextHelpFormatter)
    parser.add_argument('-S', dest='S', help='.mat with tdcslist', required=True)
    parser.add_argument('-p', dest='subpath', help='m2m-folder', required=True)
    parser.add_argument('-F', dest='fn_out', help='output name', required=True)

    parser.add_argument('--radius_surround', dest='radius_surround', default=[50.0], nargs='+')
    parser.add_argument('--pos_dir_1stsurround', dest='pos_dir_1stsurround', default=None, nargs='+')
    parser.add_argument('--N', dest='N', default=4, type=int)
    parser.add_argument('--multichannel', dest='multichannel', action='store_true')
    parser.add_argument('--phis_surround', dest='phis_surround', default=None, nargs='+')
    
    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def main():
    args = parse_arguments(sys.argv[1:])
    mat = scipy.io.loadmat(os.path.abspath(os.path.realpath(os.path.expanduser(args.S))))    
    S = sim_struct.TDCSLIST()
    S.read_mat_struct(mat)
    
    subpath = os.path.abspath(os.path.realpath(os.path.expanduser(args.subpath)))
    if not os.path.isdir(subpath):
        raise IOError('Could not find directory: {0}'.format(args.subpath))        
        
    fn_out = os.path.abspath(os.path.realpath(os.path.expanduser(args.fn_out)))
        
    radius_surround = None
    pos_dir_1stsurround = None
    phis_surround = None
    if args.radius_surround:
        radius_surround = list(map(float, args.radius_surround))
    if args.pos_dir_1stsurround:
        try:
            pos_dir_1stsurround = list(map(float, args.pos_dir_1stsurround))
        except:
            pos_dir_1stsurround = args.pos_dir_1stsurround[0]
    if args.phis_surround:
        phis_surround = list(map(float, args.phis_surround))

    S.expand_to_center_surround(subpath,
                                radius_surround = radius_surround,
                                pos_dir_1stsurround = pos_dir_1stsurround,
                                N = args.N,
                                multichannel = args.multichannel,
                                phis_surround = phis_surround)                            
    sim_struct.save_matlab_sim_struct(S, fn_out)        

        
if __name__ == '__main__':
    main()    
    
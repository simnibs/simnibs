# -*- coding: utf-8 -*-\
'''
    Updates a head model from 2.1.0 to 2.1.1
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
from simnibs import SIMNIBSDIR
from simnibs.msh.transformations import warp_coordinates
from simnibs import __version__



def parse_arguments(argv):
    parser = argparse.ArgumentParser(prog="get_eeg_positions",
                                     description=
                                     "Creates the eeg_positions "
                                     "folder. This folder is nescessary "
                                     "for using the eeg cap positions features "
                                     "in SimNIBS 2.1.1")
    parser.add_argument('-m', "--m2mpath", dest='m2mpath', required=True,
                        help="path to m2m_{subjectID} directory, created in the "
                        "segmentation")

    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def update_hm(m2m_dir):
    eeg_positions = os.path.join(m2m_dir, "eeg_positions")
    print("Transforming EEG positions")
    if not os.path.exists(eeg_positions):
        os.mkdir(eeg_positions)
    cap_file = os.path.abspath(os.path.realpath(
                    os.path.join(
                        SIMNIBSDIR,
                        'resources', 'ElectrodeCaps_MNI',
                        'EEG10-10_UI_Jurak_2007.csv')))
    cap_out = os.path.join(eeg_positions, 'EEG10-10_UI_Jurak_2007.csv')
    geo_out = os.path.join(eeg_positions, 'EEG10-10_UI_Jurak_2007.geo')
    warp_coordinates(
        cap_file, m2m_dir,
        transformation_direction='mni2subject',
        out_name=cap_out,
        out_geo=geo_out)


def main():
    args = parse_arguments(sys.argv[1:])
    m2m_dir = os.path.abspath(os.path.realpath(os.path.expanduser(args.m2mpath)))
    if not os.path.isdir(m2m_dir):
        raise IOError('Could not find directory: {0}'.format(args.m2mpath))
    update_hm(m2m_dir)


if __name__ == '__main__':
    main()

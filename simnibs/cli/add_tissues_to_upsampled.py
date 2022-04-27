# -*- coding: utf-8 -*-
'''
    Command line tool add extra tissues to the upsampled tissue labeling.
    Upsamples manual/automated segs that are in the space of the input to
    the same resolution as the upsampled tissue mask and adds them.
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
from simnibs.utils.transformations import volumetric_affine
from simnibs import __version__


def parseArguments(argv):

    parser = argparse.ArgumentParser(prog="add_tissues_to_upsampled",
                                     description="Resamples the input label data to the space "
                                     "of the target scan and adds the labels. Meant to be used "
                                     "to add tissue masks that are in the same space as the "
                                     "charm input (T1.nii.gz) to the upsampled tissue mask "
                                     "(tissue_labels_upsampled.nii.gz) for meshing. "
                                     "Ignores labels that are smaller than or equal to zero.")
    parser.add_argument("label_input", help="Input label file to be resampled")
    parser.add_argument("target_scan", help="Target label file where the iput labels will be added.")
    parser.add_argument('fn_out', help="Output file name")
    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def main():
    args = parseArguments(sys.argv[1:])
    if not os.path.isfile(args.label_input):
        raise IOError('Could not find file: {0}'.format(args.label_input))
    if not os.path.isfile(args.target_scan):
        raise IOError('Could not find file: {0}'.format(args.target_scan))
    input_scan = nibabel.load(args.label_input)
    target_scan = nibabel.load(args.target_scan)
    input_vol = input_scan.dataobj
    input_affine = input_scan.affine
    target_vol = target_scan.dataobj
    target_affine = target_scan.affine

    # Upsample the input to the target with nearest neighbor
    upsampled = volumetric_affine((input_vol, input_affine),
                                  np.eye(4), target_affine,
                                  target_vol.shape, intorder=0)

    # Check the largest label in the target and use as base label
    base_label = np.unique(target_vol).max()
    input_labels = np.unique(input_vol)
    input_labels = input_labels[input_labels > 0]
    target_array = target_scan.get_fdata()
    for l in input_labels.tolist():
        target_array[upsampled == l] = base_label + l

    if args.fn_out is None:
        split_path = os.path.split(args.target_scan)
        output_path =  split_path[0]
        filename = split_path[1]
        parts = filename.split(".")
        output_name = os.path.join(output_path, part[0]+'_added.nii.gz')
    else:
        output_name = args.fn_out

    output_image = nibabel.Nifti1Image(target_array, target_affine)
    nibabel.save(output_image, output_name)

if __name__ == '__main__':
    main()

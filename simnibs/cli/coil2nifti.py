# -*- coding: utf-8 -*-
"""
    command line tool to convert coil dipole definition ccd files to nifti1
    format. This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2021  Kristoffer H. Madsen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

"""

from itertools import chain
import os
import glob
import re
import numpy as np
import fmm3dpy
import nibabel as nib
import time

from simnibs.simulation.tms_coil.tms_coil import TmsCoil


def main():
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Convert coil files (ccd, tcd) to Nifti1 format"
    )
    parser.add_argument(
        "-i",
        "--infile",
        dest="infile",
        default=None,
        required=True,
        help="CCD file to convert",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        default=None,
        help="output filename, will default to replacing extension with .nii.gz",
    )
    parser.add_argument(
        "-f", "--force", dest="force", action="store_true", help="Force rewrite"
    )
    parser.add_argument(
        "-b",
        "--bfield",
        dest="Bfield",
        action="store_true",
        help="Write B field instead of A field",
    )

    options = parser.parse_args(sys.argv[1:])
    if os.path.isdir(options.infile):
        print(f"recursively processing CCD files in {options.infile}")
        coil_files = chain(
            glob.iglob(os.path.join(options.infile, "**", "*.ccd"), recursive=True),
            glob.iglob(os.path.join(options.infile, "**", "*.tcd"), recursive=True),
        )
        options.outfile = None
    elif os.path.isfile(options.infile):
        coil_files = (options.infile,)
    else:
        print(f"Cannot locate input file: {options.infile}")

    for coil_file in coil_files:
        if options.outfile is None:
            outfile = os.path.splitext(coil_file)[0] + ".nii.gz"
        else:
            outfile = options.outfile
        if len(glob.glob(os.path.splitext(outfile)[0] + "*")) == 0 or options.force:
            t0 = time.perf_counter()
            print(f"expanding coil file {coil_file}")

            coil = TmsCoil.from_file(coil_file)
            if coil.limits is None:
                coil.limits = np.array([[-300, 300], [-200, 200], [0, 300]])
            coil.write_nifti(outfile, b_field=options.Bfield)
            print(f"Time spend: {time.perf_counter()-t0:.0f}s")
        else:
            print(f"Nifti1 version of {coil_file} already exists")


if __name__ == "__main__":
    main()

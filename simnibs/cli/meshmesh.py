# -*- coding: utf-8 -*-
'''
    command line tool to create tetrahedral head meshes from label images
    such as 'tissue_labeling_upsampled.nii.gz' stored in the 
    'm2m_{subID}/label_prep' subfolder
    
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2020  Oula Puonti, Guilherme B Saturnino, Jesper D Nielsen, 
    Fang Cao, Axel Thielscher
    
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
import os
import re
import sys
import shutil
import textwrap
import nibabel as nib
import numpy as np
import logging

from simnibs import __version__
from simnibs import SIMNIBSDIR
from .. import utils
from simnibs.utils import cond_utils, file_finder
from simnibs.mesh_tools.mesh_io import write_msh
from simnibs.mesh_tools.meshing import create_mesh
from simnibs.utils.transformations import resample_vol, crop_vol
from simnibs.utils.settings_reader import read_ini
from simnibs.segmentation.charm_main import _read_settings_and_copy
from simnibs.utils.simnibs_logger import logger


def parseArguments(argv):
    
    usage_text = textwrap.dedent('''
                                 
CREATE HEAD MESH FROM LABEL IMAGE:
    meshmesh label_image mesh_name
            
OPTIONAL:
internally upsample label image to specific resolution (in [mm]) before meshing:
    meshmesh label_image mesh_name --voxsize_meshing 0.5
    
use custom .ini-file to control tetrahedra sizes:
    meshmesh label_image mesh_name --usesettings my_settings.ini
    ''')
    
    parser = argparse.ArgumentParser(
        prog="meshmesh",usage=usage_text)
    parser.add_argument('label_image', nargs='?', help="filename of the nifti label image")
    parser.add_argument('mesh_name', nargs='?', help="filename for the new mesh")
    parser.add_argument('--usesettings', default=None, nargs=1, metavar="settings.ini",
                        help="""ini-file with custom settings""")
    parser.add_argument('--voxsize_meshing', dest='voxsize_meshing', 
                        metavar="voxel_size_in_mm", default=None, nargs=1,
                        help="""voxel size to which the label image is internally
                        upsampled for meshing. Rule of thumb: use twice the minimum resolution 
                        unless image is already upsampled, 
                        such as tissue_labeling_upsampled.nii.gz""")
    parser.add_argument('--sizing_field', dest='fn_sizing_field', default=None, 
                        nargs=1, metavar="size_field.nii.gz", 
                        help="""optional sizing field to locally control element 
                        sizes""")
    parser.add_argument('--nthreads', type=int,  dest='nthreads', default=1,
                        help="""Number of threads to be used for the meshing.""")
    parser.add_argument('-v', '--version', action='version', version=__version__)
    parser.add_argument('--debug', help='Enable debug mode', action='store_true', required=False,
                        default=False)
    
    args=parser.parse_args(argv)
    if args.label_image is None:
        parser.print_help()
        exit()
    return args


def main():
    args = parseArguments(sys.argv[1:])
    order_upsampling = 1;  # 0 [nearest neighbor] or 1 [linear]

    # setup logger
    sub_files = file_finder.SubjectFiles(subpath=os.path.split(args.mesh_name)[0])
    _setup_logger(sub_files.charm_log)

    # load settings
    if args.usesettings is None:
        args.usesettings = os.path.join(SIMNIBSDIR, 'charm.ini')  # use standard settings
    if type(args.usesettings) == list:
        args.usesettings = args.usesettings[0]        
    # settings = read_ini(args.usesettings)['mesh']
    settings = _read_settings_and_copy(args.usesettings, sub_files.settings)
    if not settings['mesh']['skin_tag']:
        settings['mesh']['skin_tag'] = None
    if not settings['mesh']['hierarchy']:
        settings['mesh']['hierarchy'] = None

    debug_path = None
    if args.debug:
        debug_path = sub_files.subpath

    # load label image
    label_nifti = nib.load(args.label_image)
    label_image = np.squeeze(label_nifti.get_fdata().astype(np.uint16)) # Cast to uint16, otherwise meshing complains
    label_affine = label_nifti.get_qform()
    
    # load sizing field (optional) 
    sf_image = None
    if args.fn_sizing_field is not None:
        if type(args.fn_sizing_field) == list:
                args.fn_sizing_field = args.fn_sizing_field[0]
        sf_nifti = nib.load(args.fn_sizing_field)
        sf_image = np.squeeze(sf_nifti.get_fdata())
        assert sf_image.shape == label_image.shape

    # upsample (optional)
    if args.voxsize_meshing is not None:
        logger.info('upsampling label image ...')
        if type(args.voxsize_meshing) == list:
                args.voxsize_meshing = args.voxsize_meshing[0]

        args.voxsize_meshing = float(args.voxsize_meshing)
        if sf_image is not None:
            sf_image, _, _ = resample_vol(sf_image, label_affine,
                                          args.voxsize_meshing, order=order_upsampling)
        if order_upsampling == 0:
            label_image, label_affine, _ = resample_vol(label_image, label_affine,
                                                        args.voxsize_meshing, order=0)
        else:
            labels = np.unique(label_image)
            best_tag_p, new_affine, _ = resample_vol((label_image==labels[0]).astype(np.float32), 
                                                     label_affine, args.voxsize_meshing, order=1)
            best_tag = np.zeros_like(best_tag_p,dtype=np.uint16)
            best_tag[:] = labels[0]
            for i in labels[1:]:
                tag_p = resample_vol((label_image==i).astype(np.float32), 
                                     label_affine, args.voxsize_meshing, order=1)[0]
                idx = tag_p > best_tag_p
                best_tag[idx] = i
                best_tag_p[idx] = tag_p[idx]
            label_image = best_tag
            label_affine = new_affine
        
    # reduce memory consumption a bit
    if sf_image is not None:
        sf_image, _, _ = crop_vol(sf_image, label_affine,
                                  label_image>0, thickness_boundary=5)

    label_image, label_affine, _ = crop_vol(label_image, label_affine,
                                            label_image>0, thickness_boundary=5)

    # meshing
    new_mesh = create_mesh(label_image, label_affine,
                           elem_sizes=settings['mesh']['elem_sizes'],
                           smooth_size_field=settings['mesh']['smooth_size_field'],
                           skin_facet_size=settings['mesh']['skin_facet_size'],
                           facet_distances=settings['mesh']['facet_distances'],
                           optimize=settings['mesh']['optimize'],
                           remove_spikes=settings['mesh']['remove_spikes'],
                           skin_tag=settings['mesh']['skin_tag'],
                           hierarchy=settings['mesh']['hierarchy'],
                           smooth_steps=settings['mesh']['smooth_steps'],
                           sizing_field=sf_image,
                           num_threads=args.nthreads,
                           debug=args.debug,
                           debug_path=debug_path)

    # write out mesh, .opt-file with some visualization settings and .ini-file with settings
    if not args.mesh_name.lower().endswith('.msh'):
        args.mesh_name+='.msh'
    write_msh(new_mesh, args.mesh_name)
    v = new_mesh.view(cond_list=cond_utils.standard_cond())
    v.write_opt(args.mesh_name)
    try:
        shutil.copyfile(args.usesettings, args.mesh_name+'.ini')
    except shutil.SameFileError:
        pass


def _setup_logger(logfile):
    """Add FileHandler etc."""
    with open(logfile, "a") as f:
        f.write("<HTML><HEAD><TITLE>charm report</TITLE></HEAD><BODY><pre>")
        f.close()
    fh = logging.FileHandler(logfile, mode="a")
    formatter = logging.Formatter("%(levelname)s: %(message)s")
    fh.setFormatter(formatter)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    utils.simnibs_logger.register_excepthook(logger)


def _stop_logger(logfile):
    """Close down logging"""
    while logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])
    utils.simnibs_logger.unregister_excepthook()
    logging.shutdown()
    with open(logfile, "r") as f:
        logtext = f.read()

    # Explicitly remove this really annoying stuff from the log
    removetext = (
        re.escape("-\|/"),
        re.escape("Selecting intersections ... ")
        + "\d{1,2}"
        + re.escape(" %Selecting intersections ... ")
        + "\d{1,2}"
        + re.escape(" %"),
    )
    with open(logfile, "w") as f:
        for text in removetext:
            logtext = re.sub(text, "", logtext)
        f.write(logtext)
        f.write("</pre></BODY></HTML>")
        f.close()


if __name__ == '__main__':
    main()
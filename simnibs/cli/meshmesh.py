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
import sys
import shutil
import textwrap
import nibabel as nib
import numpy as np

from simnibs import __version__
from simnibs import SIMNIBSDIR
from simnibs.simulation import cond
from simnibs.mesh_tools.mesh_io import write_msh
from simnibs.mesh_tools.meshing import create_mesh
from simnibs.utils.transformations import resample_vol, crop_vol
from simnibs.utils.settings_reader import read_ini


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
    parser.add_argument('-v','--version', action='version', version=__version__)
    
    args=parser.parse_args(argv)
    if args.label_image is None:
        parser.print_help()
        exit()
    return args


def main():
    args = parseArguments(sys.argv[1:])
    
    # load settings
    if args.usesettings is None:
        args.usesettings = os.path.join(SIMNIBSDIR, 'charm.ini') # use standard settings
    if type(args.usesettings) == list:
        args.usesettings = args.usesettings[0]        
    settings = read_ini(args.usesettings)['mesh']
    if not settings['skin_tag']:
        settings['skin_tag'] = None
    if not settings['hierarchy']:
        settings['hierarchy'] = None
     
    # load label image
    label_nifti = nib.load(args.label_image)
    label_image = label_nifti.get_fdata().astype(np.uint16) # Cast to uint16, otherwise meshing complains
    label_affine = label_nifti.get_qform()
    
    # load sizing field (optional) 
    sf_image = None
    if args.fn_sizing_field is not None:
        if type(args.fn_sizing_field) == list:
                args.fn_sizing_field = args.fn_sizing_field[0]
        sf_nifti = nib.load(args.fn_sizing_field)
        sf_image = sf_nifti.get_fdata()
        assert sf_image.shape == label_image.shape
        
    # upsample (optional)
    if args.voxsize_meshing is not None:
        if type(args.voxsize_meshing) == list:
                args.voxsize_meshing = args.voxsize_meshing[0]
        args.voxsize_meshing = float(args.voxsize_meshing)
        if sf_image is not None:
            sf_image, _, _ = resample_vol(sf_image, label_affine,
                                          args.voxsize_meshing, order=0)
        label_image, label_affine, _ = resample_vol(label_image, label_affine,
                                                    args.voxsize_meshing, order=0)
    
    # reduce memory consumption a bit
    if sf_image is not None:
            sf_image, _, _ = crop_vol(sf_image, label_affine,
                                      label_image>0, thickness_boundary=5) 
    label_image, label_affine, _ = crop_vol(label_image, label_affine, 
                                            label_image>0, thickness_boundary=5) 
    
    # meshing
    new_mesh = create_mesh(label_image, label_affine,
                            elem_sizes=settings['elem_sizes'],
                            smooth_size_field=settings['smooth_size_field'],
                            skin_facet_size=settings['skin_facet_size'], 
                            facet_distances=settings['facet_distances'],
                            optimize=settings['optimize'], 
                            remove_spikes=settings['remove_spikes'], 
                            skin_tag=settings['skin_tag'],
                            remove_twins=settings['remove_twins'], 
                            hierarchy= settings['hierarchy'],
                            smooth_steps=settings['smooth_steps'],
                            sizing_field=sf_image)
    
    # write out mesh, .opt-file with some visualization settings and .ini-file with settings
    if not args.mesh_name.lower().endswith('.msh'):
        args.mesh_name+='.msh'
    write_msh(new_mesh, args.mesh_name)
    v = new_mesh.view(cond_list = cond.standard_cond())
    v.write_opt(args.mesh_name)
    try:
        shutil.copyfile(args.usesettings, args.mesh_name+'.ini')
    except shutil.SameFileError:
        pass


if __name__ == '__main__':
    main()

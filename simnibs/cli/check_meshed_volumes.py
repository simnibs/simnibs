# -*- coding: utf-8 -*-\

'''
    Checks if a head model has been correctly meshed
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

import os
import sys
import argparse
import nibabel
from simnibs.msh.hmutils import check_volumes_meshed
from simnibs import file_finder
import simnibs
from simnibs import __version__




def parse_arguments(argv):
    parser = argparse.ArgumentParser(prog="check_meshed_volumes",
                                     description=
                                     "Checks if the volumes in the head mesh "
                                     "have been meshed correctly ")
    parser.add_argument('-m', "--m2mpath", dest='m2mpath', required=True,
                        help="path to m2m_{subjectID} directory, "
                        "created in the segmentation")
    parser.add_argument("--remove", action='store_true',
                        help="whether to remove the mesh")

    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def check(m2m_folder, remove):
    files = file_finder.SubjectFiles(subpath=m2m_folder)
    # Create masks for mri2mesh
    if files.seg_type == 'mri2mesh':
        mesh = simnibs.msh.read_msh(files.fnamehead).crop_mesh(elm_type=4)
        ref = nibabel.load(files.reference_volume)
        affine = ref.affine
        n_voxels = ref.header['dim'][1:4]
        qform = ref.header.get_qform()
        ed = simnibs.msh.ElementData(mesh.elm.tag1, mesh=mesh)
        ed.to_nifti(n_voxels, affine,
                    fn=files.final_contr,
                    qform=qform,
                    method='assign')

        mask_prep_dir = os.path.join(files.subpath, 'mask_prep')
        masks = dict.fromkeys(['ventricles', 'cerebellum', 'csf', 'skull', 'skin'])
        for k in masks.keys():
            masks[k] = nibabel.load(
                os.path.join(
                    mask_prep_dir, f'{k}_raw.nii.gz')).get_data().astype(bool)

        masks['gm'] = nibabel.load(
                os.path.join(files.subpath, 'gm.nii.gz')).get_data().astype(bool)

        masks['wm'] = nibabel.load(
                os.path.join(files.subpath, 'wm.nii.gz')).get_data().astype(bool)

        mask = 8 * masks['ventricles']
        wm = (masks['wm'] + masks['cerebellum']) * (mask == 0)
        mask += 1 * wm
        gm = masks['gm'] * (mask == 0)
        mask += 2 * gm
        csf = masks['csf'] * (mask == 0)
        mask += 3 * csf
        skull = masks['skull'] * (mask == 0)
        mask += 4 * skull
        skin = masks['skin'] * (mask == 0)
        mask += 5 * skin

        ref = nibabel.load(files.reference_volume)

        img = nibabel.Nifti1Image(mask, ref.affine)
        img.header.set_xyzt_units(ref.header.get_xyzt_units()[0])
        img.header.set_qform(ref.header.get_qform(), code=2)
        img.update_header()

        nibabel.save(img, files.masks_contr)

    # Run the function
    if not check_volumes_meshed(m2m_folder):
        print('Something went wrong during the meshing.')
        if files.seg_type == 'headreco':
            print("Please delete the m2m folder and "
                  "re-run headreco increasing the '-v' argument")
            print("Type 'headreco all -h' for more informaion")
        else:
            print("Please delete the m2m folder and "
                  "re-run mri2mesh increasing the '--numvertices' argument")
            print("Type 'mri2mesh -h' for more informaion")
        if remove:
            os.remove(files.fnamehead)
    else:
        print('Mesh is fine')


def main():
    args = parse_arguments(sys.argv[1:])
    m2m_dir = os.path.abspath(os.path.realpath(os.path.expanduser(args.m2mpath)))
    if not os.path.isdir(m2m_dir):
        raise IOError('Could not find directory: {0}'.format(args.m2mpath))
    check(m2m_dir, args.remove)


if __name__ == '__main__':
    main()

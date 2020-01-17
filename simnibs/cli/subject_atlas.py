# -*- coding: utf-8 -*-\
'''
    Convert an atlas from FsAverge space to subject space

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
import nibabel

from simnibs import transformations
import simnibs.utils.file_finder as file_finder
from simnibs.utils.simnibs_logger import logger
from simnibs.mesh_tools.mesh_io import read_freesurfer_surface, read_gifti_surface
from simnibs import __version__
import textwrap

def parse_arguments(argv):

    atlasses = textwrap.dedent(
        '''
        'a2009s', 'DK40' or 'HCP_MMP1'

        'a2009s': Destrieux atlas (FreeSurfer v4.5, aparc.a2009s)
        Cite: Destrieux, C. Fischl, B. Dale, A., Halgren, E. A sulcal
        depth-based anatomical parcellation of the cerebral cortex.
        Human Brain Mapping (HBM) Congress 2009, Poster #541

        'DK40': Desikan-Killiany atlas (FreeSurfer, aparc.a2005s)
        Cite: Desikan RS, S�gonne F, Fischl B, Quinn BT, Dickerson BC,
        Blacker D, Buckner RL, Dale AM, Maguire RP, Hyman BT, Albert MS,
        Killiany RJ. An automated labeling system for subdividing the
        human cerebral cortex on MRI scans into gyral based regions of
        interest. Neuroimage. 2006 Jul 1;31(3):968-80.

        'HCP_MMP1': Human Connectome Project (HCP) Multi-Modal Parcellation
        Cite: Glasser MF, Coalson TS, Robinson EC, et al. A multi-modal
        parcellation of human cerebral cortex. Nature. 2016;536(7615):171-178.
        ''')

    def allowed_atlases(s):
        if s not in ['a2009s', 'DK40', 'HCP_MMP1']:
            raise argparse.ArgumentTypeError(atlasses)
        return s


    parser = argparse.ArgumentParser(
        prog="subject_atlas",
        description="Transforms an atlas to subject space. "
                    "Outputs two files: lh/rh.<><SUB_ID>"
    )
    parser.add_argument('-m', "--m2mpath", dest='m2mpath', required=True,
                        help="path to m2m_{subjectID} directory, created in the "
                        "segmentation")
    parser.add_argument("-a", '--atlas', dest='atlas', required=True,
                        help=atlasses, type=allowed_atlases)
    parser.add_argument('-o', '--out_folder', dest='out_folder', default='.',
                        help="Folder where output files will be saved")
    parser.add_argument('--version', action='version', version=__version__)
    return parser.parse_args(argv)


def main():
    args = parse_arguments(sys.argv[1:])
    m2m_dir = os.path.abspath(os.path.realpath(os.path.expanduser(args.m2mpath)))
    if not os.path.isdir(m2m_dir):
        raise IOError('Could not find directory: {0}'.format(args.m2mpath))
    subject_files = file_finder.SubjectFiles(subpath=m2m_dir)
    os.makedirs(args.out_folder, exist_ok=True)
    for hemi in ['lh', 'rh']:
        # I have a parallel implementation here
        fn_atlas = os.path.join(
            file_finder.templates.cat_atlases_surfaces,
            f'{hemi}.aparc_{args.atlas}.freesurfer.annot'
        )
        labels, colors, names = nibabel.freesurfer.io.read_annot(fn_atlas)
        if subject_files.seg_type == 'headreco':
            read_fun = read_gifti_surface
        elif subject_files.seg_type == 'mri2mesh':
            read_fun = read_freesurfer_surface

        if hemi == 'lh':
            ref_surf = read_gifti_surface(file_finder.templates.cat_lh_sphere_ref)
            sub_surf = read_fun(subject_files.lh_reg)
        if hemi == 'rh':
            ref_surf = read_gifti_surface(file_finder.templates.cat_rh_sphere_ref)
            sub_surf = read_fun(subject_files.rh_reg)

        labels_sub, _ = transformations._surf2surf(labels, ref_surf, sub_surf)
        fn_out = os.path.join(
            args.out_folder,
            f'{hemi}.{subject_files.subid}_{args.atlas}.annot'
        )
        logger.info('Writing: ' + fn_out)
        nibabel.freesurfer.io.write_annot(
            fn_out,
            labels_sub,
            colors,
            names,
            fill_ctab=True
        )

if __name__ == '__main__':
    main()

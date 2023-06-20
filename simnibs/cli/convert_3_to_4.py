# -*- coding: utf-8 -*-\
'''
    converts headmodels created by mri2mesh or headreco (simnibs 3) for
    use with simnibs 4
    
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 202 Axel Thielscher

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import glob
import os
import re
import shutil
import numpy as np
import nibabel as nib
import logging
import argparse
import textwrap
import sys

from simnibs import SIMNIBSDIR, __version__
from simnibs.utils import html_writer
from simnibs.utils import cond_utils
from simnibs.mesh_tools import mesh_io
from simnibs.utils import transformations
from simnibs.utils.file_finder import SubjectFiles, templates, ElectrodeCaps_MNI
from simnibs.utils.simnibs_logger import logger, register_excepthook, unregister_excepthook


class SubjectFiles_old:
    ''' Class to find files for a given subject

    Parameters
    ------------
    fnamehead: str (optional)
        Path to head mesh. Either fnamehead or subpath need to be set. subpath has
        preference when it comes to finding other files.

    subpath: str (optional)
        Path to m2m_folder. Either fnamehead or subpath need to be set. subpath has
        preference when it comes to finding other files

    Attributes
    ------------
    fnamehead: str
        Path to head mesh (.msh)

    subpath: str
        Path to m2m_folder (dir)

    subid: str
        The subject ID (eg: ernie)

    basedir: str
        Path to the folder where the m2m_subid folder is located (dir)

    tensor_file: str
        Path to the NifTi file with the tensor conductivity information (.nii.gz)

    eeg_cap_folder: str
        Path to the folder with EEG caps (dir)

    eeg_cap_1010: str
        Path to the EEG 10-10 electrode file (.csv)

    reference_volume: str
        Path to the reference subject volume (.nii.gz)

    mni2conf_nonl: str
        MNI to conform nonlinear transformation (.nii or .nii.gz)

    conf2mni_nonl: str
        Conform to MNI nonlinear tansformation (.nii or .nii.gz)

    mni2conf_6dof: str
        MNI to conform 6 DOF transfomation (.txt or .mat)

    mni2conf_12dof: str
        MNI to conform 12 DOF transfomation (.txt or .mat)

    surf_dir: str
        Directory with surfaces from CAT12/FreeSurfer segmentations (dir)

    seg_type: 'headreco' or 'mri2mesh'
        Type of segmentation (only in mri2mesh and headreco+CAT)

    lh_midgm: str
        Left hemisphere middle gray matter model (only for headreco+CAT segmentations) (.gii)

    rh_midgm: str
        Right hemisphere middle right gray matter model (only for headreco+CAT segmentations) (.gii)
    
    lh_gm: str
        Left hemisphere pial surface model (only for mri2mesh tranformations) (fs surface)

    rh_gm: str
        Right hemisphere pial surface model (only for mri2mesh tranformations) (fs surface)

    lh_wm: str
        Left hemisphere white matter surface model (only for mri2mesh tranformations) (fs surface)

    rh_wm: str
        Right hemisphere white matter surface model (only for mri2mesh tranformations) (fs surface)

    lh_reg: str
        Left hemisphere sphere registration file (only in mri2mesh and headreco+CAT) (fs surface or .gii)

    rh_reg: str
        Right hemisphere sphere registration file (only in mri2mesh and headreco+CAT) (fs surface or .gii)

    final_contr: str
        Volume mask after meshing (.nii.gz)

    masks_contr: str
        Volume mask before meshing (.nii.gz)

    T1: str
        T1 image after applying transformations
    
    Warning
    --------
    This class does not check for existance of the files
    '''

    def __init__(self, fnamehead: str = None, subpath: str = None):
        if not fnamehead and not subpath:
            raise ValueError('Either fnamehear or subpath need to be set')

        if fnamehead:
            if not fnamehead.endswith('.msh'):
                raise IOError('fnamehead must be a gmsh .msh file')
            self.fnamehead = os.path.normpath(os.path.expanduser(fnamehead))
            if not subpath:
                basedir, subid = os.path.split(self.fnamehead)
                self.subid = os.path.splitext(subid)[0]
                self.basedir = os.path.normpath(os.path.expanduser(basedir))
                self.subpath = os.path.join(self.basedir, 'm2m_' + self.subid)

        if subpath:
            self.subpath = os.path.normpath(os.path.expanduser(subpath))
            folder_name = self.subpath.split(os.sep)[-1]
            try:
                self.subid = re.search('m2m_(.+)', folder_name).group(1)
            except:
                raise IOError('Could not find subject ID from subpath. '
                              'Does the folder have the format m2m_subID?')
            self.basedir = os.path.normpath(os.path.join(self.subpath, '..'))

            if not fnamehead:
                self.fnamehead = os.path.join(self.basedir, self.subid + '.msh')

        self.tensor_file = os.path.join(
            self.basedir, 'd2c_' + self.subid,
            'dti_results_T1space', 'DTI_conf_tensor.nii.gz')

        self.eeg_cap_folder = os.path.join(self.subpath, 'eeg_positions')
        self.eeg_cap_1010 = os.path.join(self.eeg_cap_folder, 'EEG10-10_UI_Jurak_2007.csv')
        
        
        # Stuff for volume transformations
        self.reference_volume = os.path.join(self.subpath, 'T1fs_conform.nii.gz')
        mni_transf = os.path.join(self.subpath, 'toMNI')

        self.mni2conf_nonl = os.path.join(mni_transf, 'MNI2Conform_nonl.nii')
        if os.path.isfile(self.mni2conf_nonl + '.gz'):
            self.mni2conf_nonl += '.gz'

        self.conf2mni_nonl = os.path.join(mni_transf, 'Conform2MNI_nonl.nii')
        if os.path.isfile(self.conf2mni_nonl + '.gz'):
            self.conf2mni_nonl += '.gz'

        self.mni2conf_6dof = os.path.join(mni_transf, 'MNI2conform_6DOF')
        if os.path.isfile(self.mni2conf_6dof + '.txt'):
            self.mni2conf_6dof += '.txt'
        if os.path.isfile(self.mni2conf_6dof + '.mat'):
            self.mni2conf_6dof += '.mat'

        self.mni2conf_12dof = os.path.join(mni_transf, 'MNI2conform_12DOF')
        if os.path.isfile(self.mni2conf_12dof + '.txt'):
            self.mni2conf_12dof += '.txt'
        if os.path.isfile(self.mni2conf_12dof + '.mat'):
            self.mni2conf_12dof += '.mat'


        # Stuff for surface transformations
        headreco_surf_dir = os.path.join(self.subpath, 'segment', 'cat', 'surf')
        mri2mesh_surf_dir = os.path.join(self.basedir, 'fs_' + self.subid, 'surf')

        if os.path.isdir(headreco_surf_dir):
            self.surf_dir = headreco_surf_dir
            self.seg_type = 'headreco'
            self.lh_midgm = os.path.join(headreco_surf_dir, 'lh.central.T1fs_conform.gii')
            self.rh_midgm = os.path.join(headreco_surf_dir, 'rh.central.T1fs_conform.gii')
            self.lh_reg = os.path.join(headreco_surf_dir, 'lh.sphere.reg.T1fs_conform.gii')
            self.rh_reg = os.path.join(headreco_surf_dir, 'rh.sphere.reg.T1fs_conform.gii')

        elif os.path.isdir(mri2mesh_surf_dir):
            self.surf_dir = mri2mesh_surf_dir
            self.seg_type = 'mri2mesh'
            self.lh_gm = os.path.join(mri2mesh_surf_dir, 'lh.pial')
            self.lh_wm = os.path.join(mri2mesh_surf_dir, 'lh.white')
            self.rh_gm = os.path.join(mri2mesh_surf_dir, 'rh.pial')
            self.rh_wm = os.path.join(mri2mesh_surf_dir, 'rh.white')
            self.lh_reg = os.path.join(mri2mesh_surf_dir, 'lh.sphere.reg')
            self.rh_reg = os.path.join(mri2mesh_surf_dir, 'rh.sphere.reg')

        else:
            self.surf_dir = None
            self.seg_type = None

        self.final_contr = os.path.join(
            self.subpath, self.subid + '_final_contr.nii.gz')
        self.masks_contr = os.path.join(
            self.subpath, self.subid + '_masks_contr.nii.gz')
        self.T1 = os.path.join(
            self.subpath, 'T1fs_nu_conform.nii.gz')


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
    register_excepthook(logger)
    
    
def convert_old_new(subpath_old, subpath_new):
    ''' convert old head models to simnibs 4 
    '''
    # get old and new file structure
    sf_old = SubjectFiles_old(subpath = subpath_old)
    if not os.path.exists(sf_old.subpath):
        raise FileNotFoundError('Did not find folder '+sf_old.subpath)
    
    if not subpath_new:
        subpath_new = sf_old.subpath+'_v4'
    else:
        pn = os.path.split(subpath_new)[1]
        if not pn.startswith('m2m_'):
            raise IOError('The new folder name has to start with \'m2m_\'. Now it is: '+subpath_new)    
    sf_new = SubjectFiles(subpath = subpath_new)    
    
    if os.path.abspath(sf_old.subpath) == os.path.abspath(sf_new.subpath):
        raise FileExistsError('The names of the m2m-folders for the old and converted data have to differ')
    
    # create directory and start logging
    if not os.path.exists(sf_new.subpath):
        os.mkdir(sf_new.subpath)
    _setup_logger(sf_new.charm_log)  
    logger.info("Converting from "+sf_old.subpath+" to "+sf_new.subpath)
    if sf_old.seg_type is None:
        logger.error("The head model was not created with mri2mesh or headreco")
        return
    logger.info("The original head model was created by "+sf_old.seg_type)
    
    # Copying MR images and mesh
    logger.info("Copying MR images and mesh")
    # T1
    shutil.copyfile(sf_old.reference_volume, sf_new.reference_volume)
    # T2
    fn_T2old = os.path.join(sf_old.subpath,'T2_conform.nii.gz')
    if os.path.exists(fn_T2old):
        shutil.copyfile(fn_T2old, sf_new.T2_reg)
        logger.info("  T2 found ... copying")
    # MNI trafos
    if not os.path.exists(sf_new.mni_transf_folder):
        os.mkdir(sf_new.mni_transf_folder)
    # use nibabel to ensure .nii.gz ending
    nib.save(nib.load(sf_old.conf2mni_nonl), sf_new.conf2mni_nonl)
    nib.save(nib.load(sf_old.mni2conf_nonl), sf_new.mni2conf_nonl)
    # head mesh
    shutil.copyfile(sf_old.fnamehead,sf_new.fnamehead)
    # DTI tensros
    if os.path.exists(sf_old.tensor_file):
        shutil.copyfile(sf_old.tensor_file, sf_new.tensor_file)
        logger.info("  DTI tensor file found ... copying")
    
    # Copying brain surfaces
    if sf_old.surf_dir is not None:
        logger.info("Copying brain surfaces")
        if not os.path.exists(sf_new.surface_folder):
            os.mkdir(sf_new.surface_folder)
        
        if sf_old.seg_type == 'mri2mesh':
            for hemi in ['lh', 'rh']:
                wm_surface = mesh_io.read_freesurfer_surface(getattr(sf_old, hemi+'_wm'))
                gm_surface = mesh_io.read_freesurfer_surface(getattr(sf_old, hemi+'_gm'))
                middle_surf = mesh_io._middle_surface(wm_surface, gm_surface, 0.5)
                fn_out = os.path.join(sf_new.surface_folder, hemi+'.central.gii')
                mesh_io.write_gifti_surface(middle_surf,fn_out)
                
                reg_surf = mesh_io.read_freesurfer_surface(getattr(sf_old, hemi+'_reg'))
                fn_out = os.path.join(sf_new.surface_folder, hemi+'.sphere.reg.gii')
                mesh_io.write_gifti_surface(reg_surf,fn_out)
                
        elif sf_old.seg_type == 'headreco':
            for hemi in ['lh', 'rh']:
                middle_surf = mesh_io.read_gifti_surface(getattr(sf_old, hemi+'_midgm'))
                fn_out = os.path.join(sf_new.surface_folder, hemi+'.central.gii')
                mesh_io.write_gifti_surface(middle_surf,fn_out)
                
                reg_surf = mesh_io.read_gifti_surface(getattr(sf_old, hemi+'_reg'))
                fn_out = os.path.join(sf_new.surface_folder, hemi+'.sphere.reg.gii')
                mesh_io.write_gifti_surface(reg_surf,fn_out)
    
    # Transforming EEG positions
    logger.info("Transforming EEG positions")
    skin_tag = 1005
    final_mesh = mesh_io.read_msh(sf_new.fnamehead)
    idx = (final_mesh.elm.elm_type == 2) & (final_mesh.elm.tag1 == skin_tag)
    mesh = final_mesh.crop_mesh(elements=final_mesh.elm.elm_number[idx])
    
    if not os.path.exists(sf_new.eeg_cap_folder):
        os.mkdir(sf_new.eeg_cap_folder)
    
    cap_files = glob.glob(os.path.join(ElectrodeCaps_MNI, "*.csv"))
    for fn in cap_files:
        fn_out = os.path.splitext(os.path.basename(fn))[0]
        fn_out = os.path.join(sf_new.eeg_cap_folder, fn_out)
        transformations.warp_coordinates(
            fn,
            sf_new.subpath,
            transformation_direction="mni2subject",
            out_name=fn_out + ".csv",
            out_geo=fn_out + ".geo",
            mesh_in=mesh,
        )
                
    v = final_mesh.view(cond_list=cond_utils.standard_cond())
    v.write_opt(sf_new.fnamehead)
    
    # Write label image from mesh
    logger.info("Writing label image from mesh")
    MNI_template = templates.mni_volume
    mesh = final_mesh.crop_mesh(elm_type=4)
    field = mesh.elm.tag1.astype(np.uint16)
    ed = mesh_io.ElementData(field)
    ed.mesh = mesh
    ed.to_deformed_grid(
        sf_new.mni2conf_nonl,
        MNI_template,
        out=sf_new.final_labels_MNI,
        out_original=sf_new.final_labels,
        method="assign",
        reference_original=sf_new.reference_volume,
    )
    
    fn_LUT = sf_new.final_labels.rsplit(".", 2)[0] + "_LUT.txt"
    shutil.copyfile(templates.final_tissues_LUT, fn_LUT)
    
    # Write report
    logger.info("Writing report")
    fn_settings = os.path.join(SIMNIBSDIR, "charm.ini")
    shutil.copyfile(fn_settings,sf_new.settings)
    html_writer.write_report(sf_new)
    os.remove(sf_new.settings)
    
    # Stop logging ...
    logger.info("Conversion finished")
    while logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])
    unregister_excepthook()
    logging.shutdown()
    with open(sf_new.charm_log, "a") as f:
        f.write("</pre></BODY></HTML>")
        f.close()
        

def parseArguments(argv):
    usage_text = textwrap.dedent('''
                                 
Converts head models created by mri2mesh or headreco for use in simnibs 4:
    convert_3_to_4 m2m_old m2m_new
            
Note:
    m2m_new can be left out, then the new folder name will be m2m_old_v4

    ''')

    parser = argparse.ArgumentParser(
        prog="convert_3_to_4",usage=usage_text)
    parser.add_argument('m2m_old', nargs='?', help="original m2m-folder")
    parser.add_argument('m2m_new', nargs='?', help="name of converted m2m-folder")
    parser.add_argument('-v','--version', action='version', version=__version__)
    
    args=parser.parse_args(argv)
    if args.m2m_old is None:
        parser.print_help()
        exit()
    return args


def main():
    args = parseArguments(sys.argv[1:])
    convert_old_new(args.m2m_old, args.m2m_new)
    

if __name__ == '__main__':
    main()
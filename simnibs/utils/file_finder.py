# -*- coding: utf-8 -*-\
'''
    Find templates and segmentation files of SimNIBS
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2019 Guilherme B Saturnino

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
import sys
import os
import re

from .. import SIMNIBSDIR

class Templates:
    ''' Defines the Templates for file names used in SimNIBS

    Attributes
    ------------
    cat_templates_surfaces: str
        Path to the directory with the CAT12 template surfaces (dir)
    cat_atlases_surfaces: str
        Path to the directory with the CAT12 atlases surfaces (dir)
    cat_lh_cortex_ref: str
        Path to the GifTi surface with the FsAverage lh cortex template (.gii)
    cat_rh_cortex_ref: str
        Path to the GifTi surface with the FsAverage lh cortex template (.gii)
    cat_lh_sphere_ref: str
        Path to the GifTi surface with the FsAverage lh sphere template (.gii)
    cat_rh_sphere_ref: str
        Path to the GifTi surface with the FsAverage lh sphere template (.gii)
    mni_volume: str
        Path to the NifTi volume with the MNI template (T1, 1mm) (.nii.gz)
    freesurfer_templates: str
        Path to the folder with FreeSurfer templates (dir)
    fs_lh_sphere_ref: str
        Path to the fs surface file with the reference lh sphere (freesurfer surface)
    fs_rh_sphere_ref: str
        Path to the fs surface file with the reference rh sphere (freesurfer surface)
    fs_lh_cortex_ref: str
        Path to the fs surface file with the Fsavarage ls pial (freesurfer surface)
    fs_rh_cortex_ref: str
        Path to the fs surface file with the Fsavarage rs pial (freesurfer surface)
    '''
    def __init__(self):
        self._resources = os.path.join(SIMNIBSDIR, 'resources')
        self._spm12 = os.path.join(self._resources, 'spm12')
        # Cat 12 things
        self.cat_templates_surfaces = os.path.join(
            self._spm12, 'toolbox', 'cat12', 'templates_surfaces')
        self.cat_atlases_surfaces = os.path.join(
            self._spm12, 'toolbox', 'cat12', 'atlases_surfaces')
        self.cat_lh_cortex_ref = os.path.join(
            self.cat_templates_surfaces, 'lh.central.freesurfer.gii')
        self.cat_rh_cortex_ref = os.path.join(
            self.cat_templates_surfaces, 'rh.central.freesurfer.gii')
        self.cat_lh_sphere_ref = os.path.join(
            self.cat_templates_surfaces, 'lh.sphere.freesurfer.gii')
        self.cat_rh_sphere_ref = os.path.join(
            self.cat_templates_surfaces, 'rh.sphere.freesurfer.gii')
        # MNI
        self.mni_volume = os.path.join(
            self._resources, 'templates', 'MNI152_T1_1mm.nii.gz')
        # FreeSurfer
        try:
            os.environ['FREESURFER_HOME']
        except KeyError:
            self.freesurfer_templates = None
            self.fs_lh_sphere_ref = None
            self.fs_rh_sphere_ref = None
            self.fs_lh_cortex_ref = None
            self.fs_rh_cortex_ref = None
        else:
            self.freesurfer_templates = os.path.join(
                os.environ['FREESURFER_HOME'], 'subjects', 'fsaverage', 'surf')
            self.fs_lh_sphere_ref = os.path.join(self.freesurfer_templates, 'lh.sphere.reg')
            self.fs_rh_sphere_ref = os.path.join(self.freesurfer_templates, 'rh.sphere.reg')

            self.fs_lh_cortex_ref = os.path.join(self.freesurfer_templates, 'lh.pial')
            self.fs_rh_cortex_ref = os.path.join(self.freesurfer_templates, 'lh.pial')

templates = Templates()


class SubjectFiles:
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

    ref_fs: str
        Reference FreeSurfer space file (.nii.gz)

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
        self.eeg_cap_1010 = self.get_eeg_cap()

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

        self.ref_fs = os.path.join(self.subpath, 'ref_FS.nii.gz')
        headreco_suf_dir = os.path.join(self.subpath, 'segment', 'cat', 'surf')
        mri2mesh_surf_dir = os.path.join(self.basedir, 'fs_' + self.subid, 'surf')

        if os.path.isdir(headreco_suf_dir):
            self.suf_dir = headreco_suf_dir
            self.seg_type = 'headreco'
            self.lh_midgm = os.path.join(headreco_suf_dir, 'lh.central.T1fs_conform.gii')
            self.rh_midgm = os.path.join(headreco_suf_dir, 'rh.central.T1fs_conform.gii')
            self.lh_reg = os.path.join(headreco_suf_dir, 'lh.sphere.reg.T1fs_conform.gii')
            self.rh_reg = os.path.join(headreco_suf_dir, 'rh.sphere.reg.T1fs_conform.gii')

        elif os.path.isdir(mri2mesh_surf_dir):
            self.suf_dir = mri2mesh_surf_dir
            self.seg_type = 'mri2mesh'
            self.lh_gm = os.path.join(mri2mesh_surf_dir, 'lh.pial')
            self.lh_wm = os.path.join(mri2mesh_surf_dir, 'lh.white')
            self.rh_gm = os.path.join(mri2mesh_surf_dir, 'rh.pial')
            self.rh_wm = os.path.join(mri2mesh_surf_dir, 'rh.white')
            self.lh_reg = os.path.join(mri2mesh_surf_dir, 'lh.sphere.reg')
            self.rh_reg = os.path.join(mri2mesh_surf_dir, 'rh.sphere.reg')

        else:
            self.suf_dir = None
            self.seg_type = None

        self.final_contr = os.path.join(
            self.subpath, self.subid + '_final_contr.nii.gz')
        self.masks_contr = os.path.join(
            self.subpath, self.subid + '_masks_contr.nii.gz')

    def get_eeg_cap(self, cap_name: str = 'EEG10-10_UI_Jurak_2007.csv') -> str:
        ''' Gets the name of an EEG cap for this subject

        Parameters
        -----------
        cap_name: str (optional)
            Name of cap file. Default: 'EEG10-10_UI_Jurak_2007.csv'

        Returns
        --------
        cap_fn: str
            Path to the cap file

        Warning
        --------
        This does not check for existance of the file

        '''
        return os.path.join(self.eeg_cap_folder, cap_name)


def path2bin(program):
    """Return the full path to a specified program contained within the SIMNIBS
    binaries folder.

    PARAMETERS
    ----------
    program : str
        Name of the executable for which to return the full path. Must be in
        the binary directory in SIMNIBS.

    RETURNS
    ----------
    path to binary : str
        Full path to the binary.
    """
    # input must be string
    assert(type(program) is str)

    # get path to SIMNIBS
    if sys.platform == 'win32':
        p = 'win'
    elif sys.platform == 'linux':
        p = 'linux'
    elif sys.platform == 'darwin':
        p = 'osx'
    else:
        raise OSError('OS not supported!')

    path_to_binary = os.path.join(SIMNIBSDIR, 'bin', p, program)
    return f'"{path_to_binary}"'



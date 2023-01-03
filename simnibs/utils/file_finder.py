# -*- coding: utf-8 -*-\
'''
    Find templates and segmentation files of SimNIBS
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2019, 2020 Guilherme B Saturnino

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
from pathlib import Path
import sys
import os
import re
import glob
import collections
import numpy as np
import nibabel

from .. import SIMNIBSDIR

__all__ = ['templates', 'get_atlas', 'get_reference_surf', 'SubjectFiles', 'coil_models']

# TODO: add electrode mask here

HEMISPHERES = {"lh", "rh"}
VALID_HEMI = HEMISPHERES.union({"both"})

class Templates:
    ''' Defines the Templates for file names used in SimNIBS

    Attributes
    ------------
    atlases_surfaces: str
        Path to the directory with the atlases surfaces (dir)
    mni_volume: str
        Path to the NifTi volume with the MNI template (T1, 1mm) (.nii.gz)
    freesurfer_templates: str
        Path to the folder with FreeSurfer templates (dir)
    simnibs_logo: str
        Path to the SimNIBS logo stored as triangle mesh (gmsh format)
    charm_atlas_path: str
        Path to the charm atlas folders
    labeling_LUT: str
        Path to freeview color LUT for the labeling.nii.gz
    final_tissues_LUT: str
        Path to freeview color LUT for the final_tissues.nii.gz and
        label_prep/tissue_labeling_upsampled.nii.gz
    '''
    def __init__(self):
        self._resources = os.path.join(SIMNIBSDIR, 'resources')

        # atlases in fsaverage space
        self.atlases_surfaces = os.path.join(
            self._resources, 'templates', 'fsaverage_atlases')

        # MNI
        self.mni_volume = os.path.join(
            self._resources, 'templates', 'MNI152_T1_1mm.nii.gz')

        # Electrode mask of valid skin region in MNI space
        self.mni_volume_upper_head_mask = os.path.join(
            self._resources, 'templates', 'MNI152_T1_1mm_upper_head_mask.nii.gz')

        # path to fsaverage surfaces
        self.freesurfer_templates = os.path.join(
            self._resources, 'templates', 'fsaverage_surf')

        # SimNIBS logo
        self.simnibs_logo = os.path.join(self._resources, 'simnibslogo.msh')

        # labeling_LUT
        self.labeling_LUT = os.path.join(self._resources, 'labeling_FreeSurferColorLUT.txt')

        # final_tissues_LUT
        self.final_tissues_LUT = os.path.join(self._resources, 'final_tissues_FreeSurferColorLUT.txt')

        #CHARM atlas path
        self.charm_atlas_path = os.path.join(SIMNIBSDIR, 'segmentation','atlases')

        #Viewer templates
        self.html_template = os.path.join(SIMNIBSDIR, '_internal_resources', 'html', 'template.html')
        self.jquery = os.path.join(SIMNIBSDIR, '_internal_resources', 'html', 'jquery.min.js')
        self.brainsprite = os.path.join(SIMNIBSDIR, '_internal_resources', 'html', 'brainsprite.min.js')


templates = Templates()
coil_models = os.path.join(SIMNIBSDIR, 'resources', 'coil_models')
ElectrodeCaps_MNI = os.path.join(SIMNIBSDIR, 'resources', 'ElectrodeCaps_MNI')

def get_atlas(atlas_name, hemi='both'):
    ''' Loads a brain atlas based of the FreeSurfer fsaverage template

    Parameters
    -----------
    atlas_name: 'a2009s', 'DK40' or 'HCP_MMP1'
            Name of atlas to load

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

    hemi (optional): 'lh', 'rh' or 'both'
        Hemisphere to use. In the case of 'both', will assume that left hemisphere
        nodes comes before right hemisphere nodes

    Returns
    ---------
    atlas: dict
        Dictionary where atlas['region'] = roi
    '''
    if atlas_name not in ['a2009s', 'DK40', 'HCP_MMP1']:
        raise ValueError('Invalid atlas name')


    if hemi in ['lh', 'rh']:
        fn_atlas = os.path.join(
            templates.atlases_surfaces,
            f'{hemi}.aparc_{atlas_name}.freesurfer.annot'
        )
        labels, _ , names = nibabel.freesurfer.io.read_annot(fn_atlas)
        atlas = {}
        for l, name in enumerate(names):
            atlas[name.decode()] = labels == l

        return atlas
    # If both hemispheres
    elif hemi == 'both':
        atlas_lh = get_atlas(atlas_name, 'lh')
        atlas_rh = get_atlas(atlas_name, 'rh')
        atlas = {}
        pad_rh = np.zeros_like(list(atlas_rh.values())[0])
        pad_lh = np.zeros_like(list(atlas_lh.values())[0])
        for name, mask in atlas_lh.items():
            atlas[f'lh.{name}'] = np.append(mask, pad_rh)  # pad after
        for name, mask in atlas_rh.items():
            atlas[f'rh.{name}'] = np.append(pad_lh, mask)  # pad after

        return atlas
    else:
        raise ValueError('Invalid hemisphere name')


def get_reference_surf(region, surf_type='central'):
    ''' Gets the file name of a reference surface

    Parameters
    -----------
    region: 'lh', 'rh', 'lc' or 'rc'
        Name of the region of interest
    surf_type: 'central', 'sphere', 'inflated' (optional)
        Surface type. Default: central

    Returns
    --------
    fn_surf: str
        Name of surface file

    Raises
    -------
    FileNotFoundError if the specified reference surface is not found

    '''
    fn_surf = os.path.join(
        SIMNIBSDIR, 'resources', 'templates',
        'fsaverage_surf', f'{region}.{surf_type}.freesurfer.gii'
    )
    if os.path.isfile(fn_surf):
        return fn_surf
    else:
        raise FileNotFoundError('Could not find reference surface')

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

    tensor_file: str
        Path to the NifTi file with the tensor conductivity information (.nii.gz)

    eeg_cap_folder: str
        Path to the folder with EEG caps (dir)

    segmentation_folder: str
        Path to the output from the segmentation

    surface_folder: str
        Path to surfaces from middle GM reconstruction

    label_prep_folder: str
        Path to the output from upsampling

    mni_transf_folder: str
        Path to MNI transformations

    eeg_cap_1010: str
        Path to the EEG 10-10 electrode file (.csv)

    reference_volume: str
        Path to the reference subject volume (T1.nii.gz)

    mni2conf_nonl: str
        MNI to conform nonlinear transformation (.nii.gz)

    conf2mni_nonl: str
        Conform to MNI nonlinear tansformation (.nii.gz)

    mni2conf_6dof: str
        MNI to conform 6 DOF transfomation (.txt or .mat)

    mni2conf_12dof: str
        MNI to conform 12 DOF transfomation (.txt or .mat)

    final_labels_MNI: str
        Label image created from final mesh in MNI space

    ref_fs: str
        Reference FreeSurfer space file (.nii.gz)
        Now always set to True, so that a standard header is added, which seems to work

    central_surfaces: list
        List of SurfaceFile objects which containts 2 fields:
            fn: name of surface file (.gii format)
            region: 'lh', 'rh', 'lc' or 'rc'

    central_surfaces: list
        Same as above but for the pial surfaces

    sphere_reg_surfaces: list
        Same as above but for the spherical registration files

    thickness: list
        Cortical thicknesses are stored also as SurfaceFile

    regions: list
        list of region names (e.g. 'lh', 'rh') where all surfaces above are present

    T1: str
        T1 image after applying transformations

    T2_reg: str
        T2 image after rigid co-registration to T1

    T1_denoised: str
        T1 image after running filtering (optional)

    T2_reg_denoised: str
        T2 image after co-registration and filtering (optional)

    T1_bias_corrected: str
        T1 image after bias correction

    T2_reg_bias_corrected: str
        T2 image after bias correction

    labeling: str
        Output segmentation from samseg

    final_labels: str
        Label image created from final mesh

    template_coregistered: str
        Affine mapping from atlas voxel space to T1 voxel space

    T1_upsampled: str
        Bias corrected T1 upsampled to 0.5mm^3 resolution

    T2_upsampled: str
        Bias corrected T2 upsampled to 0.5mm^3 resolution

    tissue_labeling_upsampled: str
        Labeling mapped to conductivity values and reconstructed at 0.5mm^3 resolution

    head_mesh: str
        Final head mesh

    cerebrum_mask: str
        Mask indicating cerebrum and cerebellum. Needed for surfaces.

    norm_image: str
        Normalized intensity image. Needed for surfaces.

    subcortical_mask: str
        Mask of subcortical structures. Needed for surfaces

    parahippo_mask: str
        Mask of parahippocampal area. Needed for surfaces.

    hemi_mask: str
        Mask indicating left/right. Needed for surfaces.

    charm_log: str
        The charm run log (.html)

    summary_report: str
        html summarizing the run and linking to the reports (.html)

    Warning
    --------
    This class does not check for existence of the files
    '''

    def __init__(self, fnamehead: str = None, subpath: str = None):
        if not fnamehead and not subpath:
            raise ValueError('Either fnamehead or subpath need to be set')

        if fnamehead:
            if not fnamehead.endswith('.msh'):
                raise IOError('fnamehead must be a gmsh .msh file')
            self.fnamehead = os.path.normpath(os.path.expanduser(fnamehead))
            if not subpath:
                subpath, subid = os.path.split(self.fnamehead)
                self.subid = os.path.splitext(subid)[0]
                self.subpath = os.path.normpath(os.path.expanduser(subpath))
                if not re.search('m2m_(.+)', self.subpath):
                    # some mesh without m2m_folder
                    subpath = None
                    self.subpath = None

        if subpath:
            self.subpath = os.path.normpath(os.path.expanduser(subpath))
            folder_name = self.subpath.split(os.sep)[-1]
            try:
                self.subid = re.search('m2m_(.+)', folder_name).group(1)
            except:
                raise IOError('Could not find subject ID from subpath. '
                              'Does the folder have the format m2m_subID?')
            if not fnamehead:
                self.fnamehead = os.path.join(self.subpath, self.subid + '.msh')

        # otherwise SESSION._prepare fails for meshes
        # that do not have a m2m_folder
        if not self.subpath:
            self.subpath = ''

        self.tensor_file = os.path.join(self.subpath, 'DTI_coregT1_tensor.nii.gz')

        self.eeg_cap_folder = os.path.join(self.subpath, 'eeg_positions')
        self.eeg_cap_1010 = self.get_eeg_cap()

        self.segmentation_folder = os.path.join(self.subpath, 'segmentation')
        self.surface_folder = os.path.join(self.subpath, 'surfaces')
        self.label_prep_folder = os.path.join(self.subpath, 'label_prep')

        # Stuff for volume transformations
        self.reference_volume = os.path.join(self.subpath, 'T1.nii.gz')
        self.mni_transf_folder = os.path.join(self.subpath, 'toMNI')

        self.mni2conf_nonl = os.path.join(self.mni_transf_folder, 'MNI2Conform_nonl.nii.gz')
        self.conf2mni_nonl = os.path.join(self.mni_transf_folder, 'Conform2MNI_nonl.nii.gz')
        self.final_labels_MNI = os.path.join(self.mni_transf_folder, 'final_tissues_MNI.nii.gz')

        self.mni2conf_6dof = os.path.join(self.mni_transf_folder, 'MNI2conform_6DOF')
        if os.path.isfile(self.mni2conf_6dof + '.txt'):
            self.mni2conf_6dof += '.txt'
        if os.path.isfile(self.mni2conf_6dof + '.mat'):
            self.mni2conf_6dof += '.mat'

        self.mni2conf_12dof = os.path.join(self.mni_transf_folder, 'MNI2conform_12DOF')
        if os.path.isfile(self.mni2conf_12dof + '.txt'):
            self.mni2conf_12dof += '.txt'
        if os.path.isfile(self.mni2conf_12dof + '.mat'):
            self.mni2conf_12dof += '.mat'

        # Stuff for surface transformations
        # Organize the files in 3 separate lists
        SurfaceFile = collections.namedtuple('SurfaceFile', ['fn', 'region',
                                                             'subsampling'])
        # Look for all .gii files in surface_folder
        surfaces = glob.glob(os.path.join(self.surface_folder, '*.gii'))

        self.sphere_reg_surfaces = []
        #self.sphere_surfaces = []
        self.central_surfaces = []
        self.pial_surfaces = []
        for fn in surfaces:
            s = os.path.basename(fn)
            region = s[:2]
            if '.sphere.reg' in s:
                # Assumes that the filename is
                # hemi.sphere.reg[.subsampling].gii
                subsampling = _get_subsampling(s.split('.')[3:])
                self.sphere_reg_surfaces.append(
                    SurfaceFile(fn, region, subsampling)
                )
            elif '.pial.' in s:
                # Assumes that the filename is hemi.pial[.subsampling].gii
                subsampling = _get_subsampling(s.split('.')[2:])
                self.pial_surfaces.append(
                     SurfaceFile(fn, region, subsampling)
                 )
            elif '.central.' in s:
                # Assumes that the filename is hemi.central[.subsampling].gii
                subsampling = _get_subsampling(s.split('.')[2:])
                self.central_surfaces.append(
                    SurfaceFile(fn, region, subsampling)
                )

        surfaces = glob.glob(os.path.join(self.surface_folder, '*.thickness'))
        self.thickness = []
        for fn in surfaces:
            s = os.path.basename(fn)
            region = s[:2]
            self.thickness.append(
                    SurfaceFile(fn, region, None)
                )

        self.regions = sorted(
            set([s.region for s in self.sphere_reg_surfaces]) &
            #set([s.region for s in self.pial_surfaces]) &
            set([s.region for s in self.central_surfaces])
        )

        self.T2_reg = os.path.join(self.subpath, 'T2_reg.nii.gz')
        self.T1_denoised = os.path.join(self.segmentation_folder, 'T1_denoised.nii.gz')
        self.T2_reg_denoised = os.path.join(self.segmentation_folder, 'T2_reg_denoised.nii.gz')
        self.T1_bias_corrected = os.path.join(self.segmentation_folder, 'T1_bias_corrected.nii.gz')
        self.T2_bias_corrected = os.path.join(self.segmentation_folder, 'T2_bias_corrected.nii.gz')
        self.labeling = os.path.join(self.segmentation_folder, 'labeling.nii.gz')
        self.final_labels = os.path.join(self.subpath, 'final_tissues.nii.gz')
        self.template_coregistered = os.path.join(self.segmentation_folder, 'template_coregistered.mgz')
        self.T1_upsampled = os.path.join(self.label_prep_folder,'T1_upsampled.nii.gz')
        self.T2_upsampled = os.path.join(self.label_prep_folder,'T2_upsampled.nii.gz')
        self.tissue_labeling_upsampled = os.path.join(self.label_prep_folder,'tissue_labeling_upsampled.nii.gz')
        self.tissue_labeling_before_morpho = os.path.join(self.label_prep_folder, 'before_morpho.nii.gz')
        self.upper_mask = os.path.join(self.label_prep_folder, 'upper_part.nii.gz')
        self.settings = os.path.join(self.subpath, 'settings.ini')
        self.cereb_mask = os.path.join(self.surface_folder, 'cereb_mask.nii.gz')
        self.norm_image = os.path.join(self.surface_folder, 'norm_image.nii.gz')
        self.subcortical_mask = os.path.join(self.surface_folder, 'subcortical_mask.nii.gz')
        self.hemi_mask = os.path.join(self.surface_folder, 'hemi_mask.nii.gz')
        self.charm_log = os.path.join(self.subpath, 'charm_log.html')
        self.summary_report = os.path.join(self.subpath, 'charm_report.html')
        #self.ref_fs = os.path.join(self.subpath, 'ref_FS.nii.gz')
        self.ref_fs = True # when True, mesh_io.write_freesurfer_surface writes a standard header that seems to work

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
        This does not check for existence of the file

        '''
        return os.path.join(self.eeg_cap_folder, cap_name)

    def get_surface(self, region, surf_type='central', subsampling=None):
        ''' Gets the filename of a subject CAT12 surface

        Parameters
        -----------
        region: 'lh', 'rh', 'lc' or 'rc'
            Name of the region of interest
        surf_type: 'central', 'sphere_reg', 'pial', 'thickness' (optional)
            Surface type. Default: central
        subsampling : int | None
            The subsampling of the surface. None corresponds to the original
            (full resolution) surface (default = None).

        Returns
        --------
        fn_surf: str
            Name of surface file

        Raises
        -------
        FileNotFoundError if the specified reference surface is not found

        '''
        if surf_type == 'central':
            for s in self.central_surfaces:
                if s.region == region and s.subsampling == subsampling:
                    return s.fn
        elif surf_type == 'sphere_reg':
            for s in self.sphere_reg_surfaces:
                if s.region == region and s.subsampling == subsampling:
                    return s.fn
        elif surf_type == 'pial':
            for s in self.pial_surfaces:
                if s.region == region and s.subsampling == subsampling:
                    return s.fn
        elif surf_type == 'thickness':
            for s in self.thickness:
                if s.region == region and s.subsampling == subsampling:
                    return s.fn
        else:
            raise ValueError('invalid surf_type')
        raise FileNotFoundError('Could not find surface')

    #### START: new get surface methods ####
    #### Remove "v2v2_" if keep
    def v2v2_get_surface_file(self, surf, subsampling=None, hemi="both"):
        """Get surface files, e.g., central, pial, sphere, sphere.reg
        """
        subsampling, hemi = self._parse_surface_args(subsampling, hemi)
        return {h: Path(self.surface_folder) / subsampling / f"{h}.{surf}.gii" for h in hemi}


    def v2v2_get_morph_data_file(self, data, subsampling=None, hemi="both"):
        """Get morphometry data files, e.g., thickness.
        """
        subsampling, hemi = self._parse_surface_args(subsampling, hemi)
        return {h: Path(self.surface_folder) / subsampling / f"{h}.{data}" for h in hemi}


    def v2v2_read_surface(self, surf, subsampling=None, hemi="both"):
        return {h: self.v2v2_read_gifti_to_dict(f) for h, f in self.v2v2_get_surface_file(surf, subsampling, hemi).items()}


    def v2v2_read_morph_data(self, data, subsampling=None, hemi="both"):
        return {h: nibabel.freesurfer.read_morph_data(f) for h, f in self.v2v2_get_morph_data_file(data, subsampling, hemi).items()}


    def v2v2_write_surface(self, surf_dict, surf, subsampling=None):
        hemi = list(surf_dict.keys())
        hemi = "both" if len(hemi) == 2 else hemi[0]
        files = self.v2v2_get_surface_file(surf, subsampling, hemi)
        for h in surf_dict:
            if not files[h].parent.exists():
                files[h].parent.mkdir()
            self.v2v2_write_dict_to_gifti(files[h], surf_dict[h])


    def v2v2_write_morph_data(self, data_dict, data, subsampling=None):
        hemi = list(data_dict.keys())
        hemi = "both" if len(hemi) == 2 else hemi[0]
        files = self.v2v2_get_morph_data_file(data, subsampling, hemi)
        for h in data_dict:
            nibabel.freesurfer.write_morph_data(files[h], data_dict[h])


    @staticmethod
    def v2v2_write_dict_to_gifti(filename, surf):
        darrays = (
            nibabel.gifti.gifti.GiftiDataArray(surf['points'], 'pointset'),
            nibabel.gifti.gifti.GiftiDataArray(surf['tris'], 'triangle')
        )
        gii = nibabel.GiftiImage(darrays=darrays)
        gii.to_filename(filename)


    @staticmethod
    def v2v2_read_gifti_to_dict(filename):
        gii = nibabel.load(filename)
        return dict(
            points=gii.agg_data("pointset"), tris=gii.agg_data("triangle")
        )


    @staticmethod
    def _parse_surface_args(subsampling, hemi):
        assert hemi in VALID_HEMI, f"Invalid hemisphere argument. Please choose one of {VALID_HEMI}."
        hemi = HEMISPHERES if hemi == 'both' else [hemi]
        subsampling = "" if subsampling is None else str(subsampling)
        return subsampling, hemi
    #### END: new get surface methods ####


def _get_subsampling(parts):
    """Get subsampling if `parts` is of the format (subsampling, gii).
    """
    if len(parts) == 2:
        # If parts[0] cannot be interpreted as an int we will get an error
        # here. In that case, assume no subsampling has been performed.
        try:
            subsampling = int(parts[0])
        except ValueError:
            subsampling = None
    else:
        subsampling = None
    return subsampling

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

    path_to_binary = os.path.join(SIMNIBSDIR, 'external', 'bin', p, program)

    return path_to_binary



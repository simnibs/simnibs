# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 17:12:24 2020

@author: axthi
"""


import logging
import os
import shutil
import time
import subprocess
import nibabel as nib
import glob
import sys
import tempfile
import numpy as np
from functools import partial
from scipy import ndimage
import scipy.ndimage.morphology as mrph
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import affine_transform
from scipy.ndimage.measurements import label

from simnibs import SIMNIBSDIR

from .. import __version__
from . import samseg
from ._thickness import _calc_thickness
from ._cat_c_utils import sanlm
from .brain_surface import mask_from_surface
from .. import utils
from ..utils.simnibs_logger import logger
from ..utils import file_finder
from ..utils import transformations
from ..utils.transformations import resample_vol, crop_vol
from ..utils.spawn_process import spawn_process
from ..mesh_tools.meshing import create_mesh
from ..mesh_tools.mesh_io import read_gifti_surface, write_msh, read_off, write_off, Msh, ElementData
from ..simulation import cond
from ..utils import plotting

def _register_atlas_to_input_affine(T1, template_file_name,
                                    affine_mesh_collection_name, mesh_level1,
                                    mesh_level2, save_path,
                                    template_coregistered_name,
                                    init_atlas_settings, neck_tissues,
                                    visualizer,
                                    noneck,
                                    init_transform=None,
                                    world_to_world_transform_matrix=None):

    # Import the affine registration function
    scales = init_atlas_settings['affine_scales']
    thetas = init_atlas_settings['affine_rotations']
    horizontal_shifts = init_atlas_settings['affine_horizontal_shifts']
    vertical_shifts = init_atlas_settings['affine_vertical_shifts']
    thetas_rad = [theta * np.pi/180 for theta in thetas]
    neck_search_bounds = init_atlas_settings['neck_search_bounds']
    ds_factor = init_atlas_settings['dowsampling_factor_affine']

    affine = samseg.AffineWholeHead(T1, affine_mesh_collection_name,
                                    template_file_name)

    init_options = samseg.initializationOptions(pitchAngles=thetas_rad,
                                                scales=scales,
                                                scalingCenter=[0.0, -120.0, 0.0],
                                                horizontalTableShifts=horizontal_shifts,
                                                verticalTableShifts=vertical_shifts)

    image_to_image_transform, world_to_world_transform, optimization_summary =\
    affine.registerAtlas(initializationOptions=init_options,
                         targetDownsampledVoxelSpacing=ds_factor,
                         visualizer=visualizer,
                         noneck=noneck)

    affine.saveResults(T1, template_file_name, save_path,
                       template_coregistered_name, image_to_image_transform,
                       world_to_world_transform)
    
    logger.info('Template registration summary.')
    logger.info('Number of Iterations: %d, Cost: %f\n' %
                (optimization_summary['numberOfIterations'],
                 optimization_summary['cost']))

    if not noneck:
        logger.info('Adjusting neck.')
        exitcode = affine.adjust_neck(T1, template_coregistered_name, mesh_level1,
                                      mesh_level2, neck_search_bounds, neck_tissues,
                                      visualizer, downsampling_target=2.0)
        if exitcode == -1:
            file_path = os.path.split(template_coregistered_name)
            shutil.copy(mesh_level1, os.path.join(file_path[0], 'atlas_level1.txt.gz'))
            shutil.copy(mesh_level2, os.path.join(file_path[0], 'atlas_level2.txt.gz'))
        else:
            logger.info('Neck adjustment done.')
    else:
        logger.info('No neck, copying meshes over.')
        file_path = os.path.split(template_coregistered_name)
        shutil.copy(mesh_level1, os.path.join(file_path[0], 'atlas_level1.txt.gz'))
        shutil.copy(mesh_level2, os.path.join(file_path[0], 'atlas_level2.txt.gz'))

def _denoise_input_and_save(input_name, output_name):
    input_raw = nib.load(input_name)
    img = input_raw.get_data()
    # Sometimes the images have an extra unit dimension,
    # squeeze that out if it's there.
    img = img.squeeze()
    img_smoothed = sanlm(img, 3, 1)
    output_smoothed = nib.Nifti1Image(img_smoothed, input_raw.affine)
    nib.save(output_smoothed, output_name)


def _estimate_parameters(path_to_segment_folder,
                         template_coregistered_name,
                         path_to_atlas_folder, input_images,
                         segment_settings, gmm_parameters, visualizer):

    ds_targets = segment_settings['downsampling_targets']
    kernel_size = segment_settings['bias_kernel_width']
    bg_mask_sigma = segment_settings['background_mask_sigma']
    bg_mask_th = segment_settings['background_mask_threshold']
    stiffness = segment_settings['mesh_stiffness']
    covariances = segment_settings['diagonal_covariances']
    shared_gmm_parameters = samseg.io.kvlReadSharedGMMParameters(gmm_parameters)
    user_optimization_options = {'multiResolutionSpecification':
                                 [{'atlasFileName':
                                   os.path.join(path_to_segment_folder,
                                                'atlas_level1.txt.gz'),
                                   'targetDownsampledVoxelSpacing': ds_targets[0],
                                   'maximumNumberOfIterations': 100,
                                   'estimateBiasField': True},
                                  {'atlasFileName':
                                      os.path.join(path_to_segment_folder,
                                                   'atlas_level2.txt.gz'),
                                   'targetDownsampledVoxelSpacing': ds_targets[1],
                                   'maximumNumberOfIterations': 100,
                                   'estimateBiasField': True}]}

    user_model_specifications = {'atlasFileName': os.path.join(
            path_to_segment_folder, 'atlas_level2.txt.gz'),
                                 'biasFieldSmoothingKernelSize': kernel_size,
                                 'brainMaskingSmoothingSigma': bg_mask_sigma,
                                 'brainMaskingThreshold': bg_mask_th,
                                 'K': stiffness,
                                 'useDiagonalCovarianceMatrices': covariances,
                                 'sharedGMMParameters': shared_gmm_parameters}

    samseg_kwargs = dict(
        imageFileNames=input_images,
        atlasDir=path_to_atlas_folder,
        savePath=path_to_segment_folder,
        transformedTemplateFileName=template_coregistered_name,
        userModelSpecifications=user_model_specifications,
        userOptimizationOptions=user_optimization_options,
        visualizer=visualizer,
        saveHistory=False,
        saveMesh=False,
        savePosteriors=False,
        saveWarp=True)

    logger.info('Starting segmentation.')
    samsegment = samseg.SamsegWholeHead(**samseg_kwargs)
    samsegment.preProcess()
    samsegment.process()

    # Print optimization summary
    optimizationSummary = samsegment.getOptimizationSummary()
    for multiResolutionLevel, item in enumerate(optimizationSummary):
        logger.info('atlasRegistrationLevel%d %d %f\n' %
                    (multiResolutionLevel, item['numberOfIterations'],
                     item['perVoxelCost']))

    return samsegment.saveParametersAndInput()


def _post_process_segmentation(bias_corrected_image_names,
                               upsampled_image_names,
                               tissue_settings,
                               csf_factor,
                               parameters_and_inputs,
                               transformed_template_name,
                               affine_atlas,
                               before_morpho_name,
                               upper_mask):

    logger.info('Upsampling bias corrected images.')
    for input_number, bias_corrected in enumerate(bias_corrected_image_names):
        corrected_input = nib.load(bias_corrected)
        resampled_input, new_affine, orig_res = resample_vol(
                                                    corrected_input.get_data(),
                                                    corrected_input.affine,
                                                    0.5, order=1)
        upsampled = nib.Nifti1Image(resampled_input, new_affine)
        nib.save(upsampled, upsampled_image_names[input_number])

    # Next we need to reconstruct the segmentation with the upsampled data
    # and map it into the simnibs tissues
    upsampled_tissues, upper_part = samseg.simnibs_segmentation_utils.segmentUpsampled(
                                    upsampled_image_names,
                                    tissue_settings,
                                    parameters_and_inputs,
                                    transformed_template_name,
                                    affine_atlas,
                                    csf_factor)

    del parameters_and_inputs
    #Cast the upsampled image to int16  to save space
    for upsampled_image in upsampled_image_names:
        upsampled = nib.load(upsampled_image)
        upsampled.set_data_dtype(np.int16)
        nib.save(upsampled, upsampled_image)

    affine_upsampled = upsampled.affine
    upsampled_tissues_im = nib.Nifti1Image(upsampled_tissues,
                                        affine_upsampled)
    upsampled_tissues_im.set_data_dtype(np.int16)
    nib.save(upsampled_tissues_im, before_morpho_name)

    upper_part_im = nib.Nifti1Image(upper_part.astype(np.int16),affine_upsampled)
    upper_part_im.set_data_dtype(np.int16)
    nib.save(upper_part_im, upper_mask)
    del upper_part_im
    # Do morphological operations
    simnibs_tissues = tissue_settings['simnibs_tissues']
    _morphological_operations(upsampled_tissues, upper_part, simnibs_tissues)

    return upsampled_tissues


def _morphological_operations(label_img, upper_part, simnibs_tissues):
    ''' Does morphological operations to
        1. Smooth out the labeling and remove noise
        2. A CSF layer between GM and Skull and between GM and CSF
        3. Outer bone layers are compact bone
    '''
    se = ndimage.generate_binary_structure(3,3)
    se_n = ndimage.generate_binary_structure(3,1)

    # Do some clean-ups, mainly CSF and skull
    # First combine the WM and GM
    brain = (label_img == simnibs_tissues['WM']) | \
            (label_img == simnibs_tissues['GM'])
    dil = mrph.binary_dilation(_get_largest_components(
        mrph.binary_erosion(brain, se_n, 1), se_n, vol_limit=10), se_n, 1)

    unass = brain ^ dil

    #Add the CSF and open again
    csf = (label_img == simnibs_tissues['CSF'])
    brain_csf = brain | csf
    # Vol limit in voxels
    dil = mrph.binary_dilation(_get_largest_components(
        mrph.binary_erosion(brain_csf, se, 1), se_n, vol_limit=80), se, 1)
    unass |= (csf & ~dil)
    del brain, csf, dil

    #Clean the outer border of the skull
    bone = (label_img == simnibs_tissues['Compact_bone']) | \
           (label_img == simnibs_tissues['Spongy_bone'])
    veins = (label_img == simnibs_tissues['Blood'])
    air_pockets = (label_img == simnibs_tissues['Air_pockets'])
    scalp = (label_img == simnibs_tissues['Scalp'])
    muscle = (label_img == simnibs_tissues['Muscle'])
    eyes = (label_img == simnibs_tissues['Eyes'])
    #Use scalp to clean out noisy skull bits within the scalp
    skull_outer = brain_csf | bone | veins | air_pockets
    skull_outer = mrph.binary_fill_holes(skull_outer, se_n)
    num_iter = 2
    skull_outer = mrph.binary_dilation(_get_largest_components(
        mrph.binary_erosion(skull_outer, se, num_iter), se_n, vol_limit=50), se, num_iter)
    skull_inner = bone | scalp | air_pockets | muscle | eyes
    skull_inner = mrph.binary_fill_holes(skull_inner, se_n)
    skull_inner = mrph.binary_dilation(_get_largest_components(
        mrph.binary_erosion(skull_inner, se, num_iter), se_n, vol_limit=50), se, num_iter)
    dil = bone & skull_outer & skull_inner
    #Protect thin areas that would be removed by erosion
    bone_thickness = _calc_thickness(dil)
    thin_parts = bone & (bone_thickness < 3.5) & (bone_thickness > 0)
    dil |= thin_parts
    del thin_parts
    del bone_thickness

    unass = unass | (dil ^ bone)
    del bone, dil

    # Open the veins
    dil = mrph.binary_dilation(_get_largest_components(
        mrph.binary_erosion(veins, se, 1), se_n, vol_limit=10), se, 1)
    unass |= (dil ^ veins)

    # Clean the eyes
    dil = mrph.binary_dilation(_get_largest_components(
        mrph.binary_erosion(eyes, se, 1), se_n, vol_limit=10), se, 1)
    #dil = mrph.binary_opening(eyes, se, 1)
    unass |= (dil ^ eyes)

    #Clean muscles
    dil = mrph.binary_opening(muscle, se, 1)
    unass |= (dil ^ muscle)

    # And finally the scalp
    head = scalp | skull_outer | eyes | muscle
    dil = mrph.binary_dilation(_get_largest_components(
                  mrph.binary_erosion(head, se, 2), se_n, num_limit=1), se, 2)
    unass |= (scalp & ~dil)

    # Filling missing parts
    # NOTE: the labeling is uint16, so I'll code the unassigned voxels with the max value
    # which is 65535
    label_img[unass] = 65535
    #Add background to tissues
    simnibs_tissues["BG"] = 0
    _smoothfill(label_img, unass, simnibs_tissues)


    # Ensuring a CSF layer between GM and Skull and GM and Blood
    # Relabel regions in the expanded GM which are in skull or blood to CSF
    print('Ensure CSF')

    brain_gm = (label_img == simnibs_tissues['GM'])
    C_BONE = (label_img == simnibs_tissues['Compact_bone'])
    S_BONE = (label_img == simnibs_tissues['Spongy_bone'])

    brain_dilated = mrph.binary_dilation(brain_gm, se, 1)
    overlap = brain_dilated & (C_BONE | S_BONE)
    label_img[overlap] = simnibs_tissues['CSF']

    CSF_brain = brain_gm | (label_img == simnibs_tissues['CSF'])
    CSF_brain_dilated = mrph.binary_dilation(CSF_brain, se, 1)
    spongy_csf = CSF_brain_dilated & S_BONE
    upper_part = mrph.binary_erosion(upper_part,se,6)
    skin_csf = CSF_brain_dilated & (label_img == simnibs_tissues['Scalp']) & upper_part
    label_img[spongy_csf] = simnibs_tissues['Compact_bone']
    label_img[skin_csf] = simnibs_tissues['Compact_bone']

    # Ensure the outer skull label is compact bone
    C_BONE = (label_img == simnibs_tissues['Compact_bone'])
    S_BONE = (label_img == simnibs_tissues['Spongy_bone'])
    SKULL_outer = (C_BONE | S_BONE) & ~mrph.binary_erosion((C_BONE | S_BONE), se, 1)
    label_img[SKULL_outer] = simnibs_tissues['Compact_bone']
    #Relabel air pockets to air
    label_img[label_img == simnibs_tissues['Air_pockets']] = 0

def _smoothfill(label_img, unassign, simnibs_tissues):
    """Hackish way to fill unassigned voxels,
       works by smoothing the masks and binarizing
       the smoothed masks.
    """
    sum_of_unassigned = np.inf
    # for as long as the number of unassigned is changing
    # let's find the binarized version using a running
    # max index as we need to loop anyway. This way we
    # don't need to call np.argmax
    while unassign.sum() < sum_of_unassigned:
        sum_of_unassigned = unassign.sum()
        if(sum_of_unassigned == 0):
            break
        labs = 65535*np.ones_like(label_img)
        max_val = np.zeros_like(label_img, dtype=np.float32)
        for i, t in enumerate(simnibs_tissues):
            #Don't smooth WM
            vol = (label_img == simnibs_tissues[t]).astype(np.float32)
            if t == 'WM':
                cs = vol
            else:
                cs = gaussian_filter(vol, 1)

            #Check the max values and update
            max_mask = cs > max_val
            labs[max_mask] = simnibs_tissues[t]
            max_val[max_mask] = cs[max_mask]

        label_img[:] = labs[:]
        unassign = (labs == 65535)
        del labs




def _fill_missing(label_img, unassign):
    """Hackish way to fill unassigned voxels,
       works by smoothing the masks and binarizing
       the smoothed masks. Works in-place
    """

    def _get_most_frequent_label(x, y, z, label_im, pad):
        nhood = label_im[max(0,x-pad):min(label_im.shape[0],x+pad),
                         max(0,y-pad):min(label_im.shape[1],y+pad),
                         max(0,z-pad):min(label_im.shape[2],z+pad)]

        dim = nhood.shape

        return np.bincount(nhood.reshape((np.prod(dim)))).argmax()


    label_img[unassign] = 255
    inds_tmp = unassign.nonzero()
    inds_tmp = np.array(inds_tmp)
    num_unassigned = (label_img==255).sum()
    num_unassigned_new = num_unassigned
    while num_unassigned_new > 0:
        fill_map = map(partial(_get_most_frequent_label, label_im=label_img, pad=5),
                       np.nditer(inds_tmp[0,:]), np.nditer(inds_tmp[1,:]), np.nditer(inds_tmp[2,:]))
        fill_array = np.array(list(fill_map))
        label_img[inds_tmp[0,:], inds_tmp[1,:], inds_tmp[2,:]] = fill_array
        inds_tmp = inds_tmp[:,fill_array==255]
        num_unassigned_new = (label_img==255).sum()
        print('Unassigned: ' + str(num_unassigned_new))
        if num_unassigned_new == num_unassigned:
            print('Number of unassigned voxels not going down. Breaking.')
            break

        num_unassigned = num_unassigned_new

def _get_largest_components(vol, se, vol_limit=0, num_limit=-1, return_sizes=False):
    """Get the n largest components from a volume.

    PARAMETERS
    ----------
    vol : ndarray
        Image volume. A dimX x dimY x dimZ array containing the image data.
    se : ndarray
        Structuring element to use when detecting components (i.e. setting the
        connectivity which defines a component).
    n : int
        Number of (largest) components to retain.
    return_sizes : bool, optional
        Whether or not to also return the sizes (in voxels) of each component
        that was retained (default = False).

    RETURNS
    ----------
    Components : ndarray
        Binary dimX x dimY x dimZ array where entries corresponding to retained
        components are True and the remaining entries are False.
    """
    vol_lbl = label(vol, se)[0]
    labels, region_size = np.unique(vol_lbl, return_counts=True)
    labels = labels[1:]  # disregard background (label=0)
    region_size = region_size[1:]  #
    mask = region_size > vol_limit
    region_size = region_size[mask]
    labels = labels[mask]

    if num_limit == -1:
        num_limit = len(labels)

    labels = labels[np.argsort(region_size)[::-1]]
    components = np.zeros_like(vol)
    for i in labels[:num_limit]:
        components = components | (vol_lbl == i)

    if return_sizes:
        return components, region_size[labels[:num_limit]]
    else:
        return components



def _registerT1T2(fixed_image, moving_image, output_image):
    registerer = samseg.gems.KvlImageRegisterer()
    registerer.read_images(fixed_image, moving_image)
    registerer.initialize_transform()
    registerer.register()
    registerer.write_out_result(output_image)

    #The register function uses double internally
    #Let's cast to float32 and copy the header from
    # the fixed image to avoid any weird behaviour
    if os.path.exists(output_image):
        T2_reg = nib.load(output_image)
        fixed_tmp = nib.load(fixed_image)
        T2_data = T2_reg.get_data().astype(np.float32)
        T2_im = nib.Nifti1Image(T2_data,fixed_tmp.affine)
        nib.save(T2_im, output_image)


def _fillin_gm_layer(label_img, label_affine, labelorg_img, labelorg_affine, m,
                    exclusion_tissues = {"left_cerebral_wm": 2,
                                         "right_cerebral_wm": 41,
                                         "stuff_to_exclude": [4, 10, 14, 16, 24, 28, 43, 49, 60]},
                    relabel_tissues = {"GM": 2, "stuff_to_relabel": [1, 3]}):
    ''' relabels WM and CSF that intersect with the central GM surface to GM
        an exclusion mask is used to prevent relabelign in central brain regions
    '''
    # generate exclusion mask: estimate corpus callossum 
    exclude_img = mrph.binary_dilation(labelorg_img == exclusion_tissues['left_cerebral_wm'],iterations=2) 
    exclude_img *= mrph.binary_dilation(labelorg_img == exclusion_tissues['right_cerebral_wm'],iterations=2)
    # add other tissues 
    for i in exclusion_tissues['stuff_to_exclude']:
        exclude_img += labelorg_img == i
    exclude_img = mrph.binary_dilation(exclude_img,iterations=8)
    # upsample exclude_img
    iM = np.linalg.inv(labelorg_affine).dot(label_affine)
    exclude_img = affine_transform(exclude_img, iM[:3, :3], iM[:3, 3], label_img.shape, order=0)
    
    # generate voxel mask of middle GM
    mask = mask_from_surface(m.nodes[:],m.elm[:,:3]-1,label_affine,label_img.shape)
    mask = mrph.binary_dilation(mask,iterations=1) * ~mrph.binary_erosion(mask,iterations=1)
    mask[exclude_img] = 0
    
    # relabel WM and CSF parts to GM
    for i in relabel_tissues['stuff_to_relabel']:
        label_img[ (label_img == i)*mask ] = relabel_tissues['GM']
        
    return label_img


def _open_sulci(label_img, label_affine, labelorg_img, labelorg_affine, m, 
                tissue_labels = {"CSF": 3, "GM": 2, "WM": 1},
                exclusion_tissues = [17, 18, 53, 54]):
    # get thin CSF structures
    mask = mask_from_surface(m.nodes[:],m.elm[:,:3]-1,label_affine,label_img.shape)
    # mask2 = mrph.binary_dilation(mask,iterations=4)
    # mask2 = mrph.binary_erosion(mask2,iterations=5)
    mask2 = mrph.binary_dilation(mask,iterations=2)
    mask2 = mrph.binary_erosion(mask2,iterations=3)
    mask2[mask] = 0
    
    # protect hippocampi and amydalae
    exclude_img = np.zeros_like(labelorg_img,dtype=bool)
    for i in exclusion_tissues:
        exclude_img += labelorg_img == i
    iM = np.linalg.inv(labelorg_affine).dot(label_affine)
    exclude_img = affine_transform(exclude_img, iM[:3, :3], iM[:3, 3], label_img.shape, order=0)
    mask2[exclude_img] = 0
    
    # relabel GM overlapping thin CSF to CSF
    label_img[ (label_img == tissue_labels['GM'])*mask2 ] = tissue_labels['CSF']
    
    # open up remaining thin GM bridges at brain surface
    mask2 = (label_img == tissue_labels['GM']) | (label_img == tissue_labels['WM']) 
    # mask2 = mrph.binary_erosion(mask2,iterations=2)
    # mask2 = mrph.binary_dilation(mask2,iterations=2)
    #label_img[ (label_img == tissue_labels['GM'])* ~mask2 ] = tissue_labels['CSF']
    brainthickness = _calc_thickness(mask2)     
    label_img[(brainthickness <= 2.0) * (label_img == tissue_labels['GM'])] = tissue_labels['CSF']
    
    return label_img


def view(subject_dir):
    print('charm viewer not yet implemented, sorry...')


def run(subject_dir=None, T1=None, T2=None,
        registerT2=False, initatlas=False, segment=False,
        create_surfaces=False, mesh_image=False, usesettings=None,
        noneck=False, options_str=None):
    """charm pipeline

    PARAMETERS
    ----------
    subject_dir : str, mandatory
        output directory
    T1 : str
        filename of T1 image
    T2 : str
        filename of T2 image

    --> parameters to control the workflow:
    registerT2 : bool
        run T2-to-T1 registration (default = False)
    initatlas : bool
        run affine registration of atlas to input images (default = False)
    segment : bool
        run volume and surface segmentation (default = False)
    mesh_image : bool
        run tetrahedral meshing (default = False)
    --> further parameters:
    usesettings : str
        filename of alternative settings-file (default = None)
    options_str : str
        string of command line options to add to logging (default = None)
    RETURNS
    ----------
        None
    """
    # ------------------------START UP-----------------------------------------
    start = time.time()

    # make subject_dir if not existent
    if not os.path.exists(subject_dir):
        os.mkdir(subject_dir)

    # initialize subject files
    sub_files = file_finder.SubjectFiles(None, subject_dir)
    # start logging ...
    os.makedirs(sub_files.report_folder, exist_ok=True)
    logfile = os.path.join(sub_files.report_folder, "charm_log.html")
    with open(logfile, 'a') as f:
        f.write('<HTML><HEAD><TITLE>charm report</TITLE></HEAD><BODY><pre>')
        f.close()
    fh = logging.FileHandler(logfile, mode='a')
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    fh.setFormatter(formatter)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    utils.simnibs_logger.register_excepthook(logger)

    logger.info('simnibs version '+__version__)
    logger.info('charm run started: '+time.asctime())
    if options_str is not None:    
        logger.debug('options: '+options_str)
    else:
        logger.debug('options: none')

    # copy T1 (as nii.gz) if supplied
    if T1 is not None:
        if not os.path.exists(T1):
            raise FileNotFoundError(f'Could not find input T1 file: {T1}')
        else:
            #Cast to float32 and save
            T1_tmp = nib.load(T1)
            T1_tmp.set_data_dtype(np.float32)
            nib.save(T1_tmp, sub_files.reference_volume)

    if registerT2 and T2 is not None:
        if not os.path.exists(T2):
            raise FileNotFoundError(f'Could not find input T2 file: {T2}')
        else:
            _registerT1T2(T1, T2, sub_files.T2_reg)
    elif T2 is not None:
        if os.path.exists(T2):
            T2_tmp = nib.load(T2)
            T2_tmp.set_data_dtype(np.float32)
            nib.save(T2_tmp, sub_files.T2_reg)

    # Create visualization html
    if T2 is not None:
        logger.info('Creating registration visualization')
        plotting.viewer_registration(sub_files.reference_volume,
                                     sub_files.T2_reg,
                                     sub_files.reg_viewer)

    # read settings and copy settings file
    if usesettings is None:
        fn_settings = os.path.join(SIMNIBSDIR, 'charm.ini')
    else:
        if type(usesettings) == list:
            usesettings = usesettings[0]
        fn_settings = usesettings
    settings = utils.settings_reader.read_ini(fn_settings)
    try:
        shutil.copyfile(fn_settings, sub_files.settings)
    except shutil.SameFileError:
        pass
    logger.debug(settings)
    
    # -------------------------PIPELINE STEPS---------------------------------
    # TODO: denoise T1 here with the sanlm filter, T2 denoised after coreg.
    # Could be before as well but doesn't probably matter that much.
    # Denoising after has the benefit that we keep around the non-denoised
    # versions of all inputs.

    # Make the output path for segmentations
    os.makedirs(sub_files.segmentation_folder, exist_ok=True)

    denoise_settings = settings['preprocess']
    if denoise_settings['denoise'] and T1 is not None:
        logger.info('Denoising the T1 input and saving.')
        _denoise_input_and_save(T1, sub_files.T1_denoised)

    # denoise if needed
    if denoise_settings['denoise'] and os.path.exists(sub_files.T2_reg) and T2 is not None:
        logger.info('Denoising the registered T2 and saving.')
        _denoise_input_and_save(sub_files.T2_reg,
                                sub_files.T2_reg_denoised)

    # Set-up samseg related things before calling the affine registration
    # and/or segmentation
    samseg_settings = settings['samseg']

    # Specify the maximum number of threads the GEMS code will use
    # by default this is all, but can be changed in the .ini
    num_threads = samseg_settings['threads']
    if isinstance(num_threads, int) and num_threads > 0:
        samseg.setGlobalDefaultNumberOfThreads(num_threads)
        logger.info('Using %d threads, instead of all available.'
                    % num_threads)

    # TODO: Setup the visualization tool. This needs some pyqt stuff to be
    # installed. Don't know if we want to expose this in the .ini
    showFigs = False
    showMovies = False
    visualizer = samseg.initVisualizer(showFigs, showMovies)

    # Set-up atlas paths
    atlas_name = samseg_settings['atlas_name']
    logger.info('Using '+atlas_name+' as charm atlas.')
    atlas_path = os.path.join(file_finder.templates.charm_atlas_path,
                              atlas_name)
    atlas_settings = utils.settings_reader.read_ini(
            os.path.join(atlas_path, atlas_name+'.ini'))
    atlas_settings_names = atlas_settings['names']
    template_name = os.path.join(atlas_path,
                                 atlas_settings_names['template_name'])
    atlas_affine_name = os.path.join(atlas_path,
                                     atlas_settings_names['affine_atlas'])
    atlas_level1 = os.path.join(atlas_path,
                                atlas_settings_names['atlas_level1'])

    atlas_level2 = os.path.join(atlas_path,
                                    atlas_settings_names['atlas_level2'])
    if os.path.exists(sub_files.T2_reg):
        gmm_parameters = os.path.join(atlas_path,
                                      atlas_settings_names['gaussian_parameters_t2'])
    else:
        gmm_parameters = os.path.join(atlas_path,
                                      atlas_settings_names['gaussian_parameters_t1'])

    neck_tissues = atlas_settings['neck_optimization']
    neck_tissues = neck_tissues['neck_tissues']

    if initatlas:
        # initial affine registration of atlas to input images,
        # including break neck
        logger.info('Starting affine registration and neck correction.')
        init_atlas_settings = settings['initatlas']
        if denoise_settings['denoise']:
            inputT1 = sub_files.T1_denoised
        else:
            inputT1 = T1

        _register_atlas_to_input_affine(inputT1, template_name,
                                        atlas_affine_name,
                                        atlas_level1,
                                        atlas_level2,
                                        sub_files.segmentation_folder,
                                        sub_files.template_coregistered,
                                        init_atlas_settings,
                                        neck_tissues,
                                        visualizer,
                                        noneck)

        logger.info('Creating affine registration visualization')
        plotting.viewer_affine(sub_files.reference_volume,
                               sub_files.template_coregistered,
                               sub_files.affine_reg_viewer)

    if segment:
        # This part runs the segmentation, upsamples bias corrected output,
        # writes mni transforms, creates upsampled segmentation, maps tissues
        # to conductivities, runs morphological operations

        # Run the segmentation and return the class, which is needed
        # for further post-processing
        # The bias field kernel size has to be changed based on input
        input_images = []
        if denoise_settings['denoise']:
            input_images.append(sub_files.T1_denoised)
        else:
            input_images.append(T1)

        if os.path.exists(sub_files.T2_reg):
            if denoise_settings['denoise']:
                input_images.append(sub_files.T2_reg_denoised)
            else:
                input_images.append(sub_files.T2_reg)

        segment_settings = settings['segment']
        logger.info('Estimating parameters.')
        segment_parameters_and_inputs = _estimate_parameters(
                sub_files.segmentation_folder,
                sub_files.template_coregistered,
                atlas_path, input_images,
                segment_settings,
                gmm_parameters,
                visualizer)

        # Okay now the parameters have been estimated, and we can segment the
        # scan. However, we need to also do this at an upsampled resolution,
        # so first write out the bias corrected scan, and the segmentation.

        bias_corrected_image_names = [sub_files.T1_bias_corrected]
        if len(input_images) > 1:
            bias_corrected_image_names.append(sub_files.T2_bias_corrected)

        logger.info('Writing out normalized images and labelings.')
        os.makedirs(sub_files.surface_folder, exist_ok=True)
        cat_images = [sub_files.norm_image,
                      sub_files.cereb_mask,
                      sub_files.subcortical_mask,
                      sub_files.hemi_mask]

        cat_structs = atlas_settings['CAT_structures']
        tissue_settings = atlas_settings['conductivity_mapping']
        csf_factor = segment_settings['csf_factor']
        samseg.simnibs_segmentation_utils.writeBiasCorrectedImagesAndSegmentation(
                        bias_corrected_image_names,
                        sub_files.labeling,
                        segment_parameters_and_inputs,
                        tissue_settings,
                        csf_factor,
                        cat_structure_options=cat_structs,
                        cat_images=cat_images)

        fn_LUT=sub_files.labeling.rsplit('.',2)[0]+'_LUT.txt'
        shutil.copyfile(file_finder.templates.labeling_LUT, fn_LUT)
        
        # Write out MNI warps
        logger.info('Writing out MNI warps.')
        os.makedirs(sub_files.mni_transf_folder, exist_ok=True)
        samseg.simnibs_segmentation_utils.saveWarpField(
                template_name,
                sub_files.mni2conf_nonl,
                sub_files.conf2mni_nonl,
                segment_parameters_and_inputs)

        # Run post-processing
        logger.info('Post-processing segmentation')
        os.makedirs(sub_files.label_prep_folder, exist_ok=True)

        upsampled_image_names = [sub_files.T1_upsampled]
        if len(bias_corrected_image_names) > 1:
            upsampled_image_names.append(sub_files.T2_upsampled)

        cleaned_upsampled_tissues = _post_process_segmentation(
                                    bias_corrected_image_names,
                                    upsampled_image_names,
                                    tissue_settings,
                                    csf_factor,
                                    segment_parameters_and_inputs,
                                    sub_files.template_coregistered,
                                    atlas_affine_name,
                                    sub_files.tissue_labeling_before_morpho,
                                    sub_files.upper_mask)

        # Write to disk
        upsampled_image = nib.load(sub_files.T1_upsampled)
        affine_upsampled = upsampled_image.affine
        upsampled_tissues = nib.Nifti1Image(cleaned_upsampled_tissues,
                                            affine_upsampled)
        nib.save(upsampled_tissues, sub_files.tissue_labeling_upsampled)
        del cleaned_upsampled_tissues

        fn_LUT=sub_files.tissue_labeling_upsampled.rsplit('.',2)[0]+'_LUT.txt'
        shutil.copyfile(file_finder.templates.final_tissues_LUT, fn_LUT)
        
    if create_surfaces:
        # Create surfaces ala CAT12
        logger.info('Starting surface creation')
        starttime = time.time()
        fsavgDir = file_finder.Templates().freesurfer_templates
        
        surface_settings = settings['surfaces']
        nprocesses = surface_settings['processes']
        surf = surface_settings['surf']
        pial = surface_settings['pial']
        vdist = surface_settings['vdist']
        voxsize_pbt = surface_settings['voxsize_pbt']
        voxsize_refineCS = surface_settings['voxsize_refinecs']
        th_initial = surface_settings['th_initial']
        no_selfintersections = surface_settings['no_selfintersections']
        fillin_gm_from_surf = surface_settings['fillin_gm_from_surf']
        open_sulci_from_surf = surface_settings['open_sulci_from_surf']
        exclusion_tissues_fillinGM = surface_settings['exclusion_tissues_fillin_gm']
        exclusion_tissues_openCSF = surface_settings['exclusion_tissues_open_csf']

        multithreading_script = [os.path.join(SIMNIBSDIR, 'segmentation', 'run_cat_multiprocessing.py')]
        argslist = ['--Ymf', sub_files.norm_image,
                    '--Yleft_path', sub_files.hemi_mask,
                    '--Ymaskhemis_path', sub_files.cereb_mask,
                    '--surface_folder', sub_files.surface_folder,
                    '--fsavgdir', fsavgDir,
                    '--surf'] + surf + [ 
                    '--pial'] + pial + [
                    '--vdist', str(vdist[0]), str(vdist[1]),
                    '--voxsize_pbt', str(voxsize_pbt[0]), str(voxsize_pbt[1]),
                    '--voxsizeCS', str(voxsize_refineCS[0]), str(voxsize_refineCS[1]),
                    '--th_initial', str(th_initial),
                    '--no_intersect', str(no_selfintersections),
                    '--nprocesses', str(nprocesses)]

        proc = subprocess.run([sys.executable] +
                                  multithreading_script + argslist,
                                  stderr=subprocess.PIPE) # stderr: standard stream for simnibs logger
        logger.debug(proc.stderr.decode('ASCII', errors='ignore').replace('\r', ''))
        proc.check_returncode()


        # print time duration
        elapsed = time.time() - starttime
        logger.info('Total time surface creation (HH:MM:SS):')
        logger.info(time.strftime('%H:%M:%S', time.gmtime(elapsed)))
        sub_files = file_finder.SubjectFiles(None, subject_dir)
        if fillin_gm_from_surf or open_sulci_from_surf:
            logger.info('Improving GM from surfaces')
            starttime = time.time()
            # original tissue mask used for mesh generation
            shutil.copyfile(sub_files.tissue_labeling_upsampled,
                            os.path.join(sub_files.label_prep_folder, 'before_surfmorpho.nii.gz'))
            label_nii = nib.load(sub_files.tissue_labeling_upsampled)
            label_img = np.asanyarray(label_nii.dataobj)
            label_affine = label_nii.affine

            
            # orginal label mask
            label_nii = nib.load(sub_files.labeling)
            labelorg_img = np.asanyarray(label_nii.dataobj)
            labelorg_affine = label_nii.affine

            if fillin_gm_from_surf:
                # GM central surfaces
                m = Msh()
                if 'lh' in surf:
                    m = m.join_mesh(read_gifti_surface(sub_files.get_surface('lh')))
                if 'rh' in surf:
                    m = m.join_mesh(read_gifti_surface(sub_files.get_surface('rh')))
                    # fill in GM and save updated mask
                if m.nodes.nr > 0:
                    label_img = _fillin_gm_layer(label_img, label_affine, 
                                                 labelorg_img, labelorg_affine, m, 
                                                 exclusion_tissues = exclusion_tissues_fillinGM)

                    label_nii = nib.Nifti1Image(label_img, label_affine)
                    nib.save(label_nii, sub_files.tissue_labeling_upsampled)
                else:
                    logger.warning("Neither lh nor rh reconstructed. Filling in from GM surface skipped")

            if open_sulci_from_surf:
                # GM pial surfaces
                m = Msh()
                if 'lh' in pial:
                    m2 = read_gifti_surface(sub_files.get_surface('lh', surf_type='pial'))
                    # remove self-intersections using meshfix
                    with tempfile.NamedTemporaryFile(suffix='.off') as f:
                        mesh_fn = f.name
                    write_off(m2, mesh_fn)
                    cmd = [file_finder.path2bin("meshfix"), mesh_fn, '-o', mesh_fn]
                    spawn_process(cmd, lvl=logging.DEBUG)
                    m = m.join_mesh(read_off(mesh_fn))
                    if os.path.isfile(mesh_fn): os.remove(mesh_fn)
                if 'rh' in pial:
                    m2 = read_gifti_surface(sub_files.get_surface('rh', surf_type='pial'))
                    # remove self-intersections using meshfix
                    with tempfile.NamedTemporaryFile(suffix='.off') as f:
                        mesh_fn = f.name
                    write_off(m2, mesh_fn)
                    cmd = [file_finder.path2bin("meshfix"), mesh_fn, '-o', mesh_fn]
                    spawn_process(cmd, lvl=logging.DEBUG)
                    m = m.join_mesh(read_off(mesh_fn))
                    if os.path.isfile(mesh_fn): os.remove(mesh_fn)
                if m.nodes.nr > 0:
                    _open_sulci(label_img, label_affine,
                                labelorg_img, labelorg_affine, m,
                                exclusion_tissues = exclusion_tissues_openCSF)                    
                    label_nii = nib.Nifti1Image(label_img, label_affine)
                    nib.save(label_nii, sub_files.tissue_labeling_upsampled)
                else:
                    logger.warning("Neither lh nor rh pial reconstructed. Opening up of sulci skipped.")
                
            # print time duration
            elapsed = time.time() - starttime
            logger.info('Total time cost for GM imrpovements in (HH:MM:SS):')
            logger.info(time.strftime('%H:%M:%S', time.gmtime(elapsed)))

    if mesh_image:
        # create mesh from label image
        logger.info('Starting mesh')
        label_image = nib.load(sub_files.tissue_labeling_upsampled)
        label_buffer = label_image.get_fdata().astype(np.uint16) # Cast to uint16, otherwise meshing complains
        label_affine = label_image.affine
        label_buffer, label_affine, _ = crop_vol(label_buffer, label_affine, 
                                                 label_buffer>0, thickness_boundary=5) 
                                                 # reduce memory consumption a bit
                                                 
        # Read in settings for meshing
        mesh_settings = settings['mesh']
        elem_sizes = mesh_settings['elem_sizes']
        smooth_size_field = mesh_settings['smooth_size_field']
        skin_facet_size = mesh_settings['skin_facet_size']
        facet_distances = mesh_settings['facet_distances']
        optimize = mesh_settings['optimize']
        remove_spikes = mesh_settings['remove_spikes']
        skin_tag = mesh_settings['skin_tag']
        if not skin_tag:
            skin_tag = None
        remove_twins = mesh_settings['remove_twins']
        hierarchy = mesh_settings['hierarchy']
        if not hierarchy:
            hierarchy = None
        smooth_steps = mesh_settings['smooth_steps']
        
        # Meshing
        final_mesh = create_mesh(label_buffer, label_affine,
                elem_sizes=elem_sizes,
                smooth_size_field=smooth_size_field,
                skin_facet_size=skin_facet_size, 
                facet_distances=facet_distances,
                optimize=optimize, 
                remove_spikes=remove_spikes, 
                skin_tag=skin_tag,
                remove_twins=remove_twins, 
                hierarchy=hierarchy,
                smooth_steps=smooth_steps)
        
        logger.info('Writing mesh')
        write_msh(final_mesh, sub_files.fnamehead)
        v = final_mesh.view(cond_list = cond.standard_cond())
        v.write_opt(sub_files.fnamehead)
        
        logger.info('Transforming EEG positions')
        idx = (final_mesh.elm.elm_type == 2)&(final_mesh.elm.tag1 == skin_tag)
        mesh = final_mesh.crop_mesh(elements = final_mesh.elm.elm_number[idx])
        
        if not os.path.exists(sub_files.eeg_cap_folder):
            os.mkdir(sub_files.eeg_cap_folder)
                    
        cap_files = glob.glob(os.path.join(file_finder.ElectrodeCaps_MNI, '*.csv'))
        for fn in cap_files:
            fn_out = os.path.splitext(os.path.basename(fn))[0]
            fn_out = os.path.join(sub_files.eeg_cap_folder, fn_out)
            transformations.warp_coordinates(
                    fn, sub_files.subpath,
                    transformation_direction='mni2subject',
                    out_name=fn_out+'.csv',
                    out_geo=fn_out+'.geo',
                    mesh_in = mesh)
        
        logger.info('Write label image from mesh')
        MNI_template = file_finder.Templates().mni_volume
        mesh = final_mesh.crop_mesh(elm_type=4)
        field = mesh.elm.tag1.astype(np.uint16)
        ed = ElementData(field)
        ed.mesh = mesh
        ed.to_deformed_grid(sub_files.mni2conf_nonl, MNI_template, 
                            out=sub_files.final_labels_MNI,
                            out_original=sub_files.final_labels,
                            method='assign',
                            reference_original=sub_files.reference_volume)
        
        fn_LUT=sub_files.final_labels.rsplit('.',2)[0]+'_LUT.txt'
        shutil.copyfile(file_finder.templates.final_tissues_LUT, fn_LUT)

    # -------------------------TIDY UP-----------------------------------------
    # Create final seg html viewer
    # Create visualization htmls
    logger.info('Creating registration visualization')
    plotting.viewer_final(sub_files.reference_volume,
                          sub_files.final_labels,
                          sub_files.final_viewer)
    # log stopping time and total duration ...
    logger.info('charm run finished: '+time.asctime())
    logger.info('Total running time: '+utils.simnibs_logger.format_time(
                                                            time.time()-start))

    # stop logging ...
    while logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])
    utils.simnibs_logger.unregister_excepthook()
    logging.shutdown()
    with open(logfile, 'a') as f:
        f.write('</pre></BODY></HTML>')
        f.close()

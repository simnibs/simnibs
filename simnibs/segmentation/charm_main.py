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
import numpy as np
from scipy import ndimage
import scipy.ndimage.morphology as mrph
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.measurements import label
from simnibs import SIMNIBSDIR

from .. import __version__
from . import samseg
from ._cat_c_utils import sanlm
from .. import utils
from ..utils.simnibs_logger import logger
from ..utils import file_finder
from ..utils import transformations
from ..utils.transformations import resample_vol, crop_vol
from ..mesh_tools.meshing import create_mesh
from ..mesh_tools.mesh_io import write_msh, ElementData
from ..simulation import cond

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
        affine.adjust_neck(T1, template_coregistered_name, mesh_level1,
                           mesh_level2, neck_search_bounds, neck_tissues,
                           visualizer, downsampling_target=2.0)
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
                               parameters_and_inputs,
                               transformed_template_name):

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
    upsampled_tissues = samseg.simnibs_segmentation_utils.segmentUpsampled(
                        upsampled_image_names,
                        tissue_settings,
                        parameters_and_inputs,
                        transformed_template_name)

    #Cast the upsampled image to int16  to save space
    for upsampled_image in upsampled_image_names:
        upsampled = nib.load(upsampled_image)
        upsampled.set_data_dtype(np.int16)
        nib.save(upsampled, upsampled_image)
    # Do morphological operations
    simnibs_tissues = tissue_settings['simnibs_tissues']
    _morphological_operations(upsampled_tissues, simnibs_tissues)

    return upsampled_tissues


def _morphological_operations(label_img, simnibs_tissues):
    ''' Does morphological operations to
        1. Smooth out the labeling and remove noise
        2. A CSF layer between GM and Skull and between GM and CSF
        3. Outer bone layers are compact bone
    '''
    se = ndimage.generate_binary_structure(3,3)
    se_n = ndimage.generate_binary_structure(3,1)
    # Get list of masks
    tissue_masks = []
    tissue_names = []
    for tkey in simnibs_tissues:
        tmp_img = label_img == simnibs_tissues[tkey]
        tissue_masks.append(tmp_img)
        tissue_names.append(tkey)

    # Do some clean-ups, mainly CSF and skull
    # First combine the WM and GM
    brain = tissue_masks[tissue_names.index('WM')] | \
            tissue_masks[tissue_names.index('GM')]

    dil = mrph.binary_dilation(_get_n_largest_components(
        mrph.binary_erosion(brain, se, 1), se_n, 1), se, 1)

    #Remove unattached components
    unass = brain ^ dil
    tissue_masks[tissue_names.index('WM')] = \
    tissue_masks[tissue_names.index('WM')] & dil

    tissue_masks[tissue_names.index('GM')] = \
    tissue_masks[tissue_names.index('GM')] & dil

    #Add the CSF and open again
    brain_csf = brain + tissue_masks[tissue_names.index('CSF')]
    dil = mrph.binary_opening(brain_csf, se, 1)

    unass = unass + (brain_csf ^ dil)
    tissue_masks[tissue_names.index('CSF')] = \
        tissue_masks[tissue_names.index('CSF')] & dil

    #Close the skull (spongy + compact), remove noise and assign to compact
    bone = tissue_masks[tissue_names.index('Compact_bone')] | \
           tissue_masks[tissue_names.index('Spongy_bone')]

    dil = mrph.binary_fill_holes(bone,se)
    dil = mrph.binary_opening(dil,se, 1)
    dil = mrph.binary_dilation(_get_n_largest_components(
                  mrph.binary_erosion(dil, se, 1), se, 20), se, 1)

    unass = unass + (dil ^ bone)
    tissue_masks[tissue_names.index('Compact_bone')] = \
    tissue_masks[tissue_names.index('Compact_bone')] & dil

    tissue_masks[tissue_names.index('Spongy_bone')] = \
    tissue_masks[tissue_names.index('Spongy_bone')] & dil


    # Open the veins
    veins = tissue_masks[tissue_names.index('Blood')]
    #dil = mrph.binary_opening(veins, se, 1)
    dil = mrph.binary_dilation(_get_n_largest_components(
        mrph.binary_erosion(veins, se, 1), se_n, 10), se, 1)
    unass = unass + (dil ^ veins)
    tissue_masks[tissue_names.index('Blood')] = \
    tissue_masks[tissue_names.index('Blood')] & dil
    # And finally the scalp
    scalp = tissue_masks[tissue_names.index('Scalp')]
    dil = mrph.binary_dilation(_get_n_largest_components(
                  mrph.binary_erosion(scalp, se, 1), se_n, 1), se, 1)
    unass = unass + (dil ^ scalp)

    tissue_masks[tissue_names.index('Scalp')] = \
    tissue_masks[tissue_names.index('Scalp')] & dil

    # Temporarily add air here as a mask, otherwise smoothfill
    # happily dilates tissues next to air
    tissue_masks.append(label_img == 0)
    # Now run the smoothfill
    tissue_masks = _smoothfill(tissue_masks, unass)

    #Now we can get rid of the air again
    tissue_masks.pop()

    # Ensuring a CSF layer between GM and Skull and GM and Blood
    # Relabel regions in the expanded GM which are in skull or blood to CSF
    GM = tissue_masks[tissue_names.index('GM')]
    C_BONE = tissue_masks[tissue_names.index('Compact_bone')]
    S_BONE = tissue_masks[tissue_names.index('Spongy_bone')]
    BLOOD = tissue_masks[tissue_names.index('Blood')]
    GM_dilated = mrph.binary_dilation(GM, se, 1)
    GM_dilated = GM_dilated & (C_BONE | S_BONE | BLOOD)
    #CSF = tissue_masks[tissue_names.index('CSF')]
    #CSF[GM_dilated & (C_BONE | S_BONE | BLOOD)] = True
    tissue_masks[tissue_names.index('CSF')] = tissue_masks[tissue_names.index('CSF')] | GM_dilated
    tissue_masks[tissue_names.index('Compact_bone')] = tissue_masks[tissue_names.index('Compact_bone')] & \
                                                       ~GM_dilated
    tissue_masks[tissue_names.index('Spongy_bone')] = tissue_masks[tissue_names.index('Spongy_bone')] & \
                                                       ~GM_dilated
    tissue_masks[tissue_names.index('Blood')] = tissue_masks[tissue_names.index('Blood')] & \
                                                       ~GM_dilated

    # Ensure the outer skull label is compact bone
    SKULL_outer = (C_BONE | S_BONE) & ~mrph.binary_erosion((C_BONE | S_BONE), se)
    #C_BONE[SKULL_outer] = True
    tissue_masks[tissue_names.index('Compact_bone')] = tissue_masks[tissue_names.index('Compact_bone')] | SKULL_outer
    tissue_masks[tissue_names.index('Spongy_bone')] = tissue_masks[tissue_names.index('Spongy_bone')] & \
                                                       ~SKULL_outer
    # Finally dilate the Skin layer and relabel CSF to bone so we have skull everywhere
    scalp = tissue_masks[tissue_names.index('Scalp')]
    scalp = mrph.binary_dilation(scalp, se, 1)
    relabel = scalp & tissue_masks[tissue_names.index('CSF')]
    tissue_masks[tissue_names.index('Compact_bone')] = tissue_masks[tissue_names.index('Compact_bone')] | \
                                                       relabel
    tissue_masks[tissue_names.index('CSF')] = tissue_masks[tissue_names.index('CSF')] &\
                                              ~relabel

    #Relabel the tissues to produce the final label image
    for tname, tmask in zip(tissue_names, tissue_masks):
        label_img[tmask] = simnibs_tissues[tname]


def _smoothfill(vols, unassign):
    """Hackish way to fill unassigned voxels,
       works by smoothing the masks and binarizing
       the smoothed masks.
    """
    #Cast to uint8, note this is done in-place, and
    #map to 0->255
    vols = [255*v.astype(np.uint8, copy=False) for v in vols]
    unassign = unassign.copy()
    sum_of_unassigned = np.inf
    # for as long as the number of unassigned is changing
    # let's find the binarized version using a running
    # max index as we need to loop anyway. This way we
    # don't need to call np.argmax
    while unassign.sum() < sum_of_unassigned:
        sum_of_unassigned = unassign.sum()
        #Note here that unassigned voxels
        #in the end will have index 255
        #so if we at some point will have
        #255 different labels here this
        #will result in weird behavior
        inds = 255*np.ones_like(vols[0])
        max_val = np.zeros_like(vols[0],dtype=np.float32)
        for i, ivol in enumerate(vols):
            cs = gaussian_filter(ivol.astype(np.float32), 1)
            #cs = cs.astype(np.uint8, copy=False)
            cs[ivol == 255] = 255
            #Check the max values and update
            max_mask = cs > max_val
            inds[max_mask] = i
            max_val[max_mask] = cs[max_mask]


        vols = [inds==i for i in range(len(vols))]
        unassign = inds == 255

    return vols



def _get_n_largest_components(vol, se, n, return_sizes=False):
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
    labels = labels[np.argsort(region_size)[::-1]]
    components = np.any(np.array([vol_lbl == i for i in labels[:n]]),
                        axis=0)

    # if no components are found, components will be reduced to false. Replace
    # by array of appropriate size
    if components.sum() == 0:
        components = np.zeros_like(vol, dtype=bool)

    if return_sizes:
        return components, region_size[labels[:n]]
    else:
        return components

def _binarize(vols, return_empty=False):
    """Binarize a list of input volumes by finding the maximum posterior
    probability of each and assigning the voxel to this volume.

    PARAMETERS
    ----------
    vols : list
        List of filenames or nibabel image objects (describing probabilities of
        different tissue types).
    return_empty : bool
        If true, return an array containing all voxels which are not assigned
        to either of the other volumes (default: False)
    RETURNS
    ----------
    bin_vols : list
        List of ndarrays describing binarized versions of the input volumes.
    unassign : ndarray
        Array containing any unassigned voxels.
    """
    # Concatenate arrays/images
    imgs = np.concatenate(tuple([v[..., np.newaxis] for v in vols]), axis=3)
    imgs = np.concatenate((np.zeros_like(vols[0])[..., np.newaxis], imgs),
                          axis=3)

    # Find max indices
    max_idx = np.argmax(imgs, axis=3)

    # Binarize. Here vols_bin[0] contain voxels not assigned to any other
    # volume
    vols_bin = []
    for i in range(imgs.shape[-1]):
        vols_bin.append(max_idx == i)

    if return_empty:
        return vols_bin[1:] + [vols_bin[0]]
    else:
        return vols_bin[1:]

def _registerT1T2(fixed_image, moving_image, output_image):
    registerer = samseg.gems.KvlImageRegisterer()
    registerer.read_images(fixed_image, moving_image)
    registerer.initialize_transform()
    registerer.register()
    registerer.write_out_result(output_image)

    #The register function uses double internally
    #Let's cast to float32
    if os.path.exists(output_image):
        T2_reg = nib.load(output_image)
        T2_reg.set_data_dtype(np.float32)
        nib.save(T2_reg, output_image)


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

    # start logging ...
    logfile = os.path.join(subject_dir, "charm_log.html")
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
    
    # initialize subject files
    sub_files = file_finder.SubjectFiles(None, subject_dir)

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
        samseg.simnibs_segmentation_utils.writeBiasCorrectedImagesAndSegmentation(
                        bias_corrected_image_names,
                        sub_files.labeling,
                        segment_parameters_and_inputs,
                        cat_structure_options=cat_structs,
                        cat_images=cat_images)

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

        tissue_settings = atlas_settings['conductivity_mapping']
        cleaned_upsampled_tissues = _post_process_segmentation(
                                    bias_corrected_image_names,
                                    upsampled_image_names,
                                    tissue_settings,
                                    segment_parameters_and_inputs,
                                    sub_files.template_coregistered)

        # Write to disk
        upsampled_image = nib.load(sub_files.T1_upsampled)
        affine_upsampled = upsampled_image.affine
        upsampled_tissues = nib.Nifti1Image(cleaned_upsampled_tissues,
                                            affine_upsampled)
        nib.save(upsampled_tissues, sub_files.tissue_labeling_upsampled)

    if create_surfaces:
        # Create surfaces ala CAT12
        logger.info('Starting surface creation')
        starttime = time.time()
        fsavgDir = file_finder.Templates().freesurfer_templates
        
        surface_settings = settings['surfaces']
        nprocesses = surface_settings['processes']
        surf = surface_settings['surf']
        vdist = surface_settings['vdist']
        voxsize_pbt = surface_settings['voxsize_pbt']
        voxsize_refineCS = surface_settings['voxsize_refinecs']
        th_initial = surface_settings['th_initial']
        no_selfintersections = surface_settings['no_selfintersections']
        
        if sys.platform == 'win32':
            # A hack to get multithreading to work on Windows
            multithreading_script = [os.path.join(SIMNIBSDIR, 'segmentation', 'run_cat_multiprocessing.py')]
            argslist = ['--Ymf', sub_files.norm_image,
                    '--Yleft_path', sub_files.hemi_mask,
                    '--Ymaskhemis_path', sub_files.cereb_mask,
                    '--surface_folder', sub_files.surface_folder,
                    '--fsavgdir', fsavgDir,
                    '--surf'] + surf + [
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
             
        else:
            from simnibs.segmentation.run_cat_multiprocessing import run_cat_multiprocessing
            Ymf = nib.load(sub_files.norm_image)
            Yleft = nib.load(sub_files.hemi_mask)
            Yhemis = nib.load(sub_files.cereb_mask)
            
            run_cat_multiprocessing(Ymf.get_fdata(),
                                    Yleft.get_fdata(),
                                    Yhemis.get_fdata(),
                                    Ymf.affine, sub_files.surface_folder,
                                    fsavgDir, vdist, voxsize_pbt, 
                                    voxsize_refineCS, th_initial,
                                    no_selfintersections, surf, nprocesses)

        # print time duration
        elapsed = time.time() - starttime
        logger.info('Total time cost surface creation (HH:MM:SS):')
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
        field = np.zeros(mesh.elm.nr, dtype=np.int32)
        field = mesh.elm.tag1
        ed = ElementData(field)
        ed.mesh = mesh
        ed.to_deformed_grid(sub_files.mni2conf_nonl, MNI_template, 
                            out=sub_files.final_labels_MNI,
                            out_original=sub_files.final_labels,
                            method='assign',
                            reference_original=sub_files.reference_volume)

    # -------------------------TIDY UP-----------------------------------------

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


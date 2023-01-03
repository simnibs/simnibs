import logging
import os
import shutil
import nibabel as nib
import numpy as np
from functools import partial
from scipy import ndimage
import scipy.ndimage.morphology as mrph
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import affine_transform
from scipy.ndimage.measurements import label
from scipy.io import loadmat
import tempfile

from . import samseg
from ._thickness import _calc_thickness
from ._cat_c_utils import sanlm
from .brain_surface import mask_from_surface
from ..utils import file_finder
from ..utils.simnibs_logger import logger
from ..utils.transformations import resample_vol, volumetric_affine
from ..utils.spawn_process import spawn_process
from ..mesh_tools.mesh_io import read_off, write_off


def _register_atlas_to_input_affine(
    T1,
    template_file_name,
    affine_mesh_collection_name,
    mesh_level1,
    mesh_level2,
    save_path,
    template_coregistered_name,
    init_atlas_settings,
    neck_tissues,
    visualizer,
    noneck,
    init_transform=None,
    world_to_world_transform_matrix=None,
    scaling_center = [0.0, -100.0, 0.0],
    k_values = [20.0, 10.0, 5.0]
):

    # Import the affine registration function
    scales = init_atlas_settings["affine_scales"]
    thetas = init_atlas_settings["affine_rotations"]
    horizontal_shifts = init_atlas_settings["affine_horizontal_shifts"]
    vertical_shifts = init_atlas_settings["affine_vertical_shifts"]
    thetas_rad = [theta * np.pi / 180 for theta in thetas]
    neck_search_bounds = init_atlas_settings["neck_search_bounds"]
    ds_factor = init_atlas_settings["downsampling_factor_affine"]

    affine = samseg.AffineWholeHead(T1, affine_mesh_collection_name, template_file_name)

    init_options = samseg.initializationOptions(
        pitchAngles=thetas_rad,
        scales=scales,
        scalingCenter=scaling_center,
        horizontalTableShifts=horizontal_shifts,
        verticalTableShifts=vertical_shifts,
    )

    (
        image_to_image_transform,
        world_to_world_transform,
        optimization_summary,
    ) = affine.registerAtlas(
        worldToWorldTransformMatrix=world_to_world_transform_matrix,
        initTransform=init_transform,
        initializationOptions=init_options,
        targetDownsampledVoxelSpacing=ds_factor,
        visualizer=visualizer,
        noneck=noneck,
        Ks=k_values
    )

    affine.saveResults(
        T1,
        template_file_name,
        save_path,
        template_coregistered_name,
        image_to_image_transform,
        world_to_world_transform,
    )
    if world_to_world_transform_matrix is None:
        logger.info("Template registration summary.")
        logger.info(
            "Number of Iterations: %d, Cost: %f\n"
            % (optimization_summary["numberOfIterations"], optimization_summary["cost"])
        )

    if not noneck:
        logger.info("Adjusting neck.")
        exitcode = affine.adjust_neck(
            T1,
            template_coregistered_name,
            mesh_level1,
            mesh_level2,
            neck_search_bounds,
            neck_tissues,
            visualizer,
            downsampling_target=2.0,
        )
        if exitcode == -1:
            file_path = os.path.split(template_coregistered_name)
            shutil.copy(mesh_level1, os.path.join(file_path[0], "atlas_level1.txt.gz"))
            shutil.copy(mesh_level2, os.path.join(file_path[0], "atlas_level2.txt.gz"))
        else:
            logger.info("Neck adjustment done.")
    else:
        logger.info("No neck, copying meshes over.")
        file_path = os.path.split(template_coregistered_name)
        shutil.copy(mesh_level1, os.path.join(file_path[0], "atlas_level1.txt.gz"))
        shutil.copy(mesh_level2, os.path.join(file_path[0], "atlas_level2.txt.gz"))


def _denoise_input_and_save(input_name, output_name):
    input_raw = nib.load(input_name)
    img = input_raw.get_fdata()
    # Sometimes the images have an extra unit dimension,
    # squeeze that out if it's there.
    img = img.squeeze()
    img_smoothed = sanlm(img, 3, 1)
    output_smoothed = nib.Nifti1Image(img_smoothed, input_raw.affine)
    nib.save(output_smoothed, output_name)


def _init_atlas_affine(t1_scan, mni_template, affine_settings):

    registerer = samseg.gems.KvlAffineRegistration(
        affine_settings["translation_scale"],
        affine_settings["max_iter"],
        0,
        affine_settings["shrink_factors"],
        affine_settings["bg_value"],
        affine_settings["smoothing_factors"],
        affine_settings["center_of_mass"],
        affine_settings["samp_factor"],
        "b",
    )
    registerer.read_images(t1_scan, mni_template)
    registerer.initialize_transform()
    registerer.register()
    trans_mat = registerer.get_transformation_matrix()
    # registerer.write_out_result(os.path.join(path_to_segment_folder, 'mni_transformed.nii.gz'))
    logger.info(trans_mat)
    # ITK returns the matrix mapping the fixed image to the
    # moving image so let's invert it.
    return np.linalg.inv(trans_mat)


def _estimate_parameters(
    path_to_segment_folder,
    template_coregistered_name,
    path_to_atlas_folder,
    input_images,
    segment_settings,
    gmm_parameters,
    visualizer,
    user_optimization_options=None,
    user_model_specifications=None
):

    ds_targets = segment_settings["downsampling_targets"]
    kernel_size = segment_settings["bias_kernel_width"]
    bg_mask_sigma = segment_settings["background_mask_sigma"]
    bg_mask_th = segment_settings["background_mask_threshold"]
    stiffness = segment_settings["mesh_stiffness"]
    covariances = segment_settings["diagonal_covariances"]
    shared_gmm_parameters = samseg.io.kvlReadSharedGMMParameters(gmm_parameters)

    if user_optimization_options is None:
        user_optimization_options = {
            "multiResolutionSpecification": [
                {
                    "atlasFileName": os.path.join(
                        path_to_segment_folder, "atlas_level1.txt.gz"
                    ),
                    "targetDownsampledVoxelSpacing": ds_targets[0],
                    "maximumNumberOfIterations": 100,
                    "estimateBiasField": True,
                },
                {
                    "atlasFileName": os.path.join(
                        path_to_segment_folder, "atlas_level2.txt.gz"
                    ),
                    "targetDownsampledVoxelSpacing": ds_targets[1],
                    "maximumNumberOfIterations": 100,
                    "estimateBiasField": True,
                },
            ]
        }

    if user_model_specifications is None:
        user_model_specifications = {
            "atlasFileName": os.path.join(path_to_segment_folder, "atlas_level2.txt.gz"),
            "biasFieldSmoothingKernelSize": kernel_size,
            "brainMaskingSmoothingSigma": bg_mask_sigma,
            "brainMaskingThreshold": bg_mask_th,
            "K": stiffness,
            "useDiagonalCovarianceMatrices": covariances,
            "sharedGMMParameters": shared_gmm_parameters,
        }

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
        saveWarp=True,
    )

    logger.info("Starting segmentation.")
    samsegment = samseg.SamsegWholeHead(**samseg_kwargs)
    samsegment.preProcess()
    samsegment.process()

    # Print optimization summary
    optimizationSummary = samsegment.getOptimizationSummary()
    for multiResolutionLevel, item in enumerate(optimizationSummary):
        logger.info(
            "atlasRegistrationLevel%d %d %f\n"
            % (multiResolutionLevel, item["numberOfIterations"], item["perVoxelCost"])
        )

    return samsegment.saveParametersAndInput()


def _post_process_segmentation(
    bias_corrected_image_names,
    upsampled_image_names,
    tissue_settings,
    csf_factor,
    parameters_and_inputs,
    transformed_template_name,
    affine_atlas,
    before_morpho_name,
    upper_mask,
    debug=False,
):

    logger.info("Upsampling bias corrected images.")
    for input_number, bias_corrected in enumerate(bias_corrected_image_names):
        corrected_input = nib.load(bias_corrected)
        resampled_input, new_affine, orig_res = resample_vol(
            corrected_input.get_data(), corrected_input.affine, 0.5, order=1
        )
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
        csf_factor,
    )

    del parameters_and_inputs
    # Cast the upsampled image to int16  to save space
    for upsampled_image in upsampled_image_names:
        upsampled = nib.load(upsampled_image)
        upsampled.set_data_dtype(np.int16)
        nib.save(upsampled, upsampled_image)

    affine_upsampled = upsampled.affine
    if debug:
        upsampled_tissues_im = nib.Nifti1Image(upsampled_tissues, affine_upsampled)
        nib.save(upsampled_tissues_im, before_morpho_name)

    upper_part_im = nib.Nifti1Image(upper_part.astype(np.int16), affine_upsampled)
    upper_part_im.set_data_dtype(np.int16)
    nib.save(upper_part_im, upper_mask)
    del upper_part_im
    # Do morphological operations
    simnibs_tissues = tissue_settings["simnibs_tissues"]
    upsampled_tissues = _morphological_operations(
        upsampled_tissues, upper_part, simnibs_tissues
    )

    return upsampled_tissues


def _clean_brain(label_img, tissues, unass, se, se_n, vol_limit=10):
    # Do some clean-ups, mainly CSF and skull
    # First combine the WM and GM
    brain = (label_img == tissues["WM"]) | (label_img == tissues["GM"])
    dil = mrph.binary_dilation(
        _get_largest_components(
            mrph.binary_erosion(brain, se_n, 1), se, vol_limit
        ),
        se,
        1,
    )

    unass |= brain ^ dil

    # Add the CSF and open again
    csf = label_img == tissues["CSF"]
    brain_csf = brain | csf
    # Vol limit in voxels
    dil = mrph.binary_dilation(
        _get_largest_components(
            mrph.binary_erosion(brain_csf, se, 1), se_n, vol_limit=80
        ),
        se,
        1,
    )
    unass |= csf & ~dil
    del brain, csf, dil
    return brain_csf


def _get_skull(label_img, brain_csf, tissues, se_n, se, num_iter=2, vol_limit=50):
    # Clean the outer border of the skull
    bone = (label_img == tissues["Compact_bone"]) | (
        label_img == tissues["Spongy_bone"]
    )
    veins = label_img == tissues["Blood"]
    air_pockets = label_img == tissues["Air_pockets"]
    scalp = label_img == tissues["Scalp"]
    muscle = label_img == tissues["Muscle"]
    eyes = label_img == tissues["Eyes"]
    # Use scalp to clean out noisy skull bits within the scalp
    skull_outer = brain_csf | bone | veins | air_pockets
    skull_outer = mrph.binary_fill_holes(skull_outer, se_n)
    skull_outer = mrph.binary_dilation(
        _get_largest_components(
            mrph.binary_erosion(skull_outer, se, num_iter), se_n, vol_limit
        ),
        se,
        num_iter,
    )
    skull_inner = bone | scalp | air_pockets | muscle | eyes
    skull_inner = mrph.binary_fill_holes(skull_inner, se_n)
    skull_inner = mrph.binary_dilation(
        _get_largest_components(
            mrph.binary_erosion(skull_inner, se, num_iter), se_n, vol_limit
        ),
        se,
        num_iter)

    dil = bone & skull_outer & skull_inner
    return bone, skull_outer, dil


def _clean_veins(label_img, unass, tissues, se, se_n, num_iter=1, vol_limit=10):
    # Open the veins
    veins = label_img == tissues["Blood"]
    dil = mrph.binary_dilation(
        _get_largest_components(mrph.binary_erosion(veins, se, num_iter), se_n, vol_limit),
        se,
        num_iter,
    )
    unass |= dil ^ veins


def _clean_eyes(label_img, unass, tissues, se, se_n, num_iter=1, vol_limit=10):
    # Clean the eyes
    eyes = label_img == tissues["Eyes"]
    dil = mrph.binary_dilation(
        _get_largest_components(mrph.binary_erosion(eyes, se, num_iter), se_n, vol_limit),
        se,
        num_iter,
    )
    # dil = mrph.binary_opening(eyes, se, 1)
    unass |= dil ^ eyes


def _clean_muscles(label_img, unass, tissues, se, num_iter=1):
    # Clean muscles
    muscle = label_img == tissues["Muscle"]
    dil = mrph.binary_opening(muscle, se, num_iter)
    unass |= dil ^ muscle


def _clean_scalp(label_img, unass, skull_outer, tissues, se, se_n, num_iter=2, num_limit=1):
    # And finally the scalp
    scalp = label_img == tissues["Scalp"]
    eyes = label_img == tissues["Eyes"]
    muscle = label_img == tissues["Muscle"]
    head = scalp | skull_outer | eyes | muscle
    dil = mrph.binary_dilation(
        _get_largest_components(mrph.binary_erosion(head, se, num_iter), se_n, num_limit),
        se,
        num_iter,
    )
    unass |= scalp & ~dil


def _ensure_csf(label_img, tissues, upper_part, se, num_iter1=1, num_iter2=6):
    # Ensuring a CSF layer between GM and Skull and GM and Blood
    # Relabel regions in the expanded GM which are in skull or blood to CSF
    logger.info("Ensure CSF")

    brain_gm = label_img == tissues["GM"]
    C_BONE = label_img == tissues["Compact_bone"]
    S_BONE = label_img == tissues["Spongy_bone"]
    U_SKIN = (label_img == tissues["Scalp"]) & upper_part

    brain_dilated = mrph.binary_dilation(brain_gm, se, num_iter1)
    overlap = brain_dilated & (C_BONE | S_BONE | U_SKIN)
    label_img[overlap] = tissues["CSF"]

    S_BONE = label_img == tissues["Spongy_bone"]

    CSF_brain = brain_gm | (label_img == tissues["CSF"])
    CSF_brain_dilated = mrph.binary_dilation(CSF_brain, se, num_iter1)
    spongy_csf = CSF_brain_dilated & S_BONE
    upper_part = mrph.binary_erosion(upper_part, se, num_iter2)
    skin_csf = CSF_brain_dilated & (label_img == tissues["Scalp"]) & upper_part
    label_img[spongy_csf] = tissues["Compact_bone"]
    label_img[skin_csf] = tissues["Compact_bone"]


def _ensure_skull(label_img, tissues, se, num_iter=1):
    # Ensure the outer skull label is compact bone
    C_BONE = label_img == tissues["Compact_bone"]
    S_BONE = label_img == tissues["Spongy_bone"]
    SKULL_outer = (C_BONE | S_BONE) & ~mrph.binary_erosion((C_BONE | S_BONE), se, num_iter)
    label_img[SKULL_outer] = tissues["Compact_bone"]
    # Relabel air pockets to air
    label_img[label_img == tissues["Air_pockets"]] = 0



def _morphological_operations(label_img, upper_part, simnibs_tissues):
    """Does morphological operations to
    1. Smooth out the labeling and remove noise
    2. A CSF layer between GM and Skull and between GM and CSF
    3. Outer bone layers are compact bone
    """
    se = ndimage.generate_binary_structure(3, 3)
    se_n = ndimage.generate_binary_structure(3, 1)
    unass = np.zeros_like(label_img) > 0
    brain_csf = _clean_brain(label_img, simnibs_tissues, unass, se, se_n)
    bone, skull_outer, dil = _get_skull(label_img, brain_csf, simnibs_tissues, se_n, se)
    # Protect thin areas that would be removed by erosion
    bone_thickness = _calc_thickness(dil)

    thin_parts = bone & (bone_thickness < 3.5) & (bone_thickness > 0)
    dil |= thin_parts
    del thin_parts
    del bone_thickness

    unass = unass | (dil ^ bone)
    del bone, dil

    _clean_veins(label_img, unass, simnibs_tissues, se, se_n)
    _clean_eyes(label_img, unass, simnibs_tissues, se, se_n)
    _clean_muscles(label_img, unass, simnibs_tissues, se)
    _clean_scalp(label_img, unass, skull_outer, simnibs_tissues, se, se_n)

    # Filling missing parts
    # NOTE: the labeling is uint16, so I'll code the unassigned voxels with the max value
    # which is 65535
    label_unassign = 65535
    label_img[unass] = label_unassign
    # Add background to tissues
    simnibs_tissues["BG"] = 0
    label_img = label_unassigned_elements(
        label_img, label_unassign, list(simnibs_tissues.values()) + [label_unassign]
    )
    _smoothfill(label_img, unass, simnibs_tissues)

    _ensure_csf(label_img, simnibs_tissues, upper_part, se)
    _ensure_skull(label_img, simnibs_tissues, se)

    return label_img


def generate_gaussian_kernel(i, ndim, sigma=1, zero_center=False):
    """Generate a gaussian kernel of size `i` in `ndim` dimensions."""
    assert i % 2 == 1, "Window size must be odd"
    center = ndim * (np.floor(i / 2).astype(int),)
    kernel = np.zeros(ndim * (i,), dtype=np.float32)
    kernel[center] = 1
    kernel = gaussian_filter(kernel, sigma)
    kernel /= kernel[center]
    if zero_center:
        kernel[center] = 0
    return kernel


def label_unassigned_elements(
    label_arr, label_unassign, labels=None, window_size=3, ignore_labels=None
) -> np.ndarray:
    """Label unassigned elements in `label_arr`. For each unassigned element,
    find its neighbors within a certain `window_size`, weigh these according
    to their euclidean distance (Gaussian kernel) and assign the label with the
    highest weight.

    PARAMETERS
    ----------
    label_arr : ndarray
        The array
    label_unassign : int
        The label of the unassigned elements.
    labels : array-like | None
        The labels in `label_arr`. If None (default), will be inferred from the
        array using `np.unique`.
    window_size : int
        The size of the neighborhood around each element to consider when
        relabeling.
    ignore_labels : array-like | None
        Do not use these labels when labeling unassigned elements (default =
        None).

    RETURNS
    -------
    labeling : ndarray
        The relabeled array.
    """

    # Finding the indices of the unassigned voxels is slow with fortran ordered
    # arrays
    labeling = label_arr.T if label_arr.flags["F_CONTIGUOUS"] else label_arr

    # Weighting kernel and related
    pad_width = np.floor(window_size / 2).astype(int)
    kernel = generate_gaussian_kernel(window_size, labeling.ndim, zero_center=True)
    window = kernel.shape
    kernel = kernel.ravel()

    # Find the list of label (if necessary) and ensure it is unique
    labels = np.unique(labeling) if labels is None else np.unique(labels)
    labels = labels.astype(labeling.dtype)
    n_labels = labels.size
    assert label_unassign in labels, "`label_unassign` is not in the provided `labels`!"

    # Ignore desired labels and the unassigned label
    if ignore_labels is None:
        ignore_labels = [label_unassign]
    else:
        assert all(i in labels for i in ignore_labels)
        ignore_labels = ignore_labels.copy()
        ignore_labels = ignore_labels + [label_unassign]  # avoid inplace

    # Ensure continuous labels (e.g.,, [0, 1, 2, 3], not [0, 1, 10, 15])
    is_continuous_labels = labels[-1] - labels[0] + 1 == n_labels
    if is_continuous_labels:
        continuous_labels = labels
        mapped_ignore_labels = np.array(ignore_labels)
        labeling = labeling.copy()  # ensure that we do not modify the input
    else:
        continuous_labels = np.arange(n_labels)
        mapper = np.zeros(labels[-1] + 1, dtype=labeling.dtype)
        mapper[labels] = continuous_labels
        labeling = mapper[labeling]
        mapped_ignore_labels = mapper[ignore_labels]

    is_unassign = np.nonzero(labeling == mapped_ignore_labels[-1])
    while (n_unassign := is_unassign[0].size) > 0:
        # print("Number of unassigned voxels:", n_unassign)
        labeling_view = np.lib.stride_tricks.sliding_window_view(
            np.pad(labeling, pad_width, "symmetric"), window
        )
        unassign_window = labeling_view[is_unassign].reshape(-1, kernel.size)

        # Compute weights in a loop to save memory
        weights = np.array(
            [
                np.zeros(n_unassign)
                if i in mapped_ignore_labels
                else (unassign_window == i) @ kernel
                for i in continuous_labels
            ]
        )
        # The slightly faster but less memory-friendly solution
        # n_total = unassign_window.size
        # oh_enc = np.zeros((n_labels, n_total), dtype=np.uint8)
        # oh_enc[unassign_window.ravel(), np.arange(n_total)] = 1
        # oh_enc = oh_enc.reshape(n_labels, *unassign_window.shape)
        # weights = oh_enc @ kernel
        # if mapped_ignore_labels.size > 0:
        #     weights[mapped_ignore_labels] = 0

        # In the case of ties, the lowest index is returned. This is arbitrary
        # but the default behavior of argmax
        new_label = weights.argmax(0)
        valid = weights.sum(0) > 0
        invalid = ~valid

        labeling[tuple(i[valid] for i in is_unassign)] = new_label[valid]
        is_unassign = tuple(i[invalid] for i in is_unassign)

        if is_unassign[0].size == n_unassign:
            logger.warning(
                "Some elements could not be labeled (probably because they are surrounded by labels in `ignore_labels`)"
            )
            # perhaps we may want to issue the warning using the warnings
            # module instead
            # warnings.warn(
            #     "Some elements could not be labeled, probably because they are surrounded by labels in `ignore_labels`",
            #     RuntimeWarning,
            # )
            break

    # Revert array
    if not is_continuous_labels:
        labeling = labels[labeling]
    labeling = labeling.T if label_arr.flags["F_CONTIGUOUS"] else labeling

    return labeling

def _smooth(label_img, simnibs_tissues, tissues_to_smooth):
    """Smooth some of the tissues for "nicer" tissue labelings.
    """
    labs = label_img.copy()
    max_val = np.zeros_like(label_img, dtype=np.float32)
    for i, t in enumerate(simnibs_tissues):
        vol = (label_img == simnibs_tissues[t]).astype(np.float32)
        if t in tissues_to_smooth:
            cs = gaussian_filter(vol, 1)
        else:
            cs = vol

        # Check the max values and update
        max_mask = cs > max_val
        labs[max_mask] = tissues_to_smooth[t]
        max_val[max_mask] = cs[max_mask]

    label_img[:] = labs[:]
    del labs


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
        if sum_of_unassigned == 0:
            break
        labs = 65535 * np.ones_like(label_img)
        max_val = np.zeros_like(label_img, dtype=np.float32)
        for i, t in enumerate(simnibs_tissues):
            # Don't smooth WM
            vol = (label_img == simnibs_tissues[t]).astype(np.float32)
            if t == "WM":
                cs = vol
            else:
                cs = gaussian_filter(vol, 1)

            # Check the max values and update
            max_mask = cs > max_val
            labs[max_mask] = simnibs_tissues[t]
            max_val[max_mask] = cs[max_mask]

        label_img[:] = labs[:]
        unassign = labs == 65535
        del labs


def _fill_missing(label_img, unassign):
    """Hackish way to fill unassigned voxels,
    works by smoothing the masks and binarizing
    the smoothed masks. Works in-place
    """

    def _get_most_frequent_label(x, y, z, label_im, pad):
        nhood = label_im[
            max(0, x - pad) : min(label_im.shape[0], x + pad),
            max(0, y - pad) : min(label_im.shape[1], y + pad),
            max(0, z - pad) : min(label_im.shape[2], z + pad),
        ]

        dim = nhood.shape

        return np.bincount(nhood.reshape((np.prod(dim)))).argmax()

    label_img[unassign] = 255
    inds_tmp = unassign.nonzero()
    inds_tmp = np.array(inds_tmp)
    num_unassigned = (label_img == 255).sum()
    num_unassigned_new = num_unassigned
    while num_unassigned_new > 0:
        fill_map = map(
            partial(_get_most_frequent_label, label_im=label_img, pad=5),
            np.nditer(inds_tmp[0, :]),
            np.nditer(inds_tmp[1, :]),
            np.nditer(inds_tmp[2, :]),
        )
        fill_array = np.array(list(fill_map))
        label_img[inds_tmp[0, :], inds_tmp[1, :], inds_tmp[2, :]] = fill_array
        inds_tmp = inds_tmp[:, fill_array == 255]
        num_unassigned_new = (label_img == 255).sum()
        logger.info("Unassigned: " + str(num_unassigned_new))
        if num_unassigned_new == num_unassigned:
            logger.info("Number of unassigned voxels not going down. Breaking.")
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
    components = np.zeros_like(vol, dtype=bool)
    for i in labels[:num_limit]:
        components = components | (vol_lbl == i)

    if return_sizes:
        return components, region_size[labels[:num_limit]]
    else:
        return components


def _registerT1T2(fixed_image, moving_image, output_image):
    registerer = samseg.gems.KvlRigidRegistration()
    registerer.read_images(fixed_image, moving_image)
    registerer.initialize_transform()
    registerer.register()
    registerer.write_out_result(output_image)

    # The register function uses double internally
    # Let's cast to float32 and copy the header from
    # the fixed image to avoid any weird behaviour
    if os.path.exists(output_image):
        T2_reg = nib.load(output_image)
        fixed_tmp = nib.load(fixed_image)
        T2_data = T2_reg.get_fdata().astype(np.float32)
        T2_im = nib.Nifti1Image(T2_data, fixed_tmp.affine)
        nib.save(T2_im, output_image)


def _fillin_gm_layer(
    label_img,
    label_affine,
    labelorg_img,
    labelorg_affine,
    m,
    exclusion_tissues={
        "left_cerebral_wm": 2,
        "right_cerebral_wm": 41,
        "stuff_to_exclude": [4, 10, 14, 16, 24, 28, 43, 49, 60],
    },
    relabel_tissues={"GM": 2, "stuff_to_relabel": [1, 3]},
):
    """relabels WM and CSF that intersect with the central GM surface to GM
    an exclusion mask is used to prevent relabelign in central brain regions
    """
    # generate exclusion mask: estimate corpus callossum
    exclude_img = mrph.binary_dilation(
        labelorg_img == exclusion_tissues["left_cerebral_wm"], iterations=2
    )
    exclude_img *= mrph.binary_dilation(
        labelorg_img == exclusion_tissues["right_cerebral_wm"], iterations=2
    )
    # add other tissues
    for i in exclusion_tissues["stuff_to_exclude"]:
        exclude_img += labelorg_img == i
    exclude_img = mrph.binary_dilation(exclude_img, iterations=8)
    # upsample exclude_img
    iM = np.linalg.inv(labelorg_affine).dot(label_affine)
    exclude_img = affine_transform(
        exclude_img, iM[:3, :3], iM[:3, 3], label_img.shape, order=0
    )

    # generate voxel mask of middle GM
    mask = mask_from_surface(
        m.nodes[:], m.elm[:, :3] - 1, label_affine, label_img.shape
    )
    mask = mrph.binary_dilation(mask, iterations=1) * ~mrph.binary_erosion(
        mask, iterations=1
    )
    mask[exclude_img] = 0

    # relabel WM and CSF parts to GM
    for i in relabel_tissues["stuff_to_relabel"]:
        label_img[(label_img == i) * mask] = relabel_tissues["GM"]

    return label_img


def _open_sulci(
    label_img,
    label_affine,
    labelorg_img,
    labelorg_affine,
    m,
    tissue_labels={"CSF": 3, "GM": 2, "WM": 1},
    exclusion_tissues=[17, 18, 53, 54],
):
    # get thin CSF structures
    mask = mask_from_surface(
        m.nodes[:], m.elm[:, :3] - 1, label_affine, label_img.shape
    )
    # mask2 = mrph.binary_dilation(mask,iterations=4)
    # mask2 = mrph.binary_erosion(mask2,iterations=5)
    mask2 = mrph.binary_dilation(mask, iterations=2)
    mask2 = mrph.binary_erosion(mask2, iterations=3)
    mask2[mask] = 0

    # protect hippocampi and amydalae
    exclude_img = np.zeros_like(labelorg_img, dtype=bool)
    for i in exclusion_tissues:
        exclude_img += labelorg_img == i
    iM = np.linalg.inv(labelorg_affine).dot(label_affine)
    exclude_img = affine_transform(
        exclude_img, iM[:3, :3], iM[:3, 3], label_img.shape, order=0
    )
    mask2[exclude_img] = 0

    # relabel GM overlapping thin CSF to CSF
    label_img[(label_img == tissue_labels["GM"]) * mask2] = tissue_labels["CSF"]

    # open up remaining thin GM bridges at brain surface
    mask2 = (label_img == tissue_labels["GM"]) | (label_img == tissue_labels["WM"])
    # mask2 = mrph.binary_erosion(mask2,iterations=2)
    # mask2 = mrph.binary_dilation(mask2,iterations=2)
    # label_img[ (label_img == tissue_labels['GM'])* ~mask2 ] = tissue_labels['CSF']
    brainthickness = _calc_thickness(mask2)
    label_img[
        (brainthickness <= 2.0) * (label_img == tissue_labels["GM"])
    ] = tissue_labels["CSF"]

    return label_img


def _cut_and_combine_labels(fn_tissue_labeling_upsampled, fn_mni_template, 
                            fn_affine, tms_settings):
    """
    Cut away neck of tissue_labeling_upsampled.nii.gz and 
    combine some of the labels. Overwrites the original file.
    Used to create meshes optimized for TMS.

    Parameters
    ----------
    fn_tissue_labeling_upsampled : string
        filename of tissue_labeling_upsampled.nii.gz
    fn_mni_template : string
        filename of the MNI template
    fn_affine : string
        filename of the affine transformation from T1 to MNI
    tms_settings : dict
        specifies old and new labels

    Returns
    -------
    None.

    """
    # cut label image using MNI mask
    logger.info("Cutting neck region using MNI mask")
    label_image  = nib.load(fn_tissue_labeling_upsampled)
    label_buffer = np.round(label_image.get_fdata()).astype( np.uint16 )  # Cast to uint16, otherwise meshing complains
    label_affine = label_image.affine
    
    mni_image  = nib.load(fn_mni_template)
    mni_buffer = np.ones(mni_image.shape, dtype=bool)
    mni_affine = mni_image.affine
    
    trafo = loadmat(fn_affine)['worldToWorldTransformMatrix']

    upperhead = volumetric_affine((mni_buffer, mni_affine),
                                  np.linalg.inv(trafo),
                                  target_space_affine=label_affine,
                                  target_dimensions=label_image.shape,
                                  intorder=0)
    label_buffer[~upperhead] = 0
    
    # combine labels
    logger.info("Combining labels")
    for idx_old, idx_new in zip(tms_settings["old_label"], tms_settings["new_label"]):
        logger.debug("  old label: %d, new label: %d " % (idx_old, idx_new))
        label_buffer[label_buffer == idx_old] = idx_new
  
    label_image = nib.Nifti1Image(label_buffer, label_affine)
    nib.save(label_image, fn_tissue_labeling_upsampled)
    
 
def _downsample_surface(m, n_nodes):
    """
        downsample a surface using meshfix

    Parameters
    ----------
    m : simnibs.Msh 
        surface
    n_nodes : int
        target number of nodes.

    Returns
    -------
    mout : simnibs.Msh
        downsampled surface
        
    NOTE: this is a primitive wrapper around meshfix, tag1 and tag2 of the
    returned mesh will be set to 1, meshes with multiple surfaces are not 
    supported
    """
    with tempfile.NamedTemporaryFile(suffix=".off") as f:
        mesh_fn = f.name
    write_off(m, mesh_fn)
    cmd = [file_finder.path2bin("meshfix"), mesh_fn, 
           "-u", "2", "--vertices", str(n_nodes), "-o", mesh_fn]
    spawn_process(cmd, lvl=logging.DEBUG)
    
    mout = read_off(mesh_fn)
    os.remove(mesh_fn)
    if os.path.isfile("meshfix_log.txt"):
        os.remove("meshfix_log.txt")
    return mout

import charm_gems as gems
from .ProbabilisticAtlas import ProbabilisticAtlas
from .SamsegUtility import undoLogTransformAndBiasField, writeImage, maskOutBackground, logTransform, readCroppedImages
from .GMM import GMM
import numpy as np
import nibabel as nib
from simnibs.segmentation._cat_c_utils import cat_vbdist
from simnibs.segmentation.samseg.utilities import requireNumpyArray

# TODO! remove
import os

def writeBiasCorrectedImagesAndSegmentation(output_names_bias,
                                            output_name_segmentation,
                                            parameters_and_inputs,
                                            cat_structure_options=None,
                                            cat_images=None):

    # We need an init of the probabilistic segmentation class
    # to call instance methods
    probabilisticAtlas = ProbabilisticAtlas()
    # Bias correct images
    biasFields = parameters_and_inputs['biasFields']
    imageBuffers = parameters_and_inputs['imageBuffers']
    mask = parameters_and_inputs['mask']

    # Write these out
    exampleImageFileName = parameters_and_inputs['imageFileNames']
    exampleImage = gems.KvlImage(exampleImageFileName[0])
    cropping = parameters_and_inputs['cropping']
    expImageBuffers, expBiasFields = undoLogTransformAndBiasField(
                                        imageBuffers,
                                        biasFields,
                                        mask)
    for contrastNumber, out_name in enumerate(output_names_bias):
        # Bias field correct and write
        writeImage(out_name,  expImageBuffers[..., contrastNumber],
                   cropping, exampleImage)

    # Next do the segmentation, first read in the mesh
    modelSpecifications = parameters_and_inputs['modelSpecifications']
    transform_matrix = parameters_and_inputs['transform']
    transform = gems.KvlTransform(requireNumpyArray(transform_matrix))
    deformation = parameters_and_inputs['deformation']
    deformationAtlasFileName = parameters_and_inputs['deformationAtlas']

    mesh = probabilisticAtlas.getMesh(
            modelSpecifications.atlasFileName,
            transform=transform,
            K=modelSpecifications.K,
            initialDeformation=deformation,
            initialDeformationMeshCollectionFileName=deformationAtlasFileName)

    fractionsTable = parameters_and_inputs['fractionsTable']
    GMMparameters = parameters_and_inputs['GMMParameters']
    numberOfGaussiansPerClass = parameters_and_inputs['gaussiansPerClass']
    means = GMMparameters['means']
    variances = GMMparameters['variances']
    mixtureWeights = GMMparameters['mixtureWeights']
    names = parameters_and_inputs['names']
    bg_label = names.index("Background")
    FreeSurferLabels = np.array(modelSpecifications.FreeSurferLabels,
                                dtype=np.uint16)

    if cat_structure_options is not None:
        segmentation, cat_norm_scan_masks = _calculateSegmentationLoop(
                                            imageBuffers - biasFields,
                                            mask, fractionsTable, mesh,
                                            numberOfGaussiansPerClass,
                                            means, variances, mixtureWeights,
                                            FreeSurferLabels, bg_label,
                                            cat_opts=cat_structure_options)
    else:
        segmentation = _calculateSegmentationLoop(
                        imageBuffers - biasFields,
                        mask, fractionsTable, mesh,
                        numberOfGaussiansPerClass,
                        means, variances, mixtureWeights,
                        FreeSurferLabels, bg_label)
        
    writeImage(output_name_segmentation,
               segmentation, cropping, exampleImage)
    
    if cat_structure_options is not None and cat_images is not None:
        for inds, name in enumerate(cat_images):
            writeImage(name, cat_norm_scan_masks[inds], cropping, exampleImage)


def segmentUpsampled(input_bias_corrected, tissue_settings,
                     parameters_and_inputs, transformedTemplateFileName):

    # We need an init of the probabilistic segmentation class
    # to call instance methods
    probabilisticAtlas = ProbabilisticAtlas()

    # Read the input images doing the necessary cropping
    imageBuffersUpsampled, transformUpsampled, voxelSpacingUpsampled, croppingUpsampled = readCroppedImages(
                                                                                              input_bias_corrected,
                                                                                              transformedTemplateFileName)
    # Redo background masking now with the upsampled scans, note this only rasterizes
    # a single class so should be decent memory-wise
    modelSpecifications = parameters_and_inputs['modelSpecifications']
    imageBuffersUpsampled, maskUpsampled = maskOutBackground(
            imageBuffersUpsampled,
            modelSpecifications.atlasFileName,
            transformUpsampled,
            modelSpecifications.brainMaskingSmoothingSigma,
            modelSpecifications.brainMaskingThreshold,
            probabilisticAtlas)

    # Log-transform the intensities, note the scans are already
    # bias corrected so no need to remove the bias contribution

    imageBuffersUpsampled = logTransform(imageBuffersUpsampled,
                                         maskUpsampled)

    # Calculate the posteriors.
    # NOTE: the method for calculating the posteriors
    # in the gmm class rasterizes all classes in the atlas at once.
    # This is very memory-heavy if we have many classes
    # and the input resolution is high. Here we do it instead in
    # a loop to save memory.
    deformation = parameters_and_inputs['deformation']
    deformationAtlasFileName = parameters_and_inputs['deformationAtlas']
    meshUpsampled = probabilisticAtlas.getMesh(
                modelSpecifications.atlasFileName,
                transformUpsampled,
                initialDeformation=deformation,
                initialDeformationMeshCollectionFileName=deformationAtlasFileName)

    fractionsTable = parameters_and_inputs['fractionsTable']
    GMMparameters = parameters_and_inputs['GMMParameters']
    numberOfGaussiansPerClass = parameters_and_inputs['gaussiansPerClass']
    means = GMMparameters['means']
    variances = GMMparameters['variances']
    mixtureWeights = GMMparameters['mixtureWeights']
    names = parameters_and_inputs['names']
    bg_label = names.index("Background")
    FreeSurferLabels = np.array(modelSpecifications.FreeSurferLabels,
                                dtype=np.uint16)
    segmentation = _calculateSegmentationLoop(
                imageBuffersUpsampled,
                maskUpsampled, fractionsTable, meshUpsampled,
                numberOfGaussiansPerClass,
                means, variances, mixtureWeights, FreeSurferLabels, bg_label)

    simnibs_tissues = tissue_settings['simnibs_tissues']
    segmentation_tissues = tissue_settings['segmentation_tissues']

    tissue_labeling = np.zeros_like(segmentation)
    for t, label in simnibs_tissues.items():
        tissue_labeling[np.isin(segmentation,
                                FreeSurferLabels[segmentation_tissues[t]])] = label

    example_image = gems.KvlImage(input_bias_corrected[0])
    uncropped_tissue_labeling = np.zeros(
            example_image.getImageBuffer().shape, dtype=np.uint16, order='F')
    uncropped_tissue_labeling[croppingUpsampled] = tissue_labeling

    #upsampled_image = nib.load(input_bias_corrected[0])
    #affine_upsampled = upsampled_image.affine
    #upsampled_tissues = nib.Nifti1Image(uncropped_tissue_labeling,
    #                                    affine_upsampled)
    #nib.save(upsampled_tissues, './test_labeling_mem.nii.gz')
    return uncropped_tissue_labeling


def saveWarpField(template_name, warp_to_mni,
                  warp_from_mni, parameters_and_inputs):
    # Save warp field two ways: the subject space world coordinates
    # in template space, i.e., the iamge voxel coordinates in
    # physical space for every voxel in the template.
    # And the other way around, i.e., the template voxel coordinates
    # in physical space for every voxel in the image.

    # We need an init of the probabilistic segmentation class
    # to call instance methods
    probabilisticAtlas = ProbabilisticAtlas()

    # First write image -> template.
    # Get the node positions in image voxels
    modelSpecifications = parameters_and_inputs['modelSpecifications']
    transform_matrix = parameters_and_inputs['transform']
    transform = gems.KvlTransform(requireNumpyArray(transform_matrix))

    deformation = parameters_and_inputs['deformation']
    deformationAtlasFileName = parameters_and_inputs['deformationAtlas']
    nodePositions = probabilisticAtlas.getMesh(
             modelSpecifications.atlasFileName,
             transform,
             K=modelSpecifications.K,
             initialDeformation=deformation,
             initialDeformationMeshCollectionFileName=deformationAtlasFileName
         ).points

    # The image is cropped as well so the voxel coordinates
    # do not exactly match with the original image,
    # i.e., there's a shift. Let's undo that.
    cropping = parameters_and_inputs['cropping']
    nodePositions += [slc.start for slc in cropping]

    # Get mapping from voxels to world space of the image.
    imageFileNames = parameters_and_inputs['imageFileNames']
    image = nib.load(imageFileNames[0])
    imageToWorldTransformMatrix = image.affine
    image_buffer = image.get_data()
    # Transform the node positions
    nodePositionsInWorldSpace = (imageToWorldTransformMatrix @
                                 np.pad(nodePositions,
                                        ((0, 0), (0, 1)),
                                        'constant', constant_values=1).T).T
    nodePositionsInWorldSpace = nodePositionsInWorldSpace[:, 0:3]
    template = gems.KvlImage(template_name)
    templateToWorldTransformMatrix = template.transform_matrix.as_numpy_array
    template_buffer = template.getImageBuffer()

    # Rasterize the final node coordinates (in image space)
    # using the initial template mesh
    mesh = probabilisticAtlas.getMesh(
                modelSpecifications.atlasFileName,
                K=modelSpecifications.K)

    # Get node positions in template voxel space
    nodePositionsTemplate = mesh.points

    # Rasterize the coordinate values
    coordmapTemplate = mesh.rasterize_values(template_buffer.shape,
                                             nodePositionsInWorldSpace)
    # Write the warp file
    temp_header = nib.load(template_name)
    warp_image = nib.Nifti1Image(coordmapTemplate, temp_header.affine)
                                 #templateToWorldTransformMatrix)
    nib.save(warp_image, warp_to_mni)

    # Now do it the other way, i.e., template->image
    nodePositionsTemplateWorldSpace = (temp_header.affine @
                                       np.pad(nodePositionsTemplate,
                                              ((0, 0), (0, 1)),
                                              'constant', constant_values=1).T).T
    nodePositionsTemplateWorldSpace = nodePositionsTemplateWorldSpace[:, 0:3]
    # Okay get the mesh in image space
    mesh = probabilisticAtlas.getMesh(
             modelSpecifications.atlasFileName,
             transform,
             initialDeformation=deformation,
             initialDeformationMeshCollectionFileName=deformationAtlasFileName)

    imageBuffers = parameters_and_inputs['imageBuffers']
    coordmapImage = mesh.rasterize_values(imageBuffers.shape[0:-1],
                                          nodePositionsTemplateWorldSpace)
    # The image buffer is cropped so need to set
    # everything to the correct place
    uncroppedWarp = np.zeros(image_buffer.shape+(3,),
                             dtype=np.float32, order='F')

    for c in range(coordmapImage.shape[-1]):
        uncroppedMap = np.zeros(image_buffer.shape,
                                dtype=np.float32, order='F')
        uncroppedMap[cropping] = np.squeeze(coordmapImage[:, :, :, c])
        uncroppedWarp[:, :, :, c] = uncroppedMap

    # Write the warp
    warp_image = nib.Nifti1Image(uncroppedWarp,
                                 imageToWorldTransformMatrix)
    nib.save(warp_image, warp_from_mni)


def _calculateSegmentationLoop(biasCorrectedImageBuffers,
                               mask, fractionsTable, mesh,
                               numberOfGaussiansPerClass,
                               means, variances,
                               mixtureWeights, FreeSurferLabels,
                               bg_label,
                               cat_opts=None):

    data = biasCorrectedImageBuffers[mask, :]
    numberOfVoxels = data.shape[0]
    numberOfStructures = fractionsTable.shape[1]

    # We need an init of the GMM class
    # to call instance methods
    gmm = GMM(numberOfGaussiansPerClass, data.shape[1], True,
              means, variances, mixtureWeights)

    # These will store the max values and indices
    maxValues = np.zeros(biasCorrectedImageBuffers.shape[0:3],
                         dtype=np.float32)
    maxIndices = np.empty(biasCorrectedImageBuffers.shape[0:3],
                          dtype=np.uint16)

    maxIndices[:] = FreeSurferLabels[bg_label]
    print('done')

    # We need the normalizer if we are writing out the normalized images
    # for CAT
    if cat_opts is not None:
        cat_tissue_dict = cat_opts['cat_tissues']
        cat_mask_dict = cat_opts['cat_masks']
        wm_probs = np.zeros(biasCorrectedImageBuffers.shape[0:3],
                            dtype=np.float32)
        gm_probs = np.zeros(biasCorrectedImageBuffers.shape[0:3],
                            dtype=np.float32)
        csf_probs = np.zeros(biasCorrectedImageBuffers.shape[0:3],
                             dtype=np.float32)
        normalizer = np.zeros(biasCorrectedImageBuffers.shape[0:3],
                              dtype=np.float32)

    # The different structures can share their mixtures between classes
    # E.g., thalamus can be half wm and half gm. This needs to be accounted
    # for when looping. We'll first loop over the structures, i.e.,
    # everything that's in the segmentation, so that we also
    # include the correct fractions for each class (if any)
    # NOTE: this function ONLY calculates the labeling and
    # returns the normalizer (if specified). I.e., we don't
    # need to normalize to find the label with highest probability.
    # But we're not going to compute the posteriors here.
    # We're also filling zero values (outside the mask) from the prior
    print('Into the loop')
    for structureNumber in range(numberOfStructures):
        # Rasterize the current structure from the atlas
        # and cast to float from uint16
        nonNormalized = mesh.rasterize_1a(mask.shape, structureNumber)/65535.0
        prior = nonNormalized[mask]

        # Find which classes we need to look at
        # to get the fractions correct
        print(structureNumber)
        classesToLoop = [i for i, val in enumerate(fractionsTable[:, structureNumber] > 1e-10) if val]
        likelihoods = np.zeros(numberOfVoxels, dtype=np.float32)
        for classNumber in classesToLoop:
            # Compute likelihood for this class
            numberOfComponents = numberOfGaussiansPerClass[classNumber]
            fraction = fractionsTable[classNumber, structureNumber]
            for componentNumber in range(numberOfComponents):
                gaussianNumber = sum(numberOfGaussiansPerClass[:classNumber]) + componentNumber
                mean = np.expand_dims(means[gaussianNumber, :], 1)
                variance = variances[gaussianNumber, :, :]
                mixtureWeight = mixtureWeights[gaussianNumber]

                gaussianLikelihood = gmm.getGaussianLikelihoods(data,
                                                                mean,
                                                                variance)
                likelihoods += gaussianLikelihood * mixtureWeight * fraction

        # Now compute the non-normalized posterior
        nonNormalized[mask] = likelihoods*prior

        if cat_opts is not None:
            normalizer += nonNormalized
            if structureNumber in cat_tissue_dict['WM']:
                wm_probs += nonNormalized
            elif structureNumber in cat_tissue_dict['GM']:
                gm_probs += nonNormalized
            elif structureNumber in cat_tissue_dict['CSF']:
                csf_probs += nonNormalized

        if structureNumber == 0:
            # In the first iteration just save the non-normalized values
            # Indices are initialized to zero anyway
            maxValues = nonNormalized
        else:
            # Check whether we have higher values and save indices
            higher_values = nonNormalized > maxValues
            maxValues[higher_values] = nonNormalized[higher_values]
            maxIndices[higher_values] = FreeSurferLabels[structureNumber]
            

    if cat_opts is not None:
        cat_scan_and_masks = _prep_scans_and_masks([wm_probs, gm_probs,
                                                    csf_probs],
                                                   normalizer,
                                                   maxIndices,
                                                   FreeSurferLabels,
                                                   cat_mask_dict)

    if cat_opts is not None:
        return maxIndices, cat_scan_and_masks
    else:
        return maxIndices


def _prep_scans_and_masks(tissue_prob_list, normalizer, maxInds,
                          FreeSurferLabels, mask_dict):

    # Create normalized image
    eps = np.finfo(float).eps
    normalized_scan = np.zeros_like(normalizer)
    normalized_scan += tissue_prob_list[0]/(normalizer + eps)
    normalized_scan += (2.0*tissue_prob_list[1])/(3.0*(normalizer + eps))
    normalized_scan += tissue_prob_list[2]/(3.0*(normalizer + eps))

    normalized_scan = 3.0*normalized_scan
    # Create cerebrum mask
    mask_cereb = np.zeros_like(normalizer)
    for left, right in zip(mask_dict['Left_Cerebrum'], mask_dict['Right_Cerebrum']):
        mask_cereb += 1*(maxInds == FreeSurferLabels[left])
        mask_cereb += 2*(maxInds == FreeSurferLabels[right])

    for left, right in zip(mask_dict['Left_Cerebellum'], mask_dict['Right_Cerebellum']):
        mask_cereb += 3*(maxInds == FreeSurferLabels[left])
        mask_cereb += 4*(maxInds == FreeSurferLabels[right])

    mask_cereb = mask_cereb.astype(np.uint8)
    # Create subcortical mask
    mask_subcortical = np.zeros_like(normalizer)
    for struct in mask_dict['Sub_cortical']:
        mask_subcortical += 1*(maxInds == FreeSurferLabels[struct])

    # Make boolean
    mask_subcortical = mask_subcortical > 0

    # # Create parahippo mask
    # mask_parahippo = np.zeros_like(normalizer)
    # for struct in mask_dict['Parahippo']:
    #     mask_parahippo += 1*(maxInds == FreeSurferLabels[struct])

    # # Make boolean
    # mask_parahippo = mask_parahippo > 0

    # Create left/right mask by interpolating the cereb mask into the background 
    mask = 2*(mask_cereb == 1) + 2*(mask_cereb == 3) + 1*(mask_cereb == 2) + 1*(mask_cereb == 4)
    # cat_vbdist takes in a voxel size as well, need to probably use that too
    _, _, mask_hemi = cat_vbdist(mask)
    mask_hemi = mask_hemi > 1
    mask_hemi = mask_hemi.astype(np.uint8)

    # return [normalized_scan, mask_cereb, mask_subcortical,
    #         mask_parahippo, mask_hemi]
    return [normalized_scan, mask_cereb, mask_subcortical, mask_hemi]

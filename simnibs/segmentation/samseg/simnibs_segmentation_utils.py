from . import gemsbindings as gems
from .ProbabilisticAtlas import ProbabilisticAtlas
from .SamsegUtility import undoLogTransformAndBiasField, writeImage, maskOutBackground, logTransform, readCroppedImages
from .GMM import GMM
import numpy as np
import nibabel as nib


def writeBiasCorrectedImagesAndSegmentation(output_names_bias,
                                            output_name_segmentation,
                                            parameters_and_inputs):

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
    transform = parameters_and_inputs['transform']
    deformation = parameters_and_inputs['deformation']
    deformationAtlasFileName = parameters_and_inputs['deformationAtlas']

    mesh = probabilisticAtlas.getMesh(
            modelSpecifications.atlasFileName,
            transform,
            initialDeformation=deformation,
            initialDeformationMeshCollectionFileName=deformationAtlasFileName)

    fractionsTable = parameters_and_inputs['fractionsTable']
    GMMparameters = parameters_and_inputs['GMMParameters']
    numberOfGaussiansPerClass = parameters_and_inputs['gaussiansPerClass']
    means = GMMparameters['means']
    variances = GMMparameters['variances']
    mixtureWeights = GMMparameters['mixtureWeights']

    structureNumbers = _calculateSegmentationLoop(
                imageBuffers-biasFields,
                mask, fractionsTable, mesh,
                numberOfGaussiansPerClass,
                means, variances, mixtureWeights)

    FreeSurferLabels = np.array(modelSpecifications.FreeSurferLabels,
                                dtype=np.uint16)
    segmentation = FreeSurferLabels[structureNumbers]
    segmentation = segmentation.reshape(imageBuffers.shape[0:3])
    writeImage(output_name_segmentation,
               segmentation, cropping, exampleImage)


def segmentUpsampled(input_bias_corrected,
                     tissues_upsampled, tissue_settings,
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

    upsampledStructureNumbers = _calculateSegmentationLoop(
                imageBuffersUpsampled,
                maskUpsampled, fractionsTable, meshUpsampled,
                numberOfGaussiansPerClass,
                means, variances, mixtureWeights)

    FreeSurferLabels = np.array(modelSpecifications.FreeSurferLabels,
                                dtype=np.uint16)
    segmentation = FreeSurferLabels[upsampledStructureNumbers]
    segmentation = segmentation.reshape(imageBuffersUpsampled.shape[0:3])
    
    simnibs_tissues = tissue_settings['simnibs_tissues']
    segmentation_tissues = tissue_settings['segmentation_tissues']

    tissue_labeling = np.zeros_like(segmentation)
    for t, label in simnibs_tissues.items():
        tissue_labeling[np.isin(segmentation, segmentation_tissues[t])] = label

    # TODO: Writing to disk for checking
    exampleImage = gems.KvlImage(input_bias_corrected[0])

    writeImage(tissues_upsampled, tissue_labeling,
               croppingUpsampled, exampleImage)


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
    transform = parameters_and_inputs['transform']
    deformation = parameters_and_inputs['deformation']
    deformationAtlasFileName = parameters_and_inputs['deformationAtlas']
    nodePositions = probabilisticAtlas.getMesh(
             modelSpecifications.atlasFileName,
             transform,
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
    image = gems.KvlImage(imageFileNames[0])
    imageToWorldTransformMatrix = image.transform_matrix.as_numpy_array
    image_buffer = image.getImageBuffer()
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
                modelSpecifications.atlasFileName)

    # Get node positions in template voxel space
    nodePositionsTemplate = mesh.points

    # Rasterize the coordinate values
    coordmapTemplate = mesh.rasterize_values(template_buffer.shape,
                                             nodePositionsInWorldSpace)
    # Write the warp file
    warp_image = nib.Nifti1Image(coordmapTemplate,
                                 templateToWorldTransformMatrix)
    nib.save(warp_image, warp_to_mni)

    # Now do it the other way, i.e., template->image
    nodePositionsTemplateWorldSpace = (templateToWorldTransformMatrix @
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
        uncroppedMap[cropping] = coordmapImage[:, :, :, c]
        uncroppedWarp[:, :, :, c] = uncroppedMap

    # Write the warp
    warp_image = nib.Nifti1Image(coordmapImage,
                                 imageToWorldTransformMatrix)
    nib.save(warp_image, warp_from_mni)


def _calculateSegmentationLoop(biasCorrectedImageBuffers,
                               mask, fractionsTable, mesh,
                               numberOfGaussiansPerClass,
                               means, variances,
                               mixtureWeights,
                               computeNormalizer=False):

    data = biasCorrectedImageBuffers[mask, :]
    numberOfVoxels = data.shape[0]
    numberOfStructures = fractionsTable.shape[1]
    # We need an init of the GMM class
    # to call instance methods
    gmm = GMM(numberOfGaussiansPerClass, data.shape[1], True,
              means, variances, mixtureWeights)

    # These will store the max values and indices
    maxValues = np.zeros(biasCorrectedImageBuffers.shape[0:3],
                         dtype=np.float64)
    maxIndices = np.zeros(biasCorrectedImageBuffers.shape[0:3],
                          dtype=np.uint16)

    if computeNormalizer:
        normalizer = np.zeros((numberOfVoxels, 1), dtype=np.float64)

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
    for structureNumber in range(numberOfStructures):
        # Rasterize the current structure from the atlas
        # and cast to float from uint8
        nonNormalized = mesh.rasterize_1a(mask.shape, structureNumber)
        prior = nonNormalized[mask]

        # Find which classes we need to look at
        # to get the fractions correct
        print(structureNumber)
        classesToLoop = [i for i, val in enumerate(fractionsTable[:, structureNumber] > 1e-10) if val]
        likelihoods = np.zeros(numberOfVoxels)
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

        if computeNormalizer:
            normalizer += nonNormalized

        if structureNumber == 0:
            # In the first iteration just save the non-normalized values
            # Indices are intialized to zero anyway
            maxValues = nonNormalized
        else:
            # Check whether we have higher values and save indices
            higher_values = nonNormalized > maxValues
            maxValues[higher_values] = nonNormalized[higher_values]
            maxIndices[higher_values] = structureNumber

    if computeNormalizer:
        return maxIndices.flatten(), normalizer
    else:
        return maxIndices.flatten()

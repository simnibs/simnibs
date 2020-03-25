import numpy as np
import os
import pickle
from .Samseg import Samseg
from .SamsegUtility import *
import logging
from . import gems
import nibabel as nib
eps = np.finfo(float).eps


class SamsegWholeHead(Samseg):

#    @profile
    def get_segmentation_and_bias_fields(self):
        mesh = self.probabilisticAtlas.getMesh(self.modelSpecifications.atlasFileName,
                                               self.transform,
                                               initialDeformation=self.deformation,
                                               initialDeformationMeshCollectionFileName=self.deformationAtlasFileName)

        structureNumbers = self.calculateSegmentationLoop(self.imageBuffers,
                                                          self.mask,
                                                          mesh)
        # Make sure basis functions are not downsampled
        self.biasField.downSampleBasisFunctions([1, 1, 1])
        biasFields = self.biasField.getBiasFields()

        return structureNumbers, biasFields

    def writeBiasCorrectedImagesAndSegmentation(self, output_names_bias,
                                                output_name_segmentation):
        #if 1:
        structureNumbers1, biasFields = self.get_segmentation_and_bias_fields()
        #else:
        posteriors, biasFields, nodePositions, _, _ = self.segment()
        structureNumbers2 = np.array(np.argmax(posteriors, 1),
                                    dtype=np.uint32)

        expImageBuffers, expBiasFields = undoLogTransformAndBiasField(self.imageBuffers,
                                                                      biasFields,
                                                                      self.mask)

        # Write out bias corrected images
        exampleImage = gems.KvlImage(self.imageFileNames[0])
        for contrastNumber, out_name in enumerate(output_names_bias):
            # Bias field corrected image
            writeImage(out_name, expImageBuffers[..., contrastNumber],
                       self.cropping, exampleImage)

        segmentation = np.zeros(self.imageBuffers.shape[0:3], dtype=np.uint16)
        FreeSurferLabels = np.array(self.modelSpecifications.FreeSurferLabels,
                                    dtype=np.uint16)
        segmentation[self.mask] = FreeSurferLabels[structureNumbers2]
        writeImage(output_name_segmentation, segmentation,
                   self.cropping, exampleImage)

        # Write out test labeling
        segmentation = FreeSurferLabels[structureNumbers1]
        segmentation.reshape((self.imageBuffers.shape[0:3]))
        writeImage(os.path.join(self.savePath,'test_labeling.nii.gz'),
                   segmentation, self.cropping, exampleImage)

    def getOptimizationSummary(self):
        return self.optimizationSummary

    def segmentUpsampled(self, input_bias_corrected, tissue_settings):
        # Read the input images doing the necessary cropping
        imageBuffersUpsampled, transformUpsampled, voxelSpacingUpsampled, croppingUpsampled = readCroppedImages(
                                                                                              input_bias_corrected,
                                                                                              self.transformedTemplateFileName)
        # Redo background masking now with the upsampled scans, note this only rasterizes
        # a single class so should be decent memory-wise
        imageBuffersUpsampled, maskUpsampled = maskOutBackground(imageBuffersUpsampled,
                                                         self.modelSpecifications.atlasFileName,
                                                         transformUpsampled,
                                                         self.modelSpecifications.brainMaskingSmoothingSigma,
                                                         self.modelSpecifications.brainMaskingThreshold,
                                                         self.probabilisticAtlas)

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
        meshUpsampled = self.probabilisticAtlas.getMesh(self.modelSpecifications.atlasFileName, 
                                                        transformUpsampled,
                                                        initialDeformation=self.deformation,
                                                        initialDeformationMeshCollectionFileName=self.deformationAtlasFileName)
        upsampledStructureNumbers = self.calculateSegmentationLoop(imageBuffersUpsampled,
                                                               maskUpsampled,
                                                               meshUpsampled)
        
        
        FreeSurferLabels = np.array(self.modelSpecifications.FreeSurferLabels,
                                    dtype=np.uint16)
        segmentation = FreeSurferLabels[upsampledStructureNumbers]
        segmentation.reshape((imageBuffers.shape[0:3]))
        
        simnibs_tissues = tissue_settings['simnibs_tissues']
        segmentation_tissues = tissue_settings['segmentation_tissues']
        
        tissue_labeling = np.zeros_like(segmentation)
        for t, label in simnibs_tissues.items():
            tissue_labeling[np.isin(segmentation, segmentation_tissues[t])] = label
            
        # TODO: Writing to disk for checking
        exampleImage = gems.KvlImage(input_bias_corrected[0])
        writeImage(os.path.join(self.savePath,'test_labeling_upsampled.nii.gz'),
                   segmentation, croppingUpsampled, exampleImage)
        
        

    def saveMesh(self, output_name):
        mesh = self.probabilisticAtlas.getMesh(self.modelSpecifications.atlasFileName,
                                               self.transform,
                                               initialDeformation=self.deformation,
                                               initialDeformationMeshCollectionFileName=self.deformationAtlasFileName)

        estimatedNodePositions = self.probabilisticAtlas.mapPositionsFromSubjectToTemplateSpace(mesh.points,
                                                                                                self.transform)
        self.probabilisticAtlas.saveDeformedAtlas(self.modelSpecifications.atlasFileName,
                                                  output_name,
                                                  estimatedNodePositions)

    def saveHistory(self, output_path):
        history = {'input': {
                'imageFileNames': self.imageFileNames,
                'transformedTemplateFileName': self.transformedTemplateFileName,
                'modelSpecifications': self.modelSpecifications,
                'optimizationOptions': self.optimizationOptions,
                'savePath': self.savePath},
                'imageBuffers': self.imageBuffers,
                'mask': self.mask,
                'historyWithinEachMultiResolutionLevel': self.optimizationHistory,
                "labels": self.modelSpecifications.FreeSurferLabels,
                "names": self.modelSpecifications.names,
                "optimizationSummary": self.optimizationSummary}

        with open(os.path.join(output_path, 'history.p'), 'wb') as file:
            pickle.dump(history, file, protocol=pickle.HIGHEST_PROTOCOL)

    def saveWarpFieldSimNIBS(self, template_name, warp_to_mni, warp_from_mni):
        # Save warp field two ways: the subject space world coordinates
        # in template space, i.e., the iamge voxel coordinates in
        # physical space for every voxel in the template.
        # And the other way around, i.e., the template voxel coordinates
        # in physical space for every voxel in the image.

        # First write image -> template.
        # Get the node positions in image voxels
        nodePositions = self.probabilisticAtlas.getMesh(
             self.modelSpecifications.atlasFileName,
             self.transform,
             initialDeformation=self.deformation,
             initialDeformationMeshCollectionFileName=self.deformationAtlasFileName
         ).points

        # The image is cropped as well so the voxel coordinates
        # do not exactly match with the original image,
        # i.e., there's a shift. Let's undo that.
        nodePositions += [slc.start for slc in self.cropping]

        # Get mapping from voxels to world space of the image.
        image = gems.KvlImage(self.imageFileNames[0])
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
        mesh = self.probabilisticAtlas.getMesh(
                self.modelSpecifications.atlasFileName)

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
        mesh = self.probabilisticAtlas.getMesh(
             self.modelSpecifications.atlasFileName,
             self.transform,
             initialDeformation=self.deformation,
             initialDeformationMeshCollectionFileName=self.deformationAtlasFileName)

        coordmapImage = mesh.rasterize_values(self.imageBuffers.shape[0:-1],
                                              nodePositionsTemplateWorldSpace)
        # The image buffer is cropped so need to set
        # everything to the correct place
        uncroppedWarp = np.zeros(image_buffer.shape+(3,),
                                 dtype=np.float32, order='F')

        for c in range(coordmapImage.shape[-1]):
            uncroppedMap = np.zeros(image_buffer.shape,
                                    dtype=np.float32, order='F')
            uncroppedMap[self.cropping] = coordmapImage[:, :, :, c]
            uncroppedWarp[:, :, :, c] = uncroppedMap

        # Write the warp
        warp_image = nib.Nifti1Image(coordmapImage,
                                     imageToWorldTransformMatrix)
        nib.save(warp_image, warp_from_mni)

    def calculateSegmentationLoop(self, imageBuffers,
                                  mask, mesh,
                                  computeNormalizer=False):

        data = imageBuffers[mask, :]
        numberOfVoxels = data.shape[0]
        fractionsTable = self.classFractions
        numberOfStructures = fractionsTable.shape[1]

        # These will store the max values and indices
        maxValues = np.zeros(imageBuffers.shape[0:3], dtype=np.float64)
        maxIndices = np.zeros(imageBuffers.shape[0:3], dtype=np.uint16)

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
            nonNormalized = mesh.rasterize_1a(mask.shape, structureNumber)/65535.0
            prior = nonNormalized[mask]
            # Find which classes we need to look at
            # to get the fractions correct
            print(structureNumber)
            classesToLoop = [i for i, val in enumerate(fractionsTable[:, structureNumber] > 1e-10) if val] 
            likelihoods = np.zeros(numberOfVoxels)
            for classNumber in classesToLoop:
                # Compute likelihood for this class
                numberOfComponents = self.gmm.numberOfGaussiansPerClass[classNumber]
                fraction = fractionsTable[classNumber, structureNumber]
                for componentNumber in range(numberOfComponents):
                    gaussianNumber = sum(self.gmm.numberOfGaussiansPerClass[:classNumber]) + componentNumber
                    mean = np.expand_dims(self.gmm.means[gaussianNumber, :], 1)
                    variance = self.gmm.variances[gaussianNumber, :, :]
                    mixtureWeight = self.gmm.mixtureWeights[gaussianNumber]

                    gaussianLikelihood = self.gmm.getGaussianLikelihoods(data, 
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

        # Okay now we have labeling, still need to undo cropping and map to
        # label values
        # segmentation = np.zeros(imageBuffers.shape[0:3], dtype=np.uint16)
        # FreeSurferLabels = np.array(self.modelSpecifications.FreeSurferLabels,
        #                           dtype=np.uint16)
        # segmentation[mask] = maxIndices

        if computeNormalizer:
            return maxIndices.flatten(), normalizer
        else:
            return maxIndices.flatten()

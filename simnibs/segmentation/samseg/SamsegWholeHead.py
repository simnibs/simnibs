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
    
    
    def writeBiasCorrectedImagesAndSegmentation(self, output_names_bias,output_name_segmentation):
        
        posteriors, biasFields, nodePositions, _, _ = self.segment()
        
        expImageBuffers, expBiasFields = undoLogTransformAndBiasField(self.imageBuffers, biasFields, self.mask)
        
        # Write out bias corrected images
        exampleImage = gems.KvlImage(self.imageFileNames[0])
        for contrastNumber, out_name in enumerate(output_names_bias):
            # Bias field corrected image
            writeImage(out_name, expImageBuffers[..., contrastNumber],
                   self.cropping, exampleImage)
            
        # Write out segmentations
        # Convert into a crisp, winner-take-all segmentation, labeled according to the FreeSurfer labeling/naming convention
        if self.threshold is not None:
            # Figure out the structure number of the special snowflake structure
            for structureNumber, name in enumerate(self.modelSpecifications.names):
                if self.thresholdSearchString in name:
                    thresholdStructureNumber = structureNumber
                    break

            # Threshold
            print('thresholding posterior of ', self.modelSpecifications.names[thresholdStructureNumber], 'with threshold:', self.threshold)
            tmp = posteriors[:, thresholdStructureNumber].copy()
            posteriors[:, thresholdStructureNumber] = posteriors[:, thresholdStructureNumber] > self.threshold

            # Majority voting
            structureNumbers = np.array(np.argmax(posteriors, 1), dtype=np.uint32)

            # Undo thresholding in posteriors
            posteriors[:, thresholdStructureNumber] = tmp

        else:
            # Majority voting
            structureNumbers = np.array(np.argmax(posteriors, 1), dtype=np.uint32)

        segmentation = np.zeros(self.imageBuffers.shape[0:3], dtype=np.uint16)
        FreeSurferLabels = np.array(self.modelSpecifications.FreeSurferLabels, dtype=np.uint16)
        segmentation[self.mask] = FreeSurferLabels[structureNumbers]
        writeImage(output_name_segmentation, segmentation, self.cropping,exampleImage)

            
    def getOptimizationSummary(self):
        return self.optimizationSummary
    
    def segmentUpsampled(self, input_bias_corrected):
        # Read the input images doing the necessary cropping
        imageBuffersUpsampled, transformUpsampled, voxelSpacingUpsampled, croppingUpsampled = readCroppedImages(
                                                                                              input_bias_corrected,
                                                                                              self.transformedTemplateFileName)
        #Redo background masking now with the upsampled scans, note this only rasterizes
        #a single class so should be decent memory-wise
        imageBuffersUpsampled, maskUpsampled = maskOutBackground(imageBuffersUpsampled,
                                                         self.modelSpecifications.atlasFileName,
                                                         transformUpsampled,
                                                         self.modelSpecifications.brainMaskingSmoothingSigma,
                                                         self.modelSpecifications.brainMaskingThreshold,
                                                         self.probabilisticAtlas)
        
        #Log-transform the intensities, note the scans are already bias corrected so no need to remove the bias contribution
        imageBuffersUpsampled = logTransform(imageBuffersUpsampled, maskUpsampled)
        
        #Calculate the posteriors.
        #NOTE: the method for calculating the posteriors in the gmm class rasterizes
        #all classes in the atlas at once. This is very memory-heavy if we have many classes
        #and the input resolution is high. Here we do it instead in a loop to save memory.
        mesh = self.probabilisticAtlas.getMesh(self.modelSpecifications.atlasFileName, transformUpsampled,
                                               initialDeformation=self.deformation,
                                               initialDeformationMeshCollectionFileName=self.deformationAtlasFileName)
        
        self.calculatePosteriorsInLoop(upsampledImageBuffers, upsampledMask)
        
        
        
        #Do the segmentation
        
        #Write it out


    def postProcess(self):
        # =======================================================================================
        #
        # Segment the data using the estimate model parameters, and write results out
        #
        # =======================================================================================

        # OK, now that all the parameters have been estimated, try to segment the original, full resolution image
        # with all the original labels (instead of the reduced "super"-structure labels we created)
        posteriors, biasFields, nodePositions, _, _ = self.segment()

        # Write out segmentation and bias field corrected volumes
        volumesInCubicMm = writeResults(self.imageFileNames, self.savePath, self.imageBuffers, self.mask,
                                                     biasFields,
                                                     posteriors, self.modelSpecifications.FreeSurferLabels,
                                                     self.cropping,
                                                     self.targetIntensity, self.targetSearchStrings,
                                                     self.modelSpecifications.names,
                                                     self.threshold, self.thresholdSearchString,
                                                     savePosteriors=self.savePosteriors)

        logger = logging.getLogger(__name__)
            
            

        # Save the final mesh collection
        if self.saveMesh:
            logger.info('Saving the final mesh in template space')
            image_base_path, _ = os.path.splitext(self.imageFileNames[0])
            _, scanName = os.path.split(image_base_path)
            deformedAtlasFileName = os.path.join(self.savePath, scanName + '_meshCollection.txt.gz')
            self.probabilisticAtlas.saveDeformedAtlas(self.modelSpecifications.atlasFileName, deformedAtlasFileName,
                                                      nodePositions)

        # Save the history of the parameter estimation process
        if self.saveHistory:
            history = {'input': {
                'imageFileNames': self.imageFileNames,
                'transformedTemplateFileName': self.transformedTemplateFileName,
                'modelSpecifications': self.modelSpecifications,
                'optimizationOptions': self.optimizationOptions,
                'savePath': self.savePath
            }, 'imageBuffers': self.imageBuffers, 'mask': self.mask,
                'historyWithinEachMultiResolutionLevel': self.optimizationHistory,
                "labels": self.modelSpecifications.FreeSurferLabels, "names": self.modelSpecifications.names,
                "volumesInCubicMm": volumesInCubicMm, "optimizationSummary": self.optimizationSummary}
            with open(os.path.join(self.savePath, 'history.p'), 'wb') as file:
                pickle.dump(history, file, protocol=pickle.HIGHEST_PROTOCOL)

        return self.modelSpecifications.FreeSurferLabels, self.modelSpecifications.names, volumesInCubicMm, self.optimizationSummary
    
    
    def saveWarpFieldSimNIBS(self, filename):
        #Save warp field two ways: the subject space world coordinates in template space,
        #i.e., the iamge voxel coordinates in physical space for every voxel in the template.
        #And the other way around, i.e., the template voxel coordinates in physical space for
        #every voxel in the image.
        
        # First write image -> template.
        #Get the node positions in image voxels
        nodePositions = self.probabilisticAtlas.getMesh(
             self.modelSpecifications.atlasFileName,
             self.transform,
             initialDeformation=self.deformation,
             initialDeformationMeshCollectionFileName=self.deformationAtlasFileName
         ).points
        
        #The image is cropped as well so the voxel coordinates do not exactly match with the 
        #original image, i.e., there's a shift. Let's undo that.
        nodePositions += [slc.start for slc in self.cropping]
        
         # Get mapping from voxels to world space of the image.
        image = gems.KvlImage(self.imageFileNames[0])
        imageToWorldTransformMatrix = image.transform_matrix.as_numpy_array
        image_buffer = image.getImageBuffer()
        #Transform the node positions
        nodePositionsInWorldSpace = (imageToWorldTransformMatrix @ 
                                     np.pad(nodePositions, ((0, 0), (0, 1)), 
                                     'constant', constant_values=1).T).T
        nodePositionsInWorldSpace = nodePositionsInWorldSpace[:, 0:3]
        
        #TODO: Hard-coded here need to change
        atlas_path = os.path.join('C:/','Users','oupu','Documents','test_data','Atlas1_new')
        template_name = os.path.join(atlas_path,'template.nii')
        template = gems.KvlImage(template_name)
        templateToWorldTransformMatrix = template.transform_matrix.as_numpy_array
        template_buffer = template.getImageBuffer()
        
        # rasterize the final node coordinates (in image space) using the initial template mesh
        mesh = self.probabilisticAtlas.getMesh(self.modelSpecifications.atlasFileName)
        #Get node positions in template voxel space
        nodePositionsTemplate = mesh.points
        
        #Rasterize the coordinate values
        coordmapTemplate = mesh.rasterize_values(template_buffer.shape, nodePositionsInWorldSpace)
        # write the warp file
        warp_image = nib.Nifti1Image(coordmapTemplate,templateToWorldTransformMatrix)
        nib.save(warp_image,os.path.join(self.savePath,'image_to_mni_warp.nii.gz'))
        
        #Now do it the other way, i.e., template->image
        nodePositionsTemplateWorldSpace = (templateToWorldTransformMatrix @ 
                                     np.pad(nodePositionsTemplate, ((0, 0), (0, 1)), 
                                     'constant', constant_values=1).T).T
        nodePositionsTemplateWorldSpace = nodePositionsTemplateWorldSpace[:, 0:3]
        
        #Okay get the mesh in image space
        mesh = self.probabilisticAtlas.getMesh(
             self.modelSpecifications.atlasFileName,
             self.transform,
             initialDeformation=self.deformation,
             initialDeformationMeshCollectionFileName=self.deformationAtlasFileName)
        
        coordmapImage = mesh.rasterize_values(self.imageBuffers.shape[0:-1],nodePositionsTemplateWorldSpace)
        #The image buffer is cropped so need to set everything to the correct place
        uncroppedWarp = np.zeros(image_buffer.shape+(3,), dtype=np.float32, order='F')
        for c in range(coordmapImage.shape[-1]):
            uncroppedMap = np.zeros(image_buffer.shape, dtype=np.float32, order='F')
            uncroppedMap[self.cropping] = coordmapImage[:,:,:,c]
            uncroppedWarp[:,:,:,c] = uncroppedMap
            
        #Write the warp
        warp_image = nib.Nifti1Image(coordmapImage,imageToWorldTransformMatrix)
        nib.save(warp_image,os.path.join(self.savePath,'mni_to_image_warp.nii.gz'))
        
        def calculatePosteriorsInLoop(self,upsampledImageBuffers, upsampledMask, upsampledMesh):
            upsampledData = upsampledImageBuffers[upsampledMask, :]
            numberOfVoxels = upsampledData.shape[0]
            numberOfStructures = self.classFractions.shape[1]
            
            likelihoods = np.zeros((numberOfVoxels,1), dtype=np.float64)
            normalizer = np.zeros((numberOfVoxels,1), dtype=np.float64)
            
            #First calculate the normalizer in a loop
            for classNumber in range(self.numberOfClasses):
                # Compute likelihood for this class
                classLikelihoods = np.zeros(numberOfVoxels)
                numberOfComponents = self.numberOfGaussiansPerClass[classNumber]
                for componentNumber in range(numberOfComponents):
                    gaussianNumber = sum(self.numberOfGaussiansPerClass[:classNumber]) + componentNumber
                    mean = np.expand_dims(self.means[gaussianNumber, :], 1)
                    variance = self.variances[gaussianNumber, :, :]
                    mixtureWeight = self.mixtureWeights[gaussianNumber]

                    gaussianLikelihoods = self.getGaussianLikelihoods(data, mean, variance)
                    classLikelihoods += gaussianLikelihoods * mixtureWeight

                # Add contribution to the actual structures
                for structureNumber in range(numberOfStructures):
                    fraction = fractionsTable[classNumber, structureNumber]
                    if fraction < 1e-10:
                        continue
                    likelihoods[:, structureNumber] += classLikelihoods * fraction
                    
        

        
        
        

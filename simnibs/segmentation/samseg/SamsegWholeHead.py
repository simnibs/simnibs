import numpy as np
import os
import pickle
from .Samseg import Samseg
from .SamsegUtility import writeImage
import charm_gems as gems
eps = np.finfo(float).eps
from .simnibs_segmentation_utils import _calculateSegmentationLoop

class SamsegWholeHead(Samseg):

    def standard_segmentation(self):
        posteriors, _, _, _, _ = self.segment()
        structureNumbers = np.array(np.argmax(posteriors, 1), dtype=np.uint32)

        freeSurferSegmentation = np.zeros(self.imageBuffers.shape[0:3],
                                          dtype=np.uint16)
        FreeSurferLabels = np.array(self.modelSpecifications.FreeSurferLabels,
                                    dtype=np.uint16)
        freeSurferSegmentation[self.mask] = FreeSurferLabels[structureNumbers]

        # Write out various images - segmentation first
        exampleImage = gems.KvlImage(self.imageFileNames[0])
        writeImage(os.path.join(self.savePath, 'crispSegmentation.nii.gz'),
                   freeSurferSegmentation,
                   self.cropping,
                   exampleImage)

    def getOptimizationSummary(self):
        return self.optimizationSummary

    def writeMesh(self, output_name):
        mesh = self.probabilisticAtlas.getMesh(
                self.modelSpecifications.atlasFileName,
                self.transform,
                initialDeformation=self.deformation,
                initialDeformationMeshCollectionFileName=self.deformationAtlasFileName)

        estimatedNodePositions = self.probabilisticAtlas.mapPositionsFromSubjectToTemplateSpace(mesh.points,
                                                                                                self.transform)
        self.probabilisticAtlas.saveDeformedAtlas(
                self.modelSpecifications.atlasFileName,
                output_name,
                estimatedNodePositions)

    def saveHistory(self, output_path):
        history = {'input': {
                'imageFileNames': self.imageFileNames,
                'transformedTemplateFileName':
                self.transformedTemplateFileName,
                'modelSpecifications': self.modelSpecifications,
                'optimizationOptions': self.optimizationOptions,
                'savePath': self.savePath},
                'imageBuffers': self.imageBuffers,
                'mask': self.mask,
                'historyWithinEachMultiResolutionLevel':
                self.optimizationHistory,
                "labels": self.modelSpecifications.FreeSurferLabels,
                "names": self.modelSpecifications.names,
                "optimizationSummary": self.optimizationSummary.append}

        with open(os.path.join(output_path, 'history.p'), 'wb') as file:
            pickle.dump(history, file, protocol=pickle.HIGHEST_PROTOCOL)

    def saveParametersAndInput(self, output_path=None):
        # Make sure the bias field basis functions are not downsampled
        # This might happen if the downsampling factor in the last
        # resolution level is >1
        self.biasField.downSampleBasisFunctions([1, 1, 1])
        parameters_and_inputs = {'GMMParameters': {'mixtureWeights':
                                                   self.gmm.mixtureWeights,
                                                   'means': self.gmm.means,
                                                   'variances':
                                                   self.gmm.variances},
                                 'fractionsTable': self.classFractions,
                                 'gaussiansPerClass':
                                 self.gmm.numberOfGaussiansPerClass,
                                 'deformation': self.deformation,
                                 'modelSpecifications':
                                 self.modelSpecifications,
                                 'imageFileNames': self.imageFileNames,
                                 'deformationAtlas':
                                 self.deformationAtlasFileName,
                                 'imageBuffers': self.imageBuffers,
                                 'biasFields': self.biasField.getBiasFields(),
                                 'mask': self.mask,
                                 'cropping': self.cropping,
                                 'transform': self.transform.as_numpy_array,
                                 'names': self.modelSpecifications.names}

        if output_path is not None:
            with open(os.path.join(output_path, 'parameters.p'), 'wb') as file:
                pickle.dump(parameters_and_inputs,
                            file, protocol=pickle.HIGHEST_PROTOCOL)

        return parameters_and_inputs

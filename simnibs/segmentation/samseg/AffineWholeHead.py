import os
import scipy.io
import scipy.ndimage
import numpy as np
from scipy.optimize import minimize_scalar
import logging

#import freesurfer as fs
from . import gems

from .utilities import requireNumpyArray
from .figures import initVisualizer
from .SamsegUtility import readCroppedImages
eps = np.finfo(float).eps

class AffineWholeHead:
    def __init__(self,
                 scalingFactors=[[0.8,0.8,0.85],[0.9,0.9,0.9],[0.95,0.95,0.9]], 
                 thetas=[np.pi / 180 * theta for theta in [-7,-3.5,0,3.5,7]],
                 K=1e-7,
                 targetDownsampledVoxelSpacing=3.0,
                 maximalDeformationStopCriterion=0.005):
        self.scalingFactors = scalingFactors
        self.thetas = thetas
        self.targetDownsampledVoxelSpacing = targetDownsampledVoxelSpacing
        self.maximalDeformationStopCriterion = maximalDeformationStopCriterion

    def registerAtlas(self,
            imageFileName,
            meshCollectionFileName,
            templateFileName,
            savePath,
            template_coregistered_name,
            visualizer=None,
            worldToWorldTransformMatrix=None,
            initTransform=None,
            K=1e-7,
        ):
        
       
        # ------ Setup ------
        logger = logging.getLogger(__name__)
        
        # Read in image and template, as well as their coordinates in world (mm) space
        image = gems.KvlImage(imageFileName)
        imageToWorldTransformMatrix = image.transform_matrix.as_numpy_array
        template = gems.KvlImage(templateFileName)
        templateImageToWorldTransformMatrix = template.transform_matrix.as_numpy_array
        

        # Setup null visualization if necessary
        if visualizer is None:
            visualizer = initVisualizer(False, False)

        # ------ Register Image ------

        if worldToWorldTransformMatrix is not None:
            # The world-to-world transfrom is externally given, so let's just compute the corresponding image-to-image
            # transform (needed for subsequent computations) and be done
            logger.info('world-to-world transform supplied - skipping registration')
            imageToImageTransformMatrix = np.linalg.inv(imageToWorldTransformMatrix) @ worldToWorldTransformMatrix @ templateImageToWorldTransformMatrix
            optimizationSummary = None
        else:
            # The solution is not externally (secretly) given, so we need to compute it.
            logger.info('performing affine atlas registration')
            logger.info('image: %s' % imageFileName)
            logger.info('template: %s' % templateFileName)

            lineSearchMaximalDeformationIntervalStopCriterion = self.maximalDeformationStopCriterion  # Doesn't seem to matter very much

            # Initialization
            if initTransform is None:
                initialWorldToWorldTransformMatrix = np.identity(4)
            else:
                # Assume the initialization matrix is LPS2LPS
                logger.info('initializing with predifined transform')
                initialWorldToWorldTransformMatrix = initTransform

            #Calculate initial template-to-image transform
            multiplied = (initialWorldToWorldTransformMatrix @ templateImageToWorldTransformMatrix)
            initialImageToImageTransformMatrix = np.linalg.solve(imageToWorldTransformMatrix, multiplied)
            # Figure out how much to downsample (depends on voxel size)
            voxelSpacing = np.sum(imageToWorldTransformMatrix[0:3, 0:3] ** 2, axis=0) ** (1 / 2)
            downSamplingFactors = np.round(self.targetDownsampledVoxelSpacing / voxelSpacing)
            downSamplingFactors[downSamplingFactors < 1] = 1

            # Use initial transform to define the reference (rest) position of the mesh (i.e. the one
            # where the log-prior term is zero)
            mesh_collection = gems.KvlMeshCollection()
            mesh_collection.read(meshCollectionFileName)
            mesh_collection.k = K*np.prod(downSamplingFactors)
            mesh_collection.transform(gems.KvlTransform(requireNumpyArray(initialImageToImageTransformMatrix)))
            mesh = mesh_collection.reference_mesh
            
            
            # Get image data
            imageBuffer = image.getImageBuffer()
            visualizer.show(images=imageBuffer, window_id='atlas initial', title='Initial Atlas Registration')

            # Downsample
            imageBuffer = imageBuffer[
                          ::int(downSamplingFactors[0]),
                          ::int(downSamplingFactors[1]),
                          ::int(downSamplingFactors[2])]
            image = gems.KvlImage(requireNumpyArray(imageBuffer))
            
            mesh.scale(1/downSamplingFactors)
            visualizer.show(mesh=mesh, shape=imageBuffer.shape, window_id='atlas mesh', title="Atlas Mesh")
            # Get a registration cost and use it to evaluate some promising starting point proposals
            calculator = gems.KvlCostAndGradientCalculator('MutualInformation', [image], 'Affine')
            translation = self.initializeCenterOfMass(calculator,mesh,imageBuffer,initialImageToImageTransformMatrix,visualizer)
            
            
            # Initialize scales and rotation, we are trying different scales and
            # rotations here, so need to be careful that the initial positions
            # of the mesh nodes in the image define the reference position.
            # That's why let's first undo the downsampling and initial transformation
            mesh.scale(downSamplingFactors)
            mesh_collection.transform(gems.KvlTransform(requireNumpyArray(np.linalg.inv(initialImageToImageTransformMatrix))))
            
            # Remove the neck from the atlas. The reason for this is that the neck often has a non-linear
            # deformation, which throws off the affine registration.
            # The neck will be handled separately after the affine.
            alphasNeckRelabeled = self.removeNeck(mesh)
            
            initialImageToImageTransformMatrix, bestScales = self.findRotationAndScale(
                                                            mesh_collection,downSamplingFactors,
                                                            K,alphasNeckRelabeled,translation,
                                                            imageToWorldTransformMatrix,
                                                            templateImageToWorldTransformMatrix,
                                                            calculator,imageBuffer,visualizer)
            
            
            #Update image-to-image transformation matrix
            K_new = K/np.prod(bestScales)
            mesh_collection.k = K_new*np.prod(downSamplingFactors)
            mesh_collection.transform(gems.KvlTransform(requireNumpyArray(initialImageToImageTransformMatrix)))
            mesh = mesh_collection.reference_mesh
            mesh.scale(1/downSamplingFactors)
            originalNodePositions = mesh.points
            
            
            visualizer.show(mesh=mesh, images=imageBuffer, window_id='atlas mesh rotation', title="Atlas Mesh after rotation")
            # Get an optimizer, and stick the cost function into it
            optimization_parameters = {
                'Verbose': 1.0,
                'MaximalDeformationStopCriterion': self.maximalDeformationStopCriterion,
                'LineSearchMaximalDeformationIntervalStopCriterion': lineSearchMaximalDeformationIntervalStopCriterion,
                'BFGS-MaximumMemoryLength': 12.0  # Affine registration only has 12 DOF
            }
            optimizer = gems.KvlOptimizer( 'L-BFGS', mesh, calculator, optimization_parameters)

            numberOfIterations = 0
            minLogLikelihoodTimesPriors = []
            maximalDeformations = []
            visualizer.start_movie(window_id='atlas iteration', title='Atlas Registration - the movie')
            while True:
                minLogLikelihoodTimesPrior, maximalDeformation = optimizer.step_optimizer_atlas()
                minLogLikelihoodTimesPriors.append(minLogLikelihoodTimesPrior)
                maximalDeformations.append(maximalDeformation)

                if maximalDeformation == 0:
                    break
                numberOfIterations += 1
                visualizer.show(mesh=mesh, images=imageBuffer, window_id='atlas iteration', title='Atlas Registration')

            visualizer.show_movie(window_id='atlas iteration')
            nodePositions = mesh.points
            pointNumbers = [0, 110, 201, 302]
            originalY = np.vstack((np.diag(downSamplingFactors) @ originalNodePositions[pointNumbers].T, [1, 1, 1, 1]))
            Y = np.vstack((np.diag(downSamplingFactors) @ nodePositions[pointNumbers].T, [1, 1, 1, 1]))
            extraImageToImageTransformMatrix = Y @ np.linalg.inv(originalY)

            # Final result: the image-to-image (from template to image) as well as the world-to-world transform that
            # we computed (the latter would be the identity matrix if we didn't move the image at all)
            imageToImageTransformMatrix = extraImageToImageTransformMatrix @ initialImageToImageTransformMatrix
            worldToWorldTransformMatrix = imageToWorldTransformMatrix @ imageToImageTransformMatrix @ np.linalg.inv(templateImageToWorldTransformMatrix)

            optimizationSummary = {'numberOfIterations': len(minLogLikelihoodTimesPriors),
                                   'cost': minLogLikelihoodTimesPriors[-1]}

        # ------ Save Registration Results ------
        
        
        # Save the image-to-image and the world-to-world affine registration matrices
        scipy.io.savemat(os.path.join(savePath, 'coregistrationMatrices.mat'),
                         {'imageToImageTransformMatrix': imageToImageTransformMatrix,
                          'worldToWorldTransformMatrix': worldToWorldTransformMatrix } )


        # Save the coregistered template. For historical reasons, we applied the estimated transformation to the template... let's do that now
        desiredTemplateImageToWorldTransformMatrix = np.asfortranarray(imageToWorldTransformMatrix @ imageToImageTransformMatrix)
        template.write(template_coregistered_name, gems.KvlTransform(desiredTemplateImageToWorldTransformMatrix))

        return worldToWorldTransformMatrix, optimizationSummary

    def computeRotationAndScalingMatrixGuessEstimates(self,theta,scales):
        # Rotation around X-axis (direction from left to right ear)
        rotationMatrix = np.identity(4, dtype=np.double)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        rotationMatrix[1, 1] = cos_theta
        rotationMatrix[1, 2] = -sin_theta
        rotationMatrix[2, 1] = sin_theta
        rotationMatrix[2, 2] = cos_theta

        # Isotropic scaling
        scalingMatrix = np.diag([scales[0], scales[1], scales[2], 1.0])

        return rotationMatrix, scalingMatrix
    
    def initializeCenterOfMass(self,calculator,mesh,imageBuffer,initialImageToImageTransformMatrix,visualizer):
        logger = logging.getLogger(__name__)
        centerOfGravityImage = np.array(scipy.ndimage.measurements.center_of_mass(imageBuffer))
        
            
        #Calculate the average position
        nodePositionsInImage = mesh.points
        meanNodePosition = np.mean(nodePositionsInImage, axis=0)
        
        #Calculate translation
        baseTranslation = centerOfGravityImage - meanNodePosition
        
        #Move the nodes,so that the atlas more or less matches the FOV of the image
        mesh.points = nodePositionsInImage + baseTranslation
        priors = mesh.rasterize(imageBuffer.shape[0:3], -1)
        head = np.sum(priors[:,:,:,1:],axis=3)
        head = head/65535
        centerOfGravityAtlas = np.array(scipy.ndimage.measurements.center_of_mass(head))
        
        baseTranslation = baseTranslation +(centerOfGravityImage - centerOfGravityAtlas)
        mesh.points = nodePositionsInImage + baseTranslation
        
        #Check which way we should translate in image space. In mesh space
        #this would be the axial (z) direction
        zInd = np.argmax(np.abs(initialImageToImageTransformMatrix[0:3,2]))
        displacementArray = np.array([0,0,0])
        displacementArray[zInd] = 1

        # Attempt a few different initial alignments. If the initial center-of-gravity placement fails
        # try shifting the translation up and down along the Y axis by 10 voxels
        bestVerticalDisplacement = 0
        bestCost, gradient = calculator.evaluate_mesh_position_a(mesh)
        
        for verticalDisplacement in (0, -15.0/self.targetDownsampledVoxelSpacing, 15.0 / self.targetDownsampledVoxelSpacing):
            translation = baseTranslation + displacementArray*verticalDisplacement
            mesh.points = nodePositionsInImage + translation
            trialCost, trialGradient = calculator.evaluate_mesh_position_b(mesh)
            #visualizer.show(mesh=mesh, images=imageBuffer, window_id='atlas translation'+str(verticalDisplacement), title="Atlas Mesh after translation"+str(verticalDisplacement))
            
            if trialCost < bestCost:
                # This is better starting position - remember that we applied it
                bestVerticalDisplacement = verticalDisplacement
                bestCost = trialCost
                
        logger.info("Best vertical displacement is: "+str(bestVerticalDisplacement))
        translation = baseTranslation + displacementArray*bestVerticalDisplacement
        
        
        mesh.points = nodePositionsInImage + translation
        visualizer.show(mesh=mesh, images=imageBuffer, window_id='atlas mesh translation', title="Atlas Mesh after translation")
        
        mesh.points = nodePositionsInImage
        return translation
    
    def removeNeck(self,mesh):
        alphas = mesh.alphas
        spineAlphas = alphas[:,-1] #Note: the spine class is the last one in the affine atlas. Might need to change this in the future.
        mask = spineAlphas > 0.01
        
        #Get z-coordinates
        zPositions = mesh.points[:,-1]
        
        #Find the max position for spine and cut from there
        spinePosZ = zPositions[mask]
        maxInds = np.argwhere(np.argmax(spinePosZ))
        maxZCoord = zPositions[maxInds[0]]
        zMask = zPositions < maxZCoord
        
        #Create a "neck" class which includes everything in the neck
        alphasWithNeck = np.zeros((alphas.shape[0],alphas.shape[1]+1))
        alphasWithNeck[:,1:] = alphas
        alphasWithNeck[zMask,1:] = 0
        alphasWithNeck[zMask,0] = 1.0
        
        return alphasWithNeck
    
    def findRotationAndScale(self,mesh_collection,downSamplingFactors,K,alphasNeckRelabeled,translation,
                             imageToWorldTransformMatrix,templateImageToWorldTransformMatrix,calculator,
                             imageBuffer,visualizer):
        logger = logging.getLogger(__name__)
        bestTheta = 0
        bestScales = []
        bestCost = np.inf
        bestI2I = np.identity(4)
        
            
        for scales in self.scalingFactors:
            for theta in self.thetas:
                rMat, sMat = self.computeRotationAndScalingMatrixGuessEstimates(theta,scales)
                iniI2ITMP = sMat @ rMat @ templateImageToWorldTransformMatrix
                iniI2ITMP = np.linalg.solve(imageToWorldTransformMatrix,iniI2ITMP)
                iniI2ITMP[0:3,3] = iniI2ITMP[0:3,3] + (np.diag(downSamplingFactors) @ translation)
                
                #Okay now transform the mesh collection and set the K
                K_new = K/(np.prod(scales))
                mesh_collection.k = K_new*np.prod(downSamplingFactors)
                mesh_collection.transform(gems.KvlTransform(requireNumpyArray(iniI2ITMP)))
                mesh = mesh_collection.reference_mesh
                mesh.scale(1/downSamplingFactors)
                
                
                #Set alphas
                mesh.alphas = alphasNeckRelabeled
                trialCost, trialGradient = calculator.evaluate_mesh_position_b(mesh)
                #visualizer.show(mesh=mesh, images=imageBuffer, window_id='atlas mesh rotation'+str(scales[0])+str(theta), title="Atlas Mesh after rotation"+str(scales[0])+str(theta))
                if trialCost < bestCost:
                    bestCost = trialCost
                    bestTheta = theta
                    bestScales = scales
                    bestI2I = iniI2ITMP
                #Remove the transformations to set the mesh in the starting
                #Position for the next try
                mesh.scale(downSamplingFactors)
                mesh_collection.transform(gems.KvlTransform(requireNumpyArray(np.linalg.inv(iniI2ITMP))))
                        
        logger.info('Best scale: '+str(bestScales[0])+','+str(bestScales[1])+','+str(bestScales[2])+', best theta: '+str(bestTheta)+'.')
        return bestI2I, bestScales
        
    def adjust_neck(self, T1, transformed_template_file_name, mesh_level1, mesh_level2, visualizer, downsampling_target=3.0):
        logger = logging.getLogger(__name__)
        image_buffer, transformation_to_image, voxel_spacing, cropping = readCroppedImages([T1],transformed_template_file_name)
        # Figure out how much to downsample (depends on voxel size)
        downsampling_factors = np.round(downsampling_target / voxel_spacing)
        downsampling_factors[downsampling_factors < 1] = 1
        image_buffer = image_buffer[
                          ::int(downsampling_factors[0]),
                          ::int(downsampling_factors[1]),
                          ::int(downsampling_factors[2])]
        #Read in the mesh collection, and set K to be very small so that the intensity cost dominates
        K = 1e-20
        mesh_collection = gems.KvlMeshCollection()
        mesh_collection.read(mesh_level2)
        mesh_collection.k = K*np.prod(downsampling_factors)
        mesh_collection.transform(transformation_to_image)
        mesh = mesh_collection.reference_mesh
        #Downsample if needed
        mesh.scale(1/downsampling_factors)
    
        #Let's define the neck deformation based on the distance 
        alphas = mesh.alphas
        spineAlphas = alphas[:,-1] #Note: the spine class is the last one in the affine atlas. Might need to change this in the future.
        mask_spine = spineAlphas > 0.01
    
        #Let's figure out where the z-coordinate is and which way is up
        transformation_matrix = transformation_to_image.as_numpy_array
        z_dim = np.argmax(np.abs(transformation_matrix[0:3,2]))
        z_direction = transformation_matrix[z_dim,2]
        
        #Get z-coordinates
        z_positions = mesh.points[:,z_dim]
        spine_positions_z = z_positions[mask_spine]
        mask_neck = 0
        neck_pos = 0
        #The position where spine starts (from the brainstem) depends on if the direction is I->S or S->I
        if(z_direction<1): #I->S
            top_ind = np.argwhere(np.argmax(spine_positions_z))
            top_pos = spine_positions_z[top_ind[0]]
            mask_neck = z_positions<top_pos
            neck_pos = z_positions[mask_neck]
            z_dist = top_pos - neck_pos
        else: #S->I
            top_ind = np.argwhere(np.argmin(spine_positions_z))
            top_pos = spine_positions_z[top_ind[0]]
            mask_neck = z_positions>top_pos
            neck_pos = z_positions[mask_neck]
            z_dist = neck_pos - top_pos
        
        #Okay the distance from the top of the spine defines the amount of deformation in the A->P direction
        #Let's first check where the A-P direction is
        y_dim = np.argmax(np.abs(transformation_matrix[0:3,1]))
        deformation_field = np.zeros_like(mesh.points)
        deformation_field[mask_neck,y_dim] = z_dist
    
        #Okay one more trick that seems to work, only consider a subset of the structures to compute the cost
        #These are: air internal, spine, cortical bone and spongy bone
        #NOTE! This need to be read from the setup as otherwise this will fail if the atlas is changed!
        alphas_new = alphas[:,[5,6,8,9]]
        mesh.alphas = alphas_new
    
        #Let's see how it looks
        visualizer.show(mesh=mesh, images=image_buffer, window_id='Initial Neck', title="Initial Neck")
        image = gems.KvlImage(requireNumpyArray(image_buffer))
        calculator = gems.KvlCostAndGradientCalculator('MutualInformation', [image], 'Affine')
    
        #Check initial cost
        initCost, _ = calculator.evaluate_mesh_position_a(mesh)
        logger.info("Initial cost: " + str(initCost))
    
        # Optimize the linear deformation, I'll do this now naively using the Nelder-Mead Simplex algorithm
        # We should be able to compute the gradient for this simple thing, so something smarter should be used
        initial_node_positions = mesh.points
       
        res = minimize_scalar(self._neck_cost,None,None,
                              (initial_node_positions,
                               deformation_field,
                               calculator,mesh,
                               image_buffer,
                               visualizer))
        
        mesh.points = initial_node_positions - res.x*deformation_field
    
        visualizer.show(mesh=mesh, images=image_buffer, window_id='Corrected Neck', title="Corrected Neck")
    
        #Now need to write out the mesh collections with the new positions
        mesh.scale(downsampling_factors)
        mesh.alphas = alphas
        mesh_collection.transform(gems.KvlTransform(requireNumpyArray(np.linalg.inv(transformation_matrix))))
        file_path = os.path.split(transformed_template_file_name) 
        mesh_collection.write(os.path.join(file_path[0],'atlas_level2.txt'))
                          
        #Also the other mesh collection
        mesh_collection.read(mesh_level1)
        mesh_collection.k = K
        mesh_collection.transform(transformation_to_image)
        mesh = mesh_collection.reference_mesh
        
        #Get z-coordinates
        alphas = mesh.alphas
        spineAlphas = alphas[:,-1] #Note: the spine class is the last one in the affine atlas. Might need to change this in the future.
        mask_spine = spineAlphas > 0.01
        z_positions = mesh.points[:,z_dim]
        spine_positions_z = z_positions[mask_spine]
        mask_neck = 0
        neck_pos = 0
        #The position where spine starts (from the brainstem) depends on if the direction is I->S or S->I
        if(z_direction<1): #I->S
            top_ind = np.argwhere(np.argmax(spine_positions_z))
            top_pos = spine_positions_z[top_ind[0]]
            mask_neck = z_positions<top_pos
            neck_pos = z_positions[mask_neck]
            z_dist = top_pos - neck_pos
        else: #S->I
            top_ind = np.argwhere(np.argmin(spine_positions_z))
            top_pos = spine_positions_z[top_ind[0]]
            mask_neck = z_positions>top_pos
            neck_pos = z_positions[mask_neck]
            z_dist = neck_pos - top_pos
        
        deformation_field = np.zeros_like(mesh.points)
        deformation_field[mask_neck,y_dim] = z_dist
        initial_node_positions = mesh.points
        mesh.points = initial_node_positions - res.x*deformation_field
    
        mesh_collection.transform(gems.KvlTransform(requireNumpyArray(np.linalg.inv(transformation_matrix))))
        mesh_collection.write(os.path.join(file_path[0],'atlas_level1.txt'))

    def _neck_cost(self,x,initial_node_positions,deformation_field,calculator,mesh,image_buffer,visualizer):
        mesh.points = initial_node_positions - x*deformation_field
        cost, _ = calculator.evaluate_mesh_position_a(mesh)
        mesh.points = initial_node_positions
        return cost
        

        
        


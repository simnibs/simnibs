import os
import scipy.io
import scipy.ndimage
import numpy as np
from scipy.optimize import minimize_scalar
import logging

import charm_gems as gems
from .Affine import Affine
from .utilities import requireNumpyArray
from .SamsegUtility import readCroppedImages
eps = np.finfo(float).eps


class AffineWholeHead(Affine):

    def saveResults(self, image_name, template_name,
                    save_path, template_coregistered_name,
                    image_to_image_transform, world_to_world_transform):
        # ------ Save Registration Results ------
        template = gems.KvlImage(template_name)
        image = gems.KvlImage(image_name)
        image_to_world_transform = image.transform_matrix.as_numpy_array
        scipy.io.savemat(os.path.join(save_path, 'coregistrationMatrices.mat'),
                         {'imageToImageTransformMatrix': image_to_image_transform,
                          'worldToWorldTransformMatrix': world_to_world_transform})

        # Save the coregistered template. For historical reasons,
        # we applied the estimated transformation to the template...
        # let's do that now
        template_image_to_world_transform = np.asfortranarray(
            image_to_world_transform @ image_to_image_transform)
        template.write(template_coregistered_name,
                       gems.KvlTransform(template_image_to_world_transform))

    def removeNeck(self, mesh):
        alphas = mesh.alphas
        # Note: the spine class is the last one in the affine atlas.
        # Might need to change this in the future.
        spineAlphas = alphas[:, -1]
        mask = spineAlphas > 0.01

        # Get z-coordinates
        zPositions = mesh.points[:, -1]

        # Find the max position for spine and cut from there
        spinePosZ = zPositions[mask]
        maxInds = np.argwhere(np.argmax(spinePosZ))
        maxZCoord = zPositions[maxInds[0]]
        zMask = zPositions < maxZCoord

        # Create a "neck" class which includes everything in the neck
        alphasWithNeck = np.zeros((alphas.shape[0], alphas.shape[1]+1))
        alphasWithNeck[:, 1:] = alphas
        alphasWithNeck[zMask, 1:] = 0
        alphasWithNeck[zMask, 0] = 1.0

        return alphasWithNeck

    def adjust_neck(self, T1, transformed_template_file_name, mesh_level1,
                    mesh_level2, neck_bounds, neck_tissues, visualizer,
                    downsampling_target=3.0):
        logger = logging.getLogger(__name__)
        image_buffer, transformation_to_image, voxel_spacing, cropping = readCroppedImages([T1],transformed_template_file_name)
        # Figure out how much to downsample (depends on voxel size)
        downsampling_factors = np.round(downsampling_target / voxel_spacing)
        downsampling_factors[downsampling_factors < 1] = 1
        image_buffer = image_buffer[
                          ::int(downsampling_factors[0]),
                          ::int(downsampling_factors[1]),
                          ::int(downsampling_factors[2])]
        # Read in the mesh collection, and set K to be very small
        # so that the intensity cost dominates
        K = 1e-20
        mesh_collection = gems.KvlMeshCollection()
        mesh_collection.read(mesh_level2)
        mesh_collection.k = K*np.prod(downsampling_factors)
        mesh_collection.transform(transformation_to_image)
        mesh = mesh_collection.reference_mesh
        # Downsample if needed
        mesh.scale(1/downsampling_factors)

        # Let's define the neck deformation based on the distance
        alphas = mesh.alphas
        # Note: the spine class is the last one in the affine atlas.
        # Might need to change this in the future.
        spineAlphas = alphas[:, 48]
        mask_spine = spineAlphas > 0.01

        # Let's figure out where the z-coordinate is and which way is up
        transformation_matrix = transformation_to_image.as_numpy_array
        z_dim = np.argmax(np.abs(transformation_matrix[0:3, 2]))
        z_direction = transformation_matrix[z_dim, 2]

        # Get z-coordinates
        z_positions = mesh.points[:, z_dim]
        spine_positions_z = z_positions[mask_spine]
        mask_neck = 0
        neck_pos = 0
        # The position where spine starts (from the brainstem) depends on
        # if the direction is I->S or S->I
        if(z_direction < 1):  # I->S
            top_ind = np.argmax(spine_positions_z)
            top_pos = spine_positions_z[top_ind]
            mask_neck = z_positions < top_pos
            neck_pos = z_positions[mask_neck]
            z_dist = top_pos - neck_pos
        else:  # S->I
            top_ind = np.argmin(spine_positions_z)
            top_pos = spine_positions_z[top_ind]
            mask_neck = z_positions > top_pos
            neck_pos = z_positions[mask_neck]
            z_dist = neck_pos - top_pos

        # Okay the distance from the top of the spine defines the amount
        # of deformation in the A->P direction
        # Let's first check where the A-P direction is
        y_dim = np.argmax(np.abs(transformation_matrix[0:3, 1]))
        deformation_field = np.zeros_like(mesh.points)
        deformation_field[mask_neck, y_dim] = z_dist

        # Okay one more trick that seems to work, only consider a subset
        # of the structures to compute the cost
        # These are: air internal, spine, cortical bone and spongy bone
        # NOTE! This need to be read from the setup as otherwise this
        # will fail if the atlas is changed!
        alphas_new = alphas[:, neck_tissues]
        mesh.alphas = alphas_new

        # Let's see how it looks
        visualizer.show(mesh=mesh, images=image_buffer,
                        window_id='Initial Neck', title="Initial Neck")
        image = gems.KvlImage(requireNumpyArray(image_buffer))
        calculator = gems.KvlCostAndGradientCalculator('MutualInformation',
                                                       [image], 'Affine')

        # Check initial cost
        initCost, _ = calculator.evaluate_mesh_position_a(mesh)
        logger.info("Initial cost: " + str(initCost))

        # Optimize the linear deformation, I'll do this now naively using
        # the Nelder-Mead Simplex algorithm. We should be able to compute the
        # gradient for this simple thing, so something smarter should be used
        initial_node_positions = mesh.points

        res = minimize_scalar(self._neck_cost, None, neck_bounds,
                              (initial_node_positions,
                               deformation_field,
                               calculator, mesh,
                               image_buffer,
                               visualizer), method='Bounded')

        mesh.points = initial_node_positions - res.x*deformation_field

        visualizer.show(mesh=mesh, images=image_buffer,
                        window_id='Corrected Neck', title="Corrected Neck")

        # Now need to write out the mesh collections with the new positions
        mesh.scale(downsampling_factors)
        mesh.alphas = alphas
        mesh_collection.transform(gems.KvlTransform(requireNumpyArray(
            np.linalg.inv(transformation_matrix))))
        file_path = os.path.split(transformed_template_file_name)
        mesh_collection.write(os.path.join(file_path[0], 'atlas_level2.txt'))

        # Also the other mesh collection
        mesh_collection.read(mesh_level1)
        mesh_collection.k = K
        mesh_collection.transform(transformation_to_image)
        mesh = mesh_collection.reference_mesh

        # Get z-coordinates
        alphas = mesh.alphas
        # Note: the spine class is the last one in the affine atlas.
        # Might need to change this in the future.
        spineAlphas = alphas[:, 48]
        mask_spine = spineAlphas > 0.01
        z_positions = mesh.points[:, z_dim]
        spine_positions_z = z_positions[mask_spine]
        mask_neck = 0
        neck_pos = 0
        # The position where spine starts (from the brainstem) depends on
        # if the direction is I->S or S->I
        if(z_direction < 1):  # I->S
            top_ind = np.argwhere(np.argmax(spine_positions_z))
            top_pos = spine_positions_z[top_ind[0]]
            mask_neck = z_positions < top_pos
            neck_pos = z_positions[mask_neck]
            z_dist = top_pos - neck_pos
        else:  # S->I
            top_ind = np.argwhere(np.argmin(spine_positions_z))
            top_pos = spine_positions_z[top_ind[0]]
            mask_neck = z_positions > top_pos
            neck_pos = z_positions[mask_neck]
            z_dist = neck_pos - top_pos

        deformation_field = np.zeros_like(mesh.points)
        deformation_field[mask_neck, y_dim] = z_dist
        initial_node_positions = mesh.points
        mesh.points = initial_node_positions - res.x*deformation_field

        mesh_collection.transform(gems.KvlTransform(requireNumpyArray(
            np.linalg.inv(transformation_matrix))))
        mesh_collection.write(os.path.join(file_path[0], 'atlas_level1.txt'))

    def _neck_cost(self, x, initial_node_positions, deformation_field,
                   calculator, mesh, image_buffer, visualizer):
        mesh.points = initial_node_positions - x*deformation_field
        cost, _ = calculator.evaluate_mesh_position_a(mesh)
        mesh.points = initial_node_positions
        return cost

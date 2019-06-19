import copy
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools


def get_opt_grid(tms, msh, target, handle_direction_ref, radius=20, resolution_pos=1, resolution_angle=10,
                 angle_limits=[-30, 30]):
    """
    Determine the coil positions and orientations for bruteforce TMS optimization

    Parameters
    ----------
    tms : simnibs.simulation.sim_struct.TMSLIST object
        TMS simulation object instance
    msh : simnibs.msh.mesh_io.Msh object
        Simnibs mesh object
    target : list [3]
        Coordinates (x, y, z) of cortical target
    radius : float
        Size of region of interest around skin-projected cortical target, where the bruteforce simulations are
        conducted
    resolution_pos : float
        Resolution in mm of the coil positions in the region of interest
    resolution_angle :
        Resolution in deg of the coil positions in the region of interest

    Returns
    -------
    tms : simnibs.simulation.sim_struct.TMSLIST object
        TMS simulation object instance containing the coil positions and orientations to run bruteforce simulations
    """

    # project cortical target to skin surface
    tms_tmp = copy.deepcopy(tms)
    pos_target = tms_tmp.add_position()
    pos_target.centre = target
    pos_target.pos_ydir = [0, 0, 0]
    pos_target.distance = .1
    target_matsimnibs = pos_target.calc_matsimnibs(msh)
    target_skin = target_matsimnibs[0:3, 3]

    # extract ROI
    msh_surf = msh.crop_mesh(elm_type=2)
    msh_skin = msh_surf.crop_mesh([5, 1005])
    elm_center = np.mean(msh_skin.nodes.node_coord[msh_skin.elm.node_number_list[:, 0:3]-1], axis=1)
    elm_mask_roi = np.linalg.norm(elm_center - target_skin, axis=1) < 1.2*radius
    elm_center_roi = elm_center[elm_mask_roi, ]
    node_number_list_roi = msh_skin.elm.node_number_list[elm_mask_roi, 0:3]-1
    nodes_roi = msh_skin.nodes.node_coord[node_number_list_roi]
    nodes_roi = np.reshape(nodes_roi, (nodes_roi.shape[0] * 3, 3))
    nodes_roi_mean = np.mean(nodes_roi, axis=0)[np.newaxis, :]
    nodes_roi_zeromean = nodes_roi - nodes_roi_mean

    # tangential plane of target_skin point
    u, s, vh = np.linalg.svd(nodes_roi_zeromean)
    vh = vh.transpose()

    # define regular grid and rotate it to head space
    coords_plane = np.array(np.meshgrid(np.linspace(-radius, radius, 2*radius/resolution_pos + 1),
                                  np.linspace(-radius, radius, 2*radius/resolution_pos + 1))).T.reshape(-1, 2)
    coords_plane = coords_plane[np.linalg.norm(coords_plane, axis=1) <= radius]
    coords_plane = np.dot(coords_plane, vh[:, :2].transpose()) + target_skin

    # project grid-points to skin-surface
    P1 = msh_skin.nodes.node_coord[node_number_list_roi[:, 0], ]
    P2 = msh_skin.nodes.node_coord[node_number_list_roi[:, 1], ]
    P3 = msh_skin.nodes.node_coord[node_number_list_roi[:, 2], ]

    normals = np.cross(P2 - P1, P3 - P1)
    coords_mapped = np.zeros(coords_plane.shape)

    for i, c in enumerate(coords_plane):
        Q1 = c + 1e2 * vh[:, 2]
        Q2 = c - 1e2 * vh[:, 2]

        # point of intersection of infinite half-plane from triangle
        P0 = Q1[np.newaxis, :] + (np.sum((P1 - Q1) * normals, axis=1) /
                                  np.dot(Q2 - Q1, normals.transpose()))[:, np.newaxis] * (Q2 - Q1)[np.newaxis, :]

        # check if intersection points are inside triangle
        inside = np.logical_and(
            np.logical_and(np.sum(np.cross(P2 - P1, P0 - P1) * normals, axis=1) >= 0,
                                np.sum(np.cross(P3 - P2, P0 - P2) * normals, axis=1) >= 0),
                                np.sum(np.cross(P1 - P3, P0 - P3) * normals, axis=1) >= 0)

        coords_mapped[i, ] = P0[inside, ]

    # determine handle directions
    # project handle_direction_ref to reference plane
    # handle_direction_ref_proj = np.dot(handle_direction_ref, vh[:,:2])
    # handle_directions = np.tile(np.array(handle_direction_ref), (coords_mapped.shape[0], 1))

    # determine rotation matrices around z-axis of coil and rotate
    angles = np.linspace(angle_limits[0], angle_limits[1], (angle_limits[1] - angle_limits[0])/resolution_angle + 1)
    handle_directions = np.zeros((len(angles), 3))

    for i, a in enumerate(angles):
        mat_rot = np.array([[np.cos(a/180.*np.pi), -np.sin(a/180.*np.pi), 0],
                            [np.sin(a/180.*np.pi), np.cos(a/180.*np.pi), 0],
                            [0, 0, 1]])
        handle_directions[i, ] = np.dot(handle_direction_ref, mat_rot)

    # combine coil positions and orientations
    po = list(itertools.product(coords_mapped, handle_directions))

    # write coordinates in TMS object
    for c, h in po:
        pos = tms.add_position()
        pos.centre = c.tolist()
        pos.pos_ydir = h.tolist()
        pos.distance = 0.

    return tms


    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(nodes_roi[:, 0], nodes_roi[:, 1], nodes_roi[:, 2])
    # ax.quiver(target_skin[0], target_skin[1], target_skin[2],
    #           5*vh[0, 0], 5*vh[1, 0], 5*vh[2, 0], color='k')
    # ax.quiver(target_skin[0], target_skin[1], target_skin[2],
    #           5*vh[0, 1], 5*vh[1, 1], 5*vh[2, 1], color='r')
    # ax.quiver(target_skin[0], target_skin[1], target_skin[2],
    #           5 * vh[0, 2], 5 * vh[1, 2], 5 * vh[2, 2], color='g')
    #
    # ax.quiver(c[0], c[1], c[2], Q1[0], Q1[1], Q1[2], color='k')
    # ax.quiver(c[0], c[1], c[2], Q2[0], Q2[1], Q2[2], color='g')
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    # ax.quiver(c[0], c[1], c[2], 1e2 * vh[0, 2], 1e2 * vh[1, 2], 1e2 * vh[2, 2], color='k')
    # ax.quiver(c[0], c[1], c[2], -1e2 * vh[0, 2], -1e2 * vh[1, 2], -1e2 * vh[2, 2], color='kk')
    # ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], color="r")
    # ax.scatter(coords_mapped[:, 0], coords_mapped[:, 1], coords_mapped[:, 2], color="g")
    #
    #
    #
    # ax.scatter(P0[inside, ][:, 0], P0[inside, ][:, 1], P0[inside, ][:, 2], color="g")
    # ax.scatter(elm_center_roi[inside, ][:, 0], elm_center_roi[inside, ][:, 1], elm_center_roi[inside, ][:, 2], color="g")
    # ax.scatter(nodes_roi[:, 0], nodes_roi[:, 1], nodes_roi[:, 2])

    # target_skin_center, target_skin_elm_idx = msh_skin.find_closest_element(target_skin, return_index=True)
    # target_skin_nodes = msh_skin.nodes.node_coord[msh_skin.elm.node_number_list[target_skin_elm_idx][0:3]-1, ]
    # target_skin_tangents = np.array([target_skin_nodes[1, ] - target_skin_nodes[0, ],
    #                                  target_skin_nodes[2, ] - target_skin_nodes[0, ]]).transpose()


# def calc_coil_pos_from_matsimnibs(matsimnibs):
#     """ Calculate the matsimnibs matrix for TMS simulations
#
#     Parameters
#     -----------
#     center: np.ndarray
#         Position of the center of the coil, will be projected to the skin surface
#     pos_ydir: np.ndarray
#         Position of the y axis in relation to the coil center
#     distance: float
#         Distance from the center
#     skin_surface: list
#         Possible tags for the skin surface (Default: [5, 1005])
#
#     Returns
#     -------
#     matsimnibs: 2d np.ndarray
#         Matrix of the format
#         x' y' z' c
#         0  0  0  1
#         y' is the direction of the coil
#         z' is a direction normal to the coil, points inside the head
#
#     """
#     x = matsimnibs[:3, 0]
#     y = matsimnibs[:3, 1]
#     z = matsimnibs[:3, 2]
#     c = matsimnibs[:3, 3]
#
#
#     msh_surf = self.crop_mesh(elm_type=2)
#     msh_skin = msh_surf.crop_mesh(skin_surface)
#     closest = np.argmin(np.linalg.norm(msh_skin.nodes.node_coord - center, axis=1))
#     center = msh_skin.nodes.node_coord[closest]
#     # Y axis
#     y = pos_ydir - center
#     if np.isclose(np.linalg.norm(y), 0.):
#         raise ValueError('The coil Y axis reference is too close to the coil center! ')
#     y /= np.linalg.norm(y)
#     # Normal
#     normal = msh_skin.nodes_normals().value[closest]
#     if np.isclose(np.abs(y.dot(normal)), 1.):
#         raise ValueError('The coil Y axis normal to the surface! ')
#     z = -normal
#     # Orthogonalize y
#     y -= z * y.dot(z)
#     y /= np.linalg.norm(y)
#     # Determine x
#     x = np.cross(y, z)

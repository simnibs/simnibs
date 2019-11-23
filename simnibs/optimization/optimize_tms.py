"""
Functions to find optimal TMS coil positions ior a given cortical target. See examples/tms_optimization.py for
some examples on how to use these.

Written by Ole Numssen & Konstantin Weise, 2019.
Adapted by Guilherme Saturnino, 2019
"""
import numpy as np


def _create_grid(mesh, pos, distance, radius, resolution_pos):
    ''' Creates a position grid '''
    # extract ROI
    msh_surf = mesh.crop_mesh(elm_type=2)
    msh_skin = msh_surf.crop_mesh([5, 1005])
    target_skin = msh_skin.find_closest_element(pos)
    elm_center = msh_skin.elements_baricenters()[:]
    elm_mask_roi = np.linalg.norm(elm_center - target_skin, axis=1) < 1.2 * radius
    elm_center_zeromean = (
        elm_center[elm_mask_roi] -
        np.mean(elm_center[elm_mask_roi], axis=0)
    )
    msh_roi = msh_skin.crop_mesh(elements=msh_skin.elm.elm_number[elm_mask_roi])
    # tangential plane of target_skin point
    u, s, vh = np.linalg.svd(elm_center_zeromean)
    vh = vh.transpose()

    # define regular grid and rotate it to head space
    coords_plane = np.array(
        np.meshgrid(
            np.linspace(-radius, radius, int(2 * radius / resolution_pos + 1)),
            np.linspace(-radius, radius, int(2 * radius / resolution_pos + 1)),
        )
    ).T.reshape(-1, 2)
    coords_plane = coords_plane[np.linalg.norm(coords_plane, axis=1) <= radius]
    coords_plane = np.dot(coords_plane, vh[:, :2].transpose()) + target_skin

    # project grid-points to skin-surface
    coords_mapped = []
    coords_normals = []
    normals_roi = msh_roi.triangle_normals(smooth=1)
    for i, c in enumerate(coords_plane):
        # Query points inside/outside the surface
        q1 = c + 1e2 * vh[:, 2]
        q2 = c - 1e2 * vh[:, 2]
        idx, pos = msh_roi.intercept_ray(q1, q2)
        if idx is not None:
            coords_normals.append(normals_roi[idx])
            coords_mapped.append(pos)

    coords_mapped = np.array(coords_mapped)
    coords_normals = np.array(coords_normals)
    inside = np.linalg.norm(coords_mapped - target_skin, axis=1) <= radius
    coords_mapped += distance * coords_normals
    return coords_mapped[inside], coords_normals[inside]

def _rotate_system(R, angle_limits, angle_res):
    ''' Rotates the vector "y" aroud "z" between the given limits and in the given
    resolution and return rotation matrices'''
    # Define rotation matrix around Z
    n_steps = int((angle_limits[1] - angle_limits[0])/angle_res + 1)
    angles = np.linspace(angle_limits[0], angle_limits[1], n_steps)
    angles = np.deg2rad(angles[(angles > -180.1) * (angles < 180.)])
    matrices = []
    for a in angles:
        Rz = np.array((
            (np.cos(a), -np.sin(a), 0),
            (np.sin(a), np.cos(a), 0),
            (0, 0, 1),
        ))
        matrices.append(R.dot(Rz))
    return matrices


def get_opt_grid(mesh, pos, handle_direction_ref=None, distance=1., radius=20,
                 resolution_pos=1, resolution_angle=20, angle_limits=None):
    """ Determine the coil positions and orientations for bruteforce TMS optimization

    Parameters
    ----------
    mesh: simnibs.msh.mesh_io.Msh object
        Simnibs mesh object
    pos: ndarray
        Coordinates (x, y, z) of reference position
    handle_direction_ref (optinal): list of float or np.ndarray
        Vector of handle prolongation direction, in relation to "pos". (Default: do
        not select a handle direction and scan rotations from -180 to 180)
    distance: float or None
        Coil distance to skin surface [mm]. (Default: 1.)
    radius: float or None
        Radius of region of interest around the reference position, where the
        bruteforce simulations are conducted
    resolution_pos: float or None
        Resolution in mm of the coil positions in the region of interest.
    resolution_angle: float or None
        Resolution in deg of the coil positions in the region of interest (Default: 20)
    angle_limits: list of float or None
        Range of angles to get coil rotations for (Default: [-180, 180])

    Returns
    -------
    matsimnibs_list: list
        list of MATSIMNIBS matrices
    """
    # creates the spatial grid
    coords_mapped, coords_normals = _create_grid(
        mesh, pos, distance, radius, resolution_pos)
    
    # Determines the seed y direction
    if handle_direction_ref is None:
        y_seed = np.array([0., 1., 0.])
    else:
        y_seed = np.array(handle_direction_ref) - np.array(pos)
        if np.isclose(np.linalg.norm(y_seed), 0.):
            raise ValueError('The coil Y axis reference is too close to the coil center! ')
    if angle_limits is None:
        angle_limits = -180, 180

    matrices = []
    for p, n in zip(coords_mapped, coords_normals):
        z = -n
        y = y_seed - (z * y_seed.dot(z))
        y /= np.linalg.norm(y)
        x = np.cross(y, z)
        R = np.array([x, y, z]).T
        rotated = _rotate_system(R, angle_limits, resolution_angle)
        for r in rotated:
            matrices.append(
                np.vstack((
                    np.hstack((r, p[:, None])),
                    [0, 0, 0, 1]
                ))
            )

    return matrices


def plot_matsimnibs_list(matsimnibs_list, values, field_name, fn_geo):
    ''' Plots the center and the y vector of each matsimnibs matrix as a geo file

    Parameters
    -------------
    matsimnibs_list: list
        list of matsimnibs matrices
    values: array
        Value to assign to each matsimnibs
    field_name: str
        Name of field being printed
    fn_geo: str
        Name of output geo file
    '''
    with open(fn_geo, 'w') as f:
        f.write('View"' + field_name + '"{\n')
        for mat, v in zip(matsimnibs_list, values):
            c = mat[:3, 3]
            y = mat[:3, 1] * v
            f.write(
                "VP(" + ", ".join([str(i) for i in c]) + ")"
                "{" + ", ".join([str(i) for i in y]) + "};\n")
            f.write(
                "SP(" + ", ".join([str(i) for i in c]) + ")"
                "{" + str(v) + "};\n")
        f.write("};\n")


def define_target_region(mesh, target_position, target_radius, tags, elm_type=4):
    ''' Defines a target based on a position, a radius and an element tag

    Paramters
    ------------
    mesh: simnibs.mesh_io.msh
        Mesh
    target_position: array of size (3,)
        Position of target
    target_radius: float
        Size of target
    tags: array of ints
        Tag where the target is located
    elm_type: int
        Type of target element (4 for tetrahedra and 2 for triangles)

    Returns
    -------
    elements: array of size (n,)
        Numbers (1-based) of elements in the tag
    '''
    bar = mesh.elements_baricenters()[:]
    dist = np.linalg.norm(bar - target_position, axis=1)
    elm = mesh.elm.elm_number[
        (dist < target_radius) *
        np.isin(mesh.elm.tag1, tags) *
        np.isin(mesh.elm.elm_type, elm_type)
    ]
    return elm

"""
TMS coil position/orientation optimizations functions.


Authors: Ole Numssen, Konstantin Weise
"""
from simnibs import sim_struct, run_simnibs
import numpy as np


n_cpu = 10
angles = [22.5, 45, 67.5]
radius = 50
resolution = 10
brain_target = [0,-10,17]


def calc_coil_pos_from_matsimnibs(matsimnibs):
    """ Calculate the matsimnibs matrix for TMS simulations

    Parameters
    -----------
    center: np.ndarray
        Position of the center of the coil, will be projected to the skin surface
    pos_ydir: np.ndarray
        Position of the y axis in relation to the coil center
    distance: float
        Distance from the center
    skin_surface: list
        Possible tags for the skin surface (Default: [5, 1005])

    Returns
    -------
    matsimnibs: 2d np.ndarray
        Matrix of the format
        x' y' z' c
        0  0  0  1
        y' is the direction of the coil
        z' is a direction normal to the coil, points inside the head

    """
    x = matsimnibs[:3, 0]
    y = matsimnibs[:3, 1]
    z = matsimnibs[:3, 2]
    c = matsimnibs[:3, 3]


    msh_surf = self.crop_mesh(elm_type=2)
    msh_skin = msh_surf.crop_mesh(skin_surface)
    closest = np.argmin(np.linalg.norm(msh_skin.nodes.node_coord - center, axis=1))
    center = msh_skin.nodes.node_coord[closest]
    # Y axis
    y = pos_ydir - center
    if np.isclose(np.linalg.norm(y), 0.):
        raise ValueError('The coil Y axis reference is too close to the coil center! ')
    y /= np.linalg.norm(y)
    # Normal
    normal = msh_skin.nodes_normals().value[closest]
    if np.isclose(np.abs(y.dot(normal)), 1.):
        raise ValueError('The coil Y axis normal to the surface! ')
    z = -normal
    # Orthogonalize y
    y -= z * y.dot(z)
    y /= np.linalg.norm(y)
    # Determine x
    x = np.cross(y, z)


# ### General Information
S = sim_struct.SESSION()
S.fnamehead = 'ernie.msh'  # head mesh
S.pathfem = 'tms'  # Directory for the simulation

# ## Define the TMS simulation
tms = S.add_tmslist()
tms.fnamecoil = 'Magstim_70mm_Fig8.nii.gz'  # Choose a coil from the ccd-files folder

# get coil position at start position

pos = tms.add_position()
pos.centre = brain_target
pos.pos_ydir = angles[1]
pos.distance = 0


# TODO: magic here to get coil positions
pos.calc_matsimnibs(tms.mesh)
positions = dict()

# # Define the coil position
for key, val in positions.items():
    pos = tms.add_position()
    pos.centre = val['centre']
    pos.pos_ydir = val['ydir']
    pos.distance = val['dist']

# Run Simulation
run_simnibs(S, n_cpu)

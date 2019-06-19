"""
TMS coil position/orientation optimizations functions.


Authors: Ole Numssen, Konstantin Weise
"""
from simnibs.msh import mesh_io
from simnibs import sim_struct, run_simnibs
import numpy as np
import pickle


n_cpu = 10
angles = [22.5, 45, 67.5]
radius = 50
resolution = 10
brain_target = np.mean(np.array([[-27.17, -17.94, 69.94],
                                 [-28.15, -18.02, 69.24],
                                 [-27.72, -18.95, 69.82]]), axis=0).tolist()
# brain_target = [-29.2, -20.4, 95.]

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


# General Information
S = sim_struct.SESSION()
mesh = '/data/pt_01756/tmp/optim/15484.08_fixed.msh'
S.fnamehead = mesh  # head mesh
S.pathfem = '/data/pt_01756/tmp/optim'  # Directory for the simulation

# Define the TMS simulation
tms = S.add_tmslist()
tms.fnamecoil = '/data/pt_01756/coils/ccd/MagVenture_MC_B60_1233.nii.gz'  # Choose a coil from the ccd-files folder

# get coil position at start position

pos = tms.add_position()
pos.centre = brain_target
pos.pos_ydir = angles[1]
pos.distance = 1

print("prepare")
# mesh = mesh_io.read_msh(mesh)
# mesh.fix_surface_labels()
mesh = pickle.load(open("/data/pt_01756/tmp/optim/mesh.pckl", 'rb'))

print("calc matsimnibs")
mat = pos.calc_matsimnibs(mesh)  # get skull coil position for cortex target
# TODO: magic here to get coil positions
positions = dict()

# # Define the coil position
for key, val in positions.items():
    pos = tms.add_position()
    pos.centre = val['centre']
    pos.pos_ydir = val['ydir']
    pos.distance = val['dist']

# Run Simulation
run_simnibs(S, n_cpu)

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


def read_msh_from_pckl(fn, m=None):
    ''' Reads a gmsh '.msh' file

    Parameters
    ------------
    fn: str
        File name
    m: simnibs.msh.Msh (optional)
        Mesh structure to be overwritten. If unset, will create a new structure

    Returns
    --------
    msh: simnibs.msh.Msh
        Mesh structure
    '''
    print("This is the monkey-patched version.")
    assert fn.endswith('.msh')
    fn = fn[:-3] + "pckl"
    return pickle.load(open(fn, 'rb'))


from ..msh import mesh_io
mesh_io.read_msh = read_msh_from_pckl

# General Information
S = sim_struct.SESSION()
mesh = '/data/pt_01756/tmp/optim/mesh.pckl'
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

"""
Example for TMS coil position/orientation optimization.

optimize_tms_coil_pos() implements a brute-force approach to find the optimal coil position & rotation to maximize normE
at a given cortical target. The target is projected to the nearest position at the skin surface and candidate coil
positions/rotations are then distributed equally in a circular plane around that location.
The number of candiate coil positions/rotations is defined by the two resolution parameters (_pos and
angle), the radius of the circular plane, and angle_limits.

Authors: Ole Numssen, Konstantin Weise, 2019
"""

from simnibs.simulation.opt import optimize_tms_coil_pos
from simnibs.simulation.sim_struct import SESSION, TMSLIST


# define cortical target where normE should be maximal
target = [-27.17, -17.94, 69.94]
handle_direction_ref = [-2.8, 7, 8.1]  # ~45° towards fiss. long.

coil_fn = 'MagVenture_MC_B60_REF.nii.gz'
tensor_fn = "CTI_vn_tensor.nii.gz"
n_cpu = 20

# setup session with parameters
session = SESSION()
session.fnamehead = 'ernie.msh'
session.pathfem = '/data/project/optimization_left/'
session.fname_tensor = tensor_fn  # not needed for isotropic
session.open_in_gmsh = False

# setup tmslist with options
tms_list = TMSLIST()
tms_list.fnamecoil = coil_fn
tms_list.write_hdf5 = True  # write .hdf5 files to speed up read in
tms_list.remove_msh = True  # remove .msh after simulation to save storage space
tms_list.create_visuals = False  # don't create summary after simulation to speed up processing
tms_list.anisotropy_type = "vn"  # for isotropic: 'scalar'

session.add_poslist(tms_list)

# call optimization procedure with default arguments
best_positions_right = optimize_tms_coil_pos(session=session,
                                             target=target,
                                             n_cpu=n_cpu)

"""Full argument example"""
target = [27.17, -17.94, 69.94]  # x-coord in right hemisphere
handle_direction_ref = [2.8, 7, 8.1]  # adapt coil handle to right hemisphere

# these are the default values
radius = 20
resolution_pos = 1.5
resolution_angle = 15
angle_limits = [-60, 60]
distance = 1.

# setup session with parameters
session = SESSION()
session.fnamehead = 'ernie.msh'
session.pathfem = '/data/project/optimization_right/'
session.fname_tensor = tensor_fn  # not needed for isotropic
session.open_in_gmsh = False

# setup tmslist with options
tms_list = TMSLIST()
tms_list.fnamecoil = coil_fn
tms_list.write_hdf5 = True  # write .hdf5 files to speed up read in
tms_list.remove_msh = True  # remove .msh after simulation to save storage space
tms_list.create_visuals = False  # don't create summary after simulation to speed up processing
tms_list.anisotropy_type = "vn"  # for isotropic: 'scalar'

session.add_poslist(tms_list)
best_positions_left = optimize_tms_coil_pos(session=session,
                                            target=target,
                                            n_cpu=n_cpu,
                                            handle_direction_ref=handle_direction_ref,
                                            radius=radius,
                                            resolution_angle=resolution_angle,
                                            resolution_pos=resolution_pos,
                                            angle_limits=angle_limits,
                                            distance=distance)

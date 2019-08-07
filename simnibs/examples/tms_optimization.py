"""
Example for TMS coil position/orientation optimization.

optim_tms.optimize_tms_coil_pos() implements a brute-force approach to find the optimal coil position & rotation to
maximize the quantity of interest (qoi) (normE/...) at a given cortical target. The cortical target is projected to
the nearest position at the skin surface and candidate coil positions/rotations are distributed equally in a circular
plane around that location. The number of candiate coil positions/rotations is defined by the two resolution
parameters (resolution_pos and resolution_angle), the radius of the circular plane (radius), and the the range of the
coil angles (angle_limits).

The qoi is calculated for each of these candidate coil positions (whole brain) and stored in
<pathfem>/<optim_name>.hdf5 and the qoi at the cortical target is stored in  <pathfem>/<optim_name>_<datetime>.csv
with the corresponding coil positions and orientations.

nnav.simnibs2nnav() transforms simnibs positions/orientations to Localite TMS Navigator instrument markers.


See below for an examplary workflow to go from target location to optimal coil position to instrument marker file.

Authors: Ole Numssen, Konstantin Weise, 2019
"""

from simnibs.simulation.optim_tms import optimize_tms_coil_pos
from simnibs.simulation.sim_struct import TMSLIST
from simulation.optim_tms import TMSOPTIMIZATION
import numpy as np

########################################################################
# ### Minimal Example
########################################################################

# define cortical target where normE should be maximal
from utils.nnav import simnibs2nnav, write_tms_navigator_im

# Select target
target = np.array([-27.17, -17.94, 69.94])  # ~ left M1
coil_fn = 'MagVenture_MC_B60_REF.nii.gz'  # Choose a coil from the ccd-files folder

tms_optim = TMSOPTIMIZATION()
tms_optim.fnamehead = 'ernie.msh'
tms_optim.pathfem = '/data/project/optimization_left/'

poslist = TMSLIST()
poslist.fnamecoil = coil_fn  

tms_optim.add_poslist(poslist)
tms_optim.optim_name = 'm1_l'

# call optimization procedure with default arguments
best_positions_right = optimize_tms_coil_pos(tms_optim=tms_optim,
                                             target=target)
# best_positions_right['best_conds']  # is TMSOPTIMZATION object with best coil position(s)
# best_positions_right['simulations']  # is .hdf5 filename where all simulations are stored
# best_positions_right['tms_optim']  # is TMSOPTIMIZATION object with all simulations that where calculated


########################################################################
# ### Full Argument Example
########################################################################

# define cortical target where E should be maximal
target = np.array([27.17, -17.94, 69.94])  # x-coord in right hemisphere
handle_direction_ref = np.array([2.8, 7, 8.1])  # adapt coil handle to right hemisphere
handle_direction_ref += target  # handle direction reference is point, not vector
tensor_fn = "CTI_vn_tensor.nii.gz"

# these are the default values
radius = 20  # 20mm radius around target
resolution_pos = 1.5  # 1.5mm raster for coil centers
resolution_angle = 15  # 15° steps for orientation
angle_limits = [-60, 60]  # -60° to +60° around handle_direction_ref
distance = 1.  # coil distance

n_cpu = 8

# setup TMSOPTIMIZATION with global options
tms_optim = TMSOPTIMIZATION()
tms_optim.fnamehead = 'ernie.msh'
tms_optim.pathfem = '/data/project/optimization_right/'

poslist = TMSLIST()
poslist.fnamecoil = coil_fn  # Choose a coil from the ccd-files folder
tms_optim.add_poslist(tms_optim)
tms_optim.optim_name = 'm1_r'

tms_optim.fname_tensor = tensor_fn
tms_optim.anisotropy_type = 'vn'  # for isotropic: 'scalar'
tms_optim.save_fields = 'E'
tms_optim.fields = 'E'
tms_optim.qoi = 'E'


best_positions_left = optimize_tms_coil_pos(tms_optim=tms_optim,
                                            target=target,
                                            n_cpu=n_cpu,
                                            handle_direction_ref=handle_direction_ref,
                                            radius=radius,
                                            resolution_angle=resolution_angle,
                                            resolution_pos=resolution_pos,
                                            angle_limits=angle_limits,
                                            distance=distance)

# best_positions_right['best_conds']  # TMSOPTIMZATION object with best coil position(s)
# best_positions_right['simulations']  # .hdf5 filename where all simulations are stored
# best_positions_right['tms_optim']  # TMSOPTIMIZATION object with all simulations that where calculated


########################################################################
# ### Resume crashed optimization
########################################################################
tms_optim = TMSOPTIMIZATION()
tms_optim.read_mat_struct("/data/project/m1_l__20190131-120500.mat")
tms_optim.resume = True
results = optimize_tms_coil_pos(tms_optim=tms_optim)


########################################################################
# ### Write Localite TMS Navigator instrument marker for best position
########################################################################
results = optimize_tms_coil_pos(tms_optim=tms_optim)

# compute neuronavigation system coordinates for optimal coil position
fn_exp_nii = "ernie.nii"  # The T1 used in the neuronavigation
fn_conform_nii = "ernie_cornform.nii"  # The T1 used for simnibs mesh generation
# fn_exp_nii == fn_conform_nii speeds up marker generation, as no coregistration has to be performed

mat_tms_navigator = simnibs2nnav(fn_exp_nii,
                                 fn_conform_nii,
                                 results['best_conds'])

# write InstrumentMarker .xml file.
xml_fn = '/data/project/best_cond_instrument_marker.xml'
write_tms_navigator_im(mat_tms_navigator, xml_fn)

#
# At the moment (July 2019) there is no official way to import instrument marker files into the TMS Localite
# software. A workaround is to copy-paste the content of the above generated .xml file into an existing TMS Navigator
# Session instrument marker file:
# subject_name/subject_name_date/Sessions/Session_latest/InstrumentMarkers/InstrumentMarker_latest.xml
# and load the session afterwards.

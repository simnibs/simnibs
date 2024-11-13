"""
    TMS optimization using grid-search to maximize the field
    strength in a ROI defined via a mask in fsaverage space
    
    Copyright (c) 2024 SimNIBS developers. Licensed under the GPL v3.
"""
import os
import simnibs
from simnibs import opt_struct

subject_path = 'm2m_ernie'

# Define region-of-interest from mask in fsaverage space        
roi = simnibs.RegionOfInterest()
roi.subpath = subject_path
roi.method = 'surface'
roi.surface_type = 'central'
roi.mask_path = os.path.join(simnibs.SIMNIBSDIR, 'examples','utilities','P1_LH_M1_control')
roi.mask_space = 'fs_avg_lh'
# uncomment for visual control of ROI:
# roi.write_visualization('', 'roi_P1_LH_M1_control')


# Initialize structure
tms_opt = opt_struct.TMSoptimize()
tms_opt.open_in_gmsh = False # no
# Subject folder
tms_opt.subpath = subject_path
# Select output folder
tms_opt.pathfem = 'tms_optimization_region/'
# Select the coil model
tms_opt.fnamecoil = os.path.join('legacy_and_other','Magstim_70mm_Fig8.ccd')
# Select targets for the optimization. Every target will be the center of a sphere with radius tms.opt.target_size. The union of every sphere will form the target ROI
tms_opt.multiple_targets = roi.get_nodes()
# Size of every region in multiple targets
tms_opt.target_size = 3.0
tms_opt.angle_resolution = 120 # very low resolution for faster testing
tms_opt.spatial_resolution = 10
# Optional: Use the MKL PARDISO solver
# Will make the simulations much faster
# but has large (approx 12GB) memory usage
tms_opt.solver_options = 'pardiso'
# Run optimization to get optimal coil position
opt_pos=tms_opt.run(save_mat=False, cpus=1)

'''
Transform coil positions from MNI space to subject space

Please run with the command

    simnibs_python transform_coilpos.py

from the directory containing the m2m_ernie example datatset

visit www.simnibs.org for more information
'''

import simnibs


# Option 1: coil position as 4x4 matrix in MNI space (here: above left handknob)
# ------------------------------------------------------------------------------
matsimnibs_MNI = [[-0.747,   0.379,   0.546, -45.8],
                  [ 0.344,   0.923,  -0.169, -19.1],
                  [-0.568,   0.062,  -0.827,  84.0],
                  [ 0.   ,   0.   ,   0.   ,   1.0]]

# coil-skin distance in [mm]; (optional; standard: 0, i.e. the coil center is
# projected on the skin of the subject)
coil_skin_distance = 4

# convert to subject space
matsimnibs = simnibs.mni2subject_coilpos(matsimnibs_MNI, 'm2m_ernie',coil_skin_distance)

print(matsimnibs)


# Option 2: coil position defined by its coil center and unit vectors 
# for the y- and z-axes  in MNI space
# -------------------------------------------------------------------
center_MNI = [-45.8, -19.1,  84.] # corresponds to the 4th column of matsimnibs_MNI
ydir_MNI = [0.379, 0.923, 0.062] # 2nd column
zdir_MNI = [0.546, -0.169, -0.827] # 3rd column
coilpos_MNI = (center_MNI, ydir_MNI, zdir_MNI)

coil_skin_distance = 4

# convert to subject space
matsimnibs = simnibs.mni2subject_coilpos(coilpos_MNI, 'm2m_ernie', coil_skin_distance)

print(matsimnibs)


# Notes:
#   A valid matsimnibs can be obtained by, e.g.:
#    1) use the simnibs_gui to plan a coil position on the m2m_MNI152 head
#    2) save the simulation (use File-->Save in the gui, e.g. as 'myposition.mat')
#    3) load it into python: 
#       S = simnibs.sim_struct.SESSION('myposition.mat')
#    4) this will give the first coil position of the first position list:
#       import numpy as np
#       print(np.array(S.poslists[0].pos[0].matsimnibs))    
#  
#  For visualization of the converted position:
#    1) create geo-file:
#       import numpy as np
#       from simnibs.utils.transformations import write_csv_positions, csv_to_geo
#
#       write_csv_positions(
#                           'test.csv', 
#                           ['CoilPos'],
#                           matsimnibs[:3,3].reshape(1,3), 
#                           ['transformed_pos'], 
#                           [np.hstack((matsimnibs[:3,2],matsimnibs[:3,1],coil_skin_distance))]
#                           )
#       csv_to_geo('test.csv', 'test.geo')
#   2) open ernie.msh in gmsh, and use File-->Merge to load test.geo
 
"""
Example on running a SimNIBS tDCS simulation for a Nx1 montage in Python.
    
    Run with:
        simnibs_python tDCS_Nx1.py

    Copyright (c) 2020 SimNIBS developers. Licensed under the GPL v3.
"""

from simnibs import sim_struct, run_simnibs


###     SETUP
###################
S = sim_struct.SESSION()
S.subpath = 'm2m_ernie'  # m2m-folder of the subject
S.pathfem = 'tdcs_Nx1' # output directory for the simulation
S.map_to_surf = True # map to subject's middle gray matter surface (optional)

tdcs_list = S.add_tdcslist()
tdcs_list.currents = 0.001  # Current flow through center channel (A)

# define the center electrode
center = tdcs_list.add_electrode()
center.centre = 'C3'  # Place it over C3
center.shape = 'ellipse'  # round shape
center.dimensions = [30, 30]  # 30 mm diameter
center.thickness = [2, 1]  # 2 mm rubber electrodes on top of 1 mm gel layer

# when standard parameters are OK, using the following line 
# is enough to set up the surround electrodes:
# tdcs_list.expand_to_center_surround(S.subpath)


# optional parameters:
radius_surround = 60 # distance (centre-to-centre) between the centre and 
                     # surround electrodes (optional; standard: 50 mm)
                     # either a single number or an array with N entries
                     # (N: number of electrodes)             
pos_dir_1stsurround = 'C4' # a position indicating the direction in which the 
                           # first surround electrode should be placed 
                          # (optional; standard: None)                        
N = 4 # number of surround electrodes (optional; standard: 4)
multichannel = True # when set to True: Simulation of multichannel stimulator 
                     # with each suround channel receiving 1/N-th of the
                     # center channel (optional; standard: False, i.e. all 
                     # surround electrodes connected to the same channel)
phis_surround = None # Angles in degree at which the electrodes will be place 
                     # relative to the direction defined by pos_dir_1stsurround.
                     # (optional; standard: None, resulting in equal distances
                     # between surround electrodes)
                     
tdcs_list.expand_to_center_surround(S.subpath, radius_surround, N, 
                                    pos_dir_1stsurround, multichannel,
                                    phis_surround)

### RUN SIMULATION
###################
run_simnibs(S)

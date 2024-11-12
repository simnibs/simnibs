''' example script that runs a simnibs tDCS simulation to demonstrate 
    several post-processing options and the modeling of more 
    complex electrode shapes
    
    Run with:

    simnibs_python tDCS_advanced.py

    Copyright (c) 2021 SimNIBS developers. Licensed under the GPL v3.
'''
from simnibs import sim_struct, run_simnibs


# specify general parameters
S = sim_struct.SESSION()
S.subpath = 'm2m_ernie'  # m2m-folder of the subject
S.pathfem = 'tdcs_advanced_simu'  # Directory for the simulation

S.fields = 'eEjJ'  # Save the following results:
                   #  e: Electric field magnitude
                   #  E: Electric field vector
                   #  j: Current density magnitude
                   #  J: Current density vector
                   
S.map_to_surf = True    #  Map to subject's middle gray matter surface
S.map_to_fsavg = True   #  Map to FreeSurfer's FSAverage group template
S.map_to_vol = True     #  Save as nifti volume
S.map_to_MNI = True     #  Save in MNI space
S.tissues_in_niftis = [1,2,3]  # Results in the niftis will be masked 
                               # to only show WM (1), GM (2), CSF(3)
                               # (standard: only GM)
                               # To get fields everywhere: 
                               #    S.tissues_in_niftis = 'all'

S.open_in_gmsh = True  # show results in gmsh (not for the the niftis)


# specify TDCS simulation
tdcs = S.add_tdcslist()
tdcs.currents = [0.001, -0.001]  # Current flow though each channel (A)

# first electrode
el1 = tdcs.add_electrode()
el1.channelnr = 1  # Connect the electrode to the first channel
el1.centre = 'C3'  # Place it over C3
el1.pos_ydir = 'C4'  # Position on the scalp to define 
                     # the electrode orientation. The electrode's 
                     # y-axis will point from the centre to pos_ydir 
el1.shape = 'rect'  # Rectangular electrode
                    # Other shapes: 
                    #  'ellipse' (includes circular electrodes)
                    #  'custom': custom shape. In this case, define
                    #  a 2D shape using .vertices
el1.dimensions = [50, 50] # Dimension in mm
el1.thickness = 4   # 4 mm thickness
                    #  1 number: electrode with 1 (gel) layer.
                    #  2 numbers: electrode with 2 layers
                    #  (gel, conductive rubber) 
                    #  3 numbers: electrode with 3 layers
                    #  (sponge, conductive rubber, sponge) 
                    #  also .dimensions_sponge must be 
                    #  then set

# second electrode
el2 = tdcs.add_electrode()
el2.channelnr = 1
el2.centre = 'C4'
el2.pos_ydir = 'C3'
el2.shape = 'rect'
el2.dimensions = [40, 40]
el2.thickness = [3.5, 1, 3.5] # sponge electrode
el2.dimensions_sponge = [50, 70]
# define where the cable is plugged in:
plug = el2.add_plug()
plug.centre = [-10, 0] # relative to electrode center
plug.shape = 'rect'
plug.dimensions = [10, 10]

# third electrode
el3= tdcs.add_electrode()
el3.channelnr = 2 #  also connected to 2nd channel
el3.centre = 'Oz'
el3.pos_ydir = 'Pz'
el3.shape = 'custom'
el3.vertices = [[ 0, 20],
                [-40, 40],
                [-20, 0],
                [-40,-40],
                [ 0, -20],
                [40, -40],
                [20, 0],
                [40, 40]]
el3.thickness = [3, 1];
# add a 15 mm hole in the center:
hole = el3.add_hole()
hole.centre = 'Oz' 
hole.shape = 'ellipse'
hole.dimensions = [15, 15]

run_simnibs(S)
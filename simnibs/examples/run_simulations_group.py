import os
from simnibs import sim_struct, run_simnibs

# Set the subjects
subjects = ['sub01', 'sub09', 'sub10', 'sub12', 'sub15']

# Set a TDCSLIST structure with the simulation set-up
tdcslist = sim_struct.TDCSLIST()
tdcslist.currents = [0.001, -0.001]

anode = tdcslist.add_electrode()
anode.channelnr = 1
anode.centre = 'C3'
anode.pos_ydir = 'C1'
anode.shape = 'rect'
anode.dimensions = [50, 50]
anode.thickness = 4


cathode = tdcslist.add_electrode()
cathode.channelnr = 2
cathode.centre = 'AF4'
cathode.pos_ydir = 'F6'
cathode.shape = 'rect'
cathode.dimensions = [50, 70]
cathode.thickness = 4


# Run the simulation in each subject
for sub in subjects:
    # ALWAYS create a new SESSION when changing subjects
    s = sim_struct.SESSION()
    s.map_to_fsavg = True
    s.fnamehead = os.path.join(sub, sub + '.msh')
    s.pathfem = os.path.join(sub, 'bipolar')
    # Add the tdcslist we defined above
    s.add_poslist(tdcslist)
    # Run the sumulation
    run_simnibs(s)

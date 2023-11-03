"""Script for simulating the TACS challenge montage.

see @TACSchallenge on Twitter

A. Thielscher, 2020
"""
from copy import deepcopy
from simnibs import sim_struct, run_simnibs


### general Information
S = sim_struct.SESSION()
S.subpath = 'm2m_ernie'  # head mesh
S.pathfem = 'TACSchallenge'  # Directory for the simulation
S.map_to_surf=True
S.open_in_gmsh=True # open results once they are ready
                    # (set to False if you are annoyed by the popup windows)

### Define Condition A
ConditionA = S.add_tdcslist()

### Set current flow though each channel:
### 2mA peak-to-peak --> 1 mA baseline-to-peak
### The third entry is a "pseudochannel" to which electrodes 
### are connected that are not used in the simulated condition
ConditionA.currents = [0.001, -0.001, 0.0]

# Define occipital electrodes
O1 = ConditionA.add_electrode()
O1.channelnr = 1  # Connect the electrode to the first channel
O1.centre = 'O1'  # Place it over O1
O1.shape = 'ellipse'  # Elliptical / circular shape
O1.dimensions = [25, 25]  # Electrode diameter (in mm)
O1.thickness = [1.5, 1.5]  # 1.5 mm gel and 1.5 mm rubber 

ConditionA.electrode.append(deepcopy(O1))
ConditionA.electrode[-1].centre = 'Oz'

ConditionA.electrode.append(deepcopy(O1))
ConditionA.electrode[-1].centre = 'O2'

# Define return electrodes
CP5 = ConditionA.add_electrode()
CP5.channelnr = 2;
CP5.centre = 'CP5';
CP5.shape = 'ellipse';
CP5.dimensions = [50, 50];
CP5.thickness = [1.5, 1.5];

ConditionA.electrode.append(deepcopy(CP5))
ConditionA.electrode[-1].centre = 'CPz'

ConditionA.electrode.append(deepcopy(CP5))
ConditionA.electrode[-1].centre = 'CP6'

# Define retinal control electrodes
Retcontr = ConditionA.add_electrode()
Retcontr.channelnr = 3; # connect to the "pseudochannel"
Retcontr.centre = [35.4, 106.0, -40.6]; # coordinates determined in simnibs_gui
Retcontr.shape = 'ellipse';
Retcontr.dimensions = [25, 25];
Retcontr.thickness = [1.5, 1.5];

ConditionA.electrode.append(deepcopy(Retcontr))
ConditionA.electrode[-1].centre = [-33.2, 108, -40.6]; # coordinates determined in simnibs_gui


# Define Condition B
S.poslists.append(deepcopy(ConditionA))
ConditionB=S.poslists[-1]

ConditionB.electrode[0].channelnr = 3;
ConditionB.electrode[1].channelnr = 3;
ConditionB.electrode[2].channelnr = 3;
ConditionB.electrode[6].channelnr = 1;
ConditionB.electrode[7].channelnr = 1;


# Define Condition C
S.poslists.append(deepcopy(ConditionA))
ConditionC=S.poslists[-1]

ConditionC.electrode[0].channelnr = 2;
ConditionC.electrode[1].channelnr = 1;
ConditionC.electrode[2].channelnr = 2;
ConditionC.electrode[3].channelnr = 3;
ConditionC.electrode[4].channelnr = 3;
ConditionC.electrode[5].channelnr = 3;
ConditionC.electrode[6].channelnr = 3;
ConditionC.electrode[7].channelnr = 3;

run_simnibs(S)

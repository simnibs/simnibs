from simnibs import sim_struct, run_simnibs

''' General Settings '''
# Initalize a session
s = sim_struct.SESSION()
# Name of head mesh
s.fnamehead = 'ernie.msh'
# Output folder
s.pathfem = 'tutorial/'

''' TMS simulations '''
# Initialize a list of TMS simulations
tmslist = s.add_tmslist()
# Select coil
tmslist.fnamecoil = os.path.join('legacy_and_other','Magstim_70mm_Fig8.ccd')


''' First coil position '''
# Initialize a coil position
pos = tmslist.add_position()
# Select coil centre
pos.centre = 'C1'
# Select coil direction
pos.pos_ydir = 'CP1'

''' Second coil position '''
# Add another position
pos_superior = tmslist.add_position()
# Centred at C1
pos_superior.centre = 'C1'
# Pointing towards Cz
pos_superior.pos_ydir = 'Cz'

''' tDCS Simulation '''
tdcslist = s.add_tdcslist()
# Set currents
tdcslist.currents = [-1e-3, 1e-3]

''' Define cathode '''
# Initialize the cathode
cathode = tdcslist.add_electrode()
# Connect electrode to first channel (-1e-3 mA, cathode)
cathode.channelnr = 1
# Electrode dimension
cathode.dimensions = [50, 70]
# Rectangular shape
cathode.shape = 'rect'
# 5mm thickness
cathode.thickness = 5
# Electrode Position
cathode.centre = 'C3'
# Electrode direction
cathode.pos_ydir = 'Cz'

''' Define anode '''
# Add another electrode
anode = tdcslist.add_electrode()
# Assign it to the second channel
anode.channelnr = 2
# Electrode diameter
anode.dimensions = [30, 30]
# Electrode shape
anode.shape = 'ellipse'
# 5mm thickness
anode.thickness = 5
# Electrode position
anode.centre = 'C4'

''' Run Simulations '''
run_simnibs(s)

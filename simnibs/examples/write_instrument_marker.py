from simnibs.utils.nnav import get_m_localite, write_tms_navigator_im
from simnibs.simulation.sim_struct import POSITION
import numpy as np

pos = POSITION()
pos.centre = [1, 1, 1]
pos.pos_ydir = [1, 0, 0]
pos.distance = 1.
pos.calc_matsimnibs('ernie.msh')

simnibs2localite = get_m_localite('ernie.nii')
localite2simnibs = np.linalg.inv(simnibs2localite)

# construct flip matrix
m_flip = np.array([[0, 0, 1, 0],
                   [0, -1, 0, 0],
                   [1, 0, 0, 0],
                   [0, 0, 0, 1]])

# transform coil position from simnibs to neuronavigation space
matlocalite = np.dot(np.dot(localite2simnibs,
                            pos.matsimnibs.transpose([2, 0, 1])
                            ).transpose([1, 0, 2]),
                     m_flip).transpose([1, 2, 0])

# write isntrument marker file
write_tms_navigator_im(matlocalite, 'instrument_marker.xml')

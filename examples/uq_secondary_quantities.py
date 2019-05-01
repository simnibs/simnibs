import numpy as np
from simnibs import mesh_io
from simnibs.simulation import gpc


fn_hdf5 = 'tdcs_gpc/ernie_TDCS_1_gpc.hdf5'
# Read the regression object from the HDF5 file
regression = gpc.gPC_regression.read_hdf5(fn_hdf5)
# Read the mesh ROI from the HDF5 file
mesh = mesh_io.Msh.read_hdf5(fn_hdf5, 'mesh_roi')
# Define the function to be used for the expansion

def percentile_99(Es):
    # The function will receive the electric field in a format
    # N_sims x N_elm x 3
    # for each simulation, we calculate the 99th percentile
    # using mesh_io
    prc = np.zeros(Es.shape[0])
    for i, E in enumerate(Es):
        prc[i] = mesh_io.ElementData(E, mesh=mesh).get_percentiles(99)[0]

    return prc

# Calculate the gPC coefficients
gpc_coeffs = regression.expand_quantity(percentile_99)

print("99th Percentile")
print("Mean Value: ", regression.mean(gpc_coeffs))
print("Standard Deviation: ", regression.std(gpc_coeffs))

# Draw 1000 samples for the 99th percentile
samples = regression.MC_sampling(gpc_coeffs, 1000)

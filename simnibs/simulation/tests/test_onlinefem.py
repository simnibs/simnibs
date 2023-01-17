import h5py
import simnibs

# Test OnlineFEM class with TMS
########################################################################################################################
fn_mesh = "/data/pt_01756/probands/15484.08/mesh/charm_beta_coarse/m2m_15484.08/15484.08.msh"
fn_coil = "/data/pt_01756/coils/Magventure MCF-B65/Medtronic_MCF_B65_REF.nii.gz"
fn_matsimnibs = "/data/pt_01756/probands/15484.08/exp/reg_isi_05/mesh_charm_beta_coarse/matsimnibs.hdf5"
fn_roi = "/data/pt_01756/probands/15484.08/mesh/charm_beta_coarse/roi/midlayer_m1s1pmd/geo.hdf5"
fn_results = "/data/pt_01756/studies/ttf/OnlineFEM/results_TMS"
solver_options = "pardiso"
anisotropy_type = "scalar"

# load mesh
print("Loading mesh ...")
mesh = simnibs.read_msh(fn_mesh)

# setup roi
print("Initializing ROI ...")

with h5py.File(fn_roi, "r") as f:
    points = f["mesh/nodes/node_coord"][:]
    con = f["mesh/elm/triangle_number_list"][:]

roi = simnibs.RegionOfInterest(points=points, con=con, msh=mesh)

# load coil positions (in simnibs space)
print("Loading coil positions ...")
with h5py.File(fn_matsimnibs, "r") as f:
    matsimnibs = f["matsimnibs"][:]

# initialize OnlineFEM
onlinefem = simnibs.OnlineFEM(mesh=mesh, method="TMS", roi=roi, anisotropy_type=anisotropy_type,
                              solver_options=solver_options, fn_results=fn_results, useElements=True,
                              grid_spacing=None, order=1, fn_coil=fn_coil, didtmax=1, dataType=0)

e = onlinefem.update_field(matsimnibs=matsimnibs)

# Test OnlineFEM class with TES
########################################################################################################################



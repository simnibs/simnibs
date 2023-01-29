import h5py
import simnibs

# Test OnlineFEM class with TMS
########################################################################################################################
# fn_mesh = "/data/pt_01756/probands/15484.08/mesh/charm_beta_coarse/m2m_15484.08/15484.08.msh"
fn_mesh = "/home/kporzig/tmp/charm_beta_coarse/m2m_15484.08/15484.08.msh"
# fn_coil = "/data/pt_01756/coils/Magventure MCF-B65/Medtronic_MCF_B65_REF.nii.gz"
fn_coil = "/home/kporzig/tmp/Medtronic_MCF_B65_REF.nii.gz"
# fn_matsimnibs = "/data/pt_01756/probands/15484.08/exp/reg_isi_05/mesh_charm_beta_coarse/matsimnibs.hdf5"
fn_matsimnibs = "/home/kporzig/tmp/matsimnibs.hdf5"
# fn_roi = "/data/pt_01756/probands/15484.08/mesh/charm_beta_coarse/roi/midlayer_m1s1pmd/geo.hdf5"
fn_roi = "/home/kporzig/tmp/charm_beta_coarse/roi/midlayer_m1s1pmd/geo.hdf5"
# fn_results = "/data/pt_01756/studies/ttf/OnlineFEM/results_TMS"
fn_results = "/home/kporzig/tmp/OnlineFEM/e_new_iso_nodes.hdf5"

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

p1_tri = points[con[:, 0], :]
p2_tri = points[con[:, 1], :]
p3_tri = points[con[:, 2], :]

triangles_center = 1.0 / 3 * (p1_tri + p2_tri + p3_tri)

roi = simnibs.RegionOfInterest(points=triangles_center, mesh=mesh)

# load coil positions (in simnibs space)
print("Loading coil positions ...")
with h5py.File(fn_matsimnibs, "r") as f:
    matsimnibs = f["matsimnibs"][:, :, :2]

# initialize OnlineFEM
onlinefem = simnibs.OnlineFEM(mesh=mesh, method="TMS", roi=roi, anisotropy_type=anisotropy_type,
                              solver_options=solver_options, fn_results=fn_results, useElements=False,
                              fn_coil=fn_coil, dataType=0)

e = onlinefem.update_field(matsimnibs=matsimnibs, didt=1e6)

# Test OnlineFEM class with TES
########################################################################################################################



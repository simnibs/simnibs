# all of this is from the memoslap project and changed to work here
import nibabel as nib
import numpy as np
import os

from simnibs import mesh_io, opt_struct
from simnibs.utils.file_finder import SubjectFiles

from simnibs.utils import transformations
from simnibs.utils.file_finder import get_reference_surf
    

def _convert_fsavg_mask(fn_mask_fsspace, hemi, subpath):
    ''' convert mask roi from fsaverage to individual space and get positions
    '''
    assert hemi in ['lh', 'rh']

    subject_files = SubjectFiles(subpath=subpath)
    fn_sphere = get_reference_surf(hemi, 'sphere')
    fn_reg = subject_files.get_surface(hemi, 'sphere_reg')
    fn_central = subject_files.get_surface(hemi, 'central')

    surf_sphere = mesh_io.read_gifti_surface(fn_sphere)
    try:
        idx = nib.freesurfer.io.read_label(fn_mask_fsspace)
        idx_mask = np.zeros(surf_sphere.nodes.nr, dtype=np.float32)
        idx_mask[idx] = 1.
    except:
        idx_mask = nib.freesurfer.io.read_morph_data(fn_mask_fsspace)

    morph = transformations.SurfaceMorph(surf_sphere, 
                                            mesh_io.read_gifti_surface(fn_reg), 
                                            method="nearest")
    idx_mask = morph.transform(idx_mask) > 0.0001
    idx_mask = idx_mask > 0.0001
    gm_surf =  mesh_io.read_gifti_surface(fn_central)
    central_pos = gm_surf.nodes[idx_mask]

    return idx_mask, central_pos



def _convert_MNImask(fn_mask, m, subpath):
    ''' maps a mask from MNI space to mesh nodes
    '''
    subj_files = SubjectFiles(subpath = subpath)     
    # convert to subject space
    roi_buffer, roi_affine = _map_roi(subj_files, fn_mask)
    roi_buffer = roi_buffer.astype(np.uint16)
    # map on nodes
    nd = mesh_io.NodeData.from_data_grid(m, roi_buffer, roi_affine)
    nd.value = nd.value > 0
    m.add_node_field(nd,'mask')
    mask_pos = m.nodes.node_coord[nd.value, :]
    return m, mask_pos


def get_central_gm_with_mask(subpath, hemi, fn_mask_fsspace,
                             mask_type='curv', add_cerebellum = False):
    """ load bihemispheric GM and add mask as node data
        simnibs4: also cerebellum GM can be added
    
    Parameters
    ----------
    subpath : string
        m2m-folder
    hemi : string
        Defines on which hemisphere the mask is ('lh', 'rh' or 'cereb').
    fn_mask: string
        Path to the mask file
    mask_type : string, optional
        Indicates the type of mask ('curv' for masks in fsaverage space, 
        'mnivol' for masks in MNI space). The default is 'curv'.
    add_cerebellum : bool, optional
        whether to add the cerebellum central gm surface to the mesh. 
        The default is False.

    Returns
    -------
    m_surf : simnibs.mesh_io.Msh
        Central gm surfaces (labels: lh 1001, rh 1002, cerebellum 1003).
        The mask is added as nodedata field (get via m_surf.field['mask']).

    """

    assert hemi in ['lh', 'rh', 'cereb']
        
    subject_files = SubjectFiles(subpath=subpath)

    fn_cereb = None
    fn_lh_central = subject_files.get_surface('lh', 'central')
    fn_rh_central = subject_files.get_surface('rh', 'central')
    
    if add_cerebellum:
        fn_cereb = os.path.join(subject_files.surface_folder, 
                                'cerebellum.central.gii')
        if not os.path.exists(fn_cereb):
            raise FileNotFoundError(fn_cereb)

    m_surf = mesh_io.read_gifti_surface(fn_lh_central)
    m_surf.elm.tag1 = 1001 * np.ones(m_surf.elm.nr, dtype=int)
    nr_nodes_lh = m_surf.nodes.nr

    m_tmp = mesh_io.read_gifti_surface(fn_rh_central)
    m_tmp.elm.tag1 = 1002 * np.ones(m_tmp.elm.nr, dtype=int)
    nr_nodes_rh = m_tmp.nodes.nr
    m_surf = m_surf.join_mesh(m_tmp)
    
    if add_cerebellum:
        m_tmp = mesh_io.read_gifti_surface(fn_cereb)
        m_tmp.elm.tag1 = 1003 * np.ones(m_tmp.elm.nr, dtype=int)
        nr_nodes_cereb = m_tmp.nodes.nr
        m_surf = m_surf.join_mesh(m_tmp)
    
    m_surf.elm.tag2[:] = m_surf.elm.tag1

    if mask_type == 'mnivol':
        m_surf, mask_pos = _convert_MNImask(fn_mask_fsspace, m_surf, subpath)

    elif mask_type == 'curv':
        # requires surface registration to fsaverage, thus only for lh and rh
        idx_mask, _ = _convert_fsavg_mask(fn_mask_fsspace, hemi, subpath)
        
        if hemi == 'lh':
            idx_lh = idx_mask
        else:
            idx_lh = np.zeros(nr_nodes_lh,dtype=bool)
            
        if hemi == 'rh':
            idx_rh = idx_mask
        else:
            idx_rh = np.zeros(nr_nodes_rh,dtype=bool)
            
        if add_cerebellum:
            idx_cereb = np.zeros(nr_nodes_cereb,dtype=bool)
        else:
            idx_cereb = []
                
        nd=mesh_io.NodeData( np.hstack((idx_lh, idx_rh, idx_cereb)) )
        m_surf.add_node_field(nd,'mask')
        
    elif mask_type == '':
        nd=mesh_io.NodeData( np.zeros(m_surf.nodes.nr, dtype=bool) )
        m_surf.add_node_field(nd,'mask')
        
    else:
        raise ValueError(f"unknown mask_type: {mask_type}")
    
    return m_surf


subject_path = 'm2m_ernie'
hemisphere = 'lh'
mask_path = 'P1_LH_M1_control'
central_gm_with_mask = get_central_gm_with_mask(
    subpath=subject_path,
    hemi=hemisphere,
    fn_mask_fsspace=mask_path
)

coordinates_of_nodes_in_mask = central_gm_with_mask.nodes.node_coord[central_gm_with_mask.field['mask'].value > 0, :]

# Initialize structure
tms_opt = opt_struct.TMSoptimize()
tms_opt.open_in_gmsh = False # no
# Subject folder
tms_opt.subpath = subject_path
# Select output folder
tms_opt.pathfem = 'tms_optimization_region/'
# Select the coil model
tms_opt.fnamecoil = os.path.join('legacy_and_other','Magstim_70mm_Fig8.ccd')
# Select targets for the optimization. Every target will be the center of a sphere with radius tms.opt.target_size. The union of every sphere will form the target ROI
tms_opt.multiple_targets = coordinates_of_nodes_in_mask
# Size of every region in multiple targets
tms_opt.target_size = 3.0
tms_opt.angle_resolution = 120 # very low resolution for faster testing
tms_opt.spatial_resolution = 10
# Optional: Use the MKL PARDISO solver
# Will make the simulations much faster
# but has large (approx 12GB) memory usage
tms_opt.solver_options = 'pardiso'
# Run optimization to get optimal coil position
opt_pos=tms_opt.run(save_mat=False, cpus=1)

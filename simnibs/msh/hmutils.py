'''
    Utilities for Head Modeling
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2018 Jesper D Nielsen, Axel Thielscher, Guilherme Saturnino

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''


"""Suite of utilities for generating mesh files (OFF and STL) from volume masks 
(NIFTI) and vice versa. The main functions are vol2msh and msh2vol

The following functions are included:

binarize

decouple_surfaces
deform_mesh
format_time

get_connected_components
get_element_neighbors
get_large_components
get_n_largest_components
get_node_normals
get_tetrahedra_normals
get_triangle_normals
get_vx_size

in_mesh

label_points
list2numpy

make_splits
make_surface_mesh
make_volume_mask

mesh_adjust_tetrahedra_regions
mesh_load
mesh_save

msh2vol
msh_coord2vox_idx

path2bin
prepare_vols
progress_bar
ray_trace
ray_triangle_intersect
remove_duplicates
reorient2LAS
row_wise_unique

spawn_process

verts2faces
vmesh2nifti
vol2msh
vol_mask
volume_mesh_relabel
vox_idx2msh_coord

write_nifti
write_ref_fs
write_unity_xfm


EXTERNAL PROGRAMS USED
----------------------
meshfix
Gmsh

"""

# DEPENDENCIES
# =============================================================================
# Python
from glob      import glob
from itertools import islice
import logging
import nibabel as nib
import numpy as np
import os
import shutil
import scipy.ndimage.morphology as mrph
import scipy.sparse
from scipy.ndimage.filters      import gaussian_filter
from scipy.ndimage.measurements import label
from scipy.ndimage              import map_coordinates
from scipy.spatial              import cKDTree
import shlex
import struct
from subprocess import Popen, PIPE
import sys
import textwrap
from threading import Thread
from multiprocessing import Process
try:
    from Queue import Queue, Empty
except:
    from queue import Queue, Empty
import base64
if sys.version_info >= (3, ):
    base64decode = base64.decodebytes
else:
    base64decode = base64.decodestring

# simnibs
from simnibs import SIMNIBSDIR
from simnibs.msh import mesh_io
from . import _marching_cubes_lewiner_luts as mcluts
import simnibs.cython_code._marching_cubes_lewiner_cy as _marching_cubes_lewiner_cy
from ..utils.file_finder import path2bin

logger = logging.getLogger("py.warnings")
width = 72


# =============================================================================
# VOLUME UTILITIES
# =============================================================================

def prepare_vols(fin, out_format="keep", vx_size=None, img_dims=None,
                 center=False, write_vols=True):
    """Transform a set of volumes to LAS (radiological) space. Optionally,
    resample the image to a new set of voxel sizes and image dimensions. 
    Assumes that all input volumes have the same voxel size and image 
    dimensions.
    
    PARAMETERS
    ----------
    fin : list
        List of strings (filenames) or nibabel image objects.
    out_format : str
        How to write the new volumes. Options are "standard", "keep",
        "custom", "no-conform") (default = "keep").
    vx_size : array_like
        Array specifying voxel sizes (only applies if out_format = "custom").
    img_dims : array_like
        Array specifying image dimensions (only applies if out_format = 
        "custom").
    center : bool, optional
        Whether or not to center the data. The center is assumed to be
            M^(-1)*[0,0,0,1]
        where M is the affine transformation matrix between voxel coordinates
        (indices) and physical or standard space (default = False).
    write_vols : bool, optional
        Write conforming volumes to disk (default = True).
    
    RETURNS
    ----------
    vols_conform : list
        List of nibabel image objects, i.e. the volumes in the conformed space. 
    """     
    # Check inputs
    out_format = out_format.lower()
    
    if type(fin) == str:
        fin=[fin]
        
    # if input is not iterable, make it so
    try:
        next(iter(fin))
    except TypeError:
        fin = [fin]
    
    if vx_size is None:
        vx_size = []
    if img_dims is None:
        img_dims = []
        
    vols=[]
    for f in fin:
        if type(f) == str:
            try:
                vols.append(nib.load(f))
            except:
                raise IOError(f+" does not exist!") 
        elif type(f) in [nib.nifti1.Nifti1Image,nib.nifti2.Nifti2Image,nib.nifti1.Nifti1Pair,nib.nifti2.Nifti2Pair]:
            vols.append(f)
        else:
            raise ValueError("Input volumes must be specified as either a string (path to volume) or a nibabel image object.")
    
    vols_conform = []

    for v in vols:
        hdr_orig = v.header
        fname_out = add2filename(v.get_filename(), "_conform")
        fname_in = v.get_filename()
        
        if out_format == "no-conform":
            # If the files are already in the appropriate space simply copy and
            # rename the files
            log("Copying {0} to {1}", 
                [os.path.basename(fname_in),
                 os.path.basename(fname_out)])
                 
            data = v.get_data()
            M = v.affine            
        else:
            # Reorient to LAS (radiological) space
            log("Preparing {0} from {1}", 
                [os.path.basename(fname_out),
                 os.path.basename(fname_in)])
            log("-"*width)
            
            v = reorient2LAS(v)
            
            vx_size_original = get_vx_size(v)
            
            if out_format == "keep":
                vx_size  = vx_size_original
                img_dims = np.array(v.shape) 
            elif out_format == "standard":
                vx_size  = np.array([1]*3)
                img_dims = np.array([256]*3)
            elif (out_format == "custom") & (len(vx_size) == 3) & (len(img_dims) == 3):
                vx_size = np.array(vx_size)
                img_dims = np.array(img_dims)
            elif (out_format == "custom") & (len(vx_size) == 3):
                vx_size  = np.array(vx_size).astype(float)
                img_dims = np.array(v.shape)*vx_size_original/vx_size
            elif (out_format == "custom") & (len(img_dims) == 3):
                img_dims = np.array(img_dims)
                vx_size  = np.round(vx_size_original*np.array(v.shape) / img_dims.astype(np.float))
            else:
                raise KeyError("Unrecognized output format.")
            
            img_dims_int = np.ceil(img_dims).astype(np.int)
            
            log("Output format")
            print("")
            log("{:19s} {:>30s}",("Orientation","Left-Anterior-Superior (LAS)"))
            log("{:33s} [{:3.2f} {:3.2f} {:3.2f}]",("Voxel size", vx_size[0],vx_size[1],vx_size[2]))
            log("{:33s} [{:4d} {:4d} {:4d}]",("Image dimensions", img_dims_int[0],img_dims_int[1],img_dims_int[2]))
            print("")
            
            # The affine transformation matrix of the output image
            M = np.array([[-vx_size[0],         0 ,         0 ,  vx_size[0]*img_dims_int[0]/2.],
                          [         0 , vx_size[1],         0 , -vx_size[1]*img_dims_int[1]/2.],
                          [         0 ,         0 , vx_size[2],  vx_size[2]*(-img_dims_int[2]/2.+1)],
                          [         0 ,         0 ,         0 ,  1]])
            
            # Calculate the scaling factor between the new and the original
            # voxel size in each dimension to determine how densely to sample
            # the original space.
            sampling_density = np.abs( vx_size/vx_size_original.astype(np.float) )

            # Calculate the amount (in mm) by which to offset the new grid so
            # as to center it on the original data (floor it to minimize
            # smoothing effects at the expense of sampling a little
            # incorrectly).
            FOV_offset = np.round((np.array(v.shape)*vx_size_original-img_dims*vx_size)*sampling_density / 2. )

            if np.any(FOV_offset-1e-3 > 0) or np.any(img_dims*sampling_density + FOV_offset+1e-3 < v.shape):
                log("WARNING: Sampling the original volume with the current voxel size/image dimensions will reduce FOV!", level=logging.WARNING)
 
            img = v.get_data().copy()
            img_affine = v.affine.copy()
            
            if center:
                # Move to the center of the array the point corresponding to 
                # (0,0,0) in mm space (standard space or scanner space)
                t2c = np.round( (np.array(v.shape)-1)/2. - np.linalg.inv(img_affine).dot(np.array([0,0,0,1]) )[:3]).astype(np.int)              
                if np.any(t2c):
                    for i in range(3):
                        img = np.roll(img,t2c[i],axis=i)

            # Sample the new grid        
            # get grid coordinates (in mm) with the new voxel size
            grid = np.array(np.meshgrid(*tuple(map(np.arange,[0,0,0],img_dims*sampling_density,sampling_density)),indexing="ij")).reshape((3,-1))
            
            # final FOV to sample
            grid += np.array(FOV_offset)[:,np.newaxis]

            # Interpolate values of new coordinates by sampling the original
            # data matrix. Use mode="nearest" at edges to prevent coordinates JUST
            # outside from automatically being assigned a value of zero.
            # Zero pad to avoid filling of entire rows/columns
            if np.mean(sampling_density) < 1.05 and np.mean(sampling_density) > 0.95:
                spline_order = 0
            else:
                spline_order = 5
            
            data = map_coordinates( np.pad(img,1,mode="constant",constant_values=0) , grid+1 ,order=spline_order,mode="nearest").reshape(img_dims_int)
            data[data < 0] = 0
        
        # ensure that data is 16 bit unsigned integers and rescale to range
        data = (data / float(data.max()) * (2**16-1)).astype(np.uint16)
        
        # save image using the new affine transformation matrix, M
        vout = nib.Nifti1Image(data,M)
        vout.set_qform(M)
        vout.header.set_xyzt_units(*hdr_orig.get_xyzt_units())
        vout.set_filename(fname_out)
        
        vols_conform.append(vout)
        
        if write_vols:
            nib.save(vout,vout.get_filename())
         
    return vols_conform

    
def reorient2LAS(v):
    """Reorients axes so as to bring the volume in to LAS (left-anterior-
    superior) space, i.e. radialogical convention. This reorientation is not
    robust, however, in that it relies on finding the primary direction of each
    axis and assigning it to the corresponding dimension. As such, if an axis
    is more than 45 degrees off, this will fail to correctly match an array
    axis with its corresponding dimension.
    
    PARAMETERS
    ----------
    v : nibabel image object
        The image volume.
        
    RETURNS
    -------  
    vout : nibabel image object
        reoriented nibabel image object in LAS (radiological) convention.
    """    
    qform_inv = np.linalg.inv(v.get_qform())[:3,:3]
    qform_inv_sign = np.sign(qform_inv)
    
    # Get sorting to LAS. This should work except for a few "pathological"
    # cases (e.g. two 45 deg axes)
    _, LASdim = np.where((np.abs(qform_inv) == np.max(np.abs(qform_inv),axis=0)).T)
    LASsign = qform_inv_sign[LASdim,range(3)]
    LASsign[0] *= -1

    # Apply to qform
    p = np.zeros_like(v.get_qform())
    p[-1,-1] = 1
    p[LASdim,range(3)] = 1
    
    dims = np.array(v.shape)
    dims = dims[LASdim] # dimensions of the reoriented image 

    # flip axes to LAS   
    flip = np.eye(4)
    for i in range(3):
        if LASsign[i] < 0:
            flip[i,i] = -1
            flip[i,-1] = dims[i]-1
    
    # Apply to data
    data = np.transpose(v.get_data(),LASdim) # permute axes
    if LASsign[0]<0: data = data[::-1,...]   # and flip...
    if LASsign[1]<0: data = data[:,::-1,:]
    if LASsign[2]<0: data = data[...,::-1]
    
    # Make reoriented image
    vout = nib.Nifti1Image(data,affine=(v.get_qform().dot(p)).dot(flip))
    vout.set_qform((v.get_qform().dot(p)).dot(flip))
    vout.header.set_xyzt_units(*v.header.get_xyzt_units())
    
    return vout


def vol_mask(wm_raw, gm_raw, csf_raw, bone_raw, skin_raw, air_raw, prior_tissue,
             prior_eye1, prior_eye2, prior_air, prior_spinal, prior_ventricles_lateral, out_dir, usecat=False):
    """Take raw binary tissue masks, fix minor issues
    (segmentation imperfections), and create binary tissues masks which contain
    all inner tissue types as well (i.e. gray matter mask would be gray matter
    and white matter, CSF mask would be CSF, gray matter, and white matter
    etc.).
    
    PARAMETERS
    ----------
    [tissue]_raw : str
        Filename corresponding to a binary map from a
        segmentation procedure (e.g., the c* output files from the SPM 
        segmentation, )
    out_dir : str
        Directory in which to save the cleaned masks.
    usecat: Bool
        use masks for GM, WM and ventricles created by CAT; these masks
        will not be changed during volume mask creation
    
    RETURNS
    ----------
    Save the cleaned masks to the specified directory, naming them as
    "MASK_[tissue].nii.gz"   
    """    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    odtype = np.uint8 # output data type
    
    # load volumes
    vol_hdr = nib.load(wm_raw)
    vx_size = get_vx_size(vol_hdr)
    
    wm_raw   = vol_hdr.get_data()
    gm_raw   = nib.load(gm_raw).get_data()
    csf_raw  = nib.load(csf_raw).get_data()
    bone_raw = nib.load(bone_raw).get_data() 
    skin_raw = nib.load(skin_raw).get_data()
    air_raw  = nib.load(air_raw).get_data()
    
    ptissue = nib.load(prior_tissue)
    peyes = [nib.load(prior_eye1), nib.load(prior_eye2)]
    pair = nib.load(prior_air)
    pspin = nib.load(prior_spinal)
    pven = nib.load(prior_ventricles_lateral)
    
    # generate structuring elements
    se1 = mrph.generate_binary_structure(3,1) # 6-connectivity (face) neighbors
    se2 = mrph.generate_binary_structure(3,2)    
    connectivity = np.round(1/np.average(vx_size)).astype(np.int) # scaled by voxel size
    se1mm = mrph.generate_binary_structure(3, connectivity)
    
    wm_raw = wm_raw.astype(np.int)
    gm_raw = gm_raw.astype(np.int)
    csf_raw = csf_raw.astype(np.int)
    bone_raw = bone_raw.astype(np.int)
    skin_raw = skin_raw.astype(np.int)
    air_raw = air_raw.astype(np.int)
    
    alltissue_raw = wm_raw | gm_raw | csf_raw | bone_raw | skin_raw | air_raw
    alltissue_raw = mrph.binary_fill_holes(alltissue_raw)
    
    se = se1
    
    # --------------------------
    log("Cleaning bone")
    # close bone
    clo = mrph.binary_closing(np.pad( bone_raw, 1, mode="constant",
                                     constant_values=0), se)[1:-1,1:-1,1:-1]
    diff = clo-bone_raw
    bone_raw = clo
    # remove from the other masks those voxels which the closing added to bone
    # (thereby protecting GM, WM and VENTRICLES from changes in case they 
    #  were already created based on CAT results)
    if usecat: 
        bone_raw[gm_raw>0] = 0;
        diff[gm_raw>0] = 0;
    wm_raw[diff>0] = 0; gm_raw[diff>0] = 0; csf_raw[diff>0] = 0
    skin_raw[diff>0] = 0; air_raw[diff>0] = 0
    del diff
    del clo

    # --------------------------
    log("Cleaning skin")
    # Get the largest component of the eroded skin mask and dilate. Update the
    # skin mask so that it contains only voxels which are contained in both the
    # original and the "pseudo-opened" skin mask
    dil = mrph.binary_dilation(get_n_largest_components(mrph.binary_erosion(skin_raw,se,1),se,1),se,1) # dilate only largest component
    empty = skin_raw-dil
    skin_raw = dil & skin_raw
    del dil

    # --------------------------    
    log("Cleaning CSF")
    # Remove all voxels outside of the head which are labelled as CSF.
    csf_rawi = csf_raw & (ptissue.get_data() > 0.5)
    air_raw += csf_raw - csf_rawi
    csf_raw = csf_rawi
    del csf_rawi
    
    # --------------------------
    log("Cleaning WM")
    if not usecat:   
        # add spinal cord
        # add gm inside prior at 0.2 (the spinal often seems to be partly segmented
        # as gm); add the lower part of the wm segmentation; close this mask
        wm_add = np.zeros_like(wm_raw)
        idx = np.where(pspin.get_data()>0)[2].max()+1
        wm_add[...,:idx] += get_n_largest_components(wm_raw>0,se,1).astype(np.int)[...,:idx]
        wm_add[pspin.get_data()>0.2] += (gm_raw[pspin.get_data()>0.2]).astype(np.int)
        wm_add = mrph.binary_closing(np.pad(wm_add,((0,0),(0,0),(3,0)),"constant"),se,3)[...,3:].astype(np.int)
        wm_add = get_n_largest_components(wm_add,se,1).astype(np.int)
            
        wm_add = ((wm_add-wm_raw)>0).astype(np.int)
        wm_raw += wm_add
        # remove what was added from the remaining volumes
        for vol in [gm_raw, csf_raw, bone_raw, skin_raw, air_raw]:
                vol[(vol & wm_add).astype(bool)] = 0
        del wm_add
        del idx
        
        log("Labelling unassigned voxels")
        wm_raw, gm_raw, csf_raw, bone_raw, skin_raw, air_raw = \
                smoothfill([wm_raw, gm_raw, csf_raw, bone_raw, skin_raw, air_raw], empty)
    else:
        log("Labelling unassigned voxels (protecting everything inside GM mask from CAT)")
        gm_raw, csf_raw, bone_raw, skin_raw, air_raw = \
                smoothfill([gm_raw, csf_raw, bone_raw, skin_raw, air_raw], empty, protect=[0])

    skin_raw = skin_raw & alltissue_raw # to avoid "headpost" like skin structures
    del alltissue_raw                   # filled in by smoothfill in case skin binarizaition was
                                        # not covered by air binarization

    # ----------------------------------------------------
    # Work inside out to generate outer tissue surfaces        
    log("Constructing outer tissue surfaces")
    vol = {}
    
    # --------------------------
    log("White matter")
    if usecat:
        vol["wm"] = wm_raw.astype(np.int)
    else:
        wm_fill = mrph.binary_fill_holes(wm_raw)    
        vol["wm"] = get_n_largest_components(wm_fill,se,1).astype(np.int)
        del wm_fill 
    del wm_raw
    
    # --------------------------
    log("Gray matter")
    if usecat:
        vol["gm"] = gm_raw.astype(np.int)
    else:  
        gm_fill = mrph.binary_fill_holes(gm_raw | vol["wm"])
        vol["gm"] = get_n_largest_components(gm_fill,se,1).astype(np.int)
        del gm_fill
        
        # ensure a layer of GM around WM by removing from WM
        vol["wm"] = vol["wm"] & mrph.binary_erosion(vol["gm"],se)
        vol["wm"] = get_n_largest_components(vol["wm"],se,1).astype(np.int)    
    del gm_raw
    
    # --------------------------
    log("CSF and eyes")
    csf_fill = mrph.binary_fill_holes(csf_raw | mrph.binary_dilation(vol["gm"],se1mm))
    csf_fill2 = mrph.binary_dilation(get_n_largest_components(mrph.binary_erosion(csf_fill,se,2),se,1),se,2) 
    vol["csf"] = csf_fill2 & csf_fill
    vol["eyes"] = get_n_largest_components(csf_fill ^ vol["csf"],se,2)
    del csf_fill
    del csf_fill2
    
    # check eyes against prior; if nothing (or very little) add prior
    # check inside prior
    for eye in peyes:
        eyepri = eye.get_data()>0.5
        eyepercent = (vol["eyes"] & eyepri).sum() / float(eyepri.sum())
        if eyepercent < 0.25:
            vol["eyes"] += eyepri
    vol["eyes"] = vol["eyes"] > 0
    # check outside prior
    vol["eyes"][(peyes[0].get_data() < 0.05) & (peyes[1].get_data() < 0.05)] = 0
           
    smooth = 1./vx_size
    vol["csf"]  = mrph.binary_fill_holes(gaussian_filter(vol["csf"].astype(np.float),smooth) > 0.5).astype(np.int)
    vol["eyes"] = mrph.binary_fill_holes(gaussian_filter(vol["eyes"].astype(np.float),smooth) > 0.5).astype(np.int)  
    
    # --------------------------
    log("Ventricles")
    if usecat:
        vol["ventricles"] = pven.get_data().astype(np.int)
    else:
        # ventricles, CSF enclosed by GM
        # remove the largest component (should be main CSF)
        remain = get_large_components(csf_raw & ~get_n_largest_components(csf_raw,se,1),se,100)
        ven = get_n_largest_components(mrph.binary_erosion(pven.get_data(),se,1),se,1) & remain
        ven = mrph.binary_fill_holes(ven,se)

        # clean up by removing small clusters of voxels
        ven = get_large_components(ven, se, 100/np.mean(vx_size)).astype(np.int)
    
        # if overlap with some component in 'remain' then include the entire component
        remain_label = label(remain,se)[0] # perhaps use se2
        labels2use = np.unique(remain_label[ven>0])
    
        if 0 in labels2use:
                labels2use = np.delete(labels2use, 0)
        
        if len(labels2use)>0:
                for i in labels2use:
                        ven[remain_label==i] += remain[remain_label==i]
                ven = ven>0
    
        if ven.sum() > 0:
                vol["ventricles"] = ven
                # ensure by digging into WM if necessary
                vol["gm"] = (vol["gm"] | mrph.binary_dilation(vol["ventricles"],se1mm)).astype(np.int)
                vol["wm"] = vol["wm"] & mrph.binary_erosion(vol["gm"]) & ~mrph.binary_dilation(vol["ventricles"],se1mm)
                vol["wm"] = get_n_largest_components(vol["wm"],se,1).astype(np.int)
        del remain
        del ven
        del remain_label
        del labels2use
    del csf_raw
    
    # --------------------------
    log("Bone")
    bone_fill = mrph.binary_fill_holes(bone_raw | mrph.binary_dilation(vol["csf"],se1mm)) # this will also fill air cavities!
    bone_fill = mrph.binary_fill_holes(gaussian_filter(bone_fill.astype(np.float),smooth) > 0.2)    
    vol["bone"] = get_n_largest_components(bone_fill,se,1)
    vol["bone"] = mrph.binary_fill_holes(gaussian_filter(vol["bone"].astype(np.float),smooth) > 0.8).astype(np.int)
    # ensure small distance between bone and eyes
    vol["bone"][mrph.binary_dilation(vol["eyes"],se,2)] = 0
    del bone_fill
    del bone_raw

    # --------------------------
    log("Skin")
    skin_fill = mrph.binary_fill_holes(skin_raw | mrph.binary_dilation(vol["bone"], se1mm) | vol["eyes"]) # this will also fill air cavities!    
    vol["skin"] = get_n_largest_components(skin_fill,se,1)
    vol["skin"] = mrph.binary_fill_holes(gaussian_filter(vol["skin"].astype(np.float),smooth) > 0.8).astype(np.int)
    vol["skin"][mrph.binary_dilation(vol["eyes"],se,2)] = 1
    del skin_fill
    del skin_raw

    # --------------------------    
    log("Air")
    # update and clean the air mask
    vol["air"] = mrph.binary_fill_holes(gaussian_filter(((pair.get_data() > 0.1) & mrph.binary_erosion(vol["skin"],se) & air_raw).astype(np.float),smooth) > 0.3).astype(np.int)    
    
    # dilate eyes by two voxels (there should not be any air nearby) and remove
    # any air voxels within this mask; this should then be relabelled to skin
    vol["air"][mrph.binary_dilation(vol["eyes"],se,2)] = 0
    
    # clean up by removing small components (< 250 mm3)
    vol["air"] = get_large_components(vol["air"], se,threshold=250/np.mean(vx_size)).astype(np.int)
    
    # remove overlap between air outside bone and bone by removing from bone.
    # remove if more than 80 % of voxels are outside bone but there is overlap
    air_label = label(vol["air"])[0]
    fullcount=np.unique(air_label,return_counts=True)[1]
    airoutside = (vol["air"] - vol["bone"] > 0)
    for comp,count in zip(*np.unique(air_label[airoutside],return_counts=True)):  
        #print count/float(fullcount[comp])
        if (count != fullcount[comp]) & (count > fullcount[comp]*0.6):
            vol["bone"][air_label == comp] = 0
    del air_label
    del fullcount
    del airoutside

    # --------------------------    
    log("Ensuring tissue masks are properly contained within each other")
    se = se1
    if usecat:
        pairs = [["gm", "csf"], ["csf", "bone"], ["bone", "skin"]]
    else:
        pairs = [["wm","gm"], ["gm", "csf"], ["csf", "bone"], ["bone","skin"]]

    # Ensure that all voxels in inner are also in outer (inside out)
    for inner,outer in pairs:
        vol[outer] = (vol[inner] | vol[outer]).astype(np.int)
    
    # Ensure the required distance (given by the structuring element of the
    # erosion) between inner and outer by cutting from inner (outside in)
    if usecat:
        pairs = [["bone", "skin"]]
    for inner,outer in reversed(pairs):
        if inner is "wm" and outer is "gm":
            # since we do erosion/expansion of GM with se2 later, decouple with
            # this
            outer_ero = mrph.binary_erosion(vol[outer], se2)
        else:
            outer_ero = mrph.binary_erosion(vol[outer], se1)
        remove = (vol[inner]-outer_ero > 0).astype(np.int)
        del outer_ero
        if np.any(remove):
            log("# voxels removed from {:4s} : {:7d}",[inner,remove.sum()])
            vol[inner] -= remove
            vol[inner] = get_n_largest_components(vol[inner], se, 1).astype(np.int)

    if usecat:
        # make sure that bone is still at least at large as CSF to avoid side-effects
        vol["bone"] = (vol["csf"] | vol["bone"]).astype(np.int)

    # Remove possible overlaps
    remove = vol["eyes"] & vol["bone"]
    vol["eyes"] -= remove.astype(np.int)
  
    remove  = (vol["csf"] & vol["air"]) | (vol["eyes"] & vol["air"])
    vol["air"] -= remove.astype(np.int)
    del remove

    # --------------------------    
    log("Saving individual tissue masks")
    keys = vol.keys()
    if usecat: 
        keys = list(set(keys) ^ set(["ventricles", "gm", "wm"]))

    for k in keys:
        write_nifti(vol[k], os.path.join(out_dir,"MASK_"+k.upper()), vol_hdr,
                    dtype=odtype)


    # --------------------------    
    log("Saving combined tissue mask")
    # wm, eyes, air are are fine; disentangle the remaining tissues (the order
    # is important here!)
    vol["skin"] = ((vol["skin"]-vol["bone"]-vol["eyes"]-vol["air"]) > 0).astype(np.int)
    vol["bone"] = ((vol["bone"]-vol["csf"]-vol["air"]) > 0).astype(np.int)
    vol["csf"] -= vol["gm"]
    vol["gm"]  -= vol["wm"]
    if "ventricles" in vol:
        vol["gm"] -= vol["ventricles"]
    
    all_tissues = 1*vol["wm"]+2*vol["gm"]+3*vol["csf"]+4*vol["bone"]+5*vol["skin"]+6*vol["air"]+7*vol["eyes"]
    if "ventricles" in vol:
        all_tissues += 8*vol["ventricles"]

    write_nifti(all_tissues, os.path.join(out_dir,"control"), vol_hdr,
                dtype=odtype)


def smoothfill(vols, unassign, protect=None):
    """Fill using Gaussian smoothing.
    
    """
    vols = [v.astype(np.float) for v in vols]
    unassign = unassign.copy()
    if protect is None:
        protect = []
    
    try:
        next(iter(protect))
    except StopIteration:
        pass
    except TypeError:
        protect = [protect]
    noprotect = list(set(range(len(vols)))-set(protect))
    
    while unassign.sum()>0: # for as long as there are unassigned voxels
        for i in noprotect:
            cs = gaussian_filter(vols[i].astype(np.float),1)
            cs[vols[i]==1] = 1
            vols[i] = cs
            
        vols = binarize(vols, return_empty=True)
        unassign = vols.pop()
    
    return vols

   
def binarize(vols, return_empty=False):
    """Binarize a list of input volumes by finding the maximum posterior 
    probability of each and assigning the voxel to this volume.
    
    PARAMETERS
    ----------
    vols : list
        List of filenames or nibabel image objects (describing probabilities of
        different tissue types).
    return_empty : bool
        If true, return an array containing all voxels which are not assigned
        to either of the other volumes (default: False)

    RETURNS
    ----------
    bin_vols : list
        List of ndarrays describing binarized versions of the input volumes.
    unassign : ndarray
        Array containing any unassigned voxels.
    """    
    # if filenames are provided, load data
    volsi = [None]*len(vols)
    for i in range(len(vols)):
        if isinstance(vols[i], str) and os.path.isfile(vols[i]):
            volsi[i] = nib.load(vols[i]).get_data()
        elif type(vols[i]) in [nib.Nifti1Image, nib.Nifti2Image]:
            volsi[i] = vols[i].get_data()
        else:
            # assume numpy array
            volsi[i] = vols[i]

    # Concatenate arrays/images
    imgs = np.concatenate(tuple([v[...,np.newaxis] for v in volsi]), axis=3)
    imgs = np.concatenate((np.zeros_like(volsi[0])[...,np.newaxis], imgs),
                          axis=3)
    
    # Find max indices
    max_idx = np.argmax(imgs, axis=3)
    
    # Binarize. Here vols_bin[0] contain voxels not assigned to any other
    # volume
    vols_bin=[]
    for i in range(imgs.shape[-1]):
        vols_bin.append(max_idx == i)
    #vols_bin = [vb.astype(np.uint8) for vb in vols_bin]

    if return_empty:
        return vols_bin[1:]+[vols_bin[0]]
    else:
        return vols_bin[1:]

    
def get_n_largest_components(vol, se, n, return_sizes=False):
    """Get the n largest components from a volume.
    
    PARAMETERS
    ----------
    vol : ndarray
        Image volume. A dimX x dimY x dimZ array containing the image data.
    se : ndarray
        Structuring element to use when detecting components (i.e. setting the
        connectivity which defines a component).
    n : int
        Number of (largest) components to retain.
    return_sizes : bool, optional
        Whether or not to also return the sizes (in voxels) of each component
        that was retained (default = False).
    
    RETURNS
    ----------
    Components : ndarray
        Binary dimX x dimY x dimZ array where entries corresponding to retained
        components are True and the remaining entries are False.
    """    
    vol_lbl = label(vol,se)[0]
    labels,region_size = np.unique(vol_lbl,return_counts=True)
    labels = labels[1:]             # disregard background (label=0)
    region_size = region_size[1:]   #
    labels = labels[ np.argsort(region_size)[::-1]]
    components = np.any(np.array([vol_lbl == i for i in labels[:n]]),
                        axis=0)
                        
    # if no components are found, components will be reduced to false. Replace
    # by array of appropriate size
    if components.sum() == 0:
        components = np.zeros_like(vol, dtype=bool)
        
    if return_sizes:
        return components,region_size[labels[:n]]
    else:
        return components

    
def get_large_components(vol, se, threshold):
    """Get the components larger than a given threshold from a volume.
    
    PARAMETERS
    ----------
    vol : ndarray
        Image volume. A dimX x dimY x dimZ array containing the image data.
    se : ndarray
        Structuring element to use when detecting components (i.e. setting the
        connectivity which defines a component).
    threshold : float
        Only components (strictly) larger than this value (# of voxels) are
        retained.
        
    RETURNS
    ----------
    components : ndarray
        Binary dimX x dimY x dimZ array where entries corresponding to retained
        components are True and the remaining entries are False.
    """    
    vol_lbl = label(vol,se)[0]
    labels, region_size = np.unique(vol_lbl,return_counts=True)
    labels = labels[1:] # disregard background (label=0)
    region_size = region_size[1:]
    components = np.any(np.array([vol_lbl == i for i in labels[region_size > threshold]]),
                        axis=0)
                        
    # if no components are found, components will be reduced to false. Replace
    # by array of appropriate size
    if components.sum() == 0:
        components = np.zeros_like(vol, dtype=bool)
        
    return components


def get_vx_size(vol):
    """Return the voxel size of a volume.
    
    PARAMETERS
    ----------
    vol : nibabel image object
        Image volume.
    
    RETURNS
    ----------
    vx_size : numpy.ndarray
        Voxel size in all three dimensions.
    """    
    try:
        vx_size = np.linalg.norm(vol.affine[:3,:3], ord=2, axis=0)
    except AttributeError:
        raise IOError("Input must be a nibabel image object.")
        
    return vx_size



def check_volumes_meshed(m2m_folder, threshold=0.5):
    ''' Checks if the volumes have been correctly meshed

    Parameters
    -------------
    m2m_folder: str
        Path to m2m_{subID}
    threshold: float
        How much smaller is the meshed volume allowed to be
    '''
    from .. import file_finder
    files = file_finder.SubjectFiles(subpath=m2m_folder)
    after_masks = nib.load(files.final_contr).get_data()
    before_masks = nib.load(files.masks_contr).get_data()
    after_masks[(after_masks == 3) * (before_masks == 8)] = 8
    before_masks[before_masks == 6] = 0
    before_masks[before_masks == 7] = 6
    tissues = np.union1d(
        np.unique(before_masks),
        np.unique(after_masks))
    tissues = tissues[tissues != 0]
    for i in tissues:
        mask_b = before_masks[before_masks == i]
        mask_a = after_masks[after_masks == i]
        if np.sum(mask_a) < threshold * np.sum(mask_b):
            return False
    return True

# =============================================================================
# SURFACE MESH UTILITIES
# =============================================================================

def vol2msh(fnames,file_format="off",binary=True):
    """Extracts the surface of a volume and saves the resulting mesh to disk.
    Uses the marching cubes algorithm.
    
    PARAMETERS
    ----------
    fnames : list (str)
        List of files to process. If only a single file, this may be a string.
        Files may be binary or contain different labels for depending on tissue
        type.
    file_format : str
        Format in which to save mesh file. Choose OFF or STL (default = "off").
    binary : bool
        Only used together with file_format = "stl". Whether to save file as
        binary (or ascii) (default = True).
    
    RETURNS
    ----------
    The extracted mesh(es) are saved to disk appended by "_mesh" or
    "_regionN_mesh" depending on the number of regions (i.e. unique values) in
    the input volume.
    """
    return_files = []
    
    if isinstance(fnames,str):
        fnames = [fnames]
        
    for P in fnames:
        V = nib.load(P)
        Y = V.get_data()
        
        # To make the marching cubes algorithm close all surfaces, pad all
        # dimensions with one zero
        Y = np.pad(Y,1,mode="constant")
        
        # get labels of different tissues (zero assumed to be background)
        regions = np.unique(Y)[1:]
        for r in regions:
            # Extract the surface
            # Consider allow_degenerate=False and gradient_direction='ascent'
            vertices, faces, _, _ = marching_cubes(Y==r, 0.0)

            # Undo the effect of zero padding on coordinates and transform  
            vertices -= 1
            vertices = apply_affine(vertices, V.affine)
            
            # The output from marching_cubes contains several overlapping
            # triangles and duplicate vertices, however, we rely on meshfix to
            # clean that up
            
            """
            # Detect duplicate vertices
            Nv = len(vertices)
            a, uv, iv = remove_duplicates(np.sort(vertices,1), return_index=True, return_inverse=True)
            uv = np.sort(uv)
            remverts = np.delete(np.arange(Nv), uv)
                    
            # find vertices corresponding to the duplicates and replace
            q = []
            for i in remverts:
                q.append(np.where(vertices[i].sum() == vertices.sum(1))[0][:-1])
            q=np.squeeze(np.array(q))
            
            for i,ii in zip(remverts,q):
                faces[faces==i] = ii
                
                
            # Renumber faces
            y = np.zeros(Nv,dtype=int)
            y[uv] = np.arange(len(uv))
            faces = y[faces]
            
            # Remove duplicate triangles
            _, ut = remove_duplicates(np.sort(faces,1), return_index=True)
            ut = np.sort(ut)
            faces = faces[ut]            
            
            # Remove duplicate vertices
            vertices = vertices[uv]
            
            cfu = faces
            
            ###################################
            
            ####
            
            # Remove duplicate triangles
            sortarg1   = np.argsort(faces,axis=1)
            cfs        = faces[np.arange(len(faces)).T[:,np.newaxis],sortarg1]
            cfsu,Uidx  = remove_duplicates(cfs,return_index=True)
            sortarg2   = np.argsort(sortarg1[Uidx,:],axis=1)
            cfu        = cfsu[np.arange(len(cfsu)).T[:,np.newaxis],sortarg2]
            
            # Make the mesh:
            # # triangles x vertices (of triangle) x coordinates (of vertices)
            #mesh = vertices[cfu]
            mesh = vertices[faces]
            tnormals  = get_triangle_normals(mesh)
            
            # Ensure normals are pointing out
            # (1) get center of mass of the whole mesh and for each triangle
            # (2) get the vector from total COM to each triangle COM
            # (3) count how many normals are pointing away/out (positive dot
            #     product)
            # (4) If less than half of the normals point out then they probably
            #     point inwards. Hence, reverse column order of faces.
            totalCOM = np.average(mesh, axis=(0,1))
            triCOM = np.average(mesh, axis=1)
            out       = triCOM-totalCOM[np.newaxis,:]
            count_out = np.sum(np.sum(tnormals*out,axis=1)>0)
            if count_out/float(len(mesh)) < 0.5: 
                cfu = cfu[:,::-1]
            """
            # Save surface mesh
            if len(regions)>1:
                fname = add2filename(P, "_region"+str(np.where(regions==r)[0][0]+1)+"_mesh",useext=file_format)
            else:
                fname = add2filename(P, "_mesh", useext=file_format)
            
            return_files.append(fname)
            mesh_save(vertices, faces, fname, file_format, binary)
            #mesh_save(vertices, cfu, fname, file_format, binary)
            
    return return_files

        
def get_triangle_normals(mesh):
    """Get normal vectors for each triangle in the mesh.

    PARAMETERS
    ----------
    mesh : ndarray
        Array describing the surface mesh. The dimension are:
        [# of triangles] x [vertices (of triangle)] x [coordinates (of vertices)].
    
    RETURNS
    ----------
    tnormals : ndarray
        Normal vectors of each triangle in "mesh".
    """    

    tnormals = np.cross(mesh[:,1,:]-mesh[:,0,:],mesh[:,2,:]-mesh[:,0,:]).astype(np.float)
    tnormals /= np.sqrt(np.sum(tnormals**2,1))[:,np.newaxis]
    
    return tnormals


def mesh_load(fname):
    """Load a surface mesh in either .off or .stl file format. Return the 
    vertices and faces of the mesh. If .stl file, assumes only one solid, i.e. 
    only one mesh per file.
    
    PARAMETERS
    ----------
    fname : str 
        Name of the file to be read (.off or .stl file).
    
    RETURNS
    ----------
    vertices : ndarray
        Triangle vertices.
    faces : ndarray
        Triangle faces (indices into "vertices").
    """
    file_format = os.path.splitext(fname)[1].lower()
    
    if file_format == ".off":
        with open(fname, "r") as f:
            # Read header
            hdr = f.readline().rstrip("\n").lower()
            assert hdr == "off", ".off files should start with OFF"
            while hdr.lower() == "off" or hdr[0] == "#" or hdr == "\n":
                hdr = f.readline()
            hdr = [int(i) for i in hdr.split()]
            
            # Now read the data
            vertices = np.genfromtxt(islice(f,0,hdr[0]))
            faces    = np.genfromtxt(islice(f,0,hdr[1]),
                                     usecols=(1,2,3)).astype(np.uint)
        return vertices,faces
        
    elif file_format == ".stl":
        # test if ascii. If not, assume binary        
        with open(fname, "rb") as f:
            try:
                if f.readline().decode().split()[0] == "solid":
                    is_binary = False
                else:
                    is_binary = True
            except:
                is_binary = True
                
        if is_binary:
            with open(fname, "rb") as f:
                # Skip the header (80 bytes), read number of triangles (1
                # byte). The rest is the data.
                np.fromfile(f, dtype=np.uint8, count=80)         
                np.fromfile(f, dtype=np.uint32, count=1)[0]
                data = np.fromfile(f, dtype=np.uint16, count=-1)
            data = data.reshape((-1,25))[:,:24].copy().view(np.float32)
            vertices = data[:,3:].reshape(-1,3) # discard the triangle normals
            
        else:
            vertices = []
            with open(fname,"rb") as f:
                for line in f:
                    line = line.decode().lstrip().split()
                    if line[0] == "vertex":
                        vertices.append(line[1:])
            vertices = np.array(vertices, dtype=np.float)
        
        # The stl format does not contain information about the faces, hence we
        # will need to figure this out.
        faces = np.arange(len(vertices)).reshape(-1,3)

        # Remove vertice duplicates and sort rows by sum    
        sv = np.sum(vertices+vertices*(100*np.random.random(3))[np.newaxis,:],
                    axis=1)
        sv_arg = np.argsort(sv)
        sv_arg_rev = np.argsort(sv_arg) # reverse indexing for going back
        
        # Get unique rows, indices of these, and counts. Create the new indices
        # and repeat them
        u, u_idx, u_count = np.unique(sv[sv_arg],return_index=True,
                                      return_counts=True)
        repeat_idx = np.repeat(np.arange(len(u)), u_count)
        
        # Retain only unique vertices and modify faces accordingly
        vertices = vertices[sv_arg][u_idx]
        faces = repeat_idx[sv_arg_rev][faces]
        
        return vertices, faces
            
    else:
        raise IOError("Invalid file format. Only files of type .off and .stl are supported.")

        
def mesh_save(vertices,faces,fname,file_format="off",binary=True):
    """Save a surface mesh described by points in space (vertices) and indices
    into this array (faces) to an .off or .stl file.
    
    PARAMETERS
    ----------
    vertices : ndarray
        Array of vertices in the mesh.
    faces : ndarray, int
        Array describing the faces of each triangle in the mesh.
    fname : str
        Output filename.
    file_format : str, optional
        Output file format. Choose between "off" and "stl" (default = "off").
    binary : bool
        Only used when file_format="stl". Whether to save file as binary (or
        ascii) (default = True).

    RETURNS
    ----------
    Nothing, saves the surface mesh to disk.
    """
    nFaces = len(faces)
    file_format = file_format.lower()
    
    # if file format is specified in filename, use this
    if fname.split(".")[-1] in ["stl","off"]:
        file_format = fname.split(".")[-1].lower()
    else:
        fname = fname+"."+file_format
    
    if file_format == "off":
        nVertices = len(vertices)
        with open(fname, "wb") as f:
            f.write("OFF\n".encode())
            f.write("# File created by ... \n\n".encode())
            np.savetxt(f,np.array([nVertices,nFaces,0])[np.newaxis,:],fmt="%u")
            np.savetxt(f,vertices,fmt="%0.6f")
            np.savetxt(f,np.concatenate((np.repeat(faces.shape[1],nFaces)[:,np.newaxis],faces),axis=1).astype(np.uint),fmt="%u")
    
    elif file_format == "stl":
        mesh = vertices[faces]
        tnormals  = get_triangle_normals(mesh)
        data = np.concatenate((tnormals, np.reshape(mesh, [nFaces,9])),
                              axis=1).astype(np.float32)
        
        if binary:
            with open(fname, "wb") as f:
                f.write(np.zeros(80, dtype=np.uint8))
                f.write(np.uint32(nFaces))
                f.write(np.concatenate((data.astype(np.float32,order="C",copy=False).view(np.uint16),np.zeros((data.shape[0],1),dtype=np.uint16)),axis=1).reshape(-1).tobytes())
        else:
            with open(fname, "w") as f:
                f.write("solid MESH\n")
                for t in range(len(data)):
                    f.write(" facet normal {0} {1} {2}\n  outer loop\n   vertex {3} {4} {5}\n   vertex {6} {7} {8}\n   vertex {9} {10} {11}\n  endloop\n endfacet\n"\
                    .format(*data[t,:]))
                f.write("endsolid MESH\n")
    else:
        raise IOError("Invalid file format. Please choose off or stl.")


def msh2vol(fname_mesh, vol_hdr):
    """Voxelize a surface mesh.
    
    PARAMETERS
    ----------
    fname_mesh : str
        Filename of the surfacemesh to voxelize.
    vol_hdr : nibabel image object
        Image object from which to get header information (this is used to
        define the voxel space, e.g., which points to test).
    
    RETURNS
    ----------
    Y : ndarray
        Array describing the volume mask.
    """

    # Initialize image array
    img_dims = vol_hdr.shape
    Y = np.zeros(img_dims)

    # Load the surface
    vertices,faces = mesh_load(fname_mesh)
  
    if len(vertices)>0:
        vertices = apply_affine(vertices, np.linalg.inv(vol_hdr.affine))
        mesh = vertices[faces]

        # Generate points to test
        testPoints = np.array(np.meshgrid(*tuple(map(np.arange,img_dims)),
                                      indexing="ij")).reshape((3,-1)).T
      
        # Get indices of points within mesh
        inMesh,cList = in_mesh(mesh,testPoints)
        if cList.size > 0:
                log("Some voxels could not be unambiguously determined to be inside or outside the mesh:", level=logging.WARNING)
                log(str(testPoints[cList]), level=logging.WARNING)

        # Generate and write volume (assign a value of one to points inside the
        # mesh)
        Y[testPoints[inMesh,0], testPoints[inMesh,1], testPoints[inMesh,2]] = 1

    else:     
        log("Surface file was empty! Return empty volume", level=logging.WARNING)

    return Y

    
def in_mesh(mesh, testPoints, quiet=False):
    """Determines whether points are inside a mesh by tracing rays in three 
    directions (x,y,z). The idea is that if a point is inside the mesh then it
    will cross this mesh an unenven number of times in each direction
    (+/- x,y,z).
    
    PARAMETERS
    ----------
    mesh : array_like
        Array describing the triangles making up the mesh. n-by-3-by-3:
            [# triangles] x [# vertices (per triangle)] x [x,y,z of each vertex] 
    testpoints : array_like
        Points to test (if inside the mesh). m-by-3:
            [# Points] x [x,y,z]
    quiet : bool
        Whether or not to output progress to the command window.
    
    RETURNS
    ----------
    inMesh          array of booleans stating, for each point in "testpoints",
                    whether it is inside the mesh or not
    cList           list of points for which the analysis was unable to
                    determine if they were inside the mesh or not

    ACKNOWLEDGEMENTS
    ----------------
    This function is partly inspired by the MATLAB function "intriangulation"
    by Johannes Korsawe (which itself relies heavily on another MATLAB
    function, "VOXELISE", by Adam H. Aitkenhead). 
    """

    # check inputs
    
    # Only those points for which the previous analysis did not yield clear 
    # results (i.e. voxels for which the ray appeared to pass an uneven number
    # of triangles---found in cList) are passed on to subsequent analyses
    lvl = logging.DEBUG
    
    if not quiet:
        log("Making splits in z direction, tracing in x", level=lvl)
        
    inMesh,cList = ray_trace(mesh,testPoints,quiet)

    if cList.size > 0:
        if not quiet:
            log("Making splits in x direction, tracing in y", level=lvl)
        inMesh1,cList1 = ray_trace(mesh[:,:,[1,2,0]],testPoints[cList[:,np.newaxis],[1,2,0]],quiet)
        inMesh[cList[inMesh1]] = True # points determined to inside the mesh
                                      # based on the information from the x direction
        cList = cList[cList1]         # remaining indices with unclear results
        
        if cList.size > 0:
            if not quiet:
                log("Making splits in y direction, tracing in z", level=lvl)
            inMesh2,cList2 = ray_trace(mesh[:,:,[2,0,1]],testPoints[cList[:,np.newaxis],[2,0,1]],quiet)     
            inMesh[cList[inMesh2]] = True
            
            cList = cList[cList2]
        
        if cList.size > 0:
            if not quiet:
                log("Unambiguous point(s) remaining. Rotating points and mesh to make oblique splits and trace lines. Doing at most 10 iterations", level=lvl)
            for k in np.arange(1,11):
                if not quiet:
                    log("Iteration "+str(k), level=lvl)
    
                # randomly choose rotation angles (between 0 and pi)
                a = np.pi/2*np.random.random(1)[0]
                b = np.pi/2*np.random.random(1)[0]
                c = np.pi/2*np.random.random(1)[0]
                
                # make axis specific rotation matrices
                Rx = np.array([[1,0,0],[0,np.cos(a),np.sin(a)],[0,-np.sin(a),np.cos(a)]])
                Ry = np.array([[np.cos(b),0,-np.sin(b)],[0,1,0],[np.sin(b),0,np.cos(b)]])
                Rz = np.array([[np.cos(c),np.sin(c),0],[-np.sin(c),np.cos(c),0],[0,0,1]])
            
                # combined rotation matrix (rotating around all axes)
                # (R = Rz*Ry*Rx)
                R = np.dot(np.dot(Rz,Ry),Rx)
                
                # post-multiplication: when point is represented as row vector
                # pre-multiplication: when point is represented as column vector
                # this produces opposite rotations, i.e. Rv != wR, but
                # Rv == wR^T
                # Hence, here
                # np.dot(testPoints[cList],R) # rotate test points
                # np.dot(mesh,R)              # rotate surface mesh
                # which gives a CLOCKWISE rotation
                
                inMeshR,cListR = ray_trace(np.dot(mesh,R),np.dot(testPoints[cList],R),quiet=True)
                inMesh[cList[inMeshR]] = True      
                cList = cList[cListR]
                
                if cList.size == 0:
                    break
        
    return inMesh,cList


def ray_trace(mesh, testPoints, quiet=False, nsplits="auto", 
              centerOfMassLimit=10000, num_rand_points=100):
    """Determines whether points are inside the mesh by tracing rays in a given
    direction. Makes splits in the third dimension (usually z) and ray traces
    along the first dimension (x). Thus, permuting the coordinate axes of mesh
    and testPoints allows one to change this.
    If a split is sufficiently large (i.e. has many points), calculates the
    center of mass of the split and determines if it is inside the mesh. If so,
    points within a radius of this point (the radius being calculated as the
    distance of the center of mass to the nearest triangle vertice) are
    automatically determined to be within the mesh as well. If not, it simply
    goes on to test all points in the split. This is particularly useful for
    large, filled surface meshes.
    
    PARAMETERS
    ----------    
    mesh : array_like
        Array describing the triangles making up the mesh. n-by-3-by-3:
            [# triangles] x [# vertices (per triangle)] x [x,y,z of each vertex] 
    testpoints : array_like
        Points to test (if inside the mesh). m-by-3:
            [# Points] x [x,y,z]
    quiet : bool, optional
        Whether or not to output progress to the command window.
    nsplits : int, optional
        Number of splits to process the data in. By default, data are split in
        integer values between max and min values of the mesh and test points
        (default = "auto").
    centerOfMassLimit : int | float, optional
        The minimum number of points required in a split before the split center
        of mass is calculated and used as a way of quickly determining if 
        (some) points are inside or outside the mesh (default = 10000).
    num_rand_points : int, optional
        Number of randomly selected points to test if centerOfMassLimit is
        reached (default = 100).
    
    RETURNS
    ----------
    inMesh : ndarray
        Array of booleans stating, for each point in testPoints, whether it is
        inside the mesh or not.
    correction_list : ndarray
        Indices of points which might be inside the mesh but which the current
        iteration (orientation of analysis) was unable to unambiguously
        resolve.
    """
        
    # Prepare output
    inMesh = np.zeros(len(testPoints),dtype=bool)
    correction_list = []
    
    # Min and max coordinates of the mesh
    mesh_min = mesh.min(axis=(0,1))
    mesh_max = mesh.max(axis=(0,1))
     
    # Min and max coordinates of each triangle
    mesh_tri_min = mesh.min(1)
    mesh_tri_max = mesh.max(1)
    
    # Make initial bounding box; indices of points within bounding box
    bbox = np.all((testPoints >= mesh_min) & (testPoints <= mesh_max), axis=1)

    splits, split_size = make_splits(mesh, testPoints, nsplits)
    
    for s in splits:
        if not quiet:
            progress_bar(s, splits[-1], 50)
            
        # SLAB SETUP
        # splits are evenly spaced between upper and lower boundaries of
        # test points/mesh. Each split is the center +/- 0.5*split_size to
        # cover the entire space. Points within and triangles crossing this are
        # considered in a particular split.
        # ---------------------------------------------------------------------
        
        split_zmin = s-split_size/2.
        split_zmax = s+split_size/2.
  
        # get triangles in current split (i.e. all but the ones which only have
        # vertices above/below s +/- 0.5*split_size)
        stris_idx = (mesh_tri_max[:,2] >= split_zmin) & (mesh_tri_min[:,2] <= split_zmax)
        
        if stris_idx.sum() == 0:
            continue
        
        stris = mesh[stris_idx] # triangles in this split
        
        smin = stris.min(axis=(0,1))
        smax = stris.max(axis=(0,1))
        smin[2] = split_zmin
        smax[2] = split_zmax
        
        # indices of test points in the current split
        spoints_idx = np.where(bbox & np.all((testPoints >= smin) & (testPoints <= smax),
                                           axis=1))[0]

        if len(spoints_idx) == 0: # if no points to test
            continue
        
        # CHECK VOXELS
        # -----------------------------
        
        # (1) Approximate
        
        # Fast way of quickly determining a (possibly) large amount of points
        # are inside or outside. This makes circles around a few randomly
        # selected points and assigns all points within these circles at once.
        # Thus, this is most effective when the surface is not very convoluted.
        if len(spoints_idx) > centerOfMassLimit:
            # (1) Choose a number of (random) points
            center_of_mass = np.round(stris.mean((0,1))[None,:]).astype(np.int)            
            randpoints  = spoints_idx[np.random.randint(0,len(spoints_idx), num_rand_points)]
            points2test = np.concatenate((center_of_mass, testPoints[randpoints]))
            
            # (2) For each point, check if it is inside or outside the mesh
            i, c = in_mesh(stris, points2test, quiet=True)
            
            # (3) Get distance from each point to closest vertex in stris
            radii = np.min(np.sqrt(np.sum((stris[np.newaxis,...]-points2test[:,np.newaxis,np.newaxis,:])**2,axis=3)),axis=(1,2))
            
            # (4) Assign all points within a circle the same label as the
            #     center point and remove from the list to check
            inCircle = np.sqrt(np.sum((testPoints[spoints_idx][np.newaxis,...]-points2test[:,np.newaxis,:])**2,axis=2)) < np.floor(radii)[:,np.newaxis]
            inMesh[spoints_idx[np.any(inCircle[i], axis=0)]] = True
            
            spoints_idx = spoints_idx[~np.any(inCircle,axis=0)]
        
        # (2) Precise
        
        # For the (remaining) test points, make an array telling which
        # triangles which could possibly be crossed for each point. Since this
        # number may be different for different points, limit the number of
        # array columns to that of the maximum number of crossings (pad with
        # zeros)
        spoints = testPoints[spoints_idx]

        # Start from one since zero means no possibility of crossing; fix later
        stris_idx2 = np.arange(1,stris_idx.sum()+1) 

        possible_cross = stris_idx2[None,:]*((mesh_tri_max[stris_idx,1] >= spoints[:,1][:,None]) & (mesh_tri_min[stris_idx,1] <= spoints[:,1][:,None]))
        nc = (possible_cross>0).sum(1) # number of triangles which may be crossed
        mc = nc.max()
        
        if mc == 0:
            continue
        
        possible_cross = np.reshape(np.insert(possible_cross[possible_cross>0],np.cumsum(nc).repeat(mc-nc),[0]*np.sum(mc-nc)),[-1,mc])
        ok = possible_cross != 0
        
        possible_cross[possible_cross>0]-=1 # fix indexing error
        
        # Test for intersections
        ray_direction = [1, 0, 0] # trace along 1st axis
        intersect, intersect_dist = ray_triangle_intersect(stris,
                                    spoints, ray_direction,                                                           
                                    posint=possible_cross,
                                    posint_ok=ok, return_int_dists=True)
        
        intersect_dist = row_wise_unique(intersect_dist, fill_value=np.inf)
        intersect_dist = np.sort(intersect_dist, 1) # to be able to reduce array size

        nc = intersect.sum(1) # number of triangles crossed by each point
        mc = nc.max()
        if mc == 0:
            continue
        
        ncol = mc+(mc%2) # ensure equal!
        intersect_dist = intersect_dist[:,:ncol] # reduce array size
        
        # Points which cross the surface an equal number of times and have a
        # triangle pair on each side are INSIDE the surface
        inMesh[spoints_idx[(nc%2 == 0) & (np.sign(intersect_dist[:,0::2]) != np.sign(intersect_dist[:,1::2])).any(1)]] = True

        # Points which cross the surface an uneven number of times cannot be
        # determined unambiguously so save their index and test in another
        # direction
        correction_list.extend(spoints_idx[nc % 2 > 0])
       
    return inMesh, np.array(correction_list, dtype=np.uint)


def ray_triangle_intersect(triangles, ray_origin, ray_direction, posint=None,
                           posint_ok=None, plane_type="two-sided",
                           mindist2int=-np.inf, maxdist2int=np.inf,
                           return_int_points=False, return_int_dists=False, eps=1e-6):
    """Test intersection between rays and triangles.
    
    PARAMETERS
    ----------
    triangles : ndarray
        Nx3x3 where N is # of triangles each of which are described by a 3-by-3
        (vertices-by-coordinates) array.
    ray_origin : ndarray
        Points of ray origin.
    ray_direction : ndarray
        Vector describing the direction of the ray originating from ray_origin.
        This may be specified (1) one vecetor per point in ray_origin, or (2)
        a single vector in which case the ray direction is assumed to be the
        same for all origins.
    posint : ndarray, optional
        An array describing, for each point in ray_origin, which triangles the
        ray could possibly intersect. If testing a large number of points
        and/or triangles it is highly recommended to do an initial (coarse)
        check for possible point-triangle intersection pairs. Otherwise, all
        points will need to be checked against all triangles, which (1) is
        likely to be very redundant, and (2) may even result in a MemoryError
        being raised.
        Specified either as a numpy.ndarray or list of lists (in the latter
        case, posint_ok need not be specified).
    posint_ok ndarray, optional
        An array of shape posint describing which entries in posint are real
        (actual) possible intersection pairs. This array is needed since it is
        highly unlikely that all ray origin points will have the exact same
        number of possible intersections, hence the need to somehow pad the
        array. posint_ok specified which values are 'real' and which are the
        result of padding. Note, however, that this is not the case if posint
        is specified as a list of lists containing, for each ray origin, all
        possible triangles it might intersect. Lists, however, does not support
        vectorized operations, and so will need to be converted to a
        numpy.ndarray in any case.
    plane_type : {"one-sided", "two-sided"}, optional
        The point of intersection between a ray and the plane spanned by the
        triangle sides may be computed as some multiple of the ray_direction
        (from ray_origin). If one-sided, accept only intersections in the
        actual direction of the ray (positive multiple). If two-sided, accept
        all intersections no matter of there are in the actual direction of the
        ray or in the opposite direction (default = "two-sided").
    mindist2int : float
        Minimum distance for which to consider intersections
        (default = -np.inf).
    maxdist2int : float
        Maximum distance for which to consider intersections
        (default = np.inf).
    return_int_points : bool, optional
        Return the exact points where each ray crosses the plane of each
        triangle (default = False).
    return_int_dists : bool, optional
        Return the distance from ray_origin to the point of intersection with
        the plane spanned by each triangle (default = False).
    eps : float, optional
        Error tolerance (default = 1e-6).
    
    RETURNS
    ----------
    intersect : ndarray
        Array which triangles are intersected by the ray from each ray_origin.
    intersect_points : ndarray, optional
        The points where each ray crosses the plane of a particular triangle.
    intersect_distances : ndarray, optional
        The distance from ray_origin to the point of interscetion with the
        plane spanned by each triangle.
        
    NOTES
    ----------
    
    """
    # Check inputs
    ray_direction = np.array(ray_direction)
    if ray_direction.shape == ray_origin.shape:
        ray_direction = ray_direction[:,None,:]
    elif ray_direction.size == 3: # broadcast to all ray origins     
        ray_direction = ray_direction[None,None,:]
    else:
        raise ValueError("ray_direction must be either a single vector or have same dimensions as ray_origin.")
    
    try:
        assert plane_type in ["one-sided", "two-sided"]
    except AssertionError:
        raise ValueError("Choose 'one-sided' or 'two-sided'.")
    if posint is None:
        posint = []
    if posint_ok is None:
        posint_ok = []
        
    if (len(posint) > 0) & isinstance(posint,list):
        posint, posint_ok = list2numpy(posint, dtype=np.int)
        
    # Get vectors (e1, e2) spanning each triangle from point v0
    v0 = triangles[:,0,:]
    e1 = triangles[:,1,:]-v0
    e2 = triangles[:,2,:]-v0
    if len(posint)>0:
        try:
            assert posint.shape == posint_ok.shape
        except AssertionError:
            raise ValueError("posint and posint_ok must have same dimensions!")
    
        v0=v0[posint]
        e1=e1[posint]
        e2=e2[posint]
    else:
        v0=v0[None,...]
        e1=e1[None,...]
        e2=e2[None,...]
    
    # if not specified, consider all intersections possible
    if len(posint_ok) == 0:
        posint_ok = True
        
    # RAY TESTING
    # =======================
    # Implementation of the algorithm presented in
    # 
    # Moeller & Trumbore (1997). Fast, Minimum Storage Ray/Triangle
    # Intersection. Journal of Graphics Tools, 2(1):21--28.
    
    # Vector perpendicular to e2 and ray_direction
    #   P = CROSS(dir, e2)
    P = np.cross(ray_direction, e2)

    # Determinant
    #   D = DOT(P, e1)
    # D describes the relationship between triangle normal (face of the
    # triangle) and ray direction: if det>0 then the ray points towards the
    # outer face of the triangle (i.e. the ray and the triangle normal point
    # towards each other); if det<0, they point in the same direction; if
    # det=0, the ray and triangle are parallel.
    det = np.einsum("ijk,ijk->ij", P, e1)    
    inv_det = 1./(det+1e-10)
    
    # Vector from v0 to ray_origin
    #   T = O-V0
    tvec = ray_origin[:,None,:] - v0

    # 1st barycentric (e.g., x) coordinate of intersection point, i.e. where 
    # the ray intersects the plane spanned by e1 and e2.
    #   u = DOT(T, P)
    u = np.einsum("ijk,ijk->ij", P, tvec)*inv_det
    
    Q = np.cross(tvec, e1)
    
    # 2nd barycentric (e.g., y) coordinate of intersection point
    #   v = DOT(dir, Q)
    v = np.einsum("ijk,ijk->ij", ray_direction, Q)*inv_det
    
    # Distance from ray_origin to point of intersection in units of
    # ray_direction
    t = np.einsum("ijk,ijk->ij", e2, Q)*inv_det

    # CHECK CONDITIONS
    # =======================
    
    # Check that the ray crosses the plane spanned by vectors E1 and E2 in the
    # direction of the ray (one-sided) or if it crosses on either side ()
    if plane_type == "one-sided":
        posint_ok = posint_ok & (det >= eps)
    elif plane_type == "two-sided":
        posint_ok = posint_ok & (np.abs(det) >= eps)
    
    # Check that the ray actually crosses a triangle by testing if the
    # barycentric coordinates fulfills the equation
    #   T(u,v) = (1-u-v)*v0+u*v1+v*v2
    # which they do when u>=0, v>=0, (u+v)<=1.
    # Additionally, check if the intersection is within the allowed distance
    # (maxdist2int) by testing t.
    posint_ok = posint_ok & (u>=-eps) & (u<=1+eps)
    posint_ok = posint_ok & (v>=-eps) & ((u+v)<=1+eps)
    intersect = posint_ok & (t>(mindist2int+eps)) & (t<(maxdist2int+eps))

    t[~posint_ok]=np.inf
    if return_int_points:
        # The points where each ray cross the plane of each triangle
        intersect_points = ray_origin[:,np.newaxis,:]+t[...,np.newaxis]*ray_direction        
        if return_int_dists:      
            return intersect, intersect_points, t
        else:
            return intersect, intersect_points
    elif return_int_dists:
        return intersect, t
    else:
        return intersect


def make_splits(mesh, testpoints, nsplits="auto"):
    """Make splits of the data.
    
    """
    # split data along 3rd dimension and process individually
    split_min = np.floor(np.max([mesh[:,2].min(), testpoints[:,2].min()]))
    split_max = np.ceil( np.min([mesh[:,2].max(), testpoints[:,2].max()]))
   
   
    if nsplits == "auto":
        nsplits = split_max-split_min+1    
    else:
        assert isinstance(nsplits,int)
        assert nsplits > 0

    splits = np.linspace(split_min, split_max, nsplits)

    if len(splits) > 1:
        split_size = splits[1]-splits[0]
    else:
        split_size = 2*(split_max-split_min)
        
    return splits, split_size
    
       
def make_surface_mesh(fname_vol, fname_mesh, n_shells=1, min_shell_area=1,
                      vertex_density=0.5, n_smooth=1, erode_and_expand=False,
                      se_connectivity=2, mm2move=1.5, ensure_distance=1,
                      nsteps=6, remove_spikes_curv=5.0):
    """Create maximum n_shells largest surface meshes from a volume mask.
    
    PARAMETERS
    ----------
    fname_vol : str
        Filename of volume mask (nifti file).
    fname_mesh :str
        Filename of mesh to create.
    n_shells : int
        Maximum number of shells to retain.
    min_shell_area : float
        Minimum area of the shells to keep.
    vertex_density : float
        Density of vertices in the surface mesh (per mm2).
    n_smooth : int
        Number of times to apply Taubin smoothing.
    erode_and_expand : bool, optional
        Whether or not to erode the volume mask before appying marching cubes
        and then grow the resulting surface mesh (default = False).
    se_connectivity : {1,2,3}, optional
        The connectivity of the structuring element (se) used for erosion. The
        structuring element has rank 3 (default = 2).
    mm2move : float, optional
        The amount by which to move a vertice each iteration (default = 1.5).
    ensure_distance : float, optional
        Minimum distance to ensure in the mesh (default = 1).
    nsteps : int, optional
        The number of steps used to reach the desired distance to move
        (default = 6).
    remove_spikes_curv: float, optional.
        Curvature parameter for spike removal in meshfix. Default: 5.0
        
    RETURNS
    ----------
    Writes fname_mesh to disk with the number 1 appended (if multiple shells
    are retained, numbering continues as 2,3,...n_shells).
    
    NOTES
    ----------
    Uses meshfix.
    """
    assert se_connectivity in [1,2,3]
    assert isinstance(nsteps,int)
    
    lvl = logging.DEBUG
    
    mesh_format = "off"
    if fname_mesh.endswith(".off"):
        mesh_format = "off"
        fname_mesh = fname_mesh[:-(len(mesh_format)+1)]                            
    meshfix = path2bin("meshfix")
    
    # meshfix commands
    split_shells     = '"{0}" "{1}" -a 2.0 --no-clean -q --shells {2} --splitShells --shellMinArea {3} -o "{4}"'.format(meshfix,"{0}",n_shells,min_shell_area,fname_mesh)
    smooth_surface   = '"{0}" "{1}" -a 2.0 --taubin {2} --no-clean -o "{1}"'.format(meshfix,"{0}", n_smooth)
    remove_spikes    = '"{0}" "{1}" -a 2.0 --removeSpikes {2} -o {1}'.format(meshfix,"{0}", remove_spikes_curv)
    resample_surface = '"{0}" "{1}" -a 2.0 -u 5 -q --vertexDensity {2} -o "{1}"'.format(meshfix,"{0}",vertex_density)
    clean_surface    = '"{0}" "{1}" -a 2.0 -q -o "{1}"'.format(meshfix,"{0}")
    uniformremesh_surface  = '"{0}" "{1}" -a 2.0 -u 1 -q -o "{1}"'.format(meshfix,"{0}")
    save_stl = '"{0}" "{1}" --no-clean --stl -o "{2}"'.format(meshfix,"{0}", "{1}")
    #taubin_smooth_surface   = "{0} {1} -a 2.0 --smooth {2} --no-clean -o {1}".format(meshfix,"{0}", n_smooth)

    log("Generating surface mesh from {0}",os.path.basename(fname_vol), level=lvl)
    
    if erode_and_expand:
        log("Eroding", level=lvl)
        se = mrph.generate_binary_structure(3,se_connectivity)
        vol = nib.load(fname_vol)
        vole = mrph.binary_erosion(vol.get_data(), se, 1)
        fname_vol = add2filename(fname_vol, "_ERO")
        write_nifti(vole, fname_vol, vol, dtype=np.uint8)
        
        # since the erosion is specified in terms of voxel size, scale mm2move
        # by voxel size as well
        mm2move *= np.average(get_vx_size(vol))
    
    log("Creating initial mesh", level=lvl)
    input_file = vol2msh(fname_vol, mesh_format)
    
    # clean up the raw surface and split into shells, resampling to vertex
    # density and keeping only the nshells largest shells
    log("Keeping maximum {0} largest shell(s) with minimum area of {1}",(n_shells, min_shell_area), level=lvl)
    for f in input_file:
        spawn_process(split_shells.format(f))
        os.remove(f)
    
    # process each shell individually (smooth and clean twice)
    files = sorted(glob(fname_mesh+"*."+mesh_format), key=str.lower)
    for f in files:
        log("Processing {0}",os.path.basename(f), level=lvl)

        if erode_and_expand:
            log("Cleaning", level=lvl)
            spawn_process(clean_surface.format(f))

            log("Expanding", level=lvl)
            vertices,faces = mesh_load(f)
            vertices = deform_mesh(vertices, faces, mm2move, ensure_distance,
                                   nsteps, deform="expand", smooth_mesh=False,
                                   verbose=True)
            mesh_save(vertices, faces, f)

            log("Cleaning", level=lvl)
            spawn_process(clean_surface.format(f))

        spawn_process(save_stl.format(f, f[:-4]+'.stl'))
        log("Resampling to a vertex density of {0}", vertex_density, level=lvl)
        spawn_process(resample_surface.format(f))
        spawn_process(save_stl.format(f, f[:-4]+'_resampled.stl'))
        log("Cleaning", level=lvl)
        spawn_process(clean_surface.format(f))
        spawn_process(save_stl.format(f, f[:-4]+'_resampled_clean.stl'))
        log("Smoothing {0} time(s)", n_smooth, level=lvl)
        spawn_process(smooth_surface.format(f))
        spawn_process(save_stl.format(f, f[:-4]+'_resampled_clean_smoothed.stl'))
        log("Removing Spikes with curvature parameter {0}", remove_spikes_curv, level=lvl)
        spawn_process(remove_spikes.format(f))
        spawn_process(save_stl.format(f, f[:-4]+'_resampled_clean_smoothed_spike.stl'))
        log("Uniform remeshing", level=lvl)
        spawn_process(uniformremesh_surface.format(f))
        spawn_process(save_stl.format(f, f[:-4]+'_resampled_clean_smoothed_spike_remesh.stl'))
        log("Cleaning", level=lvl)        
        spawn_process(clean_surface.format(f))

    return files


def make_volume_mask(fname_mesh, fname_vol, vol_hdr, write_volume=True):    
    """Create a volume mask (nifti file) from a surface mesh and (optionally)
    save it to disk.
    
    PARAMETERS
    ----------
    fname_mesh : str
        Filename of the surface file (including file type).
    fname_vol : str
        Filename of the volume mask (to be created).
    vol_hdr : nibabel image object
        Volume object from which the header information is taken.
    write_volume : bool, optional
        Whether or not to save the volume mask to disk (default = True)
    
    RETURNS
    ----------
    vol : nibabel image object
        The volume which was created.
    """    
    lvl = logging.DEBUG
    
    log("Creating {0} from {1}",
        [os.path.basename(fname_vol), os.path.basename(fname_mesh)], level=lvl)
        
    Y = msh2vol(fname_mesh, vol_hdr)
    
    vol = write_nifti(Y, fname_vol, vol_hdr, dtype=np.uint8,
                      write_volume=write_volume, return_volume=True)
    
    return vol


def decouple_volumes(v1, v2, mode, se=None, iterations=1):
    """
    
    mode : {inner-from-outer, outer-from-inner, neighbors}
        inner-from-outer: this changes v1 by removing voxels
        outer-from-inner: this changes v2 by adding voxels
        neighbors: this changes v2 by removing voxels
    
    """
    assert mode in ["inner-from-outer","outer-from-inner","neighbors"]
    
    if isinstance(v1, str) and os.path.isfile(v1):
        v1 = nib.load(v1)
    assert isinstance(v1, nib.Nifti1Image) or isinstance(v1, nib.Nifti2Image)
    d1 = v1.get_data()
    if isinstance(v2, str) and os.path.isfile(v2):
        v2 = nib.load(v2)
    assert isinstance(v2, nib.Nifti1Image) or isinstance(v2, nib.Nifti2Image)
    d2 = v2.get_data()
    
    assert d1.ndim is d2.ndim
    
    
    if se is None:
        se = mrph.generate_binary_structure(d1.ndim,1)
    
    if mode == "inner-from-outer":
        # make v2/d2 the inner volume
        d1, d2 = d2, d1
        v1, v2 = v2, v1        
        d2 = d2 & mrph.binary_erosion(d1, se, iterations)
        
    if mode == "outer-from-inner":
        d2 = d2 | mrph.binary_dilation(d1, se, iterations)
        
    if mode == "neighbors":
        d2 = d2 & ~mrph.binary_dilation(d1, se, iterations)
    
    d2 = nib.Nifti1Image(d2, v2.affine, header=v2.header)
    d2.set_filename(v2.get_filename())
    return d2


def decouple_surfaces(surf1, surf2, mode, cut_inner=False, min_distance=0,
                      verbosity=logging.DEBUG):
    """Decouple two surfaces from each other.
    
    PARAMETERS
    ----------
    surf1 : str
        Filename of surface 1. The INNERMOST surface.
    surf2 : str
        Filename of surface 2. The OUTERMOST surface.
    mode : {"neighbor","inner-from-outer","outer-from-inner"}
        How to perform decoupling.
        
        neighbor:
            Decouples surface 1 from a neighboring one (surface 2) by pushing 
            overlapping parts of the first surface inside and---if this does
            not succeed---cutting away part of the first surface and filling
            any gaps.
            This changes surface 1.
                
        inner-from-outer:
            Decouples the inner surface (surface 1) from the outer surface
            (surface 2) by pushing overlapping parts of the inner surface
            inside and---if this does not succeed---cutting away part of
            the inner surface and filling any gaps.
            This changes surface 1.
            
        outer-from-inner:
            Decouples the outer surface (surface 2) from the inner surface
            (surface 1) by pushing overlapping parts of the outer surface
            outside and---if this does not succeed---cutting away part of the
            outer surface and fill the gaps.
            This changes surface 2 and optionally (if cut_inner = True) also
            surface 1.  
    cut_inner : bool, optional
        If True, then the parts of the inner surface (surface 1) that still
        overlaps with the outer surface (surface 2) after a few iterations are
        cut away; this helps getting rid of "spikes" in the inner surface
        facilitating the decoupling procedure, however, this comes at the price
        of changing the inner surface and so one must subsequently ensure that
        the inner surface does not overlap with any surfaces supposed to be
        contained within this surface.
        This option is only used in combination with outer-from-inner since
        surface 1 is always modified in the other modes (default = False).
    """
    assert mode in ["neighbor","inner-from-outer","outer-from-inner"], "Illegal mode '{:s}'".format(mode)
    
    meshfix = path2bin("meshfix")
    
    if mode == "neighbor":
        surf1_name = "surface 1"
        surf2_name = "surface 2"
        
        which_decouple = "outin"
        which_cut = "inner"
        
        deform = "shrink"
        
    elif mode == "inner-from-outer":
        surf1_name = "inner surface"
        surf2_name = "outer surface"
        
        which_decouple = "inin"
        which_cut = "outer"
        
        deform = "shrink"
        
    elif mode == "outer-from-inner":
        # swap surfaces
        surf1,surf2 = surf2,surf1
        
        surf1_name = "outer surface"
        surf2_name = "inner surface" 
        
        which_decouple = "outout"   
        which_cut = "inner"
        
        deform = "expand"
    else:
        raise IOError("Mode must be either 'neighbor', 'inner-from-outer', or 'outer-from-inner'.")

    surf1_file = os.path.basename(surf1.upper())
    surf2_file = os.path.basename(surf2.upper())
    
    # meshfix commands
    check_intersect = '"{0}" "{1}" "{2}" --shells 2 --no-clean --intersect'.format(meshfix,surf1,surf2)
    decouple        = '"{0}" "{1}" "{2}" -a 2.0 --shells 2 --decouple-{3} 0 -o "{1}"'.format(meshfix,surf1,surf2,which_decouple)
    cut_surf1       = '"{0}" "{1}" "{2}" -a 2.0 --shells 2 --cut-{3} 0 -o "{1}"'.format(meshfix,surf1,surf2,which_cut)
    clean1_surf1    = '"{0}" "{1}" -a 2.0 -u 1 -q -o "{1}"'.format(meshfix,surf1)
    clean2_surf1    = '"{0}" "{1}" -a 2.0 -q -o "{1}"'.format(meshfix,surf1)
    addstl = " --stl"
    if surf1.endswith("stl"):
        decouple += addstl
        cut_surf1 += addstl
        clean1_surf1 += addstl
        clean2_surf1 += addstl
        
    if mode == "outer-from-inner" and cut_inner:
        cut_surf2    = '"{0}" "{1}" "{2}" -a 2.0 --shells 2 --cut-outer 0 -o "{1}"'.format(meshfix,surf2,surf1)
        clean1_surf2 = '"{0}" "{1}" -a 2.0 -u 1 -q -o "{1}"'.format(meshfix,surf2)
        clean2_surf2 = '"{0}" "{1}" -a 2.0 -q -o "{1}"'.format(meshfix,surf2)
        if surf2.endswith("stl"):
            cut_surf2 += addstl
            clean1_surf2 += addstl
            clean2_surf2 += addstl
        
    # while there are intersections...
    # (NB. as per the convention in bash, here returncode = 0 means that there 
    # are intersections whereas returncode = 1 means that there are no
    # intersections!)    
    log("Decoupling {0} ({1}) from {2} ({3})",
        [surf1_file, surf1_name, surf2_file, surf2_name], level=verbosity)
    dist_ensured = False
    k=0
    while not dist_ensured:
        while spawn_process(check_intersect, True) == 0:
            k+=1
            log("Iteration {0}", k, level=verbosity)
            log("Decoupling", level=verbosity)
            spawn_process(decouple)
            
            for i in range(2):
                if spawn_process(check_intersect, True) == 0:
                    log("Decoupling", level=verbosity)
                    spawn_process(decouple)
                else:
                    break
            
            # if this is not working, try cutting the outer...
                
            # which surface is cut depends on the mode:
            # for the first/inner/outer surface (neighbor/inner-from-outer/
            # outer-from-inner), cut spikes overlapping with the second/outer/inner
            # surface (i.e. pointing towards surface 2/outwards/inwards) and clean
            # twice
            log("Cutting {0} ({1})", [surf1_file, surf1_name], level=verbosity)
            spawn_process(cut_surf1)
            
            log("Cleaning {0} ({1})", [surf1_file, surf1_name], level=verbosity)
            spawn_process(clean1_surf1)
            spawn_process(clean2_surf1)
            
            # if cut_inner=True and there are still intersections:
            # cut parts of the inner surface overlapping with the outer surface
            # (i.e. pointing outwards). This may alter the inner surface
            if mode == "outer-from-inner" and cut_inner:
                if spawn_process(check_intersect,True) == 0:
                    log("Cutting {0} ({1})", [surf2_file, surf2_name], level=verbosity)
                    spawn_process(cut_surf2)              
                    
                    log("Cleaning {0} ({1})", [surf2_file, surf2_name], level=verbosity) 
                    spawn_process(clean1_surf2)
                    spawn_process(clean2_surf2)
        
        # if no intersections, try to ensure minimum distance between meshes
        # if the surfaces are very complicated this may produce unexpected
        # results causing intersections, hence try this a maximum of three
        # times
        if (min_distance > 0) & (k <= 3) & (spawn_process(check_intersect, True) == 1):
            log("Attempting to ensure a minimum distance of {} between surfaces by {:s}ing {:s} ({:s})",(min_distance,deform,surf1_file,surf1_name), level=verbosity)
            v1,f1 = mesh_load(surf1)
            v2,f2 = mesh_load(surf2)
            # no smoothing since it might undo the distance ensured
            v1_new = ensure_distance(v1,f1,v2,f2, min_distance, deform)
            mesh_save(v1_new,f1,surf1)
            log("Cleaning {0} ({1})", [surf1_file, surf1_name], level=verbosity)
            spawn_process(clean1_surf1)
            spawn_process(clean2_surf1)
            # if there are no intersections
            if spawn_process(check_intersect, True) == 1:
                dist_ensured = True 
        else:
            dist_ensured = True

    log("Total iterations needed for decoupling: {0}\n", k, level=verbosity)


def deform_mesh(vertices, faces, mm2move, ensure_distance, nsteps,
                deform="expand", plane_type="one-sided", smooth_mesh=False,
                save=None, file_format="off", eps=1e-6, verbose=False):
    """Deform a mesh by either expanding it in the direction of node normals
    or shinking it in the opposite direction of the normals.
    
    PARAMETERS
    ----------
    vertices : ndarray
        Vertices describing the mesh.
    faces : ndarray
        Faces describing the mesh.    
    mm2move : float
        The amount by which to move a vertice each iteration.
    ensure_distance : float
        Minimum distance to ensure in the mesh.
    nsteps : int
        The number of steps used to reach the desired distance to move.
    deform : {"expand", "shrink"}
        Whether to expand or shink the mesh. Expand corresponds to pushing
        the nodes outwards (in the direction of the normals) whereas shrink
        will pull the nodes inwards (in the direction of minus the normals)
        (default = "expand").
    plane_type : {"one-sided", "two-sided"}, optional
        Which type of intersections to consider. If "one-sided", only
        intersections with the face of a triangle are considered (i.e. the
        direction of movement and triangle normal are opposite). If
        "two-sided", all intersections are considered (default = "one-sided").
    smooth_mesh : bool, optional
        Local smoothing of the mesh by averaging, for each vertex which has
        been moved, the coordinates of itself and all other vertices to which
        it is connected weighting the vertex itself higher than the surrounding
        vertices (default = True).
    save : {"all", "final", None}, optional
        Which steps of the process to save (default = None).
    file_format : {"off", "stl"}, optional
        Which file format to save the mesh to (default = "off").    
    eps : float, optional
        Epsilon. Error tolerance margin. For numerical stability of results
        (default = 1e-6).
    verbose : bool
        Whether or not to print to console (default = False).
    
    RETURNS
    ----------
    vertices : ndarray
        Vertices describing the expanded mesh.
    """
    
    verbosity = logging.DEBUG
    
    # check inputs
    assert deform in ["expand", "shrink"]
    assert isinstance(nsteps,int)
    assert plane_type in ["one-sided", "two-sided"]
    assert save in ["all", "final", None]
    assert file_format in ["off", "stl"]
    
    mm2move = mm2move/float(nsteps)
    n_nodes = len(vertices)
    if isinstance(mm2move, (list, tuple, np.ndarray)):
        assert len(mm2move) == len(vertices), "The length of mm2move must match that of vertices"
        ballsize = max(mm2move)
        isarray = True
        verts2consider = np.where(mm2move!=0)[0] # no need to consider others
    else:
        assert isinstance(mm2move,(int,float))
        ballsize = mm2move
        isarray = False
        verts2consider = np.arange(n_nodes)
        
    vertices = vertices.copy() # prevent modification of input "vertices"!
    v2f = verts2faces(vertices,faces)
      
    for i in range(nsteps):
        #if verbose:
        #    log("Step {0}",i+1)
        sys.stdout.flush()
        
        mesh = vertices[faces]
        
        # get node normals by averaging normals of the triangles to which the
        # node belongs
        node_normals = get_node_normals(vertices, faces, verts2consider)
        if deform == "shrink":
            node_normals *= -1
        
        # testing all nodes against all triangles is slow if the mesh is large,
        # thus, find triangles within some distance of each node and test only
        # these for intersections. This reduces # of computations dramatically.    
        barycenters = np.average(mesh, axis=1)
        bar_tree = cKDTree(barycenters)
        ver_tree = cKDTree(vertices[verts2consider])
        res = ver_tree.query_ball_tree(bar_tree, r=(2+eps)*ballsize+2*ensure_distance,
                                       p=2)
        # possible intersections
        pi, piok = list2numpy(res, dtype=np.int)

        intersect = ray_triangle_intersect(mesh, vertices[verts2consider],
                        node_normals, pi, piok, plane_type, mindist2int=0,
                        maxdist2int=ballsize+ensure_distance)
        move = ~intersect.any(1)
        verts2consider = verts2consider[move]
        
        if verbose:
            pad = str(len(str(n_nodes)))
            msg = "Moved {:"+pad+"d} of {:d} vertices"
            log(msg, (np.sum(move), n_nodes), level=verbosity)
            sys.stdout.flush()

        # update the vertices
        if isarray:
            vertices[verts2consider] += node_normals[move]*mm2move[verts2consider,None]
        else:
            vertices[verts2consider] += node_normals[move]*mm2move
        
        if smooth_mesh:
            # smooth verts2consider and neighbors
            if len(verts2consider) < len(vertices):
                f2c = [v2f[n] for n in verts2consider]
                f2c,f2cok = list2numpy(f2c,dtype=np.int)
                f2c = f2c[f2cok] # faces of verts2consider
                nn, nnok = get_element_neighbors(vertices[faces]) # get neighboring elements
                neighbors = nn[f2c][nnok[f2c]]
                # get neighboring vertices and verts2consider themselves
                v2smooth = np.unique(np.append(faces[neighbors],verts2consider)).astype(np.int)
            else:
                v2smooth = verts2consider
                
            smoo = np.zeros_like(vertices[v2smooth])
            for i,n in zip(range(len(v2smooth)),v2smooth):
                smoo[i] = np.average(vertices[faces[v2f[n]]], axis=(0,1))
            vertices[v2smooth] = smoo

        if save == "all" or (save == "final" and (i+1)==nsteps):
            pad = str(len(str(nsteps))-1)
            filename = "mesh_expand_{:0"+pad+"d}_of_{:d}."+file_format
            filename = filename.format(i+1, nsteps)
            mesh_save(vertices, faces, filename)
    
    return vertices


def list2numpy(L, pad_val=0, dtype=np.float):
    """Convert a python list of lists (the sublists being of varying length)
    to a numpy array.
    
    PARAMETERS
    ----------
    L : list
        The list of lists.
    pad_val : float, int
        The value with which to pad numpy array.
    dtype : datatype, optional
        Datatype of the output array.
        
    RETURNS
    ----------
    narr : ndarray
        L expressed as a numpy array.
    """    

    max_neighbors = len(sorted(L, key=len, reverse=True)[0])
    narr = np.array([r+[np.nan]*(max_neighbors-len(r)) for r in L])    
    ok = ~np.isnan(narr)
    narr[~ok] = pad_val    
    narr = narr.astype(dtype)
    
    return narr, ok


def verts2faces(vertices, faces, pad_val=0, array_out_type="list"):
    """Generate a mapping from vertices to faces in a mesh, i.e. for each
    vertices, which elements are it a part of.
    
    PARAMETERS
    ----------
    vertices : ndarray
        Vertices describing the mesh.
    faces : ndarray
        Faces describing the mesh.
    array_out_type : {"list", "numpy_array"}, optional
        Output type. Numpy arrays will enable vectorized operations to be
        performed on the output, however, in the case of variable number of
        elements per vertice, this will have to be padded

    RETURNS
    ----------
    v2f : {list, ndarray}
        The mapping from vertices to faces.    
    ok : ndarray
        Array describing which entries in v2f are actual faces and which are 
        "artificial". Since in a mesh, different vertices will often be part of
        different numbers of elements, some rows will have to be padded. This
        array is only returned if array_out_type is set to "numpy_array" since
        sublists of a list can be of variable length.
    """
    # Mapping from node to triangles, i.e. which nodes belongs to which
    # triangles
    v2f = [[] for i in range(len(vertices))]
    for t in range(len(faces)):  
        for n in faces[t]:
            v2f[n].append(t)
    
    if array_out_type == "list":
        return v2f        
    elif array_out_type == "numpy_array":        
        v2f, ok = list2numpy(v2f, pad_val, np.int)        
        return v2f, ok
    else:
        raise ValueError("Array output type must be list or numpy array.")    


def get_node_normals(vertices, faces, mask=None, smooth=False,
                     weigh_by_area=False):
    """Compute the normal of vertices in a mesh.
    
    PARAMETERS
    ----------
    vertices : ndarray
        Vertices describing the mesh.
    faces : ndarray
        Faces describing the mesh.
    mask : array_like
        Array of indices or boolean array describing which vertices to return
        the normals of. Returns normals for all vertices by default.
    smooth : bool, optional
        Whether or not to smooth vertice normals. A vertice normal is smoothed
        by taking into consideration also the normals of the neighboring
        vertices when computing the normal (default = False).
    weigh_by_area : bool, optional
        When computing a vertice normal, weight each triangle normal by the
        area of the triangle so that large triangles contribute more and vice
        versa (default = False).
    
    RETURNS
    ----------
    node_normals : ndarray
        The vertex normals.
    """
    verts2consider = np.arange(len(vertices))
    if not mask is None:
        verts2consider = verts2consider[mask]

    if smooth:
        assert len(verts2consider) == len(vertices),"Cannot smooth normals when using mask"
    v2f = verts2faces(vertices, faces)    
    triangle_normals = get_triangle_normals(vertices[faces])
    node_normals = []
    
    # Get the node normals
    if weigh_by_area:
        mesh = vertices[faces]
        v0 = mesh[:,0,:]
        e1 = mesh[:,1,:]-v0
        e2 = mesh[:,2,:]-v0
            
        area = np.linalg.norm(np.cross(e1,e2),axis=1)/2.
        for n in verts2consider:
            w = area[v2f[n]]
            node_normals.append(np.average(triangle_normals[v2f[n]], axis=0,
                                           weights=w/np.sum(w)))
    else:
        for n in verts2consider:
            node_normals.append(np.average(triangle_normals[v2f[n]], axis=0))
    
    # Normalize
    node_normals  = np.array(node_normals)
    node_normals /= np.sqrt(np.sum(node_normals**2,1))[:,np.newaxis]
    
    # Smooth node normals by averaging neighboring node normals
    if smooth:
        s = np.zeros_like(node_normals)
        for n in verts2consider:
            s[n] = np.average(node_normals[faces[v2f[n]]], axis=(0,1))
                
        node_normals  = s
        node_normals /= np.sqrt(np.sum(node_normals**2,1))[:,np.newaxis]
    
    return node_normals


def ensure_distance(v1, f1, v2, f2, distance, deform="shrink",
                    smooth_mesh=False, search_area=2):
    """Ensure a minimum distance between two surface meshes. To achieve this,
    the first surface is shrunk, i.e., the nodes which are closer than the
    specified distance are moved inwards (opposite direction of the node
    normal).
        
    PARAMETERS
    ----------
    v1 : array_like
        Array describing the vertices of the first mesh. This is the surface
        which is moved to ensure the required distance.
    f1 : array_like
        Array describing the faces of the first mesh. This is the surface
        which is moved to ensure the required distance.
    v2 : array_like
        Array describing the vertices of the second mesh. This surface remains
        untouched.
    f2 : array_like
        Array describing the faces of the second mesh. This surface remains
        untouched. 
    distance : float
        The minimum distance to ensure between the two meshes.
    deform : {"shrink", "expand"}
        Whether to expand or shink the (first) surface. Expand corresponds to
        pushing the nodes outwards (in the direction of the normals) whereas
        shrink will pull the nodes inwards (in the direction of minus the
        normals) to ensure the required distance. Thus, when decoupling
        neighboring surfaces one of them should be shrunk, however, when
        decoupling surfaces where one is fully contained within the other one
        could choose either to shrink the inner surface or expand the outer
        surface. Please note, that it is always the first surface that is
        modified (default = "shrink").
    search_area : float
        To reduce the number of intersections to check for, an initial KDTree
        search is done for neighbors of node within a certain distance. This
        parameter controls the size of the area in which intersections are
        tested and is a multiple of the distance parameter (e.g., if
        search_area is 2 and distance is 0.5 then intersections in a
        neighborhood of 2*0.5=1 are considered). This is a heuristic to reduce
        the computational load and ultimately constitutes a trade-off between
        efficiency and accuracy. If you find that intersections are missed, try
        increasing this parameter.
    
    RETURNS
    ----------
    v1_new : ndarray
        The updated vertices of the first surface mesh.
    """
    lvl = logging.DEBUG
    
    eps = 1e-6
    
    # Find triangles in mesh 2 to test nodes in mesh 1 against
    bar2 = v2[f2].mean(1)
    tree1 = cKDTree(v1)
    tree2 = cKDTree(bar2)
    res = tree1.query_ball_tree(tree2, r=search_area*distance, p=2)
    pi, piok = list2numpy(res, dtype=np.int)
    mask = np.where(piok.any(1))[0]
    # no nearby elements, hence nothing to test or move
    if len(mask) is 0:
        #print("surfaces not close - no vertices moved")
        return v1.copy()
        
    # Check for intersections
    normals = get_node_normals(v1, f1, mask)
    if deform == "expand":
        normals = -normals
    intersect, idists = ray_triangle_intersect(v2[f2], v1[mask], normals,
                            pi[mask], piok[mask], plane_type="one-sided",
                            mindist2int=0, maxdist2int=distance,
                            return_int_dists=True)
    irow, icol = np.where(intersect)
    
    # Shrink the first surface. Nodes are moved inwards such that the desired
    # distance between the surfaces is achieved.
    mm2move = np.zeros((len(v1), idists.shape[1]))
    mm2move[mask[irow],icol] = distance-idists[intersect]+eps
    mm2move = mm2move.max(1)    
    # if nothing to move, return
    if not (mm2move>0).any():
        #print("no need to move any vertices")
        return v1.copy()
        
    # move by maximum 0.25 per iteration
    nsteps = 2 #np.ceil(distance/0.25).astype(int)
    v1_new = deform_mesh(v1, f1, mm2move, distance, nsteps, deform,
                         plane_type="two-sided",smooth_mesh=smooth_mesh)
    log("Moved {} of {} vertices",((v1_new!=v1).sum()/3, len(v1)), level=lvl)
    return v1_new
    

# =============================================================================
# VOLUME MESH UTILITIES
# =============================================================================

def relabel_compartments(mesh, old, new):
    """Relabel elements belonging to compartment 'old' (as defined by tag1 and
    tag2 fields) to compartment 'new' of a .msh file. If 'old' and 'new' are
    lists of indices then they are zipped such that relabelling is from old[0]
    to new[0], from old[1] to new[1] etc.
    
    PARAMETERS
    ----------
    mesh : str or mesh_io msh object
        The mesh whose elements are relabelled.
    old : int or list of ints
        Old label.
    new : int or list of ints
        New label.
        
    RETURNS
    ----------
    mesh : mesh_io msh object
        The modified input object.        
    """
    if type(mesh) == str:    
        mesh = mesh_io.read_msh(mesh)
    
    # ensure iterable
    try:
        next(iter(old))
    except TypeError:
        old = [old]
    try:
        next(iter(new))
    except TypeError:
        new = [new]
    
    # renumber
    for iold, inew in zip(old,new):
        mesh.elm.tag1[mesh.elm.tag1 == iold] = inew
        mesh.elm.tag2[mesh.elm.tag2 == iold] = inew
        
    return mesh
    

def volume_mesh_relabel(mesh, surfaces, region_number):
    """Relabel tetrahedra with region_number inside of the surface(s) defined
    by surfaces and remove any such large (>100) components. The following
    steps are performed:
    
    (1) Labels of tetrahedra with region_number within surfaces are relabelled
        to 99.
    (2) Only the largest component of region_number is retained; all other are
        relabelled to 99,
    (3) Small components (<100 tetrahedra) are relabelled from 99 to
        region_number.
    (4) Large components with tag==99 are removed from the mesh.
    
    This is followed by a clean-up:
    (5) Nodes which, as a consequence of removing the above tetrahedra, are no
        longer used by any elements are removed.
    (6) Elements (triangles) which use any deleted nodes are removed.
    (7) Node and element (triangles, tetrahedra) labels are adjusted to
        accomodate the removed nodes and elements (i.e. reindexing).
    
    Note that after this procedure, the surface mesh of the volume with number
    region_number no longer accurately reflects the boundary of this volume.
    
    PARAMETERS
    ----------
    mesh : mesh_io volume mesh object or str
        Filename of the volume mesh (.msh) or mesh object.
    surfaces : str or list of strs
        An array of filenames of the air components ...
    region_number : int
        
    RETURNS
    ----------
    mesh : mesh_io volume mesh object
        The mesh with corrected tetrahedra labels.
    """
    if type(mesh) == str:    
        mesh = mesh_io.read_msh(mesh)
    if type(surfaces) == str:
        surfaces = [surfaces]
        
    lvl = logging.DEBUG
    
    # indexing in .msh files start from 1. To minimize confusion, subtract 1 to
    # get coordinates appropriate for python indexing. Readd before saving the
    # mesh to disk
    mesh.elm.node_number_list -= 1
    
    # adjust mesh such that
    log("Adjusting region numbers of tetrahedra", level=lvl)
    for f in surfaces:
        surf_verts, surf_faces = mesh_load(f)
        mesh = mesh_adjust_tetrahedra_regions(mesh, surf_verts, surf_faces, region_number, 99)
    
    # get connected skin components and relabel all but the largest component 
    # to air
    log("Getting connected components and relabelling", level=lvl)
    toi = np.where((mesh.elm.elm_type==4) & (mesh.elm.tag1==region_number))[0] # tetrahedra of interest
    
    if len(toi) > 0:    
        tetrahedra = mesh.elm.node_number_list[toi]
        tetrahedra = mesh.nodes.node_coord[tetrahedra]
    
        tetrahedra_labels = get_connected_components(tetrahedra)
        mesh.elm.tag1[toi[tetrahedra_labels>0]] = 99
     
        # relabel small air components (<100 elements) to skin
        toi = np.where((mesh.elm.elm_type==4) & (mesh.elm.tag1==99))[0]
        
        if len(toi) > 0:
            tetrahedra = mesh.elm.node_number_list[toi]
            tetrahedra = mesh.nodes.node_coord[tetrahedra]
        
            tetrahedra_labels, components, counts = get_connected_components(tetrahedra,
                                                                             return_components_and_counts=True)
            mesh.elm.tag1[toi[np.in1d(tetrahedra_labels, components[counts<100])]] = region_number
            
    # remove large air components, i.e. elements where tag1 == 99
    log("Removing large air components", level=lvl)
    rem = mesh.elm.tag1 != 99
    mesh.elm.node_number_list = mesh.elm.node_number_list[rem]
    mesh.elm.elm_type         = mesh.elm.elm_type[rem]
    mesh.elm.tag1             = mesh.elm.tag1[rem]
    
    # remove points not used by any tetrahedra
    log("Removing nodes not used by any tetrahedra", level=lvl)
    pts = np.unique(mesh.elm.node_number_list[mesh.elm.elm_type==4].ravel())    
    pts_retained = np.zeros(mesh.nodes.nr, dtype=bool)
    pts_retained[pts] = True

    mesh.nodes.node_coord      = mesh.nodes.node_coord[pts]
    
    # remove triangles which use any removed points
    log("Removing triangles using any deleted nodes", level=lvl)
    rm_triangles = ~np.all(pts_retained[mesh.elm.node_number_list], axis=1)
    rem = (mesh.elm.elm_type==2) & rm_triangles
    mesh.elm.node_number_list = mesh.elm.node_number_list[~rem]
    mesh.elm.elm_type         = mesh.elm.elm_type[~rem]
    mesh.elm.tag1             = mesh.elm.tag1[~rem]
    mesh.elm.tag2             = mesh.elm.tag1
    
    # adjust point indices of tetrahedra and triangles
    hidx      = -np.ones(pts.max()+1)
    hidx[pts] = np.arange(len(pts))
    hidx      = np.append(hidx,-1)
    
    mesh.elm.node_number_list = hidx[mesh.elm.node_number_list].astype(np.int32)
    
    # add one to all indices to conform to the .msh format (indexing starts at
    # one) before saving
    mesh.elm.node_number_list += 1 
    
    # filename
    fpath,fname = os.path.split(mesh.fn)
    mesh.fn = add2filename(mesh.fn, "_relabel")
    
    return mesh

    
def get_connected_components(elements, return_components_and_counts=False,
                             ntol=1e-6):
    """Given a set of elements (e.g., triangles, tetrahedra), find connected
    components in this set.
    
    
    PARAMETERS
    ----------
    elements : ndarray
        Array of elements (e.g., triangles or tetrahedra) described as an
        N-by-M, with N being number of elements and M being the number of
        vertices of each element (e.g., 3 and 4 for triangles and tetrahedra,
        respectively).        
    return_components_and_counts : bool, optional
        Return component labels/numbers and corresponding counts of elements
        belonging to each component. Please note that these arrays are NOT
        sorted in any way.
    ntol : float, optional
        Neighbor tolerance. This parameters controls the upper bound for when
        elements are considered neighbors, i.e. the distance between
        elements has to be smaller than this value (default = 1e-6)
    
    RETURNS
    ----------
    element_labels : ndarray
        N-by-1 vector mapping each element to a component.
    labels : ndarray, optional
        Vector of component labels/numbers.
    counts : ndarray, optional
        Vector of counts corresponding to each component label/number.
    """
    import scipy.sparse.csgraph
    elements_idx = np.arange(len(elements))

    nn, ok = get_element_neighbors(elements, ntol)
    
    # Make sparse, binary adjacency matrix and get connected components
    data = [1]*ok.sum()
    row_idx = elements_idx.repeat(ok.sum(1))
    col_idx = nn.ravel()[ok.ravel()]
    adj_matrix = scipy.sparse.csr_matrix((data, (row_idx, col_idx)), shape=(nn.shape[0],nn.shape[0])) # AT: enforce shape in case last element is not connected to any other element 
    
    element_labels = scipy.sparse.csgraph.connected_components(adj_matrix, directed=False)[1]
    
    if return_components_and_counts:
        labels, counts = np.unique(element_labels, return_counts=True)
        return element_labels, labels, counts
    else:
        return element_labels

    
def get_element_neighbors(elements, ntol=1e-6):
    """Get the neighbors of each element in elements by comparing barycenters
    of element faces (e.g., if elements are tetrahedra, the faces are
    triangles).
    
    PARAMETERS
    ----------
    elements : ndarray
        Array of elements (e.g., triangles or tetrahedra) described as an
        N-by-M, with N being number of elements and M being the number of
        vertices of each element (e.g., 3 and 4 for triangles and tetrahedra,
        respectively).
    ntol : float, optional
        Neighbor tolerance. This parameters controls the upper bound for when
        elements are considered neighbors, i.e. the distance between elements
        has to be smaller than this value (default = 1e-6).   
    
    RETURNS
    ----------
    nearest_neighbors : ndarray
        N-by-M array of indices into elements, i.e. for each face, which is its
        neighboring element.
    ok : ndarray (bool)
        This array tells, for each entry in nearest_neighbors, if this is an
        actual neighbor or not. The nearest neighbors are returned as a numpy
        ndarray of shape elements.shape for ease of interpretation and 
        efficiency (and not for example as a list of lists of [possibly]
        unequal lengths), hence this is needed.
    """
    
    elements_idx = np.arange(len(elements))
    
    # barycenters of the faces making up each element    
    barycenters = np.zeros_like(elements)
    num_nodes_per_el = elements.shape[1]
    for i in range(num_nodes_per_el):
        nodes = np.roll(np.arange(num_nodes_per_el),-i)[:-1] # nodes that make up the ith face
        barycenters[:,i,:] = np.average(elements[:,nodes,:], 1)
    
    bar_tree = cKDTree(barycenters.reshape(np.multiply(*elements.shape[:-1]),
                                           elements.shape[-1]))
    face_dist, face_idx = bar_tree.query(bar_tree.data, 2)
    
    nonself = (face_idx != np.arange(len(face_idx))[:,np.newaxis]) # get non-self-references
    
    # Distance to nearest neighbor. Neighbors having a distance shorter than
    # ntol are considered actual neighbors (i.e. sharing a face)
    face_dist = face_dist[nonself] 
    ok = face_dist < ntol
    ok = ok.reshape(elements.shape[:2])
    
    # Index of nearest neigbor. From the tree search, indices are to nearest
    # element face, however, we wish to find neighboring elements. Hence,
    # reindex.
    face_idx = face_idx[nonself]
    nearest_neighbors = elements_idx.repeat(num_nodes_per_el)[face_idx]
    nearest_neighbors = nearest_neighbors.reshape(elements.shape[:2])
    
    return nearest_neighbors, ok


def mesh_adjust_tetrahedra_regions(vmesh, surf_verts, surf_faces, old_reg_num,
                                   new_reg_num):
    """Find tetrahedra having region number old_reg_num within the surface
    spanned by the surface mesh defined by surf_verts and surf_faces and set
    the region number of these elements to new_reg_num.
    
    PARAMETERS
    ----------
    vmesh : mesh_io volume mesh object
        Object describing the volume mesh.
    surf_verts : ndarray
        Array containing the vertices of the surface mesh.
    surf_faces : ndarray        
        Array of faces in the surface mesh (i.e. an array indexing into
        surf_verts)
    old_reg_num : int
        The region within which to search for elements within the surface.
    new_reg_num : int
        Elements having region number old_reg_num and which are contained
        within the surface spanned by the surface mesh are reassigned to this
        region (number).
    
    RETURNS
    ----------
    vmesh : mesh_io mesh object
        The volume mesh with updated region numbers.
    """
    # get tetrahedra centers
    tetrahedra_idx = vmesh.elm.elm_number[
        (vmesh.elm.elm_type==4) & (vmesh.elm.tag1==old_reg_num)] - 1
    tetCenters = np.mean(vmesh.nodes.node_coord[vmesh.elm.node_number_list[tetrahedra_idx]], axis=1)
    
    # bounding box spanned by surface mesh
    surfmin = surf_verts.min(0)
    surfmax = surf_verts.max(0)
    
    # retain only tetrahedra of the specified mesh within the bounding of the 
    # surface mesh
    centerindices = np.where(np.all((tetCenters-surfmin)*(surfmax-tetCenters) > 0,axis=1))[0]
    
    inMesh,cList = in_mesh(surf_verts[surf_faces], tetCenters[centerindices,:],
                           quiet=True)
    
    # set new region numbers
    vmesh.elm.tag1[tetrahedra_idx[centerindices[inMesh]]] = new_reg_num
    
    return vmesh


def vmesh2nifti(mesh, vol, filename):
    """Given a volume mesh and an image file (e.g., .nii file), determine which
    voxels in vol are inside which elements of the mesh and label according to
    the values found in mesh.elm.tag1.
    
    PARAMETERS
    ----------
    mesh : str or mesh_io volume mesh object
        Object describing the volume mesh.
    vol : str or nibabel image object
        Reference volume. The output volume is similar to this with respect to
        image dimensions, voxel size etc. only with voxel values being replaced
        by corresponding values of the labels of the tetrahedra in the mesh.
    filename : str
        Name of the saved volume.
        
    RETURNS
    ----------
    Nothing, saves a nifti file by the name of filename to disk.
    """
    if type(mesh) == str:
        mesh = mesh_io.read_msh(mesh)
    if type(vol) == str:
        vol = nib.load(vol)
        
    img_dims = vol.shape
    
    points = np.array(np.meshgrid(*tuple(map(np.arange,img_dims)),
                                      indexing="ij")).reshape((3,-1)).T

    tetrahedra = mesh.nodes.node_coord[mesh.elm.node_number_list[mesh.elm.elm_type==4]-1]
    tetrahedra_regions = mesh.elm.tag1[mesh.elm.elm_type==4]
    
    tetrahedra = apply_affine(tetrahedra, np.linalg.inv(vol.affine))

    labels = label_points(tetrahedra, points, tetrahedra_regions)
    
    write_nifti(labels.reshape(img_dims), filename, vol, dtype=np.uint8,
                return_volume=False)


def label_points(tetrahedra, points, tetrahedra_labels=None, af=20, nsplits="auto"):
    """Find which points are inside which tetrahedra and label according to the
    values of tetrahedra_labels.
    
    PARAMETERS
    ----------
    tetrahedra : ndarray
        Array defining the tetrahedra.
    points : ndarray
        N x 3 array. Each row corresponds to a point. Columns corresponds to 
        coordinates.
    tetrahedra_labels : narray_like, optional
        Label corresponding to each tetrahedra. These values are assigned to
        points inside a given tetrahedra (default = [], in which case the index
        of the tetrahedra within which each point is found, is returned)
    af : int, optional
        Acceleration factor. The algorithm splits the tetrahedra and points
        into splits which is processed sequentially. It makes two passes when
        determining in which tetrahedra a point is located: (1) an approximate
        test in which it considers only the nearest tetrahedra of each point
        and checks these, and (2) an exact test in which remaining points are
        tested against all tetrahedra which it could possibly be located
        inside [and which were not already excluded based on (1)]. The
        acceleration factor controls how many neighbors to consider in (1). A
        small value will result in highly variable times to process each split
        (depending on the number of tetrahedra and testpoints). A high value
        will result in long computation time for splits with a large number of
        tetrahedra and/or points. An optimal value will ensure fast processing,
        even of large splits, hence decreasing the variability in computation
        time for each split (default = 20).    
    
    RETURNS
    ----------
    labels : ndarray
        Label of each point as given by tetrahedra_labels. If no labels are
        supplied, contains indices into tetrahedra.
    """
    assert isinstance(af,int)
    if tetrahedra_labels is None:
        tetrahedra_labels = []
    
    if len(tetrahedra_labels) == 0:
        tetrahedra_labels = np.arange(len(tetrahedra))
    if len(tetrahedra_labels) != len(tetrahedra):
        raise ValueError("Length of tetrahedra_labels must equal that of tetrahedra.")
    
    # prepare output
    point_labels = np.zeros(len(points))
    
    normals, tet_bars = get_tetrahedra_normals(tetrahedra,
                                               return_element_centers=True)
    
    msh_min = tetrahedra.min(axis=(0,1)) # min and max coordinates of the
    msh_max = tetrahedra.max(axis=(0,1)) # entire mesh
          
    # select only points within a rectangular bounding box of the tetrahedra
    bbox = np.all((points >= msh_min) & (points <= msh_max), axis=1)

    splits, split_size = make_splits(tetrahedra, points, nsplits)
    
    for s in splits:
        # GET TETRAHEDRA AND TEST POINTS IN SLAB
        # =======================================
        progress_bar(s, splits[-1], 50)
    
        split_zmin = s-split_size/2.
        split_zmax = s+split_size/2.
  
  
        # get tetrahedra in split (i.e. all but the ones which only have
        # vertices above/below s +/- 0.5)
        stets = ~((tetrahedra[...,2] < split_zmin).all(1) | (tetrahedra[...,2] > split_zmax).all(1))
        
        if stets.sum() == 0: # if no tetrahedra in split
            continue
        
        split_tetrahedra = tetrahedra[stets]
    
        smin = split_tetrahedra.min(axis=(0,1))
        smax = split_tetrahedra.max(axis=(0,1))
        smin[2] = split_zmin
        smax[2] = split_zmax
        
        # indices of test points in the current split
        split_bbox = np.where(bbox & np.all((points >= smin) & (points <= smax),
                                           axis=1))[0]
    
        if len(split_bbox) == 0: # if no points to test in the current split
            continue            
        
        # tetrahedra to be considered
        split_normals = normals[stets]                
        split_bars    = tet_bars[stets]
          
        # TEST POINTS AGAINST TETRAHEDRA
        # =======================================
        # Get max distance from tetrahedron barycenter to vertice. This defines
        # the upper bound for the distance between a test point and tetrahedra
        r_max = np.linalg.norm(split_tetrahedra-split_bars[:,np.newaxis,:],
                               ord=2, axis=2).max()
        
        # for each point to test: which tetrahedra are closest
        bars_tree = cKDTree(split_bars)
        
        # Make two passes, first testing points only against "af" nearest
        # neighbors and subsequently check the remaining points against all
        # tetrahedra which it could possibly be inside (as determined by r_max)
        for check in range(2):

            spoints_idx  = points[split_bbox]
            
            if check == 0: # approximate
                n_dist, n_idx = bars_tree.query(spoints_idx, af,
                                                distance_upper_bound=r_max)
                ok = n_dist != np.inf
                oko = ok.copy()
            else: # exact
                xx = bars_tree.query_ball_point(spoints_idx, r_max, p=2)
                n_idx, ok = list2numpy(xx, dtype=np.int)
       
            # For each face in the tetrahedra:
            # (1) get the vector from the first node in the ith face to the
            #     test points
            # (2) get dot product of normal of the ith face and the vector from
            #     (1)         
            # (3) if the sign of the dot product is negative then the point is
            #     on the correct side of the plane (spanned by the vectors of
            #     the ith face)
            for i in range(4):
                i2p    = spoints_idx.repeat(ok.sum(1), axis=0)-split_tetrahedra[n_idx[ok],i,:]
                #dpoint = (split_normals[n_idx[ok],i,:]*(i2p)).sum(1)
                dpoint = np.einsum("ij,ij->i",split_normals[n_idx[ok],i,:], i2p)
                ok[ok] = (dpoint <= 0)
            ok_points, ok_tetrahedra = np.where(ok)
                 
            point_labels[split_bbox[ok_points]] = tetrahedra_labels[stets][n_idx[ok_points, ok_tetrahedra]]
            
            # Get the points which had valid entries in ok originally but which
            # were not found to be inside any tetrahedra. Test these again, now
            # exact.
            if check == 0:
                split_bbox = split_bbox[(oko.sum(1)>0) ^ (ok.sum(1)>0)]
                if len(split_bbox) == 0:
                    break

    return point_labels

        
def get_tetrahedra_normals(tetrahedra, return_element_centers=False,
                           return_face_centers=False):
    """Get the vector normal to each face of all tetrahedra and ensure that
    they point outwards.

    PARAMETERS
    ----------
    tetrahedra : ndarray
        Array describing the tetrahedra with dimensions:
            [# tetrahedra] x [# vertices] x [x y z]
    return_element_centers : bool, optional
        Return barycenter of all tetrahedra.
    return_face_centers : bool, optional
        Return barycenter of all faces of all tetrahedra
        
    RETURNS
    ----------
    normals : ndarray
        Normal vectors of each face of all tetrahedra.
    tet_barycenters : ndarray, optional
        Coordinates of barycenters of all tetrahedra
    tri_barycenters : ndarray, optional
        Coordinates of barycenters of all faces of all tetrahedra        
    """    
    
    tet_barycenters = np.average(tetrahedra, axis=1)
    normals, tri_barycenters = np.zeros_like(tetrahedra), np.zeros_like(tetrahedra)
    
    # For each face of the tetrahedra, calculate normal
    # (the faces are given by [0,1,2], [1,2,3], [2,3,0], [3,0,1])
    for i in range(tetrahedra.shape[1]):
        nodes = np.roll([0,1,2,3],-i)[:3] # nodes that make up the ith face
        normals[:,i,:] = np.cross(tetrahedra[:,nodes[1],:]-tetrahedra[:,nodes[0],:],
                                  tetrahedra[:,nodes[2],:]-tetrahedra[:,nodes[0],:])
        tri_barycenters[:,i,:] = np.average(tetrahedra[:,nodes,:], 1)
    
    # Ensure normals are pointing outwards by reversing normals pointing in the
    # same direction as the vector from the face barycenter to the tetrahedral
    # barycenter (i.e. inwards)
    reverse = ((normals*(tet_barycenters[:,np.newaxis,:]-tri_barycenters)).sum(2)>=0)
    normals[reverse] *=-1
    
    # Normalize
    normals /= np.sqrt((normals*normals).sum(axis=2))[...,np.newaxis]
    
    if return_element_centers & return_face_centers:
        return normals, tet_barycenters, tri_barycenters
    elif return_element_centers:
        return normals, tet_barycenters
    elif return_face_centers:
        return normals, tri_barycenters
    else:    
        return normals
		
    
# =============================================================================
# MISCELLANEOUS
# =============================================================================

def apply_affine(points, affine):
    """Apply 4-by-4 affine transformation to points in 3D. 
    
    PARAMETERS
    ----------
    points : array_like
        Points to transform. The last axis defines the coordinates in 3D, e.g.,
        points is an n-by-3, n-by-m-by-3, etc. array where all but the last
        dimension is assumed to define points.
    affine : array_like
        The affine transformation to apply. It is right-multiplied to points.
    
    RETURNS
    ----------
    p : array_like
        Transformed points.    
    """
    
    dim = points.shape
    assert dim[-1] == 3
    assert affine.shape == (4,4)
    
    if len(dim) == 1:
        points = points[None,:]
    elif len(dim) > 2:
        points = points.reshape((-1,3))
    
    # append column of ones and apply transformation retaining original array
    # size
    p = np.concatenate((points, np.ones(len(points))[:,np.newaxis]), axis=1)
    p = p.dot(affine.T)[:,:3]
    
    if len(dim) == 1 or len(dim) > 2:
        p = p.reshape(dim)
        
    return p


def vox_idx2msh_coord(vox_idx,vox_dims,img_dims):
    """Convert LAS sorted volume voxel indices to RAS sorted surface mesh
    coordinates. Voxel dimensions may be calculated using get_vx_size or
    
        vox_dims = np.linalg.norm(M[:3,:3],ord=2,axis=0),
        
    where M is the affine transformation matrix from voxel indices to 
    coordinates (e.g., scanner or MNI). This corresponds to matrix 
    multiplication by the transformation matrix
    
        -Vx   0   0   Vx*Nx/2 
          0  Vy   0  -Vy*Ny/2 
          0   0  Vz   Vz*(-Nz/2+1) 
          0   0   0      1
      
    Thus, if vox_idx is a num_voxels-by-three array, then
    
        vertices = np.dot(np.concatenate((vox_idx,np.ones(len(vox_idx))[:,np.newaxis]),axis=1 ), M.T)[:,:3]
    
    where M now denotes the transformation matrix described above.
      
    PARAMETERS
    ----------
    vox_idx : ndarray
        Voxel indices to convert (N x 3).
    vox_dims : ndarray
        Array of voxel dimensions, [Vx,Vy,Vz].
    img_dims : ndarray
        Array of image dimensions, [Nx,Ny,Nz].
    
    RETURNS
    ----------
    vertices : ndarray
        Surface mesh coordinates.
    """    
    vox_idx  = np.array(vox_idx)
    vox_dims = np.array(vox_dims)
    img_dims = np.array(img_dims)
    
    # voxel indices (in LAS) to mesh coordinates (in RAS)
    vertices = vox_dims * ([-1,1,1]*vox_idx + ([1,-1,-1]*img_dims/2. + [0,0,1]))
    
    return vertices

    
def msh_coord2vox_idx(vertices, vox_dims, img_dims):
    """Convert RAS sorted surface mesh coordinates to LAS sorted volume voxel
    indices. Voxel dimensions may be calculated using get_vx_size or
    
        vx_dims = np.linalg.norm(M[:3,:3],ord=2,axis=0),
        
    where M is the affine transformation matrix from voxel indices to 
    coordinates (e.g., scanner or MNI). This corresponds to matrix 
    multiplication by the transformation matrix
    
        -1/Vx   0     0    Nx/2 
          0    1/Vy   0    Ny/2 
          0     0    1/Vz  Nz/2-1 
          0     0     0      1
      
    Thus, if vertices is a num_vertices-by-three array, then
    
        vox_idx = np.dot(np.concatenate((vertices,np.ones(len(vertices))[:,np.newaxis]),axis=1 ),M.T )[:,:3]
    
    where M now denotes the transformation matrix described above.
    
    PARAMETERS
    ----------
    vertices : ndarray
        Vertices describing the mesh (N x 3).
    vox_dims : ndarray
        Array of voxel dimensions, [Vx,Vy,Vz].
    img_dims : ndarray
        Array of image dimensions, [Nx,Ny,Nz].
    
    RETURNS
    ----------
    vox_idx : ndarray
        Voxel indices corresponding to the vertex coordinates.
    """
    vertices = np.array(vertices)
    vox_dims = np.array(vox_dims)
    img_dims = np.array(img_dims)
    
    # mesh coords (in RAS) to vox indices (in LAS)
    vox_idx = [-1,1,1]*vertices/vox_dims + img_dims/2. - [0,0,1]
    
    return vox_idx

    
def remove_duplicates(array_data, return_index=False, return_inverse=False):
    """Remove duplicate rows of a multi-dimensional array. Return the array 
    with the duplicates removed. If return_index is True, also return the 
    indices of array_data that result in the unique array. If return_inverse is
    True, also return the indices of the unique array that can be used to
    reconstruct array_data.

    PARAMETERS
    ----------
    array_data : ndarray
        Array from which to remove duplicate rows.
    return_index : bool, optional
        Return indices into array_data from which to create unique_array_data
        (default = False).
    return_inverse : bool, optional
        Return indices into unique_array_data that may be used to recreate
        array_data (default = False).
        
    RETURNS
    ----------
    unique_array_data : ndarray
        Array containing only unique rows of array_data.
    """ 
    unique_array_data, index_map, inverse_map = np.unique(
            array_data.view([("", array_data.dtype)] * \
                    array_data.shape[1]), return_index=True,
                    return_inverse=True)
    
    unique_array_data = unique_array_data.view(
            array_data.dtype).reshape(-1, array_data.shape[1])

    # unique returns as int64, so cast back
    index_map   = np.cast["uint32"](index_map)
    inverse_map = np.cast["uint32"](inverse_map)
    
    if return_index and return_inverse:
        return unique_array_data, index_map, inverse_map
    elif return_index:
        return unique_array_data, index_map
    elif return_inverse:
         return unique_array_data, inverse_map
    
    return unique_array_data


def row_wise_unique(arr, fill_value=0):
    """Return row-wise unique values of an array.
    
    The trick is to add a unique imaginary number to each row. That way, when 
    calling np.unique, the floats in the original array will be recognized as 
    different values if they occur in different rows, but be treated as the 
    same value if they occur in the same row.

    PARAMETERS
    ----------
    arr : ndarray
        Array in which to find unique values.

    RETURNS
    ----------
    u_arr : ndarray
        Array in which unique values are retained whereas all non-unique
        elements are set to zero.
                    
    Acknowledgements
    ----------------
    From https://stackoverflow.com/questions/26958233/numpy-row-wise-unique-elements
    (by user "unutbu").
    """    
    weight = 1j*np.linspace(0, arr.shape[1], arr.shape[0],
                            endpoint=False)      # row "weights"
    u_arr  = arr + weight[:,np.newaxis]          # add weights to rows
    u, ind = np.unique(u_arr, return_index=True) # now get unique values
    u_arr  = np.ones_like(arr)*fill_value        # initialize array
    np.put(u_arr, ind, arr.flat[ind]) # fill in unique values. Remaining values
                                     # are set to zero.
    return u_arr

def spawn_process(cmd, return_exit_status=False, new_thread=False,
                  verbose=False, shell=False, new_process=False):
    """Spawn a new process and communicate its output.
    
    PARAMETERS
    ----------
    cmd : str
        The command to be executed.
    return_exit_status : bool, optional
        Return exit status of process. Only applies if new_thread == False
        (default = False).
    new_thread : bool, optional
        By default, the child process blocks the main process while running,
        however, by setting this option to true the child progress will be
        spawned in a new thread, hence not blocking the main process (default =
        False).
    verbose : bool, optional
        Communicate stdout to console or only to logfile. Only applies if 
        new_thread == False (default = False).
    shell: bool, optional
        Whether to use shell mode. Default: False (forced to be True on Windows)
    new_process: bool, optional
        Starts command in a new process. Default: False
    RETURNS
    ----------
    p.returncode : int (optional)
        Return code of the process.
    """
    log("Running {0}",cmd, level=logging.DEBUG)
    if verbose:
        lvl = logging.INFO
    else:
        lvl = logging.DEBUG

    if new_thread or new_process:
        ON_POSIX = "posix" in sys.builtin_module_names

        def enqueue_output(out, queue):
            for line in iter(out.readline, b''):
                try:
                    queue.put(line.decode())
                except UnicodeDecodeError:
                    queue.put('Could not print line')
            out.close()

        p = Popen(cmd, stdout=PIPE,
                  bufsize=1, close_fds=ON_POSIX,
                  shell=True)

        q = Queue()
        if new_thread:
            t = Thread(target=enqueue_output, args=(p.stdout, q))
        if new_process:
            t = Process(target=enqueue_output, args=(p.stdout, q))
        t.daemon = True  # thread dies with the program
        t.start()
        return t

    else: 
        p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        for line in iter(p.stdout.readline, b''):
            # to prevent massive output to the logfile retain only the last
            # line starting with \r
            try:
                line = line.decode().split("\r")[-1]
            except UnicodeDecodeError:
                log('Could not print line', level=logging.DEBUG)
                continue
            
            # remove one \n since logger.log adds one itself
            if line[-1]=="\n":
                line = line[:-1]
            
            try:
                log(line, level=lvl)
            except KeyError: # in case line includes ".. {..} .."
                log(line.replace("{","{{").replace("}","}}"), level=lvl)
                
        if return_exit_status:
            return p.wait()

    
def format_time(running_time):
    """Format time in seconds as hours:minutes:seconds.
    
    PARAMETERS
    ----------
    running_time : float
        Time in seconds.
    
    RETURNS
    ----------
    running_time : str
        The time formatted as hours:minutes:seconds.
    """
    hrs = np.uint16(np.floor(running_time/(60.**2)))
    mts = np.uint16(np.floor(running_time/60.-hrs*60))
    sec = np.uint16(np.round(running_time-hrs*60.**2-mts*60.))

    return "{:02d}:{:02d}:{:02d}".format(hrs,mts,sec)

    
def write_unity_xfm(output_dir):
    """Write a unity transformation matrix to disk.
    
    PARAMETERS
    ----------
    output_dir : str
        The directory in which to save the generated file. 
    
    RETURNS
    ----------
    Writes the file "unity.xfm" to the specified folder.
    """
    with open(os.path.join(output_dir, "unity.xfm"),"w") as f:
        f.write(textwrap.dedent("""\
        MNI Transform File
        % tkregister2
        
        Transform_Type = Linear;
        Linear_Transform =
        1.0    0.0    0.0  0.0
        0.0    1.0    0.0  0.0
        0.0    0.0    1.0  0.0;"""))


def write_ref_fs(vol):
    """Write a FreeSurfer reference file to disk. If input format is LAS then 
    output format is LIA (from LAS to LIA: x,-z,y).
    
    PARAMETERS
    ----------
    vol : nibabel image object
        Reference volume from which the file is generated. The actual data of
        the volume is not important.
    
    RETURNS
    ----------
    Write the file "ref_FS.nii.gz" to the same folder as vol.
    """
    output_dir = os.path.dirname(vol.get_filename())
    data   = np.transpose(np.zeros_like(vol.get_data()),(0,2,1))
    # alternatively, if we wanted to write out the data
    #data   = np.transpose(vol.get_data(),(0,2,1))[:,::-1,:]
    affine = vol.affine[[0,2,1,3]]*np.array([1,1,-1,1])[:,np.newaxis]
    vx = get_vx_size(vol)
    dim = vol.shape
    affine = np.array([[-vx[0],     0 ,    0 ,  vx[0]*dim[0]/2.],
                       [    0 ,     0 , vx[1], -vx[1]*dim[1]/2.],
                       [    0 , -vx[2],    0 ,  vx[2]*dim[2]/2.],
                       [    0 ,     0 ,    0 ,                1]])
    ref = nib.Nifti1Image(data.astype(np.uint16),affine)
    ref.set_qform(affine)
    ref.header.set_xyzt_units(*vol.header.get_xyzt_units())
    nib.save(ref, os.path.join(output_dir, "ref_FS.nii.gz"))

    
def write_nifti(data, filename, reference_vol, file_format="nii.gz",
                dtype=None, write_volume=True, return_volume=False):
    """Write data as a nifti file. Compared to nibabel.save() this also sets 
    q-form and units from a refernce volume.
    
    PARAMETERS
    ----------
    data : ndarray
        The image data.
    filename : str
        Name of the file to write.
    reference_vol : nibabel image object
        Reference volume. The header information from this volume is copied and
        used when writing the new volume.
    file_format : {"nii","nii.gz"}, optional
        If filename is not specified with extension ("nii" or "nii.gz"), the
        file is saved using this format (default = "nii.gz")
    dtype : datatype, optional
        Sets the datatype of data before saving.
    write : bool, optional
        Write the data to filename (default = True).
    return_volume : bool, optional
        Return a nibabel image object containing the data and the requested 
        header information (default = False).
        
    RETURNS
    ----------
    Writes the volume to disk.
    
    v : nibabel image object
        The image which was saved.
    """
    if type(reference_vol) == str:
        reference_vol = nib.load(reference_vol)
    assert file_format == "nii" or file_format == "nii.gz"
    if dtype:
        data = data.astype(dtype)
        
    v  = nib.Nifti1Image(data, reference_vol.affine)
    v.set_qform(reference_vol.affine)
    v.header.set_xyzt_units(*reference_vol.header.get_xyzt_units())
    
    if filename.endswith(".nii") or filename.endswith(".nii.gz"):
        v.set_filename(filename)
    else:
        v.set_filename(filename+"."+file_format)
        #v.set_filename(filename.split(".")[0]+"."+file_format)
    
    if write_volume:
        nib.save(v, v.get_filename())
    if return_volume:
        return v


def progress_bar(current, end, bar_size=50, print_progress=True, symbol="#"):
    """Print a crude progress bar which measures progress in terms of number of
    loops completed.
    
    PARAMETERS
    ----------
    current : float
        Current step number in loop.
    end : float
        Number of total iterations of loop.
    bar_size : int, optional
        Length of progress bar (default = 50).
    print_progress : bool, optional
        Print progress to the command line. If false, return a string
        containing the progress bar (default = True).
    Symbol : str, optional
        Which symbol to use for the progress bar (default = "#").
        
    RETURNS
    ----------
    progress : str, optional
        Either print progress directly to the command line or return a string
        containing the progress bar.
    """
    percent_complete = current/np.float(end)*100
    percent_complete_int = np.round(percent_complete).astype(np.int)
    
    nhexs = np.round(percent_complete/(100/bar_size)).astype(np.int)

    progress = "|"+symbol*nhexs+" "*(bar_size-nhexs)+"|"
    progress += "{:3d}% Completed".format(percent_complete_int)
    
    if current == end:
        progress += "\n"
        
    if print_progress:
        sys.stdout.write(("\r"+progress))
        sys.stdout.flush()
    else:
        return progress
    
    
def log(msg, args=None, level=logging.INFO, width=72):
    """Log message with formatting.
    
    PARAMETERS
    ----------
    msg : str
        The message to log.
    args : list or str, optional
        List of arguments. Arguments are applied using .format, hence supplying
        arguments, msg should contain proper placeholders for these arguments
        (e.g., {0}, {1}, {:d} or with options such as {:06.2f}) (default = [],
        i.e. no arguments).        
    width : int, optional
        Wrapping width of msg (default = 72)
        GUILHERME: deprecated this argument, it only makes the log more confusing
    level : logging level, optional
        Specify logging level (default = logging.INFO)
        
    RETURNS
    ----------
    The message logged on logger named "logger".    
    """
    if args is None:
        args = []
    if type(args) not in [list,tuple]:
        args = [args]
        
    try:
        logger.log(level, msg.format(*args))
    # If there is an error in the message formating, logs the raw message
    except ValueError:
        logger.log(level, msg)


def add2filename(filename, text, where="end", useext=None):
    """Add text to filename, either at the front or end of the basename.
    Everything after the first punctuation mark is considered file extension.
    """
    
    where = where.lower()
    assert where in ["end","front"]
     
    directory = os.path.dirname(filename)
    basename = os.path.basename(filename).split(".")[0]
    ext = ".".join(os.path.basename(filename).split(".")[1:])

    if useext:
        ext = useext
        
    if where == "end":
        return os.path.join(directory,basename+text+"."+ext)
    elif where == "front":
        return os.path.join(directory,text+basename+"."+ext)


def read_curv(fname):
    """Read a FreeSurfer curvature file. This filetype associates one value
    with each node (e.g., gray matter thickness). Uses big-endian format.
    
    PARAMETERS
    ----------
    fname : str
        Name of the file to be read.
    
    RETURNS
    ----------
    data : ndarray
        The data in the file.
    """
    # if first three bytes equal this number, then it is the "new" curvature
    # format    
    NEW_VERSION_MAGIC_NUMBER = 16777215
    
    with open(fname,"rb") as f:
        n_vertices = read_int3(f)
        if n_vertices == NEW_VERSION_MAGIC_NUMBER:
            n_vertices, n_faces, vals_per_vertex = np.fromfile(f,">i4",3)
            data = np.fromfile(f, ">f4", n_vertices)
        else:
            read_int3(f) # number of faces
            data = np.fromfile(f, ">i4", n_vertices)/100.
    
    return data


def read_int3(f):
    """Read three bytes from the current position in the file as an integer,
    i.e. int24.
    
    PARAMETERS
    ----------
    f : file object
        The file object from which to read.
        
    RETURNS
    ----------
    val : int
        The integer read from three bytes in the file.
    """
    b1 = struct.unpack('B',f.read(1))[0]
    b2 = struct.unpack('B',f.read(1))[0]
    b3 = struct.unpack('B',f.read(1))[0]

    val = (b1 << 16) + (b2 << 8) + b3
    
    return val


"""
All the functions bellow are
Copyright (C) 2011, the scikit-image team All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

        Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
        Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
        Neither the name of skimage nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

def marching_cubes(volume, level=None, spacing=(1., 1., 1.),
                   gradient_direction='descent', step_size=1,
                   allow_degenerate=True, use_classic=False):
    ''' See the documentation for marching_cubes_lewiner '''
    # Check volume and ensure its in the format that the alg needs
    if not isinstance(volume, np.ndarray) or (volume.ndim != 3):
        raise ValueError('Input volume should be a 3D numpy array.')
    if volume.shape[0] < 2 or volume.shape[1] < 2 or volume.shape[2] < 2:
        raise ValueError("Input array must be at least 2x2x2.")
    volume = np.ascontiguousarray(volume, np.float32)  # no copy if not necessary

    # Check/convert other inputs:
    # level
    if level is None:
        level = 0.5 * (volume.min() + volume.max())
    else:
        level = float(level)
        if level < volume.min() or level > volume.max():
            raise ValueError("Surface level must be within volume data range.")
    # spacing
    if len(spacing) != 3:
        raise ValueError("`spacing` must consist of three floats.")
    # step_size
    step_size = int(step_size)
    if step_size < 1:
        raise ValueError('step_size must be at least one.')
    # use_classic
    use_classic = bool(use_classic)

    # Get LutProvider class (reuse if possible)
    L = _get_mc_luts()

    # Apply algorithm
    func = _marching_cubes_lewiner_cy.marching_cubes
    vertices, faces, normals, values = func(volume, level, L, step_size, use_classic)

    if not len(vertices):
        raise RuntimeError('No surface found at the given iso value.')

    # Output in z-y-x order, as is common in skimage
    vertices = np.fliplr(vertices)
    normals = np.fliplr(normals)

    # Finishing touches to output
    faces.shape = -1, 3
    if gradient_direction == 'descent':
        # MC implementation is right-handed, but gradient_direction is left-handed
        faces = np.fliplr(faces)
    elif not gradient_direction == 'ascent':
        raise ValueError("Incorrect input %s in `gradient_direction`, see "
                         "docstring." % (gradient_direction))
    if spacing != (1, 1, 1):
        vertices = vertices * np.r_[spacing]

    if allow_degenerate:
        return vertices, faces, normals, values
    else:
        fun = _marching_cubes_lewiner_cy.remove_degenerate_faces
        return fun(vertices, faces, normals, values)


def _to_array(args):
    shape, text = args
    byts = base64decode(text.encode('utf-8'))
    ar = np.frombuffer(byts, dtype='int8')
    ar.shape = shape
    return ar


EDGETORELATIVEPOSX = np.array([ [0,1],[1,1],[1,0],[0,0], [0,1],[1,1],[1,0],[0,0], [0,0],[1,1],[1,1],[0,0] ], 'int8')
EDGETORELATIVEPOSY = np.array([ [0,0],[0,1],[1,1],[1,0], [0,0],[0,1],[1,1],[1,0], [0,0],[0,0],[1,1],[1,1] ], 'int8')
EDGETORELATIVEPOSZ = np.array([ [0,0],[0,0],[0,0],[0,0], [1,1],[1,1],[1,1],[1,1], [0,1],[0,1],[0,1],[0,1] ], 'int8')


def _get_mc_luts():
    """ Kind of lazy obtaining of the luts.
    """
    if not hasattr(mcluts, 'THE_LUTS'):
        mcluts.THE_LUTS = _marching_cubes_lewiner_cy.LutProvider(
                EDGETORELATIVEPOSX, EDGETORELATIVEPOSY, EDGETORELATIVEPOSZ,
                _to_array(mcluts.CASESCLASSIC), _to_array(mcluts.CASES),
                _to_array(mcluts.TILING1), _to_array(mcluts.TILING2), _to_array(mcluts.TILING3_1), _to_array(mcluts.TILING3_2),
                _to_array(mcluts.TILING4_1), _to_array(mcluts.TILING4_2), _to_array(mcluts.TILING5), _to_array(mcluts.TILING6_1_1),
                _to_array(mcluts.TILING6_1_2), _to_array(mcluts.TILING6_2), _to_array(mcluts.TILING7_1),
                _to_array(mcluts.TILING7_2), _to_array(mcluts.TILING7_3), _to_array(mcluts.TILING7_4_1),
                _to_array(mcluts.TILING7_4_2), _to_array(mcluts.TILING8), _to_array(mcluts.TILING9),
                _to_array(mcluts.TILING10_1_1), _to_array(mcluts.TILING10_1_1_), _to_array(mcluts.TILING10_1_2),
                _to_array(mcluts.TILING10_2), _to_array(mcluts.TILING10_2_), _to_array(mcluts.TILING11),
                _to_array(mcluts.TILING12_1_1), _to_array(mcluts.TILING12_1_1_), _to_array(mcluts.TILING12_1_2),
                _to_array(mcluts.TILING12_2), _to_array(mcluts.TILING12_2_), _to_array(mcluts.TILING13_1),
                _to_array(mcluts.TILING13_1_), _to_array(mcluts.TILING13_2), _to_array(mcluts.TILING13_2_),
                _to_array(mcluts.TILING13_3), _to_array(mcluts.TILING13_3_), _to_array(mcluts.TILING13_4),
                _to_array(mcluts.TILING13_5_1), _to_array(mcluts.TILING13_5_2), _to_array(mcluts.TILING14),
                _to_array(mcluts.TEST3), _to_array(mcluts.TEST4), _to_array(mcluts.TEST6),
                _to_array(mcluts.TEST7), _to_array(mcluts.TEST10), _to_array(mcluts.TEST12),
                _to_array(mcluts.TEST13), _to_array(mcluts.SUBCONFIG13),
                )

    return mcluts.THE_LUTS

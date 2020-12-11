'''
    head meshing in SPM/CAT12
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2019  Jesper D Nielsen, Axel Thielscher, Guilherme Saturnino

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

'''


# DEPENDENCIES
# =============================================================================
# Python
from glob import glob
from itertools import combinations
import logging
import nibabel as nib
import numpy as np
import os
import scipy.ndimage.morphology as mrph
import shutil
from subprocess import call
import sys
import time
import traceback

# simnibs
from simnibs import SIMNIBSDIR, file_finder
from simnibs.msh import mesh_io
import simnibs.msh.hmutils as hmu
import simnibs.msh.transformations as transformations
from simnibs import __version__

# logging
logging.captureWarnings(True)
logger = logging.getLogger("py.warnings")
logger.setLevel(logging.DEBUG)


# ============================================================================
def headmodel(argv):    
    """
    
    """
    args = parse_args(argv[1:]) # get command line arguments
    # set paths and if they do not exist, create them
    pwd = os.getcwd()
    subject_dir = os.path.join(pwd, "m2m_"+args.subject_id)
    if not os.path.exists(subject_dir):
        os.mkdir(subject_dir)
    
    # log to console and file
    logfile = os.path.join(subject_dir, "headreco_log.html")
    logger.addHandler(make_handler())
    logger.addHandler(make_handler("file", logging.DEBUG, logfile ))
    with open(logfile, 'a') as f:
        f.write('<HTML><HEAD><TITLE>headreco report</TITLE></HEAD><BODY><pre>')
    width = 72
    
    # recreate command line call and log to file 
    hmu.log(" ".join(argv), level=logging.DEBUG)
    
    try:
        docat = args.cat
    except AttributeError:
        docat = True
    if args.mode == "preparecat":
        docat = True
    if docat and " " in SIMNIBSDIR:
        raise IOError(
            f'Unable to run CAT12: Found a space in the path to the SimNIBS installation: {SIMNIBSDIR} '
            f'To fix it, please reinstall SimNIBS on a path without spaces or run with headreco with --no-cat')
    # subdirectories in m2m_{subid} folder
    segment   = os.path.join(subject_dir,"segment")
    spm       = os.path.join(segment,"spm")
    spmbin    = os.path.join(spm,"binarized")
    mask_prep = os.path.join(subject_dir,"mask_prep")
    tmp       = os.path.join(subject_dir,"tmp")
    toMNI       = os.path.join(subject_dir,"toMNI")
    eeg_positions = os.path.join(subject_dir,"eeg_positions")
    templates_dir = os.path.join(subject_dir,"templates")
    if docat:
        cat    = os.path.join(segment,"cat")
        catbin = os.path.join(cat,"binarized")    
    
    # subfolders in path_simnibs containing the extended TPMs, templates and .m file(s)
    etpm_dir    = os.path.join(SIMNIBSDIR,"resources","templates")
    tspm = {}.fromkeys([os.path.basename(f) for f in glob(os.path.join(etpm_dir,"spm*.nii"))])
    
    if docat:
        tcat = {}.fromkeys([os.path.basename(f) for f in glob(os.path.join(etpm_dir,"cat*.nii"))])
    path2mfiles = [os.path.dirname(__file__),os.path.join(SIMNIBSDIR,"resources","spm12")]
    
    # time the entire meshing process
    if args.mode == "all":
        start = time.time()
    
    if args.mode in ["all","preparevols"]:
        """ -----------------------------------------------------
         this step runs SPM12 and optionally CAT12 to segment
         the input images. The detailed steps are:
           1) The input images are reoriented to LAS and
              optionally resampled. Exception: setting option
              "-d no-conform" will keep the T1 unchanged
           2) SPM12 is called to coregister the T2 to the T1,
              generate the transformations to MNI space and
              segment the images; this uses the extended tissue
              priors that cover also the neck. In addition,
              a number of MNI priors for eyes etc are warped
              to subject space for use in the later generation
              of the tissue masks
           3) CAT12 is called (optional, but recommended) to
              segment the cortex; In addition, a number of MNI
              priors are warped to subject space for use in the
              later joining of surfaces of the left and right
              hemispheres
           4) the posterior probability masks are binarized

         result files and directories:
             T1fs_conform, T2_conform: reoriented and resampled
                    input images. The T2 is coregistered to the 
                    T1.
             segment/spm/binarized: binarized posterior 
                    probability masks of the SPM segmentation
             segment/cat/surf: cortical surface reconstructions
             segment/cat/binarized: binarized posterior 
                    probability masks of GM, WM and CSF
                    of CAT segmentation
             templates: spmprior_* are created during SPM12
                    segmentation process, and are later on
                    used to, e.g. separate CSF from eyes
                        cattemplate_* are created during CAT12
                    segmentation process, and are later on
                        used for building a mid-brain-brainstem-
                        cerebellum mask that is joined with
                        the left and right cortical GM surfaces

        ----------------------------------------------------- """

        print("")
        hmu.log("="*width)
        hmu.log("PREPARING VOLUMES FOR MESHING")
        hmu.log("="*width)
        print("")
        section_start = time.time()
        
        # some files to be created
        t1conform0   = os.path.join(subject_dir,"T1fs_conform.nii.gz")
        t1conform1   = os.path.join(pwd, args.subject_id+"_T1fs_conform.nii.gz")
        t1conform_nu = os.path.join(subject_dir, "T1fs_nu_conform.nii.gz")
        mni2c6       = os.path.join(toMNI, "MNI2conform_6DOF.txt")
        mni2c12      = os.path.join(toMNI, "MNI2conform_12DOF.txt")
        
        for k in tspm.keys():
            tspm[k] = os.path.join(spm,"w"+k)
        if docat:
            for k in tcat.keys():
                tcat[k] = os.path.join(cat,"w"+k)
        
        # Rename and move image(s) to subject directory and load them
        hmu.log("Copying input files to {0}", os.path.basename(subject_dir))
        print("")
        fnames = ["T1fs","T2"]
        vols = []
        for f in args.input_files:
            fname = fnames[args.input_files.index(f)]
            dest  = os.path.join(subject_dir, fname+".nii.gz")
            nib.save(nib.load(f),dest) # move/copy and make sure it is gzipped
            vols.append(nib.load(dest))
    
        # If required, reorient to LAS and resample (T1 weighted image only)
        # (T2 image will be updated by SPM call)        
        vals = {
            "out_format":args.dimension,
            "vx_size":args.vx_size
            }
        if args.skip_coreg and len(vols)>1:
            vols_conform = hmu.prepare_vols(vols, **vals)
        else:
            vols_conform = hmu.prepare_vols(vols[0], **vals)

        # Prepare for segmentation and call SPM to segment
        hmu.log("Preparing to segment")
        print("")
        
        # make folders if they do not exist
        if not os.path.exists(segment):
            os.mkdir(segment)
        if not os.path.exists(spm):
            os.mkdir(spm)
        if not os.path.exists(spmbin):
            os.mkdir(spmbin)
        if docat and not os.path.exists(cat):
            os.mkdir(cat)
        if docat and not os.path.exists(catbin):
            os.mkdir(catbin)
        if not os.path.exists(toMNI):
            os.mkdir(toMNI)
        if not os.path.exists(templates_dir):
            os.mkdir(templates_dir)
        
        # Save non-zipped version of conformed T1
        img2sgm = []
        fspm = os.path.join(spm, os.path.basename(vols_conform[0].get_filename()).split(".")[0])+".nii"
        hmu.write_nifti(vols_conform[0].get_data(), fspm, vols_conform[0])
        img2sgm.append(fspm)
        if len(vols) > 1:
            if args.skip_coreg:
                # Save non-zipped version of T2
                fspm = os.path.join(spm, os.path.basename(vols_conform[1].get_filename()).split(".")[0])+".nii"
                hmu.write_nifti(vols_conform[1].get_data(), fspm, vols_conform[1])
                img2sgm.append(fspm)
            else:
                fspm = os.path.join(spm, os.path.basename(vols[1].get_filename()).split(".")[0])+".nii"
                hmu.write_nifti(vols[1].get_data(), fspm, vols[1])
                img2sgm.append(fspm)
        if docat:
            # save T1
            fcat = os.path.join(cat, os.path.basename(vols_conform[0].get_filename()).split(".")[0])+".nii"
            hmu.write_nifti(vols_conform[0].get_data(), fcat, vols_conform[0])
                    
        # Prepare call to MATLAB
        if len(img2sgm)==2:
            spmcall = "segment_SPM('{0}','{1}','{2}',{3},{4},{5})".format(etpm_dir, img2sgm[0], img2sgm[1], args.biasreg, args.dsf, int(args.skip_coreg))
        elif len(img2sgm)==1:
            spmcall = "segment_SPM('{0}','{1}',[],{2},{3})".format(etpm_dir, img2sgm[0], args.biasreg, args.dsf)
        if docat:
            catcall = "segment_CAT('{0}','{1}',{2})".format(fcat,etpm_dir,int(not(args.cat_print)))
            spmcall +="; "+ catcall
        
        # -nosplash -nodesktop - CAT12 will stall on using -nodisplay due to the graphical report
        if (sys.platform == 'win32'):
            cmd = "matlab -nosplash -nodesktop -wait -logfile \""+segment+"\spmcat.log\" -r "
        else:
            cmd = "matlab -nosplash -nodesktop -r "
        cmd += "\"addpath('{0}','{1}');".format(path2mfiles[0],path2mfiles[1])
        cmd += "try,{};catch ME,rethrow(ME);end,exit;\"".format(spmcall)
        
        # Call MATLAB
        hmu.log("Starting MATLAB")
        hmu.log("="*width)
        exitcode = hmu.spawn_process(cmd, verbose=True, return_exit_status=True)
        
        if exitcode is 2:
            raise RuntimeError("Segmentation using CAT12 failed. Could not find CAT12.")
            
        # Binarize posterior probability maps
        hmu.log("Binarizing tissue probability maps")
        # SPM
        fspm = sorted(glob(os.path.join(spm, "c*T1fs_conform*")),
                           key=str.lower)
        if len(fspm) == 0:
            raise OSError('Could not find tissue probability maps.\n'
                           'Probably something went wrong with the Matlab call or SPM')
        masks = hmu.binarize([nib.load(v).get_data() for v in fspm])
        
        M = nib.load(fspm[0])
        for i in range(len(masks)):
            hmu.write_nifti(masks[i], os.path.join(spmbin,
                            os.path.basename(fspm[i])[:-4]+"_bin"), M,
                            dtype=np.uint8)
        # CAT
        if docat:
            fcat = sorted(glob(os.path.join(cat,"mri", "p*.nii")), key=str.lower)
            masks = hmu.binarize([nib.load(v).get_data() for v in fcat])
            M = nib.load(fcat[0])
            for i in range(len(masks)):
                hmu.write_nifti(masks[i], os.path.join(catbin,
                                os.path.basename(fcat[i])[:-4]+"_bin"), M,
                                dtype=np.uint8)
            
        hmu.log("Creating unity.xfm and ref_FS.nii.gz")        
        hmu.write_unity_xfm(subject_dir)
        hmu.write_ref_fs(vols_conform[0])

        # move MNI-to-subject-space coregistrations to toMNI-subfolder
        hmu.log("Moving forward deformation field")
        f = glob(os.path.join(spm,"y*.nii"))
        if len(f)>1:
              raise IOError("more than one forward deformation field found ... I do not know which one to move")
        shutil.move(f[0], os.path.join(toMNI, 'MNI2Conform_nonl.nii'))
        # note: naming of warp fields according to warping of coordinates --> use opposite field for warping images

        f = glob(os.path.join(spm,"iy*.nii"))
        if len(f)>1:
              raise IOError("more than one reverse deformation field found ... I do not know which one to move")
        shutil.move(f[0], os.path.join(toMNI, 'Conform2MNI_nonl.nii'))

        hmu.log("Moving affine MNI to subject space transformations")
        shutil.move(os.path.join(spm,"MNI2conform_6DOF.txt"), mni2c6) 
        shutil.move(os.path.join(spm,"MNI2conform_12DOF.txt"), mni2c12)
        
        # copying "conformed" T1 and T2 images to m2m_{subid} folder
        hmu.log("Copying T1 image")
        shutil.copy(t1conform0, t1conform1)
        
        hmu.log("Copying bias corrected T1 image")
        BF_T1 = nib.load(os.path.join(spm, "mT1fs_conform.nii"))
        hmu.write_nifti(BF_T1.get_data(), t1conform_nu, BF_T1) 
        
        if len(img2sgm)==2:
            if args.skip_coreg:
                fname = os.path.join(spm, "mT2_conform.nii")
                dest  = os.path.join(subject_dir, "T2_conform.nii.gz")
                nib.save(nib.load(fname),dest)
            else:
                fname = os.path.join(spm, "rT2.nii")
                dest  = os.path.join(subject_dir, "T2_conform.nii.gz")
                nib.save(nib.load(fname),dest)

        # moving priors and templates (already warping into subject space) to templates subfolder
        hmu.log("Moving and compressing warped priors and templates")
        for k,src in tspm.items():
            nib.save(nib.load(src), os.path.join(templates_dir, k+".gz"))
        remove_files(tspm.values())
        if docat:
            for k,src in tcat.items():
                nib.save(nib.load(src), os.path.join(templates_dir, k+".gz"))
            remove_files(tcat.values())
        
            dartel = glob(os.path.join(cat,"wTemplate*"))[0]
            nib.save(nib.load(dartel),
                     os.path.join(templates_dir, "cattemplate_dartel.nii.gz"))
            remove_files(dartel)
        
        print("")        
        hmu.log("="*width)
        hmu.log("{:{}s}{:>10}",("Preparation time:", width-10,
                    hmu.format_time(time.time()-section_start)))
        print("")

    
    if args.mode in ["all","preparecat"] and docat:
        """ -----------------------------------------------------
         this step
         1) expands the central surfaces created by
            CAT12 to get the pial surfaces
         2) creates WM (including brainstem...), ventricles 
            and combined corpus-callosum-brainstem-cerebellum
            surfaces from the voxel segmentations
         3) joins the pial surfaces using the combined
            corpus-callosum-brainstem-cerebellum surface
         4) decouples the GM, WM and ventricle surfaces
         5) creates volume masks for those surfaces
        
         final files stored in mask_prep:
         cat_surf_gm.off, cat_surf_wm.off, cat_surf_ventricles1,2,...off
         MASK_CAT_GM.nii.gz, MASK_CAT_WM.nii.gz, MASK_CAT_VENTRICLES.nii.gz

         note: in order to convert a .off file in a .stl file for 
               visual inspection in gmsh, run:
         meshfix cat_surf_gm.off --no-clean --stl -o cat_surf_gm.stl
        
        -----------------------------------------------------"""

        print("")
        hmu.log("="*width)
        hmu.log("PREPARING SURFACES FROM CAT12")
        hmu.log("="*width)
        print("")
        section_start = time.time()
        
        if not os.path.exists(tmp):
            os.mkdir(tmp)

        if not os.path.exists(mask_prep):
            os.mkdir(mask_prep)

        vd = args.vertex_density
        
        hmu.log("Loading files...")
        print("")
        # -------
        # segmentations
        dcat, dspm = {}, {}
        dcat["wm"]  = nib.load(glob(os.path.join(catbin, "p2*"))[0])
        dcat["gm"]  = nib.load(glob(os.path.join(catbin, "p1*"))[0])
        dcat["csf"] = nib.load(glob(os.path.join(catbin, "p3*"))[0])
        dspm["wm"]  = nib.load(glob(os.path.join(spmbin, "c2*"))[0])
        dspm["gm"]  = nib.load(glob(os.path.join(spmbin, "c1*"))[0])
    
        ref = dcat["wm"]
        
        # templates
        vols = glob(os.path.join(templates_dir, "cat*.nii.gz"))
        keys = ["cerebellum", "cc", "fornix", "thalamus", "brainstem",
                "ventricles_lateral", "ventricles_third", "dartel"]
        template = {}.fromkeys(keys)
        for k in template.keys():
            template[k] = [nib.load(v) for v in vols if k in os.path.basename(v).lower()][0]
        template["spinal"] = nib.load(glob(os.path.join(templates_dir, "*spinal*"))[0])
        
        # CAT surfaces and thickness estimates
        lh = glob(os.path.join(cat, "surf", "lh.central.*.gii"))[0]
        rh = glob(os.path.join(cat, "surf", "rh.central.*.gii"))[0]
        lht = glob(os.path.join(cat, "surf", "lh.thickness.*"))[0]
        rht = glob(os.path.join(cat, "surf", "rh.thickness.*"))[0]
        # ------
        
        # Expand CAT surfaces by 0.5 thickness
        gms, gmv = {}, {}
        
        hmu.log("Expanding left hemisphere")
        gms["left"] = expand_cat_surface(lh, lht, "left.off", vertex_density=vd)
        hmu.log("Expanding right hemisphere")
        gms["right"] = expand_cat_surface(rh, rht, "right.off", vertex_density=vd)
        
        # Voxelize the expanded surfaces
        hmu.log("Voxelizing left hemisphere")
        f = os.path.join(tmp, os.path.splitext(os.path.basename(gms["left"]))[0]+".nii.gz")
        gmv["left"] = hmu.make_volume_mask(gms["left"], f, ref).get_data()
        hmu.log("Voxelizing right hemisphere")
        f = os.path.join(tmp, os.path.splitext(os.path.basename(gms["right"]))[0]+".nii.gz")
        gmv["right"]= hmu.make_volume_mask(gms["right"], f, ref).get_data()
        print("")
        
        # fs : surface file; fv : volume file
        hmu.log("Preparing to join hemispheres and extract ventricles")
        fvjoin, fvwm, fvven = make_cat_masks(dcat, dspm, template, gmv, ref, tmp)
        #fsjoin = hmu.make_surface_mesh(fvjoin, fvjoin.split(".")[0], 1, 1, vd, 1)[0]
        # AT 02.04.18: changed to avoid spikes and loops in the joining surface (messes up later decoupling from csf):
        fsjoin = hmu.make_surface_mesh(fvjoin, os.path.splitext(fvjoin)[0], 1, 1, vd, 5, erode_and_expand=True)[0]

        hmu.log("Joining hemispheres")
        fsgm = join_hemispheres(gms, fsjoin, tmp)
        print("")
     
        hmu.log("Voxelizing gray matter")
        fvgm = voxelize(fsgm, ref)[0]
        
        hmu.log("Decoupling volumes")
        vwm = hmu.decouple_volumes(fvven, fvwm, "neighbors")
        vwm = hmu.decouple_volumes(vwm, fvgm, "inner-from-outer")
        nib.save(vwm, vwm.get_filename())
        vven = hmu.decouple_volumes(fvven, fvgm, "inner-from-outer")
        
        # vven=vven&~vwm
        nib.save(vven, vven.get_filename())
         
        hmu.log("Surface extraction: white matter")
        fswm = make_surface_and_shrink(vwm.get_filename(), 1, 1, vd, 5, 0.3)[0]
        hmu.log("Surface extraction: ventricle")
        fsven = make_surface_and_shrink(vven.get_filename(), 5, 200, vd, 5, 0.5)
        
        # decouple WM, GM, and ventricles
        hmu.log("Decoupling surfaces")
        decouple_wm_gm_ventricles(fswm, fsgm, fsven)
        print("")
        
        # save GM, WM and ventricle surfaces in mask_prep subdirectory for later use; AT 01.04.18
        shutil.move(fswm, os.path.join(mask_prep, "cat_surf_wm.off"))
        shutil.move(fsgm, os.path.join(mask_prep, "cat_surf_gm.off"))
        fswm = os.path.join(mask_prep, "cat_surf_wm.off")
        fsgm = os.path.join(mask_prep, "cat_surf_gm.off")
        fsven2=[]
        for f in fsven:
            shutil.move(f, os.path.join(mask_prep, "cat_surf_" + os.path.basename(f)))
            fsven2.append(os.path.join(mask_prep, "cat_surf_" + os.path.basename(f)))
        fsven=fsven2

        # Voxelize WM and ventricles
        hmu.log("Voxelizing white matter")
        fvwm = voxelize(fswm, ref)[0]
        
        hmu.log("Voxelizing ventricles")
        #fvven = voxelize(fsven, ref, merge=True)[0]
        se = mrph.generate_binary_structure(3,1)
        cWM_ERO = nib.load(fvwm)
        cWM_ERO = mrph.binary_erosion(cWM_ERO.get_data(), se, 1)
        
        fsven_to_keep = []
        for i in range(len(fsven)):
            fvven_tmp = voxelize(fsven[i], ref)[0]
            cVEN_tmp = nib.load(fvven_tmp)
            
            if np.any(cWM_ERO & (cVEN_tmp.get_data() > 0) ):
                hmu.log("Throwing out " + fsven[i] + " because it is in WM")
                os.remove(fsven[i])
            else:
                fsven_to_keep.append(fsven[i])
                
            os.remove(fvven_tmp)
        
        fsven=fsven_to_keep
        fvven = voxelize(fsven, ref, merge=True)[0]
        print("")
         
        # move
        shutil.move(fvwm, os.path.join(mask_prep, "MASK_CAT_WM.nii.gz"))
        shutil.move(fvgm, os.path.join(mask_prep, "MASK_CAT_GM.nii.gz"))
        shutil.move(fvven, os.path.join(mask_prep, "MASK_CAT_VENTRICLES.nii.gz"))
        
        hmu.log("="*width)
        hmu.log("{:{}s}{:>10}",("CAT preparation time:", width-10,
                                hmu.format_time(time.time()-section_start))) 
        print("")

        
    if args.mode in ["all","cleanvols"]:
        """ -----------------------------------------------------
         this step cleans the binarized tissue masks
         in segment/spm/binarized and creates masks which contain
         all inner tissue types as well ("volume decoupling")
         
         when GM, WM and ventricle masks were created from the CAT
         results in the step before, they are used for volume
         decoupling the rest of the masks. They themselves are 
         not changed.
        
         final files are stored in mask_prep:
         MASK_*.nii.gz
         
        -----------------------------------------------------""" 

        print("")
        hmu.log("="*width)
        hmu.log("CLEANING SEGMENTATIONS")
        hmu.log("="*width)
        print("")
        
        section_start = time.time()
        
        if not os.path.exists(mask_prep):
            os.mkdir(mask_prep)
        
        # get priors/templates (already coregistered to individual T1)
        for k in tspm.keys():
            tspm[k] = os.path.join(templates_dir, k+".gz")
        if docat:
            for k in tcat.keys():
                tcat[k] = os.path.join(templates_dir, k+".gz")
        
        # binary tissue masks
        sgm_files = sorted(glob(os.path.join(spmbin, "c*bin.nii.gz")),
                           key=str.lower)
        # input kwargs
        vals = {
            "wm_raw"    : sgm_files[1],
            "gm_raw"    : sgm_files[0],
            "csf_raw"   : sgm_files[2],
            "bone_raw"  : sgm_files[3],
            "skin_raw"  : sgm_files[4],
            "air_raw"   : sgm_files[5],
            "out_dir"   : mask_prep,
            "prior_tissue" : [v for k,v in tspm.items() if "tissue" in k][0],
            "prior_eye1" : [v for k,v in tspm.items() if "eye1" in k][0],
            "prior_eye2" : [v for k,v in tspm.items() if "eye2" in k][0],
            "prior_air" : [v for k,v in tspm.items() if "air" in k][0],
            "prior_spinal" : [v for k,v in tspm.items() if "spinal" in k][0],
            "prior_ventricles_lateral" : [v for k,v in tspm.items() if "ventricles" in k][0]
            }
        
        if docat: 
             vals["gm_raw"]= os.path.join(mask_prep, "MASK_CAT_GM.nii.gz")
             vals["wm_raw"]= os.path.join(mask_prep, "MASK_CAT_WM.nii.gz")
             vals["prior_ventricles_lateral"]= os.path.join(mask_prep, "MASK_CAT_VENTRICLES.nii.gz")
             vals["usecat"]=True

        # save cleaned tissue masks and control image to mask_prep
        hmu.vol_mask(**vals)
        
        # rename control volume
        control_mask = os.path.join(subject_dir,
                                    args.subject_id+"_masks_contr.nii.gz")       
        shutil.move(os.path.join(mask_prep,"control.nii.gz"), control_mask)
                    
        hmu.log("="*width)
        hmu.log("{:{}s}{:>10}",("Cleaning time:", width-10,
                hmu.format_time(time.time()-section_start))) 
        print("")

        
    if args.mode in ["all","surfacemesh"]:
        """ -----------------------------------------------------
         this step creates the surface .stl files ready for 
         volume-meshing. The surfaces are build up in a 
         inner-to-outer order, i.e. the next-outer surface
         is decoupled from the next-inner one.
        
         The following steps are applied for each surface:
         1) decouple its voxel mask from the next-inner
         voxel mask created from final surface mesh (see step 4)
         2) build a surface from the updated voxel mask
         3) decouple surface from next-inner surface
         4) update its voxel mask using the decoupled surface
        
         In case surface decoupling of bone from csf fails,
         a decoupling of csf from bone is tried (i.e. 
         the order is to reversed to outer-to-inner, including
         checks for gm and wm). This was included, as "spikes"
         in csf can sometimes prevent a successful decoupling
         of bone from csf.
         
         When CAT12 was used, then the GM, WM and ventricle
         surfaces and their volume masks are only copied
         from the mask_prep folder, without any additional
         decoupling steps

         The final files (stl surfaces and corresponding
         nifti voxel masks) are stored in the m2m_{subid}
         folder
         
        -----------------------------------------------------"""

        print("")
        hmu.log("="*width)
        hmu.log("GENERATING SURFACE MESHES")
        hmu.log("="*width)
        print("")
        
        section_start = time.time()
        
        if not os.path.exists(tmp):
            os.mkdir(tmp)
        
        masks = sorted(glob(os.path.join(mask_prep, "MASK_*.nii.gz")),
                                 key=str.lower)
      
        surfaces, relabel = make_surface_meshes(masks,
            subject_dir, tmp, args.vertex_density, docat)
        
        hmu.log("="*width)
        hmu.log("{:{}s}{:>10}",("Surface meshing time:", width-10,
                    hmu.format_time(time.time()-section_start)))
        print("")

        
    if args.mode in ["all","volumemesh"]:
        """ -----------------------------------------------------
         this step creates the final volume mesh from the 
         surface .stl files
        
         The following steps are applied for each surface:
         1) a {subid}.geo file is created, and gmsh is called for
            meshing
         2) after meshing, tetrahedra corresponding to air are
            deleted and the tissue labels are updated
            to conform with the conventions used for the
            FEM calculations; in addition, triangle and tetrahedra
            orientations are checked and, if needed, updated; also
            very thin tetrahedra are made slightly thicker to improve
            the mesh quality
         3) a nifti-label-image is created from the final mesh; 
            in addition, final wm and gm volume masks are extracted
            from the mesh and converted (together with the T1) to MNI,
            to allow for visual control segmentation and
            co-registration accuracy

         The final mesh is stored in the same directory containing
         the m2m_{subid} folder
         
        -----------------------------------------------------"""

        print("")
        hmu.log("="*width)
        hmu.log("GENERATING VOLUME MESH")
        hmu.log("="*width)
        print("")
        
        section_start = time.time()


        try:
            relabel
        except NameError:
            relabel = sorted(glob(os.path.join(subject_dir,
                            "air_outer_relabel*.stl")), key=str.lower)
        try:
            surfaces
        except NameError:
            surfaces = sorted(list(set(
                glob(os.path.join(subject_dir,"*.stl")))-set(relabel)),
                key=str.lower)
        
        final_mesh = os.path.join(pwd, args.subject_id+".msh")
        final_control = os.path.join(subject_dir, args.subject_id+"_final_contr.nii.gz")    
        
        T1 = nib.load(os.path.join(subject_dir,"T1fs_conform.nii.gz"))
        
        # run Gmsh
        vmsh = make_volume_mesh(args.subject_id, surfaces, subject_dir, args.keep_air)
        ovmsh = vmsh
        
        # load volume mesh
        # using buffered reading due to ineffcient version 2 format written by gmsh
        vmsh = mesh_io.read_msh(vmsh, buffered=True)
        
        # if ventricles exist, relabel to CSF
        if any(["ventricles" in s for s in surfaces]):
            hmu.log("Relabelling ventricles to CSF")
            vmsh = hmu.relabel_compartments(vmsh, [9,1009], [3,1003])
            
        if len(relabel) > 0:
            hmu.log("Adjusting volume mesh")
            vmsh = hmu.volume_mesh_relabel(vmsh, relabel, 5)
            
        vmsh.fn = final_mesh
        hmu.log("Fix surface labels")
        vmsh.fix_surface_labels()
        hmu.log("Fix thin tetrahedra")
        vmsh.fix_thin_tetrahedra()
        hmu.log("Fix tetrahedra node ordering")
        vmsh.fix_th_node_ordering()
        hmu.log("Fix triangle node ordering")
        vmsh.fix_tr_node_ordering()
        hmu.log("Write final mesh")
        mesh_io.write_msh(vmsh, final_mesh)
        remove_files(ovmsh)

        hmu.log("Creating control volume and GM and WM masks")
        hmu.vmesh2nifti(vmsh, T1, final_control)

        hmu.log("Checking if the meshing went OK")
        if not hmu.check_volumes_meshed(subject_dir):
            os.remove(final_mesh)
            print("")
            raise RuntimeError(
                "Gmsh failed meshing one or more volumes. "
                "Please try running headreco again "
                "increasing the '-v' option "
                "(eg. -v 1.0). "
                "Type 'headreco all -h' for more informaion")

        vmsh = vmsh.crop_mesh(elm_type=4)
        mask_names = [os.path.join(subject_dir, "wm_fromMesh.nii.gz"),
                      os.path.join(subject_dir, "gm_fromMesh.nii.gz")]
        affine = T1.affine
        n_voxels = T1.header['dim'][1:4]
        for k in [1, 2]:
           field = np.zeros(vmsh.elm.nr, dtype=np.uint8)
           field[vmsh.elm.tag1 == k] = 1
           ed = mesh_io.ElementData(field)
           ed.mesh = vmsh
           ed.to_nifti(n_voxels, affine, fn=mask_names[k-1], qform=T1.header.get_qform(),
                       method='assign')

        fnames_org = [os.path.join(subject_dir,"T1fs_nu_conform.nii.gz"),
                      os.path.join(subject_dir,"wm_fromMesh.nii.gz"),
                      os.path.join(subject_dir,"gm_fromMesh.nii.gz")]
        # "MNI" will be automatically added at the end of the filenames:
        fnames_MNI = [os.path.join(subject_dir,"toMNI","T1fs_nu_nonlin.nii.gz"),
                      os.path.join(subject_dir,"toMNI","wm.nii.gz"),
                      os.path.join(subject_dir,"toMNI","gm.nii.gz")]
        for k in [0, 1, 2]:
           transformations.warp_volume(
                fnames_org[k], subject_dir, fnames_MNI[k],
                transformation_direction='subject2mni',transformation_type='nonl',
                reference=None,mask=None,labels=None,order=3, binary=k>0)

        transformations.warp_volume(
                fnames_org[0], subject_dir, os.path.join(subject_dir,"toMNI","T1fs_nu_12DOF.nii.gz"), # "MNI" will be automatically added
                transformation_direction='subject2mni',transformation_type='12dof',
                reference=None,mask=None,labels=None,order=3)

        hmu.log("Transforming EEG positions")
        if not os.path.exists(eeg_positions):
            os.mkdir(eeg_positions)
        cap_file = os.path.abspath(os.path.realpath(
                        os.path.join(
                            SIMNIBSDIR,
                            'resources', 'ElectrodeCaps_MNI',
                            'EEG10-10_UI_Jurak_2007.csv')))
        cap_out = os.path.join(eeg_positions, 'EEG10-10_UI_Jurak_2007.csv')
        geo_out = os.path.join(eeg_positions, 'EEG10-10_UI_Jurak_2007.geo')
        transformations.warp_coordinates(
            cap_file, subject_dir,
            transformation_direction='mni2subject',
            out_name=cap_out,
            out_geo=geo_out)

        hmu.log("="*width)
        hmu.log("{:{}s}{:>10}", ("Volume meshing time:", width-10,
                    hmu.format_time(time.time()-section_start)))
        print("")



        
    if args.mode == "check":
        """ -----------------------------------------------------
           this step shows some results for visual control
           of the segmentation accuracy and accuracy of the
           subject-to-MNI transformation

           if installed, freeview (from FreeSurfer package)
           is used; this is the recommended use;
           otherwise, some basic results are shown using
           the SPM viewer
        ----------------------------------------------------- """


        relabel = sorted(glob(os.path.join(subject_dir,
                                   "air_outer_relabel*.stl")),key=str.lower)
        surfaces = sorted(list(set(
                glob(os.path.join(subject_dir,"*.stl")))-set(relabel)),
                key=str.lower)
                
        T1 = os.path.join(subject_dir, "T1fs_conform.nii.gz")
        mask_initial = os.path.join(subject_dir,
                                   args.subject_id+"_masks_contr.nii.gz")
        mask_final = os.path.join(subject_dir,
                                   args.subject_id+"_final_contr.nii.gz")

        # test if freeview is installed
        usefreeview=True
        try:
            # suppress output by piping it to os.devnull
            call("freeview -h".split(), stdout=open(os.devnull, 'wb'))
        except OSError:
            usefreeview=False

        if usefreeview:                 
             visualize(args.subject_id, surfaces, T1, mask_initial, mask_final,
                       subject_dir, args.force_remake, show_mni=True)
        else:
             visualize_spm(args.subject_id,T1,mask_final,subject_dir,show_mni=True)
             args.noclean = True


    if not args.noclean:
        """ -----------------------------------------------------
           this step cleans intermediate results stored 
           in the tmp subdirectory
        ----------------------------------------------------- """
        hmu.log("Deleting temporary files")
        remove_dirs(tmp)            
        hmu.log("Done")


    if args.mode == "all":
        """ -----------------------------------------------------
           after running all steps of the pipeline, log
           required time
        ----------------------------------------------------- """

        hmu.log("="*width)
        hmu.log("{:{}s}{:>10}", ("Total running time:", width-10,
                    hmu.format_time(time.time()-start)))
        hmu.log("="*width)

    with open(logfile, 'a') as f:
        f.write('</pre></BODY></HTML>')
    # close logger handlers
        #for h in logger.handlers:
        #    h.close()


# ===========================
# HELPER FUNCTIONS
# ===========================

def expand_cat_surface(central, thickness, outname, vertex_density=0.5):
    """

    """
    lvl = logging.DEBUG
    
    # meshfix setup    
    meshfix = hmu.path2bin("meshfix")
    resample = '"{0}" "{1}" -a 2.0 -u 5 -q --vertexDensity {2} -o "{1}"'.format(meshfix, "{0}", vertex_density)
    clean    = '"{0}" "{1}" -a 2.0 -q -o "{1}"'.format(meshfix, "{0}")
    uniformremesh = '"{0}" "{1}" -a 2.0 -u 1 -q -o "{1}"'.format(meshfix, "{0}")
    
    outdir = os.path.dirname(central)
    central = nib.load(central)

    # load surfaces
    cv = central.darrays[0].data
    cf = central.darrays[1].data
    
    # load thickness
    th = hmu.read_curv(thickness)
    th /= 2
    
    hmu.log("Expanding", level=lvl)
    v1 = hmu.deform_mesh(cv,cf,th, ensure_distance=0.5, nsteps=5)
    f = os.path.join(outdir, outname)
    hmu.mesh_save(v1,cf,f)
    
    hmu.log("Cleaning", level=lvl)        
    hmu.spawn_process(clean.format(f))
    hmu.log("Resampling to a vertex density of {0}", vertex_density, level=lvl)
    hmu.spawn_process(resample.format(f))
    hmu.log("Cleaning", level=lvl)        
    hmu.spawn_process(clean.format(f))
    hmu.log("Uniform remeshing", level=lvl)
    hmu.spawn_process(uniformremesh.format(f))
    hmu.log("Cleaning", level=lvl)
    hmu.spawn_process(clean.format(f))

    return f


def make_cat_masks(cat, spm, template, gm, ref, outdir):
    """
    
    cat : wm, gm, csf
    spm : wm, gm
    template : brainstem, cc, cerebellum, dartel
    
    """
    """
    # Setup
    meshfix = hmu.path2bin("meshfix")
    
    # meshfix commands
    clean = meshfix+" {0} -a 2.0 -q -o {0}"
    erode = meshfix+" {0} -a 2.0 --dilate -0.5 -o {0}"
    """
    
    ext = "nii.gz"
    ext = "."+ext
    
    dtype = np.uint8
    se = mrph.generate_binary_structure(3,1)
    
    cat = cat.copy()
    spm = spm.copy()
    template = template.copy()
    
    # Load data
    for k,v in cat.items():
        cat[k] = v.get_data().astype(bool)
    for k,v in spm.items():
        spm[k] = v.get_data().astype(bool)
    for k,v in template.items():
        if k == "dartel":
            template[k] = (v.get_data() > 0.5)[...,1] # choose only WM
        elif k == "spinal":
            template[k] = v.get_data() > 0.25
        else:
            template[k] = v.get_data().astype(bool)
    
    # extend WM mask
    # find minimum z value of CAT template (i.e., the z value below which to
    # insert the SPM segmentation of white matter)
    hmu.log("Joining CAT and SPM white matter segmentations")
    wmx = cat["wm"].copy()
    
    # wm template
    # add WM and GM below template
    wm_add = np.zeros_like(wmx,bool)
    gm_add = np.zeros_like(wmx,bool)
    z = np.where(template["dartel"])[2]
    zcutoff = z.min()
    gm_add[...,:zcutoff] = spm["gm"][...,:zcutoff]
    wm_add[...,:zcutoff] = wmx[...,:zcutoff] | (spm["wm"][...,:zcutoff])
    wm_add = wm_add | gm_add
    
    # spinal template
    # add WM and GM inside spinal prior...
    spinmask = template["spinal"]
    wm_add[spinmask] = wm_add[spinmask] | spm["wm"][spinmask] | spm["gm"][spinmask]
    wm_add[wm_add] = hmu.get_n_largest_components((wmx | wm_add),se,1)[wm_add]
    wm_add = mrph.binary_closing(np.pad(wm_add,((0,0),(0,0),(3,3)),"constant"),se,3)[...,3:-3]    
    wmx = wmx | wm_add

    fwmx = os.path.join(outdir, "wmx"+ext)
    hmu.write_nifti(wmx, fwmx, ref, dtype=dtype)
    

    print("Preparing cerebellum")
    # we want to avoid overlap with cortical gray matter
    dmask = template["cerebellum"]
    cerebellum = np.zeros_like(dmask)
    cerebellum[dmask] = (cat["gm"] | wmx)[dmask]
    cerebellum = mrph.binary_fill_holes(cerebellum,se)
    # decouple
    nd = 2
    cerebellum = (cerebellum.astype(int)-(mrph.binary_dilation(gm["left"],se,nd)+mrph.binary_dilation(gm["right"],se,nd))) > 0 
    fcerebellum = os.path.join(outdir, "cerebellum"+ext)
    hmu.write_nifti(cerebellum, fcerebellum, ref, dtype=dtype)
    
    # make masks for joining hemispheres
    hmu.log("Preparing mask for joining hemispheres")    
    # corpus callosum
    cc = template["cc"]
    cc = mrph.binary_dilation(cc,se,2) & (wmx > 0.5)
    cc = cc | mrph.binary_dilation(cc,se,1) & cat["gm"]
    cc = hmu.get_n_largest_components(cc,se,1)
    fcc = os.path.join(outdir, "join_cc"+ext)
    hmu.write_nifti(cc, fcc, ref, dtype=dtype)
    
    # fornix
    fornix = template["fornix"]
    fornix = mrph.binary_dilation(fornix,se,2) & (wmx > 0.5)
    fornix = hmu.get_n_largest_components(fornix,se,1)
    ffornix = os.path.join(outdir, "join_fornix"+ext)
    hmu.write_nifti(fornix, ffornix, ref, dtype=dtype)
    
    # ventricles
    ventricles = (mrph.binary_dilation(template["ventricles_lateral"],se,5) | template["ventricles_third"]) & cat["csf"]
    fventriclesjoin = os.path.join(outdir, "join_ventricles"+ext)
    hmu.write_nifti(ventricles, fventriclesjoin, ref, dtype=dtype)
    
    # thalamus
    thal = template["thalamus"]
    thalamus = np.zeros_like(thal)
    thalamus[(gm["left"]+gm["right"]) > 0.5] = thal[(gm["left"]+gm["right"]) > 0.5]
    thalamus = hmu.get_n_largest_components(thalamus,se,1)
    fthalamus = os.path.join(outdir, "join_thalamus"+ext)
    hmu.write_nifti(thalamus, fthalamus, ref, dtype=dtype)
    
    # brainstem
    bstm = mrph.binary_dilation(template["brainstem"],se,5) & wmx # 2,3?
    bstm[template["brainstem"]] += cat["gm"][template["brainstem"]]
    bstm[...,:zcutoff] += mrph.binary_closing(wmx[...,:zcutoff],se,3)
    bstm = bstm | cerebellum
    bstm = hmu.get_n_largest_components(bstm,se,1)
    fbrainstem = os.path.join(outdir, "join_brainstem"+ext)
    hmu.write_nifti(bstm, fbrainstem, ref, dtype=dtype)
    
    # combine masks  
    join = mrph.binary_closing((ventricles+cc+fornix+thalamus)>0,se,2)
    join = mrph.binary_closing(join | (bstm > 0),se,1)
    join = mrph.binary_fill_holes(join,se)
    join = hmu.get_n_largest_components(join,se,1)
    fall = os.path.join(outdir, "join_all"+ext)
    vall = hmu.write_nifti(join, fall, ref, dtype=dtype, return_volume=True).get_data()
    
    hmu.log("Preparing ventricles")
    # 1) erode_once(gm | wm | vall) & csfspm > 0.5
    # 2) label components and retain those larger than X (e.g., 100/np.mean(vx_size))
    # 3) this is the ventricles?
    
    # mask to decouple WM from CSF/ventricles and leave space for a little GM?
    ventricles = mrph.binary_dilation(ventricles,se,3) & cat["csf"]# & ~wmx & ~cat["gm"] 
    ventricles = mrph.binary_fill_holes(ventricles,se)
    ventricles = mrph.binary_closing(ventricles,se,1)
    
    # for decoupling with gray matter
    ventricles = ventricles & (vall | gm["right"] | gm["left"])
    ventricles = ventricles & ~wmx
    ventricles = hmu.get_n_largest_components(ventricles,se,1)
    
    se2 = mrph.generate_binary_structure(3,2)
    vx = mrph.binary_erosion(gm["right"] | gm["left"],se,1) & cat["csf"]
    vx = vx & ~(cat["wm"] | cat["gm"])
    vx = mrph.binary_fill_holes(vx,se)
    vx = hmu.get_large_components(vx,se2,100)
    
    ventricles = mrph.binary_fill_holes(ventricles | vx, se)    
    ventricles = hmu.get_large_components(ventricles,se2,100)
    #ventricles = hmu.get_n_largest_components(ventricles,se,1)
    fven = os.path.join(outdir, "ventricles"+ext)
    hmu.write_nifti(ventricles, fven, ref, dtype=dtype)
            
    return fall, fwmx, fven


def make_surface_and_shrink(volume, nshells=1, min_shell_area=1,
                            vertex_density=0.5, nsmooth=1, shrink_amount=0.5):
    """Wrapper to make surface mesh and shrink it by some (small) amount
    """
    meshfix = hmu.path2bin("meshfix")
    clean = '"{0}" "{1}" -a 2.0 -q -o "{1}"'.format(meshfix, "{0}")
    
    surfaces = hmu.make_surface_mesh(volume, os.path.splitext(volume)[0] , nshells,
                                     min_shell_area, vertex_density, nsmooth)
    for s in surfaces:
        v,f = hmu.mesh_load(s)
        v = hmu.deform_mesh(v, f, mm2move=shrink_amount,
                            ensure_distance=0.5, nsteps=2, deform="shrink")
        hmu.mesh_save(v,f,s)
        hmu.spawn_process(clean.format(s))
    
    return surfaces

       
def join_hemispheres(gm, join, outdir):
    """
    
    """
    meshfix = hmu.path2bin("meshfix")
    
    fgmx = os.path.join(outdir, "gmx"+".off")
    # join hemispheres...
    join = '"{0}" "{1}" "{2}" -a 2.0 --shells 2 -j -o "{3}"'.format(meshfix, gm["left"], join, fgmx)
    hmu.spawn_process(join)
    join = '"{0}" "{1}" "{2}" -a 2.0 --shells 2 -j -o "{1}"'.format(meshfix, fgmx, gm["right"])
    hmu.spawn_process(join)
    
    # cleaning to get rid of putatively remaining intersections
    clean = '"{0}" "{1}" -a 2.0 -o "{1}"'.format(meshfix, fgmx)
    hmu.spawn_process(clean)

    return fgmx


def decouple_wm_gm_ventricles(wm, gm, ventricles):
    """Decouple surfaces. Modifies the files in place.
    
    """
    meshfix = hmu.path2bin("meshfix")
    
    # check intersections
    check_intersect = '"{0}" "{1}" "{2}" --shells 2 --no-clean -q --intersect'.format(meshfix, "{0}", "{1}")
    # check if any intersections with the first surface                                
    any_intersect = \
        lambda s: any([hmu.spawn_process(check_intersect.format(s1,s2), return_exit_status=True)==0 for s1,s2 in combinations(s,2)])
    
    if not isinstance(ventricles, list):
        ventricles = [ventricles]
    
    # check inputs
    assert isinstance(wm, str) and os.path.isfile(wm)
    assert isinstance(gm, str) and os.path.isfile(gm)
    for v in ventricles:
        assert isinstance(v, str) and os.path.isfile(v)
    
    veryfirstpass = True
    md = 0.5
    while veryfirstpass or any_intersect(ventricles+[wm,gm]):
        for f in ventricles[::-1]:
            fb = os.path.splitext(os.path.basename(f))[0]
            firstpass = True
            while (firstpass and (md>0)) or any_intersect([f, wm, gm]):
                hmu.log("Decoupling {} from GM".format(fb))
                hmu.decouple_surfaces(f, gm, "inner-from-outer", min_distance=md)
                hmu.log("Decoupling WM from {}".format(fb))
                hmu.decouple_surfaces(wm, f, "neighbor", min_distance=md)
                hmu.log("Decoupling GM from WM")
                hmu.decouple_surfaces(wm, gm, "outer-from-inner", cut_inner=True)
    
                firstpass = False
        veryfirstpass = False
            

def voxelize(files, ref, niiname=None, merge=False):
    """Wrapper function for voxelizing surface files (STL, OFF).
    
    """

    ext = ".nii.gz"
    dtype = np.uint8
    
    # check inputs
    if isinstance(files, str):
        files = [files]
    if len(files) is 1:
        merge = True
        
    if niiname is None:
        if merge:
            niiname = [os.path.splitext(files[0])[0]+ext]*len(files)
        else:
            niiname = [os.path.splitext(files[i]+str(i))[0]+ext for i in range(len(files))]
    elif isinstance(niiname, str):
        niiname = [niiname]

    if isinstance(niiname, list) and len(niiname) is 1:
        niiname = [niiname[0] if merge else hmu.add2filename(niiname[0],str(i)) for i in range(len(files))]
        
    if isinstance(niiname, list) and len(niiname) not in [1,len(files)]:
        raise IOError("Number of items in 'stlname' must be one or match that of 'files'.")
        
    if merge:
        msurfs = np.zeros_like(ref.get_data(), dtype)

    if len(niiname) not in [1,len(files)]:
        raise IOError("Number of items in 'niiname' must be one or match that of 'files'.")
    
    # voxelize
    vols = []
    for i in range(len(files)):
        if not (niiname[i].endswith(ext) or niiname[i].endswith(".nii")):
            niiname[i] = os.path.splitext(niiname[i])[0]+ext
        surf = hmu.make_volume_mask(files[i], niiname[i], ref, write_volume=False)  
        if merge:
            msurfs += surf.get_data()
        else:
            vols.append(niiname[i])
            nib.save(surf, niiname[i])
            
    if merge:
        vols.append(niiname[0])
        hmu.write_nifti(msurfs, niiname[0], ref, dtype=dtype)
    
    return vols


def convert2stl(files, stlname=None, merge=False):
    """Wrapper function for converting surface files to STL.
    
    """
    meshfix = hmu.path2bin("meshfix")
    write_stl = '"{0}" "{1}" --no-clean -q --stl -o "{2}"'.format(meshfix, "{0}", "{1}")
    if merge:
        add2stl = '"{0}" "{1}" "{2}" -jc --shells {3} -q --no-clean --stl -o "{1}"'.format(meshfix, "{0}", "{1}", len(files))
    
    # check inputs
    if isinstance(files, str):
        files = [files]
    if len(files) is 1:
        merge = True
        
    if stlname is None:
        if merge:
            stlname = [os.path.splitext(files[0])[0]+".stl"]*len(files)
        else:
            stlname = [os.path.splitext(files[i]+str(i))[0]+".stl" for i in range(len(files))]
    elif isinstance(stlname, str):
        stlname = [stlname]

    if isinstance(stlname, list) and len(stlname) is 1:
        stlname = [stlname[0] if merge else hmu.add2filename(stlname[0],str(i)) for i in range(len(files))]
        
    if len(stlname) not in [1,len(files)]:
        raise IOError("Number of items in 'stlname' must be one or match that of 'files'.")
    
    # convert
    stls = []
    i = 0
    for i in range(len(files)):
        if (i > 0) and merge:
            hmu.spawn_process(add2stl.format(stlname[i], files[i]))  
        else:
            hmu.spawn_process(write_stl.format(files[i], stlname[i]))
            stls.append(stlname[i])
        i += 1
    
    return stls  


def make_surface_meshes(input_files,out_dir,temp_dir,vertex_density=0.5,
                        catsurf=False):
    """Generates surface renderings from volume masks, processes them (e.g., 
    cleaning and smoothing), reconstructs (fixed) volume masks, and...
    
    
    
    Uses meshfix.
    
    PARAMETERS
    ----------
    input_files             input files or nibabel objects
    vertex_density   number of vertices per mm2 surfaces to keep

    
    RETURNS
    ----------
    sorted list of surface file (.stl) paths
    """ 
    
    if type(input_files)==str:
        input_files=[input_files]
        
    # if input is not iterable, make it        
    try:
        next(iter(input_files))
    except TypeError:
        input_files = [input_files]
    
    vols=[]
    for f in input_files:
        if type(f) == str:
            try:
                vols.append(nib.load(f))
            except:
                raise IOError(f+" does not exist!") 
        elif type(f) in [nib.nifti1.Nifti1Image,nib.nifti2.Nifti2Image,nib.nifti1.Nifti1Pair,nib.nifti2.Nifti2Pair]:
            vols.append(f)
        else:
            raise ValueError("Input volumes must be specified as either a string (path to volume) or a nibabel image object.")
    
    # get the volumes
    mask = {}
    mask = mask.fromkeys(["wm", "gm" , "csf", "eyes", "bone", "skin", "air",
                          "ventricles"])
    mask_ = list(mask.keys())
    for k in mask_:
        try:
            mask[k] = [v for v in vols if k in os.path.basename(v.get_filename()).lower()][0]
        except IndexError:
            del mask[k]
    meshfix = hmu.path2bin("meshfix")
    width = 72

    # if they do not exist, create them
    if not os.path.exists(temp_dir):
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        os.mkdir(temp_dir)
    
    se = mrph.generate_binary_structure(3,1)#np.ones((3,3,3)) # structuring element
    
    # =========================================================================
    # SURFACE GENERATION
    # =========================================================================
    """Start with white matter and work in an inside-to-outside fashion.
    """
    
    # setup input/output filenames
    off = {}
    off = off.fromkeys(mask)
    off.pop("air")
    off["cavities"] = None
    off["air_outer_relabel"] = None
    off["air_outer_decoupled"] = None
    if catsurf:
        off["ventricles"] = None
        
    nii, stl, mod = {}, {}, {}
    for k in off.keys():
        off[k] = os.path.join(temp_dir, "fixed_"+k)
        mod[k] = os.path.join(temp_dir, "mod_"+k+".nii.gz")
        nii[k] = os.path.join(out_dir, k+".nii.gz")
        stl[k] = os.path.join(out_dir, k+".stl")
    
    mdist = 0.5
    
    ref = mask["wm"] # reference volume for header information
    dtype = np.uint8 # data type
    
    # check intersections
    check_intersect = '"{0}" "{1}" "{2}" --shells 2 --no-clean -q --intersect'.format(meshfix, "{0}", "{1}")
    # check if any intersections with the first surface                                
    any_intersect_1st = \
        lambda s: any([hmu.spawn_process(check_intersect.format(s1,s2),return_exit_status=True)==0 for s1,s2 in combinations(s,2) if s[0] in (s1,s2)])
    
    # =========================================================================
    key = "wm"
    
    hmu.log("Processing {}", key.upper())
    hmu.log("="*width)   
    start = time.time()

    if catsurf:
        hmu.log("Fetching previously created files")
        assert os.path.isfile(mask[key].get_filename()), "Could not find {}".format(mask[key].get_filename())
        shutil.copy2(mask[key].get_filename(), nii[key])
        off[key] = off[key] + ".off"
        shutil.copy2(os.path.join(os.path.dirname(mask[key].get_filename()), "cat_surf_wm.off"), off[key])

        cWM = nib.load(nii[key])
    else:
        # options:
        # * n shells to keep
        # * minimum shell area
        # * vertex density
        # * n steps of smoothing
        hmu.log("Extracting surface")
        off[key] = hmu.make_surface_mesh(mask[key].get_filename(), off[key],
                                          1, 1, vertex_density, 5)[0]
        hmu.log("Voxelizing")
        cWM = hmu.make_volume_mask(off[key], nii[key], ref)
           
    print("")    
    hmu.log("{:{}s}{:>10}", ("{} processing time:".format(key.upper()),
                             width-10, hmu.format_time(time.time()-start)))
    print("")

    # =========================================================================
    key = "gm"
    
    hmu.log("Processing {}", key.upper())
    hmu.log("="*width)   
    start = time.time()
    
    if catsurf:
        hmu.log("Fetching previously created files")
        assert os.path.isfile(mask[key].get_filename()), "Could not find {}".format(mask[key].get_filename())
        shutil.copy2(mask[key].get_filename(), nii[key])
        off[key] = off[key] + ".off"
        shutil.copy2(os.path.join(os.path.dirname(mask[key].get_filename()), "cat_surf_gm.off"), off[key])

        cGM = nib.load(nii[key])
    else:       
        hmu.log("Preparing volume for meshing")
        # add wm, binarize, and fill small holes
        mGM_Y_mod = mrph.binary_fill_holes((cWM.get_data()+mask[key].get_data()) > 0)
        hmu.write_nifti(mGM_Y_mod, mod[key], ref, dtype=dtype)
        
        hmu.log("Extracting surface")
        off[key] = hmu.make_surface_mesh(mod[key], off[key], 1, 1,
                       vertex_density, 5, erode_and_expand=True)[0]
        
        hmu.log("Decoupling from WM")
        hmu.decouple_surfaces(off["wm"],off[key], "outer-from-inner",
                              cut_inner=True)
        
        hmu.log("Voxelizing")
        cGM = hmu.make_volume_mask(off[key], nii[key], ref)
    
    print("")    
    hmu.log("{:{}s}{:>10}", ("{} processing time:".format(key.upper()),
                             width-10, hmu.format_time(time.time()-start)))
    print("")

    # =========================================================================
    key = "csf"
    
    hmu.log("Processing {}", key.upper())
    hmu.log("="*width)   
    start = time.time()
    
    # add gm, binarize, and fill small holes
    hmu.log("Preparing volume for meshing")
    mCSF_Y_mod = mrph.binary_fill_holes((cGM.get_data()+mask[key].get_data()) > 0)
    hmu.write_nifti(mCSF_Y_mod, mod[key], ref, dtype=dtype)
    
    hmu.log("Extracting surface")
    off[key] = hmu.make_surface_mesh(mod[key], off[key], 1, 1, vertex_density, 5)[0] 
    
    hmu.log("Decoupling from GM")
    hmu.decouple_surfaces(off["gm"], off[key], "outer-from-inner",
                          min_distance=mdist)
    
    hmu.log("Voxelizing")                           
    cCSF = hmu.make_volume_mask(off[key], nii[key], ref)

    print("")
    hmu.log("{:{}s}{:>10}", ("{} processing time:".format(key.upper()),
                             width-10, hmu.format_time(time.time()-start)))
    print("")
    
    # =========================================================================
    key = "ventricles"
    
    if catsurf or key in mask:

        hmu.log("Processing {}", key.upper())
        hmu.log("="*width)   
        start = time.time()
        
        if catsurf:
            hmu.log("Fetching previously created files")
            assert os.path.isfile(mask[key].get_filename()), "Could not find {}".format(mask[key].get_filename())
            shutil.copy2(mask[key].get_filename(), nii[key])

            fsven=glob(os.path.join(os.path.dirname(mask[key].get_filename()), "cat_surf_ventricles*.off"))
            for i in range(len(fsven)):
                shutil.copy2(fsven[i], off[key] + str(i+1) + ".off")

            off[key] = sorted(glob(os.path.join(temp_dir, off[key] +"*.off")), key=str.lower)
            
        elif key in mask:
            hmu.log("Preparing volume for meshing")        
            modven = mask[key].get_data() & mrph.binary_erosion(cCSF.get_data(), se, 1)
            hmu.write_nifti(modven, mod[key], ref, dtype=dtype)
            
            hmu.log("Extracting surface(s)")
            off[key] = hmu.make_surface_mesh(mod[key],off[key], 9, 100,
                                             vertex_density, 5)
            hmu.log("Decoupling")
            for s1,s2 in combinations(off[key],2):
                firstpass = True
                while firstpass or any_intersect_1st([s1,s2]):
                    hmu.log("Decoupling {0} from {1}".format(
                        os.path.splitext(os.path.basename(s1))[0],
                        os.path.splitext(os.path.basename(s2))[0]))
                    hmu.decouple_surfaces(s1,s2, "neighbor", min_distance=mdist)
                    firstpass = False

            for f in off[key]:
                fb = os.path.splitext(os.path.basename(f))[0]
                firstpass = True
                while firstpass or any_intersect_1st([f,off["gm"],off["csf"],off["wm"]]):
                    hmu.log("Decoupling {} from GM", fb)
                    hmu.decouple_surfaces(f,off["gm"], "inner-from-outer", min_distance=mdist)
                    hmu.log("Decoupling {} from CSF", fb)
                    hmu.decouple_surfaces(f,off["csf"], "inner-from-outer", min_distance=mdist)
                    hmu.log("Decoupling {} from WM", fb)
                    hmu.decouple_surfaces(f,off["wm"], "neighbor", min_distance=mdist)
                    firstpass = False
            
            hmu.log("Voxelizing")
            voxelize(off[key], ref, nii[key], merge=True)

        print("")    
        hmu.log("{:{}s}{:>10}", ("{} processing time:".format(key.upper()),
                             width-10, hmu.format_time(time.time()-start)))
        print("")
    
    # =========================================================================
    key = "bone"
    
    hmu.log("Processing {}", key.upper())
    hmu.log("="*width)   
    start = time.time()
    
    hmu.log("Preparing volume for meshing")   
    # add csf, binarize, and fill small holes
    mBONE_Y_mod = mrph.binary_fill_holes((cCSF.get_data()+mask[key].get_data()) > 0)
    hmu.write_nifti(mBONE_Y_mod, mod[key], ref, dtype=dtype)
    
    hmu.log("Extracting surface")
    off[key] = hmu.make_surface_mesh(mod[key], off[key], 1, 1, vertex_density, 5)[0]
    
    hmu.log("Decoupling from CSF")
    hmu.decouple_surfaces(off["csf"], off[key], "outer-from-inner",
                              cut_inner=True, min_distance=mdist) # this may alter the CSF mask
    
    hmu.log("Voxelizing")
    cBONE = hmu.make_volume_mask(off[key], nii[key], ref)
    
    hmu.log("Ensuring that the previously created surfaces are still "+
            "decoupled")
    check_intersect = '"{0}" "{1}" "{2}" --shells 2 --no-clean --intersect'.format(meshfix, "{0}", "{1}")
    if hmu.spawn_process(check_intersect.format(off["gm"], off["csf"]), True) == 0:
        hmu.log("Decoupling GM from CSF")
        hmu.decouple_surfaces(off["gm"], off["csf"], "inner-from-outer",
                              min_distance=mdist)   
        hmu.log("Voxelizing")                        
        cGM = hmu.make_volume_mask(off["gm"], nii["gm"], ref)
        
        if hmu.spawn_process(check_intersect.format(off["wm"], off["gm"]), True) == 0:
            hmu.log("Decoupling WM from GM")
            hmu.decouple_surfaces(off["wm"], off["gm"], "inner-from-outer")
            hmu.log("Voxelizing")
            cWM = hmu.make_volume_mask(off["wm"], nii["wm"], ref)
            cGM = hmu.make_volume_mask(off["gm"], nii["gm"], ref)
    
    print("")    
    hmu.log("{:{}s}{:>10}", ("{} processing time:".format(key.upper()),
                             width-10, hmu.format_time(time.time()-start)))
    print("")
    
    # =========================================================================
    key = "skin"
    
    hmu.log("Processing {}", key.upper())
    hmu.log("="*width)   
    start = time.time()
    
    hmu.log("Preparing volume for meshing") 
    # add csf, skull, MASK_AIR, binarize, and fill small holes
    mSKIN_Y_mod = mrph.binary_fill_holes((cCSF.get_data()+cBONE.get_data()+mask["air"].get_data()+mask[key].get_data()) > 0)
    hmu.write_nifti(mSKIN_Y_mod, mod[key], ref, dtype=dtype)
    
    hmu.log("Extracting surface")
    off[key] = hmu.make_surface_mesh(mod[key], off[key], 1, 1, vertex_density, 20)[0]
    
    hmu.log("Decoupling from bone")
    hmu.decouple_surfaces(off["bone"], off[key], "outer-from-inner",
                          min_distance=mdist)
    hmu.log("Voxelizing")
    cSKIN = hmu.make_volume_mask(off[key], nii[key], ref)
    
    print("")
    hmu.log("{:{}s}{:>10}", ("{} processing time:".format(key.upper()),
                             width-10, hmu.format_time(time.time()-start)))
    print("")
    
    # =========================================================================
    key = "eyes"
    
    hmu.log("Processing {}", key.upper())
    hmu.log("="*width)   
    start = time.time()
    
    hmu.log("Extracting surfaces")
    off[key] = hmu.make_surface_mesh(mask[key].get_filename(), off[key], 2, 1,
                          vertex_density, 5)
    
    hmu.log("Decoupling")
    for f in off[key]:
        fb = os.path.splitext(os.path.basename(f))[0]
        firstpass = True
        while firstpass or any_intersect_1st([f,off["bone"],off["skin"]]):
            hmu.log("Decoupling {} from bone", fb)
            hmu.decouple_surfaces(f, off["bone"], "neighbor", min_distance=mdist)
            hmu.log("Decoupling {} from skin", fb)
            hmu.decouple_surfaces(f, off["skin"], "inner-from-outer",
                                  min_distance=mdist) 
            firstpass = False
    
    hmu.log("Voxelizing")
    voxelize(off[key], ref, nii[key], merge=True)
    
    print("")
    hmu.log("{:{}s}{:>10}", ("{} processing time:".format(key.upper()),
                             width-10, hmu.format_time(time.time()-start)))
    print("")

    # =========================================================================
    key = "cavities"
    
    hmu.log("Processing INNER AIR {}", key.upper())
    hmu.log("="*width)
    start = time.time()
    
    hmu.log("Preparing volume for meshing")
    # erode skull (by one voxel), multiply by MASK_AIR to get inner air cavities,
    # binarize, and fill small holes
    mBONE_Y_ERO = mrph.binary_erosion(np.pad( cBONE.get_data()-cCSF.get_data(),1,mode="constant",constant_values=1), se)[1:-1,1:-1,1:-1]
    mAIR_IN_Y_mod = mrph.binary_fill_holes((mBONE_Y_ERO*mask["air"].get_data()) > 0)
    hmu.write_nifti(mAIR_IN_Y_mod, mod[key], ref, dtype=dtype)
    if mAIR_IN_Y_mod.any():
        hmu.log("Extracting surface(s)")
        off[key] = hmu.make_surface_mesh(mod[key],off[key], 9, 100, vertex_density, 5)
        
        hmu.log("Decoupling")
        for f in off[key]:
            fb = os.path.splitext(os.path.basename(f))[0]
            firstpass = True
            while firstpass or any_intersect_1st([f,off["csf"],off["bone"]]):
                hmu.log("Decoupling {} from CSF", fb)
                hmu.decouple_surfaces(f,off["csf"], "neighbor", min_distance=mdist)
                hmu.log("Decoupling {} from bone", fb)
                hmu.decouple_surfaces(f,off["bone"], "inner-from-outer",
                                      min_distance=mdist)
                firstpass = False
        
        hmu.log("Voxelizing")
        voxelize(off[key], ref, nii[key], merge=True)
    else:
        del off[key]
        del stl[key]
    print("")
    hmu.log("{:{}s}{:>10}", ("INNER AIR {} processing time:".format(key.upper()),
                             width-10, hmu.format_time(time.time()-start)))
    print("")
    
    # =========================================================================
    hmu.log("Processing OUTER AIR CAVITIES")
    hmu.log("="*width)
    start = time.time()
    
    key = "air_outer_relabel"
    
    # (for later deletion of air tetrahedra)
    hmu.log("Creating non-decoupled surfaces")
    hmu.log("Preparing volume for meshing")
    # multiply skin with MASK_AIR, binarize, subtract skull to get air cavities
    # not enclosed by the skull, and fill small holes
    mAIR_OUT_Y_mod = mrph.binary_fill_holes((((cSKIN.get_data()*mask["air"].get_data())>0).astype(np.int) - (cBONE.get_data()>0)) > 0)
    if mAIR_OUT_Y_mod.any():
        hmu.write_nifti(mAIR_OUT_Y_mod, mod[key], ref, dtype=dtype)
        
        hmu.log("Extracting surface(s)")
        off[key] = hmu.make_surface_mesh(mod[key], off[key], 9, 100, vertex_density, 5)
        #off["air_outer_relabel"] = sorted(glob(off["air_outer_relabel"]+"*.off"),
        #                                  key=str.lower)
    else:
        hmu.log("Nothing to extract")
        del off[key]
        del stl[key]
    print("")
    
    key = "air_outer_decoupled"
    
    hmu.log("Creating decoupled surfaces")
    hmu.log("Preparing volume for meshing")
    mBONE_Y_DIL = mrph.binary_dilation(cBONE.get_data(), se) # structuring element=??
    mSKIN_Y_ERO = mrph.binary_erosion(np.pad(cSKIN.get_data(),1,mode="constant",constant_values=1), se)[1:-1,1:-1,1:-1]
    mAIR_OUT_Y_decoupled = mrph.binary_fill_holes((((mask["air"].get_data()*mSKIN_Y_ERO)>0).astype(np.int) - mBONE_Y_DIL) > 0)
    if mAIR_OUT_Y_decoupled.any():
        hmu.write_nifti(mAIR_OUT_Y_decoupled, mod[key], ref, dtype=dtype)
        
        hmu.log("Extracting surface(s)")
        off[key] = hmu.make_surface_mesh(mod[key],off[key], 9, 100, vertex_density, 5)

        hmu.log("Decoupling")
        for f in off[key]:
            fb = os.path.splitext(os.path.basename(f))[0]
            firstpass = True
            while firstpass or any_intersect_1st([f,off["bone"],off["skin"]]):
                hmu.log("Decoupling {} from skin", fb)
                hmu.decouple_surfaces(f,off["skin"], "inner-from-outer", min_distance=mdist)
                hmu.log("Decoupling {} from bone", fb)
                hmu.decouple_surfaces(f,off["bone"], "neighbor", min_distance=mdist)
                firstpass = False
    else:
        hmu.log("Nothing to extract")
        del off[key]
        del stl[key]
    print("")
    hmu.log("{:{}s}{:>10}", ("OUTER AIR CAVITIES processing time:",
                             width-10, hmu.format_time(time.time()-start)))
    print("")
    
    # =========================================================================
    hmu.log("Converting surface meshes to STL")
    hmu.log("="*width)
    start = time.time()
    
    # Convert to STL
    for k in off.keys():
        if k is "air_outer_relabel":
            merge = False
        else:
            merge = True
        stl[k] = convert2stl(off[k], stl[k], merge)

    hmu.log("Done")
    
    relabel = stl.pop("air_outer_relabel",[])
    #surfaces = reduce(lambda x,y: x+y, stl.values()) 
    surfaces = []
    for v in stl.values():
        surfaces += v
    
    print("")
    hmu.log("{:{}s}{:>10}", ("Conversion time:", width-10,
                             hmu.format_time(time.time()-start)))

    return surfaces, relabel


def make_volume_mesh(subject_id,input_files,out_dir, keep_air):
    """Generates a volume mesh from the given set of surface meshes using Gmsh.
    This essentially fills the surface(s) with tetrahedra. 
    Relabels
    
    Uses Gmsh.
    
    PARAMETERS
    ----------
    
    RETURNS
    ----------
    paths to ...
    a volume mesh file (.msh)
    """   
    width = 72
    
    gmsh = hmu.path2bin("gmsh")
    
    tissues = ["wm", "gm", "csf", "bone", "skin", "eyes",
               "cavities", "air_outer", "ventricles"]
    files = {}
    
    for k in tissues:
        try:
            files[k] = [path for path in input_files if k in os.path.basename(path).lower()][0]
        except IndexError: # if the file does not exist    
            continue
            
    hmu.log("Creating .geo file")
    
    # tissue numbering
    num = {"wm":         1,
           "gm":         2,
           "csf":        3,
           "bone":       4,
           "skin":       5,
           "eyes":       6,
           "cavities":   7,
           "air_outer":  8,
           "ventricles": 9}    
    
    options = [
    "//to mesh the volumes call gmsh with",
    "//gmsh -3 -bin -o {0}.msh {0}.geo",
    "Mesh.Algorithm3D=4; //1=delaunay (tetgen) and 4=frontal (netgen)",
    "Mesh.Optimize=1;",
    "Mesh.OptimizeNetgen=1;"]
    
    merge = [
    'Merge "{wm}"; // first surface',
    'Merge "{gm}";',
    'Merge "{csf}";',
    'Merge "{bone}";',
    'Merge "{skin}";']
    
    surface_loop = [
    "Surface Loop(1) = {1}; // surface number on rhs; 1: wm.stl",
    "Surface Loop(2) = {2}; // gm.stl",
    "Surface Loop(3) = {3}; // csf.stl",
    "Surface Loop(4) = {4}; // bone.stl",
    "Surface Loop(5) = {5}; // skin.stl"]
    
    volumes = [
    "Volume(1) = {1};",
    "Volume(2) = {1, 2};       // GM volume",
    "Volume(3) = {2, 3};       // CSF volume: outside GM, inside CSF",
    "Volume(4) = {3, 4};       // bone volume",
    "Volume(5) = {4, 5};       // skin volume"]
    
    physical_surf = [
    "Physical Surface(1001) = {1}; // LHS: target surface region number, RHS: surface number (i.e. from merge ...)",
    "Physical Surface(1002) = {2};",
    "Physical Surface(1003) = {3};",
    "Physical Surface(1004) = {4};",
    "Physical Surface(1005) = {5};"]
    
    physical_vol = [
    "Physical Volume(1) = {1}; // LHS: target volume region number, RHS: volume number",
    "Physical Volume(2) = {2};",
    "Physical Volume(3) = {3};",
    "Physical Volume(4) = {4};",
    "Physical Volume(5) = {5};"]
    
    # if eyes.stl / cavities.stl / air_outer.stl exist, add them and modify the
    # .geo file accordingly
    n = 6
    if "eyes" in files:
        i = num["eyes"]
        merge.append('Merge \"{eyes}\";')
        surface_loop.append("Surface Loop({0}) = {{{1}}}; // eyes.stl".format(i,n))
        volumes.append("Volume({0}) = {{{0}}};          // eyes volume".format(i))      
        physical_surf.append("Physical Surface({0}) = {{{1}}}; // eyes".format(1000+i,i))
        physical_vol.append("Physical Volume({0}) = {{{0}}}; // eyes".format(i))
        # modify skin
        volumes[4] = "Volume(5) = {{4, 5, {0}}};    // SKIN volume".format(i)
        n += 1        
        
    if "cavities" in files:
        i = num["cavities"]
        merge.append('Merge \"{cavities}\";')
        if keep_air:
            i_new = 11

        surface_loop.append("Surface Loop({0}) = {{{1}}}; // cavities.stl".format(i,n))
        # modify bone volume
        volumes[3] = "Volume(4) = {{3, 4, {0}}};    // bone volume".format(i)
        if keep_air:
            physical_surf.append("Physical Surface({0}) = {{{1}}};".format(1000+i_new,i))
            physical_vol.append("Physical Volume({0}) = {{{1}}};".format(i_new,i))
            volumes.append("Volume({0}) = {{{0}}};          // cavities volume".format(i))
        else:
            physical_surf[3] = "Physical Surface(1004) = {{4, {0}}};".format(i)
        n += 1

    if "air_outer" in files:
        # outer air : s 1007/7
        i = num["air_outer"]
        merge.append('Merge \"{air_outer}\";')

        if keep_air:
            i_new = 12

        surface_loop.append("Surface Loop({0}) = {{{1}}}; // air_outer.stl; i.e., air outside bone".format(i,n))
        # modify skin
        if "eyes" in files:
            volumes[4] = "Volume(5) = {{4, 5, 6, {0}}}; // SKIN volume".format(i)
        else:
            volumes[4] = "Volume(5) = {{4, 5, {0}}};    // SKIN volume".format(i)

        if keep_air:
            physical_surf.append("Physical Surface({0}) = {{{1}}};".format(1000 + i_new, i))
            physical_vol.append("Physical Volume({0}) = {{{1}}};".format(i_new, i))
            volumes.append("Volume({0}) = {{{0}}};          // cavities volume".format(i))

        n += 1
        
    if "ventricles" in files:
        # ventricles are 1008/8
        i = num["ventricles"]
        merge.append('Merge \"{ventricles}\";')
        surface_loop.append("Surface Loop({0}) = {{{1}}}; // ventricles.stl".format(i,n))
        volumes.append("Volume({0}) = {{{0}}};          // ventricles volume".format(i))
        physical_surf.append("Physical Surface({0}) = {{{1}}}; // ventricles".format(1000+i,i))
        physical_vol.append("Physical Volume({0}) = {{{0}}}; // ventricles".format(i))
        # modify GM volume
        volumes[1] = "Volume(2) = {{1, {0}, 2}};    // GM volume".format(i)
        
    # make strings
    options       = "\n".join(options).format(subject_id)
    merge         = "\n".join(merge).format(**files)
    surface_loop  = "\n".join(surface_loop)
    volumes       = "\n".join(volumes)
    physical_surf = "\n".join(physical_surf)
    physical_vol  = "\n".join(physical_vol)

    # write the .geo file
    with open(os.path.join(out_dir,subject_id+".geo"),"w") as f:
        f.write(("\n"*2).join([options, merge, surface_loop, volumes,
                               physical_surf, physical_vol]))

    # create volume mesh
    section_start = time.time()
    hmu.log("Creating volume mesh using Gmsh")
    hmu.log("="*width)
    
    msh = os.path.join(out_dir, subject_id+".msh")
    geo = os.path.join(out_dir, subject_id+".geo")
    make_vol_mesh = '"{0}" -3 -bin -o "{1}" "{2}"'.format(gmsh, msh, geo)
    
    # first try Frontal meshing
    # zero is exit code for succes
    attempts = 1
    while True:
        #section_start = time.time()
        exitcode = hmu.spawn_process(make_vol_mesh, return_exit_status=True,
                                     verbose=True)
        
        # if gmsh fails meshing
        if exitcode is not 0:
            # if this was first attempt, try to remesh all surfaces and see if
            # this helps
            if attempts is 1:
                meshfix = hmu.path2bin("meshfix")
                uniformremesh_surface  = '"{0}" "{1}" -a 2.0 -u 1 -q --shells 9 --stl -o "{1}"'.format(meshfix,"{0}")
                print("")
                hmu.log("Gmsh failed meshing one or more surfaces. Doing uniform remeshing and retrying")
                for f in files.values():
                    hmu.log("Remeshing {}",os.path.basename(f))
                    hmu.spawn_process(uniformremesh_surface.format(f))
                    
                # Since WM and GM may be close (particularly when using CAT12),
                # check for intersections and decouple if necessary
                hmu.log("Decoupling")
                if "ventricles" in files:
                    decouple_wm_gm_ventricles(files["wm"], files["gm"],
                                              files["ventricles"])
                else:
                    hmu.decouple_surfaces(files["wm"], files["gm"],
                                      "outer-from-inner", cut_inner=True)
                print("")
            # if gmsh still fails, try Delaunay meshing (this may change the
            # surface to some extent). Modify the .geo file  
            elif attempts is 2:
                print("")
                hmu.log("Gmsh failed meshing one or more surfaces. Trying Delaunay meshing")
                print("")
                with open(geo,"r") as f:
                    body = f.read(-1)
                body = body.replace("Mesh.Algorithm3D=4","Mesh.Algorithm3D=1")
                with open(geo,"w") as f:
                    f.write(body)
            
            # if gmsh still fails
            elif attempts is 3:
                raise RuntimeError("Gmsh failed meshing one or more surfaces in {} attempts.".format(attempts))
        else:
            print("")
            hmu.log("Meshing successful in {} attempt(s)",attempts)
            hmu.log("{:{}s}{:>10}", ("Gmsh running time:", width-10,
            hmu.format_time(time.time()-section_start)))
            print("")
            break
        attempts += 1
    
    return msh
     
    
def visualize(subject_id,surface_meshes,T1,control_mask_initial,
              control_mask_final,subject_dir,force_remake=False,show_mni=False):
    """Visualize results of the meshing in FreeView by overlaying the surfaces
    on the initial and the final (combined) tissue masks.
    
    PARAMETERS
    ----------
    subject_id : str
        subject ID.
    surface_meshes : array_like
        list of surface filenames.
    T1 : str or nibabel object
        filename of T1 image.
    control_mask_initial : str
        filename of initial (combined) tissue mask (before meshing).
    control_mask_final : str
        filename of final (combined) tissue mask (i.e., the masks recreated
        from the surface meshes).
    subject_dir : str
        path in which to find the files required by mris_transform (unity.xfm
        and ref_FS.nii.gz). This is also the path to which two FreeView meshes
        are saved.
    force_remake : boolean
        do not re-use the surface meshes created for visual control
    show_mni: boolean
        show the T1 registered to MNI space for control (taken from toMNI subdirectory)
    
    RETURNS
    ----------
    (opens FreeView to show results)
    """
    
    if type(T1)==str:
        T1 = nib.load(T1)
        
    T1_max_val = np.max(T1.get_data())
    
    fmesh_brain = os.path.join(subject_dir,"brain_contr.fsmesh")
    fmesh_head  = os.path.join(subject_dir,"head_contr.fsmesh")
    headoff     = os.path.join(subject_dir,"head.off")
    ref_fs     = os.path.join(subject_dir,"ref_FS.nii.gz")
    
    meshfix = hmu.path2bin("meshfix")
    
    # meshfix commands
    combine_meshes = '"{0}" "{1}" "{2}" -q --no-clean --shells {3} -o "{4}"'.format(meshfix, "{0}", "{1}", "{2}", "{3}")
    def combine_meshes_fsmesh(surf_a, surf_b, out):
        s_a = mesh_io.read_stl(surf_a)
        s_b = mesh_io.read_stl(surf_b)
        c = s_a.join_mesh(s_b)
        mesh_io.write_freesurfer_surface(c, out, ref_fs=ref_fs)

    visualize_with_freeview = "freeview -v "+\
          '"{0}":colormap=grayscale:grayscale=0,{1} '\
          '"{2}":colormap=nih:colorscale=0.1,6.1:opacity=0 '\
          '"{3}":colormap=nih:colorscale=0.1,6.1:opacity=0.3 '\
          '-f '\
          '"{4}":color=128,128,128:edgecolor=white '\
          '"{5}":color=255,166,133:edgecolor=gray'.format(
                                      T1.get_filename(), T1_max_val,
                                      control_mask_initial, control_mask_final,
                                      fmesh_brain, fmesh_head)
    
    stl = {}.fromkeys(["wm", "gm", "csf","eyes", "bone", "skin", "ventricles"])
    for k in stl.keys():
        try:
            stl[k] = [path for path in surface_meshes if k in os.path.basename(path).lower()][0]
        except IndexError: # if the file does not exist    
            stl[k] = ""
            
    if force_remake | (not os.path.isfile(fmesh_brain)):
        hmu.log("Combining surface meshes for visualization")
        hmu.log("Making brain mesh")
        combine_meshes_fsmesh(stl['wm'], stl['gm'], fmesh_brain)
    if force_remake | (not os.path.isfile(fmesh_head)):
        i=1
        hmu.log("Making head mesh {}".format(i))
        hmu.spawn_process(
            combine_meshes.format(stl["csf"], stl["eyes"], 10, headoff))
        i+=1
        if stl["ventricles"] is not "":
            hmu.log("Making head mesh {}".format(i))
            hmu.spawn_process(
            combine_meshes.format(headoff, stl["ventricles"], 10, headoff))
            i+=1
        hmu.log("Making head mesh {}".format(i))
        hmu.spawn_process(
            combine_meshes.format(headoff, stl["bone"], 10, headoff))
        i+=1
        hmu.log("Making head mesh {}".format(i))
        combine_meshes_fsmesh(headoff, stl['skin'], fmesh_head)
    remove_files(headoff)

    hmu.log("Opening freeview. Please wait...")
    hmu.spawn_process(visualize_with_freeview, new_thread=True)

    if show_mni:
        ref_mni = file_finder.templates.mni_volume
        if not os.path.isfile(ref_mni):
                raise IOError('Could not find reference volume in mni space at: {0}'.format(
                        ref_mni))

        T1_mni_nonlin = os.path.join(subject_dir,"toMNI","T1fs_nu_nonlin_MNI.nii.gz")
        if not os.path.isfile(T1_mni_nonlin):
                raise IOError('Could not find MNI-transformed T1: {0}'.format(T1_mni_nonlin))

        T1_mni_12dof = os.path.join(subject_dir,"toMNI","T1fs_nu_12DOF_MNI.nii.gz")
        if not os.path.isfile(T1_mni_12dof):
                raise IOError('Could not find MNI-transformed T1 (12dof): {0}'.format(T1_mni_12dof))

        T1 = nib.load(T1_mni_nonlin)
        T1_max_val = np.max(T1.get_data())
        visualize_MNI_with_freeview = 'freeview -v "{0}" "{1}":visible=0:colormap=grayscale:grayscale=0,{2} '\
                                      '"{3}":colormap=grayscale:grayscale=0,{2}'.format(ref_mni, T1_mni_12dof,
                                      T1_max_val,T1_mni_nonlin)
        hmu.spawn_process(visualize_MNI_with_freeview, new_thread=True)


def visualize_spm(subject_id,T1,control_mask_final,subject_dir,show_mni=False):
    """Visualize results of the segmentation using spm for basic quality check.
    
    PARAMETERS
    ----------
    subject_id : str
        subject ID.
    T1 : str or nibabel object
        filename of T1 image.
    control_mask_final : str
        filename of final (combined) tissue mask (i.e., the masks recreated
        from the surface meshes).
    subject_dir : str
        path in which to find the files required by mris_transform (unity.xfm
        and ref_FS.nii.gz). This is also the path to which two FreeView meshes
        are saved.
    show_mni: boolean
        show the T1 registered to MNI space for control (taken from toMNI subdirectory)
    
    RETURNS
    ----------
    (opens spm to show some results for basic quality check)
    """
              
    if type(T1)==str:
        T1 = nib.load(T1)

    # save unzipped versions
    if not os.path.isdir(os.path.join(subject_dir,"tmp")):
        os.mkdir(os.path.join(subject_dir,"tmp"))
    T1spm_name=os.path.join(subject_dir,"tmp","T1fs_conform_spm.nii")
    maskspm_name=os.path.join(subject_dir,"tmp","final_contr_spm.nii")
    nib.save(T1, T1spm_name)
    nib.save(nib.load(control_mask_final), maskspm_name)

    cmd = "matlab -nosplash -nodesktop -r "
    path2mfiles = os.path.join(SIMNIBSDIR,"resources","spm12")
    cmd += "\"addpath('{}');".format(path2mfiles)
    cmd += "try,spm_check_registration(char('{0}','{1}')); uiwait(gcf);catch ME,rethrow(ME);end,exit;\"".format(
            T1spm_name,maskspm_name)
    hmu.log("Opening spm. Please wait...")
    if (sys.platform == 'win32'):
        p = hmu.spawn_process(cmd, verbose=True, return_exit_status=True)
    else:
        p = hmu.spawn_process(cmd, new_process=True)
    
    if show_mni:
        ref_mni = file_finder.templates.mni_volume
        if not os.path.isfile(ref_mni):
                raise IOError('Could not find reference volume in mni space at: {0}'.format(
                        ref_mni))

        T1_mni_nonlin = os.path.join(subject_dir,"toMNI","T1fs_nu_nonlin_MNI.nii.gz")
        if not os.path.isfile(T1_mni_nonlin):
                raise IOError('Could not find MNI-transformed T1: {0}'.format(T1_mni_nonlin))

        ref_mni_spm_name=os.path.join(subject_dir,"tmp","MNI_template_spm.nii")
        T1_mni_nonlin_spm_name=os.path.join(subject_dir,"tmp","T1fs_nu_nonl_spm_MNI.nii")
        nib.save(nib.load(ref_mni), ref_mni_spm_name)
        nib.save(nib.load(T1_mni_nonlin), T1_mni_nonlin_spm_name)

        cmd = "matlab -nosplash -nodesktop -r "
        cmd += "\"addpath('{}');".format(path2mfiles)
        cmd += "try,spm_check_registration(char('{0}','{1}')); uiwait(gcf);catch ME,rethrow(ME);end,exit;\"".format(
                ref_mni_spm_name,T1_mni_nonlin_spm_name)

        if (sys.platform == 'win32'):
            p2 = hmu.spawn_process(cmd, verbose=True, return_exit_status=True)
        else:
            p2 = hmu.spawn_process(cmd, new_process=True)
            p2.join()
            
    if (sys.platform != 'win32'):
        p.join()


def parse_args(argv):
    """Parse arguments from command line to python function headmodel.
    
    PARAMETERS
    ----------
    argv : str
        Command line arguments.
    
    RETURNS
    ----------
    args : Namespace object
        Object with arguments ordered in fields, e.g., args.argument1,
        args.argument2 etc.
    """
    import argparse
    from argparse import RawDescriptionHelpFormatter 
    
    # kwargs definitions
    program_headmodel = {
        "prog"       : "headreco",
        "description": """Prepare a volume mesh from structural images. The
                       behavior of the program is determined by the argument
                       mode.""",
        "epilog"     : """For help on each mode run: headreco mode -h. Use
                       "headreco all -h" to get information on the standard
                        way of using it."""
        }
        
    """
    First argument should be the mode in which to run the program, i.e. what
    to do followed by the subject ID. Depending on which mode is chosen, 
    various optional arguments may be supplied. E.g.,
    
    headreco all subject_id -d custom --vx-size 0.8 0.8 0.8 file1.nii.gz
        file2.nii.gz -v 0.4 -c
    
    which would create a folder m2m_<subject_id> in the current working directory
    and run all steps necessary to generate a volume mesh. The structural
    image(s) supplied (file1 and file2) would be prepared and resampled to a
    voxel size of [0.8 0.8 0.8]. The surface meshes would be created with a
    density of 0.4 nodes per mm2 and finally all temporary files would be
    deleted.
    """
    
    # modes
    mode = {
        "metavar":"mode",
        "help"   :"Mode in which to run headreco."        
        }
    mode_all = {
        "help":"""Run all steps necessary for generating a volume mesh from
               structural images.""",
        "description":"""Run all steps necessary for generating a volume mesh
                      from structural images. This includes
                      (1) preparation and segmentation of the input volumes,
                      including an (optional, but recommended) cortical surface
                      segmentation
                      (2) cleaning of segmented volumes and surfaces
                      (3) surface meshing and decoupling, and 
                      (4) volume meshing.
                      Thus, it is equivalent to manually running preparevol,
                      (optional: preparecat), cleanvols, surfacemesh,
                      volumemesh (in this order).
                      """
        }
    mode_preparevols = {
        "help":"Prepare input files, segment and generate binarized tissue masks.",
        "description":"""
         This step runs SPM12 and optionally CAT12 to segment
         the input images. The detailed steps are:
           1) The input images are reoriented to LAS and
              optionally resampled. Exception: setting option
              "-d no-conform" will keep the T1 unchanged
           2) SPM12 is called to coregister the T2 to the T1,
              generate the transformations to MNI space and
              segment the images; this uses the extended tissue
              priors that cover also the neck. In addition,
              a number of MNI priors for eyes etc are warped
              to subject space for use in the later generation
              of the tissue masks
           3) CAT12 is called (optional, but recommended) to
              segment the cortex; In addition, a number of MNI
              priors are warped to subject space for use in the
              later joining of surfaces of the left and right
              hemispheres
           4) the posterior probability masks are binarized

         result files and directories:
             T1fs_conform, T2_conform: reoriented and resampled
                    input images. The T2 is coregistered to the 
                    T1.
             segment/spm/binarized: binarized posterior 
                    probability masks of the SPM segmentation
             segment/cat/surf: cortical surface reconstructions
             segment/cat/binarized: binarized posterior 
                    probability masks of GM, WM and CSF
                    of CAT segmentation
             templates: spmprior_* are created during SPM12
                    segmentation process, and are later on
                    used to, e.g. separate CSF from eyes
                        cattemplate_* are created during CAT12
                    segmentation process, and are later on
                        used for building a mid-brain-brainstem-
                        cerebellum mask that is joined with
                        the left and right cortical GM surfaces
               """
        }
    mode_preparecat = {
        "help" : """Process the CAT segmentation results so as to prepare them
                 for volume meshing.""",
        "description" : """
        This step
         1) expands the central surfaces created by
            CAT12 to get the pial surfaces
         2) creates WM (including brainstem...), ventricles 
            and combined corpus-callosum-brainstem-cerebellum
            surfaces from the voxel segmentations
         3) joins the pial surfaces using the combined
            corpus-callosum-brainstem-cerebellum surface
         4) decouples the GM, WM and ventricle surfaces
         5) creates volume masks for those surfaces
        
         Final files stored in mask_prep:
         cat_surf_gm.off, cat_surf_wm.off, cat_surf_ventricles1,2,...off
         MASK_CAT_GM.nii.gz, MASK_CAT_WM.nii.gz, MASK_CAT_VENTRICLES.nii.gz

         Note: in order to convert a .off file in a .stl file for 
               visual inspection in gmsh, run:
         meshfix cat_surf_gm.off --no-clean --stl -o cat_surf_gm.stl"""
        }
    mode_cleanvols = {
        "help"       :"""Clean tissue masks thus preparing them for surface
                      meshing.""",
        "description":"""
         This step cleans the binarized tissue masks
         in segment/spm/binarized and creates masks which contain
         all inner tissue types as well ("volume decoupling").
         
         When GM, WM and ventricle masks were created from the CAT
         results in the step before, they are used for volume
         decoupling the rest of the masks. They themselves are 
         kept unchanged.
        
         Final files are stored in mask_prep:
         MASK_*.nii.gz.
                      """
        } 
    mode_surfacemesh = {
        "help"       :"Generate surfaces meshes of each tissue class.",
        "description":"""
         This step creates the surface .stl files ready for 
         volume-meshing. The surfaces are build up in a 
         inner-to-outer order, i.e. the next-outer surface
         is decoupled from the next-inner one.
        
         The following steps are applied for each surface:
         1) decouple its voxel mask from the next-inner
         voxel mask created from final surface mesh (see step 4)
         2) build a surface from the updated voxel mask
         3) decouple surface from next-inner surface
         4) update its voxel mask using the decoupled surface
        
         In case surface decoupling of bone from csf fails,
         a decoupling of csf from bone is tried (i.e. 
         the order is to reversed to outer-to-inner, including
         checks for gm and wm). This was included, as "spikes"
         in csf can sometimes prevent a successful decoupling
         of bone from csf.
         
         When CAT12 was used, then the GM, WM and ventricle
         surfaces and their volume masks are only copied
         from the mask_prep folder, without any additional
         decoupling steps

         The final files (stl surfaces and corresponding
         nifti voxel masks) are stored in the m2m_{subid}
         folder.
         """
        }
    mode_volumemesh = {
        "help"       :"Generate a volume mesh from surface meshes.",
        "description":"""
         This step creates the final volume mesh from the 
         surface .stl files
        
         The following steps are applied for each surface:
         1) a {subid}.geo file is created, and gmsh is called for
            meshing
         2) after meshing, tetrahedra corresponding to air are
            deleted and the tissue labels are updated
            to conform with the conventions used for the
            FEM calculations; in addition, triangle and tetrahedra
            orientations are checked and, if needed, updated; also
            very thin tetrahedra are made slightly thicker to improve
            the mesh quality
         3) a nifti-label-image is created from the final mesh; 
            in addition, final wm and gm volume masks are extracted
            from the mesh and converted (together with the T1) to MNI,
            to allow for visual control segmentation and
            co-registration accuracy

         The final mesh is stored in the same directory containing
         the m2m_{subid} folder"""
        }
    mode_check = {
        "help":"Visually inspect the results (FreeView is strongly recommended).",
        "description":"""Visually inspect the results. Surfaces will be only
                         displayed when FreeView is installed and in the path.
                         This is the recommended usage. Otherwise, spm check reg
                         is used for some basic quality checks."""
        }    
    
    # other positional arguments
    subject_id = {
        "help":"""Subject ID. All files pertaining to a given subject is stored
               in the folder <m2m_subject_id> which initially is created in the
               current working directory."""
        }
    input_files = {
        "nargs":"+",
        "help" :"Input files, i.e. structural images"
        }
        
    # optional arguments
    cat = {
        "action"  : "store_false",
        "dest"    : "cat",
        "help"    : """Don't use CAT12 for brain segmentation."""
        }
    biasreg = {
        "type"    : float,
        "default" : 0.01,
        "help"    : """Regularization parameter for T2 bias correction in the
                    SPM segmentation. Only applies if a T2 weighted image is
                    provided. See SPM for extended documentation (default:
                    %(default)s)."""
        }
    dsf = {
        "type"    : int,
        "default" : 3,
        "help"    : """Downsampling factor in the SPM segmentation. A value of
                    zero corresponds to no downsampling. See SPM for extended
                    documentation (default: %(default)s)."""
        }
    dimension =  {
        "choices":["keep","standard","custom","no-conform"],
        "default":"no-conform",
        "help"   :"""How to prepare input image(s). keep (reorient to LAS and
                  do not resample), standard (reorient to LAS and resample to
                  voxel size of [1 1 1] mm and image dimensions [256 256 256]
                  voxels), custom (reorient to LAS and resample to a voxel size
                  as specified by --vx-size), no-conform (neither reorient to
                  LAS nor resample). (default: %(default)s)."""
        }
    vx_size = {
        "type"   :float,
        "nargs"  :3,
        "default": None,
        "metavar":("x","y","z"),
        "help"   :"Voxel size (mm)" 
        }
    vertex_density = {
        "type"   :float,
        "default":0.5,
        "help"   :"""Vertex density (nodes per mm^2) of the surface meshes.
                  Lower values will result in a coarser mesh, less disk usage,
                  and faster computations whereas the opposite is true for
                  higher values (default: %(default)s)."""
        }
    noclean = {
        "action" :"store_true",
        "default":False,
        "help"   :"""Do not remove temporary files (default: %(default)s)."""
        }
    force_remake = {
        "action" :"store_true",
        "default":False,
        "help"   :"""Force remake of all files in the current step (default:
                  %(default)s)."""
        }
    skip_coreg = {
        "action" :"store_true",
        "default":False,
        "help"   :"""Skip registration of T1 and T2 if already registered (default:
                  %(default)s)."""
        }
    cat_print = {
        "action" :"store_true",
        "default":False,
        "help"   :"""Print pdf from CAT12 (default:
                  %(default)s)."""
        }
    keep_air = {
        "action": "store_true",
        "default": False,
        "help": """Keep the surfaces and volumes of air cavities (default:
                      %(default)s)."""
    }
         
   
    # setup argument parser
    parser = argparse.ArgumentParser(**program_headmodel)
    parser.add_argument('--version', action='version', version=__version__)
    subparser_mode = parser.add_subparsers(dest="mode",**mode, required=True)

    # mode "all"    
    parser_all = subparser_mode.add_parser("all", **mode_all)
    # optional
    parser_all.add_argument("--biasreg", **biasreg)
    parser_all.add_argument("--no-cat", **cat)
    parser_all.add_argument("-d","--dimension", **dimension)     
    parser_all.add_argument("--dsf", **dsf)
    parser_all.add_argument("--vx-size", **vx_size)
    parser_all.add_argument("-v","--vertex-density", **vertex_density)
    parser_all.add_argument("--noclean", **noclean)
    parser_all.add_argument("--skip_coreg", **skip_coreg)
    parser_all.add_argument("--cat_print", **cat_print)
    parser_all.add_argument("--keep_air", **keep_air)
    # positional
    parser_all.add_argument("subject_id",**subject_id)
    parser_all.add_argument("input_files", **input_files)
    
    # mode "preparevols"   
    parser_prepvols = subparser_mode.add_parser("preparevols", **mode_preparevols)
    parser_prepvols.formatter_class=RawDescriptionHelpFormatter
    # optional
    parser_prepvols.add_argument("--no-cat", **cat)
    parser_prepvols.add_argument("--biasreg", **biasreg)
    parser_prepvols.add_argument("--dsf", **dsf)
    parser_prepvols.add_argument("-d","--dimension", **dimension)       
    parser_prepvols.add_argument("--vx-size", **vx_size)
    parser_prepvols.add_argument("--noclean", **noclean)
    parser_prepvols.add_argument("--skip_coreg", **skip_coreg)
    parser_prepvols.add_argument("--cat_print", **cat_print)
    # positional
    parser_prepvols.add_argument("subject_id",**subject_id)
    parser_prepvols.add_argument("input_files", **input_files)
    
    # mode "preparecat"   
    parser_prepcat = subparser_mode.add_parser("preparecat", **mode_preparecat)
    parser_prepcat.formatter_class=RawDescriptionHelpFormatter
    # optional
    parser_prepcat.add_argument("-v","--vertex-density", **vertex_density)
    parser_prepcat.add_argument("--noclean", **noclean)
    parser_prepcat.add_argument("--no-cat", **cat)
    # positional
    parser_prepcat.add_argument("subject_id",**subject_id)
    
    # mode "cleanvols"   
    parser_cleanvols  = subparser_mode.add_parser("cleanvols", **mode_cleanvols)
    parser_cleanvols.formatter_class=RawDescriptionHelpFormatter
    # optional
    parser_cleanvols.add_argument("--no-cat", **cat)
    parser_cleanvols.add_argument("--noclean", **noclean)
    # positional
    parser_cleanvols.add_argument("subject_id",**subject_id)

    # mode "surfacemesh"  
    parser_surfmesh = subparser_mode.add_parser("surfacemesh", **mode_surfacemesh)
    parser_surfmesh.formatter_class=RawDescriptionHelpFormatter
    # optional
    parser_surfmesh.add_argument("--no-cat", **cat)
    parser_surfmesh.add_argument("-v","--vertex-density", **vertex_density)
    parser_surfmesh.add_argument("--noclean", **noclean)
    # positional
    parser_surfmesh.add_argument("subject_id",**subject_id)
    
    # mode "volumemesh"                                                    
    parser_volmesh = subparser_mode.add_parser("volumemesh", **mode_volumemesh)
    parser_volmesh.formatter_class=RawDescriptionHelpFormatter
    # optional
    parser_volmesh.add_argument("--no-cat", **cat)
    parser_volmesh.add_argument("-nc","--noclean", **noclean)
    parser_volmesh.add_argument("--keep_air", **keep_air)
    # positional
    parser_volmesh.add_argument("subject_id",**subject_id)
    
    # mode "check"         
    parser_check = subparser_mode.add_parser("check", **mode_check)
    # optional
    parser_check.add_argument("--noclean", **noclean)
    parser_check.add_argument("-f","--force-remake", **force_remake)
    # positional
    parser_check.add_argument("subject_id",**subject_id)
    
    # get the arguments
    args = parser.parse_args(argv)
    
    return args

    
def make_handler(handler_type="stream", level=logging.INFO, filename=None):
    """Create a stream handler (logging to console) or a file handler (logging
    to a file) for a logger.

    PARAMETERS
    ----------
    handler_type : str
        The handler type, either "stream" or "file" (default = "stream").
    level : logging level
        Logging level (default = logging.INFO)
    filename : str
        The name of the file to which to write logs (only applicable if
        handler_type is "file")
        
    RETURNS
    ---------
    handler : handler object
        Handler object which can be added to a logger.
    """
    
    if handler_type == "stream":
        formatter = logging.Formatter("%(message)s")
        handler = logging.StreamHandler()

    if handler_type == "file":
        formatter = logging.Formatter(fmt="%(asctime)s [%(levelname)-5.5s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        if filename == None:
            filename = os.path.join(os.getcwd(), "log.log")
        handler = logging.FileHandler(filename, mode='a')

    handler.setLevel(level)
    handler.setFormatter(formatter)

    return handler


def log_excep(*args):
    """Log uncaught exceptions. When an exception occurs, sys.exc_info()
    returns a tuple of three variables (exception class, exception value,
    traceback). Setting
        sys.excepthook = log_excep
    will replace the standard way of handling exceptions but that of log_excep.
    log_excep takes the sys.exc_info() as input and prints the exception to 
    "logger" at level error.
    """
    logger.error("".join(traceback.format_exception(*args)).rstrip())


def remove_dirs(dirs):
    """Remove one or more directories if they exist.
    
    PARAMETERS
    ----------
    dirs : list or str
        List of directories to remove. If only one directory, a string may be
        specified.
    """
    if type(dirs) == str:
        dirs = [dirs]
        
    for d in dirs:
        try:
            shutil.rmtree(d)
        except OSError:
            pass

            
def remove_files(files):
    """Remove one or more files if they exist.
    
    PARAMETERS
    ----------
    files : list or str
        List of files to remove. If only one file, a string may be specified.
    """
    if type(files) == str:
        files = [files]
        
    for f in files:
        try:
            os.remove(f)
        except OSError:
            pass
            
# use log_excep as exception handler instead of python standard
sys.excepthook = log_excep



import logging
import os
import shutil
import time
import subprocess
import nibabel as nib
import glob
import sys
import tempfile
import re
import numpy as np

from simnibs import SIMNIBSDIR

from . import samseg
from . import charm_utils
from .. import __version__
from .. import utils
from ..utils.simnibs_logger import logger
from ..utils import file_finder
from ..utils import transformations
from ..utils.transformations import crop_vol
from ..utils.spawn_process import spawn_process
from ..mesh_tools.meshing import create_mesh
from ..mesh_tools.mesh_io import (
    read_gifti_surface,
    write_msh,
    read_off,
    write_off,
    Msh,
    ElementData,
)
from ..utils import cond_utils
from ..utils import plotting, html_writer


def run(
    subject_dir=None,
    T1=None,
    T2=None,
    registerT2=False,
    initatlas=False,
    segment=False,
    create_surfaces=False,
    mesh_image=False,
    usesettings=None,
    noneck=False,
    init_transform=None,
    use_transform=None,
    force_qform=False,
    force_sform=False,
    options_str=None,
    debug=False,
):
    """charm pipeline

    PARAMETERS
    ----------
    subject_dir : str, mandatory
        output directory
    T1 : str
        filename of T1 image
    T2 : str
        filename of T2 image

    --> parameters to control the workflow:
    registerT2 : bool
        run T2-to-T1 registration (default = False)
    initatlas : bool
        run affine registration of atlas to input images (default = False)
    segment : bool
        run volume and surface segmentation (default = False)
    mesh_image : bool
        run tetrahedral meshing (default = False)
    --> further parameters:
    use_transform: path-like
        Transformation matrix used instead of the affine registration of the
        MNI template to the subject MRI, i.e., it takes the MNI template *to*
        subject space. Supplied as a path to a space delimited .txt file
        containing a 4x4 transformation matrix (default = None, corresponding
        to the identity matrix).
    init_transform: path-like
        Transformation matrix used to initialize the affine registration of the
        MNI template to the subject MRI, i.e., it takes the MNI template *to*
        subject space. Supplied as a path to a space delimited .txt file
        containing a 4x4 transformation matrix (default = None, corresponding
        to the identity matrix).
    usesettings : str
        filename of alternative settings-file (default = None)
    options_str : str
        string of command line options to add to logging (default = None)

    RETURNS
    ----------
        None
    """
    # ------------------------START UP-----------------------------------------
    start = time.time()

    # Initialize subject directory
    if not os.path.exists(subject_dir):
        os.mkdir(subject_dir)
    sub_files = file_finder.SubjectFiles(subpath=subject_dir)

    if force_qform and force_sform:
        raise ValueError("Can't force both q- and s-forms, please use only one of the flags.")

    _setup_logger(sub_files.charm_log)
    logger.info(f"simnibs version {__version__}")
    logger.info(f"charm run started: {time.asctime()}")
    logger.debug(f"options: {options_str}")

    settings = _read_settings_and_copy(usesettings, sub_files.settings)
    denoise_settings = settings["preprocess"]
    do_denoise = denoise_settings["denoise"]
    samseg_settings = settings["samseg"]
    logger.debug(settings)

    if init_transform:
        init_transform = _read_transform(init_transform)

    if use_transform:
        use_transform = _read_transform(use_transform)

    _prepare_t1(T1, sub_files.reference_volume, force_qform, force_sform)
    _prepare_t2(T1, T2, registerT2, sub_files.T2_reg, force_qform, force_sform)

    # -------------------------PIPELINE STEPS---------------------------------
    # TODO: denoise T1 here with the sanlm filter, T2 denoised after coreg.
    # Could be before as well but doesn't probably matter that much.
    # Denoising after has the benefit that we keep around the non-denoised
    # versions of all inputs.

    # Make the output path for segmentations
    os.makedirs(sub_files.segmentation_folder, exist_ok=True)

    if do_denoise:
        _denoise_inputs(T1, T2, sub_files)

    # Set-up samseg related things before calling the affine registration
    # and/or segmentation

    # Specify the maximum number of threads the GEMS code will use
    # by default this is all, but can be changed in the .ini
    num_threads = samseg_settings["threads"]
    if isinstance(num_threads, int) and num_threads > 0:
        samseg.setGlobalDefaultNumberOfThreads(num_threads)
        logger.info("Using %d threads, instead of all available." % num_threads)

    # TODO: Setup the visualization tool. This needs some pyqt stuff to be
    # installed. Don't know if we want to expose this in the .ini
    showFigs = False
    showMovies = False
    visualizer = samseg.initVisualizer(showFigs, showMovies)

    (
        template_name,
        atlas_settings,
        atlas_path,
        atlas_level1,
        atlas_level2,
        atlas_affine_name,
        gmm_parameters,
    ) = _setup_atlas(samseg_settings, sub_files.T2_reg, usesettings)

    if initatlas:
        # initial affine registration of atlas to input images,
        # including break neck
        logger.info("Starting affine registration and neck correction.")
        inputT1 = sub_files.T1_denoised if do_denoise else T1

        if use_transform is not None:
            logger.info("Using world-to-world transform provided by user.")
            trans_mat = use_transform
        else:
            if samseg_settings["init_type"] == "atlas":
                trans_mat = None
            elif samseg_settings["init_type"] == "mni":
                mni_template = file_finder.Templates().mni_volume
                mni_settings = settings["initmni"]
                trans_mat = charm_utils._init_atlas_affine(
                    inputT1, mni_template, mni_settings
                )
            else:
                logger.info("Affine initialization type unknown. Defaulting to 'atlas'")
                trans_mat = None


        charm_utils._register_atlas_to_input_affine(
            inputT1,
            template_name,
            atlas_affine_name,
            atlas_level1,
            atlas_level2,
            sub_files.segmentation_folder,
            sub_files.template_coregistered,
            settings["initatlas"],
            atlas_settings["neck_optimization"]["neck_tissues"],
            visualizer,
            noneck,
            world_to_world_transform_matrix=trans_mat,
            init_transform=init_transform,
        )

    if segment:
        # This part runs the segmentation, upsamples bias corrected output,
        # writes mni transforms, creates upsampled segmentation, maps tissues
        # to conductivities, runs morphological operations

        # Run the segmentation and return the class, which is needed
        # for further post-processing
        # The bias field kernel size has to be changed based on input
        input_images = []
        input_images.append(sub_files.T1_denoised if do_denoise else T1)
        if os.path.exists(sub_files.T2_reg):
            input_images.append(
                sub_files.T2_reg_denoised if do_denoise else sub_files.T2_reg
            )

        segment_settings = settings["segment"]
        logger.info("Estimating parameters.")
        segment_parameters_and_inputs = charm_utils._estimate_parameters(
            sub_files.segmentation_folder,
            sub_files.template_coregistered,
            atlas_path,
            input_images,
            segment_settings,
            gmm_parameters,
            visualizer,
            gmm_params_path=sub_files.segmentation_folder if debug else None
        )

        # Okay now the parameters have been estimated, and we can segment the
        # scan. However, we need to also do this at an upsampled resolution,
        # so first write out the bias corrected scan, and the segmentation.

        bias_corrected_image_names = [sub_files.T1_bias_corrected]
        if len(input_images) > 1:
            bias_corrected_image_names.append(sub_files.T2_bias_corrected)

        logger.info("Writing out normalized images and labelings.")
        os.makedirs(sub_files.surface_folder, exist_ok=True)
        cat_images = [
            sub_files.norm_image,
            sub_files.cereb_mask,
            sub_files.subcortical_mask,
            sub_files.hemi_mask,
        ]

        cat_structs = atlas_settings["CAT_structures"]
        tissue_settings = atlas_settings["conductivity_mapping"]
        csf_factor = segment_settings["csf_factor"]
        samseg.simnibs_segmentation_utils.writeBiasCorrectedImagesAndSegmentation(
            bias_corrected_image_names,
            sub_files.labeling,
            segment_parameters_and_inputs,
            tissue_settings,
            csf_factor,
            cat_structure_options=cat_structs,
            cat_images=cat_images,
        )

        fn_LUT = sub_files.labeling.rsplit(".", 2)[0] + "_LUT.txt"
        shutil.copyfile(file_finder.templates.labeling_LUT, fn_LUT)

        # Write out MNI warps
        logger.info("Writing out MNI warps.")
        os.makedirs(sub_files.mni_transf_folder, exist_ok=True)
        samseg.simnibs_segmentation_utils.saveWarpField(
            template_name,
            sub_files.mni2conf_nonl,
            sub_files.conf2mni_nonl,
            segment_parameters_and_inputs,
        )

        # Run post-processing
        logger.info("Post-processing segmentation")
        os.makedirs(sub_files.label_prep_folder, exist_ok=True)

        upsampled_image_names = [sub_files.T1_upsampled]
        if len(bias_corrected_image_names) > 1:
            upsampled_image_names.append(sub_files.T2_upsampled)

        cleaned_upsampled_tissues = charm_utils._post_process_segmentation(
            bias_corrected_image_names,
            upsampled_image_names,
            tissue_settings,
            csf_factor,
            segment_parameters_and_inputs,
            sub_files.template_coregistered,
            atlas_affine_name,
            sub_files.tissue_labeling_before_morpho,
            sub_files.upper_mask,
            debug=debug,
        )

        # Write to disk
        upsampled_image = nib.load(sub_files.T1_upsampled)
        affine_upsampled = upsampled_image.affine
        upsampled_tissues = nib.Nifti1Image(cleaned_upsampled_tissues, affine_upsampled)
        nib.save(upsampled_tissues, sub_files.tissue_labeling_upsampled)
        del cleaned_upsampled_tissues

        fn_LUT = sub_files.tissue_labeling_upsampled.rsplit(".", 2)[0] + "_LUT.txt"
        shutil.copyfile(file_finder.templates.final_tissues_LUT, fn_LUT)

    if create_surfaces:
        # Create surfaces ala CAT12
        logger.info("Starting surface creation")
        starttime = time.time()
        fsavgDir = file_finder.Templates().freesurfer_templates

        surface_settings = settings["surfaces"]
        nprocesses = surface_settings["processes"]
        surf = surface_settings["surf"]
        pial = surface_settings["pial"]
        vdist = surface_settings["vdist"]
        voxsize_pbt = surface_settings["voxsize_pbt"]
        voxsize_refineCS = surface_settings["voxsize_refinecs"]
        th_initial = surface_settings["th_initial"]
        no_selfintersections = surface_settings["no_selfintersections"]
        fillin_gm_from_surf = surface_settings["fillin_gm_from_surf"]
        open_sulci_from_surf = surface_settings["open_sulci_from_surf"]
        exclusion_tissues_fillinGM = surface_settings["exclusion_tissues_fillin_gm"]
        exclusion_tissues_openCSF = surface_settings["exclusion_tissues_open_csf"]

        multithreading_script = [
            os.path.join(SIMNIBSDIR, "segmentation", "run_cat_multiprocessing.py")
        ]
        # fmt: off
        argslist = (
            [
                "--Ymf", sub_files.norm_image,
                "--Yleft_path", sub_files.hemi_mask,
                "--Ymaskhemis_path", sub_files.cereb_mask,
                "--surface_folder", sub_files.surface_folder,
                "--fsavgdir", fsavgDir,
                "--surf",
            ]
            + surf
            + ["--pial"]
            + pial
            + [
                "--vdist", str(vdist[0]), str(vdist[1]),
                "--voxsize_pbt", str(voxsize_pbt[0]), str(voxsize_pbt[1]),
                "--voxsizeCS", str(voxsize_refineCS[0]), str(voxsize_refineCS[1]),
                "--th_initial", str(th_initial),
                "--no_intersect", str(no_selfintersections),
                "--nprocesses", str(nprocesses),
            ]
        )
        if debug:
            argslist += ["--debug"]
        # fmt: on

        proc = subprocess.run(
            [sys.executable] + multithreading_script + argslist, stderr=subprocess.PIPE
        )  # stderr: standard stream for simnibs logger
        logger.debug(proc.stderr.decode("ASCII", errors="ignore").replace("\r", ""))
        proc.check_returncode()
        elapsed = time.time() - starttime
        logger.info("Total time surface creation (HH:MM:SS):")
        logger.info(time.strftime("%H:%M:%S", time.gmtime(elapsed)))

        sub_files = file_finder.SubjectFiles(subpath=subject_dir)
        if fillin_gm_from_surf or open_sulci_from_surf:
            logger.info("Improving GM from surfaces")
            starttime = time.time()
            # original tissue mask used for mesh generation
            if debug:
                shutil.copyfile(
                    sub_files.tissue_labeling_upsampled,
                    os.path.join(sub_files.label_prep_folder, "before_surfmorpho.nii.gz"),
                )

            label_nii = nib.load(sub_files.tissue_labeling_upsampled)
            label_img = np.asanyarray(label_nii.dataobj)
            label_affine = label_nii.affine

            # orginal label mask
            label_nii = nib.load(sub_files.labeling)
            labelorg_img = np.asanyarray(label_nii.dataobj)
            labelorg_affine = label_nii.affine

            if fillin_gm_from_surf:
                # GM central surfaces
                m = Msh()
                if "lh" in surf:
                    m = m.join_mesh(read_gifti_surface(sub_files.get_surface("lh", "central")))
                if "rh" in surf:
                    m = m.join_mesh(read_gifti_surface(sub_files.get_surface("rh", "central")))
                    # fill in GM and save updated mask
                if m.nodes.nr > 0:
                    label_img = charm_utils._fillin_gm_layer(
                        label_img,
                        label_affine,
                        labelorg_img,
                        labelorg_affine,
                        m,
                        exclusion_tissues=exclusion_tissues_fillinGM,
                    )

                    label_nii = nib.Nifti1Image(label_img, label_affine)
                    nib.save(label_nii, sub_files.tissue_labeling_upsampled)
                else:
                    logger.warning(
                        "Neither lh nor rh reconstructed. Filling in from GM surface skipped"
                    )

            if open_sulci_from_surf:
                # GM pial surfaces
                m = Msh()
                if "lh" in pial:
                    m2 = read_gifti_surface(
                        sub_files.get_surface("lh", "pial")
                    )
                    # remove self-intersections using meshfix
                    with tempfile.NamedTemporaryFile(suffix=".off") as f:
                        mesh_fn = f.name
                    write_off(m2, mesh_fn)
                    cmd = [file_finder.path2bin("meshfix"), mesh_fn, "-o", mesh_fn]
                    spawn_process(cmd, lvl=logging.DEBUG)
                    m = m.join_mesh(read_off(mesh_fn))
                    if os.path.isfile(mesh_fn):
                        os.remove(mesh_fn)
                    if os.path.isfile("meshfix_log.txt"):
                        os.remove("meshfix_log.txt")
                if "rh" in pial:
                    m2 = read_gifti_surface(
                        sub_files.get_surface("rh", "pial")
                    )
                    # remove self-intersections using meshfix
                    with tempfile.NamedTemporaryFile(suffix=".off") as f:
                        mesh_fn = f.name
                    write_off(m2, mesh_fn)
                    cmd = [file_finder.path2bin("meshfix"), mesh_fn, "-o", mesh_fn]
                    spawn_process(cmd, lvl=logging.DEBUG)
                    m = m.join_mesh(read_off(mesh_fn))
                    if os.path.isfile(mesh_fn):
                        os.remove(mesh_fn)
                    if os.path.isfile("meshfix_log.txt"):
                        os.remove("meshfix_log.txt")
                if m.nodes.nr > 0:
                    charm_utils._open_sulci(
                        label_img,
                        label_affine,
                        labelorg_img,
                        labelorg_affine,
                        m,
                        exclusion_tissues=exclusion_tissues_openCSF,
                    )
                    label_nii = nib.Nifti1Image(label_img, label_affine)
                    nib.save(label_nii, sub_files.tissue_labeling_upsampled)
                else:
                    logger.warning(
                        "Neither lh nor rh pial reconstructed. Opening up of sulci skipped."
                    )

            # print time duration
            elapsed = time.time() - starttime
            logger.info("Total time cost for GM imrpovements in (HH:MM:SS):")
            logger.info(time.strftime("%H:%M:%S", time.gmtime(elapsed)))

    if mesh_image:
        # create mesh from label image
        logger.info("Starting mesh")
        label_image = nib.load(sub_files.tissue_labeling_upsampled)
        label_buffer = np.round(label_image.get_fdata()).astype(
            np.uint16
        )  # Cast to uint16, otherwise meshing complains
        label_affine = label_image.affine
        label_buffer, label_affine, _ = crop_vol(
            label_buffer, label_affine, label_buffer > 0, thickness_boundary=5
        )
        # reduce memory consumption a bit

        # Read in settings for meshing
        mesh_settings = settings["mesh"]
        elem_sizes = mesh_settings["elem_sizes"]
        smooth_size_field = mesh_settings["smooth_size_field"]
        skin_facet_size = mesh_settings["skin_facet_size"]
        if not skin_facet_size:
            skin_facet_size = None
        facet_distances = mesh_settings["facet_distances"]
        optimize = mesh_settings["optimize"]
        remove_spikes = mesh_settings["remove_spikes"]
        skin_tag = mesh_settings["skin_tag"]
        if not skin_tag:
            skin_tag = None
        hierarchy = mesh_settings["hierarchy"]
        if not hierarchy:
            hierarchy = None
        smooth_steps = mesh_settings["smooth_steps"]
        skin_care = mesh_settings["skin_care"]
        mmg_noinsert = mesh_settings["mmg_noinsert"]

        # Meshing

        debug_path = None
        if debug:
            debug_path = sub_files.subpath

        # if num_threads is zero or less
        # set it to something fairly large
        if num_threads <= 0:
            num_threads = 32

        final_mesh = create_mesh(
            label_buffer,
            label_affine,
            elem_sizes=elem_sizes,
            smooth_size_field=smooth_size_field,
            skin_facet_size=skin_facet_size,
            facet_distances=facet_distances,
            optimize=optimize,
            remove_spikes=remove_spikes,
            skin_tag=skin_tag,
            hierarchy=hierarchy,
            smooth_steps=smooth_steps,
            skin_care=skin_care,
            num_threads=num_threads,
            mmg_noinsert=mmg_noinsert,
            debug_path=debug_path,
            debug=debug
        )

        logger.info("Writing mesh")
        write_msh(final_mesh, sub_files.fnamehead)
        v = final_mesh.view(cond_list=cond_utils.standard_cond())
        v.write_opt(sub_files.fnamehead)

        logger.info("Transforming EEG positions")
        idx = (final_mesh.elm.elm_type == 2) & (final_mesh.elm.tag1 == skin_tag)
        mesh = final_mesh.crop_mesh(elements=final_mesh.elm.elm_number[idx])

        if not os.path.exists(sub_files.eeg_cap_folder):
            os.mkdir(sub_files.eeg_cap_folder)

        cap_files = glob.glob(os.path.join(file_finder.ElectrodeCaps_MNI, "*.csv"))
        for fn in cap_files:
            fn_out = os.path.splitext(os.path.basename(fn))[0]
            fn_out = os.path.join(sub_files.eeg_cap_folder, fn_out)
            transformations.warp_coordinates(
                fn,
                sub_files.subpath,
                transformation_direction="mni2subject",
                out_name=fn_out + ".csv",
                out_geo=fn_out + ".geo",
                mesh_in=mesh,
            )

        logger.info("Write label image from mesh")
        MNI_template = file_finder.Templates().mni_volume
        mesh = final_mesh.crop_mesh(elm_type=4)
        field = mesh.elm.tag1.astype(np.uint16)
        ed = ElementData(field)
        ed.mesh = mesh
        ed.to_deformed_grid(
            sub_files.mni2conf_nonl,
            MNI_template,
            out=sub_files.final_labels_MNI,
            out_original=sub_files.final_labels,
            method="assign",
            reference_original=sub_files.reference_volume,
        )

        fn_LUT = sub_files.final_labels.rsplit(".", 2)[0] + "_LUT.txt"
        shutil.copyfile(file_finder.templates.final_tissues_LUT, fn_LUT)

    # -------------------------TIDY UP-------------------------------------
    # Create charm_report.html
    logger.info("Creating report")
    html_writer.write_report(sub_files)

    # log stopping time and total duration ...
    logger.info("charm run finished: " + time.asctime())
    logger.info(
        "Total running time: " + utils.simnibs_logger.format_time(time.time() - start)
    )

    # stop logging ...
    _stop_logger(sub_files.charm_log)


def _setup_logger(logfile):
    """Add FileHandler etc."""
    with open(logfile, "a") as f:
        f.write("<HTML><HEAD><TITLE>charm report</TITLE></HEAD><BODY><pre>")
        f.close()
    fh = logging.FileHandler(logfile, mode="a")
    formatter = logging.Formatter("%(levelname)s: %(message)s")
    fh.setFormatter(formatter)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    utils.simnibs_logger.register_excepthook(logger)


def _stop_logger(logfile):
    """Close down logging"""
    while logger.hasHandlers():
        logger.removeHandler(logger.handlers[0])
    utils.simnibs_logger.unregister_excepthook()
    logging.shutdown()
    with open(logfile, "r") as f:
        logtext = f.read()

    # Explicitly remove this really annoying stuff from the log
    removetext = (
        re.escape("-\|/"),
        re.escape("Selecting intersections ... ")
        + "\d{1,2}"
        + re.escape(" %Selecting intersections ... ")
        + "\d{1,2}"
        + re.escape(" %"),
    )
    with open(logfile, "w") as f:
        for text in removetext:
            logtext = re.sub(text, "", logtext)
        f.write(logtext)
        f.write("</pre></BODY></HTML>")
        f.close()


def _read_settings_and_copy(usesettings, fn_settingslog):
    # read settings and copy settings file
    if usesettings is None:
        fn_settings = os.path.join(SIMNIBSDIR, "charm.ini")
    else:
        if type(usesettings) == list:
            usesettings = usesettings[0]
        fn_settings = usesettings
    settings = utils.settings_reader.read_ini(fn_settings)
    try:
        shutil.copyfile(fn_settings, fn_settingslog)
    except shutil.SameFileError:
        pass
    return settings


def _prepare_t1(T1, reference_volume, force_qform, force_sform):
    # copy T1 (as nii.gz) if supplied
    if T1:
        if os.path.exists(T1):
            # Cast to float32 and save
            T1_tmp = nib.load(T1)
            T1_tmp = _check_q_and_s_form(T1_tmp, force_qform, force_sform)

            # Check for singleton dimensions and squeeze those out
            # Note: do this after the s-g-form check so affine is correct
            if (np.array(T1_tmp.shape) == 1).any():
                data_tmp = np.squeeze(T1_tmp.get_fdata())
                T1_tmp = nib.Nifti1Image(data_tmp, T1_tmp.affine)

            T1_tmp.set_data_dtype(np.float32)

            nib.save(T1_tmp, reference_volume)
        else:
            raise FileNotFoundError(f"Could not find input T1 file: {T1}")


def _prepare_t2(T1, T2, registerT2, T2_reg, force_qform, force_sform):
    if T2:
        T2_exists = os.path.exists(T2)
        if registerT2:
            if T2_exists:
                # Abuse the T2_reg filename to save a temporary
                # file in case the q and sform need to be fixed
                T2_tmp = nib.load(T2)
                T2_tmp = _check_q_and_s_form(T2_tmp, force_qform, force_sform)

                if (np.array(T2_tmp.shape) == 1).any():
                    data_tmp = np.squeeze(T2_tmp.get_fdata())
                    T2_tmp = nib.Nifti1Image(data_tmp, T2_tmp.affine)

                nib.save(T2_tmp, T2_reg)
                charm_utils._registerT1T2(T1, T2_reg, T2_reg)
            else:
                raise FileNotFoundError(f"Could not find input T2 file: {T2}")
        else:
            if T2_exists:
                T2_tmp = nib.load(T2)
                T2_tmp = _check_q_and_s_form(T2_tmp, force_qform, force_sform)

                if (np.array(T2_tmp.shape) == 1).any():
                    data_tmp = np.squeeze(T2_tmp.get_fdata())
                    T2_tmp = nib.Nifti1Image(data_tmp, T2_tmp.affine)

                T2_tmp.set_data_dtype(np.float32)
                nib.save(T2_tmp, T2_reg)


def _check_q_and_s_form(scan, force_qform=False, force_sform=False):
    if not np.array_equal(scan.get_qform(), scan.get_sform()):
        if not (force_qform or force_sform):
            raise ValueError("The qform and sform do not match. Please run charm with the --forceqform (preferred) or --forcesform option")
        elif force_qform:
            scan.set_sform(scan.get_qform())
        elif force_sform:
            # Note set_qform will strip shears silently per nibabel documentation
            scan.set_qform(scan.get_sform())
            scan.set_sform(scan.get_qform())

    return scan


def _setup_atlas(samseg_settings, T2_reg, usesettings):
    # Set-up atlas paths
    atlas_name = samseg_settings["atlas_name"]
    logger.info("Using " + atlas_name + " as charm atlas.")
    atlas_path = os.path.join(file_finder.templates.charm_atlas_path, atlas_name)
    atlas_settings = utils.settings_reader.read_ini(
        os.path.join(atlas_path, atlas_name + ".ini")
    )
    atlas_settings_names = atlas_settings["names"]
    template_name = os.path.join(atlas_path, atlas_settings_names["template_name"])
    atlas_affine_name = os.path.join(atlas_path, atlas_settings_names["affine_atlas"])
    atlas_level1 = os.path.join(atlas_path, atlas_settings_names["atlas_level1"])
    atlas_level2 = os.path.join(atlas_path, atlas_settings_names["atlas_level2"])
    custom_gmm_parameters = samseg_settings["gmm_parameter_file"]

    if type(usesettings) == list:
        usesettings = usesettings[0]

    if not usesettings or not custom_gmm_parameters:
        if os.path.exists(T2_reg):
            gmm_parameters = os.path.join(
                atlas_path, atlas_settings_names["gaussian_parameters_t2"]
            )
        else:
            gmm_parameters = os.path.join(
                atlas_path, atlas_settings_names["gaussian_parameters_t1"]
            )
    else:
        settings_dir = os.path.dirname(usesettings)
        gmm_parameters = os.path.join(settings_dir, custom_gmm_parameters)
        if not os.path.exists(gmm_parameters):
            raise FileNotFoundError(f"Could not find gmm parameter file: {gmm_parameters}")

    return (
        template_name,
        atlas_settings,
        atlas_path,
        atlas_level1,
        atlas_level2,
        atlas_affine_name,
        gmm_parameters,
    )


def _denoise_inputs(T1, T2, sub_files):
    if T1:
        logger.info("Denoising the T1 input and saving.")
        charm_utils._denoise_input_and_save(T1, sub_files.T1_denoised)
    if T2 and os.path.exists(sub_files.T2_reg):
        logger.info("Denoising the registered T2 and saving.")
        charm_utils._denoise_input_and_save(sub_files.T2_reg, sub_files.T2_reg_denoised)

def _read_transform(transform_file):
    transform = np.loadtxt(transform_file)
    assert transform.shape == (
        4,
        4,
    ), f"`transform` should have shape (4, 4), got {transform.shape}"
    # Change from RAS to LPS. ITK uses LPS internally
    RAS2LPS = np.diag([-1, -1, 1, 1])
    return RAS2LPS @ transform @ RAS2LPS

from pathlib import Path
from typing import Union

import mne
from mne.io.constants import FIFF
import numpy as np
import scipy.sparse

from simnibs.mesh_tools import mesh_io
from simnibs.utils.file_finder import SubjectFiles
from simnibs.utils.simnibs_logger import logger
from simnibs.simulation import eeg

# mne.set_log_level("warning")

# SOURCE SPACE

HEMISPHERES_MNE = dict(
    lh=FIFF.FIFFV_MNE_SURF_LEFT_HEMI, rh=FIFF.FIFFV_MNE_SURF_RIGHT_HEMI
)
MNE_LANDMARKS = ("nasion", "lpa", "rpa")
landmarks_mapper = {
    "simnibs_to_mne": {"Nz": "nasion", "Iz": "inion", "LPA": "lpa", "RPA": "rpa"},
    "mne_to_simnibs": {"nasion": "Nz", "inion": "Iz", "lpa": "LPA", "rpa": "RPA"},
}


def _get_src_sphere(src, subsampling=None):
    if isinstance(src, (Path, str)):
        sf = SubjectFiles(subpath=str(src))
        src = eeg.load_surface(sf, "sphere_reg", subsampling)
    assert isinstance(src, dict)
    return src


def make_source_morph(
    src_from, src_to, src_from_sphere, src_to_sphere, subsampling=None
):
    """
    The spherical registrations are used to create the morphing matrices for
    each hemisphere.
    The source spaces are used to determine which sources are valid (in use).
    The morphing matrices will be smoothed so as to compensate for any unused
    source positions such that all valid sources in `src_to` will be covered.

    PARAMETERS
    ----------
    src_from : mne.SourceSpaces
        Contains information about which positions are used.
    src_to : mne.SourceSpaces
        Contains information about which positions are used.
    src_from_sphere : Path | str | dict
        Path to m2m directory or a dict with `points` and `tris` describing
        the spherical registration surface.
    src_to_sphere : Path | str | dict
        See `src_from_sphere`.

    RETURNS
    -------
    SourceMorph object.
    """
    src_from_sphere = _get_src_sphere(src_from_sphere, subsampling)
    src_to_sphere = _get_src_sphere(src_to_sphere, subsampling)
    mmaps = eeg.make_morph_maps(src_from_sphere, src_to_sphere)
    return _make_source_morph(src_from, src_to, mmaps)


def setup_source_space(
    m2m_dir: Union[Path, str],
    subsampling: Union[None, int] = None,
    morph_to_fsaverage: Union[None, int] = 10
):
    """Setup a source space for use with MNE-Python.

    PARAMETERS
    ----------
    m2m_dir : Path-like
        The `m2m` directory containing the segmentation and head model results.
    subsampling : None | int
        The subsampling to use (default = None).
    morph_to_fsaverage : None | int
        Whether or not to create a mapping from subject space to the fsaverage
        template (default = 10, which constructs a morph to the fsaverage 10k
        model).

    RETURNS
    -------
    src_from : mne.SourceSpaces
        The source space object.
    morph : mne.SourceMorph | None
        The source morph object.
    """
    subjectfiles = SubjectFiles(subpath=str(m2m_dir))
    # src_from, src_from_sphere = eeg.load_subject_surfaces(subjectfiles, subsampling)
    src_from = eeg.load_surface(subjectfiles, "central", subsampling)
    src_from = make_source_spaces(src_from, subjectfiles.subid)

    # Make source morph
    if morph_to_fsaverage:
        fsavg = eeg.FsAverage(morph_to_fsaverage)
        fsavg_sphere = fsavg.get_surface("sphere")
        src_fsavg_sphere = make_source_spaces(fsavg_sphere, fsavg.name)

        morph = make_source_morph(
            src_from, src_fsavg_sphere, m2m_dir, fsavg_sphere, subsampling
        )

        # mmaps = eeg.make_morph_maps(src_from_sphere, fsavg_sphere)
        # src_from_sphere = make_source_spaces(src_from_sphere, subjectfiles.subid)
        # fsavg_sphere = make_source_spaces(fsavg_sphere, fsavg.name)
        # morph = _make_source_morph(src_from_sphere, fsavg_sphere, mmaps)
    else:
        morph = None
    return src_from, morph


def make_source_spaces(surf: dict, subject_id: Union[None,str] = None, coord_frame: str = "mri"):
    """Create an MNE-Python source space object.

    PARAMETERS
    ----------
    surf : dict
        Dictionary with entries `lh` and `rh` each being a dictionary with keys
        `points`, `tris`, and possibly `normals`.
    subject_id : str
        The id of the subject.

    RETURNS
    -------
    src : mne.SourceSpaces
        The source space object.
    """
    return mne.SourceSpaces(
        [
            make_source_space(
                **v, coord_frame=coord_frame, surf_id=k, subject_id=subject_id
            )
            for k, v in surf.items()
        ]
    )


def make_source_space(
    points: np.ndarray,
    tris: Union[None,np.ndarray] = None,
    normals: Union[None,np.ndarray] = None,
    coord_frame: str = "mri",
    surf_id: Union[None,str] = None,
    subject_id: Union[None,str] = None,
):
    """Setup a dictionary of a single hemisphere for an MNE-Python source space
    object.

    PARAMETERS
    ----------
    points : ndarray
        Source positions (in mm).
    tris : ndarray
        The triangulation of the source positions if relevant (default = None).
    normals : ndarray
        Source normals (default = None).
    coord_frame : str
        The coordinate frame of the source positions. Choose `mri` or `head`
        default = 'mri')
    surf_id : str
        Surface id, i.e., `lh`, `rh`, or None (default = None).
    subject_id = str
        Subject id (default = None).

    RETURNS
    -------
    src : dict
        The source space.
    """

    # mm -> m (input assumed to be in mm)
    rr = points * 1e-3
    n_points = rr.shape[0]

    # Just use arbitrary directions
    dummy_nn = np.zeros((n_points, 3))
    dummy_nn[:, 2] = 1.0

    if coord_frame == "mri":
        coord_frame = FIFF.FIFFV_COORD_MRI
    elif coord_frame == "head":
        coord_frame = FIFF.FIFFV_COORD_HEAD
    else:
        raise ValueError(f"coord_frame must be `mri` or `head`, got {coord_frame}")

    assert surf_id in set(("lh", "rh", None))
    if surf_id == "lh":
        surf_id = FIFF.FIFFV_MNE_SURF_LEFT_HEMI
    elif surf_id == "rh":
        surf_id = FIFF.FIFFV_MNE_SURF_RIGHT_HEMI
    elif surf_id is None:
        surf_id = FIFF.FIFFV_MNE_SURF_UNKNOWN

    # Setup as discrete source space initially
    src = dict(
        id=surf_id,
        type="discrete",
        np=n_points,
        ntri=0,
        coord_frame=coord_frame,
        rr=rr,
        nn=dummy_nn,
        tris=None,
        inuse=np.ones(n_points, dtype=int),
        nuse=n_points,
        vertno=np.arange(n_points, dtype=int),
        use_tris=None,
        nuse_tri=0,
        subject_his_id=subject_id,
    )
    src.update(
        dict(
            nearest=None,
            nearest_dist=None,
            pinfo=None,
            patch_inds=None,
            dist=None,
            dist_limit=None,
        )
    )

    # Setup as surface source space (MNE does not like surface source spaces
    # which are not lh or rh)
    if tris is not None:
        _add_surface_info(src, tris, normals)

    return src


def _add_surface_info(src, tris, nn):
    assert src["id"] in set(
        (FIFF.FIFFV_MNE_SURF_LEFT_HEMI, FIFF.FIFFV_MNE_SURF_RIGHT_HEMI)
    )
    src["type"] = "surf"
    src["tris"] = tris
    src["ntri"] = tris.shape[0]
    mne.source_space._complete_source_space_info(src)
    if nn is None:
        src["nn"] = (
            mesh_io.Msh(mesh_io.Nodes(src["rr"]), mesh_io.Elements(src["tris"] + 1))
            .nodes_normals()
            .value
        )
    else:
        src["nn"] = eeg.normalize_vectors(nn)


# SOURCE SPACE MORPHING


def _make_source_morph(
    src_from: mne.SourceSpaces,
    src_to: mne.SourceSpaces,
    mmaps: dict,
    smooth: Union[None, int] = None,
    warn: bool = True,
):
    """Make an MNE-Python source morph object.

    PARAMETERS
    ----------
    src_from : mne.SourceSpaces
        The source space to map from.
    src_to : mne.SourceSpaces
        The source space to map to.
    mmaps : dict
        The morph maps describing the morph from `src_from` to `src_to`. A dict
        of length two corresponding to left and right hemispheres should be
        provided.
    smooth : int | None
        Whether or not to smooth the morph maps. If None, will smooth until all

        (default = None)
    warn : bool
        Set warning flag.

    RETURNS
    -------
    morph : mne.SourceMorph
        The source morph object.
    """
    morph_mat = compute_morph_matrix(src_from, src_to, mmaps, smooth, warn)

    # Information required for SourceMorph
    src_data, kind, subject_from = mne.morph._get_src_data(src_from)
    subject_to = src_to[0]["subject_his_id"]
    vertices_to = [s["vertno"] for s in src_to]

    # Some defaults from mne.compute_source_morph (unused)
    zooms = "auto"
    niter_affine = (100, 100, 10)
    niter_sdr = (5, 5, 3)
    xhemi = False
    spacing = shape = affine = pre_affine = sdr_morph = None

    # fmt: off
    return mne.SourceMorph(subject_from, subject_to, kind, zooms,
                           niter_affine, niter_sdr, spacing, smooth, xhemi,
                           morph_mat, vertices_to, shape, affine,
                           pre_affine, sdr_morph, src_data, None)
    # fmt: on


def compute_morph_matrix(
    src_from: mne.SourceSpaces,
    src_to: mne.SourceSpaces,
    mmaps: dict,
    smooth: Union[None, int] = None,
    warn: bool = True,
):
    """Collect the mmaps into one sparse block matrix and apply smoothing if
    desired/needed.

    PARAMETERS
    ----------
    src_from : mne.SourceSpaces
        The source space to map from.
    src_to : mne.SourceSpaces
        The source space to map to.
    mmaps : dict of scipy.sparse.csr_matrix
        The morph maps describing the morph from `src_from` to `src_to`.
    smooth : int | None
    warn : bool
        Set warning flag.

    RETURNS
    -------
    morpher : scipy.sparse.csr_matrix

    NOTES
    -----
    This is equivalent to mne.morph._compute_morph_matrix albeit with some
    simplifications as we do not need the 'xhemi' functionality.
    """
    assert (
        len(src_from) == len(src_to) == len(mmaps)
    ), "The length of the inputs must be compatible"

    # iterate over to / block-rows of CSR matrix
    # Ensure we have the same hemispheres from src and mmaps
    morpher = []
    for i, (hemi, hemi_mne) in enumerate(HEMISPHERES_MNE.items()):
        assert src_from[i]["id"] == hemi_mne
        assert src_to[i]["id"] == hemi_mne
        morpher.append(
            mne.morph._hemi_morph(
                src_from[i]["tris"],
                src_to[i]["vertno"],
                src_from[i]["vertno"],
                smooth,
                mmaps[hemi],
                warn,
            )
        )
    return scipy.sparse.block_diag(morpher).tocsr()


# FORWARD


def make_forward(
    forward: dict,
    src: Union[Path, str, mne.SourceSpaces],
    info: Union[Path, str, mne.Info],
    trans: Union[Path, str, mne.Transform],
):
    """Create an MNE forward object from a dictionary with forward information.
    Please note that the forward solution is converted to head coordinate
    system (using `trans`).

    PARAMETERS
    ----------
    forward : dict
        Dictionary with forward solution information (as returned by
        `prepare_leadfield_for_eeg`).
    src : Path | str | mne.SourceSpaces
        Path to source space or SourceSpace object.
    info : Path | str | mne.Info
        Path to epochs, evoked, info, or raw file containing the correct info
        or an Info object.
    trans : Path | str | mne.Transform
        Path to transform or Transform object between head and MRI coordinte
        systems.

    RETURNS
    -------
    fwd : mne.Forward
        The forward object.

    NOTES
    -----
    Bad channels should *not* be included in the leadfield simulations. Should
    any such channels (as marked in the Info structure) be found in the forward
    solution, a warning will be issued.
    """
    kwargs = {
        "bem": None,
        "bem_extra": "",
        "mindist": 0,
        "n_jobs": 1,
        "meg": False,
        "eeg": True,
        "ignore_ref": True,
        "allow_bem_none": True,
    }

    if not isinstance(src, mne.SourceSpaces):
        src = mne.read_source_spaces(src)
    kwargs["src"] = src
    if not isinstance(info, mne.Info):
        kwargs["info_extra"] = str(info)
        info = mne.io.read_info(info)
    else:
        kwargs["info_extra"] = "instance of Info"
    kwargs["info"] = info
    if not isinstance(trans, mne.Transform):
        kwargs["trans"] = str(trans)
        trans = mne.read_trans(trans)
        trans = mne.transforms._ensure_trans(trans, "mri", "head")
    else:
        kwargs["trans"] = "instance of Transform"
    kwargs["mri_head_t"] = trans

    # Check the inputs
    assert all(
        s["coord_frame"] == FIFF.FIFFV_COORD_MRI for s in src
    ), f"Source space should be in {FIFF.FIFFV_COORD_MRI} coordinate system"

    # Pick the EEG channels and check
    # picks = [info["ch_names"][i] for i in mne.pick_types(info, eeg=True)]
    # info.pick_channels(picks)
    assert (
        info["ch_names"] == forward["ch_names"]
    ), f"Inconsistencies between channels in Info and leadfield {info['ch_names']} and {forward['ch_names']}"
    bads_in_forward = tuple(filter(lambda x: x in info["bads"], forward["ch_names"]))
    if any(bads_in_forward):
        logger.warning(
            "Bad channels from Info found in forward solution: " f"{bads_in_forward}"
        )

    # fmt: off
    sensors, _, _, update_kwargs, _ = mne.forward._make_forward._prepare_for_forward(**kwargs)
    # fmt: on

    # update_kwargs["info"]["command_line"] = ""

    # A note on forward solutions and coordinate systems. Below we use
    #
    #     Q[mri]      : dipole moments per position in MRI coordinates
    #     Q[head]     : dipole moments per position in HEAD coordinates
    #     I           : 3x3 identity matrix
    #     mri_head_t  : 3x3 transformation from MRI to HEAD coordinates
    #     head_mri_t  : 3x3 transformation from HEAD to MRI coordinates
    #     (i.e., without translation)
    #
    # MNE forward solutions are in head coordinates, however, they are
    # actually computed in MRI coordinates by transforming electrodes to MRI
    # coords and, instead of computing the solution along the x/y/z axes such
    # that Q[mri] = I, they use Q[mri] = head_mri_t.T as the dipole moments
    # (see mne.forward._compute_forward._prep_field_computation).
    # In this way, the solution corresponds to dipole moments Q[head] = I
    # (i.e. along x/y/z in head coordinate system).
    #
    # SimNIBS forward solutions are in the coordinate system of the volume
    # conductor model, which, in general, will be MRI coordinates since this
    # is the space of the original MRI image(s). Solutions are computed in
    # this coordinate system, i.e., using Q[mri] = I. Thus, we need to
    # transform the SimNIBS solution from Q[mri] = I to Q[head] = I.
    #
    # Since
    #
    #     head_mri_t.T @ mri_head_t.T = I
    #
    # we can convert a forward solution with Q[mri] = I to one with
    # Q[head] = I like
    #
    #     Q[head] = Q[mri] @ mri_head_t.T

    # Convert forward solution from MRI coordinates to head coordinates
    fwd = forward["data"] @ trans["trans"][:3, :3].T
    fwd = fwd.reshape(
        (forward["n_channels"], forward["n_sources"] * forward["n_orientations"])
    )
    # Transpose to (3*n_dipoles, n_sensors)
    fwd = mne.forward._make_forward._to_forward_dict(fwd.T, sensors["eeg"]["ch_names"])
    fwd.update(**update_kwargs)

    return fwd


def prepare_montage(
    fname: Union[Path, str],
    info: Union[Path, str, mne.Info],
    trans: Union[Path, str, mne.Transform],
):
    """Prepare a montage file from an MNE-Python objects. The electrode
    positions in `info` are transformed from head to MRI coordinate system.

    PARAMETERS
    ----------
    fname : str |
        The name of the file to write. If it does not end with '.csv', this
        will be added.
    info : str | mne.Info
        Path to epochs, evoked, info, or raw file containing the correct info
        or an Info object.
    trans : str | mne.Transform
        Path to transform or Transform object between head and MRI coordinte
        systems.

    RETURNS
    -------
    None
    """
    fname = Path(fname)
    if fname.suffix != ".csv":
        fname = fname.with_suffix(".csv")
    if not isinstance(info, mne.Info):
        info = mne.io.read_info(info)
    if not isinstance(trans, mne.Transform):
        trans = mne.read_trans(trans)

    # Electrodes are in head coordinates and m so transform to MRI coordinates
    # and mm
    trans = mne.transforms._ensure_trans(trans, "head", "mri")
    electrodes = {
        ch["ch_name"]: 1e3 * mne.transforms.apply_trans(trans, ch["loc"][:3])
        for ch in info["chs"]
    }

    with open(fname, "w") as f:
        # f.write("Type,X,Y,Z,ElectrodeName\n")
        for ch_name, ch_pos in electrodes.items():
            # There should be no space after ","!
            line = ",".join(["Electrode", *map(str, ch_pos), ch_name]) + "\n"
            f.write(line)

    fid_pos = mne.viz._3d._fiducial_coords(info["dig"])
    if fid_pos.size > 0:
        fid_pos = 1e3 * mne.transforms.apply_trans(trans, fid_pos)
        mapper = {
            FIFF.FIFFV_POINT_LPA: "LPA",
            FIFF.FIFFV_POINT_NASION: "Nz",
            FIFF.FIFFV_POINT_RPA: "RPA",
            FIFF.FIFFV_POINT_INION: "Iz",
        }
        fid_names = [mapper[i] for i in mne.viz._3d.FIDUCIAL_ORDER]
        with open(fname, "a") as f:
            for n, p in zip(fid_names, fid_pos):
                line = ",".join(["Fiducial", *map(str, p), n]) + "\n"
                f.write(line)


def simnibs_montage_to_mne_montage(montage, coord_frame="unknown"):
    valid_entries = {"Electrode", "ReferenceElectrode"}
    if len(montage.ch_names) > 0 and len(montage.ch_pos) > 0:
        ch_pos = {
            n: p
            for n, p, t in zip(
                montage.ch_names, 1e-3 * montage.ch_pos, montage.ch_types
            )
            if t in valid_entries
        }
    else:
        ch_pos = None

    landmarks = {k: None for k in MNE_LANDMARKS}
    if montage.landmarks:
        for k in landmarks:
            landmarks[k] = (
                1e-3 * montage.landmarks[landmarks_mapper["mne_to_simnibs"][k]]
            )

    return mne.channels.make_dig_montage(ch_pos, coord_frame=coord_frame, **landmarks)


def mne_montage_to_simnibs_montage(montage, name=None):
    d = montage.get_positions()
    ch_names = list(d["ch_pos"].keys())
    ch_pos = [p*1e3 for p in d["ch_pos"].values()]
    ch_types = ["Electrode"] * len(ch_names)
    try:
        landmarks = {
            landmarks_mapper["mne_to_simnibs"][k]: 1e3 * d[k] for k in MNE_LANDMARKS
        }
    except TypeError:
        # at least one of nasion, lpa, or rpa is None
        landmarks = None
    return eeg.Montage(name, ch_names, ch_pos, ch_types, landmarks)

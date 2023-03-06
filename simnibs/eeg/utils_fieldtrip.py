from pathlib import Path
from typing import Union

from nibabel.affines import apply_affine
import numpy as np
from scipy.io import loadmat

import simnibs
from simnibs.simulation import eeg

# EEG MONTAGE

to_m = dict(m=1, cm=1e-2, mm=1e-3)
from_m = {k: 1 / v for k, v in to_m.items()}


def prepare_montage(
    fname_montage: Union[Path, str],
    fname_info: Union[Path, str],
    fname_trans: Union[None, Path, str] = None,
):
    """Prepare SimNIBS montage file from a FieldTrip data structure containing
    and `elec` field describing the electrode configuration.

    PARAMETERS
    ----------
    fname_montage :
        Name of the SimNIBS montage file to write.
    fname_info :
        Name of a MAT file containing the electrode information in the field
        `elec`.
    fname_trans :
        Affine transformation matrix to apply to the electrode positions before
        writing to `fname_montage`.

    RETURNS
    -------
    """
    fname_info = Path(fname_info)

    info = loadmat(fname_info)["elec"][0]
    info = dict(zip(info.dtype.names, info[0]))
    info["label"] = np.array([label[0] for label in info["label"].squeeze()])
    info["chantype"] = np.array([ct[0] for ct in info["chantype"].squeeze()])
    info["unit"] = info["unit"][0]
    scale = to_m[info["unit"]] * from_m["mm"]
    info["elecpos"] *= scale
    info["chanpos"] *= scale  # unused

    if fname_trans:
        fname_trans = Path(fname_trans)
        if fname_trans.suffix == ".mat":
            trans = loadmat(fname_trans)["trans"]
        elif fname_trans.suffix == ".txt":
            trans = np.loadtxt(fname_trans)
        else:
            raise ValueError("`fname_trans` must be either a MAT or a TXT file.")

        info["elecpos"] = apply_affine(trans, info["elecpos"])
        info["chanpos"] = apply_affine(trans, info["chanpos"])

    is_eeg = info["chantype"] == "eeg"
    montage = np.column_stack(
        (
            np.broadcast_to(np.array(["Electrode"]), is_eeg.sum()),
            info["elecpos"][is_eeg],
            info["label"][is_eeg],
        )
    )
    np.savetxt(fname_montage, montage, fmt="%s", delimiter=",")


# SOURCE SPACE


def setup_source_space(
    m2m_dir: Union[Path, str],
    subsampling: Union[None, int] = None,
    morph_to_fsaverage: Union[None, int] = 10,
):
    """Setup a source space for use with FieldTrip.

    PARAMETERS
    ----------
    m2m_dir : Path-like
        The directory containing the segmentation and head model results.
    subsampling : None | int
        The subsampling to use (default = None).
    morph_to_fsaverage : None | int
        Whether or not to create a mapping from subject space to the fsaverage
        template (default = 10, which constructs a morph to the fsaverage 10k
        model).

    RETURNS
    -------
    src_from : dict
        Dictionary with the source model information.
    mmaps : dict | None
        Dictionary with scipy.sparse.csr_matrix describing the morph from
        subject to fsaverage.
    """

    subjectfiles = simnibs.utils.file_finder.SubjectFiles(subpath=str(m2m_dir))
    src_from, src_from_reg = eeg.load_subject_surfaces(subjectfiles, subsampling)
    src_from = make_sourcemodel(src_from)

    if morph_to_fsaverage:
        fsavg = eeg.FsAverage(morph_to_fsaverage)
        fsavg_sphere = fsavg.get_surface("sphere")

        mmaps = eeg.make_morph_maps(src_from_reg, fsavg_sphere)
    else:
        mmaps = None
    return src_from, mmaps


def make_sourcemodel(src: dict):
    """Create a dictionary with source model information which can be used with
    FieldTrip.

    PARAMETERS
    ----------
    src : dict
        Dictionary with entries `lh` and `rh` each being a dictionary with keys
        `points`, `tris`, and possibly `normals`.

    RETURNS
    -------
    sourcemodel : dict
        Dictionary with source model information.
    """
    hemis = ("lh", "rh")
    hemi2fieldtrip = dict(lh="CORTEX_LEFT", rh="CORTEX_RIGHT")
    assert all(h in src for h in hemis)

    pos = np.concatenate([src[h]["points"] for h in hemis])

    # Construct a composite triangulation for both hemispheres
    nrs = {h: len(src[h]["points"]) for h in hemis}
    tri_reindex = (0, nrs[hemis[0]])
    tri = np.concatenate([src[h]["tris"] + i for h, i in zip(hemis, tri_reindex)]) + 1

    # All sources are valid
    inside = np.ones(sum(nrs.values()), dtype=bool)[:, None]

    # ft_read_headshape outputs this structure
    brainstructure = np.concatenate(
        [np.full(nrs[h], i) for i, h in enumerate(hemis, start=1)]
    )
    brainstructurelabel = np.stack([hemi2fieldtrip[h] for h in hemis]).astype(object)

    sourcemodel = dict(
        pos=pos,
        tri=tri,
        unit="mm",
        inside=inside,
        brainstructure=brainstructure,
        brainstructurelabel=brainstructurelabel,
    )

    if all("normals" in src[h] for h in hemis):
        sourcemodel["normals"] = np.concatenate([src[h]["normals"] for h in hemis])

    return sourcemodel


# FORWARD


def make_forward(forward: dict, src: dict):
    """Make a forward dictionary in a FieldTrip compatible format from a
    dictionary with forward information.

    PARAMETERS
    ----------
    forward : dict
        Dictionary with forward solution information (as returned by
        `prepare_leadfield_for_eeg`).
    src : dict
        Dictionary with source space information (as returned by
        `setup_source_space`).

    RETURNS
    -------
    fwd : dict
        The forward dictionary.

    NOTES
    -----
    The leadfield of a particular electrode may be plotted in FieldTrip like so

        fwd_mat = cell2mat(fwd.leadfield);
        elec = 10; % electrode
        ori = 2;  % orientation (1,2,3 corresponding to x,y,z)
        ft_plot_mesh(fwd, 'vertexcolor', fwd_mat(elec, ori:3:end)');
    """
    # Create a cell array of matrices by filling a numpy object. Each cell is
    # [n_channels, n_orientations]
    fwd_cell = np.empty(forward["n_sources"], dtype=object)
    for i in range(forward["n_sources"]):
        fwd_cell[i] = forward["data"][:, i]

    labels = np.array(forward["ch_names"]).astype(object)[:, None]

    # The forward structure of FieldTrip
    return dict(
        pos=src["pos"],
        inside=src["inside"],
        unit=src["unit"],
        leadfield=fwd_cell,
        label=labels,
        leadfielddimord=r"{pos}_chan_ori",
        cfg=f"Created by SimNIBS {simnibs.__version__}",
    )

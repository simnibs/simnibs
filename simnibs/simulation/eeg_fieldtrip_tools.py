from pathlib import Path
from typing import Dict, Union

import numpy as np

import simnibs
from simnibs.simulation import eeg


# SOURCE SPACE


def setup_source_space(
    m2m_dir: Union[Path, str], subsampling: int = None, morph_to_fsaverage: int = 5
):
    """Setup a source space for use with FieldTrip.

    PARAMETERS
    ----------
    m2m_dir : Path-like
        The directory containing the segmentation and head model results.
    subsampling : int
        The subsampling to use (default = None).
    morph_to_fsaverage : int
        Whether or not to create a mapping from subject space to the fsaverage
        template (default = 5).

    RETURNS
    -------
    src_from : dictionary
        Dictionary with the source model information.
    mmaps : dictionary
        Dictionary with scipy.sparse.csr_matrix describing the morph from
        subject to fsaverage.
    """

    subjectfiles = simnibs.utils.file_finder.SubjectFiles(subpath=m2m_dir)
    src_from, src_from_reg = eeg.load_subject_surfaces(subjectfiles, subsampling)
    src_from = make_sourcemodel(src_from)

    if morph_to_fsaverage:
        fsavg = eeg.FsAverage(morph_to_fsaverage)
        fsavg_sphere = fsavg.get_surface("sphere")

        mmaps = eeg.make_morph_maps(src_from_reg, fsavg_sphere)
    else:
        mmaps = None
    return src_from, mmaps


def make_sourcemodel(src: Dict):
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
    tri_reindex = (0, nrs[hemis[1]])
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


def make_forward(forward: Dict, src: Dict):
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

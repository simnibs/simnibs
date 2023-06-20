from pathlib import Path
from typing import Union

import h5py
import numpy as np

import simnibs
import simnibs.eeg.utils
from simnibs.segmentation.brain_surface import subsample_surfaces
from simnibs.utils.file_finder import SubjectFiles

EEG_EXPORT_FORMATS = {"mne", "fieldtrip"}


def compute_tdcs_leadfield(
    m2m_dir: Union[Path, str],
    fem_dir: Union[Path, str],
    fname_montage: Union[Path, str],
    subsampling: Union[int, None] = None,
    point_electrodes: bool = True,
    init_kwargs: Union[dict, None] = None,
    run_kwargs: Union[dict, None] = None,
):
    """Convenience function for running a TDCS simulation. The result can be
    converted to an EEG forward solution using `make_forward`.

    PARAMETERS
    ----------
    m2m_dir : str
    fem_dir : str
    fname_montage : str
        The filename of the montage to use.
    subsampling : int | None
        The subsampling to use. If the files do not exist, they are created
        (default = None).
    point_electrodes : bool
        Whether to model electrodes as points (defaults) or mesh them onto the
        headmodel.
    init_kwargs : dict
        Kwargs used to update the initialization of `TDCSLEADFIELD`. Can be
        used to override default settings, e.g., by passing another list of
        conductivities which can be achieve by doing
        `init_kwargs=dict(cond=my_conds)` where my_conds has the same format
        as `simnibs.simulation.cond.standard_cond` (default = None).
    run_kwargs : dict
        Kwargs to pass to the run method of `TDCSLEADFIELD`. If None, will use
        dict(save_mat=False, cpus=1) (default = None).

    RETURNS
    -------
    fname : PathLike
        The name of the hdf5 file containing the simulation results.
    """
    if run_kwargs is None:
        run_kwargs = dict(save_mat=False, cpus=1)
    assert isinstance(run_kwargs, dict)

    m2m = SubjectFiles(subpath=str(m2m_dir))

    if subsampling and not all(m2m.get_surface(h, "central", subsampling=subsampling).exists() for h in m2m.hemispheres):
        _ = subsample_surfaces(m2m_dir, n_points=subsampling)

    # The paths should be strings otherwise errors might occur when writing the
    # .hdf5 file
    lf = simnibs.simulation.sim_struct.TDCSLEADFIELD()
    lf.fnamehead = str(m2m.fnamehead)
    lf.subpath = str(m2m.subpath)
    lf.pathfem = str(fem_dir)
    lf.field = "E"
    if point_electrodes:
        lf.electrode = None
    lf.eeg_cap = str(fname_montage)
    # Tissues to include results from. Default is 1006, which is the eyes,
    # however, we do not want field estimates here, we only want to interpolate
    # the results to central gray matter surface
    lf.tissues = []
    lf.interpolation = "middle gm"
    # Which tissues to interpolate from (default = 2, i.e., gray matter)
    # lf.interpolation_tissue = [2]
    lf.interpolation_subsampling = subsampling

    if init_kwargs:
        for k, v in init_kwargs.items():
            setattr(lf, k, v)

    # To use custom conductivities, do something like
    # s = cond.standard_cond()
    # s[1].value = 0.3   # gray matter
    # s[3].value = 0.006 # bone (average)
    # s[4].value = 0.3   # skin
    # lf.cond = s
    # check s[i].name/descrip to see which indices correspond to which tissues

    lf.run(**run_kwargs)

    return Path(fem_dir) / lf._lf_name()


def prepare_forward(fwd_name: Union[Path, str], apply_average_proj: bool = True):
    """Read a file (.hdf5) containing a leadfield of class TDCSLEADFIELD and
    prepare it for use with EEG. If `apply_average_proj`, the reference
    electrode is reinserted and an average reference is applied. The source
    space is unravelled to left and right hemispheres.

    Currently only supports leadfields which were interpolated to the central
    gray matter surface (exclusively, i.e., the leadfield was not read out in
    any other tissues).

    PARAMETERS
    ----------
    fwd_name : str | Path
        Path to the hdf5 file containing the leadfield.
    apply_average_proj : bool
        Apply a column-wise average reference to the leadfield matrix (default
        = True). If the data contains the actual reference electrode this can
        be turned off to retain the full rank of the data. Average projection
        reduces the data rank by 1.

    RETURNS
    -------
    forward : dict
        A dictionary containing the forward solution.

    NOTES
    -----
    See John Mosher's reply here for a discussion on rereferencing and source
    modeling rank.

        https://neuroimage.usc.edu/forums/t/reference-of-the-gain-matrix-eeg-openmeeg-bem-in-brainstorm/1248
    """
    with h5py.File(fwd_name, "r") as f:
        lf = f["mesh_leadfield"]["leadfields"]["tdcs_leadfield"]

        interp_to_gm = lf.attrs["interpolation"] == "middle gm"
        no_rois = len(lf.attrs["tissues"]) == 0

        # TODO: Support for non-'middle gm' exclusive source spaces.

        # Interpolated solutions are per point, solutions sampled in tissues
        # are per cell (tetrahedra).
        supported = interp_to_gm and no_rois
        if not supported:
            msg = (
                "This function currently only supports leadfields which "
                "have been interpolated to the central gray matter surface "
                "(exclusively)."
            )
            raise NotImplementedError(msg)

        interp_subsampling = lf.attrs["interpolation_subsampling"]
        interp_subsampling = (
            None if interp_subsampling == "None" else int(interp_subsampling)
        )

        # Get channel info
        ch_names = lf.attrs["electrode_names"].tolist()
        ch_ref = lf.attrs["reference_electrode"]

        # Get source space info
        # points = f['mesh_leadfield']['nodes']['node_coord'][:]
        # tris = f['mesh_leadfield']['elm']['node_number_list'][:, :3] - 1
        # point_start_idx = f['mesh_leadfield'].attrs['node_ix']
        # tris_start_idx = f['mesh_leadfield'].attrs['elm_ix']

        # Surface IDs (these are lh and rh)
        # interp_ids = f['mesh_leadfield'].attrs['interp_id'].tolist()

        # Get the leadfield and invert it when going from TES to 'EEG mode' due
        # to reciprocity
        lf = -lf[:]

    # Forward solution
    # Insert the reference channel and rereference to an average reference
    if apply_average_proj:
        lf = np.insert(lf, ch_names.index(ch_ref), np.zeros((1, *lf.shape[1:])), axis=0)
        lf -= lf.mean(0)
    else:
        ch_names.remove(ch_ref)
    nchan, nsrc, nori = lf.shape  # leadfield is always calculated in x, y, z
    assert len(ch_names) == nchan
    assert nori == 3

    return dict(
        data=lf,
        ch_names=ch_names,
        n_channels=nchan,
        n_sources=nsrc,
        n_orientations=nori,
        subsampling=interp_subsampling
    )


def make_forward(
    m2m_dir: Union[Path, str],
    fname_leadfield: Union[Path, str],
    out_format: str,
    info: Union[None, Path, str] = None,
    trans: Union[None, Path, str] = None,
    morph_to_fsaverage: Union[None, int] = 10,
    apply_average_proj=True,
    write: bool = False,
):
    """Create a source space object, a source morph object (to the fsaverage
    template), and a forward object (from a leadfield simulation) in the
    desired format. The outputs can be used for source analysis with the
    appropriate software.

    PARAMETERS
    ----------
    m2m_dir : Path | str
        The directory containing the segmentation and head model results. Used
        to create the source space object.
    fname_leadfield : Path | str
        Filename of the leadfield simulation to convert.
    output_format : str
        The output format of the source, morph, and forward objects. Supported
        formats are `mne` (MNE-Python) and `fieldtrip` (FieldTrip).
    info : Path | str | mne.Info
        Filename or instance of mne.Info (only used [and mandatory] when
        output_format = 'mne') (default = None).
    trans : Path | str | mne.Transform
        Filename or instance of mne.Transform (only used [and mandatory] when
        output_format = 'mne') (default = None).
    morph_to_fsaverage : None | int {10, 40, 160}
        Create a source morph object to fsaverage of the specified resolution.
        The number denotes the (approximate) number of vertices per hemisphere
        such that
            10  ->  10,242 (fsaverage 5)
            40  ->  40,962 (fsaverage 6)
            160 -> 163,842 (fsaverage 7, full resolution)
        The default is 10. Use None to omit morph.
    apply_average_proj : bool

    write : bool
        Whether or not to write the results to disk. The forward object is
        saved to the same directory as `leadfield`. The source space and source
        morph objects are saved to the surface directory in `m2m_dir`
        (default = False).

    RETURNS
    -------
    src : mne.SourceSpaces | dict
        Source space object.
    forward : mne.Forward | dict
        Forward object.
    morph : mne.SourceMorph | dict | None
        Source morph object (None if `morph_to_fsaverage` is None)

    NOTES
    -----
    This function is intended to bridge the gap between SimNIBS and software
    for neurophysiological data analysis. Specifically, the user can run a
    simulation in SimNIBS, export the results to any of the supported formats,
    and use the desired software package for source analysis. We also provide
    a mapping from subject-space to the fsaverage template for group analysis.
    """
    assert out_format in EEG_EXPORT_FORMATS

    forward = prepare_forward(fname_leadfield, apply_average_proj)
    subsampling = forward["subsampling"]

    if out_format == "mne":
        assert info is not None, "`info` must be supplied when out_format='mne'"
        assert trans is not None, "`trans` must be supplied when out_format='mne'"
        # defer import of eeg_mne_tools as it relies on mne which may not be
        # installed
        from simnibs.eeg import utils_mne

        src, morph = utils_mne.setup_source_space(
            m2m_dir, subsampling, morph_to_fsaverage
        )
        forward = utils_mne.make_forward(forward, src, info, trans)

    elif out_format == "fieldtrip":
        from simnibs.eeg import utils_fieldtrip

        src, morph = utils_fieldtrip.setup_source_space(
            m2m_dir, subsampling, morph_to_fsaverage
        )
        forward = utils_fieldtrip.make_forward(forward, src)
    else:
        raise ValueError(f"invalid `out_format` {out_format}")

    if write:
        _write_src_fwd_morph(
            src, forward, morph, fname_leadfield, subsampling, out_format
        )

    return src, forward, morph


def _write_src_fwd_morph(src, forward, morph, fname_leadfield, subsampling, out_format):
    leadfield = Path(fname_leadfield)
    stem = leadfield.stem
    stem += f"_subsampling-{subsampling:d}" if subsampling else ""
    stem += "-{:s}"
    fname_fwd = leadfield.parent / stem.format("fwd")
    fname_morph = leadfield.parent / stem.format("morph")
    fname_src = leadfield.parent / stem.format("src")

    if out_format == "mne":
        import mne

        kw = dict(overwrite=True)

        mne.write_forward_solution(fname_fwd.with_suffix(".fif"), forward, **kw)
        if morph:
            morph.save(fname_morph.with_suffix(".h5"), **kw)
        src.save(fname_src.with_suffix(".fif"), **kw)

    elif out_format == "fieldtrip":
        import scipy.io

        scipy.io.savemat(fname_fwd.with_suffix(".mat"), dict(fwd=forward))
        if morph:
            scipy.io.savemat(fname_morph.with_suffix(".mat"), dict(morph=morph))
        scipy.io.savemat(fname_src.with_suffix(".mat"), dict(src=src))

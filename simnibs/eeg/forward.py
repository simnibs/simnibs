from pathlib import Path
from typing import Union

import h5py
import numpy as np

import simnibs
import simnibs.eeg.utils
from simnibs.segmentation.brain_surface import subsample_surfaces
from simnibs.utils.file_finder import SubjectFiles

HEMISPHERES = ("lh", "rh")

EEG_EXPORT_FORMATS = {"mne", "fieldtrip"}

# MORPH MAPS

# def make_morph_map(points: np.ndarray, surf: dict, n: int = 2):
#     """Create a morph map which allows morphing of values from the nodes in
#     `surf` to `points` by linear interpolation. A morph map is a sparse matrix
#     with dimensions (n_points, n_nodes) where each row has exactly three
#     entries that sum to one. It is created by projecting each point in points
#     onto closest triangle in tris_from and determining the coefficients of each
#     point. Hence, it provides a linear interpolation from surface nodes to
#     points.

#     PARAMETERS
#     ----------
#     points : ndarray
#         The target points (i.e., the points we interpolate to).
#     surf : dict
#         The source points (i.e., the nodes we interpolate from) and the source
#         triangulation.
#     n : int
#         Number of nearest neighbors of each point in points to consider on
#         surf.

#     RETURNS
#     -------
#     mmap : scipy.sparse.csr_matrix
#         The morph map.

#     NOTES
#     -----
#     Testing all points against all triangles is expensive and inefficient, thus
#     we compute an approximation by finding, for each point in `points`, the `n`
#     nearest nodes on `surf` and the triangles to which these points belong.
#     We then test only against these triangles.
#     """

#     # Ensure on unit sphere
#     rr_ = simnibs.eeg.utils.normalize_vectors(surf["points"])
#     points_ = simnibs.eeg.utils.normalize_vectors(points)
#     n_points, d_points = points_.shape

#     surf_ = dict(points=rr_, tris=surf["tris"])

#     pttris = simnibs.transformations._get_nearest_triangles_on_surface(points_, surf_, n)

#     # Find the triangle (in surf) to which each point in points projects and
#     # get the associated weights
#     tris, weights, _, _ = simnibs.transformations._project_points_to_surface(points_, surf_, pttris)

#     rows = np.repeat(np.arange(n_points), d_points)
#     cols = surf_["tris"][tris].ravel()
#     return csr_matrix((weights.ravel(), (rows, cols)), shape=(n_points, len(rr_)))


# def make_morph_maps(src_from: dict, src_to: dict, n: int = 2):
#     """Create morph maps for left and right hemispheres.

#     PARAMETERS
#     ----------
#     src_from : dict
#         Dictionary with keys 'lh' and 'rh'
#         ...
#         corresponding to the
#         spherical registration files (e.g., [l/r]h.sphere.reg.gii).
#     src_to : dict
#         Like `src_from`.

#     RETURNS
#     -------
#     mmaps : dict
#         Dictionary of scipy.sparse.csr_matrix corresponding to the morph map
#         for left and right hemispheres.
#     """
#     morph = simnibs.transformations.SurfaceMorph("linear", n)
#     morph_maps = {}
#     for hemi in HEMISPHERES:
#         morph.fit(src_to[hemi], src_from[hemi])
#         morph_maps[hemi] = morph.morph_mat
#     return morph_maps
#     # return {
#     #     hemi: make_morph_map(src_to[hemi]["points"], src_from[hemi], n)
#     #     for hemi in HEMISPHERES
#     # }


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

    subfiles = SubjectFiles(subpath=str(m2m_dir))

    # Subsampling
    try:
        _ = [subfiles.get_surface(h, subsampling=subsampling) for h in HEMISPHERES]
    except FileNotFoundError:
        if subsampling:
            _ = subsample_surfaces(m2m_dir, n_points=subsampling)

    # The paths should be strings otherwise errors might occur when writing the
    # .hdf5 file
    lf = simnibs.simulation.sim_struct.TDCSLEADFIELD()
    lf.fnamehead = str(subfiles.fnamehead)
    lf.subpath = str(subfiles.subpath)
    lf.pathfem = str(fem_dir)
    lf.field = "E"
    if point_electrodes:
        lf.electrode = None
    # lf.solver_options = "pardiso"
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
    """Read a file (.hdf5) containing a leadfield of class TDCS_LEADFIELD and
    prepare it for use with EEG. The leadfield has the reference electrode
    reinserted and is rereferenced to an average reference. The source space is
    unravelled to left and right hemispheres.

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
            None if interp_subsampling == "None" else interp_subsampling
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
        # interp_ids = interp_ids,
    )

    # # Source space
    # points = np.vsplit(points, [point_start_idx[1]])
    # tris = np.vsplit(tris, [tris_start_idx[1]])
    # tris[1] -= point_start_idx[1] # fix indexing

    # src = dict()
    # for i, p, t in zip(interp_ids, points, tris):
    #     src[i] = dict(points=p, tris=t, nr=len(p))

    # # Normals
    # if interp_subsampling is not None and m2m_dir is not None:
    #     sf = SubjectFiles(subpath=m2m_dir)
    #     for i in interp_ids:
    #         src[i]['n'] = np.loadtxt(sf.get_surface(i, 'central',
    #                                                  interp_subsampling))

    # return forward


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


# def project_points_to_surface_example():
#     """Visualize a projection of random points on two triangles in 2D."""
#     import matplotlib.pyplot as plt
#     import matplotlib.patches as mp
#     import matplotlib.lines as ml

#     mpts = np.array([[1, 1], [2, 2], [-2, 3], [0, -1]], dtype=np.float64)
#     mtris = np.array([[0, 1, 2], [0, 3, 2]])
#     m = mpts[mtris]

#     n = 10
#     points = (np.random.random((2 * n)).reshape(n, 2) - 0.5) * 8
#     surf = dict(points=mpts, tris=mtris)
#     pttris = np.repeat(np.array([0, 1]), n).reshape(2, n).T

#     tris, weights, projs, dists = project_points_to_surface(points, surf, pttris)

#     fig, ax = plt.subplots(1, 1, figsize=(10, 10))
#     for t in m:
#         ax.add_patch(mp.Polygon(t[:, :2], alpha=0.3, edgecolor="k"))
#     for i in range(n):
#         ax.add_artist(
#             ml.Line2D([points[i, 0], projs[i, 0]], [points[i, 1], projs[i, 1]])
#         )
#         ax.add_patch(mp.Circle(points[i], dists[i], alpha=0.5, ec="k", fc="none"))
#     ax.scatter(*points.T, color="r")
#     ax.scatter(*projs.T, color="g")
#     ax.set_xlim((-4, 4))
#     ax.set_ylim((-4, 4))
#     ax.grid(True, alpha=0.2)

#     fig.show()

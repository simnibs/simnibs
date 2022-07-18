from pathlib import Path
from typing import (
    Dict,
    Iterable,
    List,
    Tuple,
    Union,
)  # replace Dict with dict when using python >3.9

import h5py
import nibabel as nib
import numpy as np
import pyvista as pv
import scipy.ndimage as ndi
from scipy.optimize import least_squares
from scipy.sparse import coo_matrix, csr_matrix
from scipy.spatial.ckdtree import cKDTree

import simnibs
from simnibs.utils.file_finder import SubjectFiles

HEMISPHERES = ("lh", "rh")
MORPH_DATA = {"area", "curv", "sulc", "thickness"}
GEOMETRY = {"inflated", "pial", "sphere", "white"}
ANNOT = {
    "aparc",
    "aparc.a2005s",
    "aparc.a2009s",
    "oasis.chubs",
    "PALS_B12_Brodmann",
    "PALS_B12_Lobes",
    "PALS_B12_OrbitoFrontal",
    "PALS_B12_Visuotopic",
    "Yeo2011_7Networks_N1000",
    "Yeo2011_17Networks_N1000",
}
# LABEL = set((""))

# UTILITIES


def load_surface(
    sub_files: SubjectFiles, surf: str, subsampling: Union[int, None] = None
):
    """Load subject-specific surface files into a dictionary. Also load the
    corresponding *normals.txt file if present (and relevant).

    PARARMETERS
    -----------
    sub_files : simnibs.utils.file_finder.SubjectFiles
        SubjectFiles object.
    surf : str
        The surface type to load (e.g., 'central').
    subsampling : int | None
        The subsampling to load (default = None).
    """
    surfaces = {}
    for hemi in sub_files.regions:
        fname = sub_files.get_surface(hemi, surf, subsampling)
        gii = nib.load(fname)
        surfaces[hemi] = dict(
            points=gii.agg_data("pointset"), tris=gii.agg_data("triangle")
        )
        if subsampling is not None and surf == "central":
            fname_nn = "".join((fname.rstrip("gii"), "normals.txt"))
            surfaces[hemi]["normals"] = np.loadtxt(fname_nn)
    return surfaces


def load_subject_surfaces(sub_files: SubjectFiles, subsampling=None):
    """Load subject-specific surface files and fsaverage into dictionaries."""
    # When a leadfield simulation is run, the hemisphere surfaces are appended
    # according to the order in SubjectFiles.regions, hence we load them in the
    # same way here to ensure that the order is the same.
    return (
        load_surface(sub_files, "central", subsampling),
        load_surface(sub_files, "sphere_reg", subsampling),
    )


class FsAverage(object):
    def __init__(self, subdivision: Union[int, str] = 7):
        """A class that provides access to the FsAverage surface files included
        with SimNIBS.

        PARAMETERS
        ----------
        subdivision : int
            The FsAverage subdivision factor (default = 7 which is the full
            resolution fsaverage template). Other available factors are 5 and 6.

        """
        mapper = {5: "fsaverage5", 6: "fsaverage6", 7: "fsaverage"}
        if isinstance(subdivision, str):
            subdivision = subdivision.lstrip("fsaverage")
            if len(subdivision) == 0:
                subdivision = 7
            elif len(subdivision) == 1:
                subdivision = int(subdivision)
            else:
                raise ValueError

        try:
            self.name = mapper[subdivision]
        except KeyError as e:
            raise ValueError(
                f"{subdivision} is not a valid FsAverage subdivision factor."
            ) from e

        self.path = (
            Path(simnibs.SIMNIBSDIR)
            / "resources"
            / "templates"
            / "freesurfer"
            / self.name
        )
        self.subpaths = {
            "surface": self.path / "surf",
            "morph_data": self.path / "surf",
            "annot": self.path / "label",
        }

    def _get_files(self, what, path):
        return {h: self.subpaths[path] / ".".join((h, what)) for h in HEMISPHERES}

    def get_surface(self, surface):
        assert surface in GEOMETRY, f"{surface} is not a valid surface."
        files = self._get_files(surface, "surface")
        return {
            k: dict(zip(("points", "tris"), nib.freesurfer.read_geometry(v)))
            for k, v in files.items()
        }

    def get_central_surface(self):
        central = self.get_surface("white")
        pial = self.get_surface("pial")
        for h, v in central.items():
            v["points"] += pial[h]["points"]
            v["points"] *= 0.5
        return central

    def get_morph_data(self, data):
        assert data in MORPH_DATA, f"{data} is not a valid morph data type."
        files = self._get_files(data, "morph_data")
        return {k: nib.freesurfer.read_morph_data(v) for k, v in files.items()}

    def get_annotation(self, annot):
        # only aparc and aparc.a2009s are available for fsaverage5 and fsaverage6
        assert annot in ANNOT, f"{annot} is not a valid annotation."

        files = self._get_files(f"{annot}.annot", "annot")
        keys = ("labels", "ctab", "names")
        return {
            k: dict(zip(keys, nib.freesurfer.read_annot(v))) for k, v in files.items()
        }


def read_montage(fname: Union[Path, str]):
    ch_types, x, y, z, ch_names = np.genfromtxt(
        fname, dtype=None, delimiter=",", encoding=None, unpack=True,
    )
    if not isinstance(ch_names.dtype, str):
        ch_names = ch_names.astype(str)
    return ch_types, np.column_stack((x, y, z)), ch_names


def load_landmarks(name: Union[Path, str]):
    ch_types, ch_pos, ch_names = read_montage(
        simnibs.utils.file_finder.get_landmarks(name)
    )
    assert all(i == "Fiducial" for i in ch_types)
    return dict(zip(ch_names, ch_pos))


def make_montage(name: Union[Path, str]):
    """Load an EEG montage.

    If `name` is not one of the standard montages included with SimNIBS, it is
    assumed to be the filename of a CSV file defining an EEG montage in the
    same format as those included in simnibs.
    """
    ch_types, ch_pos, ch_names = read_montage(
        simnibs.utils.file_finder.get_montage(name)
    )
    lm = ch_types == "Fiducial"
    if any(lm):
        landmarks = dict(zip(ch_names[lm], ch_pos[lm]))
        not_fids = ~lm
        ch_types = ch_types[not_fids]
        ch_pos = ch_pos[not_fids]
        ch_names = ch_names[not_fids]
    else:
        landmarks = None
    assert all(np.isin(ch_types, ("Electrode", "ReferenceElectrode")))
    return Montage(name, ch_names, ch_pos, ch_types, landmarks)


def apply_trans(trans, arr):
    return np.squeeze(np.atleast_2d(arr) @ trans[:3, :3].T + trans[:3, 3])


class Montage:
    def __init__(
        self,
        name: Union[None, Path, str] = None,
        ch_names: Union[None, List[str], Tuple[str]] = None,
        ch_pos: Union[None, List, Tuple, np.ndarray] = None,
        ch_types: Union[None, List[str], Tuple[str], str] = None,
        landmarks: Union[None, Dict] = None,
        headpoints: Union[None, List, Tuple, np.ndarray] = None,
    ):
        self.name = name
        self.ch_names = [] if ch_names is None else ch_names
        self.ch_pos = np.array([]).reshape(0, 3) if ch_pos is None else np.array(ch_pos)
        self.n_channels = self.ch_pos.shape[0]
        assert self.n_channels == len(self.ch_names)

        if isinstance(ch_types, str):
            ch_types = [ch_types] * self.n_channels
        self.ch_types = [] if ch_types is None else ch_types

        self.landmarks = landmarks or {}
        if self.landmarks:
            self.landmarks = {k: np.asarray(v) for k, v in self.landmarks.items()}
        self.headpoints = (
            np.array([]).reshape(0, 3) if headpoints is None else np.array(headpoints)
        )
        self.n_headpoints = self.headpoints.shape[0]

    def add_landmarks(self, name="Fiducials"):
        """The name of the fiducial file to read fiducials from..."""
        self.landmarks = load_landmarks(name)

    def get_landmark_names(self):
        return list(self.landmarks.keys())

    def get_landmark_pos(self, names: Union[Iterable, None] = None):
        names = names or self.get_landmark_names()
        return np.array([self.landmarks[n] for n in names]).reshape(-1, 3)

    def apply_trans(self, trans):
        """Apply affine transformation."""
        self.ch_pos = apply_trans(trans, self.ch_pos)
        self.headpoints = apply_trans(trans, self.headpoints)
        if self.landmarks:
            for k, v in self.landmarks.items():
                self.landmarks[k] = apply_trans(trans, v)

    def apply_deform(self, deform):
        """Apply a nonlinear deformation field, e.g., MNI2Conform_nonl or
        Conform2MNI_nonl.

        `deform` is a 4D array where the last dimension encodes the coordinates
        corresponding to a particular voxel.
        """
        # `from` coords -> `from` voxels
        self.apply_trans(np.linalg.inv(deform.affine))
        # `from` voxels -> `to` coords
        coords = deform.get_fdata()
        self.ch_pos = self.map_voxels_to_coords(coords, self.ch_pos)
        self.headpoints = self.map_voxels_to_coords(coords, self.headpoints)
        if self.landmarks:
            for k, v in self.landmarks.items():
                self.landmarks[k] = self.map_voxels_to_coords(coords, v)

    @staticmethod
    def map_voxels_to_coords(coords, voxels):
        """
        E.g., if voxels are in 3d then `coords` have shape (i, j, k, 3).
        """
        voxels = np.atleast_2d(voxels)
        ndim = voxels.shape[1]
        assert coords.shape[-1] == ndim
        return np.stack(
            [ndi.map_coordinates(coords[..., i], voxels.T) for i in range(ndim)]
        ).T.squeeze()

    def project_to_surface(self, surf, subset=None):
        n = 2

        split_idx = np.cumsum((self.n_channels, self.n_headpoints))
        points = np.concatenate((self.ch_pos, self.headpoints, self.get_landmark_pos()))

        pttris = get_nearest_triangles_on_surface(points, surf, n, subset)
        _, _, projs, _ = project_points_to_surface(points, surf, pttris)

        self.ch_pos, self.headpoints, proj_landmarks = np.split(projs, split_idx)
        if self.landmarks:
            for k, v in zip(self.get_landmark_names(), proj_landmarks):
                self.landmarks[k] = v

    def fit_to(self, other, **kwargs):
        """Fit the landmarks in `self` to those in `other` in a least squares
        sense, i.e., the transformation matrix moves points from `self` to
        `other`.

        kwargs to pass to `fit_matched_points_analytical`

        """
        assert isinstance(other, Montage)
        assert self.landmarks and other.landmarks
        names = set(self.get_landmark_names()).intersection(
            set(other.get_landmark_names())
        )
        return fit_matched_points_analytical(
            self.get_landmark_pos(names), other.get_landmark_pos(names), **kwargs
        )

    def write(self, filename):
        kw = {"fmt": "%s", "delimiter": ","}
        with open(filename, "w") as f:
            arr = np.column_stack((self.ch_types, self.ch_pos, self.ch_names))
            np.savetxt(f, arr, **kw)
            if self.landmarks:
                n = len(self.landmarks)
                ch_types = ["Fiducial"] * n
                ch_pos = list(self.landmarks.values())
                ch_names = list(self.landmarks.keys())
                arr = np.column_stack((ch_types, ch_pos, ch_names))
                np.savetxt(f, arr, **kw)


# fsavg = FsAverage(5)
# fsavg.get_central_surface()
# fsavg.get_surface("pial")
# fsavg.get_morph_data("thickness")
# fsavg.get_annotation("Yeo2011_17Networks_N1000")

# montage = make_standard_montage("easycap-m10")
# montage.add_landmarks()


def rotation(ax: float, ay: float, az: float) -> np.ndarray:
    """Rotation"""
    rx = np.array(
        [
            [1, 0, 0, 0],
            [0, np.cos(ax), -np.sin(ax), 0],
            [0, np.sin(ax), np.cos(ax), 0],
            [0, 0, 0, 1],
        ]
    )
    ry = np.array(
        [
            [np.cos(ay), 0, np.sin(ay), 0],
            [0, 1, 0, 0],
            [-np.sin(ay), 0, np.cos(ay), 0],
            [0, 0, 0, 1],
        ]
    )
    rz = np.array(
        [
            [np.cos(az), -np.sin(az), 0, 0],
            [np.sin(az), np.cos(az), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
    return rz @ ry @ rx


def scaling(sx: float, sy: float, sz: float):
    return np.array([[sx, 0, 0, 0], [0, sy, 0, 0], [0, 0, sz, 0], [0, 0, 0, 1]])


def translation(tx: float, ty: float, tz: float):
    return np.array([[1, 0, 0, tx], [0, 1, 0, ty], [0, 0, 1, tz], [0, 0, 0, 1]])


def affine_transformation(
    rx: float = 0.0,
    ry: float = 0.0,
    rz: float = 0.0,
    tx: float = 0.0,
    ty: float = 0.0,
    tz: float = 0.0,
    sx: float = 1.0,
    sy: float = 1.0,
    sz: float = 1.0,
):
    """Build an affine transformation matrix from a set of parameters. Applies
    scaling, rotation, and translation in that order.
    """
    s = scaling(sx, sy, sz)
    r = rotation(rx, ry, rz)
    t = translation(tx, ty, tz)
    return t @ r @ s


def fit_matched_points_analytical(
    src_pts, tgt_pts, scale: bool = False, out: str = "trans"
):
    """
    Rigid body transformation with optional scaling.

    PARAMETERS
    ----------
    scale : bool
        If True, estimate a uniform scaling parameter for all axes (default = False).


    REFERENCES
    ----------
    * P.J. Besl and N.D. McKay, A Method for  Registration of 3-D Shapes, IEEE
        Trans. Patt. Anal. Machine Intell., 14, 239 - 255, 1992.
    * Horn, Closed-form solution of absolute orientation using unit
        quaternions, J Opt. Soc. Amer. A vol 4 no 4 pp 629-642, Apr. 1987.
    * Implementation in mne.transforms._fit_matched_points

    """
    assert out in {"params", "trans"}

    # To be consistent with the notation in the references
    p = src_pts
    x = tgt_pts

    n_p = p.shape[0]

    mu_p = p.mean(0)
    mu_x = x.mean(0)
    p_ = p - mu_p
    x_ = x - mu_x

    # Get the optimal rotation
    sigma_px = p_.T @ x_ / n_p  # cross covariance
    Aij = sigma_px - sigma_px.T
    delta = Aij[(1, 2, 0), (2, 0, 1)]

    trace_sigma_px = np.trace(sigma_px)
    Q = np.empty((4, 4))
    Q[0, 0] = trace_sigma_px
    Q[0, 1:] = delta
    Q[1:, 0] = delta
    Q[1:, 1:] = sigma_px + sigma_px.T
    idx = np.arange(1, 4)
    Q[idx, idx] -= trace_sigma_px

    _, v = np.linalg.eigh(Q)

    # q0 of the unit quaternion should be >= 0
    qr = v[1:, -1]  # eigh sorts eigenvalues in ascending order
    if (q0_sign := np.sign(v[0, -1])) != 0:
        qr *= q0_sign

    # Get the scaling
    qs = x_.std() / p_.std() if scale else 1

    # Get the transformation
    qt = mu_x - qs * quaternion_to_rotation(qr) @ mu_p

    params = np.array([*quaternion_to_euler_angles(qr), *qt, qs, qs, qs])
    if out == "params":
        return params
    elif out == "trans":
        return affine_transformation(*params)


def quaternion_real_part(q):
    """Get q0 of q = [q0, q1, q2, q3] where `q` is a unit quaternion."""
    return np.sqrt(np.maximum(0, 1 - np.sum(q ** 2)))


def quaternion_to_rotation(q):
    """Rotation matrix corresponding to the rotation around the axis defined by
    `q`.

    PARAMETERS
    ----------
    q :
        Unit quaternion.

    RETURNS
    -------
    The rotation matrix associated with q.
    """
    q1, q2, q3 = q
    q0 = quaternion_real_part(q)
    qq0 = q0 ** 2
    qq1, qq2, qq3 = q ** 2

    q0q1_2 = 2 * q0 * q1
    q0q2_2 = 2 * q0 * q2
    q0q3_2 = 2 * q0 * q3
    q1q2_2 = 2 * q1 * q2
    q1q3_2 = 2 * q1 * q3
    q2q3_2 = 2 * q2 * q3

    # homogeneous expression
    return np.array(
        [
            [qq0 + qq1 - qq2 - qq3, q1q2_2 - q0q3_2, q1q3_2 + q0q2_2],
            [q1q2_2 + q0q3_2, qq0 - qq1 + qq2 - qq3, q2q3_2 - q0q1_2],
            [q1q3_2 - q0q2_2, q2q3_2 + q0q1_2, qq0 - qq1 - qq2 + qq3],
        ]
    )


def quaternion_to_euler_angles(q):
    """Rotation angles around each (x,y,z) axis."""
    q0 = quaternion_real_part(q)
    q1, q2, q3 = q
    return np.array(
        [
            np.arctan2(2 * (q0 * q1 + q2 * q3), 1 - 2 * (q1 * q1 + q2 * q2)),
            np.arcsin(2 * (q0 * q2 - q1 * q3)),
            np.arctan2(2 * (q0 * q3 + q1 * q2), 1 - 2 * (q2 * q2 + q3 * q3)),
        ]
    )


def fit_matched_points_generic(
    src_pts,
    tgt_pts,
    rotate: bool = True,
    translate: bool = True,
    scale: bool = False,
    init_params=None,
    return_trans: bool = True,
):
    """Find an affine transformation which matches the points in `src_pts`
    to those in `tgt_pts` in a least squares sense.


    PARAMETERS
    ----------


    RETURNS
    -------
    The parameters (rotation, translation, scaling) or as an affine transformation
    matrix.
    """
    assert src_pts.shape == tgt_pts.shape
    assert any((rotate, translate, scale))

    if not rotate and translate:
        raise NotImplementedError(
            "Only implemented with rotation/translation or scaling/rotation/translation"
        )

    apply_params = lambda p, points: apply_trans(affine_transformation(*p), points)
    error_func = lambda p: np.ravel(tgt_pts - apply_params(p, src_pts))

    if init_params is None:
        # Identity
        init_params = []
        if rotate:
            init_params += [0, 0, 0]
        if translate:
            init_params += [0, 0, 0]
        if scale:
            init_params += [1, 1, 1]
    else:
        assert len(init_params) == 3 * rotate + 3 * translate + 3 * scale

    p = least_squares(error_func, init_params, method="lm").x

    return affine_transformation(*p) if return_trans else p


def normalize_vectors(v: np.ndarray):
    """Normalize the rows of `v` avoiding RuntimeWarning due to division by
    zero.

    PARAMETERS
    ----------
    v : ndarray
        Array with vectors in rows.

    RETURNS
    -------
    v (normalized)
    """
    assert v.ndim == 2
    size = np.linalg.norm(v, axis=1, keepdims=True)
    return np.divide(v, size, where=size != 0)


# MORPH MAPS


def get_triangle_neighbors(tris: np.ndarray, nr: Union[int, None] = None):
    """For each point get its neighboring triangles (i.e., the triangles to
    which it belongs).

    PARAMETERS
    ----------
    tris : ndarray
        Array describing a triangulation with size (n, 3) where n is the number
        of triangles.
    nr : int
        Number of points. If None, it is inferred from `tris` as tris.max()+1
        (default = None).

    RETURNS
    -------
    pttris : ndarray
        Array of arrays where pttris[i] are the neighboring triangles of the
        ith point.
    """
    n_tris, n_dims = tris.shape
    nr = tris.max() + 1 if nr is None else nr

    rows = tris.ravel()
    cols = np.repeat(np.arange(n_tris), n_dims)
    data = np.ones_like(rows)
    csr = coo_matrix((data, (rows, cols)), shape=(nr, n_tris)).tocsr()
    return np.array(np.split(csr.indices, csr.indptr[1:-1]), dtype=object)


def sliced_argmin(x: np.ndarray, indptr: np.ndarray):
    """Perform argmin on slices of x.

    PARAMETERS
    ----------
    x : 1-d array
        The array to perform argmin on.
    indptr : 1-d array-like
        The indices of the slices. The ith slice is indptr[i]:indptr[i+1].

    RETURNS
    -------
    res : 1-d array
        The indices (into x!) corresponding to the minimum values in each
        chunk.
    """
    assert x.ndim == 1
    return np.array([x[i:j].argmin() + i for i, j in zip(indptr[:-1], indptr[1:])])


def project_points_to_surface(
    points: np.ndarray,
    surf: Dict,
    pttris: Union[List, np.ndarray],
    return_all: bool = False,
):
    """Project each point in `points` to the closest point on the surface
    described by `surf` restricted to the triangles in `pttris`.

    PARAMETERS
    ----------
    points : ndarray
        Array with shape (n, d) where n is the number of points are d is the
        dimension.
    surf : dict
        Dictionary with keys 'points' and 'tris' describing the surface mesh on
        which the points are to be projected.
    pttris : ndarray | list
        If a ragged/nested array, the ith entry contains the triangles against
        which the ith point will be tested.
    return_all : bool
        Whether to return all projection results (i.e., the projection of a
        point on each of the triangles which it was tested against) or only the
        projection on the closest triangle.

    RETURNS
    -------
    tris : ndarray
        The index of the triangle onto which a point was projected.
    weights : ndarray
        The linear interpolation weights resulting in the projection of a point
        onto a particular triangle.
    projs :
        The coordinates of the projection of a point on a triangle.
    dists :
        The distance of a point to its projection on a triangle.

    NOTES
    -----
    The cost function to be minimized is the squared distance between a point
    P and a triangle T

        Q(s,t) = |P - T(s,t)|**2 =
               = a*s**2 + 2*b*s*t + c*t**2 + 2*d*s + 2*e*t + f

    The gradient

        Q'(s,t) = 2(a*s + b*t + d, b*s + c*t + e)

    is set equal to (0,0) to find (s,t).

    REFERENCES
    ----------
    https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf

    """

    npttris = list(map(len, pttris))
    pttris = np.concatenate(pttris)

    m = surf["points"][surf["tris"]]
    v0 = m[:, 0]  # Origin of the triangle
    e0 = m[:, 1] - v0  # s coordinate axis
    e1 = m[:, 2] - v0  # t coordinate axis

    # Vector from point to triangle origin (if reverse, the negative
    # determinant must be used)
    rep_points = np.repeat(points, npttris, axis=0)
    w = v0[pttris] - rep_points

    a = np.sum(e0 ** 2, 1)[pttris]
    b = np.sum(e0 * e1, 1)[pttris]
    c = np.sum(e1 ** 2, 1)[pttris]
    d = np.sum(e0[pttris] * w, 1)
    e = np.sum(e1[pttris] * w, 1)
    # f = np.sum(w**2, 1)

    # s,t are so far unnormalized!
    s = b * e - c * d
    t = b * d - a * e
    det = a * c - b ** 2

    # Project points (s,t) to the closest points on the triangle (s',t')
    sp, tp = np.zeros_like(s), np.zeros_like(t)

    # We do not need to check a point against all edges/interior of a triangle.
    #
    #          t
    #     \ R2|
    #      \  |
    #       \ |
    #        \|
    #         \
    #         |\
    #         | \
    #     R3  |  \  R1
    #         |R0 \
    #    _____|____\______ s
    #         |     \
    #     R4  | R5   \  R6
    #
    # The code below is equivalent to the following if/else structure
    #
    # if s + t <= 1:
    #     if s < 0:
    #         if t < 0:
    #             region 4
    #         else:
    #             region 3
    #     elif t < 0:
    #         region 5
    #     else:
    #         region 0
    # else:
    #     if s < 0:
    #         region 2
    #     elif t < 0
    #         region 6
    #     else:
    #         region 1

    # Conditions
    st_l1 = s + t <= det
    s_l0 = s < 0
    t_l0 = t < 0

    # Region 0 (inside triangle)
    i = np.flatnonzero(st_l1 & ~s_l0 & ~t_l0)
    deti = det[i]
    sp[i] = s[i] / deti
    tp[i] = t[i] / deti

    # Region 1
    # The idea is to substitute the constraints on s and t into F(s,t) and
    # solve, e.g., here we are in region 1 and have Q(s,t) = Q(s,1-s) = F(s)
    # since in this case, for a point to be on the triangle, s+t must be 1
    # meaning that t = 1-s.
    i = np.flatnonzero(~st_l1 & ~s_l0 & ~t_l0)
    aa, bb, cc, dd, ee = a[i], b[i], c[i], d[i], e[i]
    numer = cc + ee - (bb + dd)
    denom = aa - 2 * bb + cc
    sp[i] = np.clip(numer / denom, 0, 1)
    tp[i] = 1 - sp[i]

    # Region 2
    i = np.flatnonzero(~st_l1 & s_l0)  # ~t_l0
    aa, bb, cc, dd, ee = a[i], b[i], c[i], d[i], e[i]
    tmp0 = bb + dd
    tmp1 = cc + ee
    j = tmp1 > tmp0
    j_ = ~j
    k, k_ = i[j], i[j_]
    numer = tmp1[j] - tmp0[j]
    denom = aa[j] - 2 * bb[j] + cc[j]
    sp[k] = np.clip(numer / denom, 0, 1)
    tp[k] = 1 - sp[k]
    sp[k_] = 0
    tp[k_] = np.clip(-ee[j_] / cc[j_], 0, 1)

    # Region 3
    i = np.flatnonzero(st_l1 & s_l0 & ~t_l0)
    cc, ee = c[i], e[i]
    sp[i] = 0
    tp[i] = np.clip(-ee / cc, 0, 1)

    # Region 4
    i = np.flatnonzero(st_l1 & s_l0 & t_l0)
    aa, cc, dd, ee = a[i], c[i], d[i], e[i]
    j = dd < 0
    j_ = ~j
    k, k_ = i[j], i[j_]
    sp[k] = np.clip(-dd[j] / aa[j], 0, 1)
    tp[k] = 0
    sp[k_] = 0
    tp[k_] = np.clip(-ee[j_] / cc[j_], 0, 1)

    # Region 5
    i = np.flatnonzero(st_l1 & ~s_l0 & t_l0)
    aa, dd = a[i], d[i]
    tp[i] = 0
    sp[i] = np.clip(-dd / aa, 0, 1)

    # Region 6
    i = np.flatnonzero(~st_l1 & t_l0)  # ~s_l0
    aa, bb, cc, dd, ee = a[i], b[i], c[i], d[i], e[i]
    tmp0 = bb + ee
    tmp1 = aa + dd
    j = tmp1 > tmp0
    j_ = ~j
    k, k_ = i[j], i[j_]
    numer = tmp1[j] - tmp0[j]
    denom = aa[j] - 2 * bb[j] + cc[j]
    tp[k] = np.clip(numer / denom, 0, 1)
    sp[k] = 1 - tp[k]
    tp[k_] = 0
    sp[k_] = np.clip(-dd[j_] / aa[j_], 0, 1)

    # Distance from original point to its projection on the triangle
    projs = v0[pttris] + sp[:, None] * e0[pttris] + tp[:, None] * e1[pttris]
    dists = np.linalg.norm(rep_points - projs, axis=1)
    weights = np.column_stack((1 - sp - tp, sp, tp))

    if return_all:
        tris = pttris
    else:
        # Find the closest projection
        indptr = [0] + np.cumsum(npttris).tolist()
        i = sliced_argmin(dists, indptr)
        tris = pttris[i]
        weights = weights[i]
        projs = projs[i]
        dists = dists[i]

    return tris, weights, projs, dists


def make_morph_map(points: np.ndarray, surf: Dict, n: int = 1):
    """Create a morph map which allows morphing of values from the nodes in
    `surf` to `points` by linear interpolation. A morph map is a sparse matrix
    with dimensions (np, n_nodes) where each row has exactly three
    entries that sum to one. It is created by projecting each point in points
    onto closest triangle in tris_from and determining the coefficients of each
    point. Hence, it provides a linear interpolation from surface nodes to
    points.

    PARAMETERS
    ----------
    points : ndarray
        The target points (i.e., the points we interpolate to).
    surf : dict
        The source points (i.e., the nodes we interpolate from) and the source triangulation.
    n : int
        Number of nearest neighbors of each point in points to consider on surf
        (default = 1).

    RETURNS
    -------
    mmap : scipy.sparse.csr_matrix
        The morph map.

    NOTES
    -----
    Testing all points against all triangles is expensive and inefficient, thus
    we compute an approximation by finding, for each point in `points`, the `n`
    nearest nodes on `surf` and the triangles to which these points belong.
    We then test only against these triangles.
    """

    # Ensure on unit sphere
    rr_ = normalize_vectors(surf["points"])
    points_ = normalize_vectors(points)
    n_points, d_points = points_.shape

    surf_ = dict(points=rr_, tris=surf["tris"])

    pttris = get_nearest_triangles_on_surface(points_, surf_, n)

    # Find the triangle (in surf) to which each point in points projects and
    # get the associated weights
    tris, weights, _, _ = project_points_to_surface(points_, surf_, pttris)

    rows = np.repeat(np.arange(n_points), d_points)
    cols = surf_["tris"][tris].ravel()
    return csr_matrix((weights.ravel(), (rows, cols)), shape=(n_points, len(rr_)))


def get_nearest_triangles_on_surface(
    points: np.ndarray, surf: dict, n: int = 1, subset=None
):
    """For each point in `points` get the `n` nearest nodes on `surf` and
    return the triangles to which these nodes belong.

    subset :
        use only a subset of the vertices in surf. should be *indices* not a
        boolean mask.

    RETURNS
    -------
    pttris : list
        Point to triangle mapping.
    """
    assert isinstance(n, int) and n >= 1

    surf_points = surf["points"] if subset is None else surf["points"][subset]
    tree = cKDTree(surf_points)
    _, ix = tree.query(points, n)

    pttris = get_triangle_neighbors(surf["tris"], len(surf["points"]))
    pttris = pttris[ix] if subset is None else pttris[subset[ix]]
    if n > 1:
        pttris = list(map(lambda x: np.unique(np.concatenate(x)), pttris))

    return pttris


def make_morph_maps(src_from: Dict, src_to: Dict, n: int = 1):
    """Create morph maps for left and right hemispheres.

    PARAMETERS
    ----------
    src_from : dict
        Dictionary with keys 'lh' and 'rh' each being a dictionary with keys
        'points' (points) and 'tris' (triangulation) corresponding to the spherical
        registration files (e.g., [l/r]h.sphere.reg.gii).
    src_to : dict
        Like `src_from`.

    RETURNS
    -------
    mmaps : List
        List of scipy.sparse.csr_matrix corresponding to the morph map for left
        and right hemispheres, respectively.
    """
    return [
        make_morph_map(src_to[hemi]["points"], src_from[hemi], n)
        for hemi in HEMISPHERES
    ]


def compute_tdcs_leadfield(
    m2m_dir: Union[Path, str],
    fem_dir: Union[Path, str],
    fname_montage: Union[Path, str],
    subsampling: Union[int, None] = None,
    init_kwargs: Union[Dict, None] = None,
    run_kwargs: Union[Dict, None] = None,
):
    """Convenience function for running a TDCS simulation. The result can be
    converted to an EEG forward solution using `prepare_for_inverse`.

    PARAMETERS
    ----------
    m2m_dir : str
    fem_dir : str
    fname_montage : str
        The filename of the montage to use.
    subsampling : int | None
        The subsampling to use. If the files do not exist, they are created
        (default = None).
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
    if (
        subsampling
        and len((subfiles.get_surface(h, subsampling=subsampling) for h in HEMISPHERES))
        == 0
    ):
        _ = simnibs.segmentation.brain_surface.subsample_surfaces(
            m2m_dir, n_points=subsampling
        )

    # The paths should be strings otherwise errors might occur when writing the
    # .hdf5 file
    lf = simnibs.simulation.sim_struct.TDCSLEADFIELD()
    lf.fnamehead = str(subfiles.fnamehead)
    lf.subpath = str(subfiles.subpath)
    lf.pathfem = str(fem_dir)
    lf.field = "E"
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
    fwd_name : str | path-like
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
    lf = np.insert(lf, ch_names.index(ch_ref), np.zeros((1, *lf.shape[1:])), axis=0)
    ch_names.remove(ch_ref)
    if apply_average_proj:
        lf -= lf.mean(0)
    nchan, nsrc, nori = lf.shape  # leadfield is always calculated in x, y, z
    assert len(ch_names) == nchan

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


def prepare_for_inverse(
    m2m_dir: Union[Path, str],
    fname_leadfield: Union[Path, str],
    out_format: str,
    info: Union[None, Path, str] = None,
    trans: Union[None, Path, str] = None,
    morph_to_fsaverage: int = 5,
    apply_average_proj=True,
    write: bool = False,
):
    """Create a source space object, a source morph object (to the fsaverage
    template), and a forward object (from a leadfield simulation) in the
    desired format. The outputs can be used for source analysis with the
    appropriate software.

    PARAMETERS
    ----------
    m2m_dir : Path-like
        The directory containing the segmentation and head model results. Used
        to create the source space object.
    fname_leadfield : str
        Filename of the leadfield simulation to convert.
    output_format : str
        The output format of the source, morph, and forward objects. Supported
        formats are `mne` (MNE-Python) and `fieldtrip` (FieldTrip).
    info : str | path-like | mne.Info
        Filename or instance of mne.Info (only used [and mandatory] when
        output_format = 'mne') (default = None).
    trans : str | path-like | mne.Transform
        Filename or instance of mne.Transform (only used [and mandatory] when
        output_format = 'mne') (default = None).
    morph_to_fsaverage : int {5,6,7}
        Create a source morph object to fsaverage. Choose from {5, 6, 7} to
        morph to fsaverage5/6/7 (n = 10,242/40,962/163,842 vertices pe
        hemisphere) where 7 is the standard (full resolution) fsaverage
         template (default = 5). Use None to omit morph.
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
    ALLOWED_EEG_FORMATS = ("mne", "fieldtrip")
    assert out_format in ALLOWED_EEG_FORMATS

    forward = prepare_forward(fname_leadfield, apply_average_proj)
    subsampling = forward["subsampling"]

    if out_format == "mne":
        assert info is not None, "`info` must be supplied when out_format='mne'"
        assert trans is not None, "`trans` must be supplied when out_format='mne'"
        # defer import of eeg_mne_tools as it relies on mne which may not be
        # installed
        from simnibs.simulation import eeg_mne_tools

        src, morph = eeg_mne_tools.setup_source_space(
            m2m_dir, subsampling, morph_to_fsaverage
        )
        forward = eeg_mne_tools.make_forward(forward, src, info, trans)

    elif out_format == "fieldtrip":
        from simnibs.simulation import eeg_fieldtrip_tools

        src, morph = eeg_fieldtrip_tools.setup_source_space(
            m2m_dir, subsampling, morph_to_fsaverage
        )
        forward = eeg_fieldtrip_tools.make_forward(forward, src)
    else:
        raise ValueError(f"invalid `out_format` {out_format}")

    if write:
        _write_src_fwd_morph(
            src, forward, morph, fname_leadfield, subsampling, out_format
        )

    return src, forward, morph


def _write_src_fwd_morph(src, forward, morph, fname_leadfield, subsampling, out_format):
    leadfield = Path(fname_leadfield)
    stem = f"{leadfield.stem}_subsampling-{subsampling:d}-{{:s}}"
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

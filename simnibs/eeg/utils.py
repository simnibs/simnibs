from pathlib import Path
from typing import (
    Iterable,
    Union,
)

import nibabel as nib
from nibabel.affines import apply_affine
import numpy as np
import scipy.ndimage as ndi
from scipy.optimize import least_squares

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
VALID_FSAVERAGE_RESOLUTIONS = {10, 40, 160}
ALLOWED_EEG_EXPORT_FORMATS = {"mne", "fieldtrip"}

class FsAverage:
    def __init__(self, resolution: int = 160):
        """A class that provides access to the FsAverage surface files included
        with SimNIBS.

        fsavg40 = FsAverage(40)
        fsavg40.get_central_surface()
        fsavg40.get_surface("pial")
        fsavg40.get_morph_data("thickness")
        fsavg40.get_annotation("Yeo2011_17Networks_N1000")

        PARAMETERS
        ----------
        resolution : int
            The FsAverage resolution factor. Available resolutions are
            160, 40, 10 denoting the approximate number of thousands of
            vertices per hemisphere, e.g., 160 has ~160,000 vertices per
            hemisphere and correspondings to the default `fsaverage` template
            (default = 160).
        """
        if resolution in VALID_FSAVERAGE_RESOLUTIONS:
            self.name = f"fsaverage{resolution:d}k"
        else:
            raise ValueError(f"{resolution} is not a valid resolution for FsAverage. Please choose from {VALID_FSAVERAGE_RESOLUTIONS}")

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


class Montage:


    def __init__(
        self,
        name: Union[None, Path, str] = None,
        ch_names: Union[None, list[str], tuple[str]] = None,
        ch_pos: Union[None, list, tuple, np.ndarray] = None,
        ch_types: Union[None, list[str], tuple[str], str] = None,
        landmarks: Union[None, dict] = None,
        headpoints: Union[None, list, tuple, np.ndarray] = None,
        from_file: Union[None, Path, str] = None,
    ):
        """
        montage = Montage(from_file="easycap-m10")
        # montage = make_standard_montage("easycap-m10")
        montage.add_landmarks()

        """
        if isinstance(from_file, (Path, str)):
            self = make_montage(from_file)
            return
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
        """The name of the file to read fiducials from..."""
        self.landmarks = load_landmarks(name)

    def get_landmark_names(self):
        return list(self.landmarks.keys())

    def get_landmark_pos(self, names: Union[Iterable, None] = None):
        names = names or self.get_landmark_names()
        return np.array([self.landmarks[n] for n in names]).reshape(-1, 3)

    def apply_trans(self, trans):
        """Apply affine transformation."""
        self.ch_pos = apply_affine(trans, self.ch_pos)
        self.headpoints = apply_affine(trans, self.headpoints)
        if self.landmarks:
            for k, v in self.landmarks.items():
                self.landmarks[k] = apply_affine(trans, v).squeeze()

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
        pttris = simnibs.transformations._get_nearest_triangles_on_surface(points, surf, n, subset)
        _, _, projs, _ = simnibs.transformations._project_points_to_surface(points, surf, pttris)

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


# Point fitting

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


def affine_from_params(
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
    out : str
        Whether to return an affine transformation ("trans") or the raw
        parameters ("params") (default = "trans").

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
        return affine_from_params(*params)


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
    src_pts: ndarray

    tgt_pts: ndarray

    rotate: bool = True

    translate: bool = True

    scale: bool = False

    init_params: =None

    return_trans: bool = True


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

    apply_params = lambda p, points: apply_affine(affine_from_params(*p), points)
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

    return affine_from_params(*p) if return_trans else p


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

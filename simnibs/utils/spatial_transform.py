from nibabel.affines import apply_affine
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation


def rotation(
    ax: float, ay: float, az: float, ref_frame="extrinsic", degrees=False
) -> np.ndarray:
    """Construct a rotation matrix from Euler angles.

    PARAMETERS
    ----------
    ax, ay, az : float
        Euler angles.
    ref_frame : str
        Whether to interpret the Euler angles as extrinsic (default) or
        intrinsic rotations.
    degrees : bool
        Interpret angles as degrees.

    RETURNS
    -------
    m : ndarray
        Rotation matrix.

    REFERENCES
    ----------
    https://en.wikipedia.org/wiki/Euler_angles
    https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
    """
    m = np.identity(4)
    seq = "xyz"
    if ref_frame == "intrinsic":
        seq = seq.upper()
    m[:3,:3] = Rotation.from_euler(seq, [ax, ay, az], degrees).as_matrix()
    return m

def scaling(sx: float, sy: float, sz: float):
    return np.array([[sx, 0, 0, 0], [0, sy, 0, 0], [0, 0, sz, 0], [0, 0, 0, 1]])


def translation(tx: float, ty: float, tz: float):
    return np.array([[1, 0, 0, tx], [0, 1, 0, ty], [0, 0, 1, tz], [0, 0, 0, 1]])


def matrix_from_params(
    rx: float = 0.0,
    ry: float = 0.0,
    rz: float = 0.0,
    tx: float = 0.0,
    ty: float = 0.0,
    tz: float = 0.0,
    sx: float = 1.0,
    sy: float = 1.0,
    sz: float = 1.0,
    ref_frame = "extrinsic",
):
    """Build an affine transformation matrix from a set of parameters. Applies
    scaling, rotation, and translation in that order.
    """
    s = scaling(sx, sy, sz)
    r = rotation(rx, ry, rz, ref_frame)
    t = translation(tx, ty, tz)
    return t @ r @ s


def fit_matched_points_analytical(
    src_pts, tgt_pts, scale: bool = False, return_matrix: bool = True,
):
    """Estimate a transformation (rigid body with optional scaling) which moves
    `src_pts` to fit `tgt_pts` in a least squares sense.

    PARAMETERS
    ----------
    src_pts : ndarray
        Source points.
    tgt_pts : ndarray
        Target points.
    scale : bool
        If True, estimate a uniform scaling parameter for all axes (default = False).
    return_matrix : bool
        If True, return an affine transformation matrix. Otherwise, return the
        raw parameters (default = True).

    REFERENCES
    ----------
    * P.J. Besl and N.D. McKay, A Method for  Registration of 3-D Shapes, IEEE
        Trans. Patt. Anal. Machine Intell., 14, 239 - 255, 1992.
    * Horn, Closed-form solution of absolute orientation using unit
        quaternions, J Opt. Soc. Amer. A vol 4 no 4 pp 629-642, Apr. 1987.
    * Implementation in mne.transforms._fit_matched_points

    """
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

    # scipy.spatial.transform.Rotation expects quaternions of the form
    # [v0, v1, v2, s] where s is the scalar (real) part and v0-v2 is the vector
    # (imaginary) part
    v = np.roll(v, -1, 0)
    qr = v[:, -1]
    if (q0_sign := np.sign(qr[-1])) != 0:
        qr *= q0_sign
    rotation = Rotation.from_quat(qr)

    # Get the scaling
    qs = x_.std() / p_.std() if scale else 1

    # Get the translation
    qt = mu_x - qs * rotation.as_matrix() @ mu_p

    # get the extrinsic rotation parameters
    params = np.array([*rotation.as_euler("xyz"), *qt, qs, qs, qs])

    return matrix_from_params(*params) if return_matrix else params


def fit_matched_points_generic(
    src_pts,
    tgt_pts,
    rotate: bool = True,
    translate: bool = True,
    scale: bool = False,
    init_params=None,
    return_matrix: bool = True,
):
    """Estimate a transformation (rigid body with optional scaling) which moves
    `src_pts` to fit `tgt_pts` in a least squares sense.

    PARAMETERS
    ----------
    src_pts: ndarray
        Source points.
    tgt_pts: ndarray
        Target points.
    rotate: bool
        Whether to allow rotation or not (default = True).
    translate: bool
        Whether to allow translation or not (default = True).
    scale: bool
        Whether to allow scaling or not (default = False).
    init_params: None | list
        Initial parameters. If None, will use zero translation/rotation and
        scaling of 1 in all directions. If list, a list of nine elements.
    return_matrix: bool


    RETURNS
    -------
    The parameters (rotation, translation, scaling) or as an affine
    transformation matrix.
    """
    assert src_pts.shape == tgt_pts.shape
    assert any((rotate, translate, scale))

    if not rotate and translate:
        raise NotImplementedError(
            "Only implemented with rotation/translation or scaling/rotation/translation"
        )

    apply_params = lambda p, points: apply_affine(matrix_from_params(*p), points)
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

    return matrix_from_params(*p) if return_matrix else p

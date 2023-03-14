import numpy as np
from scipy.optimize import least_squares

# replace with functionality from scipy.spatial.transform?


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
    return np.sqrt(np.maximum(0, 1 - np.sum(q**2)))


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
    qq0 = q0**2
    qq1, qq2, qq3 = q**2

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

    rotate: bool

    translate: bool

    scale: bool

    init_params: bool

    return_trans: bool


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

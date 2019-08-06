"""
Functions to find optimal TMS coil positions for a given cortical target. See examples/tms_optimization.py for
some examples on how to use these.

Written by Ole Numssen & Konstantin Weise, 2019.
"""
import os
import glob
import numpy as np
import multiprocessing
import datetime
import logging
import h5py
import itertools
import pandas as pd
import scipy.io
import copy

import simnibs
from ..msh import mesh_io
from ..simulation import SESSION, SimuList, TMSLIST, save_matlab_sim_struct, FEMSystem, coil, calc_fields
from ..utils.matlab_read import try_to_read_matlab_field, remove_None
from ..utils.file_finder import SubjectFiles
from ..utils.simnibs_logger import logger

try:
    import pyfempp
except ImportError:
    print("pyfempp not found. Some visualizations will not work.")
    pyfempp = False


def eval_optim(simulations, target, tms_optim, res_fn, qoi="normE", elmtype='elms'):
    """
    Evaluation of optimization results.

    Parameters
    ----------
    simulations: str
        .hdf5 filename of simulation.

    target: list of float [3]
        Cortical target position in subject space [x,y,z].

    tms_optim: simnibs.simulation.sim_struct.TMSOPTIMIZATION() or str
        TMSOPTIMIZATION() object or filename of such.

    res_fn: str
        Filename for results .csv file.

    qoi: str
        Quantity of interest. Which field type to find optimum for. If multidimensional (like 'E'), the mean computed.

    elmtype: string
        Read field from 'nodes' or 'elms'.

    Returns
    --------
    simnibs.simulation.sim_struct.TMSOPTIMIZATION()
        TMSOPTIMIZATION() object with optimal conditions in TMSOPTIMIZATION.optimlist

    """
    # load tms_list if filename
    if type(tms_optim) == str:
        tms_list_new = TMSOPTIMIZATION()
        tms_list_new.read_mat_struct(tms_optim)
        tms_optim = tms_list_new

    logger.info("Evaluating simulations to find optimal coil position")

    # check if all simulations were completed
    if tms_optim.get_unfinished_sim_idx():
        logger.warn("Found unfinished simulations in {}. "
                    "Set TMSOPTIMZATION.resume = True and run optimization again.".format(tms_optim.fname_hdf5))

    # get coil distance
    assert all(pos.distance == tms_optim.optimlist.pos[0].distance for pos in tms_optim.optimlist.pos), \
        "Different coil distances found."
    distance = tms_optim.distance

    if elmtype == 'elms':
        path = '/elmdata/'
    elif elmtype == 'nodes':
        path = '/nodedata/'
    else:
        raise ValueError('elmtype {} is unknown. Use "elms" for elmdata or "nodes" for nodedata'.format(elmtype))

    path += qoi
    with h5py.File(simulations, 'r') as f:
        res = f[path]

        assert res.shape[-1] == len(tms_optim.optimlist.pos), \
            'Simulation number ({}) does not fit to TMSLIST positions ({})'.format(len(simulations),
                                                                                   len(tms_optim.optimlist.pos))
        assert res.shape[-1] > 1, 'Only one simulation found. No optimization possible.'

        if not tms_optim.target or tms_optim.target != target or not tms_optim.target_idx:
            # read first mesh to get index of target element
            try:
                # try to read mesh geom from .hdf5
                mesh = mesh_io.Msh()
                mesh = mesh.read_hdf5(simulations, load_data=False)
                _, idx = mesh.find_closest_element(target, return_index=True)
            except (ValueError, KeyError):
                mesh = mesh_io.read_msh(tms_optim.fnamehead)
                _, idx = mesh.find_closest_element(target, return_index=True)
        else:
            idx = tms_optim.target_idx

        qoi_in_target = np.linalg.norm(res[idx - 1], axis=0)
        assert qoi_in_target.shape == (len(tms_optim.optimlist.pos),)

        # build dictionary to easily create pandas afterwards
        d = dict()
        for i in range(qoi_in_target.shape[-1]):
            centre = tms_optim.optimlist.pos[i].centre
            if not centre:
                centre = [np.nan, np.nan, np.nan]

            ydir = tms_optim.optimlist.pos[i].pos_ydir
            if not ydir:
                ydir = [np.nan, np.nan, np.nan]

            d[i] = [i, qoi_in_target[i]] + centre + ydir + [idx - 1] + target
        data = pd.DataFrame.from_dict(d)
        data = data.transpose()
        data.columns = ['sim_nr'] + [qoi] + ['x', 'y', 'z'] + \
                       ['handle_x', 'handle_y', 'handle_z',
                        'idx_hdf5',
                        'target_x', 'target_y', 'target_z']
        logger.info("Saving {} for all simulations in target element {} to {}".format(qoi, idx, res_fn))
        data.to_csv(res_fn, index=False)  # save normE from all simulations at target elm to .csv

        if pyfempp:
            # write stimsite xdmf with qoi
            stimsite_vis_fn = os.path.splitext(res_fn)[0] + ".hdf5"
            pyfempp.create_stimsite_from_tmslist(stimsite_vis_fn,
                                                 tms_optim.optimlist, qoi,
                                                 qoi_in_target, overwrite=False)

        # get best coil position/orientation from the results
        data[qoi] = data[qoi].astype(float)
        best_cond = data.iloc[[data[qoi].idxmax()]]
        tms_optim = copy.copy(tms_optim)
        tms_optim.optimlist.pos = []  # use original tms_list and replace pos with best conditions
        if not (best_cond.count(0) == 1).all():
            simnibs.utils.simnibs_logger.logger.warn("Multiple optimal coil pos/rot found!")

        for _, row in best_cond.iterrows():
            best_pos = simnibs.simulation.sim_struct.POSITION()

            best_pos.distance = distance
            best_pos.centre = row[['x', 'y', 'z']].values.tolist()
            best_pos.pos_ydir = row[['handle_x', 'handle_y', 'handle_z']].values.tolist()
            best_pos.calc_matsimnibs(msh=tms_optim.fnamehead, log=False)
            logger.info("Optimal coil position for target element {}: \n{}".format(idx, best_pos.matsimnibs))
            tms_optim.optimlist.add_position(best_pos)

    return tms_optim


def get_opt_grid(tms_optim, mesh,
                 target=None, handle_direction_ref=None, distance=1., radius=20,
                 resolution_pos=1, resolution_angle=10, angle_limits=None):
    """
    Determine the coil positions and orientations for bruteforce TMS optimization

    Parameters
    ----------
    tms_optim: simnibs.simulation.sim_struct.TMSOPTIMIZATION object
        TMS simulation object instance
    mesh: simnibs.msh.mesh_io.Msh object
        Simnibs mesh object.
    target: list of float or np.ndarray or None
        Coordinates (x, y, z) of cortical target.
    handle_direction_ref: list of float or np.ndarray
        Vector of handle prolongation direction.
    distance: float or None
        Coil distance to skin surface [mm]. Default: 1.
    radius: float or None
        Radius of region of interest around skin-projected cortical target.
    resolution_pos: float or None
        Resolution in mm of the coil positions in the region of interest.
    resolution_angle: float or None
        Resolution in deg of the coil positions in the region of interest. Default: 20.
    angle_limits: list of float or None
        Range of angles to get coil rotations for (Default: [-60, 60]).

    Returns
    -------
    tms_optim: simnibs.simulation.sim_struct.TMSOPTIMIZATION
        TMS simulation object instance.
    """
    if not angle_limits:
        angle_limits = [-60, 60]

    if not tms_optim.angle_limits:
        tms_optim.angle_limits = angle_limits
    if not tms_optim.handle_direction_ref:
        tms_optim.handle_direction_ref = handle_direction_ref
    if not tms_optim.resolution_angle:
        tms_optim.resolution_angle = resolution_angle
    if not tms_optim.resolution_pos:
        tms_optim.resolution_pos = resolution_pos
    if not tms_optim.target:
        tms_optim.target = target
    if not tms_optim.radius:
        tms_optim.radius = radius

    # project cortical target to skin surface
    tms_tmp = copy.deepcopy(tms_optim.optimlist)
    pos_target = tms_tmp.add_position()
    pos_target.centre = tms_optim.target
    pos_target.pos_ydir = tms_optim.handle_direction_ref
    pos_target.distance = .1
    target_matsimnibs = pos_target.calc_matsimnibs(mesh, log=False)
    tms_optim.target_coil_matsim = target_matsimnibs

    _, tms_optim.target_idx = mesh.find_closest_element(tms_optim.target, return_index=True)
    tms_optim.target = tms_optim.target
    target_skin = target_matsimnibs[0:3, 3]

    # extract ROI
    msh_surf = mesh.crop_mesh(elm_type=2)
    msh_skin = msh_surf.crop_mesh([5, 1005])
    elm_center = np.mean(msh_skin.nodes.node_coord[msh_skin.elm.node_number_list[:, 0:3] - 1], axis=1)
    elm_mask_roi = np.linalg.norm(elm_center - target_skin, axis=1) < 1.2 * tms_optim.radius
    # elm_center_roi = elm_center[elm_mask_roi, ]
    node_number_list_roi = msh_skin.elm.node_number_list[elm_mask_roi, 0:3] - 1
    nodes_roi = msh_skin.nodes.node_coord[node_number_list_roi]
    nodes_roi = np.reshape(nodes_roi, (nodes_roi.shape[0] * 3, 3))
    nodes_roi_mean = np.mean(nodes_roi, axis=0)[np.newaxis, :]
    nodes_roi_zeromean = nodes_roi - nodes_roi_mean

    # tangential plane of target_skin point
    u, s, vh = np.linalg.svd(nodes_roi_zeromean)
    vh = vh.transpose()

    # define regular grid and rotate it to head space
    coords_plane = np.array(np.meshgrid(np.linspace(-tms_optim.radius,
                                                    tms_optim.radius,
                                                    int(2 * tms_optim.radius / tms_optim.resolution_pos + 1)),
                                        np.linspace(-tms_optim.radius,
                                                    tms_optim.radius,
                                                    int(2 * tms_optim.radius / tms_optim.resolution_pos + 1)))
                            ).T.reshape(-1, 2)
    coords_plane = coords_plane[np.linalg.norm(coords_plane, axis=1) <= tms_optim.radius]
    coords_plane = np.dot(coords_plane, vh[:, :2].transpose()) + target_skin

    # project grid-points to skin-surface
    p1 = msh_skin.nodes.node_coord[node_number_list_roi[:, 0]]
    p2 = msh_skin.nodes.node_coord[node_number_list_roi[:, 1]]
    p3 = msh_skin.nodes.node_coord[node_number_list_roi[:, 2]]

    normals = np.cross(p2 - p1, p3 - p1)
    # coords_mapped = np.zeros(coords_plane.shape)
    coords_mapped = []

    for i, c in enumerate(coords_plane):
        q1 = c + 1e2 * vh[:, 2]
        q2 = c - 1e2 * vh[:, 2]

        # point of intersection of infinite half-plane from triangle
        p0 = q1[np.newaxis, :] + (np.sum((p1 - q1) * normals, axis=1) /
                                  np.dot(q2 - q1, normals.transpose()))[:, np.newaxis] * (q2 - q1)[np.newaxis, :]

        # check if intersection points are inside triangle
        inside = np.logical_and(
            np.logical_and(np.sum(np.cross(p2 - p1, p0 - p1) * normals, axis=1) >= 0,
                           np.sum(np.cross(p3 - p2, p0 - p2) * normals, axis=1) >= 0),
            np.sum(np.cross(p1 - p3, p0 - p3) * normals, axis=1) >= 0)
        if np.sum(inside) == 1:
            coords_mapped.append(np.squeeze(p0[inside]))
        else:
            logger.debug(f"Cannot map coord for radius: {tms_optim.radius}, plane:: {i}. "
                         f"Try a smaller radius (or > 0).")
        # coords_mapped[i] = p0[inside]

    # determine rotation matrices around z-axis of coil and rotate
    angles = np.linspace(tms_optim.angle_limits[0],
                         tms_optim.
                         angle_limits[1],
                         int((tms_optim.angle_limits[1] - tms_optim.angle_limits[0]) / tms_optim.resolution_angle + 1))
    handle_directions = np.zeros((len(angles), 3))

    for i, a in enumerate(angles):
        x = target_matsimnibs[0:3, 0]
        y = target_matsimnibs[0:3, 1]
        handle_directions[i] = np.cos(a / 180. * np.pi) * y + np.sin(a / 180. * np.pi) * -x

    # combine coil positions and orientations
    po = list(itertools.product(coords_mapped, handle_directions))
    n_pos, starttime = len(po), 0
    assert n_pos, "Cannot create any positions from given arguments. "

    n_zeros = len(str(n_pos))

    # store coordinates in TMS object
    for i, val in enumerate(po):
        c, h = val
        pos = tms_optim.optimlist.add_position()
        pos.centre = c.tolist()
        pos.pos_ydir = (h + c).tolist()

        if not i % 100:  # on every 100th iteration print status
            if i:
                duration = datetime.datetime.now() - starttime
                try:
                    time_left = ((n_pos - i) / 100.) * duration
                    time_left = datetime.timedelta(seconds=time_left.seconds)
                except OverflowError:
                    time_left = 'unknown'
                simnibs.utils.simnibs_logger.logger.info("Determining coil position "
                                                         "{0:0>{3}}/{1} (time left: {2}).".format(i,
                                                                                                  n_pos,
                                                                                                  time_left,
                                                                                                  n_zeros))
            else:
                simnibs.utils.simnibs_logger.logger.info("Determining coil position "
                                                         "{0:0>{2}}/{1}.".format(i,
                                                                                 n_pos,
                                                                                 n_zeros))
            starttime = datetime.datetime.now()

        pos.matsimnibs = None
        pos.distance = distance
        pos.matsimnibs = pos.calc_matsimnibs(mesh, log=False, msh_surf=msh_surf)

    return tms_optim


def optimize_tms_coil_pos(tms_optim, target=None,
                          handle_direction_ref=None, radius=20, angle_limits=None,
                          resolution_pos=1.5, resolution_angle=15, distance=1., n_cpu=1):
    """Compute optimal TMS coil position/rotation for maximum field at cortical target.

    This implements a brute-force approach to find the optimal coil position & rotation to maximize normE
    at a given cortical target. The target is projected to the nearest position at the skin surface and candidate coil
    positions/rotations are then distributed equally in a circular plane around that location.
    The number of candiate coil positions/rotations is defined by the two resolution parameters (_pos and
    angle), the radius of the circular plane, and angle_limits.

    All simulations are stored in one .hdf5 files.

    Parameters
    ----------
    target: np.ndarray or list of float
        XYZ coordinates to optimize coil position/rotation for.
    tms_optim: simnibs.simulation.stim_struct.TMSOPTIMIZATION
        Session object with 1 TMSLIST() as poslist.
        If session.poslists[0].pos is empty, is is filled according to parameters below.
    handle_direction_ref: list of num or np.ndarray
        Vector of coil handle prolongation direction. Defaults to left M1 45° ([-2.8, 7, 8.1]).
    radius: float
        Radius of region of interest around skin-projected cortical target, where the bruteforce simulations are
        conducted (Default: 20).
    resolution_pos: float
        Resolution in mm of the coil positions in the region of interest (Default: 1.5).
    resolution_angle: float (Default: 15)
        Resolution in deg of the coil positions in the region of interest.
    angle_limits: list of int or np.ndarray
        Range of angles to get coil rotations for. Default: [-60, 60].
    distance: float
        Coil distance to scalp (Default: 1.).
    n_cpu: Int
        Number of cores to use.

    Returns
    -------
    file
        all simulations are stored in tms_optim.fname_hdf5

    file
        |tms_optim.qoi| at tms_optim.target element for all simulations
        is saved to tms_optim.pathfem/tms_optim.optim_name + tms_optim.timestr + .csv

    file
        TMSOPTIMIZATION object used for simulation in tms_optim.pathfem/tms_optim.optim_name + tms_optim.time_str + .mat

    file
        TMSOPTIMIZATION object with optimal coil positions in
        tms_optim.pathfem/tms_optim.optim_name + tms_optim.time_str + optimal_coil_pos.mat

    dict
        'best_conds': simnibs.sim_struct.TMSLIST for optimal coil position(s)
        'simulations': simulation filenames,
        'tms_optim': simnibs.sim_struct.SESSION object used for simulation
    """

    # some argument checks, defaults, type handling
    if angle_limits is None:
        angle_limits = [-60, 60]
    if handle_direction_ref is None:
        handle_direction_ref = np.array([-2.8, 7, 8.1])

    assert tms_optim.optimlist, "Please provide simnibs.simulation.sim_struct.TMSOPTIMIZATION with TMSLIST optimlist."

    if type(target) == np.ndarray:
        target = target.tolist()
    if type(handle_direction_ref) == np.ndarray:
        handle_direction_ref = handle_direction_ref.tolist()
    if type(resolution_angle) == np.ndarray:
        resolution_angle = resolution_angle.tolist()
    if type(angle_limits) == np.ndarray:
        angle_limits = angle_limits.tolist()

    tms_optim.angle_limits = angle_limits
    tms_optim.handle_direction_ref = handle_direction_ref
    tms_optim.resolution_angle = resolution_angle
    tms_optim.resolution_pos = resolution_pos
    tms_optim.target = target
    tms_optim.radius = radius

    assert tms_optim.radius >= tms_optim.resolution_pos
    if tms_optim.optimlist.pos and \
            (tms_optim.target or tms_optim.handle_direction_ref or
             tms_optim.resolution_angle or tms_optim.angle_limits or tms_optim.distance) and \
            not tms_optim.resume:
        simnibs.utils.simnibs_logger.logger.warn("tms_list positions provided. Ignoring other provided arguments!")

    if not tms_optim.optimlist.pos and not tms_optim.target:
        raise ValueError("Provide either target or tms_optim.optimlist.pos .")

    # setup tmslist if not provided
    if not tms_optim.optimlist.pos:

        mesh = mesh_io.read_msh(tms_optim.fnamehead)
        mesh.fix_surface_labels()

        # update session poslist with changed tms_list
        tms_optim = get_opt_grid(tms_optim=tms_optim,
                                 mesh=mesh)

    else:
        simnibs.utils.simnibs_logger.logger.info("Optimization: using positions from provided TMSLIST.pos .")

    if pyfempp:
        pyfempp.create_stimsite_from_tmslist(os.path.join(tms_optim.pathfem,
                                                          tms_optim.optim_name + "_coil_positions.hdf5"),
                                             tms_optim.optimlist, overwrite=True)

    # run simulations
    simnibs.utils.simnibs_logger.logger.info("Starting optimization: "
                                             "{} FEMs on {} cpus".format(len(tms_optim.optimlist.pos), n_cpu))

    if tms_optim.resume:
        # create copy of tms_optim, as resuming changes tms_optim.optimlist.pos
        # tms_optim_org = copy.copy(tms_optim)
        optimlist_org = copy.copy(tms_optim.optimlist.pos)
        simulations_fn = tms_optim.run(allow_multiple_runs=True, cpus=n_cpu)
        tms_optim.optimlist.opt = optimlist_org
    else:
        simulations_fn = tms_optim.run(allow_multiple_runs=True, cpus=n_cpu)

    # evaluate all simulations
    tms_list_optim = eval_optim(simulations_fn,
                                target,
                                tms_optim,
                                os.path.join(tms_optim.pathfem,
                                             tms_optim.optim_name + '_' + tms_optim.time_str + '.csv'),
                                qoi=tms_optim.qoi)

    save_matlab_sim_struct(tms_list_optim,
                           os.path.join(tms_list_optim.pathfem,
                                        f"{tms_list_optim.optim_name}_{tms_list_optim.time_str}_optimal_coil_pos.mat"
                                        ))
    logging.shutdown()

    # return best position
    return {'best_conds': tms_list_optim,
            'simulations': simulations_fn,
            'tms_optim': tms_optim}


class TMSOPTIMIZATION(SESSION):
    """
    Defines the structure to optimize the TMS coil position/rotation for a given cortical target.

    Use one TMSOPTIMIZATION per optimization target. That means:
    - only one TMSLIST is accepted, with
        - same mesh & anisotropy setup for all positions
        - same coil for all positions
        - same intensity for all positions

    Attributes
    ----------
    angle_limits: list of float or None
        Range of angles to get coil rotations for. Default: [-60, 60].
    anisotropy_type: str
        Type of anisotropy for simulation. ('scalar', 'vn', or 'mc').
    compress_hdf5: str or None
        Compression used for .hdf5. 'gzip' or None.
    cond: list
        List of COND structures with conductivity information.
    didt: float
        Change of current in coil, in A/s.
    distance: float or None
        Coil distance to skin surface [mm]. Default: 1.
    fname_hdf5: str
        Filename for the qoi results.
    fnamecoil: str
        Filename of coil file.
    handle_direction_ref: list of float or np.ndarray[float]
        Vector of handle prolongation direction.
    resolution_angle: float or None
        Resolution in deg of the coil positions in the region of interest. Default: 20.
    resolution_pos: float or None
        Resolution in mm of the coil positions in the region of interest.
    mesh: simnibs.mesh_io.Msh
        Mesh where the simulations will be performed.
    n_sim: int
        Number of simulations.
    optim_name: str
        Filename prefix for result files.
    optimlist: sim_struct.TMSLIST()
        List of coil positions/rotations.
    qoi: str
        Quality of interest to maximize, a member of self.save_fields. 'E', 'normE', ...
    radius: float or None
        Radius of region of interest around skin-projected cortical target.
    resume: bool
        Resume an unfinished optimization. Default: False.
    save_fields: str
        Which fields to store in <optim_name>.hdf5. 'normE' | 'E' | list of both.
    target: ndarray
        Target location (X,Y,Z) at which the maximum qoi is determined. target.shape == (3,)
    target_coil_matsim: np.ndarray[float]
        Matsimnibs for start coil position/orientation according to self.target.
    target_idx: int
        Nearest headmesh element index according to target.
    write_mesh_geom: bool
        Add mesh geometry to <optim_name>.hdf5. Default: True.
    """

    def __init__(self, matlab_struct=None):
        SESSION.__init__(self, matlab_struct)

        # SESSION fields
        del self.poslists
        self.open_in_gmsh = False

        # new fields
        self.anisotropy_type = 'scalar'
        self.compress_hdf5 = None  # None | 'gzip'
        self.cond = None
        self.didt = None
        self.distance = None
        self.fname_hdf5 = None
        self.fnamecoil = None
        self.mesh = None
        self.n_sim = None
        self.optim_name = ''
        self.optimlist = None
        self.qoi = 'E'
        self.resume = False
        self.save_fields = None
        self.write_mesh_geom = True

        # target information
        self.angle_limits = None
        self.handle_direction_ref = None
        self.radius = None
        self.resolution_angle = None
        self.resolution_pos = None
        self.target = None
        self.target_coil_matsim = None
        self.target_idx = None

    def _set_logger(self, fname_prefix=None, summary=False):
        """Set logfile name to self.optim_name"""
        if not fname_prefix:
            fname_prefix = self.optim_name
        super()._set_logger(fname_prefix=fname_prefix)

    def add_poslist(self, pl):
        """ Adds a SimList object to the poslist variable

        Parameters:
        ----------------
        pl: sim_struct.SimuList
            SimuList object
        """
        if not isinstance(pl, SimuList):
            raise TypeError('Elements in poslist must be subclasses of SimuList')

        # we only want one poslist for the optimization
        self.optimlist = pl
        self._prepared = False

    def read_mat_struct(self, mat):
        if type(mat) == str:
            mat = scipy.io.loadmat(mat)

        super().read_mat_struct(mat)

        self.anisotropy_type = try_to_read_matlab_field(mat, 'anisotropy_type', str, 'scalar')
        self.cond = try_to_read_matlab_field(mat, 'cond', str, self.date)
        self.compress_hdf5 = try_to_read_matlab_field(mat, 'compress_hdf5', str)
        self.didt = try_to_read_matlab_field(mat, 'didt', float)
        self.distance = try_to_read_matlab_field(mat, 'distance', float)
        self.fnamecoil = try_to_read_matlab_field(mat, 'fnamecoil', str)
        self.fname_hdf5 = try_to_read_matlab_field(mat, 'hdf5_fn', str)
        if len(mat['optimlist']) > 0:
            for PL in mat['optimlist'][0]:
                if PL['type'][0] == 'TMSLIST':
                    self.add_poslist(TMSLIST(PL))
                else:
                    raise IOError(
                        "poslist type is not of type TMSLIST")
        self.n_sim = try_to_read_matlab_field(mat, 'n_sim', int)
        self.qoi = try_to_read_matlab_field(mat, 'qoi', str)

        # target information
        self.target = try_to_read_matlab_field(mat, 'target', np.array)
        self.target_idx = try_to_read_matlab_field(mat, 'target_idx', int)
        self.target_coil_matsim = try_to_read_matlab_field(mat, 'target_coil_matsim', np.array)
        self.angle_limits = try_to_read_matlab_field(mat, 'angle_limits', list)
        self.handle_direction_ref = try_to_read_matlab_field(mat, 'handle_direction_ref', list)
        self.radius = try_to_read_matlab_field(mat, 'radius', float)
        self.resolution_angle = try_to_read_matlab_field(mat, 'resolution_angle', float)
        self.resolution_pos = try_to_read_matlab_field(mat, 'resolution_pos', float)

        self.save_fields = try_to_read_matlab_field(mat, 'save_fields', list)  # 'normE', 'E', both
        self.save_fields = [f.strip(' ') for f in self.save_fields]
        self.write_mesh_geom = try_to_read_matlab_field(mat, 'write_mesh_geom', bool)

        if self.compress_hdf5 == '':
            self.compress_hdf5 = None
        self.optim_name = try_to_read_matlab_field(mat, 'optim_name', str, '')

    def remove_poslist(self, number=None):
        """Removes the specified poslist

        Parameters:
        ------------------------------
        number: int
        indice of postist to be removed
        """

        if number:
            logger.warn("TMSOPTIMISATION.remove_poslist() does not eat number= argument.")
        self.optimlist = None
        self._prepared = False

    def clear_poslist(self):
        """ Removes all poslists
        """
        self.remove_poslist()
        self._prepared = False

    def sim_struct2mat(self):
        """ Creates a dictionary to save a matlab structure with scipy.io.savemat()

        Returns
        --------------------
        dict
            Dictionaty for usage with scipy.io.savemat
        """
        mat = {}
        mat['type'] = 'SESSION'
        mat['date'] = remove_None(self.date)
        mat['volfn'] = remove_None(self.volfn)
        mat['subpath'] = remove_None(self.subpath)
        mat['eeg_cap'] = remove_None(self.eeg_cap)
        mat['fnamehead'] = remove_None(self.fnamehead)
        mat['pathfem'] = remove_None(self.pathfem)
        mat['fname_tensor'] = remove_None(self.fname_tensor)
        mat['vol'] = self.vol.sim_struct2mat()
        mat['map_to_vol'] = remove_None(self.map_to_vol)
        mat['map_to_MNI'] = remove_None(self.map_to_MNI)
        mat['map_to_fsavg'] = remove_None(self.map_to_fsavg)
        mat['map_to_surf'] = remove_None(self.map_to_surf)
        mat['fields'] = remove_None(self.fields)
        mat['fiducials'] = self.fiducials.sim_struct2mat()

        mat['poslist'] = []
        # mat['cond'] = remove_None(self.cond)
        mat['hdf5_fn'] = remove_None(self.fname_hdf5)
        mat['didt'] = remove_None(self.didt)
        mat['distance'] = remove_None(self.distance)
        mat['fnamecoil'] = remove_None(self.fnamecoil)
        mat['qoi'] = remove_None(self.qoi)
        mat['save_fields'] = [remove_None(self.save_fields)]  # as list to cope with matlab_read
        mat['n_sim'] = remove_None(self.n_sim)
        mat['write_mesh_geom'] = remove_None(self.write_mesh_geom)
        mat['compress_hdf5'] = remove_None(self.compress_hdf5)
        mat['anisotropy_type'] = remove_None(self.anisotropy_type)
        mat['optim_name'] = remove_None(self.optim_name)

        if self.optimlist:
            mat['optimlist'] = self.optimlist.sim_struct2mat()

        # target information
        mat['angle_limits'] = remove_None(self.angle_limits)
        mat['handle_direction_ref'] = remove_None(self.handle_direction_ref)
        mat['radius'] = remove_None(self.radius)
        mat['resolution_angle'] = remove_None(self.resolution_angle)
        mat['resolution_pos'] = remove_None(self.resolution_pos)
        mat['target'] = remove_None(self.target)
        mat['target_coil_matsim'] = [remove_None(self.target_coil_matsim)]  # as list to cope with matlab_read
        mat['target_idx'] = remove_None(self.target_idx)

        return mat

    def _prepare(self):
        """
        Prepares session for optimizations.

        - relative paths are made absolute,
        - empty fields are set to default values,
        - check if required fields exist
        - check if TMSLIST parameters are set correctly (e.g. no different diDt within poslist, ...)
        - non-changing fields from the single poslist are stored in self


        """
        if self._prepared:
            logger.warn('TMSOPTIMIZATION already was prepared. Repeating all steps now.')
            self._prepared = False
        self.fnamehead = os.path.abspath(os.path.expanduser(self.fnamehead))
        if not os.path.isfile(self.fnamehead):
            raise IOError('Cannot locate head mesh file: %s' % self.fnamehead)

        sub_files = SubjectFiles(self.fnamehead, self.subpath)
        self.fnamehead = sub_files.fnamehead
        self.subpath = sub_files.subpath

        if not os.path.isdir(self.subpath):
            logger.warning('Cannot locate subjects m2m folder')
            logger.warning('some postprocessing options might fail')
            self.subpath = None

        if not self.fname_tensor:
            self.fname_tensor = sub_files.tensor_file

        self.fname_hdf5 = self.optim_name + '.hdf5'

        if self.eeg_cap:
            raise NotImplementedError('eeg_cap is not tested for TMS optimization.')

        logger.info('Head Mesh: {0}'.format(self.fnamehead))
        logger.info('Subject Path: {0}'.format(self.subpath))
        self.pathfem = os.path.abspath(os.path.expanduser(self.pathfem))
        logger.info('Simulation Folder: {0}'.format(self.pathfem))

        assert self.optimlist, "No TMSLIST as optimlist provided."
        if os.path.isfile(self.fnamehead):
            mesh = mesh_io.read_msh(self.fnamehead)
            mesh.fix_surface_labels()
            self.optimlist.postprocess = self.fields
            self.optimlist.fn_tensor_nifti = self.fname_tensor

            self.optimlist._prepare()
            if not self.optimlist.mesh:
                self.optimlist.mesh = mesh

            # do some checks
            assert len(self.optimlist.pos)
            assert np.all([pos.didt == self.optimlist.pos[0].didt for pos in self.optimlist.pos]), \
                "For TMSOPTIMIZATION, use the same didt for all positions."

            for pos in self.optimlist.pos:
                if pos.name:
                    pos.name = int(pos.name)
            assert np.all([not pos.name or type(pos.name) == int for pos in self.optimlist.pos]), \
                "Provide no positions names or integer vales"

            # if no pos names are given, use index
            if not self.optimlist.pos[0].name:
                for i, pos in enumerate(self.optimlist.pos):
                    pos.name = str(i)

            assert np.all([pos.distance == self.optimlist.pos[0].distance for pos in self.optimlist.pos]), \
                "For TMSOPTIMMIZATION, use the same coil distance for all positions."

            # store some fields in self as they don't change over the optimization positions
            self.didt = self.optimlist.pos[0].didt
            self.distance = self.optimlist.pos[0].distance
            self.mesh = self.optimlist.mesh
            self.fnamecoil = self.optimlist.fnamecoil
            self.n_sim = len(self.optimlist.pos)

            # store fields from TMSOPTIMIZATION in TMSLIST
            self.optimlist.anisotropy_type = self.anisotropy_type
        else:
            raise IOError(
                'Could not find head mesh file: {0}'.format(self.fnamehead))

        # hdf5 fn stuff
        if not self.fname_hdf5:
            self.fname_hdf5 = 'optim_results_{}.hdf5'.format(self.time_str)
        if not os.path.isabs(self.fname_hdf5):
            self.fname_hdf5 = os.path.join(self.pathfem, self.fname_hdf5)
        if not self.fname_hdf5.endswith('.hdf5'):
            self.fname_hdf5 += '.hdf5'

        # build list of fields to save from self.fields value
        if not self.save_fields:
            if self.fields != 'eE':
                logger.warn("Only 'eE' (that is normE and E) are tested for optimization.")
            self.save_fields = [field for field in self.fields]

            # replace e by normE
            try:
                idx = self.save_fields.index('e')
                self.save_fields[idx] = 'normE'
            except ValueError:
                pass

        # hdf5 file init
        if os.path.exists(self.fname_hdf5) and not self.resume:
            backup_fn = self.fname_hdf5[:-5] + '_' + self.time_str + '_backup' + '.hdf5'
            logger.warn("{} already exists. Moving to {}.".format(self.fname_hdf5, backup_fn))
            os.rename(self.fname_hdf5, backup_fn)
        if self.write_mesh_geom and not self.resume:
            self.mesh.write_hdf5(self.fname_hdf5, compression='gzip')

        # resuming unfinished optimization tasks
        if self.resume:
            self.optimlist = self.get_resume_optimlist()
        self._prepared = True

    def get_resume_optimlist(self):
        """
        Build TMSLIST of unfinished simulations.

        Returns
        --------
        TMSLIST
        """
        # get unfinished indices
        unfinished = self.get_unfinished_sim_idx()
        if not unfinished:
            logger.warn("No unfinished simulations found. Start with resume=False to start new optimization.")
        resume_tms_list = copy.copy(self.optimlist)
        resume_tms_list.pos = []
        for i in unfinished:
            resume_tms_list.add_position(self.optimlist.pos[i])

        return resume_tms_list

    def get_unfinished_sim_idx(self):
        """Checks self.hdf5_fn to find not processed simulations.

        Returns
        -------
        set
            Contains the .hdf5/data/field indices that are not processed.
            May be empty.
        """

        paths = ['/nodedata/', 'elmdata/']
        unfinished = set()
        try:
            with h5py.File(self.fname_hdf5, 'r') as f:
                for field, path in itertools.product(self.save_fields, paths):
                    try:
                        g = f[path + field]
                        unfinished = unfinished.union(np.where(g[-1, :, :] == 0)[1].tolist())
                    except KeyError:
                        pass
        except OSError:
            raise OSError("Cannot open {}. "
                          "Start optimization with resume=False to start from scratch.".format(self.fname_hdf5))

        if unfinished:
            logger.info("Found {} (out of {}) unfinished simulations.".format(len(unfinished), self.n_sim))
        return unfinished

    def run(self, cpus=1, allow_multiple_runs=True, save_mat=True):
        """ Run simulations in the current optimization.

        Parameters
        -----------
        cpus: int (optional)
            Number of cpus to use. Not necessarily will use all cpus. Default: 1
        allow_multiple_runs: bool (optinal)
            Whether to allow multiple runs in one folder. Default: False
        save_mat: bool (optional)
            Whether to save the .mat file of this structure

        Returns
        ---------
        Writes the optimization results.

        """
        if not self.optim_name:
            print("optim_name not set. Defaulting to 'optimization'")
            self.optim_name = "optimization"

        self._set_logger()
        if not self._prepared:
            self._prepare()
        dir_name_sim = os.path.abspath(os.path.expanduser(self.pathfem))

        if not self.optimlist.pos:
            logger.info("No simulations left.")
            return self.fname_hdf5

        # some directory checks
        if os.path.isdir(dir_name_sim):
            g = glob.glob(os.path.join(dir_name_sim, '*optimization*.mat'))
            if g and not allow_multiple_runs:
                raise IOError('Found already existing simulation results in directory.'
                              ' Please run the simulation in a new directory or delete'
                              ' the simnibs_simulation*.mat files from the folder : {0}'.format(dir_name_sim))

            logger.info('Saving optimization results in {0}'.format(dir_name_sim))
        else:
            logger.info('Running optimization results in new directory {0}'.format(dir_name_sim))
            os.makedirs(dir_name_sim)

        if save_mat:
            save_matlab_sim_struct(self,
                                   os.path.join(dir_name_sim,
                                                self.optim_name + '_{}.mat'.format(self.time_str)))

        self.cond2elmdata(force=True)

        cpus = np.min((cpus, self.n_sim))

        # set up solver
        solver = FEMSystem.tms(self.mesh, self.cond)
        solver.lock = multiprocessing.Lock()  # for the thread-save .hdf5 writes

        if cpus > multiprocessing.cpu_count():
            logger.warn(
                "{} parallel jobs requested, but only {} cores found.".format(cpus, multiprocessing.cpu_count()))

        def ec(e):
            """Error callback for pool.apply_async()"""
            raise e

        starttime = datetime.datetime.now()
        if cpus == 1:
            _set_up_global_solver(solver)
            # simulate each position
            for pos in self.optimlist.pos:
                run_optim_sim(pos,
                              self.mesh,
                              self.fnamecoil,
                              self.didt,
                              self.fname_hdf5,
                              self.fields,
                              self.cond,
                              self.save_fields,
                              self.n_sim,
                              self.compress_hdf5)
            _finalize_global_solver()

        else:
            with multiprocessing.Pool(processes=cpus,
                                      initializer=_set_up_global_solver,
                                      initargs=(solver,)) as pool:
                sims = []
                for pos in self.optimlist.pos:
                    sims.append(pool.apply_async(run_optim_sim, (pos,
                                                                 self.mesh,
                                                                 self.fnamecoil,
                                                                 self.didt,
                                                                 self.fname_hdf5,
                                                                 self.fields,
                                                                 self.cond,
                                                                 self.save_fields,
                                                                 self.n_sim,
                                                                 self.compress_hdf5), error_callback=ec))
                pool.close()
                pool.join()

        duration = datetime.datetime.now() - starttime
        logger.info('=====================================')
        logger.info("SimNIBS finished running simulations: #{}, {} ({}/simulation)".format(self.n_sim, duration,
                                                                                           duration / self.n_sim))

        logger.info('Simulation results stored in {}'.format(self.fname_hdf5))
        return self.fname_hdf5

    def cond2elmdata(self, force=False):
        """
        Wrapper method for optimlist.cond2elmdata().

        Parameters:
        -----------
        force: bool (false)
            Set force to recompute conductivity information, even if already present.
        """
        assert isinstance(self.optimlist, TMSLIST), "No poslist present."
        if not self.cond or force:
            self.cond = self.optimlist.cond2elmdata()
        else:
            logger.info("Reusing already processed conductivity information.")


def run_optim_sim(pos, mesh, fnamecoil, didt, hdf5_fn, field_names, cond, fields_to_save, n_sim, compression=None):
    """
    Runs a single simulation for TMSOPTIMIZATION

    Parameters:
    -----------
    pos: simnibs.simulation.mesh_io.POSITION
        Position object with matsimnibs.
    mesh: simnibs.simulation.mesh_io.MSH
        Head mesh.
    fnamecoil: str
        Coil filename. Endswith .ccd or .nii or .nii.gz
    didt: int
        Intensity used for simulation.
    hdf5_fn: str
        Filename of .hdf5 file. Path must exist, file may.
    field_names: list of str
        Field types to be computed. E.g. ['eE']
    cond
    save_fields: list of str
        Field names that shall be saved in hdf5. Subset of fields. E.g. ['eE']
    n_sim: int
        Number of simulations in total.
    compression: str or None (Default: None)
        Compression used for .hdf5. 'gzip' or None.
    """
    global tms_global_solver

    # compute dadt
    dadt = coil.set_up_tms_dAdt(
        mesh,
        fnamecoil,
        pos.matsimnibs,
        didt=didt)
    if isinstance(dadt, mesh_io.NodeData):
        dadt = dadt.node_data2elm_data()  # 5s
    dadt.field_name = 'dAdt'

    # compute fields
    b = tms_global_solver.assemble_tms_rhs(dadt)  # 9s
    v = tms_global_solver.solve(b)
    v = mesh_io.NodeData(v, name='v', mesh=mesh)
    mesh.nodedata = [v]
    v.mesh = mesh
    computed_fields = calc_fields(v, field_names, cond=cond, dadt=dadt)  # 10s

    # write results thread-save into hdf5
    tms_global_solver.lock.acquire()
    add_optim_fields_to_hdf5(computed_fields, pos,
                             hdf5_fn=hdf5_fn,
                             fields_to_save=fields_to_save,
                             n_sim=n_sim,
                             compression=compression)  # 5s

    tms_global_solver.lock.release()

    del v


def add_optim_fields_to_hdf5(data, pos, hdf5_fn, fields_to_save, n_sim, compression=None):
    """
    Adds a single simulation (fields and position information) from optimization to hdf5 file.

    For each field in fields_to_save, one array is created with (field.shape,n_sim). The data of this simulation is
    stored in field[::,pos.name], so pos.name should be the simulation's index number.

    Parameters:
    -----------
    pos: simnibs.simulation.sim_struct.POSITION
        Position object that according to simulation
    data: simnibs.msh.mesh_io.Msh
        Msh with simnibs.msh.mesh_io.ElementData in data.fields[]
    hdf5_fn: string
        Filename of .hdf5 file. Defaults to self.hdf5_fn
    fields_to_save: list of string
        Fields to pick from data.fields
    compression: str or None (Default: None)
        Compression used for .hdf5. 'gzip' or None.
    """

    sim_nr = int(pos.name)
    logger.debug("Writing results of position {} to .hdf5 .".format(sim_nr))

    with h5py.File(hdf5_fn, 'a') as f:
        # write elmdata
        for elmdata in data.elmdata:
            field = elmdata.field_name
            logger.debug("Starting field {} of position {} to .hdf5 done.".format(field, sim_nr))

            # only write out fields from fields_to_save
            if field in fields_to_save:
                try:
                    g = f.create_group('/elmdata')
                except ValueError:
                    g = f['/elmdata']

                if field not in g.keys():
                    # create initial zero-filled dataset with the size field X n_optimizations
                    try:
                        shape = (elmdata.value.shape[0],
                                 elmdata.value.shape[1],
                                 n_sim)  # 3-dimensional data (E)
                    except IndexError:
                        shape = (elmdata.value.shape[0],
                                 1,
                                 n_sim)  # 1-dimensional data (normE)
                    g.create_dataset(field, shape, compression=compression, chunks=shape[:2] + (1,))

                try:
                    # add this simulation's field to the array
                    g[field][:, :, sim_nr] = elmdata.value[:]
                except TypeError:
                    g[field][:, 0, sim_nr] = elmdata.value[:]

        # write nodedata
        for nodedata in data.nodedata:
            field = nodedata.field_name
            if field in fields_to_save:
                try:
                    g = f.create_group('/nodedata')
                except ValueError:
                    g = f['/nodedata']

                if field not in g.keys():
                    # create initial zero-filled dataset with the size field X n_optimizations
                    try:
                        shape = (nodedata.value.shape[0],
                                 nodedata.value.shape[1],
                                 n_sim)  # 3-dimensional data (E)
                    except IndexError:
                        shape = (nodedata.value.shape[0],
                                 1,
                                 n_sim)  # 1-dimensional data (normE)
                    g.create_dataset(field, shape, compression=compression, chunks=shape[:2] + (1,))

                try:
                    # add this simulation's field to the array
                    g[field][:, :, sim_nr] = nodedata.value[:]
                except TypeError:
                    g[field][:, 0, sim_nr] = nodedata.value[:]
        # add coil information
        try:
            g = f.create_group('/coil')
        except ValueError:
            g = f['/coil']
        if 'matsimnibs' not in g.keys():
            # create initial zero-filled dataset with the size field X n_optimizations
            shape = (4, 4, n_sim)
            g.create_dataset('matsimnibs', shape, compression="gzip")
        g['matsimnibs'][:, :, sim_nr] = pos.matsimnibs

        logger.debug("Writing results of position {} to .hdf5 done.".format(sim_nr))


def _set_up_global_solver(S):
    global tms_global_solver
    tms_global_solver = S


def _finalize_global_solver():
    global tms_global_solver
    del tms_global_solver

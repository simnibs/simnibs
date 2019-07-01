import copy
import datetime
import gc
import glob
import logging
import h5py
import os
import pandas as pd
import pickle
import numpy as np
import itertools
from simnibs import msh
from ..simulation import TMSLIST, sim_struct
from ..utils.simnibs_logger import logger
from simnibs.msh.mesh_io import _find_mesh_version, _read_msh_2, _read_msh_4, Msh, read_msh

try:
    import pyfempp
except ImportError:
    print("pyfempp not found")
    pyfempp = False


def eval_optim(simulations, target, tms_list, res_fn):
    """
    Evaluation of optimization results.

    Parameters
    ----------
    simulations: list of str
        list of filenames of simulations (.hdf5 files)

    target: list of float [3]
        cortical target position in subject space [x,y,z]

    tms_list: simnibs.sim_struct.TMSLIST object
        TMS simulation object

    res_fn: str
        filename for results .csv file
    """
    # load tms_list if filename
    if type(tms_list) == str:
        tms_list_new = TMSLIST()
        tms_list_new.read_mat_struct(tms_list)
        tms_list = tms_list_new

    # get coil distance
    assert all(pos.distance == tms_list.pos[0].distance for pos in tms_list.pos), "Different coil distances found."
    distance = tms_list.pos[0].distance

    assert len(simulations) == len(tms_list.pos), \
        'Simulation number ({}) does not fit to TMSLIST positions ({})'.format(len(simulations), len(tms_list.pos))

    # read first mesh to get index of target element
    first_res_mesh = Msh()
    first_res_mesh = first_res_mesh.read_hdf5(simulations[0])
    _, idx = first_res_mesh.find_closest_element(target, return_index=True)

    # read normE at from target element from each simulation
    # p = Pool(n_cpu)
    # f = partial(get_target_e_from_mesh, target=target)
    # f = partial(get_target_e_from_hdf5, idx=idx - 1, verbose=True)  # .hdf5 indexing is +1 to .msh

    # read fields in a pooled manner
    # target_e = p.map(f, simulations)

    # read data
    target_e = [get_target_e_from_hdf5(sim, idx=idx - 1, verbose=False) for sim in simulations]

    # build dictionary to easily create pandas afterwards
    d = dict()
    for i, simulation in enumerate(simulations):
        d[i] = [target_e[i]] + tms_list.pos[i].centre + tms_list.pos[i].pos_ydir + [simulation, idx - 1] + target
    data = pd.DataFrame.from_dict(d)
    data = data.transpose()
    data.columns = ['normE', 'x', 'y', 'z',
                    'handle_x', 'handle_y', 'handle_z',
                    'fn', 'idx_hdf5',
                    'target_x', 'target_y', 'target_z']
    data.to_csv(res_fn)  # save normE from all simulations at target elm to .csv

    # get best coil position/orientation from the results
    data['normE'] = data['normE'].astype(float)
    best_cond = data.iloc[[data['normE'].idxmax()]]
    tms_list = copy.copy(tms_list)
    tms_list.pos = []  # use original tms_list and replace pos with best conditions
    if not (best_cond.count(0) == 1).all():
        logger.warn("Multiple optimal coil pos/rot found!")

    for _, row in best_cond.iterrows():
        best_pos = sim_struct.POSITION()
        best_pos.distance = distance
        best_pos.centre = row[['x', 'y', 'z']].values.tolist()
        best_pos.pos_ydir = row[['handle_x', 'handle_y', 'handle_z']].values.tolist()

        tms_list.add_position(best_pos)

    return tms_list


def get_target_e_from_mesh(mesh_fn, target, na_val=-1, verbose=False):
    """
    Worker function to pool.map() .msh reading.

    Parameters
    ----------
    mesh_fn: basestring
        location of .msh file
    target: list of int
        target location to extract normE from
    verbose: bool (Default: False)
        print some verbosity information

    Returns
    -------
    float:
        normE from element nearest to target location

    """
    if verbose:
        print("Reading target E from {}".format(os.path.basename(mesh_fn)))
        t = datetime.datetime.now()
    try:
        sim_msh = read_msh(mesh_fn)
        elm, idx = sim_msh.find_closest_element(target, return_index=True)
        e = sim_msh.field['normE'][idx]
    except (OSError, IOError):
        logger.warn("Cannot read {}".format(mesh_fn))
        e = na_val
    gc.collect()
    if verbose:
        logger.info("Reading target E from {} done ({}s).".format(os.path.basename(mesh_fn),
                                                                  datetime.datetime.now() - t))
    return e


def get_target_e_from_hdf5(mesh_fn, idx, na_val=np.float64(-1), verbose=False):
    """
    Worker function to pool.map() .hdf5 reading.

    Remember, the hdf5 indexing is +1 towards .msh indexing

    Parameters
    ----------
    mesh_fn: basestring
        location of .msh file
    idx: int
        element index to extract normE from
    verbose: bool (Default: False)
        print some verbosity information

    Returns
    -------
    float:
        normE from idx element
    """
    try:
        e = h5py.File(mesh_fn, 'r')['/elmdata/normE'][idx]
    except (OSError, IOError):
        logger.warn("Cannot read {}".format(mesh_fn))
        e = na_val
    if verbose:
        logger.info("Reading target E from {} done.".format(os.path.basename(mesh_fn)))
    return e


def read_msh_from_pckl(fn, m=None):
    """ Reads a gmsh '.msh' file

    Parameters
    ------------
    fn: str
        File name
    m: simnibs.msh.Msh (optional)
        Mesh structure to be overwritten. If unset, will create a new structure

    Returns
    --------
    msh: simnibs.msh.Msh
        Mesh structure
    """
    print("This is the monkey-patched version of read_msh().")
    assert fn.endswith('.msh')
    fn_pckl = fn[:-3] + "pckl"
    try:
        return pickle.load(open(fn_pckl, 'rb'))

    except FileNotFoundError:
        print("Falling back to original read_msh for {}".format(fn))
        if m is None:
            m = msh.Msh()

        fn = os.path.expanduser(fn)

        if not os.path.isfile(fn):
            raise IOError(fn + ' not found')

        version_number = _find_mesh_version(fn)
        if version_number == 2:
            m = _read_msh_2(fn, m)

        elif version_number == 4:
            m = _read_msh_4(fn, m)

        else:
            raise IOError('Unrecgnized Mesh file version : {}'.format(version_number))

        return m


def get_opt_grid(tms_list, mesh, target, handle_direction_ref, distance=1., radius=20, resolution_pos=1, resolution_angle=10,
                 angle_limits=None):
    """
    Determine the coil positions and orientations for bruteforce TMS optimization

    Parameters
    ----------
    tms_list: simnibs.simulation.sim_struct.TMSLIST object
        TMS simulation object instance
    mesh: simnibs.msh.mesh_io.Msh object
        Simnibs mesh object
    target: np.ndarray or list of float
        Coordinates (x, y, z) of cortical target
    handle_direction_ref: list of float or np.ndarray
        Vector of handle prolongation direction
    distance: float (Default: 1.)
        Coil distance to skin surface [mm].
    radius: float
        Radius of region of interest around skin-projected cortical target, where the bruteforce simulations are
        conducted
    resolution_pos: float
        Resolution in mm of the coil positions in the region of interest
    resolution_angle: float (Default: 20)
        Resolution in deg of the coil positions in the region of interest
    angle_limits: list of float (Default: [-30, 30])
        Range of angles to get coil rotations for.

    Returns
    -------
    tms: simnibs.simulation.sim_struct.TMSLIST object
        TMS simulation object instance containing the coil positions and orientations to run bruteforce simulations
    """
    if not angle_limits:
        angle_limits = [-30, 30]

    # project cortical target to skin surface
    tms_tmp = copy.deepcopy(tms_list)
    pos_target = tms_tmp.add_position()
    pos_target.centre = target
    pos_target.pos_ydir = handle_direction_ref
    pos_target.distance = .1
    target_matsimnibs = pos_target.calc_matsimnibs(mesh, log=False)
    target_skin = target_matsimnibs[0:3, 3]

    # extract ROI
    msh_surf = mesh.crop_mesh(elm_type=2)
    msh_skin = msh_surf.crop_mesh([5, 1005])
    elm_center = np.mean(msh_skin.nodes.node_coord[msh_skin.elm.node_number_list[:, 0:3] - 1], axis=1)
    elm_mask_roi = np.linalg.norm(elm_center - target_skin, axis=1) < 1.2 * radius
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
    coords_plane = np.array(np.meshgrid(np.linspace(-radius, radius, int(2 * radius / resolution_pos + 1)),
                                        np.linspace(-radius, radius, int(2 * radius / resolution_pos + 1)))).T.reshape(
        -1, 2)
    coords_plane = coords_plane[np.linalg.norm(coords_plane, axis=1) <= radius]
    coords_plane = np.dot(coords_plane, vh[:, :2].transpose()) + target_skin

    # project grid-points to skin-surface
    p1 = msh_skin.nodes.node_coord[node_number_list_roi[:, 0],]
    p2 = msh_skin.nodes.node_coord[node_number_list_roi[:, 1],]
    p3 = msh_skin.nodes.node_coord[node_number_list_roi[:, 2],]

    normals = np.cross(p2 - p1, p3 - p1)
    coords_mapped = np.zeros(coords_plane.shape)

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

        coords_mapped[i,] = p0[inside,]

    # determine rotation matrices around z-axis of coil and rotate
    angles = np.linspace(angle_limits[0],
                         angle_limits[1],
                         int((angle_limits[1] - angle_limits[0]) / resolution_angle + 1))
    handle_directions = np.zeros((len(angles), 3))

    for i, a in enumerate(angles):
        x = target_matsimnibs[0:3, 0]
        y = target_matsimnibs[0:3, 1]
        handle_directions[is] = np.cos(a / 180. * np.pi) * y + np.sin(a / 180. * np.pi) * -x

    # combine coil positions and orientations
    po = list(itertools.product(coords_mapped, handle_directions))
    n_pos = len(tms_list.pos)
    logger.info("Determining {} coil positions/rotations for optimization.".format(n_pos))

    # write coordinates in TMS object
    for i, val in enumerate(po):
        c, h = val
        pos = tms_list.add_position()
        pos.centre = c.tolist()
        pos.pos_ydir = h.tolist() + c.tolist()
        pos.distance = 0.

        if not i % 100:  # on every 100th interation print status
            if i:
                duration = datetime.datetime.now() - starttime
                try:
                    time_left = ((n_pos - i) / 100.) * duration
                    time_left = datetime.timedelta(seconds=time_left.seconds)
                except OverflowError:
                    time_left = 'unknown'
                logger.info("Creating matsimnibs {0:0>4}/{1} (time left: {2}).".format(i, n_pos, time_left))
            else:
                logger.info("Creating matsimnibs {0:0>4}/{1}.".format(i, n_pos))
            starttime = datetime.datetime.now()

        # pos.pos_ydir = pos.pos_ydir
        # pos.centre = pos.centre
        pos.matsimnibs = None
        pos.distance = distance
        pos.matsimnibs = pos.calc_matsimnibs(mesh, log=False, msh_surf=msh_surf)

    # get matsimnibs for each position

    return tms_list


def optimize_tms_coil_pos(session, target=None,
                          handle_direction_ref=None, radius=20, angle_limits=None,
                          resolution_pos=1.5, resolution_angle=15, distance=1.,
                          n_cpu=8):
    """Compute optimal TMS coil position/rotation for maximum field at cortical target.

    This implements a brute-force approach to find the optimal coil position & rotation to maximize normE
    at a given cortical target. The target is projected to the nearest position at the skin surface and candidate coil
    positions/rotations are then distributed equally in a circular plane around that location.
    The number of candiate coil positions/rotations is defined by the two resolution parameters (_pos and
    angle), the radius of the circular plane, and angle_limits.

    All simulations are stored in opt_folder/simulations/ as .hdf5 files, and results.csv .

    Parameters
    ----------
    target: np.ndarray or list of float
        XYZ coordinates to optimize coil position/rotation for
    session: simnibs.simulation.stim_struct.SESSION
        Session object with 1 TMSLIST() as poslist.
        If session.poslists[0].pos is empty, is is filled according to parameters below.
    handle_direction_ref: list of num or np.ndarray (Default: [-2.8, 7, 8.1])
        Vector of coil handle prolongation direction. Defaults to left M1 45°.
    radius: float (Default: 20)
        Radius of region of interest around skin-projected cortical target, where the bruteforce simulations are
        conducted
    resolution_pos: float (Default: 1.5)
        Resolution in mm of the coil positions in the region of interest
    resolution_angle: float (Default: 15)
        Resolution in deg of the coil positions in the region of interest
    angle_limits: list of int or np.ndarray (Default: [-60, 60])
        Range of angles to get coil rotations for. Default: [-30, 30]
    distance: float (Default: 1.)
        Coil distance to scalp
    n_cpu: Int
        Number of cores to use

    Returns
    -------
    files
        all simulations are stored in session.pathfem/simulations/ (.hdf, .geo)

    file
        normE at target element for all simulations is saved to opt_folder/optim.csv

    file
        session object used for simulation in opt_folder/optim_session.mat

    file
        tmslist object with optimal coil positions

    dict
        'best_conds': simnibs.sim_struct.TMSLIST for optimal coil position(s)
        'simulations': simulation filenames,
        'session': simnibs.sim_struct.SESSION object used for simulation
    """

    # some argument checks, defaults, type handling
    if angle_limits is None:
        angle_limits = [-60, 60]
    if handle_direction_ref is None:
        handle_direction_ref = [-2.8, 7, 8.1]

    assert len(session.poslists) == 1, "Please provide session with 1 TMS poslist."

    if type(target) == np.ndarray:
        target = target.tolist()
    if type(handle_direction_ref) == np.ndarray:
        handle_direction_ref = handle_direction_ref.tolist()
    if type(resolution_angle) == np.ndarray:
        resolution_angle = resolution_angle.tolist()
    if type(angle_limits) == np.ndarray:
        angle_limits = angle_limits.tolist()

    tms_list = session.poslists[0]
    if tms_list.pos and (
            target or handle_direction_ref or resolution_angle or angle_limits or distance):
        logger.warn("tms_list positions provided. Ignoring other provided arguments!")

    if not tms_list.pos and not (target):
        raise ValueError("Provide either target or tms_list.pos .")

    opt_folder = session.pathfem
    sim_folder = opt_folder + "/simulations/"
    session.pathfem = sim_folder

    if not os.path.exists(sim_folder):
        logger.info("Creating folder {}".format(sim_folder))
        os.makedirs(sim_folder)
    if glob.glob(sim_folder + '*.msh'):
        raise FileExistsError(".msh files present in {}. Quitting.".format(sim_folder))

    # setup tmslist if not provided
    if not tms_list.pos:

        # mesh = mesh_io.read_msh(mesh)  # TODO: remove this
        # mesh.fix_surface_labels()
        mesh = pickle.load(open(session.fnamehead[:-3] + "pckl", 'rb'))

        tms_list = get_opt_grid(tms_list=tms_list,
                                mesh=mesh,
                                target=target,
                                handle_direction_ref=handle_direction_ref,
                                radius=radius,
                                angle_limits=angle_limits,
                                resolution_pos=resolution_pos,
                                resolution_angle=resolution_angle)

        # get msh_surf to speed up things a bit


    else:
        logger.info("Optimization: using positions from provided TMSLIST.pos .")

    if pyfempp:
        pyfempp.create_stimsite_from_tmslist(opt_folder + "/coil_positions.hdf5",
                                             tms_list, overwrite=True)

    tms_list.create_visuals = False
    tms_list.remove_msh = True
    tms_list.open_in_gmsh = False

    # update session poslist with changed tms_list
    session.poslists[0] = tms_list

    # save session with positions to disk
    tmslist_optim_all = opt_folder + "/optim_session.mat"
    sim_struct.save_matlab_sim_struct(session, tmslist_optim_all)

    # run simulations
    logger.info("Starting optimization: {} FEMs on {} cpus".format(len(tms_list.pos), n_cpu))
    out = session.run(allow_multiple_runs=True, cpus=n_cpu)

    simulations = [fn[:-3] + 'hdf5' for fn in out]

    # evaluate all simulations
    tms_list_optim = eval_optim(simulations, target, tms_list, opt_folder + '/optim.csv')

    # save session with positions to disk
    tmslist_optim_all = opt_folder + "/optim_best_conds_tmslist.mat"
    sim_struct.save_matlab_sim_struct(tms_list_optim, tmslist_optim_all)

    # return best position
    logging.shutdown()

    return {'best_conds': tms_list_optim,
            'simulations': simulations,
            'session': session}

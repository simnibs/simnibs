import copy
import datetime
import gc
import glob
import logging
import h5py
import os
import time
from functools import partial
from multiprocessing.pool import Pool
import pandas as pd
import pickle
import numpy as np
import itertools
import pyfempp
from simnibs import msh
from ..simulation import TMSLIST, sim_struct, save_matlab_sim_struct
from ..utils.simnibs_logger import logger
from simnibs.msh.mesh_io import _find_mesh_version, _read_msh_2, _read_msh_4, Msh, read_msh


def get_target_e_from_mesh(mesh_fn, target, verbose=False):
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
    print("Reading target E from {}".format(os.path.basename(mesh_fn)))
    t = datetime.datetime.now()
    sim_msh = read_msh(mesh_fn)
    elm, idx = sim_msh.find_closest_element(target, return_index=True)
    e = sim_msh.field['normE'][idx]
    gc.collect()
    if verbose:
        logger.info("Reading target E from {} done ({}s).".format(os.path.basename(mesh_fn),
                                                                  datetime.datetime.now() - t))
    return e


def get_target_e_from_hdf5(mesh_fn, idx, verbose=False):
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
    e = h5py.File(mesh_fn, 'r')['/elmdata/normE'][idx]
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


def get_opt_grid(tms, msh, target, handle_direction_ref, radius=20, resolution_pos=1, resolution_angle=10,
                 angle_limits=None):
    """
    Determine the coil positions and orientations for bruteforce TMS optimization

    Parameters
    ----------
    tms: simnibs.simulation.sim_struct.TMSLIST object
        TMS simulation object instance
    msh: simnibs.msh.mesh_io.Msh object
        Simnibs mesh object
    target: np.ndarray or list of float
        Coordinates (x, y, z) of cortical target
    handle_direction_ref: list of float or np.ndarray
        Vector of handle prolongation direction
    radius: float
        Radius of region of interest around skin-projected cortical target, where the bruteforce simulations are
        conducted
    resolution_pos: float
        Resolution in mm of the coil positions in the region of interest
    resolution_angle: int or float
        Resolution in deg of the coil positions in the region of interest
    angle_limits: list of int
        Range of angles to get coil rotations for. Default: [-30, 30]

    Returns
    -------
    tms: simnibs.simulation.sim_struct.TMSLIST object
        TMS simulation object instance containing the coil positions and orientations to run bruteforce simulations
    """
    if not angle_limits:
        angle_limits = [-30, 30]

    # project cortical target to skin surface
    tms_tmp = copy.deepcopy(tms)
    pos_target = tms_tmp.add_position()
    pos_target.centre = target
    pos_target.pos_ydir = [0, 0, 0]
    pos_target.distance = .1
    target_matsimnibs = pos_target.calc_matsimnibs(msh, log=False)
    target_skin = target_matsimnibs[0:3, 3]

    # extract ROI
    msh_surf = msh.crop_mesh(elm_type=2)
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
    p1 = msh_skin.nodes.node_coord[node_number_list_roi[:, 0], ]
    p2 = msh_skin.nodes.node_coord[node_number_list_roi[:, 1], ]
    p3 = msh_skin.nodes.node_coord[node_number_list_roi[:, 2], ]

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

        coords_mapped[i, ] = p0[inside, ]

    # determine handle directions
    # project handle_direction_ref to reference plane
    # handle_direction_ref_proj = np.dot(handle_direction_ref, vh[:,:2])
    # handle_directions = np.tile(np.array(handle_direction_ref), (coords_mapped.shape[0], 1))

    # determine rotation matrices around z-axis of coil and rotate
    angles = np.linspace(angle_limits[0],
                         angle_limits[1],
                         int((angle_limits[1] - angle_limits[0]) / resolution_angle + 1))
    handle_directions = np.zeros((len(angles), 3))

    for i, a in enumerate(angles):
        mat_rot = np.array([[np.cos(a / 180. * np.pi), -np.sin(a / 180. * np.pi), 0],
                            [np.sin(a / 180. * np.pi), np.cos(a / 180. * np.pi), 0],
                            [0, 0, 1]])
        handle_directions[i, ] = np.dot(handle_direction_ref, mat_rot)

    # combine coil positions and orientations
    po = list(itertools.product(coords_mapped, handle_directions))

    # write coordinates in TMS object
    for c, h in po:
        pos = tms.add_position()
        pos.centre = c.tolist()
        pos.pos_ydir = h.tolist()
        pos.distance = 0.

    return tms

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(nodes_roi[:, 0], nodes_roi[:, 1], nodes_roi[:, 2])
    # ax.quiver(target_skin[0], target_skin[1], target_skin[2],
    #           5*vh[0, 0], 5*vh[1, 0], 5*vh[2, 0], color='k')
    # ax.quiver(target_skin[0], target_skin[1], target_skin[2],
    #           5*vh[0, 1], 5*vh[1, 1], 5*vh[2, 1], color='r')
    # ax.quiver(target_skin[0], target_skin[1], target_skin[2],
    #           5 * vh[0, 2], 5 * vh[1, 2], 5 * vh[2, 2], color='g')
    #
    # ax.quiver(c[0], c[1], c[2], Q1[0], Q1[1], Q1[2], color='k')
    # ax.quiver(c[0], c[1], c[2], Q2[0], Q2[1], Q2[2], color='g')
    #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    # ax.quiver(c[0], c[1], c[2], 1e2 * vh[0, 2], 1e2 * vh[1, 2], 1e2 * vh[2, 2], color='k')
    # ax.quiver(c[0], c[1], c[2], -1e2 * vh[0, 2], -1e2 * vh[1, 2], -1e2 * vh[2, 2], color='kk')
    # ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], color="r")
    # ax.scatter(coords_mapped[:, 0], coords_mapped[:, 1], coords_mapped[:, 2], color="g")
    #
    #
    #
    # ax.scatter(P0[inside, ][:, 0], P0[inside, ][:, 1], P0[inside, ][:, 2], color="g")
    # ax.scatter(elm_center_roi[inside, ][:, 0],
    #            elm_center_roi[inside, ][:, 1], elm_center_roi[inside, ][:, 2], color="g")
    # ax.scatter(nodes_roi[:, 0], nodes_roi[:, 1], nodes_roi[:, 2])

    # target_skin_center, target_skin_elm_idx = msh_skin.find_closest_element(target_skin, return_index=True)
    # target_skin_nodes = msh_skin.nodes.node_coord[msh_skin.elm.node_number_list[target_skin_elm_idx][0:3]-1, ]
    # target_skin_tangents = np.array([target_skin_nodes[1, ] - target_skin_nodes[0, ],
    #                                  target_skin_nodes[2, ] - target_skin_nodes[0, ]]).transpose()


def optimize_tms_coil_pos(session,
                          target=None, coil_fn=None, tms_list=None,
                          handle_direction_ref=None, radius=20, angle_limits=None,
                          resolution_pos=1.5, resolution_angle=15, distance=1.,
                          n_cpu=8):
    """Compute optimal TMS coil position/rotation for maximum field at cortical target.

    This implements a brute-force approach to find the optimal coil position & rotation to maximize normE
    at given target. The cortical target is projected to the nearest position at the skin surface and candidate coil
    positions/rotations are then distributed equally in a circular plance around that location.
    The number of candiate coil positions/rotations is defined by the two resolution parameters (_pos and
    angle), the radius of the circular plane, and angle_limits.

    All simulations are stored in opt_folder/simulations/ as .hdf5 files, and results.csv .

    Parameters
    ----------
    target: np.ndarray or list of float
        XYZ coordinates to optimize coil position/rotation for
    session: simnibs.simulation.stim_struct.SESSION
        Session object
    tms_list: simnibs.stimulation.sim_struct.TMSLIST (optional)
        Filled tmslist to simulate. If provided, arguments below are ignored.
    coil_fn: basestring
        Location of coil nifti-file
    handle_direction_ref: list of num or np.ndarray
        Vector of handle prolongation direction
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
        all simulations are stored in opt_folder/simulations/ (.hdf, .geo)

    file
        normE at target element for all simulations is saved to opt_folder/optim.csv

    file
        tms_poslist.pckl in opt_folder/

    file
        tmslist for best coil position/orientation(s) is saved to opt_folder/optim.mat

    simnibs.sim_struct.TMSLIST
        tmslist for best coil position(s)
    """

    # some argument checks, defaults, type handling
    if angle_limits is None:
        angle_limits = [-60, 60]
    if handle_direction_ref is None:
        handle_direction_ref = [22.5, 45, 67.5]

    if type(target) == np.ndarray:
        target = target.tolist()
    if type(handle_direction_ref) == np.ndarray:
        handle_direction_ref = handle_direction_ref.tolist()
    if type(resolution_angle) == np.ndarray:
        resolution_angle = resolution_angle.tolist()
    if type(angle_limits) == np.ndarray:
        angle_limits = angle_limits.tolist()

    if tms_list is not None and (
            target or coil_fn or handle_direction_ref or resolution_angle or angle_limits or distance):
        logger.warn("tms_list provided. Ignoring other provided arguments!")

    if not tms_list and not (target and coil_fn):
        raise ValueError("Provide either tms_list or (target, coil_fn, and handle_direction_ref).")

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

        logger.info("Determining coil positions and orientations for optimization.")

        tms_list = get_opt_grid(tms=tms_list,
                                msh=mesh,
                                target=target,
                                handle_direction_ref=handle_direction_ref,
                                radius=radius,
                                angle_limits=angle_limits,
                                resolution_pos=resolution_pos)

        # get msh_surf to speed up things a bit
        msh_surf = mesh.crop_mesh(elm_type=2)

        # get matsimnibs for each position
        n_pos = len(tms_list.pos)
        logger.info("Determining {} coil positions/rotations for optimization.".format(n_pos))
        for i, pos in enumerate(tms_list.pos):  # pos = tms.pos[0]

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

            pos.pos_ydir = pos.pos_ydir
            pos.centre = pos.centre
            pos.matsimnibs = None
            pos.distance = distance
            pos.matsimnibs = pos.calc_matsimnibs(mesh, log=False, msh_surf=msh_surf)

        # store poslist as matfile
        save_matlab_sim_struct(opt_folder + '/tms_poslist.mat')
        # pickle.dump(tms_list, open(opt_folder + '/tms_poslist.pckl', 'wb'))
    else:
        logger.info("Optimization: Calculation positions from provided TMSLIST.pos.")

    # TODO: remove pyfempp
    pyfempp.create_stimsite_from_tmslist(opt_folder + "/coil_positions.hdf5",
                                         tms_list, overwrite=True)
    tms_list.create_visuals = False
    tms_list.remove_msh = True
    tms_list.open_in_gmsh = False

    session.add_poslist(tms_list)

    # run simulations
    logger.info("Starting optimization: {} FEMs on {} cpus".format(len(tms_list.pos), n_cpu))
    out = session.run(allow_multiple_runs=True, cpus=n_cpu)

    simulations = [fn[:-3] + 'hdf5' for fn in out]

    # read first mesh to get index of target element
    first_res_mesh = Msh()
    first_res_mesh = first_res_mesh.read_hdf5(simulations[0])
    _, idx = first_res_mesh.find_closest_element(target, return_index=True)

    # read normE at from target element from each simulation
    p = Pool(n_cpu)
    # f = partial(get_target_e_from_mesh, target=target)
    f = partial(get_target_e_from_hdf5, idx=idx-1, verbose=True)  # .hdf5 indexing is +1 to .msh

    # read fields in a pooled manner
    target_e = p.map(f, simulations)

    # build dictionary to easily create pandas afterwards
    d = dict()
    for i, simulation in enumerate(simulations):
        d[i] = [target_e[i]] + tms_list.pos[i].centre + tms_list.pos[i].pos_ydir
    data = pd.DataFrame.from_dict(d)
    data = data.transpose()
    data.columns = ['normE', 'x', 'y', 'z', 'handle_x', 'handle_y', 'handle_z']
    data.to_csv(opt_folder + '/optim.csv')  # save normE from all simulations at target elm to .csv

    # get best coil position/orientation from the results
    best_cond = data.iloc[[data['normE'].idxmax()]]
    tms_list.pos = []  # use original tms_list and replace pos with best conditions
    if not (best_cond.count(0) == 1).all():
        logger.warn("Multiple optimal coil pos/rot found!")

    for _, row in best_cond.iterrows():
        best_pos = sim_struct.POSITION()
        best_pos.distance = distance
        best_pos.centre = row[['x', 'y', 'z']].values.tolist()
        best_pos.pos_ydir = row[['handle_x', 'handle_y', 'handle_z']].values.tolist()

        tms_list.add_position(best_pos)

    # save results TMSLIST to disk
    sim_struct.save_matlab_sim_struct(tms_list, opt_folder + "/optim.mat")

    # return best position
    logging.shutdown()

    return tms_list

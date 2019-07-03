import os
import copy
import datetime
import logging
import h5py
import pickle
import itertools
import pandas as pd
import numpy as np
import simnibs

from .. import msh
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
        tms_list_new = simnibs.simulation.TMSOPTIMIZATION()
        tms_list_new.read_mat_struct(tms_optim)
        tms_optim = tms_list_new

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
            'Simulation number ({}) does not fit to TMSLIST positions ({})'.format(len(simulations), len(tms_optim.pos))
        assert res.shape[-1] > 1, 'Only one simulation found. No optimization possible.'

        # read first mesh to get index of target element
        try:
            # try to read mesh geom from .hdf5
            mesh = msh.mesh_io.Msh()
            mesh = mesh.read_hdf5(simulations, load_data=False)
            _, idx = mesh.find_closest_element(target, return_index=True)
        except (ValueError, KeyError):
            mesh = msh.mesh_io.read_msh(tms_optim.fnamehead)
            _, idx = mesh.find_closest_element(target, return_index=True)

        qoi_in_target = np.linalg.norm(res[idx - 1], axis=0)
        assert qoi_in_target.shape == (len(tms_optim.optimlist.pos),)

        # build dictionary to easily create pandas afterwards
        d = dict()
        for i in range(qoi_in_target.shape[-1]):
            d[i] = [i, qoi_in_target[i]] + tms_optim.optimlist.pos[i].centre + tms_optim.optimlist.pos[i].pos_ydir + [
                idx - 1] + target
        data = pd.DataFrame.from_dict(d)
        data = data.transpose()
        data.columns = ['sim_nr'] + [qoi] + ['x', 'y', 'z'] + \
                       ['handle_x', 'handle_y', 'handle_z',
                        'idx_hdf5',
                        'target_x', 'target_y', 'target_z']
        data.to_csv(res_fn, index=False)  # save normE from all simulations at target elm to .csv

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

            tms_optim.optimlist.add_position(best_pos)

    return tms_optim


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

        version_number = simnibs.msh.msh.mesh_io_find_mesh_version(fn)
        if version_number == 2:
            m = simnibs.simulation.mesh_io._read_msh_2(fn, m)

        elif version_number == 4:
            m = simnibs.simulation.mesh_io._read_msh_4(fn, m)

        else:
            raise IOError('Unrecgnized Mesh file version : {}'.format(version_number))

        return m


def get_opt_grid(tms_list, mesh, target, handle_direction_ref, distance=1., radius=20, resolution_pos=1,
                 resolution_angle=10,
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
    coords_plane = np.array(np.meshgrid(np.linspace(-radius,
                                                    radius,
                                                    int(2 * radius / resolution_pos + 1)),
                                        np.linspace(-radius,
                                                    radius,
                                                    int(2 * radius / resolution_pos + 1)))).T.reshape(-1, 2)
    coords_plane = coords_plane[np.linalg.norm(coords_plane, axis=1) <= radius]
    coords_plane = np.dot(coords_plane, vh[:, :2].transpose()) + target_skin

    # project grid-points to skin-surface
    p1 = msh_skin.nodes.node_coord[node_number_list_roi[:, 0]]
    p2 = msh_skin.nodes.node_coord[node_number_list_roi[:, 1]]
    p3 = msh_skin.nodes.node_coord[node_number_list_roi[:, 2]]

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

        coords_mapped[i] = p0[inside]

    # determine rotation matrices around z-axis of coil and rotate
    angles = np.linspace(angle_limits[0],
                         angle_limits[1],
                         int((angle_limits[1] - angle_limits[0]) / resolution_angle + 1))
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

    # write coordinates in TMS object
    for i, val in enumerate(po):
        c, h = val
        pos = tms_list.add_position()
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

    return tms_list


def optimize_tms_coil_pos(tms_optim, target=None,
                          handle_direction_ref=None, radius=20, angle_limits=None,
                          resolution_pos=1.5, resolution_angle=15, distance=1., n_cpu=8):
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
    tms_optim: simnibs.simulation.stim_struct.SESSION
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

    assert tms_optim.optimlist, "Please provide simnibs.simulation.sim_struct.TMSOPTIMIZATION with TMSLIST optimlist."

    if type(target) == np.ndarray:
        target = target.tolist()
    if type(handle_direction_ref) == np.ndarray:
        handle_direction_ref = handle_direction_ref.tolist()
    if type(resolution_angle) == np.ndarray:
        resolution_angle = resolution_angle.tolist()
    if type(angle_limits) == np.ndarray:
        angle_limits = angle_limits.tolist()

    # tms_list = session.poslists[0]
    if tms_optim.optimlist.pos and (
            target or handle_direction_ref or resolution_angle or angle_limits or distance):
        simnibs.utils.simnibs_logger.logger.warn("tms_list positions provided. Ignoring other provided arguments!")

    if not tms_optim.optimlist.pos and not target:
        raise ValueError("Provide either target or tms_list.pos .")

    # setup tmslist if not provided
    if not tms_optim.optimlist.pos:

        # mesh = msh.mesh_io.read_msh(mesh)  # TODO: remove this
        # mesh.fix_surface_labels()
        mesh = pickle.load(open(tms_optim.fnamehead[:-3] + "pckl", 'rb'))

        tms_list = get_opt_grid(tms_list=tms_optim.optimlist,
                                mesh=mesh,
                                target=target,
                                handle_direction_ref=handle_direction_ref,
                                radius=radius,
                                angle_limits=angle_limits,
                                resolution_pos=resolution_pos,
                                resolution_angle=resolution_angle)

        # update session poslist with changed tms_list
        tms_optim.add_poslist(tms_list)

    else:
        simnibs.utils.simnibs_logger.logger.info("Optimization: using positions from provided TMSLIST.pos .")

    if pyfempp:
        pyfempp.create_stimsite_from_tmslist(tms_optim.pathfem + "/coil_positions.hdf5",
                                             tms_optim.optimlist, overwrite=True)

    # run simulations
    simnibs.utils.simnibs_logger.logger.info("Starting optimization: "
                                             "{} FEMs on {} cpus".format(len(tms_optim.optimlist.pos), n_cpu))
    simulations_fn = tms_optim.run(allow_multiple_runs=True, cpus=n_cpu)

    # evaluate all simulations
    tms_list_optim = eval_optim(simulations_fn,
                                target,
                                tms_optim,
                                tms_optim.pathfem + '/optim_' + tms_optim.time_str + '.csv',
                                qoi=tms_optim.qoi)

    # return best position
    logging.shutdown()

    return {'best_conds': tms_list_optim,
            'simulations': simulations_fn,
            'tms_optim': tms_optim}

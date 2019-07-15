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
from ..simulation import save_matlab_sim_struct

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

        if not tms_optim.target or tms_optim.target != target:
            # read first mesh to get index of target element
            try:
                # try to read mesh geom from .hdf5
                mesh = msh.mesh_io.Msh()
                mesh = mesh.read_hdf5(simulations, load_data=False)
                _, idx = mesh.find_closest_element(target, return_index=True)
            except (ValueError, KeyError):
                mesh = msh.mesh_io.read_msh(tms_optim.fnamehead)
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
        Simnibs mesh object
    target: list of float or np.ndarray or None
        Coordinates (x, y, z) of cortical target
    handle_direction_ref: list of float or np.ndarray
        Vector of handle prolongation direction
    distance: float or None
        Coil distance to skin surface [mm]. (Default: 1.)
    radius: float or None
        Radius of region of interest around skin-projected cortical target, where the bruteforce simulations are
        conducted
    resolution_pos: float or None
        Resolution in mm of the coil positions in the region of interest.
    resolution_angle: float or None
        Resolution in deg of the coil positions in the region of interest (Default: 20)
    angle_limits: list of float or None
        Range of angles to get coil rotations for (Default: [-60, 60]).

    Returns
    -------
    tms_optim: simnibs.simulation.sim_struct.TMSOPTIMIZATION object
        TMS simulation object instance
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

        mesh = msh.mesh_io.read_msh(tms_optim.fnamehead)  # TODO: remove this
        mesh.fix_surface_labels()
        # mesh = pickle.load(open(tms_optim.fnamehead[:-3] + "pckl", 'rb'))

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

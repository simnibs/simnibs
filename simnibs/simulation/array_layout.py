import os
import copy
import numpy as np
import simnibs
from . current_estimator import CurrentEstimator
from . sim_struct import SESSION
from .. utils.ellipsoid import subject2ellipsoid, ellipsoid2subject
from .. mesh_tools.surface import Surface
from .. utils.ellipsoid import Ellipsoid
from .. utils.file_finder import Templates


class ElectrodeInitializer():
    """
    Helper class to catch electrode properties before initialization
    """
    def __init__(self):
        """
        Initializes ElectrodeInitializer class instance
        """

        # general properties
        self.type = None        # "ElectrodeArrayPair" or "CircularArray"
        self.current = None

        # CircularArray properties
        self.radius_inner = None
        self.radius_inner_bounds = None
        self.radius_outer = None
        self.radius_outer_bounds = None
        self.distance = None
        self.distance_bounds = None
        self.n_outer = None
        self.n_outer_bounds = None

        # ElectrodeArrayPair properties
        self._center = None
        self._radius = None
        self.radius_bounds = None
        self._length_x = None
        self.length_x_bounds = None
        self._length_y = None
        self.length_y_bounds = None

        # simulation properties
        self.dirichlet_correction = None
        self.dirichlet_correction_detailed = None
        self.current_estimator_method = None

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        if type(value) is int or type(value) is float:
            value = [value]

        self._radius = value

    @property
    def length_x(self):
        return self._length_x

    @length_x.setter
    def length_x(self, value):
        if type(value) is int or type(value) is float:
            value = [value]

        self._length_x = value

    @property
    def length_y(self):
        return self._length_y

    @length_y.setter
    def length_y(self, value):
        if type(value) is int or type(value) is float:
            value = [value]

        self._length_y = value

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, value):
        if type(value) is list:
            value = np.array(value)

        if value.ndim == 1:
            value = value[np.newaxis, :]

        self._center = value

    def initialize(self):
        """
        Run electrode initialization and return final electrode

        Returns
        -------
        electrode : ElectrodeArrayPair or CircularArray class instance
            Electrode
        """
        if self.type == "CircularArray":
            electrode = CircularArray(radius_inner=self.radius_inner,
                                      distance=self.distance,
                                      n_outer=self.n_outer,
                                      radius_outer=self.radius_outer,
                                      current_estimator_method=self.current_estimator_method,
                                      dirichlet_correction=self.dirichlet_correction,
                                      dirichlet_correction_detailed=self.dirichlet_correction_detailed,
                                      current=self.current,
                                      distance_bounds=self.distance_bounds,
                                      radius_inner_bounds=self.radius_inner_bounds,
                                      radius_outer_bounds=self.radius_outer_bounds,
                                      n_outer_bounds=self.n_outer_bounds)

        elif self.type == "ElectrodeArrayPair":
            electrode = ElectrodeArrayPair(center=self.center,
                                           radius=self.radius,
                                           length_x=self.length_x,
                                           length_y=self.length_y,
                                           dirichlet_correction=self.dirichlet_correction,
                                           dirichlet_correction_detailed=self.dirichlet_correction_detailed,
                                           current=self.current,
                                           radius_bounds=self.radius_bounds,
                                           length_x_bounds=self.length_x_bounds,
                                           length_y_bounds=self.length_y_bounds)

        else:
            electrode = None

        return electrode


class ElectrodeMaster():
    """
    ElectrodeMaster class containing attributes and methods for ElectrodeArrayPair and CircularArray
    """

    def __init__(self, dirichlet_correction, dirichlet_correction_detailed, current_outlier_correction,
                 current_estimator):
        """
        Initializes ElectrodeMaster class
        """
        # global nodal arrays (set by compile_node_arrays)
        self.node_channel_id = None
        self.node_array_id = None
        self.node_ele_id = None
        self.node_coords = None
        self.node_idx = None
        self.node_current = None
        self.node_current_sign = None
        self.node_voltage = None
        self.node_voltage_sign = None
        self.node_area = None

        self.electrode_arrays = None

        self.dirichlet_correction = dirichlet_correction
        self.dirichlet_correction_detailed = dirichlet_correction_detailed
        self.current_outlier_correction = current_outlier_correction
        self.current_estimator = current_estimator

    def estimate_currents(self, electrode_pos):
        """
        Estimate electrode currents for fake Dirichlet BC using the CurrentEstimator class.
        Writes current in self.current.

        Parameters
        ----------
        electrode_pos : list of np.ndarray of float [n_array_free][n_pos x 3]
            Positions and orientations of ElectrodeArrayPair or CircularArray
        """
        if self.current_estimator.current is None:
            return
        return self.current_estimator.estimate_current(electrode_pos=np.hstack(electrode_pos))

    def compile_node_arrays(self):
        """
        Gathers all information from the nodes of the electrode arrays and the containing electrodes
        and collect them in global arrays
        """
        _node_channel_id = []
        _node_array_id = []
        _node_ele_id = []
        _node_coords = []
        _node_idx = []
        _node_current = []
        _node_current_sign = []
        _node_voltage = []
        _node_voltage_sign = []
        _node_area = []

        for i_array, _electrode_array in enumerate(self.electrode_arrays):
            for i_ele, _ele in enumerate(_electrode_array.electrodes):
                if _ele.node_coords is None:
                    self.node_channel_id = None
                    self.node_array_id = None
                    self.node_ele_id = None
                    self.node_coords = None
                    self.node_idx = None
                    self.node_current = None
                    self.node_current_sign = None
                    self.node_voltage = None
                    self.node_voltage_sign = None
                    self.node_area = None
                    return

                _node_channel_id.append(_ele.channel_id * np.ones(_ele.n_nodes))
                _node_array_id.append(i_array * np.ones(_ele.n_nodes))
                _node_ele_id.append(_ele.ele_id * np.ones(_ele.n_nodes))
                _node_coords.append(_ele.node_coords)
                _node_idx.append(_ele.node_idx)
                _node_current.append(_ele.node_current)
                _node_voltage.append(_ele.node_voltage)
                _node_area.append(_ele.node_area)

                if _ele.ele_voltage_sign is not None:
                    _node_voltage_sign.append(_ele.ele_voltage_sign * np.ones(_ele.n_nodes))
                else:
                    _node_voltage_sign = None

                if _ele.ele_current_sign is not None:
                    _node_current_sign.append(_ele.ele_current_sign * np.ones(_ele.n_nodes))
                else:
                    _node_current_sign = None

        self.node_channel_id = np.hstack(_node_channel_id)
        self.node_array_id = np.hstack(_node_array_id)
        self.node_ele_id = np.hstack(_node_ele_id)
        self.node_coords = np.vstack(_node_coords)
        self.node_idx = np.hstack(_node_idx)
        self.node_current = np.hstack(_node_current)
        self.node_voltage = np.hstack(_node_voltage)
        self.node_area = np.hstack(_node_area)

        if _node_voltage_sign is not None:
            self.node_voltage_sign = np.hstack(_node_voltage_sign)

        if _node_current_sign is not None:
            self.node_current_sign = np.hstack(_node_current_sign)

    def update_electrode_from_posmat(self, mesh, fn_electrode_mask=None):
        """
        Takes information from the posmat matrices and writes global and local node array information to electrode
        instances of the electrode arrays.

        Parameters
        ----------
        mesh : Msh Object
            Head mesh
        fn_electrode_mask : str, optional, default: Templates().mni_volume_upper_head_mask
            Filename of mask to use to define valid skin region, where the electrodes can be applied
        """
        # get some paths
        subpath = os.path.split(mesh.fn)[0]

        if fn_electrode_mask is None:
            ff_templates = Templates()
            fn_electrode_mask = ff_templates.mni_volume_upper_head_mask

        # relabel internal air
        mesh_relabel = simnibs.opt_struct.relabel_internal_air(m=mesh,
                                                               subpath=subpath,
                                                               label_skin=1005,
                                                               label_new=1099,
                                                               label_internal_air=501)

        # make final skin surface including some additional distance
        skin_surface_ellipsoid = Surface(mesh=mesh_relabel, labels=1005)
        skin_surface = copy.deepcopy(skin_surface_ellipsoid)
        skin_surface_ellipsoid = simnibs.opt_struct.valid_skin_region(skin_surface=skin_surface_ellipsoid,
                                                                      fn_electrode_mask=fn_electrode_mask,
                                                                      mesh=mesh_relabel,
                                                                      additional_distance=0)

        # get mapping between skin_surface node indices and global mesh nodes
        node_idx_msh = np.where(np.isin(mesh.nodes.node_coord, skin_surface.nodes).all(axis=1))[0]

        # fit optimal ellipsoid to valid skin points
        ellipsoid = Ellipsoid()
        ellipsoid.fit(points=skin_surface_ellipsoid.nodes)

        alpha_array = []
        start_array = []
        distance_array = []

        for i_array, _electrode_array in enumerate(self.electrode_arrays):
            # transform electrode position in posmat from subject space to ellipsoid space
            coords_base_eli_eli = subject2ellipsoid(coords=_electrode_array.posmat[:3, 3],
                                                    normals=_electrode_array.posmat[:3, 2],
                                                    ellipsoid=ellipsoid)

            # transform electrode center position to cartesian coordinates
            coords_base_eli_cart = ellipsoid.ellipsoid2cartesian(coords=coords_base_eli_eli,
                                                                 norm=False,
                                                                 return_normal=False)

            q1 = ellipsoid.get_geodesic_destination(start=np.tile(coords_base_eli_cart, (3,1)),
                                                    distance=np.array([0.1, 0.1, 0.1]),
                                                    alpha=np.array([0, -10/180. * np.pi, 10/180. * np.pi]),
                                                    n_steps=10)
            q = q1[0, :] - coords_base_eli_cart
            q /= np.linalg.norm(q)
            q = q.flatten()

            if (np.linalg.norm(q1[1, :] - (coords_base_eli_cart + _electrode_array.posmat[:3, 1])) <
                    np.linalg.norm(q1[2, :] - (coords_base_eli_cart + _electrode_array.posmat[:3, 1]))):
                q_sign = -1
            else:
                q_sign = +1

            # calculate angle between vector of constant lambda and electrode direction
            alpha_array.append(q_sign * np.arccos(np.dot(_electrode_array.posmat[:3, 1], q)) + _electrode_array.angle)

            # save starting point (center of array on ellipsoid)
            start_array.append(coords_base_eli_cart)

            # save electrode distances
            distance_array.append(_electrode_array.distance)

        electrode_array_idx = np.hstack([i_array * np.ones(_electrode_array.n_ele)
                                         for i_array, _electrode_array in
                                         enumerate(self.electrode_arrays)])

        # combine data from different electrode_arrays to run geodesic distance calculation only once
        start = np.vstack([np.tile(start_array[i_array], (_electrode_array.n_ele, 1))
                           for i_array, _electrode_array in enumerate(self.electrode_arrays)])

        alpha = np.hstack(alpha_array)
        distance = np.hstack(distance_array)

        for i_a, _alpha in enumerate(alpha):
            if _alpha > np.pi:
                alpha[i_a] = _alpha - 2 * np.pi
            elif _alpha < -np.pi:
                alpha[i_a] = _alpha + 2 * np.pi

        if not (distance == 0.).all():
            electrode_coords_eli_cart = ellipsoid.get_geodesic_destination(start=start,
                                                                           distance=distance,
                                                                           alpha=alpha,
                                                                           n_steps=400)
        else:
            electrode_coords_eli_cart = start

        n = ellipsoid.get_normal(coords=electrode_coords_eli_cart)

        # transform to ellipsoidal coordinates
        electrode_coords_eli_eli = ellipsoid.cartesian2ellipsoid(coords=electrode_coords_eli_cart)

        # project coordinates to subject
        tmp_arrays = []
        i_ele = 0
        for i_array, _electrode_array in enumerate(self.electrode_arrays):
            ele_idx, tmp = ellipsoid2subject(coords=electrode_coords_eli_eli[electrode_array_idx == i_array, :],
                                             ellipsoid=ellipsoid,
                                             surface=skin_surface)
            tmp_arrays.append(tmp)

            if len(ele_idx) != len(alpha[electrode_array_idx == i_array]):
                raise AssertionError("Electrode position invalid (can not project all electrodes on skin surface)")

            electrode_coords_subject = np.vstack(tmp_arrays)

            # loop over electrodes and determine node indices
            for _electrode in _electrode_array.electrodes:
                if _electrode.type == "spherical":
                    # mask with a sphere
                    mask = np.linalg.norm(
                        skin_surface.nodes - electrode_coords_subject[i_ele, :],
                        axis=1) < _electrode.radius

                    # save position of electrode in subject space to posmat field
                    _electrode.posmat[:3, 3] = electrode_coords_subject[i_ele, :]

                elif _electrode.type == "rectangular":
                    cx_local = np.cross(n[i_ele, :], _electrode_array.posmat[:3, 1])

                    # rotate skin nodes to normalized electrode space
                    rotmat = np.array([[cx_local[0], _electrode_array.posmat[0, 1], n[i_ele, 0]],
                                       [cx_local[1], _electrode_array.posmat[1, 1], n[i_ele, 1]],
                                       [cx_local[2], _electrode_array.posmat[2, 1], n[i_ele, 2]]])
                    center = np.array([electrode_coords_subject[i_ele, 0],
                                       electrode_coords_subject[i_ele, 1],
                                       electrode_coords_subject[i_ele, 2]])

                    # save position of electrode in subject space to posmat field
                    _electrode.posmat = np.vstack(
                        (np.hstack((rotmat, center[:, np.newaxis])), np.array([0, 0, 0, 1])))

                    skin_nodes_rotated = (skin_surface.nodes - center) @ rotmat

                    # mask with a box
                    mask_x = np.logical_and(skin_nodes_rotated[:, 0] > -_electrode.length_x / 2,
                                            skin_nodes_rotated[:, 0] < +_electrode.length_x / 2)
                    mask_y = np.logical_and(skin_nodes_rotated[:, 1] > -_electrode.length_y / 2,
                                            skin_nodes_rotated[:, 1] < +_electrode.length_y / 2)
                    mask_z = np.logical_and(skin_nodes_rotated[:, 2] > -30,
                                            skin_nodes_rotated[:, 2] < +30)
                    mask = np.logical_and(np.logical_and(mask_x, mask_y), mask_z)
                else:
                    raise AssertionError("Electrodes have to be either 'spherical' or 'rectangular'")

                # node areas
                _electrode.node_area = skin_surface.nodes_areas[mask]

                # total effective area of all nodes
                _electrode.area_skin = _electrode.node_area_total

                # electrode position is invalid if it overlaps with invalid skin region and area is not "complete"
                # if _electrode.area_skin < 0.90 * _electrode.area:
                #     raise AssertionError("Electrode position invalid (partly overlaps with invalid skin region)")

                # save node indices (referring to global mesh)
                _electrode.node_idx = node_idx_msh[mask]

                # save node coords (refering to global mesh)
                _electrode.node_coords = skin_surface.nodes[mask]

                i_ele += 1

        self.compile_node_arrays()

    def update_electrode_from_node_arrays(self):
        """
        Updates information from the node arrays and writes information to the electrode instances
        of the electrode arrays (inverse of compile_node_arrays)
        """
        for _electrode_array in self.electrode_arrays:
            for _ele in _electrode_array.electrodes:
                mask = (_ele.ele_id == self.node_ele_id) * (_ele.channel_id == self.node_channel_id)

                _ele.node_idx = self.node_idx[mask]
                _ele.node_area = self.node_area[mask]
                _ele.node_coords = self.node_coords[mask]

                if self.node_voltage[0] is not None:
                    _ele.node_voltage = self.node_voltage[mask]

                if self.node_current[0] is not None:
                    _ele.node_current = self.node_current[mask]

    def export_node_coords(self, fn_out):
        """
        Export node coordinates and node currents on subject skin surface to .txt file.
        The first 3 columns are the x, y, and z coordinates and the last column is the coil current.

        Parameters
        ----------
        fn_out : str
            Filename of output .txt file.
        """
        np.savetxt(fn_out, np.hstack((self.node_coords, self.node_current[:, np.newaxis])))

    def apply_current_outlier_correction(self):
        """
        Perform current outlier correction. Modifies self.node_current in place.
        """
        for channel_id in np.unique(self.node_channel_id):
            mask = (channel_id == self.node_channel_id).flatten()
            node_current_abs = np.abs(self.node_current[mask])
            c_min = np.mean(node_current_abs) - 8 * np.std(node_current_abs)
            c_max = np.mean(node_current_abs) + 8 * np.std(node_current_abs)
            outlier_mask = np.logical_not(np.logical_and(c_min < node_current_abs, c_max > node_current_abs))
            outlier_idx = np.where(outlier_mask)[0]

            print(f"Current outlier correction: corrected {len(outlier_idx)} outlier(s) in channel {int(channel_id)}.")
            for _idx in outlier_idx:
                # find 8 next neighbors
                neighbor_idx = np.argsort(
                    np.linalg.norm(self.node_coords[mask, :] - self.node_coords[mask, :][_idx, :], axis=1))[1:9]

                # filter out outliers in the neighbors
                neighbor_idx_valid = np.array([i for i in neighbor_idx if i not in outlier_idx])

                # set new current to maximal current of one of the neighbors
                current_outlier = node_current_abs[_idx]

                if len(neighbor_idx_valid) > 0:
                    current_new = np.max(node_current_abs[neighbor_idx_valid])
                else:
                    current_new = np.mean(node_current_abs)

                node_current_abs[_idx] = current_new

                # distribute the remaining current to all other nodes
                current_diff = current_outlier - current_new
                node_current_abs += current_diff * self.node_area[mask] / np.sum(self.node_area[mask])

            # write new current with correct sign into array again
            self.node_current[mask] = self.node_current_sign[mask] * node_current_abs

        self.update_electrode_from_node_arrays()


class Electrode():
    """
    Electrode class.
    Contains a single circular or rectangular electrode.
    ElectrodeArrays are built from Electrode instances.

    Parameters
    ----------
    channel_id : int
        Channel identifier, indicates connected electrodes (same identifier)
    ele_id : int
        Unique electrode identifier in whole setup
    center : np.ndarray of float [3]
        Center of electrode (in mm)
    radius : float, optional, default: None
        Radius of circular electrode, None in case of rectangular electrodes (in mm)
    length_x : float, optional, default: None
        Electrode extension in x-direction (rectangular electrode), None in case of circular electrodes (in mm)
    length_y : float, optional, default: None
        Electrode extension in y-direction (rectangular electrode), None in case of circular electrodes (in mm)
    ele_current : float, optional, default: None
        Current assigned to electrode
    ele_voltage : float, optional, default: None
        Voltage assigned to electrode
    node_current : np.ndarray of float [n_points]
        Associated current of nodes
    node_voltage : float, optional, default: None
        Associated voltages of nodes
    node_area : np.ndarray of float [n_points]
        Associated area of the nodes
    node_idx : np.ndarray of int [n_nodes]
        Node indices assigned to electrode on subject skin surface (referring to global mesh)
    node_coords : np.ndarray of int [n_nodes x 3]
        Node coordinates assigned to electrode on subject skin surface

    Attributes
    ----------
    channel_id : int
        Channel identifier, indicates connected electrodes (same identifier)
    center : np.ndarray of float [3]
        Center of electrode
    radius : float or None
        Radius of circular electrode
    length_x : float or None
        Electrode extension in x-direction (rectangular electrode)
    length_y : float or None
        Electrode extension in y-direction (rectangular electrode)
    posmat_norm : np.ndarray of float [4x4]
        Position matrix [4x4] of electrode in normalized space (in mm)
        np.array([1, 0, 0, center[0],
                  0, 1, 0, center[1],
                  0, 0, 1, center[2],
                  0, 0, 0, 1])
    posmat : np.ndarray of float [4x4]
        Position matrix [4x4] of electrode in head coordinate system (in mm)
    area : float
        Electrode area (in mm²)
    node_area : np.ndarray of float [n_points]
        Associated area of the nodes
    node_current : np.ndarray of float [n_points]
        Associated current of nodes
    node_idx : np.ndarray of int [n_nodes]
        Node indices assigned to electrode on subject skin surface (referring to global mesh)
    node_coords : np.ndarray of int [n_nodes x 3]
        Node coordinates assigned to electrode on subject skin surface
    node_voltage : float, optional, default: None
        Associated voltages of nodes
    type : str
        Type of electrode ("spherical", "rectangular")
    """

    def __init__(self, channel_id, ele_id, center, radius=None, length_x=None, length_y=None,
                 ele_current=None, ele_voltage=None,
                 node_current=None, node_voltage=None, node_area=None, node_idx=None, node_coords=None):

        if (ele_voltage is not None and ele_current is not None) or \
                (node_voltage is not None and node_current is not None):
            raise AssertionError("Define either voltage or current to electrode.")

        if ele_current is not None and node_current is not None:
            raise AssertionError("Define current either node wise or electrode wise")

        if ele_voltage is not None and node_voltage is not None:
            raise AssertionError("Define voltage either node wise or electrode wise")

        # geometrical properties
        ################################################################################################################
        if radius is None:
            radius = 0

        if length_x is None:
            length_x = 0

        if length_y is None:
            length_y = 0

        self.node_area = node_area
        self.node_coords = node_coords
        self.node_idx = node_idx
        self.channel_id = channel_id
        self.ele_id = ele_id
        self.center = center
        self.radius = radius
        self.length_x = length_x
        self.length_y = length_y
        self.transmat = None
        self.posmat_norm = np.array([[1, 0, 0, center[0]],
                                     [0, 1, 0, center[1]],
                                     [0, 0, 1, 0],
                                     [0, 0, 0, 1]])
        self.posmat = copy.deepcopy(self.posmat_norm)

        if radius != 0 and (length_x != 0 or length_y != 0):
            raise AssertionError("Define either radius for circular electrode or "
                                 "length_x and length_y for rectangular electrode")

        if self.radius is not None and radius != 0:
            self.area = np.pi*radius**2
            self.type = "spherical"
        else:
            self.area = length_x * length_y
            self.type = "rectangular"

        # source across whole electrode
        ################################################################################################################
        # current
        self.ele_current = ele_current
        self.ele_current_init = ele_current

        if self.ele_current is not None:
            self.ele_current_sign = np.sign(self.ele_current)
        else:
            self.ele_current_sign = None

        # voltage
        self.ele_voltage = ele_voltage
        self.ele_voltage_init = ele_voltage

        if self.ele_voltage is not None:
            self.ele_voltage_sign = np.sign(self.ele_voltage)
        else:
            self.ele_voltage_sign = None

        # source for individual nodes
        ################################################################################################################
        # raise error if only parts of node properties are provided
        if 1 <= np.sum([p is None for p in [node_area, node_coords, node_idx]]) < 3:
            raise AssertionError("Provide node_area, node_coords and node_idx when initializing an "
                                 "electrode on the node level")

        # current
        if node_current is not None:
            self.node_current_sign = np.sign(node_current)
            self.node_current = node_current
            self.node_current_init = node_current
        else:
            self.node_current_sign = None
            self.node_current = None
            self.node_current_init = None

        # voltage
        if node_voltage is not None:
            self.node_voltage_sign = np.sign(node_voltage)
            self.node_voltage = node_voltage
            self.node_voltage_init = node_voltage
        else:
            self.node_voltage_sign = None
            self.node_voltage = None
            self.node_voltage_init = None

    @property
    def ele_current(self):
        return self._ele_current

    @ele_current.setter
    def ele_current(self, value):
        self._ele_current = value

        # set node currents according to the node area
        if self.node_area is not None:
            self._node_current = value * self.node_area / self.node_area_total

    @property
    def node_current(self):
        return self._node_current

    @node_current.setter
    def node_current(self, value):
        if value is not None:
            assert len(value) == self.n_nodes, "Number of node currents does not match total number of nodes!"

            self._node_current = value

            # set ele current according to the node area
            if self.node_area is not None and len(self.node_area) == len(value):
                self._ele_current = np.sum(value)
        else:
            self._node_current = None

    @property
    def ele_voltage(self):
        return self._ele_voltage

    @ele_voltage.setter
    def ele_voltage(self, value):
        self._ele_voltage = value

        # set also node voltages
        if self.node_area is not None:
            self._node_voltage = value * np.ones(self.n_nodes)

    @property
    def node_voltage(self):
        return self._node_voltage

    @node_voltage.setter
    def node_voltage(self, value):
        if value is not None:
            assert len(value) == self.n_nodes, "Number of node voltages does not match total number of nodes!"

            self._node_voltage = value

            # also set ele voltage (average voltage over all nodes)
            self._ele_voltage = np.mean(value)  # * self.node_area / self.node_area_total)
        else:
            self._node_voltage = None

    @property
    def node_area(self):
        return self._node_area

    @node_area.setter
    def node_area(self, value):
        self._node_area = value

        # set also total node area and total number of nodes
        if value is not None:
            self.node_area_total = np.sum(value)
            self.n_nodes = len(value)

            # update node currents
            if self._node_current is None or len(self.node_current) != len(value):
                self._node_current = self.ele_current * self.node_area / self.node_area_total

    @property
    def node_coords(self):
        return self._node_coords

    @node_coords.setter
    def node_coords(self, value):
        self._node_coords = value

        # set also total number of nodes
        if value is not None:
            self.n_nodes = value.shape[0]

    @property
    def node_idx(self):
        return self._node_idx

    @node_idx.setter
    def node_idx(self, value):
        self._node_idx = value

        # set also total number of nodes
        if value is not None:
            self.n_nodes = len(value)

    def transform(self, transmat):
        """
        Transforms the electrode position matrix self.posmat_norm by the given 4x4 transformation matrix
        and saves new position in self.posmat.

        Parameters
        ----------
        transmat : nparray of float [4x4]
            Transformation matrix to place the electrode (in mm)
        """
        self.transmat = transmat
        self.posmat = self.posmat_norm @ transmat


class ElectrodeArray():
    """ Array of physically connected electrodes. Can only be moved together. Contains multiple Electrode instances.

    Circular and rectangular electrodes can be mixed.

    Example: Array with 2 circular electrode of radius 1 and 2 and one rectangular electrode of size 2x2
    radius   = [1,    None, 2   ]
    length_x = [None, 2,    None]
    length_y = [None, 2,    None]

    Parameters
    ----------
    channel_id : np.ndarray of int [n_ele]
        Channel identifiers, indicates connected electrodes (same identifier), in case of mixed definitions (circular
        and rectangular), the IDs of the circular electrodes are coming first followed by the rectangular ones.
    ele_id : np.ndarray of int [n_ele]
        Unique electrode identifier of whole setup
    center : np.ndarray of float [n_ele, 2]
        Center coordinates (x, y) of electrodes in normalized space (in mm)
    radius : np.ndarray of float [n_ele_circ] or None
        Radii of circular electrodes (in mm)
    length_x : np.ndarray of float [n_ele_rect] or None
        Electrode extensions in x-direction (rectangular electrodes)
    length_y : np.ndarray of float [n_ele_rect] or None
        Electrode extensions in y-direction (rectangular electrodes)
    current : np.ndarray of float [n_ele]
        Current through electrodes (in A)

    Attributes
    ----------
    channel_id : np.ndarray of int [n_ele]
        Channel identifiers, indicates connected electrodes (same identifier)
    array_center : np.ndarray of float [3]
        Center of the whole electrode array
    center : np.ndarray of float [n_ele, 3]
        Center of the individual electrodes (in mm)
    radius : np.ndarray of float [n_ele_circ]
        Radii of circular electrodes (in mm)
    length_x : np.ndarray of float [n_ele_rect] or None
        Electrode extensions in x-direction (rectangular electrodes)
    length_y : np.ndarray of float [n_ele_rect] or None
        Electrode extensions in y-direction (rectangular electrodes)
    n_ele : int
        Total number of electrodes (n_ele_circ + n_ele_rect)
    electrodes : list of Electrode instances
        Electrodes in the array
    type : list of str
        List containing the type of electrodes ("circ" and/or "rect")
    posmat_norm : np.ndarray of float [4x4]
        Position matrix [4x4] of electrode array in normalized space (in mm)
        np.array([1, 0, 0, center[0],
                  0, 1, 0, center[1],
                  0, 0, 1, center[2],
                  0, 0, 0, 1])
    posmat : np.ndarray of float [4x4]
        Position matrix [4x4] of electrode array in head coordinate system (in mm)
    distance : np.ndarray of float [n_ele]
        Euclidean distances of electrodes to array center (0, 0, 0)
    angle : np.ndarray of float [n_ele]
        Angle between connection line electrodes <-> array center and principal direction axis of array (0, 1, 0)
    electrode_pos : np.ndarray of float [3]
        Position and orientation of electrode array in ellipsoidal coordinates [beta, lambda, alpha]
    optimize_alpha : bool
        Flag to indicate if rotation angle of electrode array is also optimized. Generally yes but in case of single
        circular electrodes, the angle will not be optimized.
    """

    def __init__(self, channel_id, center, radius=None, length_x=None, length_y=None, current=None, ele_id=None):
        self.channel_id = channel_id
        self.array_center = np.array([0, 0])
        self.center = center
        self.radius = radius
        self.length_x = length_x
        self.length_y = length_y
        self.n_ele = len(channel_id)
        self.distance = np.zeros(self.n_ele)
        self.angle = np.zeros(self.n_ele)
        self.electrodes = []
        self.posmat_norm = np.array([[1, 0, 0, 0],
                                     [0, 1, 0, 0],
                                     [0, 0, 1, 0],
                                     [0, 0, 0, 1]])
        self.posmat = copy.deepcopy(self.posmat_norm)
        self.transmat = None
        self.electrode_pos = None

        if ele_id is None:
            self.ele_id = np.arange(self.n_ele)
        else:
            self.ele_id = ele_id

        if current is None:
            self.current = 1/self.n_ele * np.ones(self.n_ele)
        else:
            self.current = current

        assert center.shape[0] == self.n_ele, "Number of elements in 'center' does not match number of electrodes!"
        assert len(radius) == self.n_ele, "Number of elements in 'radius' does not match number of electrodes!"
        assert len(length_x) == self.n_ele, "Number of elements in 'length_x' does not match number of electrodes!"
        assert len(length_y) == self.n_ele, "Number of elements in 'length_y' does not match number of electrodes!"
        assert len(self.current) == self.n_ele, "Number of elements in 'current' does not match number of electrodes!"

        for i_ele in range(self.n_ele):
            self.distance[i_ele] = np.linalg.norm(center[i_ele, :] - self.array_center)

            if self.distance[i_ele] == 0:
                self.angle[i_ele] = 0
            else:
                if (center[i_ele, 0] > 0 and center[i_ele, 1] > 0) or (center[i_ele, 0] < 0 and center[i_ele, 1] > 0):
                    self.angle[i_ele] = np.arcsin(center[i_ele, 0] / self.distance[i_ele])

                elif center[i_ele, 0] < 0 and center[i_ele, 1] < 0:
                    self.angle[i_ele] = - np.pi - np.arcsin(center[i_ele, 0] / self.distance[i_ele])

                elif center[i_ele, 0] > 0 and center[i_ele, 1] < 0:
                    self.angle[i_ele] = np.pi - np.arcsin(center[i_ele, 0] / self.distance[i_ele])

                elif center[i_ele, 0] == 0 and center[i_ele, 1] > 0:
                    self.angle[i_ele] = 0

                elif center[i_ele, 0] == 0 and center[i_ele, 1] < 0:
                    self.angle[i_ele] = -np.pi

                elif center[i_ele, 0] > 0 and center[i_ele, 1] == 0:
                    self.angle[i_ele] = np.pi/2

                elif center[i_ele, 0] < 0 and center[i_ele, 1] == 0:
                    self.angle[i_ele] = -np.pi/2

                elif center[i_ele, 0] == 0 and center[i_ele, 1] == 0:
                    self.angle[i_ele] = 0

            self.electrodes.append(Electrode(channel_id=self.channel_id[i_ele],
                                             ele_id=self.ele_id[i_ele],
                                             center=self.center[i_ele, :],
                                             radius=self.radius[i_ele],
                                             length_x=self.length_x[i_ele],
                                             length_y=self.length_y[i_ele],
                                             ele_current=self.current[i_ele]))

        # check if only one spherical electrode is present, if yes, do not optimize rotation angle
        if len(self.electrodes) == 1 and self.electrodes[0].type == "spherical":
            self.optimize_alpha = False
        else:
            self.optimize_alpha = True

    def transform(self, transmat):
        """
        Transforms electrode array configuration according to given 4x4 transformation matrix

        Parameters
        ----------
        transmat : nparray of float [4x4]
            Transformation matrix to place the electrode
        """
        self.transmat = transmat

        for i in range(self.n_ele):
            self.electrodes[i].transform(transmat=transmat)

    def plot(self, fn_plot=None, show=True):
        """
        Plot electrode array

        Parameters
        ----------
        show : bool
            Show plot
        fn_plot : str
            Filename of output.png file
        """
        import matplotlib
        import matplotlib.pyplot as plt

        try:
            matplotlib.use('Qt5Agg')
        except:
            pass

        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']

        if len(colors) < np.max(self.channel_id):
            colors = colors * int(np.ceil(colors/np.max(self.channel_id)))

        plt.ioff()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        i_ele = 0

        for cen, cha, rad, l_x, l_y in zip(self.center, self.channel_id, self.radius, self.length_x, self.length_y):
            if rad != 0:
                obj = plt.Circle(cen[:2],
                                 self.radius[i_ele],
                                 facecolor=colors[cha], edgecolor='k')

            elif l_x != 0:
                obj = plt.Rectangle(xy=(cen[0]-self.length_x[i_ele]/2,
                                        cen[1]-self.length_y[i_ele]/2),
                                    width=self.length_x[i_ele],
                                    height=self.length_y[i_ele],
                                    facecolor=colors[cha], edgecolor='k')
            i_ele += 1

            ax.add_artist(obj)
            plt.annotate(f"Ch: {cha}", xy=cen[:2], xytext=(-12, 0), va='center',
            xycoords='data', textcoords='offset points')

        # plt.tight_layout()
        plt.grid()
        ax.set_axisbelow(True)
        ax.set_xlim((np.min(self.center-1.5*np.max(self.radius)), np.max(self.center+1.5*np.max(self.radius))))
        ax.set_ylim((np.min(self.center-1.5*np.max(self.radius)), np.max(self.center+1.5*np.max(self.radius))))
        ax.set_aspect('equal', 'box')
        ax.set_xlabel("x in mm")
        ax.set_ylabel("y in mm")

        if show:
            plt.show()

        if fn_plot:
            plt.savefig(fn_plot)

        plt.close()


class ElectrodeArrayPair(ElectrodeMaster):
    """
    Symmetric pair of electrode arrays with n_x * n_y electrodes each.
    Contains 2 ElectrodeArray instances with n_x * n_y Electrode instances each.
    Circular and rectangular electrodes can be mixed.

    Example: Array with 2 circular electrode of radius 1 and 2 and one rectangular electrode of size 2x2
    radius   = [1,    None, 2   ]
    length_x = [None, 2,    None]
    length_y = [None, 2,    None]

    Parameters
    ----------
    center : np.ndarray of float [n_ele/2 x 2]
        Center positions of electrodes in normalized x/y plane (z=0). Values will be copied for second array.
    radius : np.ndarray of float [n_ele/2] or None
        Radii of circular electrodes in array (in mm). Array will be copied for second electrode array.
    length_x : np.ndarray of float [n_ele/2] or None
        x-dimension of rectangular electrodes in array (in mm). Values will be copied for second array.
    length_y : np.ndarray of float [n_ele/2] or None
        y-dimension of rectangular electrodes in array (in mm). Values will be copied for second array.
    radius_bounds : list or np.ndarray of float [2]
        (min, max) values of electrode radius for optimization (applied to all electrodes in array).
    length_x_bounds : list or np.ndarray of float [2]
        (min, max) values of electrode x-dimension for optimization (applied to all electrodes in array).
    length_y_bounds : list or np.ndarray of float [2]
        (min, max) values of electrode y-dimension for optimization (applied to all electrodes in array).
    current : np.ndarray of float [n_ele] or None
        Current through electrodes. Has to sum up to zero net current.
    current_estimator_method : str, optional, default: "linear"
        Method to estimate the electrode currents:
        - "linear": linear regression
        - "gpc": generalized polynomial chaos
    dirichlet_correction : bool, optional, default: True
        If electrodes are connected to the same channel, they have to have the same voltage. This is ensured by
        setting this flag.
    dirichlet_correction_detailed : bool, optional, default: False
        Apply detailed Dirichlet correction such that every node current is optimized separately to match the equal
        voltage constraint of an electrode (recommended for large electrodes as in regular TES applications)
    current_outlier_correction : bool, optional, default: False
        Apply current outlier correction after node-wise dirichlet approximation (dirichlet_correction_detailed=True).

    Attributes
    ----------
    n_ele : int
        Total number of electrodes including all channels and arrays
    center : np.ndarray of float [n_ele/2, 3]
        Center of electrodes (in mm). Same for second array, so we save it only once.
    radius : list of float
        Radii of electrodes (in mm). Same for second array, so we save it only once.
    electrode_arrays : list of ElectrodeArray instances [2]
        Two ElectrodeArray instances
    """

    def __init__(self, center, radius=None, length_x=None, length_y=None, current=None,
                 radius_bounds=None, length_x_bounds=None, length_y_bounds=None,
                 current_estimator_method=None, dirichlet_correction=True, dirichlet_correction_detailed=False,
                 current_outlier_correction=False):

        if type(center) is list:
            center = np.array(center)

        if radius is None:
            radius = np.zeros(center.shape[0])

        if length_x is None:
            assert length_y is None, "length_x provided but not length_y"

        if length_y is None:
            assert length_x is None, "length_y provided but not length_x"

        if length_x is None and length_y is None:
            length_x = np.zeros(center.shape[0])
            length_y = np.zeros(center.shape[0])

        # check free and fixed parameters
        if radius_bounds is not None:
            assert len(radius_bounds) == 2, ("Dimension mismatch! "
                                             "Please provide [min, max] values for radius_bounds.")
            self.radius_free = True
            self.radius_bounds = radius_bounds
            self.radius = np.array([np.mean(radius_bounds)] * center.shape[0])
        else:
            self.radius_free = False
            self.radius_bounds = [0, 0]
            self.radius = radius

        if length_x_bounds is not None:
            assert len(length_x_bounds) == 2, ("Dimension mismatch! "
                                               "Please provide [min, max] values for length_x.")
            self.length_x_free = True
            self.length_x_bounds = length_x_bounds
            self.length_x = np.array([np.mean(length_x_bounds)] * center.shape[0])
        else:
            self.length_x_free = False
            self.length_x_bounds = [0, 0]
            self.length_x = length_x

        if length_y_bounds is not None:
            assert len(length_y_bounds) == 2, ("Dimension mismatch! "
                                               "Please provide [min, max] values for length_y.")
            self.length_y_free = True
            self.length_y_bounds = length_y_bounds
            self.length_y = np.array([np.mean(length_y_bounds)] * center.shape[0])
        else:
            self.length_y_free = False
            self.length_y_bounds = [0, 0]
            self.length_y = length_y

        self.center = center
        self.n_ele = center.shape[0] * 2                                                    # total number of electrodes
        self.channel_id = np.hstack([[i for _ in range(self.center.shape[0])] for i in range(2)])    # array of all channel_ids
        self.channel_id_unique = np.unique(self.channel_id)
        self.n_channel = len(self.channel_id_unique)
        self.ele_id = np.arange(self.n_ele)

        # order of free parameters
        self.free_geometry = np.array([self.radius_free,
                                       self.length_x_free,
                                       self.length_y_free])
        self.any_free_geometry = self.free_geometry.any()
        self.geo_para_bounds = np.vstack((self.radius_bounds,
                                          self.length_x_bounds,
                                          self.length_y_bounds))
        self.geo_para_mean = np.hstack((self.radius,
                                        self.length_x,
                                        self.length_y))

        if current is None:
            self.current = np.hstack((np.array([1/(self.n_ele/2.) for _ in range(int(self.n_ele/2))]),
                                      np.array([-1/(self.n_ele/2.) for _ in range(int(self.n_ele/2))])))
        else:
            self.current = np.array(current)

        if np.abs(np.sum(self.current)) > 1e-12:
            raise AssertionError("Please check electrode currents. They do not sum up to 0. (atol = 1e-12)")

        # number of electrodes per channel [n_channel]
        self.n_ele_per_channel = np.array([np.sum(self.channel_id == i) for i in np.unique(self.channel_id)])

        # total current entering domain (read from first channel)
        self.current_total = np.sum(self.current[self.channel_id == self.channel_id[0]])

        # determine mean current of electrodes for each channel [n_channel]
        self.current_mean = self.current_total / self.n_ele_per_channel
        self.current_mean[1] *= -1

        # total current of each channel (here we only have 2)
        self.current_channel = np.array([self.current_total, -self.current_total])

        # initialize current estimator for fake Dirichlet BC
        self.current_estimator_method = current_estimator_method
        if current_estimator_method is None or current_estimator_method == "" or (self.n_ele_per_channel == 1).all():
            current_estimator = None
        else:
            current_estimator = CurrentEstimator(method=current_estimator_method,
                                                 channel_id=self.channel_id,
                                                 ele_id=self.ele_id,
                                                 current_sign=np.sign(self.current),
                                                 current_total=self.current_total)

        super(ElectrodeArrayPair, self).__init__(dirichlet_correction=dirichlet_correction,
                                                 dirichlet_correction_detailed=dirichlet_correction_detailed,
                                                 current_outlier_correction=current_outlier_correction,
                                                 current_estimator=current_estimator)

        # create two ElectrodeArray instance where all electrodes of the first array have channel_id=0 (common
        # connection) and all electrodes of the second array have channel_id=1
        self.electrode_arrays = [ElectrodeArray(center=self.center,
                                                radius=self.radius,
                                                length_x=self.length_x,
                                                length_y=self.length_y,
                                                channel_id=self.channel_id[self.channel_id == self.channel_id_unique[i]],
                                                ele_id=self.ele_id[self.channel_id == self.channel_id_unique[i]],
                                                current=self.current[self.channel_id == self.channel_id_unique[i]]) for i in range(2)]

        # compile node arrays
        self.compile_node_arrays()

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        if not hasattr(value, "__len__"):
            value = np.array([value])
        self._radius = value

    @property
    def length_x(self):
        return self._length_x

    @length_x.setter
    def length_x(self, value):
        if not hasattr(value, "__len__"):
            value = np.array([value])
        self._length_x = value

    @property
    def length_y(self):
        return self._length_y

    @length_y.setter
    def length_y(self, value):
        if not hasattr(value, "__len__"):
            value = np.array([value])
        self._length_y = value

    def set_geometrical_parameters_optimization(self, params):
        """
        Sets the geometrical parameters of electrode and update geometry.
        Only free parameters during optimization have to be provided in given order:
        radius, length_x, length_y

        Parameters
        ----------
        params : np.ndarray of float [n_free_geometry]
            Parameters in given order (radius, length_x, length_y) but only free ones
        """
        assert np.sum(self.free_geometry) == len(params), "Number of free parameters in electrode does not match " \
                                                          "number of provided parameters!"

        # the order of the parameters is fixed
        i_para = 0
        if self.radius_free:
            self.radius = params[i_para]
            i_para += 1

        if self.length_x_free:
            self.length_x = params[i_para]
            i_para += 1

        if self.length_y_free:
            self.length_y = params[i_para]
            i_para += 1

        # update geometry
        self.update_geometry()

    def update_geometry(self, center=None, radius=None, length_x=None, length_y=None):
        """
        Modifies the geometry of the electrode setup (e.g. for optimization). old parameters are kept if not changed.

        Parameters
        ----------
        center : np.ndarray of float [n_ele x 3]
            Center positions of electrodes in normalized x/y plane (z=0)
        radius : np.ndarray of float [n_ele_circ] or None
            Radii of circular electrodes (in mm)
        length_x : np.ndarray of float [n_ele_rect] or None
            Electrode extensions in x-direction (rectangular electrodes)
        length_y : np.ndarray of float [n_ele_rect] or None
            Electrode extensions in y-direction (rectangular electrodes)
        """
        if center is not None:
            self.center = center

        if radius is not None:
            if not hasattr(radius, "__len__"):
                radius = np.array([radius])
            self.radius = radius

        if length_x is not None:
            if not hasattr(length_x, "__len__"):
                length_x = np.array([length_x])
            self.length_x = length_x

        if length_y is not None:
            if not hasattr(length_y, "__len__"):
                length_y = np.array([length_y])
            self.length_y = length_y

        self.electrode_arrays = [ElectrodeArray(center=self.center,
                                                radius=self.radius,
                                                length_x=self.length_x,
                                                length_y=self.length_y,
                                                channel_id=self.channel_id[self.channel_id == self.channel_id_unique[i]],
                                                ele_id=self.ele_id[self.channel_id == self.channel_id_unique[i]],
                                                current=self.current[self.channel_id == self.channel_id_unique[i]]) for i in range(2)]


class CircularArray(ElectrodeMaster):
    """
    Generates a circular electrode array with one center electrode and n_outer equally spaced electrodes.
    Generates one ElectrodeArray instance because it can only be moved together.

    Parameters
    ----------
    radius_inner : float
        Radius of inner electrodes or (min, max) values for optimization.
    distance : float
        Distance between inner and outer electrodes (from center) (in mm)
    n_outer : int
        Number of outer electrodes or (min, max) values for optimization.
    radius_outer : float
        Radius of outer electrodes, default: same radius as inner electrodes or (min, max) values for optimization.
    radius_inner_bounds : list or np.ndarray of float [2], optional, default: None
        (min, max) values of radius of inner electrodes for optimization.
    distance_bounds : list or np.ndarray of float [2], optional, default: None
        (min, max) values of distance between inner and outer electrodes (from center) (in mm) for optimization.
    n_outer_bounds : list or np.ndarray of int [2], optional, default: None
        (min, max) values of number of outer electrodes for optimization.
    radius_outer_bounds : list or np.ndarray of float [2], optional, default: None
        (min, max) values radius of outer electrodes for optimization.
    current : np.ndarray of float [n_ele]
        Current through electrodes (First entry is central electrode).
        Has to sum up to zero net current.
    current_estimator_method : str, optional, default: "linear"
        Method to estimate the electrode currents:
        - "linear": linear regression
        - "gpc": generalized polynomial chaos
    dirichlet_correction : bool, optional, default: True
        If electrodes are connected to the same channel, they have to have the same voltage. This is ensured by
        setting this flag.
    dirichlet_correction_detailed : bool, optional, default: False
        Apply detailed Dirichlet correction such that every node current is optimized separately to match the equal
        voltage constraint of an electrode (recommended for large electrodes as in regular TES applications)
    current_outlier_correction : bool, optional, default: False
        Apply current outlier correction after node-wise dirichlet approximation (dirichlet_correction_detailed=True).

    Attributes
    ----------
    n_ele : int
        Number of electrodes
    center : np.ndarray of float [n_ele, 3]
        Center of electrodes (in mm)
    radius : list of float
        Radii of electrodes (in mm)
    electrode_arrays : list of ElectrodeArray instances [1]
        One ElectrodeArray instance containing the Electrode instances
    """
    def __init__(self,
                 radius_inner=None, distance=None, n_outer=None, radius_outer=None,
                 radius_inner_bounds=None, distance_bounds=None, n_outer_bounds=None, radius_outer_bounds=None,
                 current=None, current_estimator_method=None,
                 dirichlet_correction=True, dirichlet_correction_detailed=False, current_outlier_correction=False):

        if type(current) is int or type(current) is float:
            current = [current]

        if radius_outer is None:
            radius_outer = radius_inner

        # check free and fixed parameters
        if radius_inner_bounds is not None:
            assert len(radius_inner_bounds) == 2, ("Dimension mismatch! "
                                                   "Please provide [min, max] values for radius_inner_bounds.")
            self.radius_inner_free = True
            self.radius_inner_bounds = radius_inner_bounds
            self.radius_inner = np.mean(radius_inner_bounds)
        else:
            self.radius_inner_free = False
            self.radius_inner_bounds = [radius_inner, radius_inner]
            self.radius_inner = radius_inner

        if distance_bounds is not None:
            assert len(distance_bounds) == 2, ("Dimension mismatch! "
                                               "Please provide [min, max] values for distance_bounds.")
            self.distance_free = True
            self.distance_bounds = distance_bounds
            self.distance = np.mean(distance_bounds)
        else:
            self.distance_free = False
            self.distance_bounds = [distance, distance]
            self.distance = distance

        if n_outer_bounds is not None:
            assert len(n_outer_bounds) == 2, ("Dimension mismatch! "
                                              "Please provide [min, max] values for n_outer_bounds.")
            self.n_outer_free = True
            self.n_outer_bounds = n_outer_bounds
            self.n_outer = int(np.round(np.mean(n_outer_bounds)))
        else:
            self.n_outer_free = False
            self.n_outer_bounds = [n_outer, n_outer]
            self.n_outer = n_outer

        if radius_outer_bounds is not None:
            assert len(radius_outer_bounds) == 2, ("Dimension mismatch! "
                                                   "Please provide [min, max] values for radius_outer_bounds.")
            self.radius_outer_free = True
            self.radius_outer_bounds = radius_outer_bounds
            self.radius_outer = np.mean(radius_outer_bounds)
        else:
            self.radius_outer_free = False
            self.radius_outer_bounds = [radius_outer, radius_outer]
            self.radius_outer = radius_outer

        # order of free parameters
        self.free_geometry = np.array([self.radius_inner_free,
                                       self.distance_free,
                                       self.n_outer_free,
                                       self.radius_outer_free])
        self.any_free_geometry = self.free_geometry.any()
        self.geo_para_bounds = np.vstack((self.radius_inner_bounds,
                                          self.distance_bounds,
                                          self.n_outer_bounds,
                                          self.radius_outer_bounds))
        self.geo_para_mean = np.hstack((self.radius_inner,
                                        self.distance,
                                        self.n_outer,
                                        self.radius_outer))
        self.n_ele = self.n_outer + 1
        self.channel_id = np.array([0] + [1] * self.n_outer)
        self.channel_id_unique = np.unique(self.channel_id)
        self.n_channel = len(self.channel_id_unique)

        # global electrode arrays (set by compile_electrode_arrays)
        self.ele_id = np.arange(self.n_ele)
        # self.ele_channel_id = self.channel_id

        if current is None:
            self.current = np.hstack((1, -1/(self.n_ele-1) * np.ones(self.n_ele-1)))
        elif self.free_geometry[2]:
            self.current = np.hstack((current[0], -current[0]/(self.n_ele-1) * np.ones(self.n_ele-1)))
        else:
            self.current = current

        if len(self.current) != self.n_ele:
            self.current = np.hstack((self.current[0], -self.current[0] / (self.n_ele - 1) * np.ones(self.n_ele - 1)))

        if np.abs(np.sum(self.current)) > 1e-12:
            raise AssertionError("Please check electrode currents. They do not sum up to 0. (atol = 1e-12)")

        # number of electrodes per channel [n_channel]
        self.n_ele_per_channel = np.array([np.sum(self.channel_id == i) for i in np.unique(self.channel_id)])

        # total current entering domain (read from center electrode)
        self.current_total = self.current[0]

        # determine mean current of electrodes for each channel [n_channel]
        self.current_mean = self.current_total / self.n_ele_per_channel
        self.current_mean[1] *= -1

        # total current of each channel (here we only have 2)
        self.current_channel = np.array([self.current_total, -self.current_total])

        self.radius = np.append(np.array([self.radius_inner]), self.radius_outer*np.ones(self.n_outer))
        self.center = np.array([[0., 0.]])
        self.length_x = np.zeros(self.n_ele)
        self.length_y = np.zeros(self.n_ele)

        for i in range(self.n_outer):
            self.center = np.vstack((self.center, np.array([[0., 0.]])))
            self.center[-1, 0] = np.cos(i*(2*np.pi)/self.n_outer+np.pi/2) * self.distance
            self.center[-1, 1] = np.sin(i*(2*np.pi)/self.n_outer+np.pi/2) * self.distance

            if np.isclose(self.center[-1, 0], 0):
                self.center[-1, 0] = 0.

            if np.isclose(self.center[-1, 1], 0):
                self.center[-1, 1] = 0.

        # initialize current estimator for fake Dirichlet BC
        self.current_estimator_method = current_estimator_method
        if current_estimator_method is None or current_estimator_method == "" or (self.n_ele_per_channel == 1).all():
            current_estimator = None
        else:
            current_estimator = CurrentEstimator(method=current_estimator_method,
                                                 channel_id=self.channel_id,
                                                 ele_id=self.ele_id,
                                                 current_sign=np.sign(self.current),
                                                 current_total=self.current_total)

        super(CircularArray, self).__init__(dirichlet_correction=dirichlet_correction,
                                            dirichlet_correction_detailed=dirichlet_correction_detailed,
                                            current_outlier_correction=current_outlier_correction,
                                            current_estimator=current_estimator)

        # initialize freely movable electrode arrays
        self.electrode_arrays = [ElectrodeArray(channel_id=self.channel_id,
                                                ele_id=self.ele_id,
                                                center=self.center,
                                                radius=self.radius,
                                                length_x=self.length_x,
                                                length_y=self.length_y,
                                                current=self.current)]

        # compile node arrays
        self.compile_node_arrays()

    def set_geometrical_parameters_optimization(self, params):
        """
        Sets the geometrical parameters of electrode and update geometry.
        Only free parameters during optimization have to be provided in given order:
        radius_inner, distance, n_outer, radius_outer

        Parameters
        ----------
        params : np.ndarray of float [n_free_geometry]
            Parameters in given order (radius_inner, distance, n_outer, radius_outer) but only free ones
        """
        assert np.sum(self.free_geometry) == len(params), "Number of free parameters in electrode does not match " \
                                                          "number of provided parameters!"

        # the order of the parameters is fixed
        i_para = 0
        if self.radius_inner_free:
            self.radius_inner = params[i_para]
            i_para += 1

        if self.distance_free:
            self.distance = params[i_para]
            i_para += 1

        if self.n_outer_free:
            self.n_outer = int(np.round(params[i_para]))
            i_para += 1

        if self.radius_outer_free:
            self.radius_outer = params[i_para]
            i_para += 1

        # update geometry
        self.update_geometry()

    def update_geometry(self, radius_inner=None, distance=None, n_outer=None,  radius_outer=None):
        """
        Modifies the geometry of the electrode setup (e.g. for optimization). old parameters are kept if not changed.

        Parameters
        ----------
        radius_inner : float
            Radius of inner electrodes
        distance : float
            Distance between inner and outer electrodes (from center) (in mm)
        n_outer : int
            Number of outer electrodes
        radius_outer : float
            Radius of outer electrodes, default: same radius as inner electrodes
        """
        if radius_inner is not None:
            self.radius_inner = radius_inner

        if distance is not None:
            self.distance = distance

        if n_outer is not None:
            self.n_outer = int(np.round(n_outer))

        if radius_outer is not None:
            self.radius_outer = radius_outer

        self.n_ele = self.n_outer + 1
        self.ele_id = np.arange(self.n_ele)
        self.current = np.hstack((self.current_total, -self.current_total / (self.n_ele - 1) * np.ones(self.n_ele - 1)))

        if self.radius_outer is None:
            self.radius_outer = self.radius_inner
        else:
            self.radius_outer = self.radius_outer

        self.radius = np.append(np.array([self.radius_inner]), self.radius_outer*np.ones(self.n_outer))
        self.center = np.array([[0., 0.]])
        self.length_x = np.zeros(self.n_ele)
        self.length_y = np.zeros(self.n_ele)

        for i in range(self.n_outer):
            self.center = np.vstack((self.center, np.array([[0., 0.]])))
            self.center[-1, 0] = np.cos(i*(2*np.pi)/self.n_outer+np.pi/2) * self.distance
            self.center[-1, 1] = np.sin(i*(2*np.pi)/self.n_outer+np.pi/2) * self.distance

        self.electrode_arrays = [ElectrodeArray(channel_id=np.array([0] + [1] * self.n_outer),
                                                ele_id=self.ele_id,
                                                center=self.center,
                                                radius=self.radius,
                                                length_x=self.length_x,
                                                length_y=self.length_y,
                                                current=self.current)]


class ElectrodeArrayPairOpt(ElectrodeArrayPair):
    """
    Symmetric pair of electrode arrays with n_x * n_y electrodes each.
    Contains 2 ElectrodeArray instances with n_x * n_y Electrode instances each.
    Contains only circular electrodes.

    Parameters
    ----------
    n_ele_x : int
        Number of electrodes in x-direction
    n_ele_y : int
        Number of electrodes in y-direction
    separation_distance : float
        Distance between individual electrodes in the array from center to center in mm.
    radius : float
        Radius of electrodes in mm. All the electrodes in the array will have the radius.
    current_estimator_method : str, optional, default: "linear"
        Method to estimate the electrode currents:
        - "linear": linear regression
        - "gpc": generalized polynomial chaos
    dirichlet_correction : bool, optional, default: True
        If electrodes are connected to the same channel, they have to have the same voltage. This is ensured by
        setting this flag.
    dirichlet_correction_detailed : bool, optional, default: False
        Apply detailed Dirichlet correction such that every node current is optimized separately to match the equal
        voltage constraint of an electrode (recommended for large electrodes as in regular TES applications)
    current_outlier_correction : bool, optional, default: False
        Apply current outlier correction after node-wise dirichlet approximation (dirichlet_correction_detailed=True).

    Attributes
    ----------

    """

    def __init__(self, n_ele_x, n_ele_y, separation_distance, radius, current_estimator_method=None,
                 dirichlet_correction=True, dirichlet_correction_detailed=False, current_outlier_correction=False):
        """
        Initializes ElectrodeArrayPairOpt
        """
        # check free and fixed parameters
        if isinstance(n_ele_x, np.ndarray) or isinstance(n_ele_x, list):
            self.OPT_n_ele_x_free = True
            self.OPT_n_ele_x_bounds = n_ele_x
            self.OPT_n_ele_x = int(np.round(np.mean(n_ele_x)))

            # disabling current estimator because number of electrodes changes during optimization
            current_estimator_method = None
        else:
            self.OPT_n_ele_x_free = False
            self.OPT_n_ele_x_bounds = [n_ele_x, n_ele_x]
            self.OPT_n_ele_x = int(np.round(n_ele_x))

        if isinstance(n_ele_y, np.ndarray) or isinstance(n_ele_y, list):
            self.OPT_n_ele_y_free = True
            self.OPT_n_ele_y_bounds = n_ele_y
            self.OPT_n_ele_y = int(np.round(np.mean(n_ele_y)))

            # disabling current estimator because number of electrodes changes during optimization
            current_estimator_method = None
        else:
            self.OPT_n_ele_y_free = False
            self.OPT_n_ele_y_bounds = [n_ele_y, n_ele_y]
            self.OPT_n_ele_y = int(np.round(n_ele_y))

        if isinstance(separation_distance, np.ndarray) or isinstance(separation_distance, list):
            self.OPT_separation_distance_free = True
            self.OPT_separation_distance_bounds = separation_distance
            self.OPT_separation_distance = np.mean(separation_distance)
        else:
            self.OPT_separation_distance_free = False
            self.OPT_separation_distance_bounds = [separation_distance, separation_distance]
            self.OPT_separation_distance = separation_distance

        if isinstance(radius, np.ndarray) or isinstance(radius, list):
            self.OPT_radius_free = True
            self.OPT_radius_bounds = radius
            self.OPT_radius = np.mean(radius)
        else:
            self.OPT_radius_free = False
            self.OPT_radius_bounds = [radius, radius]
            self.OPT_radius = radius

        self.current_estimator_method = current_estimator_method
        self.dirichlet_correction = dirichlet_correction
        self.dirichlet_correction_detailed = dirichlet_correction_detailed
        self.current_outlier_correction = current_outlier_correction

        n_ele = self.OPT_n_ele_x * self.OPT_n_ele_y
        extent_x = 2 * self.OPT_radius + (self.OPT_n_ele_x - 1) * self.OPT_separation_distance
        extent_y = 2 * self.OPT_radius + (self.OPT_n_ele_y - 1) * self.OPT_separation_distance
        ele_pos_x = np.linspace(-extent_x / 2 + self.OPT_radius, extent_x / 2 - self.OPT_radius, self.OPT_n_ele_x)
        ele_pos_y = np.linspace(-extent_y / 2 + self.OPT_radius, extent_y / 2 - self.OPT_radius, self.OPT_n_ele_y)
        center = np.hstack((np.tile(ele_pos_x, (self.OPT_n_ele_y, 1)).flatten()[:, np.newaxis],
                            np.tile(ele_pos_y[:, np.newaxis], (1, self.OPT_n_ele_x)).flatten()[:, np.newaxis],
                            np.zeros((n_ele, 1))))

        # Initialize parent class (ElectrodeArrayPair)
        super(ElectrodeArrayPairOpt, self).__init__(center=center,
                                                    radius=self.OPT_radius * np.ones(n_ele),
                                                    length_x=np.zeros(n_ele),
                                                    length_y=np.zeros(n_ele),
                                                    current_estimator_method=current_estimator_method,
                                                    dirichlet_correction=dirichlet_correction,
                                                    dirichlet_correction_detailed=dirichlet_correction_detailed)

        # order of free parameters
        self.free_geometry = np.array([self.OPT_n_ele_x_free,
                                       self.OPT_n_ele_y_free,
                                       self.OPT_separation_distance_free,
                                       self.OPT_radius_free])
        self.any_free_geometry = self.free_geometry.any()
        self.geo_para_bounds = np.vstack((self.OPT_n_ele_x_bounds,
                                          self.OPT_n_ele_y_bounds,
                                          self.OPT_separation_distance_bounds,
                                          self.OPT_radius_bounds))
        self.geo_para_mean = np.hstack((self.OPT_n_ele_x,
                                        self.OPT_n_ele_y,
                                        self.OPT_separation_distance,
                                        self.OPT_radius))

    def update_geometry(self, n_ele_x=None, n_ele_y=None, separation_distance=None, radius=None):
        """
        Modifies the geometry of the electrode setup (e.g. for optimization). old parameters are kept if not changed.

        Parameters
        ----------
        n_ele_x : int, optional, default: None
            Number of electrodes in x-direction
        n_ele_y : int, optional, default: None
            Number of electrodes in y-direction
        separation_distance : float, optional, default: None
            Distance between individual electrodes in the array from center to center in mm.
        radius : float, optional, default: None
            Radius of electrodes in mm. All the electrodes in the array will have the radius.
        """

        if n_ele_x is not None:
            self.OPT_n_ele_x = int(np.round(n_ele_x))

        if n_ele_y is not None:
            self.OPT_n_ele_y = int(np.round(n_ele_y))

        if separation_distance is not None:
            self.OPT_separation_distance = separation_distance

        if radius is not None:
            self.OPT_radius = radius

        n_ele = self.OPT_n_ele_x * self.OPT_n_ele_y
        extent_x = 2 * self.OPT_radius + (self.OPT_n_ele_x - 1) * self.OPT_separation_distance
        extent_y = 2 * self.OPT_radius + (self.OPT_n_ele_y - 1) * self.OPT_separation_distance
        ele_pos_x = np.linspace(-extent_x / 2 + self.OPT_radius, extent_x / 2 - self.OPT_radius, self.OPT_n_ele_x)
        ele_pos_y = np.linspace(-extent_y / 2 + self.OPT_radius, extent_y / 2 - self.OPT_radius, self.OPT_n_ele_y)
        center = np.hstack((np.tile(ele_pos_x, (self.OPT_n_ele_y, 1)).flatten()[:, np.newaxis],
                            np.tile(ele_pos_y[:, np.newaxis], (1, self.OPT_n_ele_x)).flatten()[:, np.newaxis],
                            np.zeros((n_ele, 1))))

        # Initialize parent class (ElectrodeArrayPair)
        super(ElectrodeArrayPairOpt, self).__init__(center=center,
                                                    radius=self.OPT_radius * np.ones(n_ele),
                                                    length_x=np.zeros(n_ele),
                                                    length_y=np.zeros(n_ele),
                                                    current_estimator_method=self.current_estimator_method,
                                                    dirichlet_correction=self.dirichlet_correction,
                                                    dirichlet_correction_detailed=self.dirichlet_correction_detailed)

    def set_geometrical_parameters_optimization(self, params):
        """
        Sets the geometrical parameters of electrode and updates geometry.
        Only free parameters during optimization have to be provided in given order:
        n_ele_x, n_ele_y, separation_distance, radius

        Parameters
        ----------
        params : np.ndarray of float [n_free_geometry]
            Parameters in given order (n_ele_x, n_ele_y, separation_distance, radius) but only free ones
        """
        assert np.sum(self.free_geometry) == len(params), "Number of free parameters in electrode does not match " \
                                                          "number of provided parameters!"

        # the order of the parameters is fixed
        i_para = 0
        if self.OPT_n_ele_x_free:
            self.OPT_n_ele_x = int(np.round(params[i_para]))
            i_para += 1

        if self.OPT_n_ele_y_free:
            self.OPT_n_ele_y = int(np.round(params[i_para]))
            i_para += 1

        if self.OPT_separation_distance_free:
            self.OPT_separation_distance = params[i_para]
            i_para += 1

        if self.OPT_radius_free:
            self.OPT_radius = params[i_para]
            i_para += 1

        # update geometry
        self.update_geometry()


def create_tdcs_session_from_array(ElectrodeArray, fnamehead, pathfem, thickness=None,
                                   plug_center=None, plug_dimensions=None, rubber_size=None,
                                   sigma_rubber=None, sigma_saline=None):
    """
    Create a sim_struct.SESSION including a TDCSLIST object with ELECTRODE instances for regular TDCS
    simulations including electrode meshing etc. for reference simulations.
    Returns the SESSION object "s", which can be run with run_simnibs(s)

    Parameters
    ----------
    ElectrodeArray : CircularArray or Electrode object

    fnamehead : str
        Path to headmodel .msh file
    pathfem : str
        Output folder of simulation data
    thickness : list of float [1-3]
        Can have up to 3 arguments. 1st argument is the lower part of the electrode, 2nd arrgument is the rubber;
        3rd is the upper part. Examples:
        - thickness = [1] simple electrode with 1 mm gel
        - thickness = [1, 2] 1 mm gel and 2 mm rubber (good conductor)
        - thickness = [3.5, 1, 3.5] 3.5 mm sponge, 1 mm rubber (good conductor), 3.5 mm sponge
        If a sponge electrode is defined, the size of the sponge is the size of the given electrode. The
        size of the rubber inside the two sponge layers has to be specified separately using the rubber_size parameter.
    plug_center : list of float [2]
        If provided, a plug will be added to the electrode where the current is injected. The location can be provided
        in relative coordinates in mm with respect to the center of the electrode.
    plug_dimensions : list of float [2]
        Size of the plug in mm. If not provided, the default value of [10, 10] will be used.
    rubber_size : list of float [2]
        Size of the rubber in mm in case of sponge electrodes with len(thickness) = 3

    Returns
    -------
    s : simnibs.sim_struct.SESSION() object
        SESSION object, which includes a TDCSLIST and the corresponding ELECTRODE instances.
        Can be run with simnibs.run_simnibs(s).
    """

    if thickness is None:
        thickness = [1, 2]

    # Initalize a session
    s = SESSION()
    # Name of head mesh
    s.fnamehead = fnamehead
    # Output folder
    s.pathfem = pathfem
    # disable visualization
    s.open_in_gmsh = False

    # Initialize a tDCS simulation
    tdcslist = s.add_tdcslist()

    if sigma_rubber is not None:
        tdcslist.cond[99].value = sigma_rubber

    if sigma_saline is not None:
        tdcslist.cond[499].value = sigma_saline

    # Set currents
    tdcslist.currents = ElectrodeArray.current_channel

    # Initialize the electrodes
    for i_array, _electrode_array in enumerate(ElectrodeArray.electrode_arrays):
        for i_ele, _electrode in enumerate(_electrode_array.electrodes):
            # add new electrode
            electrode = tdcslist.add_electrode()

            if _electrode.type == "spherical":
                # Circular shape
                electrode.shape = 'ellipse'

                # Electrode (rubber) and Sponge dimension
                if len(thickness) == 3:
                    electrode.dimensions_sponge = [2 * _electrode.radius, 2 * _electrode.radius]
                    electrode.dimensions = rubber_size[i_ele]
                else:
                    electrode.dimensions = [2*_electrode.radius, 2*_electrode.radius]

            elif _electrode.type == "rectangular":
                # Rectangular shape
                electrode.shape = 'rect'

                # Electrode (rubber) and Sponge dimension
                if len(thickness) == 3:
                    electrode.dimensions_sponge = [_electrode.length_x, _electrode.length_y]
                    electrode.dimensions = rubber_size[i_ele]
                else:
                    electrode.dimensions = [_electrode.length_x, _electrode.length_y]

            else:
                raise AssertionError("Electrodes have to be either 'spherical' or 'rectangular'")

            # add plug
            if plug_center is not None:
                plug = electrode.add_plug()
                plug.centre = plug_center
                plug.shape = electrode.shape
                plug.dimensions = [10, 10]

                if plug_dimensions is not None:
                    plug.dimensions = plug_dimensions

            # Connect electrode to first channel (-1e-3 mA, cathode)
            electrode.channelnr = _electrode.channel_id

            # electrode thickness
            electrode.thickness = thickness  # 1: gel; 2: gel + conductive rubber; 3: sponge + rubber + sponge

            # Electrode Position
            electrode.centre = _electrode.posmat[:3, 3]

            # Electrode direction
            electrode.pos_ydir = electrode.centre + 20*_electrode.posmat[:3, 1]

    return s

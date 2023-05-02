import copy
import numpy as np
from .current_estimator import CurrentEstimator


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
    current : float, optional, default: None
        Current assigned to electrode
    voltage : float, optional, default: None
        Voltage assigned to electrode

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
    current : float, optional, default: None
        Current assigned to electrode
    voltage : float, optional, default: None
        Voltage assigned to electrode
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
                                     [0, 0, 1, center[2]],
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
                self._ele_current = np.sum(value * self.node_area / self.node_area_total)
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
            if self.node_area is not None:
                self._ele_voltage = np.sum(value * self.node_area / self.node_area_total)
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
    center : np.ndarray of float [n_ele, 3]
        Center coordinates (x,y,z) of electrodes in normalized space (in mm)
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
    """

    def __init__(self, channel_id, ele_id, center, radius=None, length_x=None, length_y=None, current=None):
        self.channel_id = channel_id
        self.ele_id = ele_id
        self.array_center = np.array([0, 0, 0])
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

        # apply Dirichlet correction step in fast TES calculations if electrodes share the same channel
        if self.n_ele != len(np.unique(self.channel_id)):
            self.dirichlet_correction = True
        else:
            self.dirichlet_correction = False

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
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt

        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']

        if len(colors) < np.max(self.channel_id):
            colors = colors * int(np.ceil(colors/np.max(self.channel_id)))

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


class ElectrodeArrayPair():
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
    center : np.ndarray of float [n_ele/2 x 3]
        Center positions of electrodes in normalized x/y plane (z=0). Array will be copied for second array.
    radius : np.ndarray of float [n_ele/2] or list [2] or None
        Radii of single circular electrodes in array (in mm). Array will be copied for second array.
        If a list is passed, it has to be of length 2 and the parameter will be included in the optimization
        varying between the provided values [min, max].
        (i.e. optimization of geometrical parameters can only performed for standard TES with 2 electrodes)
    length_x : np.ndarray of float [n_ele/2] or list [2] or None
        Extensions in x-direction of single rectangular electrodes in array. Array will be copied for second array.
        If a list is passed, it has to be of length 2 and the parameter will be included in the optimization
        varying between the provided values [min, max].
        (i.e. optimization of geometrical parameters can only performed for standard TES with 2 electrodes)
    length_y : np.ndarray of float [n_ele/2] or list [2] or None
        Extensions in y-direction of single rectangular electrodes in array. Array will be copied for second array.
        If a list is passed, it has to be of length 2 and the parameter will be included in the optimization
        varying between the provided values [min, max].
        (i.e. optimization of geometrical parameters can only performed for standard TES with 2 electrodes)
    current : np.ndarray of float [n_ele] or None
        Current through electrodes. Has to sum up to zero net current.
    current_estimator_method : str, optional, default: "linear"
        Method to estimate the electrode currents:
        - "linear": linear regression
        - "gpc": generalized polynomial chaos
    dirichlet_correction_detailed : bool, optional, default: False
        Apply detailed Dirichlet correction such that every node current is optimized separately to match the equal
        voltage constraint of an electrode (recommended for large electrodes as in regular TES applications)

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
                 current_estimator_method="gpc", dirichlet_correction_detailed=False):

        # check free and fixed parameters
        if isinstance(radius, list):
            assert len(radius) == 2, "radius is provided as list, optimizing this parameter but length != 2 [min, max]"
            self.radius_free = True
            self.radius_bounds = radius
            self.radius = np.array([np.mean(radius)])
        else:
            self.radius_free = False
            self.radius_bounds = [0, 0]
            self.radius = radius

        if isinstance(length_x, list):
            assert len(length_x) == 2, "length_x is provided as list, optimizing this parameter but length != 2 [min, max]"
            self.length_x_free = True
            self.length_x_bounds = length_x
            self.length_x = np.array([np.mean(length_x)])
        else:
            self.length_x_free = False
            self.length_x_bounds = [0, 0]
            self.length_x = length_x

        if isinstance(length_y, list):
            assert len(length_y) == 2, "length_y is provided as list, optimizing this parameter but length != 2 [min, max]"
            self.length_y_free = True
            self.length_y_bounds = length_y
            self.length_y = np.array([np.mean(length_y)])
        else:
            self.length_y_free = False
            self.length_y_bounds = [0, 0]
            self.length_y = length_y

        self.center = center
        self.n_ele = len(radius) * 2                                                        # total number of electrodes
        self.channel_id = np.hstack([[i for _ in range(len(center))] for i in range(2)])    # array of all channel_ids
        self.channel_id_unique = np.unique(self.channel_id)
        self.n_channel = len(self.channel_id_unique)
        self.dirichlet_correction_detailed = dirichlet_correction_detailed
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

        if len(self.channel_id) != len(self.channel_id_unique):
            self.dirichlet_correction = True
        else:
            self.dirichlet_correction = False

        if current is None:
            self.current = np.hstack((np.array([1/(self.n_ele/2.) for _ in range(int(self.n_ele/2))]),
                                      np.array([-1/(self.n_ele/2.) for _ in range(int(self.n_ele/2))])))
        else:
            self.current = current

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
        if current_estimator_method is None or current_estimator_method == "" or (self.n_ele_per_channel == 1).all():
            self.current_estimator = None
        else:
            self.current_estimator = CurrentEstimator(method=current_estimator_method,
                                                      channel_id=self.channel_id,
                                                      ele_id=self.ele_id,
                                                      current_sign=np.sign(self.current),
                                                      current_total=self.current_total)

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


class CircularArray():
    """
    Generates a circular electrode array with one center electrode and n_outer equally spaced electrodes.
    Generates one ElectrodeArray instance because it can only be moved together.

    Parameters
    ----------
    radius_inner : float or np.ndarray of float [2]
        Radius of inner electrodes or (min, max) values for optimization.
    distance : float or np.ndarray of float [2]
        Distance between inner and outer electrodes (from center) (in mm)
    n_outer : int or np.ndarray of int [2]
        Number of outer electrodes or (min, max) values for optimization.
    radius_outer : float or np.ndarray of float [2]
        Radius of outer electrodes, default: same radius as inner electrodes or (min, max) values for optimization.
    current : np.ndarray of float [n_ele]
        Current through electrodes (First entry is central electrode).
        Has to sum up to zero net current.
    current_estimator_method : str, optional, default: "linear"
        Method to estimate the electrode currents:
        - "linear": linear regression
        - "gpc": generalized polynomial chaos
    dirichlet_correction_detailed : bool, optional, default: False
        Apply detailed Dirichlet correction such that every node current is optimized separately to match the equal
        voltage constraint of an electrode (recommended for large electrodes as in regular TES applications)

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
                 radius_inner, distance, n_outer=4, radius_outer=None,
                 current=None, current_estimator_method="gpc", dirichlet_correction_detailed=False):

        if radius_outer is None:
            radius_outer = radius_inner

        # check free and fixed parameters
        if isinstance(radius_inner, np.ndarray) or isinstance(radius_inner, list):
            self.radius_inner_free = True
            self.radius_inner_bounds = radius_inner
            self.radius_inner = np.mean(radius_inner)
        else:
            self.radius_inner_free = False
            self.radius_inner_bounds = [radius_inner, radius_inner]
            self.radius_inner = radius_inner

        if isinstance(distance, np.ndarray) or isinstance(distance, list):
            self.distance_free = True
            self.distance_bounds = distance
            self.distance = np.mean(distance)
        else:
            self.distance_free = False
            self.distance_bounds = [distance, distance]
            self.distance = distance

        if isinstance(n_outer, np.ndarray) or isinstance(n_outer, list):
            self.n_outer_free = True
            self.n_outer_bounds = n_outer
            self.n_outer = int(np.round(np.mean(n_outer)))
        else:
            self.n_outer_free = False
            self.n_outer_bounds = [n_outer, n_outer]
            self.n_outer = n_outer

        if isinstance(radius_outer, np.ndarray) or isinstance(radius_outer, list):
            self.radius_outer_free = True
            self.radius_outer_bounds = radius_outer
            self.radius_outer = np.mean(radius_outer)
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
        self.dirichlet_correction = True
        self.dirichlet_correction_detailed = dirichlet_correction_detailed

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

        # global electrode arrays (set by compile_electrode_arrays)
        self.ele_id = np.arange(self.n_ele)
        # self.ele_channel_id = self.channel_id

        if current is None:
            self.current = np.hstack((1, -1/(self.n_ele-1) * np.ones(self.n_ele-1)))
        else:
            self.current = current

        if np.abs(np.sum(self.current)) > 1e-12:
            raise AssertionError("Please check electrode currents. They do not sum up to 0. (atol = 1e-12)")

        # number of electrodes per channel [n_channel]
        self.n_ele_per_channel =  np.array([np.sum(self.channel_id == i) for i in np.unique(self.channel_id)])

        # total current entering domain (read from center electrode)
        self.current_total = self.current[0]

        # determine mean current of electrodes for each channel [n_channel]
        self.current_mean = self.current_total / self.n_ele_per_channel
        self.current_mean[1] *= -1

        # total current of each channel (here we only have 2)
        self.current_channel = np.array([self.current_total, -self.current_total])

        self.radius = np.append(np.array([self.radius_inner]), self.radius_outer*np.ones(self.n_outer))
        self.center = np.array([[0., 0., 0.]])
        self.length_x = np.zeros(self.n_ele)
        self.length_y = np.zeros(self.n_ele)

        for i in range(self.n_outer):
            self.center = np.vstack((self.center, np.array([[0., 0., 0.]])))
            self.center[-1, 0] = np.cos(i*(2*np.pi)/self.n_outer+np.pi/2) * self.distance
            self.center[-1, 1] = np.sin(i*(2*np.pi)/self.n_outer+np.pi/2) * self.distance

            if np.isclose(self.center[-1, 0], 0):
                self.center[-1, 0] = 0.

            if np.isclose(self.center[-1, 1], 0):
                self.center[-1, 1] = 0.

        # initialize current estimator for fake Dirichlet BC
        if current_estimator_method is None or current_estimator_method == "" or (self.n_ele_per_channel == 1).all():
            self.current_estimator = None
        else:
            self.current_estimator = CurrentEstimator(method=current_estimator_method,
                                                      channel_id=self.channel_id,
                                                      ele_id=self.ele_id,
                                                      current_sign=np.sign(self.current),
                                                      current_total=self.current_total)

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
        self.current = np.hstack((1, -1 / (self.n_ele - 1) * np.ones(self.n_ele - 1)))

        if self.radius_outer is None:
            self.radius_outer = self.radius_inner
        else:
            self.radius_outer = self.radius_outer

        self.radius = np.append(np.array([self.radius_inner]), self.radius_outer*np.ones(self.n_outer))
        self.center = np.array([[0., 0., 0.]])
        self.length_x = np.zeros(self.n_ele)
        self.length_y = np.zeros(self.n_ele)

        for i in range(self.n_outer):
            self.center = np.vstack((self.center, np.array([[0., 0., 0.]])))
            self.center[-1, 0] = np.cos(i*(2*np.pi)/self.n_outer+np.pi/2) * self.distance
            self.center[-1, 1] = np.sin(i*(2*np.pi)/self.n_outer+np.pi/2) * self.distance

        self.electrode_arrays = [ElectrodeArray(channel_id=np.array([0] + [1] * self.n_outer),
                                                ele_id=self.ele_id,
                                                center=self.center,
                                                radius=self.radius,
                                                length_x=self.length_x,
                                                length_y=self.length_y,
                                                current=self.current)]

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

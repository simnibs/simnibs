import copy
import numpy as np


class Electrode():
    """ Electrode class. Contains a single electrode. ElectrodeArrays are built from Electrode instances

    Parameters
    ----------
    channel_id : int
        Channel identifier, indicates connected electrodes (same identifier)
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
        Associated area of the points
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

    def __init__(self, channel_id, center, radius=None, length_x=None, length_y=None, current=None, voltage=None):
        if radius is None:
            radius = 0

        if length_x is None:
            length_x = 0

        if length_y is None:
            length_y = 0

        self.channel_id = channel_id
        self.center = center
        self.radius = radius
        self.length_x = length_x
        self.length_y = length_y
        self.current = current
        self.voltage = voltage
        self.transmat = None
        self.node_area = None
        self.node_idx = None
        self.node_coords = None
        self.n_nodes = None
        self.posmat_norm = np.array([[1, 0, 0, center[0]],
                                     [0, 1, 0, center[1]],
                                     [0, 0, 1, center[2]],
                                     [0, 0, 0, 1]])
        self.posmat = copy.deepcopy(self.posmat_norm)

        if self.voltage is not None and self.current is not None:
            raise AssertionError("Define either voltage or current to electrode.")

        if radius != 0 and (length_x != 0 or length_y != 0):
            raise AssertionError("Define either radius for circular electrode or "
                                 "length_x and length_y for rectangular electrode")

        if self.radius is not None and radius != 0:
            self.area = np.pi*radius**2
            self.type = "spherical"
        else:
            self.area = length_x * length_y
            self.type = "rectangular"

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
    """

    def __init__(self, channel_id, center, radius=None, length_x=None, length_y=None, current=None):
        self.channel_id = channel_id
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
                                             center=self.center[i_ele, :],
                                             radius=self.radius[i_ele],
                                             length_x=self.length_x[i_ele],
                                             length_y=self.length_y[i_ele],
                                             current=self.current[i_ele]))

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
            xycoords = 'data', textcoords = 'offset points')

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
    radius : np.ndarray of float [n_ele/2] or None
        Radii of single circular electrodes in array (in mm). Array will be copied for second array.
    length_x : np.ndarray of float [n_ele/2] or None
        Extensions in x-direction of single rectangular electrodes in array. Array will be copied for second array.
    length_y : np.ndarray of float [n_ele/2] or None
        Extensions in y-direction of single rectangular electrodes in array. Array will be copied for second array.
    current : np.ndarray of float [n_ele] or None
        Current through electrodes. Has to sum up to zero net current.

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

    def __init__(self, center, radius=None, length_x=None, length_y=None, current=None):
        self.radius = radius
        self.length_x = length_x
        self.length_y = length_y
        self.center = center
        self.n_ele = len(center) * 2                                                        # total number of electrodes
        self.channel_id = np.hstack([[i for _ in range(len(center))] for i in range(2)])    # array of all channel_ids
        self.n_channel = np.unique(self.channel_id)

        if current is None:
            self.current = np.hstack((np.array([1/(self.n_ele/2.) for _ in range(int(self.n_ele/2))]),
                                      np.array([-1/(self.n_ele/2.) for _ in range(int(self.n_ele/2))])))
        else:
            self.current = current

        if np.abs(np.sum(self.current)) > 1e-12:
            raise AssertionError("Please check electrode currents. They do not sum up to 0. (atol = 1e-12)")

        # number of electrodes per channel [n_channel]
        self.n_ele_per_channel = [np.sum(self.channel_id == i) for i in np.unique(self.channel_id)]

        # total current entering domain (read from first channel)
        self.current_total = np.sum(self.current[self.channel_id == self.channel_id[0]])

        # determine mean current of electrodes for each channel [n_channel]
        self.current_mean = self.current_total / self.n_ele_per_channel
        self.current_mean[1] *= -1

        # create two ElectrodeArray instance where all electrodes of the first array have channel_id=0 (common
        # connection) and all electrodes of the second array have channel_id=1
        self.electrode_arrays = [ElectrodeArray(center=self.center,
                                                radius=self.radius,
                                                length_x=self.length_x,
                                                length_y=self.length_y,
                                                channel_id=self.channel_id[i],
                                                current=self.current[i]) for i in range(2)]

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
            self.radius = radius

        if length_x is not None:
            self.length_x = length_x

        if length_y is not None:
            self.length_y = length_y

        self.electrode_arrays = [ElectrodeArray(center=self.center,
                                                radius=self.radius,
                                                length_x=self.length_x,
                                                length_y=self.length_y,
                                                channel_id=[i for _ in range(self.n_ele)]) for i in range(2)]


class CircularArray():
    """
    Generates a circular electrode array with one center electrode and n_outer equally spaced electrodes.
    Generates one ElectrodeArray instance because it can only be moved together.

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
    current : np.ndarray of float [n_ele]
        Current through electrodes (First entry is central electrode).
        Has to sum up to zero net current.

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
    def __init__(self, radius_inner, distance, n_outer=4, radius_outer=None, current=None):
        self.radius_inner = radius_inner
        self.distance = distance
        self.n_ele = n_outer + 1
        self.n_outer = n_outer
        self.channel_id = np.array([0] + [1] * n_outer)
        self.n_channel = len(np.unique(self.channel_id))
        self.dirichlet_correction = True

        if radius_outer is None:
            self.radius_outer = radius_inner
        else:
            self.radius_outer = radius_outer

        if current is None:
            self.current = np.hstack((1, -1/(self.n_ele-1) * np.ones(self.n_ele-1)))
        else:
            self.current = current

        if np.abs(np.sum(self.current)) > 1e-12:
            raise AssertionError("Please check electrode currents. They do not sum up to 0. (atol = 1e-12)")

        # number of electrodes per channel [n_channel]
        self.n_ele_per_channel = [np.sum(self.channel_id == i) for i in np.unique(self.channel_id)]

        # total current entering domain (read from center electrode)
        self.current_total = self.current[0]

        # determine mean current of electrodes for each channel [n_channel]
        self.current_mean = self.current_total / self.n_ele_per_channel
        self.current_mean[1] *= -1

        self.radius = np.append(np.array([self.radius_inner]), self.radius_outer*np.ones(n_outer))
        self.center = np.array([[0., 0., 0.]])
        self.length_x = np.zeros(self.n_ele)
        self.length_y = np.zeros(self.n_ele)

        for i in range(n_outer):
            self.center = np.vstack((self.center, np.array([[0., 0., 0.]])))
            self.center[-1, 0] = np.cos(i*(2*np.pi)/n_outer+np.pi/2) * distance
            self.center[-1, 1] = np.sin(i*(2*np.pi)/n_outer+np.pi/2) * distance

            if np.isclose(self.center[-1, 0], 0):
                self.center[-1, 0] = 0.

            if np.isclose(self.center[-1, 1], 0):
                self.center[-1, 1] = 0.

        self.electrode_arrays = [ElectrodeArray(channel_id=self.channel_id,
                                                center=self.center,
                                                radius=self.radius,
                                                length_x=self.length_x,
                                                length_y=self.length_y,
                                                current=self.current)]

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
            self.n_outer = n_outer

        if radius_outer is not None:
            self.radius_outer = radius_outer

        self.n_ele = self.n_outer + 1

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
                                                center=self.center,
                                                radius=self.radius,
                                                length_x=self.length_x,
                                                length_y=self.length_y)]

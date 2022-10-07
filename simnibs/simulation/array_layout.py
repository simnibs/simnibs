import numpy as np


class ElectrodeArray():
    """ Array of physically connected electrodes. Can only be moved together

    Parameters
    ----------

    Attributes
    ----------

    """

    def __init__(self, electrode_center, electrode_radii, channel_id):
        self.channel_id = electrode_center
        self.electrode_center = electrode_radii
        self.electrode_radii = channel_id

    def transform(self, transmat):
        """
        Transforms electrode array according to given 4x4 transformation matrix

        Parameters
        ----------
        transmat : nparray of float [4x4]
            Transformation matrix to place the array
        """
        pass

class ArrayLayout():
    """ Array Layout class for TES

    Parameters
    ----------

    Attributes
    ----------

    """

    def __init__(self, electrode_arrays):
        self.electrode_arrays = electrode_arrays


class TTF9Pair(ArrayLayout):
    """
    Pair of 3x3 electrodes used for TTF

    Parameters
    ----------

    Attributes
    ----------

    """

    def __init__(self, n_x=3, n_y=3, X=100, Y=50, r=10):

        electrode_center = np.array(np.meshgrid(np.linspace(-X/2, X/2, n_x),
                                                np.linspace(-Y/2, Y/2, n_y))).T.reshape(-1, 2)
        electrode_radii = [r for _ in range(n_x*n_y)]
        electrode_arrays = [0 for _ in range(2)]

        for i in range(2):
            electrode_arrays[i] = ElectrodeArray(electrode_center=electrode_center,
                                                 electrode_radii=electrode_radii,
                                                 channel_id=[i for _ in range(n_x*n_y)])


        super(TTF9Pair, self).__init__(electrode_arrays=electrode_arrays)


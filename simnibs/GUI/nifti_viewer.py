''' Nifti viewer Based on the OrthoSlicer3D from nibabel
Modified to handle overlays
'''

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

from nibabel.affines import voxel_sizes
from nibabel.orientations import aff2axcodes, axcodes2ornt
import simnibs
'''
The NiftiViewer class is based on the OthoSlicer3D class from nibabel, and has the following license

The MIT License

Copyright (c) 2009-2019 Matthew Brett <matthew.brett@gmail.com>
Copyright (c) 2010-2013 Stephan Gerhard <git@unidesign.ch>
Copyright (c) 2006-2014 Michael Hanke <michael.hanke@gmail.com>
Copyright (c) 2011 Christian Haselgrove <christian.haselgrove@umassmed.edu>
Copyright (c) 2010-2011 Jarrod Millman <jarrod.millman@gmail.com>
Copyright (c) 2011-2019 Yaroslav Halchenko <debian@onerussian.com>
Copyright (c) 2015-2019 Chris Markiewicz <effigies@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


Modifications by Guilherme Saturnino, 2019
'''


class NiftiViewer(object):
    def __init__(self, volume, affine=None, title=None, cmap='gray', clim=None, alpha=1.):
        """
        Parameters
        ----------
        volume : array-like
            The data that will be displayed by the slicer. Should have 3
            dimensions.
        affine : array-like or None, optional
            Affine transform for the data. This is used to determine
            how the data should be sliced for plotting into the sagittal,
            coronal, and axial view axes. If None, identity is assumed.
            The aspect ratio of the data are inferred from the affine
            transform.
        title : str or None, optional
            The title to display. Can be None (default) to display no
            title.
        cmap: matplotlib colormap, optional
            Colormap to use for ploting. Default: 'gray'
        clim: [min, max] or None
            Limits to use for plotting. Default: 1 and 99th percentiles
        alpha: float
            Transparency value
        """
        # Use these late imports of matplotlib so that we have some hope that
        # the test functions are the first to set the matplotlib backend. The
        # tests set the backend to something that doesn't require a display.
        self._title = title
        self._closed = False
        self._cross = True

        volume = np.asanyarray(volume)
        if volume.ndim < 3:
            raise ValueError('volume must have at least 3 dimensions')
        if np.iscomplexobj(volume):
            raise TypeError("Complex data not supported")
        affine = np.array(affine, float) if affine is not None else np.eye(4)
        if affine.shape != (4, 4):
            raise ValueError('affine must be a 4x4 matrix')
        # determine our orientation
        self._affine = affine
        codes = axcodes2ornt(aff2axcodes(self._affine))
        self._order = np.argsort([c[0] for c in codes])
        self._flips = np.array([c[1] < 0 for c in codes])[self._order]
        self._flips = list(self._flips) + [False]  # add volume dim
        self._scalers = voxel_sizes(self._affine)
        self._inv_affine = np.linalg.inv(affine)
        # current volume info
        self._volume_dims = volume.shape[3:]
        if len(self._volume_dims) > 0:
            raise NotImplementedError('Cannot handle 4-D Datasets')
        self._volumes = []

        # ^ +---------+   ^ +---------+
        # | |         |   | |         |
        #   |   Sag   |     |   Cor   |
        # S |    0    |   S |    1    |
        #   |         |     |         |
        #   |         |     |         |
        #   +---------+     +---------+
        #        A  -->
        # ^ +---------+
        # | |         |
        #   |  Axial  |
        # A |    2    |
        #   |         |
        #   |         |
        #   +---------+
        #   <--  R
        fig, axes = plt.subplots(2, 2)
        fig.set_size_inches((8, 8), forward=True)
        self._axes = [axes[0, 0], axes[0, 1], axes[1, 0]]
        plt.tight_layout(pad=0.1)
        fig.delaxes(axes[1, 1])
        if self._title is not None:
            fig.canvas.set_window_title(str(title))

        # Start midway through each axis, idx is current slice number
        self._ims, self._data_idx = list(), list()

        # set up axis crosshairs
        self._crosshairs = [None] * 3
        r = [self._scalers[self._order[2]] / self._scalers[self._order[1]],
             self._scalers[self._order[2]] / self._scalers[self._order[0]],
             self._scalers[self._order[1]] / self._scalers[self._order[0]]]
        self._sizes = [volume.shape[order] for order in self._order]
        for ii, xax, yax, ratio, label in zip([0, 1, 2], [1, 0, 0], [2, 2, 1],
                                              r, ('SAIP', 'SRIL', 'ARPL')):
            ax = self._axes[ii]
            vert = ax.plot([0] * 2, [-0.5, self._sizes[yax] - 0.5],
                           color=(0, 1, 0), linestyle='-')[0]
            horiz = ax.plot([-0.5, self._sizes[xax] - 0.5], [0] * 2,
                            color=(0, 1, 0), linestyle='-')[0]
            self._crosshairs[ii] = dict(vert=vert, horiz=horiz)
            # add text labels (top, right, bottom, left)
            lims = [0, self._sizes[xax], 0, self._sizes[yax]]
            bump = 0.01
            poss = [[lims[1] / 2., lims[3]],
                    [(1 + bump) * lims[1], lims[3] / 2.],
                    [lims[1] / 2., 0],
                    [lims[0] - bump * lims[1], lims[3] / 2.]]
            anchors = [['center', 'bottom'], ['left', 'center'],
                       ['center', 'top'], ['right', 'center']]
            for pos, anchor, lab in zip(poss, anchors, label):
                ax.text(pos[0], pos[1], lab,
                        horizontalalignment=anchor[0],
                        verticalalignment=anchor[1])
            ax.axis(lims)
            ax.set_aspect(ratio)
            ax.patch.set_visible(False)
            ax.set_frame_on(False)
            ax.axes.get_yaxis().set_visible(False)
            ax.axes.get_xaxis().set_visible(False)
            self._data_idx.append(0)
        self._data_idx.append(-1)  # volume

        self._figs = set([a.figure for a in self._axes])
        for fig in self._figs:
            fig.canvas.mpl_connect('scroll_event', self._on_scroll)
            fig.canvas.mpl_connect('motion_notify_event', self._on_mouse)
            fig.canvas.mpl_connect('button_press_event', self._on_mouse)

        # actually set data meaningfully
        self.add_overlay(volume, cmap=cmap, clim=clim, alpha=alpha, draw=False)
        self._position = np.zeros(4)
        self._position[3] = 1.  # convenience for affine multiplication
        self._changing = False  # keep track of status to avoid loops
        plt.draw()
        for fig in self._figs:
            fig.canvas.draw_idle()
            fig.canvas.draw()
        plt.pause(1e-3) # give a little bit of time for the renderer (needed on MacOS)
        self._set_position(0., 0., 0.)
        self._draw()

    def add_overlay(self, vol, cmap='gray', clim=None, alpha=1., draw=True):
        '''
        volume : array-like
            The data that will be overlaid. Must have the same dimensions as the original
            plot
        cmap: matplotlib colormap, optional
            Colormap to use for ploting. Default: 'gray'
        clim: [min, max] or None
            Limits to use for plotting. Default: 1 and 99th percentiles
        alpha: float
            transparency value
        '''
        if len(self._volumes) > 0 and vol.shape != self._volumes[0].shape:
            raise ValueError('Cannot add overlay, different shape')
        # add volume
        self._volumes.append(vol)
        if clim is None:
            clim = np.percentile(vol, (1., 99.))
        # create new images
        ims = []
        for ii, xax, yax in zip([0, 1, 2], [1, 0, 0], [2, 2, 1]):
            d = np.zeros((self._sizes[yax], self._sizes[xax]))
            im = self._axes[ii].imshow(
                d, aspect=1, cmap=cmap, clim=clim, alpha=alpha,
                interpolation='nearest', origin='lower'
            )
            ims.append(im)
        self._ims.append(ims)
        if draw:
            self.set_position()

    def __repr__(self):
        title = '' if self._title is None else ('%s ' % self._title)
        vol = '' if self.n_volumes <= 1 else (', %s' % self.n_volumes)
        r = ('<%s: %s(%s, %s, %s%s)>'
             % (self.__class__.__name__, title, self._sizes[0], self._sizes[1],
                self._sizes[2], vol))
        return r

    # User-level functions ###################################################
    def show(self):
        """Show the slicer in blocking mode; convenience for ``plt.show()``
        """
        plt.show()

    def close(self):
        """Close the viewer figures
        """
        for f in self._figs:
            plt.close(f)

    def draw(self):
        """Redraw the current image"""
        for fig in self._figs:
            fig.canvas.draw()

    @property
    def position(self):
        """The current coordinates"""
        return self._position[:3].copy()

    @property
    def figs(self):
        """A tuple of the figure(s) containing the axes"""
        return tuple(self._figs)

    def set_position(self, x=None, y=None, z=None):
        """Set current displayed slice indices
        Parameters
        ----------
        x : float | None
            X coordinate to use. If None, do not change.
        y : float | None
            Y coordinate to use. If None, do not change.
        z : float | None
            Z coordinate to use. If None, do not change.
        """
        self._set_position(x, y, z)
        self._draw()

    def _set_position(self, x, y, z):
        """Set the plot data using a physical position"""
        # deal with volume first
        if self._changing:
            return
        self._changing = True
        x = self._position[0] if x is None else float(x)
        y = self._position[1] if y is None else float(y)
        z = self._position[2] if z is None else float(z)

        # deal with slicing appropriately
        self._position[:3] = [x, y, z]
        idxs = np.dot(self._inv_affine, self._position)[:3]
        for ii, (size, idx) in enumerate(zip(self._sizes, idxs)):
            self._data_idx[ii] = max(min(int(round(idx)), size - 1), 0)
        for ii in range(3):
            # sagittal: get to S/A
            # coronal: get to S/L
            # axial: get to A/L
            for jj, vol in enumerate(self._volumes):
                data = np.rollaxis(vol, axis=self._order[ii])[self._data_idx[ii]]
                xax = [1, 0, 0][ii]
                yax = [2, 2, 1][ii]
                if self._order[xax] < self._order[yax]:
                    data = data.T
                if self._flips[xax]:
                    data = data[:, ::-1]
                if self._flips[yax]:
                    data = data[::-1]
                self._ims[jj][ii].set_data(data)
            # deal with crosshairs
            loc = self._data_idx[ii]
            if self._flips[ii]:
                loc = self._sizes[ii] - loc
            loc = [loc] * 2
            if ii == 0:
                self._crosshairs[2]['vert'].set_xdata(loc)
                self._crosshairs[1]['vert'].set_xdata(loc)
            elif ii == 1:
                self._crosshairs[2]['horiz'].set_ydata(loc)
                self._crosshairs[0]['vert'].set_xdata(loc)
            else:  # ii == 2
                self._crosshairs[1]['horiz'].set_ydata(loc)
                self._crosshairs[0]['horiz'].set_ydata(loc)

            self._changing = False

    # Matplotlib handlers ####################################################
    def _in_axis(self, event):
        """Return axis index if within one of our axes, else None"""
        if getattr(event, 'inaxes') is None:
            return None
        for ii, ax in enumerate(self._axes):
            if event.inaxes is ax:
                return ii

    def _on_scroll(self, event):
        """Handle mpl scroll wheel event"""
        assert event.button in ('up', 'down')
        ii = self._in_axis(event)
        if ii is None:
            return
        assert ii in range(4)
        dv = 10. if event.key is not None and 'control' in event.key else 1.
        dv *= 1. if event.button == 'up' else -1.
        dv *= -1 if self._flips[ii] else 1
        val = self._data_idx[ii] + dv
        coords = [self._data_idx[k] for k in range(3)] + [1.]
        coords[ii] = val
        self._set_position(*np.dot(self._affine, coords)[:3])
        self._draw()

    def _on_mouse(self, event):
        """Handle mpl mouse move and button press events"""
        if event.button != 1:  # only enabled while dragging
            return
        ii = self._in_axis(event)
        if ii is None:
            return
        else:
            # translate click xdata/ydata to physical position
            xax, yax = [[1, 2], [0, 2], [0, 1]][ii]
            x, y = event.xdata, event.ydata
            x = self._sizes[xax] - x if self._flips[xax] else x
            y = self._sizes[yax] - y if self._flips[yax] else y
            idxs = [None, None, None, 1.]
            idxs[xax] = x
            idxs[yax] = y
            idxs[ii] = self._data_idx[ii]
            self._set_position(*np.dot(self._affine, idxs)[:3])
        self._draw()

    def _draw(self):
        """Update all four (or three) plots"""
        if self._closed:  # make sure we don't draw when we shouldn't
            return
        for ii in range(3):
            ax = self._axes[ii]
            for im in self._ims:
                ax.draw_artist(im[ii])
            if self._cross:
                for line in self._crosshairs[ii].values():
                    ax.draw_artist(line)
            ax.figure.canvas.blit(ax.bbox)

def check_segmentation(fn_subject):
    from scipy import ndimage
    import matplotlib.pylab as pl
    from matplotlib.colors import ListedColormap
    files = simnibs.SubjectFiles(fn_subject + '.msh')
    T1 = nib.load(files.T1)
    masks = nib.load(files.final_contr).get_fdata()
    lines = np.linalg.norm(np.gradient(masks), axis=0) > 0
    print(lines.shape)
    viewer = NiftiViewer(T1.get_fdata(), T1.affine)
    cmap = pl.cm.jet
    my_cmap = cmap(np.arange(cmap.N))
    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)
    my_cmap = ListedColormap(my_cmap)
    viewer.add_overlay(lines, cmap=my_cmap)
    viewer.show()

if __name__ == '__main__':
    #from simnibs import templates
    #image = nib.load(templates.mni_volume)
    #viewer = NiftiViewer(image.get_fdata(), image.affine)
    #overlay = np.random.random(image.get_fdata().shape)
    #viewer.add_overlay(overlay, cmap='viridis', alpha=0.5)
    #viewer.show()
    check_segmentation('/Users/gbs/simnibs_examples/ernie/ernie')
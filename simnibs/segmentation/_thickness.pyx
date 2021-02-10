#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: language_level=3

from scipy import ndimage as ndi
import numpy as np
import cython

import nibabel

cimport numpy as np


cdef inline np.float_t float_max(np.float_t a, np.float_t b) nogil: return a if a >= b else b
cdef inline np.uint32_t int_max(np.int32_t a, np.int32_t b) nogil: return a if a >= b else b
cdef inline np.uint32_t int_min(np.int32_t a, np.int32_t b) nogil: return a if a <= b else b


def _calc_thickness(label_img):
    ''' Calculates the thichkess of each layer in a 3D binary image'''
    thickness = np.zeros_like(label_img, dtype=np.float32)
    for t in np.unique(label_img):
        if t == 0:
            continue
        else:
            thickness += _thickness_3d_binary_image(
                (label_img == t).astype(np.uint8)
            )
    # If for some reason a voxel had unasigned thickness
    thickness[np.isinf(thickness)] = 1
    return thickness

def _thickness_3d_binary_image(image):
    ''' Calculate thicness in a 3D binary image '''
    thickness = np.zeros_like(image, dtype=np.float32)
    thickness[image > 0] = np.inf
    # Calculate thickness per-slice along 3 different cuts
    # and take the smallest one
    for i in range(3):
        thickness_ax = np.zeros(
            image.swapaxes(0, i).shape, dtype=np.float32
        )
        for j, slice_ in enumerate(image.swapaxes(0, i)):
            if not np.any(slice_):
                continue
            thick_slice = _thickness_slice(slice_)
            thickness_ax[j, ...] = thick_slice
        thickness_ax = thickness_ax.swapaxes(0, i)
        thickness_ax[(thickness_ax < 1e-3) * (image > 0)] = np.inf
        thickness = np.min([thickness, thickness_ax], axis=0)

    return thickness

def _thickness_slice(np.ndarray[np.uint8_t, ndim=2] slice_):
    """ Based on Hilderbrand and Ruegsegger, J. of Microscopy, 1997 """

    if not np.any(slice_):
        return np.zeros_like(slice_, dtype=np.float32)

    ma, distances = medial_axis(
        slice_, return_distance=True
    )

    cdef np.ndarray[np.float32_t, ndim=2] radius_sq = distances.astype(np.float32) ** 2
    cdef np.ndarray[np.int32_t, ndim=2] r = distances.astype(np.int32) + 1
    cdef np.ndarray[np.uint8_t, ndim=2] medial_axis_mask = ma.astype(np.uint8)


    cdef np.ndarray[np.float32_t, ndim=2] thickness = np.zeros_like(
        radius_sq, dtype=np.float32
    )

    cdef np.int32_t x, y
    cdef np.int32_t x_ma, y_ma
    cdef np.int32_t sx = slice_.shape[0]
    cdef np.int32_t sy = slice_.shape[1]

    cdef np.float32_t distance
    cdef np.float_t eps = 1e-3


    # For each voxel in the medial axis
    for x_ma in range(sx):
        for y_ma in range(sy):
            if medial_axis_mask[x_ma, y_ma]:
                # look into the l_infty ball centered around the
                # medial axis voxel
                for x in range(int_max(x_ma - r[x_ma, y_ma], 0),
                               int_min(x_ma + r[x_ma, y_ma] + 1, sx)):
                    for y in range(int_max(y_ma - r[x_ma, y_ma], 0),
                                   int_min(y_ma + r[x_ma, y_ma] + 1, sy)):
                    # if is in the mask
                        if slice_[x, y]:
                            # and in the L2 ball
                            if (x_ma - x)**2 + (y_ma - y)**2 < radius_sq[x_ma, y_ma] + eps:
                                # update thickness
                                thickness[x, y] = float_max(
                                    radius_sq[x_ma, y_ma],
                                    thickness[x, y]
                                )

    return np.sqrt(thickness)

''' All code below has the following license

    Copyright (C) 2019, the scikit-image team
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in
        the documentation and/or other materials provided with the
        distribution.
     3. Neither the name of skimage nor the names of its contributors may be
        used to endorse or promote products derived from this software without
        specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
    IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
    IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
'''

#GBS: Added this function to accelerate medial axis transformation
_medial_axis_lookup_table = None

def _build_medial_axis_lookup_table():
    #
    # Build lookup table - three conditions
    # 1. Keep only positive pixels (center_is_foreground array).
    # AND
    # 2. Keep if removing the pixel results in a different connectivity
    # (if the number of connected components is different with and
    # without the central pixel)
    # OR
    # 3. Keep if # pixels in neighbourhood is 2 or less
    # Note that table is independent of image
    _eight_connect = ndi.generate_binary_structure(2, 2) # GBS: Modified
    global _medial_axis_lookup_table
    center_is_foreground = (np.arange(512) & 2**4).astype(bool)
    _medial_axis_lookup_table = (
        center_is_foreground  # condition 1.
                &
        (np.array([ndi.label(_pattern_of(index), _eight_connect)[1] !=
                   ndi.label(_pattern_of(index & ~ 2**4),
                                _eight_connect)[1]
                   for index in range(512)])  # condition 2 GBS: This is the slowest part
            |
        np.array([np.sum(_pattern_of(index)) < 3 for index in range(512)])) # condition 3
        )



def medial_axis(image, mask=None, return_distance=False):
    """
    Compute the medial axis transform of a binary image
    Parameters
    ----------
    image : binary ndarray, shape (M, N)
        The image of the shape to be skeletonized.
    mask : binary ndarray, shape (M, N), optional
        If a mask is given, only those elements in `image` with a true
        value in `mask` are used for computing the medial axis.
    return_distance : bool, optional
        If true, the distance transform is returned as well as the skeleton.
    Returns
    -------
    out : ndarray of bools
        Medial axis transform of the image
    dist : ndarray of ints, optional
        Distance transform of the image (only returned if `return_distance`
        is True)
    See also
    --------
    skeletonize
    Notes
    -----
    This algorithm computes the medial axis transform of an image
    as the ridges of its distance transform.
    The different steps of the algorithm are as follows
     * A lookup table is used, that assigns 0 or 1 to each configuration of
       the 3x3 binary square, whether the central pixel should be removed
       or kept. We want a point to be removed if it has more than one neighbor
       and if removing it does not change the number of connected components.
     * The distance transform to the background is computed, as well as
       the cornerness of the pixel.
     * The foreground (value of 1) points are ordered by
       the distance transform, then the cornerness.
     * A cython function is called to reduce the image to its skeleton. It
       processes pixels in the order determined at the previous step, and
       removes or maintains a pixel according to the lookup table. Because
       of the ordering, it is possible to process all pixels in only one
       pass.
    Examples
    --------
    >>> square = np.zeros((7, 7), dtype=np.uint8)
    >>> square[1:-1, 2:-2] = 1
    >>> square
    array([[0, 0, 0, 0, 0, 0, 0],
           [0, 0, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 0, 0],
           [0, 0, 1, 1, 1, 0, 0],
           [0, 0, 0, 0, 0, 0, 0]], dtype=uint8)
    >>> medial_axis(square).astype(np.uint8)
    array([[0, 0, 0, 0, 0, 0, 0],
           [0, 0, 1, 0, 1, 0, 0],
           [0, 0, 0, 1, 0, 0, 0],
           [0, 0, 0, 1, 0, 0, 0],
           [0, 0, 0, 1, 0, 0, 0],
           [0, 0, 1, 0, 1, 0, 0],
           [0, 0, 0, 0, 0, 0, 0]], dtype=uint8)

    """
    if mask is None:
        masked_image = image.astype(np.bool)
    else:
        masked_image = image.astype(bool).copy()
        masked_image[~mask] = False

    # Build lookup table
    # GBS: as this part is independent of the image
    # and takes a lot of time, I moved it into a global
    # variable
    if _medial_axis_lookup_table is None:
        _build_medial_axis_lookup_table()
    table = np.ascontiguousarray(_medial_axis_lookup_table, dtype=np.uint8)

    # Build distance transform
    distance = ndi.distance_transform_edt(masked_image)
    if return_distance:
        store_distance = distance.copy()


    # Corners
    # The processing order along the edge is critical to the shape of the
    # resulting skeleton: if you process a corner first, that corner will
    # be eroded and the skeleton will miss the arm from that corner. Pixels
    # with fewer neighbors are more "cornery" and should be processed last.
    # We use a cornerness_table lookup table where the score of a
    # configuration is the number of background (0-value) pixels in the
    # 3x3 neighbourhood
    cornerness_table = np.array([9 - np.sum(_pattern_of(index))
                                 for index in range(512)])
    corner_score = _table_lookup(masked_image, cornerness_table)

    # Define arrays for inner loop
    i, j = np.mgrid[0:image.shape[0], 0:image.shape[1]]
    result = masked_image.copy()
    distance = distance[result]
    i = np.ascontiguousarray(i[result], dtype=np.intp)
    j = np.ascontiguousarray(j[result], dtype=np.intp)
    result = np.ascontiguousarray(result, np.uint8)

    # Determine the order in which pixels are processed.
    # We use a random # for tiebreaking. Assign each pixel in the image a
    # predictable, random # so that masking doesn't affect arbitrary choices
    # of skeletons
    #
    generator = np.random.RandomState(0)
    tiebreaker = generator.permutation(np.arange(masked_image.sum()))
    order = np.lexsort((tiebreaker,
                        corner_score[masked_image],
                        distance))
    order = np.ascontiguousarray(order, dtype=np.int32)

    # Remove pixels not belonging to the medial axis
    _skeletonize_loop(result, i, j, order, table)

    result = result.astype(bool)
    if mask is not None:
        result[~mask] = image[~mask]
    if return_distance:
        return result, store_distance
    else:
        return result

def _pattern_of(index):
    """
    Return the pattern represented by an index value
    Byte decomposition of index
    """
    return np.array([[index & 2**0, index & 2**1, index & 2**2],
                     [index & 2**3, index & 2**4, index & 2**5],
                     [index & 2**6, index & 2**7, index & 2**8]], bool)


def _table_lookup(image, table):
    """
    Perform a morphological transform on an image, directed by its
    neighbors
    Parameters
    ----------
    image : ndarray
        A binary image
    table : ndarray
        A 512-element table giving the transform of each pixel given
        the values of that pixel and its 8-connected neighbors.
    border_value : bool
        The value of pixels beyond the border of the image.
    Returns
    -------
    result : ndarray of same shape as `image`
        Transformed image
    Notes
    -----
    The pixels are numbered like this::
      0 1 2
      3 4 5
      6 7 8
    The index at a pixel is the sum of 2**<pixel-number> for pixels
    that evaluate to true.
    """
    #
    # We accumulate into the indexer to get the index into the table
    # at each point in the image
    #
    if image.shape[0] < 3 or image.shape[1] < 3:
        image = image.astype(bool)
        indexer = np.zeros(image.shape, int)
        indexer[1:, 1:]   += image[:-1, :-1] * 2**0
        indexer[1:, :]    += image[:-1, :] * 2**1
        indexer[1:, :-1]  += image[:-1, 1:] * 2**2

        indexer[:, 1:]    += image[:, :-1] * 2**3
        indexer[:, :]     += image[:, :] * 2**4
        indexer[:, :-1]   += image[:, 1:] * 2**5

        indexer[:-1, 1:]  += image[1:, :-1] * 2**6
        indexer[:-1, :]   += image[1:, :] * 2**7
        indexer[:-1, :-1] += image[1:, 1:] * 2**8
    else:
        indexer = _table_lookup_index(np.ascontiguousarray(image, np.uint8))
    image = table[indexer]
    return image


def _table_lookup_index(np.uint8_t[:, ::1] image):
    """
    Return an index into a table per pixel of a binary image
    Take the sum of true neighborhood pixel values where the neighborhood
    looks like this::
       1   2   4
       8  16  32
      64 128 256
    This code could be replaced by a convolution with the kernel::
      256 128 64
       32  16  8
        4   2  1
    but this runs about twice as fast because of inlining and the
    hardwired kernel.
    """
    cdef:
        Py_ssize_t[:, ::1] indexer
        Py_ssize_t *p_indexer
        np.uint8_t *p_image
        Py_ssize_t i_stride
        Py_ssize_t i_shape
        Py_ssize_t j_shape
        Py_ssize_t i
        Py_ssize_t j
        Py_ssize_t offset

    i_shape   = image.shape[0]
    j_shape   = image.shape[1]
    indexer = np.zeros((i_shape, j_shape), dtype=np.intp)
    p_indexer = &indexer[0, 0]
    p_image   = &image[0, 0]
    i_stride  = image.strides[0]
    assert i_shape >= 3 and j_shape >= 3, \
        "Please use the slow method for arrays < 3x3"
    with nogil:
        for i in range(1, i_shape-1):
            offset = i_stride* i + 1
            for j in range(1, j_shape - 1):
                if p_image[offset]:
                    p_indexer[offset + i_stride + 1] += 1
                    p_indexer[offset + i_stride] += 2
                    p_indexer[offset + i_stride - 1] += 4
                    p_indexer[offset + 1] += 8
                    p_indexer[offset] += 16
                    p_indexer[offset - 1] += 32
                    p_indexer[offset - i_stride + 1] += 64
                    p_indexer[offset - i_stride] += 128
                    p_indexer[offset - i_stride - 1] += 256
                offset += 1
        #
        # Do the corner cases (literally)
        #
        if image[0, 0]:
            indexer[0, 0] += 16
            indexer[0, 1] += 8
            indexer[1, 0] += 2
            indexer[1, 1] += 1

        if image[0, j_shape - 1]:
            indexer[0, j_shape - 2] += 32
            indexer[0, j_shape - 1] += 16
            indexer[1, j_shape - 2] += 4
            indexer[1, j_shape - 1] += 2

        if image[i_shape - 1, 0]:
            indexer[i_shape - 2, 0] += 128
            indexer[i_shape - 2, 1] += 64
            indexer[i_shape - 1, 0] += 16
            indexer[i_shape - 1, 1] += 8

        if image[i_shape - 1, j_shape - 1]:
            indexer[i_shape - 2, j_shape - 2] += 256
            indexer[i_shape - 2, j_shape - 1] += 128
            indexer[i_shape - 1, j_shape - 2] += 32
            indexer[i_shape - 1, j_shape - 1] += 16
        #
        # Do the edges
        #
        for j in range(1, j_shape - 1):
            if image[0, j]:
                indexer[0, j - 1] += 32
                indexer[0, j] += 16
                indexer[0, j + 1] += 8
                indexer[1, j - 1] += 4
                indexer[1, j] += 2
                indexer[1, j + 1] += 1
            if image[i_shape - 1, j]:
                indexer[i_shape - 2, j - 1] += 256
                indexer[i_shape - 2, j] += 128
                indexer[i_shape - 2, j + 1] += 64
                indexer[i_shape - 1, j - 1] += 32
                indexer[i_shape - 1, j] += 16
                indexer[i_shape - 1, j + 1] += 8

        for i in range(1, i_shape - 1):
            if image[i, 0]:
                indexer[i - 1, 0] += 128
                indexer[i, 0] += 16
                indexer[i + 1, 0] += 2
                indexer[i - 1, 1] += 64
                indexer[i, 1] += 8
                indexer[i + 1, 1] += 1
            if image[i, j_shape - 1]:
                indexer[i - 1, j_shape - 2] += 256
                indexer[i, j_shape - 2] += 32
                indexer[i + 1, j_shape - 2] += 4
                indexer[i - 1, j_shape - 1] += 128
                indexer[i, j_shape - 1] += 16
                indexer[i + 1, j_shape - 1] += 2
    return np.asarray(indexer)


def _skeletonize_loop(np.uint8_t[:, ::1] result,
                      Py_ssize_t[::1] i, Py_ssize_t[::1] j,
                      np.int32_t[::1] order, np.uint8_t[::1] table):
    """
    Inner loop of skeletonize function
    Parameters
    ----------
    result : ndarray of uint8
        On input, the image to be skeletonized, on output the skeletonized
        image.
    i, j : ndarrays
        The coordinates of each foreground pixel in the image
    order : ndarray
        The index of each pixel, in the order of processing (order[0] is
        the first pixel to process, etc.)
    table : ndarray
        The 512-element lookup table of values after transformation
        (whether to keep or not each configuration in a binary 3x3 array)
    Notes
    -----
    The loop determines whether each pixel in the image can be removed without
    changing the Euler number of the image. The pixels are ordered by
    increasing distance from the background which means a point nearer to
    the quench-line of the brushfire will be evaluated later than a
    point closer to the edge.
    Note that the neighbourhood of a pixel may evolve before the loop
    arrives at this pixel. This is why it is possible to compute the
    skeleton in only one pass, thanks to an adapted ordering of the
    pixels.
    """
    cdef:
        Py_ssize_t accumulator
        Py_ssize_t index, order_index
        Py_ssize_t ii, jj
        Py_ssize_t rows = result.shape[0]
        Py_ssize_t cols = result.shape[1]

    with nogil:
        for index in range(order.shape[0]):
            accumulator = 16
            order_index = order[index]
            ii = i[order_index]
            jj = j[order_index]
            # Compute the configuration around the pixel
            if ii > 0:
                if jj > 0 and result[ii - 1, jj - 1]:
                    accumulator += 1
                if result[ii - 1, jj]:
                    accumulator += 2
                if jj < cols - 1 and result[ii - 1, jj + 1]:
                        accumulator += 4
            if jj > 0 and result[ii, jj - 1]:
                accumulator += 8
            if jj < cols - 1 and result[ii, jj + 1]:
                accumulator += 32
            if ii < rows - 1:
                if jj > 0 and result[ii + 1, jj - 1]:
                    accumulator += 64
                if result[ii + 1, jj]:
                    accumulator += 128
                if jj < cols - 1 and result[ii + 1, jj + 1]:
                    accumulator += 256
            # Assign the value of table corresponding to the configuration
            result[ii, jj] = table[accumulator]


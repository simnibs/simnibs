#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Copyright (C) 2022  Kristoffer H. Madsen, Oula Puonti

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

    The included brainsprite.js uses the MIT License (MIT):

    Copyright (c) 2016-2020 Pierre Bellec and the Brainsprite contributors

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
"""

import nibabel as nib
import numpy as np
import scipy.ndimage as ndi
from io import BytesIO
from base64 import b64encode
import json
import os
from PIL import Image
from string import Template

def mosaic(data):
    '''
    Make 2D mosaic of 3D input data

    Parameters
    ----------
    data : Numpy array
        Input data
    Returns
    -------
    mosaic : 2D numpy array
        If each sagittal slice is nz (height) x ny (width) pixels, the sprite
        size is (M x nz) x (N x ny), where M and N are computed to be roughly
        equal. All slices are pasted together row by row, from top left to
        bottom right. The last row is completed with empty slices.
    '''
    data = np.ascontiguousarray(data)
    nx, ny, nz = data.shape
    nrows = int(np.ceil(np.sqrt(nx)))
    ncolumns = int(np.ceil(nx / float(nrows)))

    mosaic = np.zeros((nrows * nz, ncolumns * ny))
    indrow, indcol = np.where(np.ones((nrows, ncolumns)))

    for xx in range(nx):
        for i in range(len(data)):
           mosaic[(indrow[xx] * nz):((indrow[xx] + 1) * nz), (indcol[xx] * ny):
                ((indcol[xx] + 1) * ny)] = data[xx, :, ::-1].transpose()
    return mosaic


def getRASdata(nii, returncoords=False):
    '''
    Return near RAS orientation data without reslicing

    Parameters
    ----------
    nii : nibabel image instance
        input volume.
    returncoords : boolean
        return meshgrid with sampled coordinates in world space.
    Returns
    -------
    data : ndarray
        data in near RAS space.
    aff : corresponding affine transformation matrix to world space
    coord : ndarray (only if returncoords=True)
    world space coordinates sampled if requested.
    '''

    orient = (('R','A','S'),('L','P','I'))

    axcodes = nib.aff2axcodes(nii.affine)
    #-1 will indicate axis flip
    flip = [1, 1, 1]
    order = [0,1,2]
    for i in range(3):
        order[i] = np.nonzero([ax==orient[0][i] for ax in axcodes])[0]
        if len(order[i])==0:
            order[i] = np.nonzero([ax==orient[1][i] for ax in axcodes])[0][0]
            flip[i] = -1
        else:
            order[i]=order[i][0]
    data = np.transpose(nii.dataobj,order)[::flip[0],::flip[1],::flip[2]]
    aff = np.zeros((4,4))
    for i in range(3):
        aff[order[i],i]=flip[i]
    aff[3,3]=1
    aff=aff@nii.affine
    if returncoords:
        voxcoords = np.meshgrid(*[np.arange(s) for s in nii.shape[:3]],
                             indexing='ij')
        worldcoords = aff[:3,:3]@np.array(voxcoords).reshape(3,-1)+aff[:3,3,None]
        return (data, aff, worldcoords)
    else:
        return (data, aff)

def resample_world_iso(nii, maxsize=256, order=1):
    corners = np.array(np.meshgrid(
        *[np.linspace(0,i,2) for i in nii.shape]+[1,],indexing='ij')
        ).reshape(4,-1)
    corners = (nii.affine@corners)[:3]
    bb = (np.min(corners, axis=1), np.max(corners, axis=1))
    res = np.max(np.diff(bb, axis=0) / maxsize)
    coord = np.meshgrid(*[np.arange(bb[0][i], bb[1][i], res)
                          for i in range(3)], indexing='ij')
    data = resample_at_coords(nii, coord)
    M = np.diag(res*np.ones(4))
    M[3, 3] = 1
    M[:3, 3] = bb[0]
    return (data, M, coord)

def resample_at_coords(nii, coord, output_shape, order=1):
    '''
    Resample data at given world space coordinates using map_coordinates

    Parameters
    ----------
    nii : nibabel image
        input image.
    coord : ndarray
        world space coordinates to sample.
    order : int
        interpolation order.
    Returns
    -------
    data : ndarray
        resampled data.

    '''
    iM = np.linalg.inv(nii.affine)
    vox_coord = iM[:3, :3] @ np.array(coord).reshape(3, -1) + iM[:3, 3][:, None]
    data = ndi.map_coordinates(nii.dataobj, vox_coord,
                                   order=order).reshape(
                                       coord[0].shape)
    print(coord.shape)
    return data.reshape(output_shape)

def niis_to_mosaic(niis, interpolation_order, maxsize=256):
    '''
    Convert volumes to mosaic

    Parameters
    ----------
    niis : list of nibabel volumes
        input images.
    maxsize : int, optional
        maximum size of slice in voxels. The default is 256.

    Returns
    -------
    list of ndarrays
        mosaic'ed data.

    '''
    if np.max(niis[0].shape) > maxsize:
        data0, affine, coord = resample_world_iso(niis[0], order=interpolation_order[0], maxsize=maxsize)
        data = np.zeros((len(niis),) + data0.shape)
        data[0] = data0
        for i in range(1,len(niis)):
            data[i] = resample_at_coords(niis[i], coord, order=interpolation_order[i])
    else:
        all_equal = True
        for i in range(1,len(niis)):
            if not np.allclose(niis[i-1].affine, niis[i].affine):
                all_equal = False
                break
        if all_equal:
            data = [getRASdata(nii) for nii in niis]
            affine = data[0][1]
            data = [d[0] for d in data]
        else:
            data0, affine, coord = getRASdata(niis[0], returncoords=True)
            data = np.zeros((len(niis),)+(data0.shape))
            data[0] = data0
            for i in range(1,len(niis)):
                data[i] = resample_at_coords(niis[i], coord, data0.shape, order=interpolation_order[i])
    imgs = [mosaic(d) for d in data]
    return (imgs, data[0].shape, affine)

def _get_im_data(im, filename=None):
    if filename is None:
        sprite = BytesIO()
        im.save(sprite,format='WebP', lossless=True, method=6)
        sprite.seek(0)
        data = b64encode(sprite.read()).decode('utf-8')
        sprite.close()
        return data
    else:
        im.save(filename, format='WebP', lossless=True, method=6)

def form_brainsprite(imgs, shape, affine, template, js_query_path,
                     brain_sprite_path, names=None, imgPath=None,
                     select=(0,1)):
    if imgPath is None:
        embed = True
    else:
        embed = False

    cfont = '#FFFFFF'
    cbg = '#000000'
    params = {'canvas': '3Dviewer',
              'sprite': 'spriteImg',
              'overlay': {'sprite': 'overlayImg',
                           'nbSlice': {'X': shape[0],
                                       'Y': shape[1],
                                       'Z': shape[2]},
                           'opacity': 1},
              'nbSlice': {'X': shape[0],
                          'Y': shape[1],
                          'Z': shape[2]},
              'colorBackground': cbg,
              'colorFont': cfont,
              'crosshair': True,
              'affine': affine.tolist(),
              'flagCoordinates': True,
              'title': None,
              'flagValue': False,
              'numSlice': {'X': shape[0]//2,
                           'Y': shape[1]//2,
                           'Z': shape[2]//2}}

    if names is None:
        names = [f'{i+1}' for i in range(len(imgs))]

    jsonv = dict.fromkeys(['params', 'js_jquery', 'js_brainsprite'])
    #convert to b64 encoded image data
    data = []
    dropd1 = []
    dropd2 = []
    for i,name in enumerate(names):
        if embed:
            data.append(''.join(
                (f'<img id="{name}" class="hidden" src="data:image/webp;base64,',
                _get_im_data(imgs[i]), f'" alt="{name}" />\n')))
        else:
            fn = os.path.join(imgPath,f'{name}.webp')
            _get_im_data(imgs[i],fn)
            data.append(''.join(
                (f'<img id="{name}" class="hidden" src="{fn}" alt="{name}" />\n')))
        if i==select[0]:
            dropd1.append(f'<option selected value="img{i}">{name}</option>\n')
        else:
            dropd1.append(f'<option value="img{i}">{name}</option>\n')
        if i==select[1]:
            dropd2.append(f'<option selected value="img{i}">{name}</option>\n')
        else:
            dropd2.append(f'<option value="img{i}">{name}</option>\n')

    jsonv['data_base64'] = ''.join(data)
    jsonv['dropdown1'] = ''.join(dropd1)
    jsonv['dropdown2'] = ''.join(dropd2)
    jsonv['params'] = params
    jsonv['params'] = json.dumps(jsonv['params'])

    with open(js_query_path) as f:
       jsonv['js_jquery'] = f.read()
    with open(brain_sprite_path) as f:
       jsonv['js_brainsprite'] = f.read()

    with open(template, 'rb') as f:
        template = Template(f.read().decode('utf-8'))
    template=template.safe_substitute(jsonv)
    return template

def write_viewer(images, color_maps, interpolation_order, viewer,
                 template, jquery, brainsprite, names=None, imgPath=None,
                 select=(0,1)):
    import time
    t0=time.time()
    imgs,shape,affine=niis_to_mosaic(images, interpolation_order, maxsize=256)
    print(f'time to form overlays: {time.time()-t0:.1f}s')
    im=[]
    for img, cmap, order in zip(imgs, color_maps, interpolation_order):
        if order == 0:
            img = np.uint8(img)
        im.append(Image.fromarray(np.uint8(cmap(img)*255)))
    t0=time.time()
    doc = form_brainsprite(im, shape, affine, template, jquery, brainsprite,
                           names=names, imgPath = imgPath)
    print(f'time to compress images and create brainsprite: {time.time()-t0:.1f}s')
    with open(viewer,'wb') as f:
        f.write(doc.encode('utf-8'))

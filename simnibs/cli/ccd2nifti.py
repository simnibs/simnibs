# -*- coding: utf-8 -*-
'''
    command line tool to convert coil dipole definition ccd files to nifti1
    format. This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2021  Kristoffer H. Madsen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

'''

import os
import glob
import numpy as np
import fmm3dpy
import simnibs.simulation.coil_numpy as coil_numpy
import nibabel as nib
import re

def parseccd(ccd_file):
    '''
    Parse ccd file, and return intended bounding box, resolution and maxdIdt
        if availible
    
    Parameters
    ----------
    ccd_file : string
        ccd file to parse
    '''
    
    d_position, d_moment = coil_numpy.read_ccd(ccd_file)
    f = open(ccd_file,'r')
    data = f.readline()
    bb = []
    for dim in ('x','y','z'):
        a = re.search('%s=[-_.\d,]*'%dim, data)
        if a is None:
            bb.append(None)
        else:
            a = np.fromstring(data[a.start()+2:a.end()],sep=',')
            if len(a)<2:
                bb.append(-np.abs(a),np.abs(a))
            else:
                bb.append(a)
    res = []
    a = re.search('%s=[-_.\d,]*'%'resolution', data)
    if a is None:
        res.append(None)
    else:
        a = np.fromstring(data[a.start()+11:a.end()],sep=',')
        if len(a)<3:
            for i in range(len(a),3):
                a = np.concatenate((a, (a[i-1],)))
        res = a
    a = re.search('%s=[-_.\d,]*'%'maxdIdt', data)
    if a is None:
        maxdIdt = None
    else:
        maxdIdt = np.fromstring(data[a.start()+8:a.end()],sep=',')      
    return d_position, d_moment, bb, res, maxdIdt


def ccd2nifti(ccdfn, resolution=None, boundingbox=None, dIdt=None, eps=1e-3):
    '''
    Convert CCD coil dipole files to nifti1 format

    Parameters
    ----------
    ccdfn : string
        CCD file.
    resolution : ndarray, optional
        Resolution (dx,dy,dz). The default is (3,3,3).
    boundingbox : ndarray, optional
        Bounding box ((xmin,xmax),(ymin,ymax),(zmin,zmax)). 
        The default is ((-300, 300), (-200, 200), (0, 300)).
    dIdt : float, optional
        Maximum dIdt. The default is None.
    eps : float, optional
        precision for fmm3dpy. The default is 1e-3.

    Returns
    -------
    nii : Nifti1Volume
        Nifti1Volume containing dA/dt field.

    '''
    
    #read and parse ccd file
    d_position, d_moment, res, bb, maxdIdt = parseccd(ccdfn)
    if boundingbox is None:
        bb = boundingbox
    if resolution is None:
        res = resolution
    if res is None:
        res = np.array((3., 3., 3.))
    if bb is None:
        bb = np.array(((-300, 300), (-200, 200), (0, 300)))
    if dIdt is None:
        dIdt = maxdIdt
    #create grid
    eps = np.spacing(1e4)
    x = np.arange(bb[0][0], bb[0][1] + eps, res[0]) #xgrid
    y = np.arange(bb[1][0], bb[1][1] + eps, res[1]) #ygrid
    z = np.arange(bb[2][0], bb[2][1] + eps, res[2]) #zgrid
    xyz = np.meshgrid(x, y, z, indexing='ij') #grid
    #reshape to 2D
    xyz = np.array(xyz).reshape((3, len(x) * len(y) * len(z)))
    xyz *= 1.0e-3 #from mm to SI (meters)
    #use fmm3dpy to calculate expansion fast
    A = np.zeros((xyz.shape[1], 3), dtype=float)
    out = [
        fmm3dpy.lfmm3d(
            eps=eps,
            sources=d_position.T,
            charges=d_m,
            targets=xyz,
            pgt=2
        )
        for d_m in d_moment.T
    ]
    A[:, 0] = (out[1].gradtarg[2] - out[2].gradtarg[1])
    A[:, 1] = (out[2].gradtarg[0] - out[0].gradtarg[2])
    A[:, 2] = (out[0].gradtarg[1] - out[1].gradtarg[0])
    A *= -1e-7
 
    A = A.reshape((len(x), len(y), len(z), 3))
    #header info
    hdr = nib.Nifti1Header()
    hdr.set_data_dtype(np.float32)
    hdr.set_xyzt_units('mm','unknown')
    #affine matrix
    M = np.identity(4) * np.array((res[0], res[1], res[2], 1))
    M[0, 3] = x[0]
    M[1, 3] = y[0]
    M[2, 3] = z[0]
    #create nifti1 volume
    nii = nib.Nifti1Image(A, M, hdr)
    #set maxdIdt if availible
    if not dIdt is None:
        nii.header['descrip'] = 'maxdIdt=%f' % maxdIdt
    return nii

def convert_recursive(indir, force=False):
    '''
    Converts all ccd files within input directory recursively

    Parameters
    ----------
    indir : string
        input directory.
    force : boolean, optional
        Force overwriting nifti1 files. The default is False.

    '''
    ccd_files = glob.iglob(os.path.join(indir, '**', '*.ccd'),
                          recursive=True)
    for ccd_file in ccd_files:
        if len(glob.glob(os.path.splitext(ccd_file)[0]
                         + '.nii*')) == 0 or force:
            print('Converting %s to nifti1 format' % ccd_file)
            nii = ccd2nifti(ccd_file)
            nii.to_filename(os.path.splitext(ccd_file)[0] + '.nii.gz')
        else:
            print('Nifti1 version of %s already exists' % ccd_file)

if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Convert CCD files to Nifti1 format.')
    parser.add_argument('-i', '--infile', dest='infile', default=None, required=True,
                        help='CCD file to convert')
    parser.add_argument('-o', '--outfile', dest='outfile', default=None,
                    help='sum the integers (default: find the max)')
    options = parser.parse_args(sys.argv[1:])
    if os.path.isfile(options.infile):
        nii = ccd2nifti(options.infile)
        if options.outfile is None:
            options.outfile = os.path.splitext(options.infile)[0] + '.nii.gz'
        nii.to_filename(options.outfile)
    else:
        print('Cannot locate input file: %s' % options.infile)

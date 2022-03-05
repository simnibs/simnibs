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
import time
from simnibs.simulation.coil_numpy import parseccd

def writeccd(fn, mpos, m, info=None, extra=None):
    N=m.shape[0]
    f=open(fn,'w')
    f.write('# %s version 1.0;'%os.path.split(fn)[1])
    f.write('# %s version 1.1;'%os.path.split(fn)[1])
    if not info is None:
        for i,key in enumerate(info.keys()):
            f.write(f'{key}={info[key]};')
    if not extra is None:
	    f.write(f'{extra};')
    f.write('\n')
    f.write('%i\n'%N)
    f.write('# centers and weighted directions of the elements (magnetic dipoles)\n')
    for i in range(N):
        f.write('%.15e %.15e %.15e '%tuple(mpos[i]))
        f.write('%.15e %.15e %.15e\n'%tuple(m[i]))
    f.close()

def rescale_ccd(ccd_file, outname=None):
    d_position, d_moment, bb, res, info = parseccd(ccd_file)
    try:
        scale = np.fromstring(info['dIdtmax'],sep=',')[0] / \
            np.fromstring(info['dIdtstim'],sep=',')[0]
    except:
        raise ValueError('cannot find dIdtstim and/or dIdtmax field(s) in ccd file')
    #rescale dipole moment
    d_moment *= scale
    #set maxdIdt and remove dIdtstim
    info['dIdtmax'] = info['dIdtstim']
    info.pop('dIdtstim')

    if outname is None:
        inname = os.path.split(ccd_file)
        outname = os.path.join(inname[0],
                               os.path.splitext(inname[1])[0] +
                               '_rescaled' + '.ccd')
    writeccd(outname, d_position, d_moment, info)

def _lfmm3d(charges, eps, sources, targets, nd=1):
    '''
    Wrapper function for fmm3dpy.lfmm3d
    Parameters
    ----------
    charges : ndarray
        Charges
    eps : float
        precision
    sources : ndarray
        Sources.
    targets : ndarray
        Targets

    Returns
    -------
    fmm3dpy output object

    '''
    return fmm3dpy.lfmm3d(eps=eps, sources=sources, charges=charges,
                          targets=targets, pgt=2, nd=nd)

def ccd2nifti(ccdfn, info={}, eps=1e-3, Bfield=False):
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
    d_position, d_moment, boundingbox, resolution, info = parseccd(ccdfn)
    if not boundingbox[0] is None:
        bb = boundingbox
    else:
	    bb = np.array(((-300, 300), (-200, 200), (0, 300)))
    if not resolution[0] is None:
        res = resolution
    else:
        res = np.array((3., 3., 3.))
    #create grid
    dx = np.spacing(1e4)
    x = np.arange(bb[0][0], bb[0][1] + dx, res[0]) #xgrid
    y = np.arange(bb[1][0], bb[1][1] + dx, res[1]) #ygrid
    z = np.arange(bb[2][0], bb[2][1] + dx, res[2]) #zgrid
    xyz = np.meshgrid(x, y, z, indexing='ij') #grid
    #reshape to 2D
    xyz = np.array(xyz).reshape((3, len(x) * len(y) * len(z)))
    xyz *= 1.0e-3 #from mm to SI (meters)
    if Bfield:
        A = B_from_dipoles(d_moment, d_position, xyz.T)
    else:
        A = A_from_dipoles(d_moment, d_position, xyz.T)
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
    #set dIdtmax if availible
    try:
        dstr = f"dIdtmax={info['dIdtmax']}"
        try:
            dstr += f";coilname={info['coilname']}"
        except:
            pass
        try:
            dstr += f";stimulator={info['stimulator']}"
        except:
            pass
        try:
            dstr += f";brand={info['brand']}";
        except:
            pass
        nii.header['descrip'] = dstr
    except:
        print('no information on dIdtmax found omitting from nii file.')
    return nii

def A_from_dipoles(d_moment, d_position, target_positions, eps=1e-3, direct='auto'):
    #if set to auto use direct methods if # dipoles less than 300
    if direct=='auto':
        if d_moment.shape[0]<300:
            direct = True
        else:
            direct = False
    if direct is True:
        out = fmm3dpy.l3ddir(charges=d_moment.T, sources=d_position.T,
                  targets=target_positions.T, nd=3, pgt=2)
    elif direct is False:
        #use fmm3dpy to calculate expansion fast
        out = fmm3dpy.lfmm3d(charges=d_moment.T, eps=eps, sources=d_position.T,
                  targets=target_positions.T, nd=3, pgt=2)
    else:
        print('Error: direct flag needs to be either "auto", True or False')
    A = np.empty((target_positions.shape[0], 3), dtype=float)
    #calculate curl
    A[:, 0] = (out.gradtarg[1][2] - out.gradtarg[2][1])
    A[:, 1] = (out.gradtarg[2][0] - out.gradtarg[0][2])
    A[:, 2] = (out.gradtarg[0][1] - out.gradtarg[1][0])
    #scale
    A *= -1e-7
    return A

def B_from_dipoles(d_moment, d_position, target_positions, eps=1e-3, direct='auto'):
    #if set to auto use direct methods if # dipoles less than 300
    if direct=='auto':
        if d_moment.shape[0]<300:
            direct = True
        else:
            direct = False
    if direct is True:
        out = fmm3dpy.l3ddir(dipvec=d_moment.T, sources=d_position.T,
                  targets=target_positions.T, nd=1, pgt=2)
    elif direct is False:
        out = fmm3dpy.lfmm3d(dipvec=d_moment.T, eps=eps, sources=d_position.T,
                  targets=target_positions.T, nd=1, pgt=2)
    else:
        print('Error: direct flag needs to be either "auto", True or False')
    B = out.gradtarg.T
    B *= -1e-7
    return B

def main():
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Convert CCD files to Nifti1 format')
    parser.add_argument('-i', '--infile', dest='infile', default=None, required=True,
                        help='CCD file to convert')
    parser.add_argument('-o', '--outfile', dest='outfile', default=None,
                    help='output filename, will default to replacing extension with .nii.gz')
    parser.add_argument('-r', '--rescale', dest='rescale', action='store_true',
                        help='Rescale CCD file according to stimulator reported'
                           'dI/dt (writes new ccd file - with suffix _rescaled)')
    parser.add_argument('-f', '--force', dest='force', action='store_true',
                        help='Force rewrite')
    parser.add_argument('-b', '--bfield', dest='Bfield', action='store_true',
                        help='Write B field instead of A field')

    options = parser.parse_args(sys.argv[1:])
    if os.path.isdir(options.infile):
        print(f'recursively processing CCD files in {options.infile}')
        ccd_files = glob.iglob(os.path.join(options.infile, '**', '*.ccd'),
                          recursive=True)
        options.outfile = None
    elif os.path.isfile(options.infile):
        ccd_files = (options.infile,)
    else:
        print(f'Cannot locate input file: {options.infile}')
    for ccdfile in ccd_files:
        if options.rescale:
            try:
                rescale_ccd(ccdfile,options.outfile)
                print(f'Successfully rescaled {ccdfile}')
            except:
                print(f'Rescaling {ccdfile} failed, check if dIdtstim exists'
                      'and that only one stimulator value is present')
        else:
            if options.outfile is None:
                outfile = os.path.splitext(ccdfile)[0] + '.nii.gz'
            else:
                outfile=options.outfile
            if len(glob.glob(os.path.splitext(outfile)[0] + '*')) == 0 or options.force:
                t0 = time.perf_counter()
                print(f'expanding CCD file {ccdfile}')
                nii = ccd2nifti(ccdfile, Bfield=options.Bfield)
                nii.to_filename(outfile)
                print(f'Time spend: {time.perf_counter()-t0:.0f}s')
            else:
                print(f'Nifti1 version of {ccdfile} already exists')

if __name__ == '__main__':
    main()

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
# import multiprocessing
# import functools

def parseccd(ccd_file):
    '''
    Parse ccd file, and return intended bounding box, resolution and dIdtmax
        and other fields if available
    
    Parameters
    ----------
    ccd_file : string
        ccd file to parse
    
    On ccd file format version 1.1:
    1) First line is a header line escaped with # which can contain any number
       of variables in the form variable=value, a few of these variables are 
       reserved as can be seen below, they are separated by semicolons (;).
    2) The second line is the contains the number of dipoles expected
       (this is actually not used in practice).
    3) Third line contains a header text excaped by #, typically:
       # centers and weighted directions of the elements (magnetic dipoles)
    4-end) Remaining lines are space separated dipole positions and dipole
       moments in a string number format readable by numpy. E.g. each line must contain 
       six values: x y z mx my mz, where the first three are x,y and positions
       in meters, and the remaining three are dipole moments in x,y and z direction
       in Coulumb * meter per 1 A/s input current.
       an example could be:
       0 0 0 0 0 1.0e-03
       indicating a dipole at position 0,0,0 in z direction with strength
       0.001 C*m*s/A
     
    The variables are used to encode additonal optional information in text, specifically:
    dIdtmax=147,100
        which would indicate a max dI/dt (at 100% MSO) of 146.9 A/microsecond
        for first stimulator and 100 A/microsecond for the second stimulator. 
    dIdtstim=162
        Indicating the max dI/dt reported on the stimulation display of 162, 
        typically this is used to create a rescaled version of the ccd file 
        such that the stimulator reported dI/dt max can be used directly.
        This is currently only supported for one stimulator.
    stimulator=Model name 1,Model name 2
    brand=Brand name 1,Brand name 2
    coilname=name of coil
        Indicates the name of the coil for display purposes
    Some variables are used for expansion in to nifti1 format:
    x=-300,300
        Indicates that the ccd file should be expanded into a FOV from 
        x=-300mm to x=300mm, this could also be indicated as x=300
    y=-300,300
        The same for y
    z=-200,200
        The same for z
    resolution=3,3,3
        Indicates that the resolution should be 3mm in x,y and z directions,
        this could also be given as resolution=3
    
    The below is an example header line:
    #Test CCD file;dIdtmax=162;x=300;y=300;z=200;resolution=3;stimulator=MagProX100;brand=MagVenture;
        
        
    '''
    def parseField(info, field):
        try:
            a = np.fromstring(info[field],sep=',')
        except:
            a = None
        return a
    
    d_position, d_moment = coil_numpy.read_ccd(ccd_file)
    
    #reopen to read header
    f = open(ccd_file,'r')
    data = f.readline()
    fields = re.findall('(\w*=[^;]*)', data + ';')
    labels = [f.split('=')[0] for f in fields]
    values = [f.split('=')[1] for f in fields]
    info = {}
    for i,label in enumerate(labels):
        info[label]=values[i]
    
    #parse bounding box for nii
    bb = []
    for dim in ('x','y','z'):
        a = parseField(info, dim)
        if a is None:
            bb.append(None)
        else:
            if len(a)<2:
                bb.append((-np.abs(a),np.abs(a)))
            else:
                bb.append(a)
    
    #parse resolution
    res = []
    a = parseField(info, 'resolution')
    if a is None:
        res.append(None)
    else:
        if len(a)<3:
            for i in range(len(a),3):
                a = np.concatenate((a, (a[i-1],)))
        res = a
        
    
    return d_position, d_moment, bb, res, info

def writeccd(fn, mpos, m, info=None, extra=None):
    N=m.shape[0]
    f=open(fn,'w')
<<<<<<< HEAD
    f.write('# %s version 1.0;'%os.path.split(fn)[1])
=======
    f.write('# %s version 1.1;'%os.path.split(fn)[1])
>>>>>>> 7d1d615cac139da8265fb72baaf78dffe935ab58
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
    
def _lfmm3d(charges, eps, sources, targets):
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
                          targets=targets, pgt=2)

def ccd2nifti(ccdfn, info={}, eps=1e-3):
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
<<<<<<< HEAD
	    bb = np.array(((-300, 300), (-200, 200), (0, 300)))
=======
        bb = np.array(((-300, 300), (-200, 200), (0, 300)))
>>>>>>> 7d1d615cac139da8265fb72baaf78dffe935ab58
    if not resolution[0] is None:
        res = resolution
    else:
        res = np.array((3., 3., 3.))
<<<<<<< HEAD

=======
    
>>>>>>> 7d1d615cac139da8265fb72baaf78dffe935ab58
    #create grid
    eps = np.spacing(1e4)
    x = np.arange(bb[0][0], bb[0][1] + eps, res[0]) #xgrid
    y = np.arange(bb[1][0], bb[1][1] + eps, res[1]) #ygrid
    z = np.arange(bb[2][0], bb[2][1] + eps, res[2]) #zgrid
    xyz = np.meshgrid(x, y, z, indexing='ij') #grid
    #reshape to 2D
    xyz = np.array(xyz).reshape((3, len(x) * len(y) * len(z)))
    xyz *= 1.0e-3 #from mm to SI (meters)
    A = np.empty((xyz.shape[1], 3), dtype=float)
    
    #use fmm3dpy to calculate expansion fast
    #create a pool of workers to process each dimension separately
    #to accomplish the same but using multiprocessing:
    #pool = multiprocessing.Pool()
    #f = functools.partial(_lfmm3d, eps=eps, sources=d_position.T, targets=xyz)
    #out = pool.map(f, d_moment.T)
    
    #non-parallel version (note, fmm3d is already somewhat parallel)
    out = [_lfmm3d(charges=d_m, eps=eps, sources=d_position.T, targets=xyz)
        for d_m in d_moment.T]

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
            if len(glob.glob(os.path.splitext(outfile)[0]
                         + '.nii*')) == 0 or options.force:
                t0 = time.perf_counter()
                print(f'expanding CCD file {ccdfile}')
                nii = ccd2nifti(ccdfile)
                nii.to_filename(outfile)
                print(f'Time spend: {time.perf_counter()-t0:.0f}s')
            else:
                print(f'Nifti1 version of {ccdfile} already exists')        

if __name__ == '__main__':
    main()
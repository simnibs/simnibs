# -*- coding: utf-8 -*-\
'''
    Convert vtk text files to simNIBS mesh objects
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2022 Kristoffer H Madsen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import numpy as np
from simnibs.mesh_tools import mesh_io
import gzip

def read_itk_tetrahedron(fn):
    with gzip.open(fn,'rt') as fid:
        lines=fid.read().split('\n')
    keys = ('Number of points:','Number of cells:','Number of labels:')
    vals = [0 for key in keys]
    j = 0
    i = 0
    for line in lines:
        while j<len(keys):
            if line.find(keys[j])>=0:
                #print(line[line.find(keys[j])+len(keys[j]):])
                vals[j] = int(line[line.find(keys[j])+len(keys[j]):])
                break
            j += 1
        if line.find('Reference position')>=0:
            break
        i += 1
    vertices = np.fromstring(''.join(lines[i+1:i+1+vals[0]]), sep=' ', dtype='float32').reshape(vals[0], 4)
    i += vals[0] + 1
    
    i += next(j for j,line in enumerate(lines[i+1:]) if line.startswith('Cells:'))
    
    ii = next(j for j,line in enumerate(lines[i+1:]) if line.startswith('Point parameters:'))
    tet = np.fromstring(' '.join([' '.join(line.split()[2:]) 
                                for line in lines[i+1:i+1+ii] if line.find('TETRAHEDRON')>=0]
                               ), dtype=int, sep=' ').reshape(-1,4)

    i += ii
    vertidx = vertices[:, 0].astype('int')
    tet = np.searchsorted(vertidx, tet)
    vertices = vertices[:, 1:]
    
    nodedata = np.array([np.array(line.split()[:-4], dtype='float32') for line in lines[i+2:i+2+vals[0]]])
    
    assert np.allclose(nodedata[:,0],vertidx)
    return vertices, tet, nodedata[:,1:]

def itk_to_msh(fn):
    vertices,tetrahedra,nodedata=read_itk_tetrahedron(fn)
    msh = mesh_io.Msh()
    msh.elm = mesh_io.Elements(tetrahedra=tetrahedra + 1)
    msh.nodes = mesh_io.Nodes(vertices)	
    msh.add_node_field(nodedata,'tissue_probability')
    return msh
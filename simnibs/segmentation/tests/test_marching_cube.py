# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 20:28:04 2020

@author: axthi
"""

import pytest
import numpy as np
from simnibs.segmentation.marching_cube import marching_cube

def two_spheres_image():
    """
    create image with 2 spheres (WM is 1, GM boundary is 0)

    Returns
    -------
    R_sp1 : 3d image (WM is 1, GM boundary is 0)
    vox2mm : affine trafo
    center_sp1 : center coords of sphere 1
    radius_sp1 : radius of sphere 1
    center_sp2 : ditto
    radius_sp2 : ditto
    """  

    imgsize=[101, 51, 71]
    voxsize=0.5
    
    center_sp1=[10, 0, 0]
    radius_sp1=10
    cortexth_sp1=4
    
    center_sp2=[-10, 0, 0]
    radius_sp2=5
    cortexth_sp2=3
    
    vox2mm=np.float32(([voxsize, 0, 0, -voxsize*(imgsize[0]-1)/2],
                       [0, voxsize, 0, -voxsize*(imgsize[1]-1)/2],
                       [0, 0, voxsize, -voxsize*(imgsize[2]-1)/2],
                       [0, 0, 0, 1]))
    
    G = np.meshgrid(np.arange(imgsize[0]), np.arange(imgsize[1]), 
                    np.arange(imgsize[2]), indexing='ij')
    
    xyz_mm=np.inner(vox2mm, np.vstack((G[0].flatten(),
                                       G[1].flatten(),
                                       G[2].flatten(),
                                       np.ones_like(G[0].flatten()))).T
                    )[0:3]
    
    R_sp1 = np.linalg.norm(xyz_mm.T-center_sp1, axis=1)
    R_sp1=(radius_sp1-R_sp1)/cortexth_sp1 + 0.5 # add 0.5 to center central surface at given radius
    R_sp1[R_sp1>1]=1
    R_sp1[R_sp1<0]=0
    
    R_sp2 = np.linalg.norm(xyz_mm.T-center_sp2, axis=1)
    R_sp2=(radius_sp2-R_sp2)/cortexth_sp2 + 0.5  
    R_sp2[R_sp2>1]=1
    R_sp2[R_sp2<0]=0
    
    R_sp1+=R_sp2
    R_sp1=R_sp1.reshape(imgsize)
       
    return R_sp1, vox2mm, center_sp1, radius_sp1, center_sp2, radius_sp2


class TestMarchingCube:
    def test_with_spheres(self): 
        testdata = two_spheres_image()
        img = testdata[0]
        vox2mm = testdata[1]
        center_sp1 = testdata[2] 
        radius_sp1 = testdata[3] 
        
        _, EC = marching_cube(img, vox2mm, level=0.5, step_size=2, only_largest_component=False)
        assert EC==4
        
        CS, EC = marching_cube(img, vox2mm, level=0.5, step_size=2, only_largest_component=True)
        assert EC==2
        assert np.allclose(np.linalg.norm(CS.nodes.node_coord-center_sp1, axis=1), radius_sp1, atol=2)

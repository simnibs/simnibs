# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 16:28:21 2020

@author: oupu
"""

from simnibs.segmentation.brain_surface import createCS, expandCS
from simnibs.mesh_tools.mesh_io import read_gifti_surface, read_curv, write_gifti_surface
import functools
import multiprocessing
import os


def expandCS_wrapper(actualsurf, surffolder):
    
    Pcentral = os.path.join(surffolder,actualsurf+'.central.gii')
    Ppial = os.path.join(surffolder,actualsurf+'.pial.gii')
    Pthickness = os.path.join(surffolder,actualsurf+'.thickness')
    
    m = read_gifti_surface(Pcentral)
    thickness = read_curv(Pthickness)
    m.nodes.node_coord = expandCS(m.nodes[:], m.elm[:,:3]-1, thickness/2, 
                                  ensure_distance=0.2, nsteps=5,
                                  deform="expand", smooth_mesh=True, 
                                  skip_lastsmooth=True, smooth_mm2move=True, 
                                  despike_nonmove=True, fix_faceflips=True,
                                  actualsurf=actualsurf)
    write_gifti_surface(m, Ppial)
    return Ppial


def run_cat_multiprocessing(Ymf, Yleft, Ymaskhemis,
                            vox2mm, surface_folder, fsavgDir, vdist,
                            voxsize_pbt, voxsize_refineCS, th_initial,
                            no_selfintersections, surf, pial = [], nprocesses = 0):

    Pcentral_all = []
    Pspherereg_all = []
    Pthick_all = []
    Ppial_all = []
    EC_all = []
    defect_size_all = []

    processes = len(surf)
    if nprocesses > 0:
        processes=min(nprocesses,processes)
        
    with multiprocessing.Pool(processes=processes) as pool:
        partial_create_cs = functools.partial(
            createCS, Ymf, Yleft, Ymaskhemis, vox2mm,
            surffolder=surface_folder, fsavgDir=fsavgDir, vdist=vdist,
            voxsize_pbt=voxsize_pbt, voxsize_refineCS=voxsize_refineCS,
            th_initial=th_initial, no_selfintersections=no_selfintersections)

        # call pool.map to run in parallel
        results = pool.map(partial_create_cs, surf)
        for r in results:
            Pcentral_all.append(r[0])
            Pspherereg_all.append(r[1])
            Pthick_all.append(r[2])
            EC_all.append(r[3])
            defect_size_all.append(r[4])

        if len(pial) > 0:
            assert all(elem in surf for elem in pial)
            partial_expand_cs = functools.partial(expandCS_wrapper, surffolder=surface_folder)

            # call pool.map to run in parallel
            results = pool.map(partial_expand_cs, pial)
            for r in results:
                Ppial_all.append(r[0])



if __name__ == '__main__':
    import argparse
    import nibabel as nib
   
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--Ymf_path", nargs="+", type=str)
    argument_parser.add_argument('--Yleft_path', nargs='+', type=str)
    argument_parser.add_argument('--Ymaskhemis_path', nargs='+', type=str)
    argument_parser.add_argument('--surface_folder', nargs='+', type=str)
    argument_parser.add_argument('--fsavgdir', nargs='+', type=str)
    argument_parser.add_argument('--vdist', nargs='+', type=float,
                                 default=[1.0, 0.75])
    argument_parser.add_argument('--voxsize_pbt', nargs='+', type=float,
                                 default=[0.5, 0.25])
    argument_parser.add_argument('--voxsizeCS', nargs='+', type=float,
                                 default=[0.75, 0.5])
    argument_parser.add_argument('--th_initial', nargs='+', type=float,
                                 default=[0.714])
    argument_parser.add_argument('--no_intersect', nargs='+', type=bool, default=[True])
    argument_parser.add_argument('--surf', nargs='+',
                                 default=['lh', 'rh', 'lc', 'rc'])
    argument_parser.add_argument('--pial', nargs='+',
                                 default=['lh', 'rh', 'lc', 'rc'])
    argument_parser.add_argument('--nprocesses', nargs='+', type=int,
                                 default=[0])
    parsed = argument_parser.parse_args()
    
    Ymf = nib.load(parsed.Ymf_path[0])
    Yleft = nib.load(parsed.Yleft_path[0])
    Yhemis = nib.load(parsed.Ymaskhemis_path[0])
    surf_folder = parsed.surface_folder[0]
    fsavgdir = parsed.fsavgdir[0]
    vdist = parsed.vdist
    voxsize_pbt = parsed.voxsize_pbt
    voxsize_refineCS = parsed.voxsizeCS
    th_initial = parsed.th_initial[0]
    no_selfintersections = parsed.no_intersect[0]
    surf = parsed.surf
    pial = parsed.pial
    nprocesses = parsed.nprocesses[0]
    
    run_cat_multiprocessing(Ymf.get_fdata(), Yleft.get_fdata(), Yhemis.get_fdata(),
                            Ymf.affine, surf_folder, fsavgdir, vdist,
                            voxsize_pbt, voxsize_refineCS, th_initial,
                            no_selfintersections, surf, pial, nprocesses)
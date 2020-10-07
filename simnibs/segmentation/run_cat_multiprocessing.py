# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 16:28:21 2020

@author: oupu
"""

from simnibs.segmentation.brain_surface import createCS
import functools
import multiprocessing
import argparse
import nibabel as nib


def run_cat_multiprocessing(Ymf, Yleft, Ymaskhemis, Ymaskparahipp,
                            vox2mm, surface_folder, fsavgDir, vdist,
                            voxsize_pbt, voxsize_refineCS, th_initial,
                            no_selfintersections, add_parahipp,
                            close_parahipp, surf):

    Pcentral_all = []
    Pspherereg_all = []
    Pthick_all = []
    EC_all = []
    defect_size_all = []

    with multiprocessing.Pool(processes=len(surf)) as pool:
        partial_create_cs = functools.partial(
            createCS, Ymf, Yleft, Ymaskhemis, Ymaskparahipp, vox2mm,
            surffolder=surface_folder, fsavgDir=fsavgDir, vdist=vdist,
            voxsize_pbt=voxsize_pbt, voxsize_refineCS=voxsize_refineCS,
            th_initial=th_initial, no_selfintersections=no_selfintersections,
            add_parahipp=add_parahipp, close_parahipp=close_parahipp)

        # call pool.map to run in parallel
        results = pool.map(partial_create_cs, surf)
        for r in results:
            Pcentral_all.append(r[0])
            Pspherereg_all.append(r[1])
            Pthick_all.append(r[2])
            EC_all.append(r[3])
            defect_size_all.append(r[4])


if __name__ == '__main__':
    argument_parser = argparse.ArgumentParser()
    argument_parser.add_argument("--Ymf_path", nargs="+", type=str)
    argument_parser.add_argument('--Yleft_path', nargs='+', type=str)
    argument_parser.add_argument('--Ymaskhemis_path', nargs='+', type=str)
    argument_parser.add_argument('--Ymaskparahipp_path', nargs='+', type=str)
    argument_parser.add_argument('--surface_folder', nargs='+', type=str)
    argument_parser.add_argument('--fsavgdir', nargs='+', type=str)
    argument_parser.add_argument('--vdist', nargs='+', type=float,
                                 default=[1.0, 0.75])
    argument_parser.add_argument('--voxsize_pbt', nargs='+', type=float,
                                 default=[0.5, 0.25])
    argument_parser.add_argument('--voxsizeCS', nargs='+', type=float,
                                 default=[0.75, 0.5])
    argument_parser.add_argument('--th_initial', nargs='+', type=float,
                                 default=[0.5])
    argument_parser.add_argument('--no_intersect', nargs='+', type=bool, default=[True])
    argument_parser.add_argument('--add_parahipp', nargs='+', type=bool, default=[False])
    argument_parser.add_argument('--close_parahipp', nargs='+', type=bool,
                                 default=[False])
    argument_parser.add_argument('--surf', nargs='+',
                                 default=['lh', 'rh', 'lc', 'rc'])

    parsed = argument_parser.parse_args()
    Ymf = nib.load(parsed.Ymf_path[0])
    Yleft = nib.load(parsed.Yleft_path[0])
    Yhemis = nib.load(parsed.Ymaskhemis_path[0])
    Yparahipp = nib.load(parsed.Ymaskparahipp_path[0])
    surf_folder = parsed.surface_folder[0]
    fsavgdir = parsed.fsavgdir[0]
    vdist = parsed.vdist
    voxsize_pbt = parsed.voxsize_pbt
    voxsize_refineCS = parsed.voxsizeCS
    th_initial = parsed.th_initial[0]
    no_selfintersections = parsed.no_intersect[0]
    add_parahipp = parsed.add_parahipp[0]
    close_parahipp = parsed.close_parahipp[0]
    surf = parsed.surf

    vox2mm = Ymf.affine
    run_cat_multiprocessing(Ymf.get_fdata(), Yleft.get_fdata(),
                            Yhemis.get_fdata(), Yparahipp.get_fdata(),
                            vox2mm, surf_folder, fsavgdir, vdist,
                            voxsize_pbt, voxsize_refineCS, th_initial,
                            no_selfintersections, add_parahipp,
                            close_parahipp, surf)

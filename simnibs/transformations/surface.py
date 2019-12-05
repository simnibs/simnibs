'''
    Surface-based transformations for SimNIBS
'''
'''
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.
    Copyright (C) 2020 Kristoffer H Madsen, Guilherme B Saturnino, Axel Thielscher

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
import os
from functools import partial

import numpy as np
import scipy.ndimage
import scipy.spatial
import nibabel as nib

from ..utils.simnibs_logger import logger
from .. import SIMNIBSDIR
from ..utils.file_finder import templates, SubjectFiles, get_atlas

__all__ = [
    'subject_atlas',
    'middle_gm_interpolation'
]

def get_surface_names_from_folder_structure(m2m_folder):
    # Subject name
    sub_files = SubjectFiles(subpath=m2m_folder)
    if not os.path.isdir(sub_files.subpath):
        raise IOError('The given m2m folder name does not correspond to a directory')

    def look_up(f):
        if os.path.isfile(f) or os.path.isdir(f):
            return f
        else:
            raise IOError('Could not find file or directory: {0}'.format(f))

    names = {}

    if sub_files.seg_type == 'headreco':
        names['surf_dir'] = look_up(sub_files.surf_dir)
        names['lh_midgm'] = look_up(sub_files.lh_midgm)
        names['rh_midgm'] = look_up(sub_files.rh_midgm)
        names['lh_reg'] = look_up(sub_files.lh_reg)
        names['rh_reg'] = look_up(sub_files.rh_reg)
        names['lh_sphere_ref'] = look_up(templates.cat_lh_sphere_ref)
        names['rh_sphere_ref'] = look_up(templates.cat_rh_sphere_ref)
        names['ref_fs'] = look_up(sub_files.ref_fs)
        names['lh_cortex_ref'] = look_up(templates.cat_lh_cortex_ref)
        names['rh_cortex_ref'] = look_up(templates.cat_rh_cortex_ref)

    elif sub_files.seg_type == 'mri2mesh':
        names['subj_id'] = sub_files.subid
        names['surf_dir'] = look_up(sub_files.surf_dir)
        names['lh_gm'] = look_up(sub_files.lh_gm)
        names['lh_wm'] = look_up(sub_files.lh_wm)
        names['rh_gm'] = look_up(sub_files.rh_gm)
        names['rh_wm'] = look_up(sub_files.rh_wm)
        names['lh_reg'] = look_up(sub_files.lh_reg)
        names['rh_reg'] = look_up(sub_files.rh_reg)
        names['lh_sphere_ref'] = look_up(templates.fs_lh_sphere_ref)
        names['rh_sphere_ref'] = look_up(templates.fs_rh_sphere_ref)
        names['ref_fs'] = look_up(sub_files.ref_fs)
        names['lh_cortex_ref'] = look_up(templates.fs_lh_cortex_ref)
        names['rh_cortex_ref'] = look_up(templates.fs_rh_cortex_ref)

    else:
        raise IOError('Could not find surface files in m2m folder. SPM-only segmentation?')

    return names, sub_files.seg_type


def _surf2surf(field, in_surf, out_surf, kdtree=None):
    ''' Interpolates the field defined in in_vertices to the field defined in
    out_vertices using nearest neighbour '''
    if kdtree is None:
        # Normalize the radius of the input sphere
        in_v = np.copy(in_surf.nodes.node_coord)
        in_v /= np.average(np.linalg.norm(in_v, axis=1))
        kdtree = scipy.spatial.cKDTree(in_v)
    out_v = np.copy(out_surf.nodes.node_coord)
    # Normalize the radius of the output sphere
    out_v /= np.average(np.linalg.norm(out_v, axis=1))
    _, closest = kdtree.query(out_v)
    return field[closest], kdtree


def middle_gm_interpolation(mesh_fn, m2m_folder, out_folder, out_fsaverage=None,
                            depth=0.5, quantities=['norm', 'normal', 'tangent','angle'],
                            fields=None, open_in_gmsh=False):
    ''' Interpolates the vector fieds in the middle gray matter surface

    Parameters
    -----------
    mesh_fn: str
        String with file name to mesh
    m2m_folder: str
        Path to the m2m_{subject_id} folder, generated during the segmantation
    out_folder: str
        Name of output folder. Output files will be written to this folder
    out_fsaverage: str (optional)
        Name of output folder for transformed files. In not set, fsaverage transormation
        will not be performed
    depth: float (optional)
        The distance bewteen grey and white matter where the
        interpolation should be done. p = depth * wm + (1 - depth) * gm (default: .5)
    quantities: list with the elements {norm, normal, tangent, angle}
        Quantites to be calculated from vector field
    fields: list of strings (optional)
        Fields to be transformed. Default: all fields
    open_in_gmsh: bool
        If true, opens a Gmsh window with the interpolated fields
    '''
    from . import mesh_io
    m2m_folder = os.path.abspath(os.path.normpath(m2m_folder))
    names, segtype = get_surface_names_from_folder_structure(m2m_folder)
    if depth < 0. or depth > 1.:
        raise ValueError('Invalid depth value. Should be between 0 and 1')

    if any([q not in ['norm', 'normal', 'tangent', 'angle'] for q in quantities]):
        raise ValueError('Invalid quanty in {0}'.format(quantities))

    def calc_quantities(nd, quantities):
        d = dict.fromkeys(quantities)
        for q in quantities:
            if q == 'norm':
                d[q] = nd.norm()
            elif q == 'normal':
                d[q] = nd.normal()
                d[q].value *= -1
            elif q == 'tangent':
                d[q] = nd.tangent()
            elif q == 'angle':
                d[q] = nd.angle()
            else:
                raise ValueError('Invalid quantity: {0}'.format(q))
        return d

    m = mesh_io.read_msh(mesh_fn)
    subdir, sim_name = os.path.split(mesh_fn)
    sim_name = '.' + os.path.splitext(sim_name)[0]
    # Crio out GM
    m = m.crop_mesh(2)
    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)
    out_folder = os.path.abspath(os.path.normpath(out_folder))

    if out_fsaverage is not None and not os.path.isdir(out_fsaverage):
        os.mkdir(out_fsaverage)
    if out_fsaverage is not None:
        out_fsaverage = os.path.abspath(os.path.normpath(out_fsaverage))

    middle_surf = {}
    reg_surf = {}
    ref_surf = {}
    avg_surf = {}
    # Load and write furfaces
    if segtype == 'mri2mesh':
        for hemi in ['lh', 'rh']:
            wm_surface = mesh_io.read_freesurfer_surface(names[hemi + '_wm'])
            gm_surface = mesh_io.read_freesurfer_surface(names[hemi + '_gm'])
            middle_surf[hemi] = mesh_io._middle_surface(wm_surface, gm_surface, depth)
            mesh_io.write_freesurfer_surface(
                middle_surf[hemi],
                os.path.join(out_folder, hemi + '.central'),
                names['ref_fs'])
    elif segtype == 'headreco':
        for hemi in ['lh', 'rh']:
            middle_surf[hemi] = mesh_io.read_gifti_surface(names[hemi + '_midgm'])
            mesh_io.write_freesurfer_surface(
                middle_surf[hemi],
                os.path.join(out_folder, hemi + '.central'),
                names['ref_fs'])
    # Load average space things
    if out_fsaverage:
        for hemi in ['lh', 'rh']:
            if segtype == 'headreco':
                reg_surf[hemi] = \
                    mesh_io.read_gifti_surface(names[hemi+'_reg'])
                ref_surf[hemi] = \
                    mesh_io.read_gifti_surface(names[hemi+'_sphere_ref'])
                avg_surf[hemi] = \
                    mesh_io.read_gifti_surface(names[hemi+'_cortex_ref'])
            if segtype == 'mri2mesh':
                reg_surf[hemi] = \
                    mesh_io.read_freesurfer_surface(names[hemi+'_reg'])
                ref_surf[hemi] = \
                    mesh_io.read_freesurfer_surface(names[hemi+'_sphere_ref'])
                avg_surf[hemi] = \
                    mesh_io.read_freesurfer_surface(names[hemi+'_cortex_ref'])

    names_subj = []
    names_fsavg = []
    kdtree = {'lh': None, 'rh': None}
    h = []
    for name, data in m.field.items():
        for hemi in ['lh', 'rh']:
            if fields is None or name in fields:
                # Interpolate to middle gm
                data = data.as_nodedata()
                interpolated = data.interpolate_to_surface(middle_surf[hemi])
                # For vector quantities, calculate quantities (normal, norm, ...)
                if data.nr_comp == 3:
                    q = calc_quantities(interpolated, quantities)
                    for q_name, q_data in q.items():
                        out_subj = os.path.join(
                            out_folder, hemi + sim_name + '.central.' + name + '.' + q_name)
                        mesh_io.write_curv(
                            out_subj,
                            q_data.value,
                            middle_surf[hemi].elm.nr)
                        names_subj.append(out_subj)
                        middle_surf[hemi].add_node_field(q_data, name + '_' + q_name)
                        h.append(hemi)
                        # Interpolate to fsavg
                        if out_fsaverage is not None:
                            q_transformed, kdtree[hemi] = _surf2surf(
                                q_data.value,
                                reg_surf[hemi],
                                ref_surf[hemi],
                                kdtree[hemi])
                            out_avg = os.path.join(
                                          out_fsaverage,
                                          hemi + sim_name + '.fsavg.'
                                          + name + '.' + q_name)
                            mesh_io.write_curv(
                                out_avg,
                                q_transformed,
                                ref_surf[hemi].elm.nr)
                            avg_surf[hemi].add_node_field(q_transformed, name + '_' + q_name)
                            names_fsavg.append(out_avg)

                # For scalar quantities
                elif data.nr_comp == 1:
                    field_name = name[-1]
                    q_name = name[:-1]
                    if field_name in m.field.keys() and q_name in quantities:
                        # If we have an equivalent quantity being calculated, skip
                        pass
                    else:
                        out_subj = os.path.join(
                            out_folder, hemi + sim_name + '.central.' + name)
                        mesh_io.write_curv(
                            out_subj,
                            interpolated.value.squeeze(),
                            middle_surf[hemi].elm.nr)
                        names_subj.append(out_subj)
                        h.append(hemi)
                        middle_surf[hemi].add_node_field(interpolated, name)
                        if out_fsaverage is not None:
                            f_transformed, kdtree[hemi] = _surf2surf(
                                interpolated.value.squeeze(),
                                reg_surf[hemi],
                                ref_surf[hemi],
                                kdtree[hemi])
                            out_avg = os.path.join(
                                out_fsaverage,
                                hemi + sim_name + '.fsavg.' + name)
                            mesh_io.write_curv(
                                out_avg, f_transformed, ref_surf[hemi].elm.nr)
                            names_fsavg.append(out_avg)
                            avg_surf[hemi].add_node_field(f_transformed, name)


    # Join surfaces, fields and open in gmsh
    def join_and_write(surfs, fn_out, open_in_gmsh):
        mesh = surfs['lh'].join_mesh(surfs['rh'])
        mesh.elm.tag1 = 1002 * np.ones(mesh.elm.nr, dtype=int)
        mesh.elm.tag2 = 1002 * np.ones(mesh.elm.nr, dtype=int)
        mesh.nodedata = []
        mesh.elmdata = []
        for k in surfs['lh'].field.keys():
            mesh.add_node_field(
                np.append(surfs['lh'].field[k].value,
                          surfs['rh'].field[k].value),
                k)
        v = mesh.view(visible_fields=list(surfs['lh'].field.keys())[0])
        v.write_opt(fn_out)
        mesh_io.write_msh(mesh, fn_out)
        if open_in_gmsh:
            mesh_io.open_in_gmsh(fn_out, True)

    join_and_write(
        middle_surf,
        os.path.join(out_folder, sim_name[1:] + '_central.msh'),
        open_in_gmsh)
    if out_fsaverage:
        join_and_write(
            avg_surf,
            os.path.join(out_fsaverage, sim_name[1:] + '_fsavg.msh'),
            open_in_gmsh)


def subject_atlas(atlas_name, m2m_dir, hemi='both'):
    ''' Loads a brain atlas based of the FreeSurfer fsaverage template

    Parameters
    -----------
    atlas_name: 'a2009s', 'DK40' or 'HCP_MMP1'
            Name of atlas to load

            'a2009s': Destrieux atlas (FreeSurfer v4.5, aparc.a2009s)
            Cite: Destrieux, C. Fischl, B. Dale, A., Halgren, E. A sulcal
            depth-based anatomical parcellation of the cerebral cortex.
            Human Brain Mapping (HBM) Congress 2009, Poster #541

            'DK40': Desikan-Killiany atlas (FreeSurfer, aparc.a2005s)
            Cite: Desikan RS, S�gonne F, Fischl B, Quinn BT, Dickerson BC,
            Blacker D, Buckner RL, Dale AM, Maguire RP, Hyman BT, Albert MS,
            Killiany RJ. An automated labeling system for subdividing the
            human cerebral cortex on MRI scans into gyral based regions of
            interest. Neuroimage. 2006 Jul 1;31(3):968-80.

            'HCP_MMP1': Human Connectome Project (HCP) Multi-Modal Parcellation
            Cite: Glasser MF, Coalson TS, Robinson EC, et al. A multi-modal
            parcellation of human cerebral cortex. Nature. 2016;536(7615):171-178.

    m2m_folder: str
        Path to the m2m_{subject_id} folder, generated during the segmantation
 
    hemi (optional): 'lh', 'rh' or 'both'
        Hemisphere to use. In the case of 'both', will assume that left hemisphere
        nodes comes before right hemisphere nodes

    Returns
    ---------
    atlas: dict
        Dictionary where atlas['region'] = roi
    '''
    from .mesh_io import read_gifti_surface, read_freesurfer_surface
    if atlas_name not in ['a2009s', 'DK40', 'HCP_MMP1']:
        raise ValueError('Invalid atlas name')

    subject_files = SubjectFiles(subpath=m2m_dir)

    if hemi in ['lh', 'rh']:
        fn_atlas = os.path.join(
            templates.cat_atlases_surfaces,
            f'{hemi}.aparc_{atlas_name}.freesurfer.annot'
        )
        labels, _ , names = nib.freesurfer.io.read_annot(fn_atlas)
        if subject_files.seg_type == 'headreco':
            read_fun = read_gifti_surface
        elif subject_files.seg_type == 'mri2mesh':
            read_fun = read_freesurfer_surface

        if hemi == 'lh':
            labels_sub, _ = _surf2surf(
                labels,
                read_gifti_surface(templates.cat_lh_sphere_ref),
                read_fun(subject_files.lh_reg)
                ) 
        if hemi == 'rh':
            labels_sub, _ = _surf2surf(
                labels,
                read_gifti_surface(templates.cat_rh_sphere_ref),
                read_fun(subject_files.rh_reg)
                ) 
        atlas = {}
        for l, name in enumerate(names):
            atlas[name.decode()] = labels_sub == l

        return atlas

    # If both hemispheres
    elif hemi == 'both':
        atlas_lh = subject_atlas(atlas_name, m2m_dir, 'lh')
        atlas_rh = subject_atlas(atlas_name, m2m_dir, 'rh')
        atlas = {}
        pad_rh = np.zeros_like(list(atlas_rh.values())[0])
        pad_lh = np.zeros_like(list(atlas_lh.values())[0])
        for name, mask in atlas_lh.items():
            atlas[f'lh.{name}'] = np.append(mask, pad_rh)  # pad after
        for name, mask in atlas_rh.items():
            atlas[f'rh.{name}'] = np.append(pad_lh, mask)  # pad after

        return atlas
    else:
        raise ValueError('Invalid hemisphere name')


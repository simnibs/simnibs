# -*- coding: utf-8 -*-\
'''
    IO functions for Gmsh .msh files
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2013-2018 Andre Antunes, Guilherme B Saturnino, Kristoffer H Madsen

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
from __future__ import division
from __future__ import print_function
import os
import struct
import copy
import datetime
import warnings
import gc
import hashlib
import tempfile
import subprocess
import threading

from functools import partial
import numpy as np
import scipy.spatial
import scipy.ndimage
import scipy.sparse
import scipy.interpolate
import nibabel
import h5py

from .transformations import nifti_transform
from . import gmsh_view
from ..utils.file_finder import path2bin, templates
import simnibs.cython_code.cython_msh as cython_msh


__all__ = [
    'read_msh',
    'write_msh',
    'read_freesurfer_surface',
    'write_freesurfer_surface',
    'read_gifti_surface',
    'read_curv',
    'write_curv',
    'read_stl',
    'write_geo_spheres',
    'write_geo_text',
    'Msh',
    'Nodes',
    'Elements',
    'ElementData',
    'NodeData'
]
# =============================================================================
# CLASSES
# =============================================================================

class Nodes:
    """class to handle the node information:

    Parameters
    -----------------------
    node_coord (optional): (Nx3) ndarray
        Coordinates of the nodes

    Attributes
    ----------------------
    node_coord: (Nx3) ndarray
        Coordinates of the nodes

    nr: property
        Number of nodes

    Examples
    -----------------------------------
     >>> nodes = Nodes(np.array([[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
     array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
     >>> nodes.node_number
     array([1, 2, 3])
     >>> nodes[1]
     array([1, 0, 0])
     >>> nodes[array(True, False, True)]
     array([[1, 0, 0], [0, 0, 1]])
    """

    def __init__(self, node_coord=None):
        # gmsh fields
        self.node_coord = np.array([], dtype='float64')
        if node_coord is not None:
            self.node_coord = node_coord

    @property
    def nr(self):
        ''' Number of nodes '''
        return self.node_coord.shape[0]

    @property
    def node_number(self):
        ''' Node numbers (1, ..., nr) '''
        return np.array(range(1, self.nr + 1), dtype='int32')

    def find_closest_node(self, querry_points, return_index=False):
        """ Finds the closest node to each point in p

        Parameters
        --------------------------------
        querry_points: (Nx3) ndarray
            List of points (x,y,z) to which the closes node in the mesh should be found

        return_index: (optional) bool
        Whether to return the index of the closes nodes, default: False

        Returns
        -------------------------------
        coords: Nx3 array of floats
            coordinates of the closest points

        indexes: Nx1 array of ints
            Indices of the nodes in the mesh

        --------------------------------------
        The indices are in the mesh listing, that starts at one!
       """
        if len(self.node_coord) == 0:
            raise ValueError('Mesh has no nodes defined')

        kd_tree = scipy.spatial.cKDTree(self.node_coord)
        _, indexes = kd_tree.query(querry_points)
        coords = self.node_coord[indexes, :]

        if return_index:
            return (coords, indexes + 1)

        else:
            return coords

    def __eq__(self, other):
        try:
            return self.__dict__ == other.__dict__
        except AttributeError:
            return False

    def __getitem__(self, index):
        return _getitem_one_indexed(self.node_coord, index)

    def __str__(self):
        return str(self.node_coord)


class Elements:
    """ Mesh elements.

    Can only handle triangles and tetrahedra!

    Parameters
    --------------------------
    triangles (optional): (Nx3) ndarray
        List of nodes composing each triangle
    tetrahedra(optional): (Nx3) ndarray
        List of nodes composing each tetrahedra


    Attributes
    ----------------------------------
    elm_number: (Nx1) ndarray
          element ID (u from 1 till nr)
    elm_type: (Nx1) ndarray
        elm-type (2=triangle, 4=tetrahedron, etc)
    tag1: (Nx1) ndarray
        first tag for each element
    tag2: (Nx1) ndarray
        second tag for each elenent
    node_number_list: (Nx4) ndarray
        4xnumber_of_element matrix of the nodes that constitute the element.
        For the triangles, the fourth element = -1
    nr: int
        Number or elemets


    Notes
    -------------------------
    Node and element count starts at 1!

    """

    def __init__(self, triangles=None, tetrahedra=None):
        # gmsh fields
        self.elm_type = np.array([], 'int8')
        self.tag1 = np.array([], dtype='int16')
        self.tag2 = np.array([], dtype='int16')
        self.node_number_list = np.array([], dtype='int32')

        if triangles is not None:
            assert triangles.shape[1] == 3
            assert np.all(triangles > 0), "Node count should start at 1"
            self.node_number_list = np.zeros(
                (triangles.shape[0], 4), dtype='int32')
            self.node_number_list[:, :3] = triangles.astype('int32')
            self.node_number_list[:, 3] = -1
            self.elm_type = np.ones((self.nr,), dtype='int32') * 2

        if tetrahedra is not None:
            assert tetrahedra.shape[1] == 4
            assert np.all(tetrahedra > 0), "Node count should start at 1"
            if len(self.node_number_list) == 0:
                self.node_number_list = tetrahedra.astype('int32')
                self.elm_type = np.ones((self.nr,), dtype='int32') * 4
            else:
                self.node_number_list = np.vstack(
                    (self.node_number_list, tetrahedra.astype('int32')))
                self.elm_type = np.append(
                    self.elm_type, np.ones((self.nr,), dtype='int32') * 4)

        if len(self.node_number_list) > 0:
            self.tag1 = np.ones((self.nr,), dtype='int32')
            self.tag2 = np.ones((self.nr,), dtype='int32')

    @property
    def nr(self):
        ''' Number of elements '''
        return self.node_number_list.shape[0]

    @property
    def triangles(self):
        ''' Triangle element numbers '''
        return self.elm_number[self.elm_type == 2]

    @property
    def tetrahedra(self):
        ''' Tetrahedra element numbers '''
        return self.elm_number[self.elm_type == 4]

    @property
    def elm_number(self):
        ''' Element numbers (1, ..., nr) '''
        return np.arange(1, self.nr + 1, dtype='int32')

    def find_all_elements_with_node(self, node_nr):
        """ Finds all elements that have a given node

        Parameters
        -----------------
        node_nr: int
            number of node

        Returns
        ---------------
        elm_nr: np.ndarray
            array with indices of element numbers

        """
        elm_with_node = np.any(
            np.isin(self.node_number_list, node_nr),
            axis=1)
        return self.elm_number[elm_with_node]

    def find_neighbouring_nodes(self, node_nr):
        """ Finds the nodes that share an element with the specified node

        Parameters
        -----------------------------
        node_nr: int
            number of query node (mesh listing)
        Returns
        ------------------------------
        all_neighbours: np.ndarray
            list of all nodes what share an element with the node

        Example
        ----------------------------
        >>> elm = msh.Elements()
        >>> elm.node_number_list = np.array([[1, 2, 3, 4], [3, 2, 5, 7], [1, 8, 2, 6]])
        >>> elm.find_neighbouring_nodes(1))
        array([2, 3, 4, 8, 6])

        """
        elm_with_node = np.array(
            sum(self.node_number_list[:, i] == node_nr for i in range(4)),
            dtype=bool)
        all_neighbours = self.node_number_list[elm_with_node, :].reshape(-1)
        all_neighbours = np.unique(all_neighbours)
        all_neighbours = all_neighbours[all_neighbours != node_nr]
        all_neighbours = all_neighbours[all_neighbours != 0]
        return np.array(all_neighbours, dtype=int)

    def get_faces(self, tetrahedra_indexes=None):
        ''' Creates a list of nodes in each face and a list of faces in each tetrahedra

        Parameters
        ----------------
        tetrahedra_indexes: np.ndarray
            Indices of the tetrehedra where the faces are to be determined (default: all
            tetrahedra)

        Returns
        -------------
        faces: np.ndarray
            List of nodes in faces, in arbitrary order
        th_faces: np.ndarray
            List of faces in each tetrahedra, starts at 0, order=((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3))
        face_adjacency_list: np.ndarray
            List of tetrahedron adjacent to each face, filled with -1 if a face is in a
            single tetrahedron. Not in the normal element ordering, but only in the order
            the tetrahedra are presented
        '''
        if tetrahedra_indexes is None:
            tetrahedra_indexes = self.tetrahedra
        th = self[tetrahedra_indexes]
        faces = th[:, [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]]
        faces = faces.reshape(-1, 3)

        # Try to hash a few times
        warn = True
        for i in range(3):
            hash_array = _hash_rows(faces, mult=1000003 + i)
            unique, idx, inv, count = np.unique(hash_array, return_index=True,
                                                return_inverse=True, return_counts=True)
            if np.all(count <= 2):
                warn = False
                break

        if warn:
            warnings.warn('Invalid Mesh: Found a face with more than 2 adjacent'
                          ' tetrahedra!')

        faces = faces[idx]
        face_adjacency_list = -np.ones((len(unique), 2), dtype=int)
        face_adjacency_list[:, 0] = idx // 4

        # Remove the faces already seen from consideration
        # Second round in order to make adjacency list
        # create a new array with a mask in the elements already seen
        mask = unique[-1] + 1
        hash_array_masked = np.copy(hash_array)
        hash_array_masked[idx] = mask
        # make another array, where we delete the elements we have already seen
        hash_array_reduced = np.delete(hash_array, idx)
        # Finds where each element of the second array is in the first array
        # (https://stackoverflow.com/a/8251668)
        hash_array_masked_sort = hash_array_masked.argsort()
        hash_array_repeated_pos = hash_array_masked_sort[
            np.searchsorted(hash_array_masked[hash_array_masked_sort], hash_array_reduced)]
        # Now find the index of the face corresponding to each element in the
        # hash_array_reduced
        faces_repeated = np.searchsorted(unique, hash_array_reduced)
        # Finally, fill out the second column in the adjacency list
        face_adjacency_list[faces_repeated, 1] = hash_array_repeated_pos // 4

        return faces, inv.reshape(-1, 4), face_adjacency_list

    def get_outside_faces(self, tetrahedra_indexes=None):
        ''' Creates a list of nodes in each face that are in the outer volume

        Parameters
        ----------------
        tetrahedra_indexes: np.ndarray
            Indices of the tetrehedra where the outer volume is to be determined (default: all
            tetrahedra)
        Returns
        -------------
        faces: np.ndarray
            List of nodes in faces, in arbitrary order
        '''
        if tetrahedra_indexes is None:
            tetrahedra_indexes = self.tetrahedra
        th = self[tetrahedra_indexes]
        faces = th[:, [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]]
        faces = faces.reshape(-1, 3)
        warn = True
        for i in range(3):
            hash_array = _hash_rows(faces, mult=1000003 + i)
            unique, idx, inv, count = np.unique(hash_array, return_index=True,
                                                return_inverse=True, return_counts=True)
            if np.all(count <= 2):
                warn = False
                break

        if warn:
            warnings.warn('Invalid Mesh: Found a face with more than 2 adjacent'
                          ' tetrahedra!')

        outside_faces = faces[idx[count == 1]]
        return outside_faces

    def nodes_with_tag(self, tags):
        ''' Gets all nodes indexes that are part of at least one element with the given
        tags

        Parameters
        -----------
        tags: list
            Integer tags to search

        Returns
        -------------
        nodes: ndarray of integer
            Indexes of nodes with given tag
        '''
        nodes = np.unique(self[np.isin(self.tag1, tags)].reshape(-1))
        nodes = nodes[nodes > 0]
        return nodes


    def __getitem__(self, index):
        return _getitem_one_indexed(self.node_number_list, index)

    def __eq__(self, other):
        try:
            return self.__dict__ == other.__dict__
        except AttributeError:
            return False

    def __str__(self):
        s = ''
        s += 'Nr elements: {0}\n'.format(self.nr)
        s += 'elm types: {0}\n'.format(self.elm_type)
        s += 'tags: {0}\n'.format(self.tag1)
        s += 'node list: {0}'.format(self.node_number_list)
        return s

class Msh:
    """class to handle the meshes.
    Gatters Nodes, Elements and Data

    Parameters
    -------------------------
    nodes: (optional) simnibs.msh.Nodes
        Nodes structure

    elements: (optional) simnibs.msh.Elements()
        Elements structure

    fn: str (optional)
        Name of ".msh" file to be read.
        Overides nodes and elements

    Attributes
    -------------------------
    nodes: simnibs.msh.Nodes
        a Nodes field
    elm: simnibs.msh.Elements
        A Elements field
    nodedata: simnibs.msh.NodeData
        list of NodeData filds
    elmdata: simnibs.msh.ElementData
        list of ElementData fields
    fn: str
        name of file
    binary: bool
        wheather or not the mesh was in binary format
   """

    def __init__(self, nodes=None, elements=None, fn=None):
        self.nodes = Nodes()
        self.elm = Elements()
        self.nodedata = []
        self.elmdata = []
        self.fn = ''  # file name to save msh

        if nodes is not None:
            self.nodes = nodes
        if elements is not None:
            self.elm = elements

        if fn is not None:
            self = read_msh(fn, m=self)

    @property
    def field(self):
        '''Dictionary of fields indexed by their name'''
        return dict(
            [(data.field_name, data) for data in self.nodedata + self.elmdata])

    def write(self, out_fn):
        ''' Writes out the mesh as a ".msh" file

        Parameters
        ---------------
        out_fn: str
            Name of output file
        '''
        write_msh(self, out_fn)

    def crop_mesh(self, tags=None, elm_type=None, nodes=None, elements=None):
        """ Crops the specified tags from the mesh
        Generates a new mesh, with only the specified tags
        The nodes are also reordered

        Parameters
        ---------------------
        tags:(optinal) int or list
            list of tags to be cropped, default: all

        elm_type: (optional) list of int
            list of element types to be croped (2 for triangles, 4 for tetrahedra), default: all

        nodes: (optional) list of ints
            List of nodes to be cropped, returns the minimal mesh containing elements
            wih at least one of the given nodes

        elements: (optional) list of ints
            List of elements to be cropped

        Returns
        ---------------------
        simnibs.msh.Msh
            Mesh with only the specified tags/elements/nodes

        Raises
        -----------------------
            ValueError, if the tag and elm_type combination is not foud

        Notes
        -----------
        If more than one (tags, elm_type, nodes, elements) is selected, they are joined by an OR
        operation
        """
        if tags is None and elm_type is None and nodes is None and elements is None:
            raise ValueError("At least one type of crop must be specified")

        elm_keep = np.zeros((self.elm.nr, ), dtype=bool)

        if tags is not None:
            elm_keep += np.in1d(self.elm.tag1, tags)

        if elm_type is not None:
            elm_keep += np.in1d(self.elm.elm_type, elm_type)

        if nodes is not None:
            elm_keep += np.any(np.in1d(self.elm.node_number_list, nodes).reshape(-1, 4), axis=1)

        if elements is not None:
            elm_keep += np.in1d(self.elm.elm_number, elements)

        if not np.any(elm_keep):
            raise ValueError("Could not find any element to crop!")

        idx = np.where(elm_keep)[0]
        nr_elements = len(idx)
        unique_nodes = np.unique(self.elm.node_number_list[idx, :].reshape(-1))
        unique_nodes = unique_nodes[unique_nodes != -1]
        if unique_nodes[0] == 0:
            unique_nodes = np.delete(unique_nodes, 0)
        nr_unique = np.size(unique_nodes)

        # creates a dictionary
        nodes_dict = np.zeros(self.nodes.nr + 1, dtype='int')
        nodes_dict[unique_nodes] = np.arange(1, 1 + nr_unique, dtype='int32')

        # Gets the new node numbers
        node_number_list = nodes_dict[self.elm.node_number_list[idx, :]]

        # and the positions in appropriate order
        node_coord = self.nodes.node_coord[unique_nodes - 1]

        # gerenates new mesh
        cropped = Msh()

        cropped.elm.tag1 = np.copy(self.elm.tag1[idx])
        cropped.elm.tag2 = np.copy(self.elm.tag2[idx])
        cropped.elm.elm_type = np.copy(self.elm.elm_type[idx])
        cropped.elm.node_number_list = np.copy(node_number_list)
        cropped.elm.node_number_list[cropped.elm.elm_type == 2, 3] = -1

        cropped.nodes.node_coord = np.copy(node_coord)

        cropped.nodedata = copy.deepcopy(self.nodedata)

        for nd in cropped.nodedata:
            nd.mesh = cropped
            if nd.nr_comp == 1:
                nd.value = np.copy(nd.value[unique_nodes - 1])
            else:
                nd.value = np.copy(nd.value[unique_nodes - 1, :])

        for ed in self.elmdata:
            cropped.elmdata.append(
                ElementData(ed.value[idx],
                            ed.field_name,
                            mesh=cropped))

        return cropped

    def join_mesh(self, other):
        ''' Join the current mesh with another

        Parameters
        -----------
        other: simnibs.msh.Msh
            Mesh to be joined

        Returns
        --------
        joined: simnibs.msh.Msh
            Mesh with joined nodes and elements
        '''
        joined = copy.deepcopy(self)
        joined.elmdata = []
        joined.nodedata = []
        other = copy.deepcopy(other)

        joined.nodes.node_coord = np.vstack([joined.nodes.node_coord,
                                             other.nodes.node_coord])
        other_node_number_list = other.elm.node_number_list + self.nodes.nr
        joined.elm.node_number_list = np.vstack([joined.elm.node_number_list,
                                                 other_node_number_list])
        joined.elm.tag1 = np.hstack([joined.elm.tag1, other.elm.tag1])
        joined.elm.tag2 = np.hstack([joined.elm.tag2, other.elm.tag2])
        joined.elm.elm_type = np.hstack([joined.elm.elm_type, other.elm.elm_type])
        triangles = np.where(joined.elm.elm_type == 2)[0]
        tetrahedra = np.where(joined.elm.elm_type != 2)[0]
        new_elm_order = np.hstack([triangles, tetrahedra])
        joined.elm.node_number_list = joined.elm.node_number_list[new_elm_order]
        joined.elm.tag1 = joined.elm.tag1[new_elm_order]
        joined.elm.tag2 = joined.elm.tag2[new_elm_order]
        joined.elm.elm_type = joined.elm.elm_type[new_elm_order]


        for nd in self.nodedata:
            assert len(nd.value) == self.nodes.nr
            pad_length = [(0, other.nodes.nr)] + [(0, 0)] * (nd.value.ndim - 1)
            new_values = np.pad(nd.value.astype(float), pad_length, 'constant', constant_values=np.nan)
            joined.nodedata.append(NodeData(new_values, nd.field_name))

        for ed in self.elmdata:
            assert len(ed.value) == self.elm.nr
            pad_length = [(0, other.elm.nr)] + [(0, 0)] * (ed.value.ndim - 1)
            new_values = np.pad(ed.value.astype(float), pad_length, 'constant', constant_values=np.nan)
            joined.elmdata.append(ElementData(new_values[new_elm_order], ed.field_name))

        for nd in other.nodedata:
            assert len(nd.value) == other.nodes.nr
            pad_length = [(self.nodes.nr, 0)] + [(0, 0)] * (nd.value.ndim - 1)
            new_values = np.pad(nd.value.astype(float), pad_length, 'constant', constant_values=np.nan)
            joined.nodedata.append(NodeData(new_values, nd.field_name))

        for ed in other.elmdata:
            assert len(ed.value) == other.elm.nr
            pad_length = [(self.elm.nr, 0)] + [(0, 0)] * (ed.value.ndim - 1)
            new_values = np.pad(ed.value.astype(float), pad_length, 'constant', constant_values=np.nan)
            joined.elmdata.append(ElementData(new_values[new_elm_order], ed.field_name))

        return joined

    def remove_from_mesh(self, tags):
        """ Removes the specified tags from the mesh
        Generates a new mesh, with the specified tags removed
        The nodes are also reordered

        Parameters
        ---------------------
        tags: int or list
            list of tags to be removed

        Returns
        ---------------------
        simnibs.msh.Msh
            Mesh without the specified tags

        """
        remove = np.in1d(self.elm.tag1, tags)
        keep = self.elm.elm_number[~remove]
        mesh = self.crop_mesh(elements=keep)
        return mesh

    def elements_baricenters(self):
        """ Calculates the baricenter of the elements

        Returns
        ------------
        baricenters: ElementData
            ElementData with the baricentes of the elements

        """
        bar = ElementData(
            np.zeros((self.elm.nr, 3), dtype=float),
            'baricenter', mesh=self)
        bar.mesh = self
        th_indexes = self.elm.tetrahedra
        tr_indexes = self.elm.triangles

        if len(th_indexes) > 0:
            bar[th_indexes] = np.average(
                self.nodes[self.elm[th_indexes]],
                axis=1)

        if len(tr_indexes) > 0:
            bar[tr_indexes] = np.average(
                self.nodes[self.elm[tr_indexes][:, :3]],
                axis=1)

        return bar

    def elements_volumes_and_areas(self):
        """ Calculates the volumes of tetrahedra and areas of triangles

        Returns
        ----------
        v: simnibs.Msh.ElementData
            Volume/areas of tetrahedra/triangles

        Note
        ------
            In the mesh's unit (normally mm)
        """
        vol = ElementData(np.zeros(self.elm.nr, dtype=float),
                          'volumes_and_areas')
        th_indexes = self.elm.tetrahedra
        tr_indexes = self.elm.triangles

        node_tr = self.nodes[self.elm[tr_indexes, :3]]
        sideA = node_tr[:, 1] - node_tr[:, 0]

        sideB = node_tr[:, 2] - node_tr[:, 0]

        n = np.cross(sideA, sideB)

        vol[tr_indexes] = np.linalg.norm(n, axis=1) * 0.5

        node_th = self.nodes[self.elm[th_indexes]]
        M = node_th[:, 1:] - node_th[:, 0, None]
        vol[th_indexes] = np.abs(np.linalg.det(M)) / 6.

        return vol

    def find_closest_element(self, querry_points, return_index=False,
                             elements_of_interest=None,
                             k=1, return_distance=False):
        """ Finds the closest node to each point in p

        Parameters
        --------------------------------
        querry_points: (Nx3) ndarray
            List of points (x,y,z) to which the closest element in the mesh should be found

        return_index: (optional) bool
            Whether to return the index of the closes nodes, default=False

        elements_of_interest: (opional) list
            list of element indices that are of interest

        k: (optional) int
            number of nearest neighbourt to return

        Returns
        -------------------------------
        coords: Nx3 array of floats
            coordinates of the baricenter of the closest element

        indexes: Nx1 array of ints
            Indice of the closest elements

        Notes
        --------------------------------------
        The indices are in the mesh listing, that starts at one!

        """
        if len(self.elm.node_number_list) == 0:
            raise ValueError('Mesh has no elements defined')

        baricenters = self.elements_baricenters()
        if elements_of_interest is not None:
            bar = baricenters[elements_of_interest]
        else:
            elements_of_interest = baricenters.elm_number
            bar = baricenters.value
        kd_tree = scipy.spatial.cKDTree(bar)
        d, indexes = kd_tree.query(querry_points, k=k)
        indexes = elements_of_interest[indexes]
        coords = baricenters[indexes]

        if return_distance and return_index:
            return coords, indexes, d

        elif return_index:
            return coords, indexes

        elif return_distance:
            return coords, d

        else:
            return coords

    def elm_node_coords(self, elm_nr=None, tag=None, elm_type=None):
        """ Returns the position of each of the element's nodes

        Arguments
        -----------------------------
        elm_nr: (optional) array of ints
            Elements to return, default: Return all elements
        tag: (optional) array of ints
            Only return elements with specified tag. default: all tags
        elm_type: (optional) array of ints
            Only return elements of specified type. default: all

        Returns
        -----------------------------
        Nx4x3 ndarray
            Array with node position of every element
            For triangles, the fourth coordinates are 0,0,0
        """
        elements_to_return = np.ones((self.elm.nr, ), dtype=bool)

        if elm_nr is not None:
            elements_to_return[elm_nr] = True

        if elm_type is not None:
            elements_to_return = np.logical_and(
                elements_to_return,
                np.in1d(self.elm.elm_type, elm_type))

        if tag is not None:
            elements_to_return = np.logical_and(
                elements_to_return,
                np.in1d(self.elm.tag1, tag))

        tmp_node_coord = np.vstack((self.nodes.node_coord, [0, 0, 0]))

        elm_node_coords = \
            tmp_node_coord[self.elm.node_number_list[elements_to_return, :] - 1]

        return elm_node_coords

    def write_hdf5(self, hdf5_fn, path='./'):
        """ Writes a HDF5 file with mesh information

        Parameters
        -----------
        hdf5_fn: str
            file name of hdf5 file
        path: str
            path in the hdf5 file where the mesh should be saved
        """
        with h5py.File(hdf5_fn, 'a') as f:
            try:
                g = f.create_group(path)
            except ValueError:
                g = f[path]
            if 'elm' in g.keys():
                del g['elm']
            if 'nodes' in g.keys():
                del g['nodes']
            if 'fields' in g.keys():
                del g['fields']
            g.attrs['fn'] = self.fn
            elm = g.create_group('elm')
            for key, value in vars(self.elm).items():
                elm.create_dataset(key, data=value)
            node = g.create_group('nodes')
            for key, value in vars(self.nodes).items():
                node.create_dataset(key, data=value)
            elmdata = g.create_group('elmdata')
            for d in self.elmdata:
                elmdata.create_dataset(d.field_name, data=d.value)
            nodedata = g.create_group('nodedata')
            for d in self.nodedata:
                nodedata.create_dataset(d.field_name, data=d.value)

    @classmethod
    def read_hdf5(self, hdf5_fn, path='./'):
        """ Reads mesh information from an hdf5 file

        Parameters
        ----------
        hdf5_fn: str
            file name of hdf5 file
        path: str
            path in the hdf5 file where the mesh is saved
        """
        import h5py
        self = self()
        with h5py.File(hdf5_fn, 'r') as f:
            g = f[path]
            try:
                self.fn = g.attrs['fn']
            except KeyError:
                pass
            for key, value in self.elm.__dict__.items():
                setattr(self.elm, key, np.squeeze(np.array(g['elm'][key])))
            for key, value in self.nodes.__dict__.items():
                try:
                    setattr(self.nodes, key, np.squeeze(np.array(g['nodes'][key])))
                except KeyError:
                    pass
            try:
                for field_name, field in g['elmdata'].items():
                    self.elmdata.append(
                        ElementData(np.squeeze(np.array(field)), field_name, mesh=self))
            except KeyError:
                pass

            try:
                for field_name, field in g['nodedata'].items():
                    self.nodedata.append(
                        NodeData(np.squeeze(np.array(field)), field_name, mesh=self))
            except KeyError:
                pass

        return self

    def open_in_gmsh(self):
        ''' Opens the mesh in gmsh '''
        f, fn_tmp = tempfile.mkstemp(suffix='.msh')
        os.close(f)
        open_in_gmsh(fn_tmp)
        os.remove(fn_tmp)

    def __eq__(self, other):
        try:
            return self.__dict__ == other.__dict__
        except AttributeError:
            return False

    def tetrahedra_quality(self, tetrahedra='all'):
        """ calculates the quality measures of the tetrahedra

        Parameters
        ------------
        tetrahedra: np.ndarray
            tags of the tetrahedra where the quality parameters are to be calculated

        Returns
        ----------
        measures: dict
            dictionary with ElementData with measures
        """
        if tetrahedra is 'all':
            tetrahedra = self.elm.tetrahedra
        if not np.all(np.in1d(tetrahedra, self.elm.tetrahedra)):
            raise ValueError('No tetrahedra with element number'
                             '{0}'.format(tetrahedra[
                                 np.logical_not(
                                     np.in1d(tetrahedra, self.elm.tetrahedra))]))
        M = self.nodes[self.elm[tetrahedra]]
        measures = {}
        V = self.elements_volumes_and_areas()[tetrahedra]
        E = np.array([
            M[:, 0] - M[:, 1],
            M[:, 0] - M[:, 2],
            M[:, 0] - M[:, 3],
            M[:, 1] - M[:, 2],
            M[:, 1] - M[:, 3],
            M[:, 2] - M[:, 3]])
        E = np.swapaxes(E, 0, 1)
        S = np.linalg.norm(E, axis=2)
        # calculate the circunstribed radius
        a = S[:, 0] * S[:, 5]
        b = S[:, 1] * S[:, 4]
        c = S[:, 2] * S[:, 3]
        s = 0.5 * (a + b + c)
        delta = np.sqrt(s * (s - a) * (s - b) * (s - c))
        CR = delta / (6 * V)
        # radius or inscribed sphere
        SA = np.linalg.norm([
            np.cross(E[:, 0], E[:, 1]),
            np.cross(E[:, 0], E[:, 2]),
            np.cross(E[:, 1], E[:, 2]),
            np.cross(E[:, 3], E[:, 4])], axis=2) * 0.5
        SA = np.swapaxes(SA, 0, 1)
        IR = 3 * V / np.sum(SA, axis=1)
        measures['beta'] = ElementData(CR / IR, 'beta')

        gamma = (np.sum(S * S, axis=1) / 6.) ** (3./2.) / V
        measures['gamma'] = ElementData(gamma, 'gamma')
        return measures

    def triangle_normals(self):
        """ Calculates the normals of triangles

        Returns
        --------
        normals: ElementData
            normals of triangles, zero at the tetrahedra

        """
        normals = ElementData(np.zeros((self.elm.nr, 3), dtype=float),
                              'normals')
        tr_indexes = self.elm.triangles
        node_tr = self.nodes[self.elm[tr_indexes, :3]]
        sideA = node_tr[:, 1] - node_tr[:, 0]

        sideB = node_tr[:, 2] - node_tr[:, 0]
        n = np.cross(sideA, sideB)
        normals[tr_indexes] = n / np.linalg.norm(n, axis=1)[:, None]
        return normals

    def nodes_volumes_or_areas(self):
        ''' Return the volume (volume mesh) if area (surface mesh) of all nodes
        Only works for ordered values of mesh and node indices

        Returns
        -------------
        nd: NodeData
            NodeData structure with the volume or area of each node
        '''
        nd = np.zeros(self.nodes.nr)
        if len(self.elm.tetrahedra) > 0:
            name = 'volumes'
            volumes = self.elements_volumes_and_areas()[self.elm.tetrahedra]
            th_nodes = self.elm[self.elm.tetrahedra] - 1
            for i in range(4):
                nd[:np.max(th_nodes[:, i]) + 1] += \
                    np.bincount(th_nodes[:, i], volumes / 4.)

        elif len(self.elm.triangles) > 0:
            name = 'areas'
            areas = self.elements_volumes_and_areas()[self.elm.triangles]
            tr_nodes = self.elm[self.elm.triangles] - 1
            for i in range(3):
                nd[:np.max(tr_nodes[:, i]) + 1] += \
                    np.bincount(tr_nodes[:, i], areas / 3.)

        return NodeData(nd, name)

    def nodes_areas(self):
        ''' Areas for all nodes in a surface

        Returns
        ---------
        nd: NodeData
            NodeData structure with normals for each node

        '''
        areas = self.elements_volumes_and_areas()[self.elm.triangles]
        triangle_nodes = self.elm[self.elm.triangles] - 1
        triangle_nodes = triangle_nodes[:, :3]
        nd = np.bincount(triangle_nodes.reshape(-1),
                         np.repeat(areas / 3., 3), self.nodes.nr) 

        return NodeData(nd, 'areas')


    def nodes_normals(self, smooth=0):
        ''' Normals for all nodes in a surface

        Parameters
        ------------
        smooth: int (optional)
            Number of smoothing cycles to perform. Default: 0

        Returns
        ---------
        nd: NodeData
            NodeData structure with normals for each node

        '''
        nodes = np.unique(self.elm[self.elm.triangles, :3])
        elements = self.elm.triangles

        nd = np.zeros((self.nodes.nr, 3))

        node_tr = self.nodes[self.elm.node_number_list[elements - 1, :3]]
        sideA = node_tr[:, 1] - node_tr[:, 0]

        sideB = node_tr[:, 2] - node_tr[:, 0]
        normals = np.cross(sideA, sideB)

        triangle_nodes = self.elm.node_number_list[elements - 1, :3] - 1
        for s in range(smooth + 1):
            for i in range(3):
                nd[:, i] = \
                    np.bincount(triangle_nodes.reshape(-1),
                                np.repeat(normals[:, i], 3),
                                self.nodes.nr)

            normals = np.sum(nd[self.elm.node_number_list[elements - 1, :3] - 1], axis=1)
            normals /= np.linalg.norm(normals, axis=1)[:, None]

        nd[nodes - 1] = nd[nodes-1] / \
            np.linalg.norm(nd[nodes-1], axis=1)[:, None]
 
        return NodeData(nd, 'normals')

    def find_tetrahedron_with_points(self, points, compute_baricentric=True):
        ''' Finds the tetrahedron that contains each of the described points using a
        stochastic walk algorithm

        Parameters
        ----------------
        points: Nx3 ndarray
            List of points to be queried

        compute_baricenters: bool
            Wether or not to compute baricentric coordinates of the points

        Returns
        ----------------
        th_with_points: ndarray
            List with the tetrahedron that contains each point. If the point is outside
            the mesh, the value will be -1
        baricentric: n x 4 ndarray (if compute_baricentric == True)
            Baricentric coordinates of point. If the point is outside, a list of zeros

        References
        ------------------
        Devillers, Olivier, Sylvain Pion, and Monique Teillaud. "Walking in a
        triangulation." International Journal of Foundations of Computer Science 13.02
        (2002): 181-199.
        '''
        #
        th_indices = self.elm.tetrahedra
        th_nodes = self.nodes[self.elm[th_indices]]
        # Reduce the number of elements
        points_max = np.max(points, axis=0)
        points_min = np.min(points, axis=0)
        th_max = np.max(th_nodes, axis=1)
        th_min = np.min(th_nodes, axis=1)
        slack = (points_max - points_min) * .05
        th_in_box = np.where(
            np.all(
                (th_min <= points_max + slack) * (th_max >= points_min - slack), axis=1))[0]
        th_indices = th_indices[th_in_box]
        # if all the points are outside the bounding box
        if len(th_indices) == 0: 
            if compute_baricentric:
                return -np.ones(len(points), dtype=int), np.zeros_like(points)
            else:
                return -np.ones(len(points), dtype=int)
        th_nodes = th_nodes[th_in_box]

        # Calculate a few things we will use later
        faces, th_faces, adjacency_list = self.elm.get_faces(th_indices)
        # Find initial positions
        th_baricenters = np.average(th_nodes, axis=1)
        kdtree = scipy.spatial.cKDTree(th_baricenters)

        # Starting position for walking algorithm: the closest baricenter
        _, closest_th = kdtree.query(points)
        pts = np.array(points, dtype=np.float)
        th_nodes = np.array(th_nodes, dtype=np.float)
        closest_th = np.array(closest_th, dtype=np.int)
        th_faces = np.array(th_faces, dtype=np.int)
        adjacency_list = np.array(adjacency_list, dtype=np.int)
        th_with_points = cython_msh.find_tetrahedron_with_points(
            pts, th_nodes, closest_th, th_faces, adjacency_list)

        # calculate baricentric coordinates
        inside = th_with_points != -1
        if compute_baricentric:
            M = np.transpose(th_nodes[th_with_points[inside], :3, :3] -
                             th_nodes[th_with_points[inside], 3, None, :], (0, 2, 1))
            baricentric = np.zeros((len(points), 4), dtype=float)
            baricentric[inside, :3] = np.linalg.solve(
                M, points[inside] - th_nodes[th_with_points[inside], 3, :])
            baricentric[inside, 3] = 1 - np.sum(baricentric[inside], axis=1)

        # Return indices
        th_with_points[inside] = \
            th_indices[th_with_points[inside]]

        if compute_baricentric:
            return th_with_points, baricentric
        else:
            return th_with_points

        '''
        th_with_points = -np.ones(points.shape[0], dtype=int)
        # We can now start the walking algorithm
        face_order = range(4)
        for i, p, t in zip(range(points.shape[0]), points, closest_th):
            previous_t = False
            end = False
            while not end:
                np.random.shuffle(face_order)
                end = True
                for f in face_order:
                    orient = orientation(p, t, f)
                    if (orient < 0 and
                            not np.any(adjacency_list[th_faces[t, f]] == previous_t)):
                        adj = adjacency_list[th_faces[t, f]]
                        previous_t = t
                        t = adj[adj != t]
                        if t == -1:
                            end = True
                        else:
                            end = False
                        break
            th_with_points[i] = th_indices[t] if t != -1 else t
        '''
    def test_inside_volume(self, points):
        ''' Tests if points are iside the volume using the Möller–Trumbore intersection
        algorithm

        Notice: This algorithm is vulnerable to degenerate cases (example: when the ray
        goes right through an edge
        Parameters
        ----------------
        points: Nx3 ndarray
            List of points to be queried

        Returns
        ----------------
        inside: ndarray
            array of booleans

        '''
        triangles = self.nodes[self.elm.get_outside_faces()]
        quantities = cython_msh.calc_quantities_for_test_point_in_triangle(triangles)
        points = np.array(points, dtype=np.float)
        inside = cython_msh.test_point_in_triangle(points, *quantities)

        return inside

    def _check_th_node_ordering(self):
        th_nodes = self.nodes[self.elm[self.elm.tetrahedra]]
        th_baricenters = np.average(th_nodes, axis=1)
        face_points = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
        for t, b in enumerate(th_baricenters):
            for i in range(4):
                d = th_nodes[t, face_points[i]] - b
                assert np.linalg.det(d) > 0, \
                    'Found a face pointing the wrong way!'

    def fix_th_node_ordering(self):
        ''' Fixes the node ordering of tetrahedra in-place '''
        th = self.nodes[self.elm.node_number_list[self.elm.elm_type == 4, :]]
        M = th[:, 1:] - th[:, 0, None]
        switch = self.elm.elm_type == 4
        switch[switch] = np.linalg.det(M) < 0
        tmp = np.copy(self.elm.node_number_list[switch, 1])
        self.elm.node_number_list[switch, 1] = self.elm.node_number_list[switch, 0]
        self.elm.node_number_list[switch, 0] = tmp
        del tmp
        gc.collect()

    def fix_tr_node_ordering(self):
        ''' Fixes the node ordering of the triangles in-place '''
        corresponding = self.find_corresponding_tetrahedra()
        triangles = np.where(self.elm.elm_type == 2)[0]

        triangles = triangles[corresponding != -1]
        corresponding = corresponding[corresponding != -1]

        normals = self.triangle_normals().value[triangles]
        baricenters = self.elements_baricenters().value
        pos_bar = baricenters[corresponding - 1] - baricenters[triangles]

        dotp = np.einsum('ij, ij -> i', normals, pos_bar)
        switch = triangles[dotp > 0]

        tmp = np.copy(self.elm.node_number_list[switch, 1])
        self.elm.node_number_list[switch, 1] = self.elm.node_number_list[switch, 0]
        self.elm.node_number_list[switch, 0] = tmp
        del tmp
        gc.collect()


    def fix_surface_labels(self):
        ''' Fixels labels of surfaces '''
        change = (self.elm.elm_type == 2) * (self.elm.tag1 < 1000)
        self.elm.tag1[change] += 1000
        change = (self.elm.elm_type == 2) * (self.elm.tag2 < 1000)
        self.elm.tag2[change] += 1000

    def compact_ordering(self, node_number):
        ''' Changes the node and element ordering so that it goes from 1 to nr_nodes
        
        Parameters
        --------------
        node_number: N_nodes x 1 ndarray of inte
            Node numbering in the original mesh
        '''
        rel = int(-9999) * np.ones((np.max(node_number) + 1), dtype=int)
        rel[node_number] = np.arange(1, self.nodes.nr + 1, dtype=int)
        self.elm.node_number_list = rel[self.elm.node_number_list]
        self.elm.node_number_list[self.elm.elm_type == 2, 3] = -1

    def prepare_surface_tags(self):
        triangles = self.elm.elm_type == 2
        surf_tags = np.unique(self.elm.tag1[triangles])
        for s in surf_tags:
            if s < 1000:
                self.elm.tag1[triangles *
                              (self.elm.tag1 == s)] += 1000

    def find_corresponding_tetrahedra(self):
        ''' Finds the tetrahedra corresponding to each triangle

        Returns
        ---------
        corresponding_th_indices: ndarray of ints
           List of the element indices of the tetrahedra corresponding to each triangle.
           Note: This is in mesh ordering (starts at 1), -1 if there's no corresponding

        '''
        # Look into the cache
        node_nr_list_hash = hashlib.sha1(
            np.hstack((self.elm.tag1[:, None], self.elm.node_number_list))).hexdigest()
        try:
            if self._correspondance_node_nr_list_hash == node_nr_list_hash:
                return self._corresponding_tetrahedra
            else:
                raise AttributeError

        except AttributeError:
            pass

        # If could not find correspondence in cache
        tr_tags = np.unique(self.elm.tag1[self.elm.elm_type == 2])
        corresponding_th_indices = -np.ones(len(self.elm.triangles), dtype=int)
        for t in tr_tags:
            # look into tetrahedra with tags t, 1000-t
            to_crop = [t]
            to_crop.append(t - 1000)
            if t >= 1100:  # Electrodes
                to_crop.append(t - 600)
            if t >= 2000:
                to_crop.append(t - 2000)
                to_crop.append(t - 1600)
            # Select triangles and tetrahedra with tags
            th_of_interest = np.where((self.elm.elm_type == 4) *
                                      np.in1d(self.elm.tag1, to_crop))[0]
            if len(th_of_interest) == 0:
                continue
            tr_of_interest = np.where((self.elm.elm_type == 2) *
                                      (self.elm.tag1 == t))[0]

            th = self.elm.node_number_list[th_of_interest]
            faces = th[:, [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]].reshape(-1, 3)
            faces = faces.reshape(-1, 3)
            faces_hash_array = _hash_rows(faces)

            tr = self.elm.node_number_list[tr_of_interest, :3]
            tr_hash_array = _hash_rows(tr)
            # This will check if all triangles have a corresponding face
            has_tetra = np.in1d(tr_hash_array, faces_hash_array)
            faces_argsort = faces_hash_array.argsort()
            faces_hash_array = faces_hash_array[faces_argsort]
            tr_search = np.searchsorted(faces_hash_array, tr_hash_array[has_tetra])
            # Find the indice
            index = faces_argsort[tr_search] // 4

            # Put the values in corresponding_th_indices
            position = np.searchsorted(self.elm.triangles,
                                       self.elm.elm_number[tr_of_interest[has_tetra]])
            corresponding_th_indices[position] = \
                self.elm.elm_number[th_of_interest[index]]

        self._correspondance_node_nr_list_hash = node_nr_list_hash
        self._corresponding_tetrahedra = corresponding_th_indices
        gc.collect()
        return corresponding_th_indices

    def fix_thin_tetrahedra(self, n_iter=7, step_length=.05):
        ''' Optimize the locations of the points by moving them towards the center
        of their patch. This is done iterativally for all points for a number of
        iterations'''
        vol_before = self.elements_volumes_and_areas()[self.elm.tetrahedra].sum()
        boundary_faces = []
        vol_tags = np.unique(self.elm.tag1)
        for t in vol_tags:
            th_indexes = self.elm.elm_number[
                (self.elm.tag1 == t) * (self.elm.elm_type == 4)]
            boundary_faces.append(
                self.elm.get_outside_faces(
                    tetrahedra_indexes=th_indexes))
        boundary_faces = np.vstack(boundary_faces)
        boundary_nodes = np.unique(boundary_faces) - 1
        mean_bar = np.zeros_like(self.nodes.node_coord)
        nodes_bk = np.copy(self.nodes.node_coord)
        new_nodes = self.nodes.node_coord
        tetrahedra = self.elm[self.elm.tetrahedra] - 1
        k = np.bincount(tetrahedra.reshape(-1), minlength=self.nodes.nr)
        for n in range(n_iter):
            bar = np.mean(new_nodes[tetrahedra], axis=1)
            for i in range(3):
                mean_bar[:, i] = np.bincount(tetrahedra.reshape(-1),
                                             weights=np.repeat(bar[:, i], 4),
                                             minlength=self.nodes.nr)
            mean_bar /= k[:, None]
            new_nodes += step_length * (mean_bar - new_nodes)
            new_nodes[boundary_nodes] = nodes_bk[boundary_nodes]
        self.nodes.node_coord = new_nodes
        vol_after = self.elements_volumes_and_areas()[self.elm.tetrahedra].sum()
        # Cheap way of comparing the valitidy of the new trianglulation (assuming the one
        # in the input is valid)
        if not np.isclose(vol_before, vol_after):
            self.nodes.node_coord = nodes_bk

    def calc_matsimnibs(self, center, pos_ydir, distance, skin_surface=[5, 1005]):
        ''' Calculate the matsimnibs matrix for TMS simulations

        Parameters
        -----------
        center: np.ndarray
            Position of the center of the coil, will be projected to the skin surface
        pos_ydir: np.ndarray
            Position of the y axis in relation to the coil center
        distance: float
            Distance from the center
        skin_surface: list
            Possible tags for the skin surface (Default: [5, 1005])

        Returns
        -------
        matsimnibs: 2d np.ndarray
            Matrix of the format
            x' y' z' c
            0  0  0  1
            y' is the direction of the coil
            z' is a direction normal to the coil, points inside the head

        '''
        msh_surf = self.crop_mesh(elm_type=2)
        msh_skin = msh_surf.crop_mesh(skin_surface)
        closest = np.argmin(np.linalg.norm(msh_skin.nodes.node_coord - center, axis=1))
        center = msh_skin.nodes.node_coord[closest]
        # Y axis
        y = pos_ydir - center
        if np.isclose(np.linalg.norm(y), 0.):
            raise ValueError('The coil Y axis reference is too close to the coil center! ')
        y /= np.linalg.norm(y)
        #Normal
        normal = msh_skin.nodes_normals().value[closest]
        if np.isclose(np.abs(y.dot(normal)), 1.):
            raise ValueError('The coil Y axis normal to the surface! ')
        z = -normal
        #Orthogonalize y
        y -= z * y.dot(z)
        y /= np.linalg.norm(y)
        # Determine x
        x = np.cross(y, z)
        # Determine C
        c = center + distance * normal
        # matsimnibs
        matsimnibs = np.zeros((4, 4), dtype=float)
        matsimnibs[:3, 0] = x
        matsimnibs[:3, 1] = y
        matsimnibs[:3, 2] = z
        matsimnibs[:3, 3] = c
        matsimnibs[3, 3] = 1
        return matsimnibs

    def elm2node_matrix(self, th_indices=None):
        ''' Calculates a sparse matrix to tranform from ElementData to NodeData
        Uses Superconvergent patch recovery for volumetric data. Will not work well for discontinuous fields (like E, if
        several tissues are used)

        Returns
        ---------
        M: scipy.sparse.csc
            Sparse matrix, where M.dot(elm_data) = node_data, elm_data is a vector with
            element data, and node_data is a vector with node data.
        th_indices: np.ndarray (optional)
            Indices of the tetrahedra to be considered in the volume. Default: use all
            tetrahedra

        References
        -----------
            Zienkiewicz, Olgierd Cecil, and Jian
            Zhong Zhu. "The superconvergent patch recovery and a posteriori error
            estimates. Part 1: The recovery technique." International Journal for
            Numerical Methods in Engineering 33.7 (1992): 1331-1364.

        '''
        if th_indices is not None:
            th_indices = np.intersect1d(self.elm.tetrahedra, th_indices)
        else:
            th_indices = self.elm.tetrahedra

        if len(th_indices) == 0:
            raise ValueError("Can only create elm2node matrices for tetrahedral meshes")

        # Get the point in the outside surface
        points_outside = np.unique(self.elm.get_outside_faces(th_indices))
        outside_points_mask = np.in1d(
            self.elm[th_indices], points_outside).reshape(-1, 4)
        # Get indices starting at zero
        points_outside = points_outside - 1
        th_nodes = self.elm[th_indices] - 1
        # This will map all points outside to the first position
        masked_th_nodes = np.copy(th_nodes)
        masked_th_nodes[outside_points_mask] = -1
        masked_th_nodes += 1

        # Calculates the quantities needed for the superconvergent patch recovery
        baricenters = self.elements_baricenters()[th_indices]
        volumes = self.elements_volumes_and_areas()[th_indices]
        baricenters = np.hstack(
            [np.ones((baricenters.shape[0], 1)), baricenters])

        # NOTICE: Bellow, I will add everything in the outer nodes to the
        # first row. I will remove it afterwards
        A = np.zeros((self.nodes.nr + 1, 4, 4))
        for i in range(4):
            for j in range(i, 4):
                A[:, i, j] = np.bincount(
                    masked_th_nodes.reshape(-1),
                    np.repeat(baricenters[:, i], 4) *
                    np.repeat(baricenters[:, j], 4),
                    minlength=self.nodes.nr + 1)

        # This here only ensures we can invert
        outside = np.isclose(A[:, 0, 0], 0)
        for i in range(4):
            A[outside, i, i] = 1

        A[:, 1, 0] = A[:, 0, 1]
        A[:, 2, 0] = A[:, 0, 2]
        A[:, 3, 0] = A[:, 0, 3]
        A[:, 2, 1] = A[:, 1, 2]
        A[:, 3, 1] = A[:, 1, 3]
        A[:, 3, 2] = A[:, 2, 3]

        Ainv = np.linalg.inv(A)

        M = scipy.sparse.csr_matrix((self.nodes.nr + 1, self.elm.nr))
        node_pos = np.hstack(
            [np.ones((self.nodes.nr, 1)), self.nodes.node_coord])
        # Added a dummy to the first position
        node_pos = np.vstack([np.ones((1, 4)), node_pos])
        for i in range(4):
            M += scipy.sparse.csr_matrix(
                (np.einsum(
                    'bi, bij, bj -> b',
                    node_pos[masked_th_nodes[:, i]],
                    Ainv[masked_th_nodes[:, i]],
                    baricenters),
                 (masked_th_nodes[:, i], th_indices - 1)),
                shape=M.shape)

        # Assigns the average value to the points in the outside surface
        masked_th_nodes = np.copy(th_nodes)
        masked_th_nodes[~outside_points_mask] = -1
        masked_th_nodes += 1

        for i in range(4):
            M += scipy.sparse.csr_matrix(
                 (volumes, (masked_th_nodes[:, i], th_indices - 1)),
                 shape=M.shape)

        M = M[1:]

        node_vols = np.bincount(
            th_nodes.reshape(-1),
            np.repeat(volumes, 4),
            minlength=self.nodes.nr+1)

        normalization = np.ones(self.nodes.nr)
        normalization[points_outside] = 1 / node_vols[points_outside]

        D = scipy.sparse.dia_matrix(
            (normalization, 0), shape=(self.nodes.nr, self.nodes.nr))
        M = D.dot(M)

        return M

    def interp_matrix(self, pos, out_fill=np.nan, th_indices=None, element_wise=False):
        ''' Calculates a matrix to perform interpolation
        y = M.dot(x)

        Where x is node-wise data (if element_wise=False, default) or element-wise data
        otherwise

        Parameters
        ----------
        pos: (N x 3) np.ndarray
            Positions where interpolation is to be performed

        out_fill: float or 'nearest'
            How to fill values for positions outside the volume

        th_indices: np.ndarray (optional)
            Indices of the tetrahedra to be considered in the volume. Default: use all
            tetrahedra

        element_wise: bool (optional)
            Wether to do interpolations one element-wise data. Default=False

        Returns
        -------
        M: scipy.sparse.csc
            Sparse matrix, interpolation represented by dot product

        '''
        if len(self.elm.tetrahedra) == 0:
            raise ValueError("Can only create interpolation matrices for tetrahedral meshes")

        th_with_points, bar = self.find_tetrahedron_with_points(
            pos, compute_baricentric=True)
        if th_indices is not None:
            th_with_points[~np.isin(th_with_points, th_indices)] = -1
        inside = th_with_points != -1
        pos_nr = np.arange(len(pos))
        th_nodes = self.elm[th_with_points[inside]]

        M = scipy.sparse.csc_matrix((len(pos), self.nodes.nr))
        for i in range(4):
            M += scipy.sparse.csc_matrix(
                (bar[inside, i],
                 (pos_nr[inside], th_nodes[:, i] - 1)),
                shape=M.shape)

        if out_fill != 'nearest':
            v = out_fill * np.ones(np.sum(~inside))
            M += scipy.sparse.csc_matrix(
                (v, (pos_nr[~inside], np.zeros(np.sum(~inside)))),
                shape=M.shape)

        else:
            _, nearest = self.nodes.find_closest_node(
                    pos[~inside], return_index=True)
            M += scipy.sparse.csc_matrix(
                (np.ones(np.sum(~inside)), (pos_nr[~inside], nearest-1)),
                shape=M.shape)

        if element_wise:
            if th_indices is None:
                M = M.dot(self.elm2node_matrix())
            else:
                M = M.dot(self.elm2node_matrix(th_indices))

        return M


    def view(self,
             visible_tags=None,
             visible_fields=[]):
        ''' Visualize mesh in Gmsh

        Parameters
        ------------
        visible_tags: list (optional)
            List of tags to be visible. Default: all tags
        visible_fields: list (optional) or 'all'
            Name of visible fields or 'all' to view all fields. Default: no fields visible

        Returns
        --------
        vis: simnibs.gmsh_view.Visualization
            Visualization object

        Example
        ---------
        >>> mesh = simnibs.msh.read_msh('ernie.msh')
        >>> vis = mesh.view()
        >>> vis.show()
        '''
        vis = gmsh_view.Visualization(self)
        if visible_tags is not None:
            vis.visibility = visible_tags
        vis.View = []
        if visible_fields == 'all':
            visible_fields = self.field.keys()

        for i, d in enumerate(self.nodedata + self.elmdata):
            vis.View.append(
                d.view_options(
                   visible=d.field_name in visible_fields,
                   visible_tags=visible_tags,
                   idx=i))

        if len(visible_fields) > 0:
            vis.Mesh.SurfaceEdges = 0
            vis.Mesh.SurfaceFaces = 0
            vis.Mesh.VolumeEdges = 0
            vis.Mesh.VolumeFaces = 0

        return vis

    def add_node_field(self, field, field_name):
        ''' Adds field defined in the nodes

        Parameters
        ----------
        field: np.ndarray or simnibs.NodeData
            Value of field in the mesh nodes. Should have the shape
            - (n_nodes,) or (n_nodes, 1) for scalar fields
            - (n_nodes, 3) for vector fields
            - (n_nodes, 9) for tensors

        field_name: str
            Name for the field.

        Returns
        --------
        nd: simnibs.NodeData
            NodeData class with the input field

        '''
        if isinstance(field, NodeData):
            if field.nr != self.nodes.nr:
                raise ValueError('Number of data points in the field '
                                 'and of mesh nodes do not match')
            field.field_name = field_name
            field.mesh = self
            self.nodedata.append(field)

            return field

        else:
            if self.nodes.nr != field.shape[0]:
                raise ValueError('Number of data points in the field '
                                 'and of mesh nodes do not match')

            nd = NodeData(field, field_name, self)
            self.nodedata.append(nd)
            return nd

    def add_element_field(self, field, field_name):
        ''' Adds field defined in the elements

        Parameters
        ----------
        field: np.ndarray or simnibs.ElementData
            Value of field in the mesh elements. Should have the shape
            - (n_elm,) or (n_elm, 1) for scalar fields
            - (n_elm, 3) for vector fields
            - (n_elm, 9) for tensors

        field_name: str
            Name for the field.

        Returns
        --------
        ed: simnibs.ElementData
            ElementData class with the input field

        '''
        if isinstance(field, ElementData):
            if field.nr != self.elm.nr:
                raise ValueError('Number of data points in the field '
                                 'and number of mesh elements do not match')
            field.field_name = field_name
            field.mesh = self
            self.elmdata.append(field)

            return field

        else:
            if self.elm.nr != field.shape[0]:
                raise ValueError('Number of data points in the field '
                                 'and number of mesh elements do not match')
            ed = ElementData(field, field_name, self)
            self.elmdata.append(ed)
            return ed

    def fields_summary(self, roi=None, fields=None,
                       percentiles=(99.9, 99, 95),
                       focality_cutoffs=(75, 50)):
        ''' Creates a text summaty of the field

        Parameters
        ------------
        roi: list (optional)
            Regions of interest, in tissue tags. Default: Whole mesh
        fields: list (optional)
            Fields for which to calculate the summary
        percentiles: ndarray (optinal)
            Field percentiles to be printed. Default: (99.9, 99, 95)
        focality_cutoffs: ndarray (optional)
            Cuttofs for focality calculations. Default: (75, 50)
        '''
        if roi is None:
            mesh = self
        else:
            try:
                mesh = self.crop_mesh(roi)
            except ValueError:
                warnings.warn(f"Could not find any element with tags {roi}")
                return f"Could not find any element with tags {roi}"

        if fields is None:
            fields = self.field.keys()

        units = []
        for f in fields:
            if f in ['E', 'normE', 'D', 'g']:
                units.append(' V/m')
            elif f in ['J', 'normJ']:
                units.append(' A/m²')
            elif f == 'v':
                units.append(' V')
            else:
                units.append('')

        if 2 in mesh.elm.elm_type:
            units_mesh = ' mm²'
        elif 4 in mesh.elm.elm_type:
            units_mesh = ' mm³'
        if np.all(np.isin([2, 4], mesh.elm.elm_type)):
            warnings.warn("Can't report Field summary in meshes with volumes and surfaces")
            return ''

        percentiles_table = [['Field'] + [f'{p:.1f}%' for p in percentiles]]
        focality_table = [['Field'] + [f'{f:.1f}%' for f in focality_cutoffs]]
        for fn, u in zip(fields, units):
            f = mesh.field[fn]
            prc = f.get_percentiles(percentiles)
            percentiles_table.append([fn] + [f'{p:.2e}{u}' for p in prc])
            focality = f.get_focality(focality_cutoffs, 99.9)
            focality_table.append([fn] + [f'{fv:.2e}{units_mesh}' for fv in focality])
        
        def format_table(table):
            entry_sizes = np.array([[len(e) for e in row] for row in table])
            col_sizes = np.max(entry_sizes, axis=0)
            align_string = ['{:<' + str(cs) + '}' for cs in col_sizes]
            t = ''
            for i, row in enumerate(table):
                t += '|'
                t += ' |'.join(a_s.format(col) for a_s, col in zip(align_string, row))
                t += ' |\n'
                if i == 0:
                    t += '|' + '|'.join((cs+1) * '-' for cs in col_sizes) + '|\n'
            return t

        string = ''
        string += 'Field Percentiles\n'
        string += '-----------------\n'
        string += 'Top percentiles of the field (or field norm for vector fields)\n'
        string += format_table(percentiles_table)
        string += '\n'
        string += 'Field Focality\n'
        string += '---------------\n'
        string += 'Mesh volume or area with a field >= X% of the 99.9th percentile\n'
        string += format_table(focality_table)

        return string



class Data(object):
    """Store data in elements or nodes

    Parameters
    -----------------------
    value: np.ndarray
        Value of field in nodes

    field_name: str (optional)
        name of field. Default: empty

    mesh: simnibs.msh.Msh (optional)
        Mesh where the field is define. Required for several methods

    Attributes
    --------------
    value: ndarray
        Value of field in nodes
    field_name: str
        name of field
    nr: property
        number of data points
    nr_comp: property
        number of dimensions per data point (1 for scalars, 3 for vectors)

    """

    def __init__(self, value, name='', mesh=None):
        self.field_name = name
        self.value = value
        self.mesh = mesh

        if value.ndim > 2:
            raise ValueError('Can only hadle 1 and 2 dimensional fields '
                             'Tensors should be given as a Nx9 array')

        if self.nr_comp > self.nr:
            warnings.warn('Second axis larger than the first '
                          'Field is probably transposed')

    @property
    def type(self):
        '''NodeData of ElementData'''
        return self.__class__.__name__

    @property
    def nr(self):
        '''Number of data entries'''
        return self.value.shape[0]

    @property
    def nr_comp(self):
        '''Number of field components'''
        try:
            return self.value.shape[1]
        except IndexError:
            return 1

    def interpolate_to_surface(self, surface, out_fill='nearest'):
        ''' Interpolates the field in the nodes of a given surface
        The interpolation occurs in the tetrahedra!

        Parameters
        -----------
        surface: Msh
            Mesh structure with triangles only
        out_fill: float or 'nearest' (optional)
            Value to be assigned to the points in the surface outside the volume.

        Returns
        ---------
        node_data: NodeData
            Node data structure with the interpolated field
        '''
        interp = self.interpolate_scattered(surface.nodes.node_coord, out_fill=out_fill)
        return NodeData(interp, name=self.field_name, mesh=surface)

    def to_nifti(self, n_voxels, affine, fn=None, units='mm', qform=None,
                 method='linear', continuous=False):
        ''' Transforms the data in a nifti file

        Parameters
        -----------
        n_voxels: list of ints
            Number of vexels in each dimension
        affine: 4x4 ndarray
            Transformation of voxel space into xyz. This sets the sform
        fn: str (optional)
            String with file name to be used, if the result is to be saved
        units: str (optional)
            Units to be set in the NifTI header. Default: mm
        qform: 4x4 ndarray (optional)
            Header qform. Default: set the same as the affine
        method: {'assign' or 'linear'} (Optional)
            If 'assign', gives to each voxel the value of the element that contains
            it. If linear, first assign fields to nodes, and then perform
            baricentric interpolatiom. Only for ElementData input. Default: linear
        continuous: bool
            Wether fields is continuous across tissue boundaries. Changes the
            behaviour of the function only if method == 'linear'. Default: False

        Returns
        ---------
        img: nibabel.Nifti1Pair
            Image object with the field interpolated in the voxels
        '''
        data = self.interpolate_to_grid(n_voxels, affine, method=method,
                                        continuous=continuous)
        if data.dtype == np.bool or data.dtype == bool:
            data = data.astype(np.uint8)
        img = nibabel.Nifti1Pair(data, affine)
        img.header.set_xyzt_units(units)
        if qform is not None:
            img.set_qform(qform)
        else:
            img.set_qform(affine)
        img.update_header()
        del data
        if fn is not None:
            nibabel.save(img, fn)
        else:
            return img

    def to_deformed_grid(self, warp, reference, out=None,
                         out_original=None, tags=None, order=1,
                         method='linear', continuous=False,
                         inverse_warp=None, reference_original=None,
                         binary=False):
        ''' Interpolates field to a grid and apply non-linear interpolation

        We first interpolate to a grid and then apply the transformation in order to
        avoid problems from deformed triangles

        Parameters
        ------------
        warp: str
            Name of file with the transformation. Can either be a nifti file
            where each voxel corresponds to coordinates (x, y, z) in the original space
            or an affine transformation defined from the target space to the original
            space. In the later case, the name must finish in ".mat", and it will be read
            with the numpy loadtxt file
        ref: str
            Name of reference file. The output will be in the same space
            as the reference (same affine transfomation and same dimensions)
        out: str (optional)
            If not None, the result will be written to this file as a nifti
        out_original: str (optional)
            If not None, the volume in the original grid be written to this file as a nifti
        tags: list (optional)
            Mesh tags to be transformed. Defaut: transform the entire mesh
        order: int (optional)
            Interpolation order to be used. Default: 1
        method: {'assign' or 'linear'} (Optional)
            Method for gridding the data.
            If 'assign', gives to each voxel the value of the element that contains
            it. If linear, first assign fields to nodes, and then perform
            baricentric interpolatiom. Only for ElementData input. Default: linear
        continuous: bool
            Wether fields is continuous across tissue boundaries. Changes the
            behaviour of the function only if method == 'linear'. Default: False
        inverse_warp: str
            Name of nifti file with inverse the transformation. Used to rotate vectors to the
            target space in the case of non-liner transformations. If the transformation is
            linear, the inverse matrix is used.
        reference_original: str
            Name of nifti file with reference in the original space. Used to determine
            the dimensions and affine transformation for the initial griding

        Returns:
        --------
        img: nibabel.Nifti1Pair
            Nibabel image object with tranformed field
        '''
        self._test_msh()
        # Figure out a good space where to grid the data
        if tags is not None:
            msh = copy.deepcopy(self.mesh)
            if isinstance(self, ElementData):
                msh.nodedata = []
                msh.elmdata = [self]
                msh = msh.crop_mesh(tags)
                data = msh.elmdata[0]
            elif isinstance(self, NodeData):
                msh.nodedata = [self]
                msh.elmdata = []
                msh = msh.crop_mesh(tags)
                data = msh.nodedata[0]
        else:
            data = self
        # Get the affine transformation and the dimensions of the image in original space
        if reference_original is None:
            reference_nifti = nibabel.load(reference)
        else:
            reference_nifti = nibabel.load(reference_original)
        affine = reference_nifti.affine
        dimensions = reference_nifti.shape[:3]
        image = data.interpolate_to_grid(dimensions, affine, method=method,
                                         continuous=continuous)
        if len(image.shape) == 3:
            image = image[..., None]

        if out_original is not None:
            img = nibabel.Nifti1Pair(image, affine)
            img.header.set_xyzt_units('mm')
            img.set_qform(affine)
            nibabel.save(img, out_original)

        img = nifti_transform(
            (image, affine),
            warp, reference, out=out, order=order,
            inverse_warp=inverse_warp, binary=binary)

        del image
        gc.collect()
        return img

    def _norm(self):
        ''' simple norm of the field '''
        if self.nr_comp == 1:
            return np.abs(self.value).reshape(-1)
        else:
            return np.linalg.norm(self.value, axis=1)

    def _weights(self, roi=slice(None)):
        ''' Area / volume of each nodes or element '''
        if isinstance(self, NodeData):
            return self.mesh.nodes_volumes_or_areas()[roi]

        elif isinstance(self, ElementData):
            return self.mesh.elements_volumes_and_areas()[roi]

        else:
            raise NotImplementedError


    def mean_field_norm(self):
        ''' Calculates V*w/sum(w)
        Where V is the norm of the field, and w is the volume or area of the mesh where
        the field is defined. This can be used as a focality metric. It should give out
        small values when the field is focal and

        Returns
        ----------
        eff_area: float
            Area or volume of mesh, weighted by the field
        '''
        self._test_msh()
        if np.all(np.isin([2, 4], self.mesh.elm.elm_type)):
            warnings.warn('Calculating effective volume/area of fields in meshes with'
                          ' triangles and tetrahedra can give misleading results')

        norm = self._norm()
        weights = self._weights()

        return np.sum(norm * weights) / np.sum(weights)

    def get_percentiles(self, percentile=[99.9], roi=None):
        ''' Get percentiles of field (or field norm, if a vector field)

        Parameters
        ------------
        percentile: ndarray (optional)
            Percentiles of interest, between 0 and 100. Defaut: 99.9

        roi: ndarray (optinal)
            Region of interest in terms of element/node indices. Default: the whole mesh

        Returnts
        ----------
        f_p: ndarray
            Field at the given percentiles
        '''
        self._test_msh()
        if roi is None:
            roi = slice(None)

        if self.nr_comp > 1:
            v = np.linalg.norm(self[roi], axis=1)
        else:
            v = np.squeeze(self[roi])
        s = np.argsort(v)
        v = v[s]
        weights = self._weights(roi)[s]
        weights = np.cumsum(weights)
        weights /= weights[-1]
        perc = np.array(percentile, dtype=float) / 100
        perc = perc.reshape(-1)
        closest = np.zeros(perc.shape, dtype=int)
        for i, p in enumerate(perc):
            closest[i] = np.argmin(np.abs(weights - p))

        return v[closest]

    def get_focality(self, cuttofs=[50, 70], peak_percentile=99.9):
        ''' Caluclates field focality as the area/volume of the mesh experiencing a field
        norm of above (cut_off% of the field peak). peak_percentile gives what is the
        field peak

        Parameters
        ------------
        cuttofs: ndarray (optional)
            Percentage of the peak value for the cut_off, between 0 and 100. Default: [50, 70]

        peak_percentile: float (optional)
            Percentile to be used to calculate peak value. Default: 99.9
        
        Returns
        ---------
        focality: ndarray
            Area/volume exceeding the cuttof of the peak value
        '''
        self._test_msh()
        if np.all(np.isin([2, 4], self.mesh.elm.elm_type)):
            warnings.warn('Calculating focality of fields in meshes with'
                          ' triangles and tetrahedra can give misleading results')

        norm = self._norm()
        s = np.argsort(norm)
        norm = norm[s]
        weights = self._weights()[s]
        weights_norm = np.cumsum(weights)
        weights_norm = weights_norm * 100 / weights_norm[-1]
        peak_value = norm[np.argmin(np.abs(weights_norm - peak_percentile))]

        co = np.array(cuttofs, dtype=float).reshape(-1) / 100
        focality = np.zeros(co.shape, dtype=int)
        for i, c in enumerate(co):
            focality[i] = np.sum(weights[norm > c * peak_value])

        return focality

    def summary(self, percentiles=(99.9, 99, 95), focality_cutoffs=(75, 50), units=None):
        ''' Creates a text summaty of the field

        Parameters
        ------------
        percentiles: ndarray (optinal)
            Field percentiles to be printed. Default: (99.9, 99, 95)
        focality_cutoffs: ndarray (optional)
            Cuttofs for focality calculations. Default: (75, 50)
        units: str or None
            Name of field units or automatically determine from name
        '''
        if units is None:
            if self.field_name in ['E', 'normE', 'D', 'g']:
                units = 'V/m'
            elif self.field_name in ['J', 'normJ']:
                units = 'A/m²'
            elif self.field_name == 'v':
                units = 'V'
            else:
                units = ''
        if units:
            units = ' ' + units

        if 2 in self.mesh.elm.elm_type:
            units_mesh = ' mm²'
        elif 4 in self.mesh.elm.elm_type:
            units_mesh = ' mm³'
        if np.all(np.isin([2, 4], self.mesh.elm.elm_type)):
            warnings.warn('Field summary in meshes with'
                          ' triangles and tetrahedra can give misleading results')
        norm = self._norm()
        weights = self._weights()
        mean_norm = np.sum(norm * weights)/np.sum(weights)
        prc = self.get_percentiles(percentiles)
        focality = self.get_focality(focality_cutoffs, percentiles[-1])
        string = f'Field: {self.field_name}\n'
        string += 'Peak Values:\n'
        n_spaces = len(f'{percentiles[0]:.2e}{units} ') - 5
        string += (n_spaces*' ' + '|').join(f'{p:.1f}%' for p in percentiles) + '\n'
        string += ' |'.join(f'{p:.2e}{units}' for p in prc) + '\n'
        string += 'Focality:\n'
        n_spaces = len(f'{focality[0]:.2e}{units_mesh} ') - 5
        string += (n_spaces*' ' + '|').join(f'{fc:.1f}%' for fc in focality_cutoffs) + '\n'
        string += ' |'.join(f'{fv:.2e}{units_mesh}' for fv in focality) + '\n'

        string += f'Mean Field:\n{mean_norm:.2f}{units}'
        return string

    def view_options(self, v_range='auto',
                     percentile=False,
                     visible=True,
                     visible_tags=None,
                     saturate=True,
                     idx=None):
        ''' Generates a View object with visualization opions

        Parameters
        -----------
        v_range: [min, max] or 'auto' (optional)
            Range of the values to be displayed. Defaut: will automatically try to find a
            good range (0.1 - 99.9% percentile for positive or vector data, symetrical
            for two-sided data)
        percentile: bool (optional)
            Whether the values in v_range arre given in percentile (between 0 and 100).
            Default: False
        visible: bool (optional)
            Whether to turn on visualization of this field. Default: True
        visible_tags: list (optional)
            List of tags to be visible. Will also be used to calulate field ranges in
            'auto' mode or if percentile=True. Default: all tags
        saturate: bool (optional)
            Whether to saturate values. Default: true
        idx: int (optional)
            Index of this field in the mesh
        Returns
        ----------
        view: gmsh_visualization.View
            view object
        '''
        self._test_msh()
        view = gmsh_view.View(indx=idx)
        if isinstance(self, NodeData):
            view.GlyphLocation = 2
            if visible_tags is not None:
                roi = self.mesh.elm.nodes_with_tag(visible_tags)
            else:
                roi = None
        elif isinstance(self, ElementData):
            view.GlyphLocation = 1
            if visible_tags is not None:
                roi = np.isin(self.mesh.elm.tag1, visible_tags)
            else:
                roi = None
        else:
            raise NotImplementedError

        if self.nr_comp == 1 and self.value.min() * self.value.max() < 0:
            two_sided = True
        else:
            two_sided = False
        if v_range == 'auto':
            if two_sided:
                prc = self.get_percentiles([0.1, 99.9], roi)
                max_ = np.max(np.abs(prc))
                min_ = -max_
                view.ColormapNumber = 10

            else:
                # all positive or vector
                if self.nr_comp > 1 or self.value.min() > -1e-10:
                    min_ = 0
                    max_ = self.get_percentiles(99.9, roi)[0]

                # All negative
                else:
                    max_ = 0
                    min_ = self.get_percentiles(0.1, roi)[0]

            # if the min and the max are close together (eg. masks)
            if np.isclose(min_, max_, atol=1e-12):
                min_ = self.value.min()
                max_ = self.value.max()

        else:
            if percentile is True:
                min_, max_ = self.get_percentiles([min_, max_], roi)
            else:
                min_, max_ = v_range

        view.CustomMax = float(max_)
        view.CustomMin = float(min_)
        view.RangeType = 2
        view.Visible = int(visible)
        view.SaturateValues = int(saturate)

        return view

    @property
    def indexing_nr(self):
        ''' Nodes or element numbers '''
        raise Exception('indexing_nr is not defined')

    def __eq__(self, other):
        try:
            return self.__dict__ == other.__dict__
        except AttributeError:
            return False

    def __getitem__(self, index):
        return _getitem_one_indexed(self.value, index)

    def __setitem__(self, index, item):
        index = _fix_indexing_one(index)
        self.value[index] = item

    def __mul__(self, other):
        cp = copy.copy(self)
        cp.value = self.value.__mul__(other)
        return cp

    def __neg__(self):
        cp = copy.copy(self)
        cp.value = self.value.__neg__()
        return cp

    def __sub__(self, other):
        cp = copy.copy(self)
        cp.value = self.value.__sub__(other)
        return cp

    def __add__(self, other):
        cp = copy.copy(self)
        cp.value = self.value.__add__(other)
        return cp

    def __str__(self):
        return self.field_name + '\n' + self.value.__str__()

    def __truediv__(self, other):
        cp = copy.copy(self)
        cp.value = self.value.__truediv__(other)
        return cp

    def __div__(self, other):
        cp = copy.copy(self)
        cp.value = self.value.__div__(other)
        return cp

    def __pow__(self, p):
        cp = copy.copy(self)
        cp.value = self.value.__pow__(p)
        return cp

    def _test_msh(self):
        if self.mesh is None:
            raise ValueError('Cannot evaluate function if .mesh property is not '
                             'assigned')

    def write_hdf5(self, hdf5_fn, path='./'):
        ''' Writes the field to an hdf5 file

        Parameters
        -----------
        hdf5_fn: str
            file name of hdf5 file
        path: str
            path in the hdf5 file where the field should be saved

        '''
        with h5py.File(hdf5_fn, 'a') as f:
            try:
                g = f.create_group(path)
            except ValueError:
                g = f[path]
            try:
                g.create_dataset(self.field_name, data=self.value)
                # If the dataset already exists, overwrite it
            except RuntimeError:
                del g[self.field_name]
                g.create_dataset(self.field_name, data=self.value)

    @classmethod
    def read_hdf5_data_matrix_row(cls, leadfield_fn, field_name, row):
        """ Reads a row of an hdf5 data matrix and store it as Data

        Parameters
        -----------
        leadfield_fn: str
            Name of file with leadfield
        field_name: str
            name of field
        row: int
            number of the row to be read

        Returns
        ---------
        data: Data()
            instance with the fields
        """
        import h5py
        with h5py.File(leadfield_fn, 'r') as f:
            value = f[field_name][row]
        return cls(value, field_name)

class ElementData(Data):
    """
    Data (scalar, vector or tensor) defined in mesh elements.

    Parameters
    -----------------------
    value: ndarray
        Value of field in elements. Should have the shape
            - (n_elm,) or (n_elm, 1) for scalar fields
            - (n_elm, 3) for vector fields
            - (n_elm, 9) for tensors

    field_name: str (optional)
        name of field. Default: ''

    mesh: simnibs.msh.Msh (optional)
        Mesh structure where the field is defined. Required for many methods


    Attributes
    --------------------------
    value: ndarray
        Value of field in elements
    field_name: str
        name of field
    elm_number: ndarray
        index of elements
    nr: property
        number of data points
    nr_comp: property
        number of dimensions per data point (1 for scalars, 3 for vectors)
    """

    def __init__(self, value, name='', mesh=None):
        Data.__init__(self, value=value, name=name, mesh=mesh)

    @property
    def elm_number(self):
        '''Element numbers (1, ..., nr)'''
        return np.arange(1, self.nr + 1, dtype='int32')

    @property
    def indexing_nr(self):
        '''Same as elm_number'''
        return self.elm_number

    def elm_data2node_data(self):
        """Transforms an ElementData field into a NodeData field using Superconvergent
        patch recovery
        For volumetric data. Will not work well for discontinuous fields (like E, if
        several tissues are used)

        Returns
        ----------------------
        simnibs.NodeData
            Structure with NodeData

        References
        -------------------
            Zienkiewicz, Olgierd Cecil, and Jian
            Zhong Zhu. "The superconvergent patch recovery and a posteriori error
            estimates. Part 1: The recovery technique." International Journal for
            Numerical Methods in Engineering 33.7 (1992): 1331-1364.

        """
        self._test_msh()
        msh = self.mesh
        if self.nr != msh.elm.nr:
            raise ValueError("The number of data points in the data "
                             "structure should be equal to the number of elements in the mesh")
        nd = np.zeros((msh.nodes.nr, self.nr_comp))

        if len(msh.elm.tetrahedra) == 0:
            raise ValueError("Can only transform volume data")

        value = np.atleast_2d(self[msh.elm.tetrahedra])
        if value.shape[0] < value.shape[1]:
            value = value.T

        # get all nodes used in tetrahedra, creates the NodeData structure
        uq = np.unique(msh.elm[msh.elm.tetrahedra])
        nd = NodeData(np.zeros((len(uq), self.nr_comp)), self.field_name, mesh=msh)

        # Get the point in the outside surface
        points_outside = np.unique(msh.elm.get_outside_faces())
        outside_points_mask = np.in1d(msh.elm[msh.elm.tetrahedra],
                                      points_outside).reshape(-1, 4)
        masked_th_nodes = np.copy(msh.elm[msh.elm.tetrahedra])
        masked_th_nodes[outside_points_mask] = -1

        # Calculates the quantities needed for the superconvergent patch recovery
        uq_in, th_nodes = np.unique(masked_th_nodes, return_inverse=True)
        baricenters = msh.elements_baricenters()[msh.elm.tetrahedra]
        volumes = msh.elements_volumes_and_areas()[msh.elm.tetrahedra]
        baricenters = np.hstack([np.ones((baricenters.shape[0], 1)), baricenters])

        A = np.empty((len(uq_in), 4, 4))
        b = np.empty((len(uq_in), 4, self.nr_comp), self.value.dtype)
        for i in range(4):
            for j in range(i, 4):
                A[:, i, j] = np.bincount(th_nodes.reshape(-1),
                                         np.repeat(baricenters[:, i], 4) *
                                         np.repeat(baricenters[:, j], 4))
        A[:, 1, 0] = A[:, 0, 1]
        A[:, 2, 0] = A[:, 0, 2]
        A[:, 3, 0] = A[:, 0, 3]
        A[:, 2, 1] = A[:, 1, 2]
        A[:, 3, 1] = A[:, 1, 3]
        A[:, 3, 2] = A[:, 2, 3]

        for j in range(self.nr_comp):
            for i in range(4):
                b[:, i, j] = np.bincount(th_nodes.reshape(-1),
                                         np.repeat(baricenters[:, i], 4) *
                                         np.repeat(value[:, j], 4))

        a = np.linalg.solve(A[uq_in != -1], b[uq_in != -1])
        p = np.hstack([np.ones((np.sum(uq_in != -1), 1)), msh.nodes[uq_in[uq_in != -1]]])
        f = np.einsum('ij, ijk -> ik', p, a)
        nd[uq_in[uq_in != -1]] = f

        # Assigns the average value to the points in the outside surface
        masked_th_nodes = np.copy(msh.elm[msh.elm.tetrahedra])
        masked_th_nodes[~outside_points_mask] = -1
        uq_out, th_nodes_out = np.unique(masked_th_nodes,
                                         return_inverse=True)
        sum_vals = np.empty((len(uq_out), self.nr_comp), self.value.dtype)
        for j in range(self.nr_comp):
            sum_vals[:, j] = np.bincount(th_nodes_out.reshape(-1),
                                         np.repeat(value[:, j], 4) *
                                         np.repeat(volumes, 4))

        sum_vols = np.bincount(th_nodes_out.reshape(-1), np.repeat(volumes, 4))

        nd[uq_out[uq_out != -1]] = (sum_vals/sum_vols[:, None])[uq_out != -1]

        nd.value = np.squeeze(nd.value)
        return nd

    def as_nodedata(self):
        ''' Converts the current ElementData instance to NodaData
        For more information see the elm_data2node_data method
        '''
        return self.elm_data2node_data()

    def interpolate_scattered(self, points, out_fill=np.nan, method='linear',
                              continuous=False, squeeze=True):
        ''' Interpolates the ElementData into the points by finding the element
        containing the point and assigning the value in it

        Parameters
        -------
        points: Nx3 ndarray
            List of points where we want to interpolate
        out_fill: float
            Value to be goven to points outside the volume, if 'nearest' assigns the
            nearest value (default: NaN)
        method: {'assign' or 'linear'} (Optional)
            If 'assign', gives to each voxel the value of the element that contains
            it. If linear, first assign fields to nodes, and then perform
            baricentric interpolatiom. Default: linear
        continuous: bool
            Wether fields is continuous across tissue boundaries. Changes the
            behaviour of the function only if method == 'linear'. Default: False
        squeeze: bool
            Wether to squeeze the output. Default: True

        Returns:
        -------
        f: np.ndarray
            Value of function in the points
        '''
        self._test_msh()
        msh = self.mesh
        if len(self.value.shape) > 1:
            f = np.zeros((points.shape[0], self.nr_comp), self.value.dtype)
        else:
            f = np.zeros((points.shape[0], ), self.value.dtype)
        th_with_points = \
            msh.find_tetrahedron_with_points(points, compute_baricentric=False)
        inside = th_with_points != -1

        if method == 'assign':
            f[inside] = self[th_with_points[inside]]
            if out_fill == 'nearest':
                _, nearest = msh.find_closest_element(
                    points[~inside], return_index=True)

                f[~inside] = self[nearest]

            else:
                f[~inside] = out_fill

        elif method == 'linear':
            if continuous:
                nd = self.elm_data2node_data()
                f = nd.interpolate_scattered(points, out_fill=out_fill, squeeze=False)
            else:
                # if all points are outside
                if not np.any(inside):
                    tags = [np.unique(msh.elm.tag1[msh.elm.elm_type==4])]
                # if there are points inside
                else:
                    tags = np.unique(msh.elm.tag1[th_with_points[inside] - 1])
                msh_copy = copy.deepcopy(msh)
                msh_copy.elmdata = [ElementData(self.value, mesh=msh_copy)]
                # create a list of fields at each tag
                field_at_tags = []
                for t in tags:
                    msh_tag = msh_copy.crop_mesh(tags=t)
                    nd = msh_tag.elmdata[0].elm_data2node_data()
                    field_at_tags.append(nd.interpolate_scattered(points, out_fill=np.nan, squeeze=False))
                    del msh_tag
                    del nd
                # Join the list of field values
                del msh_copy
                f = field_at_tags[0]
                count = np.ones(len(f), dtype=int)
                for f_t in field_at_tags[1:]:
                    # find where f and f_t are unasigned
                    unasigned_f = np.isnan(f)
                    if unasigned_f.ndim == 2:
                        unasigned_f = np.any(unasigned_f, axis=1)
                    unasigned_ft = np.isnan(f_t)
                    if unasigned_ft.ndim == 2:
                        unasigned_ft = np.any(unasigned_ft, axis=1)
                    # Assign to f unassigned values
                    f[unasigned_f] = f_t[unasigned_f]
                    # if for some reason a value is in 2 tissues, calculate the average
                    f[~unasigned_f * ~unasigned_ft] += f_t[~unasigned_f * ~unasigned_ft]
                    count[~unasigned_f * ~unasigned_ft] += 1
                del field_at_tags
                if f.ndim == 2:
                    f /= count[:, None]
                else:
                    f /= count

                # Finally, fill in the unassigned values
                unasigned_f = np.isnan(f)
                if unasigned_f.ndim == 2:
                    unasigned_f = np.any(unasigned_f, axis=1)
                if out_fill == 'nearest':
                    _, nearest = msh.find_closest_element(points[unasigned_f],
                                                          return_index=True)
                    f[unasigned_f] = self[nearest]
                else:
                    f[unasigned_f] = out_fill
        else:
            raise ValueError('Invalid interpolation method!')

        if squeeze:
            f = np.squeeze(f)
        return f

    def interpolate_to_grid(self, n_voxels, affine, method='linear', continuous=False):
        ''' Interpolates the ElementData into a grid.
            finds which tetrahedra contais the given voxel and
            assign the value of the tetrahedra to the voxel.

        Parameters
        -----------
        n_voxels: list or tuple
            number of voxels in x, y, and z directions
        affine: ndarray
            A 4x4 matrix specifying the transformation from voxels to xyz
        method: {'assign' or 'linear'} (Optional)
            If 'assign', gives to each voxel the value of the element that contains
            it. If linear, first assign fields to nodes, and then perform
            baricentric interpolatiom. Default: linear
        continuous: bool
            Wether fields is continuous across tissue boundaries. Changes the
            behaviour of the function only if method == 'linear'. Default: False

        Returns
        --------
        image: ndarray
            An (n_voxels[0], n_voxels[1], n_voxels[2], nr_comp) matrix with
            interpolated values. If nr_comp == 1, the last dimension is squeezed out
        '''

        msh = self.mesh
        self._test_msh()
        if self.nr != msh.elm.nr:
            raise ValueError('Invalid Mesh! Mesh should have the same number of elements'
                             'as the number of data points')
        if len(n_voxels) != 3:
            raise ValueError('n_voxels should have length = 3')
        if affine.shape != (4, 4):
            raise ValueError('Affine should be a 4x4 matrix')

        msh_th = msh.crop_mesh(elm_type=4)
        msh_th.elmdata = []
        msh_th.nodedata = []
        v = np.atleast_2d(self.value)
        if v.shape[0] < v.shape[1]:
            v = v.T
        v = v[msh.elm.elm_type == 4]

        if method == 'assign':
            nd = np.hstack([msh_th.nodes.node_coord, np.ones((msh_th.nodes.nr, 1))])
            inv_affine = np.linalg.inv(affine)
            nd = inv_affine.dot(nd.T).T[:, :3]

            # initialize image
            image = np.zeros([n_voxels[0], n_voxels[1], n_voxels[2], self.nr_comp], dtype=np.float)
            field = v.astype(np.float)
            image = cython_msh.interp_grid(
                np.array(n_voxels, dtype=np.int), field, nd.astype(np.float),
                (msh_th.elm.node_number_list - 1).astype(np.int))
            image = image.astype(self.value.dtype)
            if self.nr_comp == 1:
                image = np.squeeze(image, axis=3)

        elif method == 'linear':
            if continuous:
                nd = self.elm_data2node_data()
                image = nd.interpolate_to_grid(n_voxels, affine)

            else:
                if self.nr_comp != 1:
                    image = np.zeros(list(n_voxels) + [self.nr_comp], dtype=float)
                else:
                    image = np.zeros(list(n_voxels), dtype=float)
                # Interpolate each tag separetelly
                tags = np.unique(msh_th.elm.tag1)
                msh_th.elmdata = [ElementData(v, mesh=msh_th)]
                for t in tags:
                    msh_tag = msh_th.crop_mesh(tags=t)
                    nd = msh_tag.elmdata[0].elm_data2node_data()
                    image += nd.interpolate_to_grid(n_voxels, affine)
                    del msh_tag
                    del nd
                    gc.collect()
        else:
            raise ValueError('Invalid interpolation method!')

        del msh_th
        del v

        gc.collect()
        return image

    def assign_triangle_values(self):
        ''' In-place Assigns field value at triangle as the same as the one of the tetrahedra with
        a similar tag to it.

        '''
        self._test_msh()
        msh = self.mesh
        corrensponding = msh.find_corresponding_tetrahedra()
        self[msh.elm.triangles[corrensponding >= 0]] =\
                self[corrensponding[corrensponding >= 0]]

    def calc_flux(self, triangles=None):
        ''' Calculates the flux of a vectorial field

        Parameters
        -----------
        triangles: ndarray of ints (optional)
            Triangles where the flux should be calculated. Default: All

        Returns
        -------
        flux: float
            total field crossing the surface

        '''
        self._test_msh()
        msh = self.mesh
        if triangles is None:
            triangles = msh.elm.triangles
        normals = msh.triangle_normals()
        areas = msh.elements_volumes_and_areas()
        flux = np.sum(areas[triangles] *
                      np.sum(normals[triangles] * self[triangles], axis=1))
        return flux

    def norm(self, ord=2):
        ''' Calculate the norm of the field

        Parameters
        ------------
        ord: float
            Order of norm. Default: 2 (euclidian norm)

        Returns
        -----------
        norm: NodeData
            NodeData field with the norm the field
        '''
        if len(self.value.shape) == 1:
            ed = ElementData(np.abs(self.value),
                             name='norm' + self.field_name,
                             mesh=self.mesh)
        else:
            ed = ElementData(np.linalg.norm(self.value, axis=1, ord=ord),
                             name='norm' + self.field_name,
                             mesh=self.mesh)
        return ed



    @classmethod
    def from_data_grid(cls, mesh, data_grid, affine, field_name='', **kwargs):
        ''' Defines an ElementData field form a mesh and gridded data

        Parameters
        ---------
        mesh: Msh()
            Mesh structure where the field is to be interpolated
        data_grid: ndarray
            Array of 3 or 4 dimensions with data
        affine: 4x4 ndarray
            Array describing the affine transformation from the data grid to the mesh
            space
        kwargs: see the scipy.ndimage.map_coordinates documentation
        '''
        assert len(data_grid.shape) in [3, 4], \
                'The data grid must have 3 or 4 dimensions'
        bar = mesh.elements_baricenters().value.T
        iM = np.linalg.inv(affine)
        coords = iM[:3, :3].dot(bar) + iM[:3, 3, None]
        f = partial(
            scipy.ndimage.map_coordinates, coordinates=coords,
            output=data_grid.dtype, **kwargs)
        if len(data_grid.shape) == 4:
            indim = data_grid.shape[3]
            outdata = np.array(
                [f(data_grid[..., i]) for i in range(indim)]).T
        elif len(data_grid.shape) == 3:
            outdata = f(data_grid)

        ed = cls(outdata, name=field_name, mesh=mesh)
        return ed

    def append_to_mesh(self, fn, mode='binary'):
        """Appends this ElementData fields to a file

        Parameters
        ----------
            fn: str
                file name
            mode: binary or ascii
                mode in which to write
        """
        with open(fn, 'ab') as f:
            f.write(b'$ElementData\n')
            # string tags
            f.write((str(1) + '\n').encode('ascii'))
            f.write(('"' + self.field_name + '"\n').encode('ascii'))

            f.write((str(1) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))

            f.write((str(4) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))
            f.write((str(self.nr_comp) + '\n').encode('ascii'))
            f.write((str(self.nr) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))

            if mode == 'ascii':
                for ii in range(self.nr):
                    f.write((str(self.elm_number[ii]) + ' ' +
                             str(self.value[ii]).translate(None, '[](),') +
                             '\n').encode('ascii'))

            elif mode == 'binary':

                elm_number = self.elm_number.astype('int32')
                value = self.value.astype('float64')
                try:
                    value.shape[1]
                except IndexError:
                    value = value[:, np.newaxis]
                m = elm_number[:, np.newaxis]
                for i in range(self.nr_comp):
                    m = np.concatenate((m,
                                        value[:, i].astype('float64').view('int32').reshape(-1, 2)),
                                       axis=1)
                f.write(m.tostring())

            else:
                raise IOError("invalid mode:", mode)

            f.write(b'$EndElementData\n')

    def write(self, fn):
        """Writes this ElementData fields to a file with field information only
        This file needs to be merged with a mesh for visualization

        Parameters
        ---------
            fn: str
                file name
            mode: binary or ascii
                mode in which to write
        """
        with open(fn, 'wb') as f:
            f.write(b'$MeshFormat\n2.2 1 8\n')
            f.write(struct.pack('i', 1))
            f.write(b'\n$EndMeshFormat\n')
            f.write(b'$ElementData\n')
            # string tags
            f.write((str(1) + '\n').encode('ascii'))
            f.write(('"' + self.field_name + '"\n').encode('ascii'))

            f.write((str(1) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))

            f.write((str(4) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))
            f.write((str(self.nr_comp) + '\n').encode('ascii'))
            f.write((str(self.nr) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))

            elm_number = self.elm_number.astype('int32')
            value = self.value.astype('float64')
            try:
                value.shape[1]
            except IndexError:
                value = value[:, np.newaxis]
            m = elm_number[:, np.newaxis]
            for i in range(self.nr_comp):
                m = np.concatenate((m,
                                    value[:, i].astype('float64').view('int32').reshape(-1, 2)),
                                   axis=1)
            f.write(m.tostring())

            f.write(b'$EndElementData\n')



class NodeData(Data):
    """
    Data (scalar, vector or tensor) defined in mesh nodes.

    Parameters
    ----------
    value: ndarray
        Value of field in nodes. Should have the shape
            - (n_nodes,) or (n_nodes, 1) for scalar fields
            - (n_nodes, 3) for vector fields
            - (n_nodes, 9) for tensors

    field_name: str (optional)
        name of field. Default: ''

    mesh: simnibs.msh.Msh (optinal)
        Mesh where the field is defined. Required for many methods

    Attributes
    --------------------------
    value: ndarray
        Value of field in elements
    field_name: str
        name of field
    node_number: ndarray
        index of elements
    nr: property
        number of data points
    nr_comp: property
        number of dimensions per data point (1 for scalars, 3 for vectors)
    """

    def __init__(self, value, name='', mesh=None):
        Data.__init__(self, value=value, name=name, mesh=mesh)

    @property
    def node_number(self):
        ''' Node numbers (1, ..., nr)'''
        return np.arange(1, self.nr + 1, dtype='int32')

    @property
    def indexing_nr(self):
        ''' Same as node_numbers'''
        return self.node_number

    def as_nodedata(self):
        return self

    def node_data2elm_data(self):
        """Transforms an ElementData field into a NodeData field
        the value in the element is the average of the value in the nodes

        Returns
        --------
        simnibs.ElementData
            structure with field value interpolated at element centers

        """
        self._test_msh()
        msh = self.mesh
        if (self.nr != msh.nodes.nr):
            raise ValueError(
                "The number of data points in the data structure should be"
                "equal to the number of elements in the mesh")

        triangles = np.where(msh.elm.elm_type == 2)[0]
        tetrahedra = np.where(msh.elm.elm_type == 4)[0]

        if self.nr_comp == 1:
            elm_data = np.zeros((msh.elm.nr,), dtype=float)
        else:
            elm_data = np.zeros((msh.elm.nr, self.nr_comp), dtype=float)

        if len(triangles) > 0:
            elm_data[triangles] = \
                np.average(self.value[msh.elm.node_number_list[
                           triangles, :3] - 1], axis=1)
        if len(tetrahedra) > 0:
            elm_data[tetrahedra] = \
                np.average(self.value[msh.elm.node_number_list[
                           tetrahedra, :4] - 1], axis=1)

        return ElementData(elm_data, self.field_name, mesh=msh)

    def gradient(self):
        ''' Calculates the gradient of a field in the middle of the tetrahedra

        Parameters
        ------
        mesh: simnibs.Msh
            A mesh with the geometrical information

        Returns
        -----
        grad: simnibs.ElementData
            An ElementData field with gradient in the middle of each tetrahedra
        '''
        self._test_msh()
        mesh = self.mesh
        if self.nr_comp != 1:
            raise ValueError('can only take gradient of scalar fields')
        if mesh.nodes.nr != self.nr:
            raise ValueError('mesh must have the same number of nodes as the NodeData')

        # Tetrahedra gradients
        elm_node_coords = mesh.nodes[mesh.elm[mesh.elm.tetrahedra]]

        tetra_matrices = elm_node_coords[:, 1:4, :] - \
            elm_node_coords[:, 0, :][:, None]

        dif_between_tetra_nodes = self[mesh.elm[mesh.elm.tetrahedra][:, 1:4]] - \
            self[mesh.elm[mesh.elm.tetrahedra][:, 0]][:, None]

        th_grad = np.linalg.solve(tetra_matrices, dif_between_tetra_nodes)

        gradient = np.zeros((mesh.elm.nr, 3), dtype=float)
        gradient[mesh.elm.elm_type == 4] = th_grad
        gradient[mesh.elm.elm_type == 2] = 0.

        return ElementData(gradient, 'grad_' + self.field_name, mesh)

    def calc_flux(self, nodes=None):
        ''' Calculates the flux of a vector field though the given nodes

        Parameters
        -----------
        nodes: list of ints
            List of node indices where to calculate flux. Default: all nodes in triangle
            surface

        Returns
        --------
        flux: float
            Total fux through all surfcaces
        '''
        self._test_msh()
        msh = self.mesh
        if msh.nodes.nr != self.nr:
            raise ValueError('Mesh should have {0} nodes'.format(self.nr))
        if nodes is None:
            nodes = np.unique(msh.elm[msh.elm.triangles][:, :3])
        assert self.nr_comp == 3, 'Can only calculate flux of vectors'
        normals = msh.nodes_normals()
        areas = msh.nodes_areas()
        flux = np.sum(
            np.sum(normals[nodes] * self[nodes], axis=1) * areas[nodes])
        return flux

    def interpolate_scattered(self, points, out_fill=np.nan, squeeze=True):
        ''' Interpolates the NodeaData into the points by finding the element
        containing the point and performing linear interpolation inside the element

        Parameters
        ------
        points: Nx3 ndarray
            List of points where we want to interpolate
        out_fill: float
            Value to be goven to points outside the volume. If 'nearest', will assign the
            value of th nearest node. (default: NaN)
        squeeze: bool
            Wether to squeeze the output. Default: True
        Returns:
        ----
        f: np.ndarray
            Value of function in the points
        '''
        self._test_msh()
        msh = self.mesh
        if len(self.value.shape) > 1:
            f = np.zeros((points.shape[0], self.nr_comp), self.value.dtype)
        else:
            f = np.zeros((points.shape[0], ), self.value.dtype)

        th_with_points, bar = \
            msh.find_tetrahedron_with_points(points, compute_baricentric=True)
        inside = th_with_points != -1
        if np.any(inside) and len(self.value.shape) == 1:
            f[inside] = np.einsum('ik, ik -> i',
                                  self[msh.elm[th_with_points[inside]]],
                                  bar[inside])
        elif np.any(inside):
            f[inside] = np.einsum('ikj, ik -> ij',
                                  self[msh.elm[th_with_points[inside]]],
                                  bar[inside])

        if out_fill == 'nearest':
            _, nearest = msh.nodes.find_closest_node(points[~inside],
                                                     return_index=True)
            f[~inside] = self[nearest]
        else:
            f[~inside] = out_fill

        if squeeze:
            f = np.squeeze(f)

        return f

    def interpolate_to_grid(self, n_voxels, affine, **kwargs):
        ''' Interpolates the NodeData into a grid.
            finds which tetrahedra contais the given voxel and
            performs linear interpolation inside the voxel

        The kwargs is ony to have the same interface as the ElementData version

        Parameters
        ------
            n_voxels: list or tuple
                number of voxels in x, y, and z directions
            affine: ndarray
                A 4x4 matrix specifying the transformation from voxels to xyz
        Returns
        ----
            image: ndarray
                An (n_voxels[0], n_voxels[1], n_voxels[2], nr_comp) matrix with
                interpolated values. If nr_comp == 1, the last dimension is squeezed out
        '''

        msh = self.mesh
        self._test_msh()
        if len(n_voxels) != 3:
            raise ValueError('n_voxels should have length = 3')
        if affine.shape != (4, 4):
            raise ValueError('Affine should be a 4x4 matrix')

        v = np.atleast_2d(self.value)
        if v.shape[0] < v.shape[1]:
            v = v.T
        msh_th = msh.crop_mesh(elm_type=4)
        inv_affine = np.linalg.inv(affine)
        nd = np.hstack([msh_th.nodes.node_coord, np.ones((msh_th.nodes.nr, 1))])
        nd = inv_affine.dot(nd.T).T[:, :3]

        # initialize image
        image = np.zeros([n_voxels[0], n_voxels[1], n_voxels[2], self.nr_comp], dtype=float)
        field = v.astype(float)
        if v.shape[0] != msh_th.nodes.nr:
            raise ValueError('Number of data points in the structure does not match '
                             'the number of nodes present in the volume-only mesh')
        image = cython_msh.interp_grid(
            np.array(n_voxels, dtype=np.int), field, nd,
            msh_th.elm.node_number_list - 1)
        image = image.astype(self.value.dtype)
        if self.nr_comp == 1:
            image = np.squeeze(image, axis=3)
        del nd
        del msh_th
        del field
        gc.collect()
        return image

    def norm(self, ord=2):
        ''' Calculate the norm of the field

        Parameters
        ------------
        ord: float 
            Order of norm. Default: 2 (euclidian norm)

        Returns
        ---------
        norm: NodeData
            NodeData field with the norm the field
        '''
        if len(self.value.shape) == 1:
            nd = NodeData(np.abs(self.value),
                          name='norm' + self.field_name,
                          mesh=self.mesh)
        else:
            nd = NodeData(np.linalg.norm(self.value, axis=1, ord=ord),
                          name='norm' + self.field_name,
                          mesh=self.mesh)
        return nd

    def normal(self, fill=np.nan):
        ''' Calculate the normal component of the field in the mesh surfaces

        Parameters
        -----------
        fill: float (optional)
            Value to be used when node is not in surface (Default: NaN)

        Returns
        ---------
        normal: NodeData
            NodeData field with the normal the field where a surface is defined and the
            fill value where it's not
        '''
        self._test_msh()
        assert self.nr_comp == 3, 'Normals are only defined for vector fields'
        normal = NodeData(fill * np.ones(self.nr, dtype=self.value.dtype),
                          name='normal' + self.field_name,
                          mesh=self.mesh)
        nodes_in_surface = self.mesh.elm[self.mesh.elm.triangles, :3]
        nodes_in_surface = np.unique(nodes_in_surface)
        nodes_normals = self.mesh.nodes_normals()
        normal[nodes_in_surface] = np.sum(self[nodes_in_surface] * nodes_normals[nodes_in_surface], axis=1)
        return normal

    def angle(self, fill=np.nan):
        ''' Calculate the angle between the field and the surface normal
        
        Parameters
        -------------
        fill: float (optional)
            Value to be used when node is not in surface

        Returns
        -----
        angle: NodeData
            NodeData field with the angles the field where a surface is defined and the
            fill value where it's not
        '''
        self._test_msh()
        assert self.nr_comp == 3, 'angles are only defined for vector fields'
        angle = NodeData(fill * np.ones(self.nr, dtype=self.value.dtype),
                          name='angle' + self.field_name,
                          mesh=self.mesh)

        nodes_in_surface = self.mesh.elm[self.mesh.elm.triangles, :3]
        nodes_in_surface = np.unique(nodes_in_surface)
        nodes_normals = self.mesh.nodes_normals()
        normal = np.sum(self[nodes_in_surface] * nodes_normals[nodes_in_surface], axis=1)
        norm = np.linalg.norm(self[nodes_in_surface], axis=1)
        tan = np.sqrt(norm ** 2 - normal ** 2)
        #angle[nodes_in_surface] = np.arccos(normal/norm)
        angle[nodes_in_surface] = np.arctan2(tan, normal)
        return angle


    def tangent(self, fill=np.nan):
        ''' Calculate the tangent component of the field in the surfaces
        
        Parameters
        -----------
        fill: float (optional)
            Value to be used when node is not in surface

        Returns
        -----
        tangent: NodeData
            NodeData field with the tangent component of the field where a surface is defined and the
            fill value where it's not
        '''
        self._test_msh()
        assert self.nr_comp == 3, 'angles are only defined for vector fields'
        tangent = NodeData(fill * np.ones(self.nr, dtype=self.value.dtype),
                          name='tangent' + self.field_name,
                          mesh=self.mesh)
        nodes_in_surface = self.mesh.elm[self.mesh.elm.triangles, :3]
        nodes_in_surface = np.unique(nodes_in_surface)
        nodes_normals = self.mesh.nodes_normals()
        normal = np.sum(self[nodes_in_surface] * nodes_normals[nodes_in_surface], axis=1)
        norm = np.linalg.norm(self[nodes_in_surface], axis=1)
        tangent[nodes_in_surface] = np.sqrt(norm ** 2 - normal ** 2)
        return tangent



    def append_to_mesh(self, fn, mode='binary'):
        """Appends this NodeData fields to a file

        Parameters
        ------------
            fn: str
                file name
            mode: binary or ascii
                mode in which to write
        """
        with open(fn, 'ab') as f:
            f.write(b'$NodeData\n')
            # string tags
            f.write((str(1) + '\n').encode('ascii'))
            f.write(('"' + self.field_name + '"\n').encode('ascii'))

            f.write((str(1) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))

            f.write((str(4) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))
            f.write((str(self.nr_comp) + '\n').encode('ascii'))
            f.write((str(self.nr) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))

            if mode == 'ascii':
                for ii in range(self.nr):
                    f.write((
                        str(self.node_number[ii]) + ' ' +
                        str(self.value[ii]).translate(None, '[](),') +
                        '\n').encode('ascii'))

            elif mode == 'binary':
                value = self.value.astype('float64')
                try:
                    value.shape[1]
                except IndexError:
                    value = value[:, np.newaxis]

                m = self.node_number[:, np.newaxis].astype('int32')
                for i in range(self.nr_comp):
                    m = np.concatenate((m,
                                        value[:, i].astype('float64').view('int32').reshape(-1, 2)),
                                       axis=1)

                f.write(m.tostring())
            else:
                raise IOError("invalid mode:", mode)

            f.write(b'$EndNodeData\n')


    def write(self, fn):
        """Writes this NodeData field to a file with field information only
        This file needs to be merged with a mesh for visualization

        Parameters
        -----------
            fn: str
                file name
        """
        with open(fn, 'wb') as f:
            f.write(b'$MeshFormat\n2.2 1 8\n')
            f.write(struct.pack('i', 1))
            f.write(b'\n$EndMeshFormat\n')
            f.write(b'$NodeData\n')
            f.write((str(1) + '\n').encode('ascii'))
            f.write(('"' + self.field_name + '"\n').encode('ascii'))

            f.write((str(1) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))

            f.write((str(4) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))
            f.write((str(self.nr_comp) + '\n').encode('ascii'))
            f.write((str(self.nr) + '\n').encode('ascii'))
            f.write((str(0) + '\n').encode('ascii'))

            value = self.value.astype('float64')
            try:
                value.shape[1]
            except IndexError:
                value = value[:, np.newaxis]

            m = self.node_number[:, np.newaxis].astype('int32')
            for i in range(self.nr_comp):
                m = np.concatenate((m,
                                    value[:, i].astype('float64').view('int32').reshape(-1, 2)),
                                   axis=1)

            f.write(m.tostring())

            f.write(b'$EndNodeData\n')


def read_msh(fn, m=None):
    ''' Reads a gmsh '.msh' file

    Parameters
    ------------
    fn: str
        File name
    m: simnibs.msh.Msh (optional)
        Mesh structure to be overwritten. If unset, will create a new structure

    Returns
    --------
    msh: simnibs.msh.Msh
        Mesh structure
    '''
    if m is None:
        m = Msh()

    fn = os.path.expanduser(fn)

    if not os.path.isfile(fn):
        raise IOError(fn + ' not found')

    version_number = _find_mesh_version(fn)
    if version_number == 2:
        m = _read_msh_2(fn, m)

    elif version_number == 4:
        m = _read_msh_4(fn, m)

    else:
        raise IOError('Unrecgnized Mesh file version : {}'.format(version_number))

    return m


def _find_mesh_version(fn):
    if not os.path.isfile(fn):
        raise IOError(fn + ' not found')

    # file open
    with open(fn, 'rb') as f:
        # check 1st line
        first_line = f.readline().decode()
        if first_line != '$MeshFormat\n':
            raise IOError(fn, "must start with $MeshFormat")

        # parse 2nd line
        version_number, file_type, data_size = f.readline().decode().split()
        version_number = int(version_number[0])
        file_type = int(file_type)
        data_size = int(data_size)
    return version_number


def _read_msh_2(fn, m):
    m.fn = fn

    # file open
    with open(fn, 'rb') as f:
        # check 1st line
        first_line = f.readline()
        if first_line != b'$MeshFormat\n':
            raise IOError(fn, "must start with $MeshFormat")

        # parse 2nd line
        version_number, file_type, data_size = f.readline().decode().strip().split()
        version_number = int(version_number[0])
        file_type = int(file_type)
        data_size = int(data_size)

        if version_number != 2:
            raise IOError("Can only handle v2 meshes")

        if file_type == 1:
            binary = True
        elif file_type == 0:
            binary = False
        else:
            raise IOError("File_type not recognized: {0}".format(file_type))

        if data_size != 8:
            raise IOError(
                "data_size should be double (8), i'm reading: {0}".format(data_size))

        # read next byte, if binary, to check for endianness
        if binary:
            endianness = struct.unpack('i', f.readline()[:4])[0]
        else:
            endianness = 1

        if endianness != 1:
            raise IOError("endianness is not 1, is the endian order wrong?")

        # read 3rd line
        if f.readline() != b'$EndMeshFormat\n':
            raise IOError(fn + " expected $EndMeshFormat")

        # read 4th line
        if f.readline() != b'$Nodes\n':
            raise IOError(fn + " expected $Nodes")

        # read 5th line with number of nodes
        try:
            node_nr = int(f.readline().decode().strip())
        except:
            raise IOError(fn + " something wrong with Line 5 - should be a number")

        # read all nodes
        if binary:
            # 0.02s to read binary.msh
            dt = np.dtype([
                ('id', np.int32, 1),
                ('coord', np.float64, 3)])

            temp = np.fromfile(f, dtype=dt, count=node_nr)
            node_number = np.copy(temp['id'])
            node_coord = np.copy(temp['coord'])

            # sometimes there's a line feed here, sometimes there is not...
            LF_byte = f.read(1)  # read \n
            if not ord(LF_byte) == 10:
                # if there was not a LF, go back 1 byte from the current file
                # position
                f.seek(-1, 1)

        else:
            # nodes has 4 entries: [node_ID x y z]
            node_number = np.empty(node_nr, dtype='int32')
            # array Nx3 for (x,y,z) coordinates of the nodes
            node_coord = np.empty(3 * node_nr, dtype='float64')
            for ii in range(node_nr):
                line = f.readline().decode().strip().split()
                node_number[ii] = line[0]
                # it's faster to use a linear array and than reshape
                node_coord[3 * ii] = line[1]
                node_coord[3 * ii + 1] = line[2]
                node_coord[3 * ii + 2] = line[3]
            node_coord = node_coord.reshape((node_nr, 3))

        if not np.all(node_number == np.arange(1, node_nr + 1)):
            warnings.warn("Mesh file with discontinuos nodes, things can fail"
                          " unexpectedly")
        m.nodes.node_coord = node_coord

        if f.readline() != b'$EndNodes\n':
            raise IOError(fn + " expected $EndNodes after reading " +
                          str(node_nr) + " nodes")

        # read all elements
        if f.readline() != b'$Elements\n':
            raise IOError(fn, "expected line with $Elements")

        try:
            elm_nr = int(f.readline().decode().strip())
        except:
            raise IOError(
                fn + " something wrong when reading number of elements (line after $Elements)"
                "- should be a number")

        if binary:
            current_element = 0

            elm_number = np.empty(elm_nr, dtype='int32')
            m.elm.elm_type = np.empty(elm_nr, dtype='int32')
            m.elm.tag1 = np.empty(elm_nr, dtype='int32')
            m.elm.tag2 = np.empty(elm_nr, dtype='int32')
            m.elm.node_number_list = -np.ones((elm_nr, 4), dtype='int32')
            read = np.ones(elm_nr, dtype=bool)

            nr_nodes_elm = [None, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9,
                            10, 27, 18, 14, 1, 8, 20, 15, 13]
            while current_element < elm_nr:
                elm_type, nr, _ = np.fromfile(f, 'int32', 3)
                if elm_type == 2:
                    tmp = np.fromfile(f, 'int32', nr * 6).reshape(-1, 6)

                    m.elm.elm_type[current_element:current_element+nr] = \
                        2 * np.ones(nr, 'int32')
                    elm_number[current_element:current_element+nr] = tmp[:, 0]
                    m.elm.tag1[current_element:current_element+nr] = tmp[:, 1]
                    m.elm.tag2[current_element:current_element+nr] = tmp[:, 2]
                    m.elm.node_number_list[current_element:current_element+nr, :3] = tmp[:, 3:]
                    read[current_element:current_element+nr] = 1

                elif elm_type == 4:
                    tmp = np.fromfile(f, 'int32', nr * 7).reshape(-1, 7)

                    m.elm.elm_type[current_element:current_element+nr] = \
                        4 * np.ones(nr, 'int32')
                    elm_number[current_element:current_element+nr] = tmp[:, 0]
                    m.elm.tag1[current_element:current_element+nr] = tmp[:, 1]
                    m.elm.tag2[current_element:current_element+nr] = tmp[:, 2]
                    m.elm.node_number_list[current_element:current_element+nr] = tmp[:, 3:]
                    read[current_element:current_element+nr] = 1

                else:
                    warnings.warn('element of type {0} '
                                  'cannot be read, ignoring it'.format(elm_type))
                    np.fromfile(f, 'int32', nr * (3 + nr_nodes_elm[elm_type]))
                    read[current_element:current_element+nr] = 0
                current_element += nr

            elm_number = elm_number[read]
            m.elm.elm_type = m.elm.elm_type[read]
            m.elm.tag1 = m.elm.tag1[read]
            m.elm.tag2 = m.elm.tag2[read]
            m.elm.node_number_list = m.elm.node_number_list[read]

        else:

            elm_number = np.empty(elm_nr, dtype='int32')
            m.elm.elm_type = np.empty(elm_nr, dtype='int32')
            m.elm.tag1 = np.empty(elm_nr, dtype='int32')
            m.elm.tag2 = np.empty(elm_nr, dtype='int32')
            m.elm.node_number_list = -np.ones((elm_nr, 4), dtype='int32')
            read = np.ones(elm_nr, dtype=bool)

            for ii in range(elm_nr):
                line = f.readline().decode().strip().split()
                if line[1] == '2':
                    elm_number[ii] = line[0]
                    m.elm.elm_type[ii] = line[1]
                    m.elm.tag1[ii] = line[3]
                    m.elm.tag2[ii] = line[4]
                    m.elm.node_number_list[ii, :3] = [int(i) for i in line[5:]]
                elif line[1] == '4':
                    elm_number[ii] = line[0]
                    m.elm.elm_type[ii] = line[1]
                    m.elm.tag1[ii] = line[3]
                    m.elm.tag2[ii] = line[4]
                    m.elm.node_number_list[ii] = [int(i) for i in line[5:]]
                else:
                    read[ii] = 0
                    warnings.warn('element of type {0} '
                                  'cannot be read, ignoring it'.format(elm_type))

            elm_number = elm_number[read]
            m.elm.elm_type = m.elm.elm_type[read]
            m.elm.tag1 = m.elm.tag1[read]
            m.elm.tag2 = m.elm.tag2[read]
            m.elm.node_number_list = m.elm.node_number_list[read]

        elm_nr_changed = False
        if not np.all(elm_number == np.arange(1, m.elm.nr + 1)):
            warnings.warn('Changing element numbering')
            elm_nr_changed = True
            elm_number = np.arange(1, m.elm.nr + 1)

        line = f.readline()
        if b'$EndElements' not in line:
            line = f.readline()
            if b'$EndElements' not in line:
                raise IOError(fn + " expected $EndElements after reading " +
                              str(m.elm.nr) + " elements. Read " + line)

        # read the header in the beginning of a data section
        def parse_Data():
            section = f.readline()
            if section == b'':
                return 'EOF', '', 0, 0
            # string tags
            number_of_string_tags = int(f.readline().decode('ascii'))
            assert number_of_string_tags == 1, "Invalid Mesh File: invalid number of string tags"
            name = f.readline().decode('ascii').strip().strip('"')
            # real tags
            number_of_real_tags = int(f.readline().decode('ascii'))
            assert number_of_real_tags == 1, "Invalid Mesh File: invalid number of real tags"
            f.readline()
            # integer tags
            number_of_integer_tags = int(f.readline().decode().strip())  # usually 3 or 4
            integer_tags = [int(f.readline().decode().strip())
                            for i in range(number_of_integer_tags)]
            nr = integer_tags[2]
            nr_comp = integer_tags[1]
            return section.strip(), name, nr, nr_comp

        def read_NodeData(t, name, nr, nr_comp, m):
            data = NodeData(np.empty((nr, nr_comp)), name=name, mesh=m)
            if binary:
                dt = np.dtype([
                    ('id', np.int32, 1),
                    ('values', np.float64, nr_comp)])

                temp = np.fromfile(f, dtype=dt, count=nr)
                node_number = np.copy(temp['id'])
                data.value = np.copy(temp['values'])
            else:
                node_number = np.empty(nr, dtype='int32')
                data.value = np.empty((nr, nr_comp), dtype='float64')
                for ii in range(nr):
                    line = f.readline().decode().strip().split()
                    node_number[ii] = int(line[0])
                    data.value[ii, :] = [float(v) for v in line[1:]]

            if f.readline() != b'$EndNodeData\n':
                raise IOError(fn + " expected $EndNodeData after reading " +
                              str(nr) + " lines in $NodeData")

            if np.any(node_number != m.nodes.node_number):
                raise IOError("Can't read NodeData field: "
                              "it does not have one data point per node")

            return data

        def read_ElementData(t, name, nr, nr_comp, m):
            if elm_nr_changed or not np.all(read):
                raise IOError('Could not read ElementData: '
                              'Element ordering not compact or invalid element type')
            data = ElementData(np.empty((nr, nr_comp)), name=name, mesh=m)
            if binary:
                dt = np.dtype([
                    ('id', np.int32, 1),
                    ('values', np.float64, nr_comp)])

                temp = np.fromfile(f, dtype=dt, count=nr)
                elm_number = np.copy(temp['id'])
                data.value = np.copy(temp['values'])

            else:
                elm_number = np.empty(nr, dtype='int32')
                data.value = np.empty([nr, nr_comp], dtype='float64')

                for ii in range(nr):
                    line = f.readline().decode().strip().split()
                    elm_number[ii] = int(line[0])
                    data.value[ii, :] = [float(jj) for jj in line[1:]]

            if f.readline() != b'$EndElementData\n':
                raise IOError(fn + " expected $EndElementData after reading " +
                              str(nr) + " lines in $ElementData")

            if np.any(elm_number != m.elm.elm_number):
                raise IOError("Can't read ElementData field: "
                              "it does not have one data point per element")

            return data

        # read sections recursively
        def read_next_section():
            t, name, nr, nr_comp = parse_Data()
            if t == 'EOF':
                return
            elif t == b'$NodeData':
                m.nodedata.append(read_NodeData(t, name, nr, nr_comp, m))
            elif t == b'$ElementData':
                m.elmdata.append(read_ElementData(t, name, nr, nr_comp, m))
            else:
                raise IOError("Can't recognize section name:" + t)

            read_next_section()
            return

        read_next_section()
    m.compact_ordering(node_number)
    return m


def _read_msh_4(fn, m):
    m.fn = fn
    # file open
    with open(fn, 'rb') as f:
        # check 1st line
        first_line = f.readline()
        if first_line != b'$MeshFormat\n':
            raise IOError(fn, "must start with $MeshFormat")

        # parse 2nd line
        version_number, file_type, data_size = f.readline().decode().strip().split()
        version_number = int(version_number[0])
        file_type = int(file_type)
        data_size = int(data_size)

        if version_number != 4:
            raise IOError("Can only handle v4 meshes")

        if file_type == 1:
            binary = True
        elif file_type == 0:
            binary = False
        else:
            raise IOError("File_type not recognized: {0}".format(file_type))

        if data_size != 8:
            raise IOError(
                "data_size should be double (8), i'm reading: {0}".format(data_size))

        # read next byte, if binary, to check for endianness
        if binary:
            endianness = struct.unpack('i', f.readline()[:4])[0]
        else:
            endianness = 1

        if endianness != 1:
            raise IOError("endianness is not 1, is the endian order wrong?")

        # read 3rd line
        if f.readline() != b'$EndMeshFormat\n':
            raise IOError(fn + " expected $EndMeshFormat")

        # Skip  everyting untill nodes
        while f.readline() != b'$Nodes\n':
            continue

        # Number of nodes and number of blocks
        if binary:
            entity_blocks, node_nr = struct.unpack('LL', f.read(struct.calcsize('LL')))
        else:
            line = f.readline().strip()
            entity_blocks, node_nr = line.decode().split()
            entity_blocks = int(entity_blocks)
            node_nr = int(node_nr)
        n_read = 0
        node_number = np.zeros(node_nr, dtype=np.int32)
        node_coord = np.zeros((node_nr, 3), dtype=np.float)
        for block in range(entity_blocks):
            if binary:
                _, _, parametric, n_in_block = struct.unpack(
                    'iiii', f.read(struct.calcsize('iiii')))
                # We need to read 4 extra bytes here
                f.read(4)
            else:
                _, _, parametric, n_in_block = f.readline().decode().strip().split()
                parametric = int(parametric)
                n_in_block = int(n_in_block)
            if parametric:
                raise IOError("Can't read parametric entity!")
            if binary:
                dt = np.dtype([
                    ('id', np.int32, 1),
                    ('coord', np.float64, 3)])
                temp = np.fromfile(f, dtype=dt, count=n_in_block)
                node_nbr_block = temp['id']
                node_coord_block = temp['coord']
            else:
                node_nbr_block = np.zeros(n_in_block, dtype=np.int)
                node_coord_block = np.zeros(3 * n_in_block, dtype=np.float)
                for i in range(n_in_block):
                    line = f.readline().decode().strip().split()
                    node_nbr_block[i] = line[0]
                    node_coord_block[3*i] = line[1]
                    node_coord_block[3*i + 1] = line[2]
                    node_coord_block[3*i + 2] = line[3]
                node_coord_block = node_coord_block.reshape(-1, 3)

            node_number[n_read:n_read+n_in_block] = node_nbr_block
            node_coord[n_read:n_read+n_in_block, :] = node_coord_block
            n_read += n_in_block

        m.nodes.node_coord = node_coord

        line = f.readline()
        if b'$EndNodes' not in line:
            line = f.readline()
            if b'$EndNodes' not in line:
                raise IOError(fn + " expected $EndNodes after reading " +
                              str(m.noedes.nr) + " nodes. Read " + line)

        # Read Elements
        if f.readline() != b'$Elements\n':
            raise IOError(fn, "expected line with $Elements")
        if binary:
            entity_blocks, elm_nr = struct.unpack('LL', f.read(struct.calcsize('LL')))
        else:
            line = f.readline().strip()
            entity_blocks, elm_nr = line.decode().split()
            entity_blocks = int(entity_blocks)
            elm_nr = int(elm_nr)
        n_read = 0
        elm_number = np.zeros(elm_nr, dtype=np.int32)
        m.elm.elm_type = np.zeros(elm_nr, dtype=np.int32)
        m.elm.tag1 = np.zeros(elm_nr, dtype=np.int32)
        m.elm.node_number_list = -np.ones((elm_nr, 4), dtype=np.int32)
        read = np.ones(elm_nr, dtype=bool)
        for block in range(entity_blocks):
            if binary:
                tag, _, elm_type, n_in_block = struct.unpack(
                    'iiii', f.read(struct.calcsize('iiii')))
                f.read(4)
            else:
                tag, _, elm_type, n_in_block = f.readline().decode().strip().split()
                tag = int(tag)
                elm_type = int(elm_type)
                n_in_block = int(n_in_block)

            if elm_type == 2:
                nr_nodes_elm = 3
            elif elm_type == 4:
                nr_nodes_elm = 4
            else:
                warnings.warn(
                    "Can't read element type: {}. Ignoring it".format(elm_type))
                continue

            if binary:
                dt = np.dtype([
                    ('id', np.int32, 1),
                    ('nodes', np.int32, nr_nodes_elm)])
                temp = np.fromfile(f, dtype=dt, count=n_in_block)
                elm_nbr_block = temp['id']
                elm_node_block = temp['nodes']

            else:
                elm_nbr_block = np.zeros(n_in_block, dtype=np.int32)
                elm_node_block = -np.ones(nr_nodes_elm * n_in_block, dtype=np.int32)
                for i in range(n_in_block):
                    line = f.readline().decode().strip().split()
                    elm_nbr_block[i] = int(line[0])
                    elm_node_block[nr_nodes_elm*i:nr_nodes_elm*(i+1)] = [int(l) for l in line[1:]]
                elm_node_block = elm_node_block.reshape(-1, nr_nodes_elm)

            elm_number[n_read:n_read+n_in_block] = elm_nbr_block
            m.elm.node_number_list[n_read:n_read+n_in_block, :nr_nodes_elm] = elm_node_block
            m.elm.tag1[n_read:n_read+n_in_block] = tag
            m.elm.elm_type[n_read:n_read+n_in_block] = elm_type
            read[n_read:n_read+n_in_block] = True
            n_read += n_in_block

        elm_number = elm_number[read]
        m.elm.node_number_list = m.elm.node_number_list[read]
        m.elm.tag1 = m.elm.tag1[read]
        m.elm.elm_type = m.elm.elm_type[read]

        order = np.argsort(elm_number)
        elm_number = elm_number[order]
        m.elm.node_number_list = m.elm.node_number_list[order]
        m.elm.tag1 = m.elm.tag1[order]
        m.elm.elm_type = m.elm.elm_type[order]
        m.elm.tag2 = m.elm.tag1

        line = f.readline()
        if b'$EndElements' not in line:
            line = f.readline()
            if b'$EndElements' not in line:
                raise IOError(fn + " expected $EndElements after reading " +
                              str(m.elm.nr) + " elements. Read " + line)

        # read the header in the beginning of a data section
        def parse_Data():
            section = f.readline()
            if section == b'':
                return 'EOF', '', 0, 0
            # string tags
            number_of_string_tags = int(f.readline().decode('ascii'))
            assert number_of_string_tags == 1, "Invalid Mesh File: invalid number of string tags"
            name = f.readline().decode('ascii').strip().strip('"')
            # real tags
            number_of_real_tags = int(f.readline().decode('ascii'))
            assert number_of_real_tags == 1, "Invalid Mesh File: invalid number of real tags"
            f.readline()
            # integer tags
            number_of_integer_tags = int(f.readline().decode('ascii'))  # usually 3 or 4
            integer_tags = [int(f.readline().decode('ascii'))
                            for i in range(number_of_integer_tags)]
            nr = integer_tags[2]
            nr_comp = integer_tags[1]
            return section.strip(), name, nr, nr_comp

        def read_NodeData(t, name, nr, nr_comp, m):
            data = NodeData(np.empty((nr, nr_comp)), name=name, mesh=m)
            if binary:
                dt = np.dtype([
                    ('id', np.int32, 1),
                    ('values', np.float64, nr_comp)])

                temp = np.fromfile(f, dtype=dt, count=nr)
                node_number = np.copy(temp['id'])
                data.value = np.copy(temp['values'])
            else:
                node_number = np.empty(nr, dtype='int32')
                data.value = np.empty((nr, nr_comp), dtype='float64')
                for ii in range(nr):
                    line = f.readline().decode('ascii').split()
                    node_number[ii] = int(line[0])
                    data.value[ii, :] = [float(v) for v in line[1:]]

            if f.readline() != b'$EndNodeData\n':
                raise IOError(fn + " expected $EndNodeData after reading " +
                              str(nr) + " lines in $NodeData")

            if np.any(node_number != m.nodes.elm_number):
                raise IOError("Can't read NodeData field: "
                              "it does not have one data point per node")

            return data

        def read_ElementData(t, name, nr, nr_comp, m):
            if not np.all(read):
                raise IOError('Could not read ElementData: '
                              'Element ordering not compact or invalid element type')
            data = ElementData(np.empty((nr, nr_comp)), name=name, mesh=m)
            if binary:
                dt = np.dtype([
                    ('id', np.int32, 1),
                    ('values', np.float64, nr_comp)])

                temp = np.fromfile(f, dtype=dt, count=nr)
                elm_number = np.copy(temp['id'])
                data.value = np.copy(temp['values'])

            else:
                elm_number = np.empty(nr, dtype='int32')
                data.value = np.empty([nr, nr_comp], dtype='float64')

                for ii in range(nr):
                    line = f.readline().decode('ascii').split()
                    elm_number[ii] = int(line[0])
                    data.value[ii, :] = [float(jj) for jj in line[1:]]

            if f.readline() != b'$EndElementData\n':
                raise IOError(fn + " expected $EndElementData after reading " +
                              str(nr) + " lines in $ElementData")

            if np.any(elm_number != m.elm.elm_number):
                raise IOError("Can't read ElementData field: "
                              "it does not have one data point per element")

            return data

        # read sections recursively
        def read_next_section():
            t, name, nr, nr_comp = parse_Data()
            if t == 'EOF':
                return
            elif t == b'$NodeData':
                m.nodedata.append(read_NodeData(t, name, nr, nr_comp, m))
            elif t == b'$ElementData':
                m.elmdata.append(read_ElementData(t, name, nr, nr_comp, m))
            else:
                raise IOError("Can't recognize section name:" + t)

            read_next_section()
            return

        read_next_section()
    m.compact_ordering(node_number)
    return m


# write msh to mesh file
def write_msh(msh, file_name=None, mode='binary'):
    ''' Writes a gmsh 'msh' file

    Parameters
    ------------
    msh: simnibs.msh.Msh
        Mesh structure
    file_name: str (optional)
        Name of file to be writte. Default: msh.fn
    mode: 'binary' or 'ascii':
        The mode in which the file should be read
    '''
    if file_name is not None:
        msh.fn = file_name

    fn = msh.fn

    if fn[0] == '~':
        fn = os.path.expanduser(fn)

    if mode not in ['ascii', 'binary']:
        raise ValueError("Only 'ascii' and 'binary' are allowed")

    with open(fn, 'wb') as f:
        if mode == 'ascii':
            f.write(b'$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')

        elif mode == 'binary':
            f.write(b'$MeshFormat\n2.2 1 8\n')
            f.write(struct.pack('i', 1))
            f.write(b'\n$EndMeshFormat\n')

        # write nodes
        f.write(b'$Nodes\n')
        f.write('{0}\n'.format(msh.nodes.nr).encode('ascii'))

        if mode == 'ascii':
            for ii in range(msh.nodes.nr):
                f.write((str(msh.nodes.node_number[ii]) + ' ' +
                         str(msh.nodes.node_coord[ii][0]) + ' ' +
                         str(msh.nodes.node_coord[ii][1]) + ' ' +
                         str(msh.nodes.node_coord[ii][2]) + '\n').encode('ascii'))

        elif mode == 'binary':
            node_number = msh.nodes.node_number.astype('int32')
            node_coord = msh.nodes.node_coord.astype('float64')
            m = node_number[:, np.newaxis]
            for i in range(3):
                nc = node_coord[:, i].astype('float64')
                m = np.concatenate((m,
                                    nc.view('int32').reshape(-1, 2)), axis=1)
            f.write(m.tostring())
        f.write(b'$EndNodes\n')
        # write elements
        f.write(b'$Elements\n')
        f.write((str(msh.elm.nr) + '\n').encode('ascii'))

        if mode == 'ascii':
            for ii in range(msh.elm.nr):
                line = str(msh.elm.elm_number[ii]) + ' ' + \
                    str(msh.elm.elm_type[ii]) + ' ' + str(2) + ' ' +\
                    str(msh.elm.tag1[ii]) + ' ' + str(msh.elm.tag2[ii]) + ' '

                if msh.elm.elm_type[ii] == 2:
                    line += str(msh.elm.node_number_list[ii, :3]
                                ).translate(None, '[](),') + '\n'
                elif msh.elm.elm_type[ii] == 4:
                    line += str(msh.elm.node_number_list[ii, :]
                                ).translate(None, '[](),') + '\n'
                else:
                    raise IOError(
                        "ERROR: cant write meshes with elements of type",
                        msh.elm.elm_type[ii])

                f.write(line.encode('ascii'))

        elif mode == 'binary':
            triangles = np.where(msh.elm.elm_type == 2)[0]
            if len(triangles > 0):
                triangles_header = np.array((2, len(triangles), 2), 'int32')
                triangles_number = msh.elm.elm_number[triangles].astype('int32')
                triangles_tag1 = msh.elm.tag1[triangles].astype('int32')
                triangles_tag2 = msh.elm.tag2[triangles].astype('int32')
                triangles_node_list = msh.elm.node_number_list[
                    triangles, :3].astype('int32')
                f.write(triangles_header.tostring())
                f.write(np.concatenate((triangles_number[:, np.newaxis],
                                        triangles_tag1[:, np.newaxis],
                                        triangles_tag2[:, np.newaxis],
                                        triangles_node_list), axis=1).tostring())

            tetra = np.where(msh.elm.elm_type == 4)[0]
            if len(tetra > 0):
                tetra_header = np.array((4, len(tetra), 2), 'int32')
                tetra_number = msh.elm.elm_number[tetra].astype('int32')
                tetra_tag1 = msh.elm.tag1[tetra].astype('int32')
                tetra_tag2 = msh.elm.tag2[tetra].astype('int32')
                tetra_node_list = msh.elm.node_number_list[tetra].astype('int32')

                f.write(tetra_header.tostring())
                f.write(np.concatenate((tetra_number[:, np.newaxis],
                                        tetra_tag1[:, np.newaxis],
                                        tetra_tag2[:, np.newaxis],
                                        tetra_node_list), axis=1).tostring())

        f.write(b'$EndElements\n')

    # write nodeData, if existent
    for nd in msh.nodedata:
        nd.append_to_mesh(fn, mode)

    for eD in msh.elmdata:
        eD.append_to_mesh(fn, mode)


'''
# Adds 1000 to the label of triangles, if less than 100
def create_surface_labels(msh):
    triangles = np.where(msh.elm.elm_type == 2)[0]
    triangles = np.where(msh.elm.tag1[triangles] < 1000)[0]
    msh.elm.tag1[triangles] += 1000
    msh.elm.tag2[triangles] += 1000
    return msh
'''

def read_res_file(fn_res, fn_pre=None):
    """ Reads a .res file

    Parameters
    ------------
    fn_res: str
        name of res file

    fn_pre: str
        name of pre file, if additional information is needed. Default: do not read .pre file

    Returns
    --------
    ndarray
        values
    """
    with open(fn_res, 'rb') as f:
        f.readline()  # skip first line
        check, type_of_file = f.readline().decode('ascii').strip('\n').split(' ')
        if check != '1.1':
            raise IOError('Unexpected value in res!')

        if type_of_file == '0':
            v = np.loadtxt(f, comments='$', skiprows=3, usecols=[
                           0], delimiter=' ', dtype=np.float)

        elif type_of_file == '1':
            f.readline()
            f.readline()
            f.readline()
            s = b''
            for line in f:
                if line == b'$EndSolution\n':
                    break
                else:
                    s += line
            s = s.strip()
            cols = np.frombuffer(s, dtype=np.float)
            v = cols[::2]

        else:
            raise IOError(
                'Do not recognize file type: %s for res file' % type_of_file)

    if fn_pre is not None:
        with open(fn_pre, 'br') as f:
            assert "Resolution" in f.readline().decode('ascii'), 'Invalid .pre file'
            main_resolution_number, number_of_dofdata = f.readline().decode('ascii')[:-1].split()
            main_resolution_number = int(main_resolution_number)
            number_of_dofdata = int(number_of_dofdata)
            assert f.readline() == b"$EndResolution\n", 'Invalid .pre file'
            assert b"$DofData" in f.readline(), 'Invalid .pre file'
            [f.readline() for i in range(4)]
            number_of_any_dof, number_of_dof = f.readline().decode('ascii')[:-1].split()
            number_of_any_dof = int(number_of_any_dof)
            number_of_dof = int(number_of_dof)

        dof_type = np.loadtxt(fn_pre, comments='$',
                              skiprows=9, usecols=[3], delimiter=' ',
                              dtype=int)

        dof_data = np.loadtxt(fn_pre, comments='$',
                              skiprows=9, usecols=[4], delimiter=' ',
                              dtype=float)

        vals = np.zeros(number_of_any_dof, dtype=float)
        vals[dof_type == 2] = dof_data[dof_type == 2]
        vals[dof_type == 1] = v
        v = vals

    return v


def write_geo_spheres(positions, fn, values=None, name="", mode='bw'):
    """ Writes a .geo file with spheres in specified positions

    Parameters
    ------------
    positions: nx3 ndarray:
        position of spheres
    fn: str
        name of file to be written
    values: nx1 ndarray  (optional)
        values to be assigned to the spheres. Default: 0
    name: str (optional)
        Name of the view
    mode: str (optional)
        Mode in which open the file. Default: 'bw'

    """
    if values is None:
        values = np.zeros((len(positions), ))

    if len(values) != len(positions):
        raise ValueError(
            'The length of the vector of positions is different from the'
            ' length of the vector of values')

    with open(fn, mode) as f:
        f.write(('View"' + name + '"{\n').encode('ascii'))
        for p, v in zip(positions, values):
            f.write(("SP(" + ", ".join([str(i) for i in p]) +
                     "){" + str(v) + "};\n").encode('ascii'))
        f.write(b"};\n")


def write_geo_text(positions, text, fn, name="", mode='bw'):
    """ Writes a .geo file with text in specified positions

    Parameters
    -----------
    positions: nx3 ndarray:
        position of spheres
    text: list of strings
        Text to be printed in each position
    fn: str
        name of file to be written
    mode: str
        Mode in which open the file. Default: 'bw'

    """
    if len(positions) != len(text):
        raise ValueError('The length of the vector of positions is different from the' +
                         'length of the list of text strings')

    with open(fn, mode) as f:
        f.write(('View"' + name + '"{\n').encode('ascii'))
        for p, t in zip(positions, text):
            f.write(("T3(" + ", ".join([str(i) for i in p]) +
                     ', TextAttributes("FontSize", "24")'
                     '){"' + t + '"};\n').encode('ascii'))
        f.write(b"};\n")


def write_geo_triangles(triangles, nodes, fn,  values=None, name="", mode="bw"):
    """ Writes a .geo file with triangles and field

    Parameters
    -----------
    triangles: tx3 ndarray
        Indices of the nodes composing each triangle
    nodes: nx3 ndarray:
        position of nodes
    fn: str
        name of file to be written
    values(optional): ndarray
        values to be assigned to the triangles. Default: 0
    name(optional): str
        Name of the view
    mode(oprional): str
        Mode in which open the file. Default: 'bw'
    """

    if values is None:
        values = np.zeros((len(triangles), 3))

    if values.ndim == 1:
        values = np.repeat(values[:, None], 3, axis=1)

    if len(values) != len(triangles):
        raise ValueError(
            'The length of the vector of triangles is different from the'
            ' length of the vector of values')

    with open(fn, mode) as f:
        f.write(('View"' + name + '"{\n').encode('ascii'))
        for t, v in zip(triangles, values):
            string = "ST("
            string += ", ".join([str(c) for c in nodes[t[0]]]) + ","
            string += ", ".join([str(c) for c in nodes[t[1]]]) + ","
            string += ", ".join([str(c) for c in nodes[t[2]]]) + ")"
            string += "{" + ",".join([str(vi) for vi in v]) + "};\n"
            f.write(string.encode('ascii'))
        f.write(b"};\n")



def read_freesurfer_surface(fn):
    ''' Function to read FreeSurfer surface files, based on the read_surf.m script by
    Bruce Fischl

    Parameters
    ------------
    fn: str
        File name

    Returns
    --------
    msh: Msh()
        Mesh structure

    '''
    TRIANGLE_FILE_MAGIC_NUMBER = 16777214
    QUAD_FILE_MAGIC_NUMBER = 16777215

    def read_3byte_integer(f):
        n = 0
        for i in range(3):
            b = f.read(1)
            n += struct.unpack("B", b)[0] << 8 * (2 - i)
        return n

    with open(fn, 'rb') as f:
        # Read magic number as a 3 byte integer
        magic = read_3byte_integer(f)
        if magic == QUAD_FILE_MAGIC_NUMBER:
            raise IOError('Quad files not supported!')

        elif magic == TRIANGLE_FILE_MAGIC_NUMBER:
            f.readline()
            f.readline()
            # notice that the format uses big-endian machine format
            vnum = struct.unpack(">i", f.read(4))[0]
            fnum = struct.unpack(">i", f.read(4))[0]
            vertex_coords = np.fromfile(f, np.dtype('>f'), vnum * 3)
            vertex_coords = vertex_coords.reshape((-1, 3)).astype(np.float64)
            faces = np.fromfile(f, np.dtype('>i'), fnum * 3)
            faces = faces.reshape((-1, 3)).astype(np.int64)
        else:
            raise IOError('Invalid magic number')

    msh = Msh()
    msh.elm = Elements(triangles=faces + 1)
    msh.nodes = Nodes(vertex_coords)
    return msh

def write_freesurfer_surface(msh, fn, ref_fs=None):
    ''' Writes a FreeSurfer surface
    Only the surfaces (triangles) are writen to the FreeSurfer surface file

    Parameters
    -----------
    msh: Msh()
        Mesh structure

    fn: str
        output file name

    ref_fs: str
        Name of ref_fs file, Used to set the tail
    '''
    def write_3byte_integer(f, n):
        b1 = struct.pack('B', (n >> 16) & 255)
        b2 = struct.pack('B', (n >> 8) & 255)
        b3 = struct.pack('B', n & 255)
        f.write(b1)
        f.write(b2)
        f.write(b3)

    TRIANGLE_FILE_MAGIC_NUMBER = 16777214
    m = msh.crop_mesh(elm_type=2)
    faces = m.elm.node_number_list[:, :3] - 1
    vertices = m.nodes.node_coord
    vnum = m.nodes.nr
    fnum = m.elm.nr

    if ref_fs is not None:
        ref_vol = nibabel.load(ref_fs)
        affine = ref_vol.affine.copy()[:3, :3]
        volume = ref_vol.header['dim'][1:4]
        voxelsize = np.sqrt(np.sum(affine ** 2, axis=0))
        affine /= voxelsize[None, :]
        affine = np.linalg.inv(affine)
        write_tail = True

    else:
        write_tail = False

    with open(fn, 'wb') as f:
        write_3byte_integer(f, TRIANGLE_FILE_MAGIC_NUMBER)
        f.write('Created by {0} on {1} with SimNIBS\n\n'.format(
            os.getenv('USER'), str(datetime.datetime.now())).encode('ascii'))
        f.write(struct.pack('>i', vnum))
        f.write(struct.pack('>i', fnum))
        f.write(vertices.reshape(-1).astype('>f').tostring())
        f.write(faces.reshape(-1).astype('>i').tostring())
        if write_tail:
            f.write(b'\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x14')
            f.write(b'valid = 1  # volume info valid\n')
            f.write('filename = {0}\n'.format(ref_fs).encode('ascii'))
            f.write('volume = {0:d} {1:d} {2:d}\n'.format(*volume).encode('ascii'))
            f.write('voxelsize = {0:.15e} {1:.15e} {2:.15e}\n'.format(*voxelsize).encode('ascii'))
            f.write('xras   = {0:.15e} {1:.15e} {2:.15e}\n'.format(*affine[0, :]).encode('ascii'))
            f.write('yras   = {0:.15e} {1:.15e} {2:.15e}\n'.format(*affine[1, :]).encode('ascii'))
            f.write('zras   = {0:.15e} {1:.15e} {2:.15e}\n'.format(*affine[2, :]).encode('ascii'))
            f.write('cras   = {0:.15e} {1:.15e} {2:.15e}\n'.format(0, 0, 0).encode('ascii'))


def read_gifti_surface(fn):
    ''' Reads a gifti surface

    Parameters
    -----------
    fn: str
        File name

    Returns
    ---------
    msh: Msh()
        mesh structure with geometrical information
    '''
    s = nibabel.load(fn)
    faces = s.get_arrays_from_intent('NIFTI_INTENT_TRIANGLE')[0].data
    nodes = s.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data
    msh = Msh()
    msh.elm = Elements(triangles=np.array(faces + 1, dtype=int))
    msh.nodes = Nodes(np.array(nodes, dtype=float))
    return msh


def read_curv(fn):
    ''' Reads a freesurfer .curv file

    Parameters
    ------------
    fn: str
        File name

    Returns
    ---------
    curv: np.ndarray
        array with file informatio
    '''
    NEW_VERSION_MAGIC_NUMBER = 16777215

    def read_3byte_integer(f):
        b = f.read(3)
        n = struct.unpack('>i', b'\x00' + b)[0]
        return n

    with open(fn, 'rb') as f:
        # Read magic number as a 3 byte integer
        magic = read_3byte_integer(f)
        if magic == NEW_VERSION_MAGIC_NUMBER:
            vnum = struct.unpack(">i", f.read(4))[0]
            fnum = struct.unpack(">i", f.read(4))[0]
            vals_per_vertex = struct.unpack(">i", f.read(4))[0]
            curv = np.fromfile(f, np.dtype('>f'), vnum)

        else:
            fnum = read_3byte_integer(f)
            curv = np.fromfile(f, np.dtype('>f'), magic) / 100.
    return curv


def write_curv(fn, curv, fnum):
    ''' Writes a freesurfer .curv file

    Parameters
    ------------
    fn: str
        File name to be written
    curv: ndaray
        Data array to be written
    fnum: int
        Number of faces in the mesh
    '''
    def write_3byte_integer(f, n):
        b1 = struct.pack('B', (n >> 16) & 255)
        b2 = struct.pack('B', (n >> 8) & 255)
        b3 = struct.pack('B', (n & 255))
        f.write(b1)
        f.write(b2)
        f.write(b3)


    NEW_VERSION_MAGIC_NUMBER = 16777215
    vnum = len(curv)
    with open(fn, 'wb') as f:
        write_3byte_integer(f, NEW_VERSION_MAGIC_NUMBER)
        f.write(struct.pack(">i", int(vnum)))
        f.write(struct.pack('>i', int(fnum)))
        f.write(struct.pack('>i', 1))
        f.write(curv.astype('>f').tostring())


def _middle_surface(wm_surface, gm_surface, depth):
    interp_surface = copy.deepcopy(wm_surface)
    interp_surface.nodes.node_coord = \
        depth * wm_surface.nodes.node_coord + \
        (1 - depth) * gm_surface.nodes.node_coord
    return interp_surface


def read_stl(fn):
    from .hmutils import mesh_load
    vertices, faces = mesh_load(fn)
    msh = Msh()
    msh.elm = Elements(triangles=faces + 1)
    msh.nodes = Nodes(vertices)

    return msh

def open_in_gmsh(fn, new_thread=False):
    ''' Opens the mesh in gmsh

    Parameters
    ------------
    fn: str
        Name of mesh file
    new_thread: bool
        Wether to open gmsh in a new thread. Defaut: False

    '''
    gmsh_bin = path2bin('gmsh')
    if new_thread:
        t = threading.Thread(target=subprocess.run,
                             args=([gmsh_bin, fn], ),
                             kwargs={'check': True})
        t.daemon = False  # thread dies with the program
        t.start()
    else:
        subprocess.run([gmsh_bin, fn], check=True)


def _hash_rows(array, mult=1000003, dtype=np.uint64):
    # Code based on python's tupleobject hasing
    # Generate hash
    mult = dtype(mult)
    sorted = np.sort(array, axis=1).astype(dtype)
    hash_array = 0x345678 * np.ones(array.shape[0], dtype=dtype)
    for i in reversed(range(array.shape[1])):
        hash_array = np.bitwise_xor(hash_array, sorted[:, i]) * mult
        mult += dtype(82520 + i + i)
    hash_array += dtype(97531)
    hash_array[hash_array == -1] = -2
    return hash_array



def _fix_indexing_one(index):
    '''Fix indexing to allow getting and setting items with one-idexed arrays'''
    def fix_slice(slc):
        start = slc.start
        stop = slc.stop
        step = slc.step
        if start is not None:
            if start > 0:
                start -= 1
            elif start == 0:
                raise IndexError('Cant get item 0 in one-indexed array')
            else:
                raise IndexError('Cant get negative slices in one-indexed array')

        if stop is not None:
            if stop > 0:
                stop -= 1
            elif stop == 0:
                raise IndexError('Cant get item 0 in one-indexed array')
            else:
                raise IndexError('Cant get negative slices in one-indexed array')

        if step is not None and step < 0:
            raise IndexError('Cant get negative slices in one-indexed array')

        return slice(start, stop, step)

    def fix_index_array(idx_array):
        idx_array = np.array(idx_array)
        if idx_array.dtype == bool:
            return idx_array
        else:
            if np.any(idx_array == 0):
                raise IndexError('Cant get item 0 in one-indexed array')
            idx_array[idx_array > 0] -= 1
            return idx_array

    def is_integer(index):
        answer = isinstance(index, int)
        answer += isinstance(index, np.int)
        answer += isinstance(index, np.intc)
        answer += isinstance(index, np.intp)
        answer += isinstance(index, np.int8)
        answer += isinstance(index, np.int16)
        answer += isinstance(index, np.int32)
        answer += isinstance(index, np.int64)
        answer += isinstance(index, np.uint8)
        answer += isinstance(index, np.uint16)
        answer += isinstance(index, np.uint32)
        answer += isinstance(index, np.uint64)
        return answer

    if is_integer(index):
        if index > 0:
            return index - 1
        if index < 0:
            return index
        if index == 0:
            raise IndexError('Cant get item 0 in one-indexed array')

    elif isinstance(index, slice):
        return fix_slice(index)

    elif isinstance(index, list) or isinstance(index, np.ndarray):
        return fix_index_array(index)

    elif isinstance(index, tuple):
        index = list(index)
        if isinstance(index[0], list) or isinstance(index[0], np.ndarray):
            index[0] = fix_index_array(index[0])
        if isinstance(index[0], slice):
            index[0] = fix_slice(index[0])
        return tuple(index)

    else:
        return index

def _getitem_one_indexed(array, index):
    index = _fix_indexing_one(index)
    return array.__getitem__(index)

class _GetitemTester():
    def __init__(self, array):
        self.array = array

    def __getitem__(self, index):
        return _getitem_one_indexed(self.array, index)

#!/usr/bin/python2.7 -u
# -*- coding: utf-8 -*-\
'''
   Structures for the GUI head model
   This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

   Copyright (C) 2015-2018 Guilherme B Saturnino

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
import gc


class Surface():

    def __init__(self, mesh, labels):
        """
        Class for working with mesh surfaces

        Parameters:
        ----------------
        mesh - gmsh mesh structure 
        labels - mesh labels with the surfaces of interest

        Notes:
        -------------------

        """

        # @var nodes
        # the nodes in the surface
        self.nodes = []

        # @var tr_nodes
        # The nodes in each triangle
        self. tr_nodes = []

        # @var tr_normals
        # the triangle normals
        self.tr_normals = []

        # @var tr_areas
        # Triangle areas
        self.tr_areas = []

        # @var tr_centers
        # The baricenter of each triangle
        self.tr_centers = []

        # @var node_normals
        # normal of a node. defined as the average of the average of the normals
        # that contain such node
        self.nodes_normals = []

        self.nodes_areas = []

        # @var values
        # field values
        self.values = []

        # @var surf2msh_triangles
        # list of mesh indices for surface triangles
        self.surf2msh_triangles = []

        # @var surf2msh_nodes
        # list of inidices for nodes
        self.surf2msh_nodes = []

        # nodes os a new array, of numpy arrays
        # it's filled up in the order the nodes appear in the mesh triangles
        # The nodes_dict makes the link between the gmsh numbering and the new
        # numbering
        nodes_dict = np.zeros(mesh.nodes.nr + 1, dtype='int')

        # Sees if labels is an array or just a number, it must be an array
        try:
            _ = labels[0]
        except TypeError:
            labels = [labels]

        # We need to create a dictionary, linking the "new node numbers" to the
        # old ones
        label_triangles = []
        self.surf2msh_triangles = np.empty(0, dtype='int')

        for label in labels:
            label_triangles.append(np.where(np.logical_and(
                mesh.elm.elm_type == 2, mesh.elm.tag1 == label))[0])

        if len(label_triangles) == 0:
            raise ValueError('No triangles with tags: ' + str(labels) + ' found')

        for ii in range(len(labels)):
            self.surf2msh_triangles = np.union1d(
                self.surf2msh_triangles, label_triangles[ii])

        # Creates a unique list of all triangle nodes.
        self.surf2msh_nodes = np.unique(
            mesh.elm.node_number_list[self.surf2msh_triangles, :3].reshape(-1))

        nr_unique = np.size(self.surf2msh_nodes)
        np_range = np.array(range(nr_unique), dtype='int')
        #nodes_dict[old_index] = new_index
        nodes_dict[self.surf2msh_nodes] = np_range

        # Gets the new node numbers
        self.tr_nodes = nodes_dict[
            mesh.elm.node_number_list[self.surf2msh_triangles, :3]]

        # and the positions in appropriate order
        self.nodes = mesh.nodes.node_coord[self.surf2msh_nodes - 1]
        self.nodes_normals = np.zeros(
            [np.size(self.surf2msh_nodes), 3], dtype='float')

        self.nodes_areas = np.zeros(
            [np.size(self.surf2msh_nodes), ], dtype='float')

        # positions of each triangle's node
        tr_pos = self.nodes[self.tr_nodes]

        # triangles normals and areas
        self.tr_normals = np.cross(
            tr_pos[:, 1] - tr_pos[:, 0], tr_pos[:, 2] - tr_pos[:, 0])
        self.tr_areas = 0.5 * np.linalg.norm(self.tr_normals, axis=1)
        self.tr_normals /= 2 * self.tr_areas[:, np.newaxis]

        # defining the normal of a node as being the average of the normals of
        # surrounding triangles
        for ii in range(len(self.surf2msh_triangles)):
            self.nodes_normals[self.tr_nodes[ii]
                               ] += self.tr_normals[ii][np.newaxis, :]
            self.nodes_areas[self.tr_nodes[ii]] += 1. / 3. * self.tr_areas[ii]

        self.nodes_normals /= np.linalg.norm(self.nodes_normals,
                                             axis=1)[:, np.newaxis]

        # out if the normals point inside or outside, calculates the center of
        # the mesh
        self.center = self.nodes.sum(axis=0) / self.nodes.shape[0]
        node2center = np.subtract(self.nodes, self.center)
        dotProducts = np.sum(node2center * self.nodes_normals, 1)
        # if the sum of the dot product is smaller than 0, the normals point
        # inwards. correct that
        if np.sum(dotProducts) < 0:
            self.nodes_normals *= -1
            self.tr_normals *= -1

        self.tr_centers = np.sum(self.nodes[self.tr_nodes], 1) / 3

        mesh = None
        gc.collect()

    def interceptRay(self, Near, Far, return_index=False):
        # Find point in surface in the near-far line thats nearest to the "near point"
        # based on http://geomalgorithms.com/a06-_intersect-2.html
        delta = 1e-9
        P1 = np.array(Near, float)
        P0 = np.array(Far, float)
        intersect_point = None
        intersect_normal = None
        triangle_index = None
        # nescessary when dealing with concave surfaces (like GM surface)

        V0 = self.nodes[self.tr_nodes[:, 0]]
        V1 = self.nodes[self.tr_nodes[:, 1]]
        V2 = self.nodes[self.tr_nodes[:, 2]]
        perp_project = np.sum(self.tr_normals * (P1 - P0)[None, :], axis=1)
        r = np.sum(self.tr_normals * np.subtract(V0, P0), 1)
        r = np.divide(r, perp_project)
        Pr = P0 + (P1 - P0) * r[:, np.newaxis]
        w = Pr - V0
        u = V1 - V0
        v = V2 - V0
        len2u = np.sum(u * u, 1)
        len2v = np.sum(v * v, 1)
        UdotV = np.sum(u * v, 1)
        WdotV = np.sum(w * v, 1)
        WdotU = np.sum(w * u, 1)
        denominator = (UdotV)**2 - (len2u) * (len2v)
        s = (UdotV * WdotV - len2v * WdotU) / denominator
        t = (UdotV * WdotU - len2u * WdotV) / denominator

        triangles = np.arange(len(self.tr_nodes), dtype='int')
        candidates = triangles[(perp_project > delta) * (s > -delta) * (
            t > -delta) * (s + t < 1 + delta) * (r > delta) * (r < 1 + delta)]

        if len(candidates) == 0:
            if return_index:
                return None, None, None
            else:
                return None, None

        else:
            # most of the time, there will only be one candidate
            triangle_index = candidates[np.argmin(r[candidates])]
            intersect_point = V0[triangle_index, :] + s[triangle_index] * \
                u[triangle_index, :] + t[triangle_index] * v[triangle_index, :]
            intersect_normal = self.tr_normals[triangle_index]
            if return_index:
                return triangle_index, intersect_point, intersect_normal
            else:
                return intersect_point, intersect_normal

    def projectPoint(self, point, reference=None, smooth=False):
        f = 1.0
        point_proj = None
        normal = None
        if reference is None:
            node_index = np.argmin(
                np.linalg.norm(point - self.nodes, axis=1))
            point_proj = self.nodes[node_index]
            normal = self.nodes_normals[node_index]
        while point_proj is None and f < 10:
            P1 = reference + f * (point - reference)
            P0 = reference
            triangle_index, point_proj, normal = self.interceptRay(
                P1, P0, True)
            f += 0.2

        return point_proj, normal

    def getTangentCoordinates(self, normal):
        # orthogonalization with the basis
        if normal[2] > 0.95:
            u = np.array([1, 0, 0], 'float')
            v = np.array([0, 1, 0], 'float')
            return u, v
        bases = []
        bases.append(np.array([1, 0, 0], 'float'))
        bases.append(np.array([0, 1, 0], 'float'))
        bases.append(np.array([0, 0, 1], 'float'))
        tangents = []
        # orthogonalize using the 2 smallest components
        dirs = [0, 1, 2]
        dirs = [x for (y, x) in sorted(zip(list(abs(normal)), dirs))]
        tangents.append(bases[dirs[0]])
        tangents[0] = tangents[0] - tangents[0].dot(normal) * normal
        tangents[0] = tangents[0] / np.linalg.norm(tangents[0])

        tangents.append(bases[dirs[1]])
        tangents[1] = tangents[1] - tangents[1].dot(normal) *\
            normal - tangents[1].dot(tangents[0]) * tangents[0]
            
        tangents[1] = tangents[1] / np.linalg.norm(tangents[1])

        # gets the standard vectors
        if abs(tangents[0][2]) >= abs(tangents[1][2]):
            u = -(tangents[1][2] / tangents[0][2]) * tangents[0] + tangents[1]
            u = u / np.linalg.norm(u)
            v = tangents[0] - tangents[0].dot(u) * u
            v = v / np.linalg.norm(v)

        if abs(tangents[0][2]) < abs(tangents[1][2]):
            u = tangents[0] - (tangents[0][2] / tangents[1][2]) * tangents[1]
            u = u / np.linalg.norm(u)
            v = tangents[1] - tangents[1].dot(u) * u
            v = v / np.linalg.norm(v)

        if v[2] < 0:
            v = -v
        if u[0] < 0:
            u = -u

        return u, v

    # Smooths the normal of the triangle using the normals of its nodes
    def smoothNormal(self, triangle_index):
        triangle = self.tr_nodes[triangle_index]
        new_normal = np.sum(self.nodes_normals[triangle], 0)
        new_normal /= np.linalg.norm(new_normal)

        return new_normal

    # Given 2 points, calculates a transformation matrix for a coordinate system which has:
    #self.center in p1
    # y axis going from p1 to p2
    def calculateMatSimnibs(self, p1, p2, skin_distance=4):
        p1_np = np.array(p1, 'float')
        p2_np = np.array(p2, 'float')
        center, z_axis = self.projectPoint(p1_np, smooth=True)
        z_axis = -1 * z_axis
        center = center - skin_distance * z_axis

        # y direction: user input, orthogonalized
        y_axis = p2_np - p1_np
        y_axis = y_axis - y_axis.dot(z_axis) * z_axis
        if np.linalg.norm(y_axis) > 1e-6:
            y_axis /= np.linalg.norm(y_axis)
        else:
            y_axis = np.random.rand(3)
            y_axis = y_axis - y_axis.dot(z_axis) * z_axis
            y_axis /= np.linalg.norm(y_axis)


        x_axis = np.cross(y_axis, z_axis)

        matsimnibs = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

        (matsimnibs[0][0], matsimnibs[1][0],
         matsimnibs[2][0]) = x_axis.tolist()
        (matsimnibs[0][1], matsimnibs[1][1],
         matsimnibs[2][1]) = y_axis.tolist()
        (matsimnibs[0][2], matsimnibs[1][2],
         matsimnibs[2][2]) = z_axis.tolist()
        (matsimnibs[0][3], matsimnibs[1][3],
         matsimnibs[2][3]) = center.tolist()
        matsimnibs[3][3] = 1

        return matsimnibs

    def calculateDistance(self, p1):
        p1_np = np.array(p1, 'float')
        center, z_axis = self.projectPoint(p1_np, smooth=True)
        return np.linalg.norm(p1_np - center)

    # Documentation for a method.
    # @param point 3-dimensional np vector
    # Calculates the closest triangle to a point, based on the triangle's
    # center
    def findClosestTriangle2Point(self, point):
        return np.argmin(np.sum(np.square(self.tr_centers - point), 1))

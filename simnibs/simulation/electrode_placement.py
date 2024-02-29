'''
    Electrode placement in meshes
    This program is part of the SimNIBS package.
    Please check on www.simnibs.org how to cite our work in publications.

    Copyright (C) 2018 Guilherme Saturnino

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


import warnings
import copy
import numpy as np
import scipy.spatial

from simnibs.utils.mesh_element_properties import ElementTags
from ..mesh_tools.mesh_io import _hash_rows
from ..utils.transformations import project_points_on_surface


def _remove_unconnected_triangles(mesh, roi_triangles, center,
                                  roi_nodes, triangles):
    cc = mesh.elm.connected_components(roi_triangles)
    if len(cc) > 1:
        _, center_index = mesh.find_closest_element(center, return_index=True,
                                                    elements_of_interest=roi_triangles)
        for comps in cc:
            if center_index in comps:
                roi_triangles = comps
        
        roi_nodes, triangles = \
        np.unique(mesh.elm[roi_triangles, :3], return_inverse=True)
        triangles = triangles.reshape(-1,3)
        
    return roi_nodes, triangles, roi_triangles
        
def _get_nodes_in_surface(mesh, surface_tags):
    tr_of_interest = (mesh.elm.elm_type == 2) * (np.in1d(mesh.elm.tag1, surface_tags))
    nodes_in_surface = np.unique(mesh.elm.node_number_list[tr_of_interest, :3])
    return nodes_in_surface

def _get_transform(center, y_axis, mesh, mesh_surface=[ElementTags.SCALP, ElementTags.SCALP_TH_SURFACE], y_type='relative',
                   nodes_roi=None):
    ''' Finds the transformation to make  '''
    center = np.array(center, dtype=float)

    c = project_points_on_surface(mesh, center, surface_tags = mesh_surface).flatten()

    if nodes_roi is None:
        nodes_in_surface = _get_nodes_in_surface(mesh, mesh_surface)
        kd_tree = scipy.spatial.cKDTree(mesh.nodes[nodes_in_surface])
        _, center_idx = kd_tree.query(center)
        normal = mesh.nodes_normals().value[nodes_in_surface - 1][center_idx]
    else:
        normal = np.nanmean(
            mesh.nodes_normals().value[nodes_roi - 1],
            axis=0
        )

    z_axis = normal / np.linalg.norm(normal)
    normal = z_axis

    if y_axis is None:
        if np.abs(normal[1]) < .8:
            y_axis = np.array([0., 1., 0.])
        else:
            y_axis = np.array([1., 0., 0.])
        y_type = 'absolute'
    y_axis = np.array(y_axis, dtype=float)

    if y_type == 'relative' and np.allclose(y_axis, center):
        raise ValueError('Y reference is at the same place as the electrode center')

    # Fixes the coordinate system if y_ref and the center are too close
    if y_type == 'relative':
        alpha = 0
        nr_iter = 0
        y_ref = c
        nodes_in_surface = _get_nodes_in_surface(mesh, mesh_surface)
        kd_tree = scipy.spatial.cKDTree(mesh.nodes[nodes_in_surface])
        while np.linalg.norm(c - y_ref) < 1e-5:
            _, y_ref_idx = kd_tree.query(y_axis + alpha * (y_axis - center))
            y_ref = mesh.nodes.node_coord[nodes_in_surface - 1][y_ref_idx]
            alpha += 1.
            nr_iter += 1
            if nr_iter == 10:
                raise ValueError('Could not define coordinate system: '
                                 'Y reference perpendicular to surface?')
        y_axis = (y_ref - c) - normal * normal.dot(y_ref - c)
        y_axis /= np.linalg.norm(y_axis)

    elif y_type == 'absolute':
        y_axis -= normal * normal.dot(y_axis)
        if np.isclose(np.linalg.norm(y_axis), 0):
            raise ValueError('Could not define coordinate system: '
                             'Y axis perpendicular to surface?')
        y_axis /= np.linalg.norm(y_axis)

    x_axis = np.cross(y_axis, -z_axis)  #  For coherence with the GUI
    affine = np.zeros((4, 4), dtype=float)
    affine[0, :3] = x_axis
    affine[1, :3] = y_axis
    affine[2, :3] = z_axis
    affine[:3, 3] = -affine[:3, :3].dot(c)
    affine[3, 3] = 1
    return affine, c

def _get_roi(center, radius, mesh, mesh_surface=[5, 1005], min_cos=.1):
    ''' Defines the region of interest of a given radius. Around a given center. Only
    node with normals pointing in the given direction are
    considered. Returns a list of triangles with at least 1 element in the ROI and their adjacent tetrahedra'''
    center = np.array(center, dtype=float)
    nodes_in_surface = _get_nodes_in_surface(mesh, mesh_surface)
    if len(nodes_in_surface) == 0:
        raise ValueError('Could not find surface {0} in mesh'.format(mesh_surface))

    distances = np.linalg.norm(center - mesh.nodes[nodes_in_surface], axis=1)
    normals = mesh.nodes_normals()[nodes_in_surface]
    center_normal = normals[np.argmin(distances)]
    
    # determine average edge length in roi by sampling
    avg_l = np.average(np.linalg.norm(mesh.nodes[mesh.elm[mesh.elm.triangles, 0][nodes_in_surface]] - \
                                      mesh.nodes[mesh.elm[mesh.elm.triangles, 1][nodes_in_surface]], axis=1))
    
    nodes_in_roi_ext = nodes_in_surface[(distances <= radius + 2*avg_l) *
                            (center_normal.dot(normals.T) > min_cos)]
    nodes_in_roi = nodes_in_surface[(distances <= radius) *
                            (center_normal.dot(normals.T) > min_cos)]

    tr_in_roi = np.any(
        np.in1d(mesh.elm[mesh.elm.triangles, :3], nodes_in_roi).reshape(-1, 3), axis=1)

    th_in_roi = np.any(
        np.in1d(mesh.elm[mesh.elm.tetrahedra], nodes_in_roi_ext).reshape(-1, 4), axis=1)

    roi_tr_nodes, roi_triangles_reordering = \
        np.unique(mesh.elm[mesh.elm.triangles[tr_in_roi], :3], return_inverse=True)

    roi_th_nodes, roi_tetrahedra_reordering = \
        np.unique(mesh.elm[mesh.elm.tetrahedra[th_in_roi]], return_inverse=True)

    return mesh.elm.triangles[tr_in_roi], roi_tr_nodes, roi_triangles_reordering.reshape(-1,3), \
        mesh.elm.tetrahedra[th_in_roi], roi_th_nodes, roi_tetrahedra_reordering.reshape(-1,4)

def _point_inside_polygon(vertices, points, tol=1e-3):
    '''uses the line-tracing algorithm Known issues: if the point is right in the edge,
    the behaviour becomes unstable '''
    vertices = np.array(vertices, dtype=float)
    vertices = np.vstack([vertices, vertices[0]])
    # Cound intersections along x and y directions
    intersection_counts = np.zeros((len(points)), dtype=int)
    s_x = np.maximum(np.max(vertices[:, 0]), np.max(points[:, 0])) -\
        np.minimum(np.min(vertices[:, 0]), np.min(points[:, 0]))
    s_y = np.maximum(np.max(vertices[:, 1]), np.max(points[:, 1])) -\
        np.minimum(np.min(vertices[:, 1]), np.min(points[:, 1]))
    r = 1.5 + np.random.rand(2)
    p = [r[0] * s_x, r[1] * s_y]
    for v in range(len(vertices) - 1):
        intersection_points, inside = \
            _line_line_intersect(vertices[v], vertices[v+1],
                                 points.T, (points + p).T)

        intersection_counts += inside

    return intersection_counts % 2 == 1

def _line_line_intersect(end_point_a1, end_point_a2, end_point_b1, end_point_b2, tol=1e-5):
    ''' Here end_points_a is 1-dimensional, while end_points_b might be 2D
    https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line'''
    x1 = end_point_a1[0]
    x2 = end_point_a2[0]
    x3 = end_point_b1[0]
    x4 = end_point_b2[0]
    y1 = end_point_a1[1]
    y2 = end_point_a2[1]
    y3 = end_point_b1[1]
    y4 = end_point_b2[1]
    den = (x1-x2) * (y3-y4) - (y1-y2) * (x3-x4)
    try:
        p = np.zeros((2, len(x3)), dtype=float)
        inside = np.zeros(len(x3), dtype=bool)
        paralel = np.abs(den) < tol
        p[:, paralel] = np.nan
        inside[paralel] = False

        x3 = x3[~paralel]
        x4 = x4[~paralel]
        y3 = y3[~paralel]
        y4 = y4[~paralel]
        den = den[~paralel]

        p[0, ~paralel] = ((x1*y2 - y1*x2) * (x3 - x4) - (x1-x2) * (x3*y4 - y3*x4)) / den
        p[1, ~paralel] = ((x1*y2 - y1*x2) * (y3 - y4) - (y1-y2) * (x3*y4 - y3*x4)) / den

        inside[~paralel] = p[0, ~paralel] >= np.minimum(x1, x2) - tol
        inside[~paralel] *= p[0, ~paralel] >= np.min(np.vstack([x3, x4]), axis=0) - tol

        inside[~paralel] *= p[0, ~paralel] <= np.maximum(x1, x2) + tol
        inside[~paralel] *= p[0, ~paralel] <= np.max(np.vstack([x3, x4]), axis=0) + tol

        inside[~paralel] *= p[1, ~paralel] >= np.minimum(y1, y2) - tol
        inside[~paralel] *= p[1, ~paralel] >= np.min(np.vstack([y3, y4]), axis=0) - tol

        inside[~paralel] *= p[1, ~paralel] <= np.maximum(y1, y2) + tol
        inside[~paralel] *= p[1, ~paralel] <= np.max(np.vstack([y3, y4]), axis=0) + tol

        return p, inside


    except TypeError:
        p = np.zeros((2), dtype=float)
        inside = False
        if np.abs(den) < tol:
            return [None, None], False

        Px = ((x1*y2 - y1*x2) * (x3 - x4) - (x1-x2) * (x3*y4 - y3*x4)) / den
        Py = ((x1*y2 - y1*x2) * (y3 - y4) - (y1-y2) * (x3*y4 - y3*x4)) / den

        inside = (Px >= min(x1, x2, x3, x4) - tol) and (Px < max(x1, x2, x3, x4) + tol)
        inside *= (Py >= min(y1, y2, y3, y4) - tol) and (Py < max(y1, y2, y3, y4) + tol)
        return np.array([Px, Py], dtype=float), inside

def _apply_affine(affine, p):
    p_inc = np.hstack([p, np.ones((len(p), 1))])
    return affine.dot(p_inc.T).T[:, :3]

def _point_line_distance(end_point_1, end_point_2, points, tol=1e-5):
    ''' Distance of points to line. Notice: Here we condider a line segment!'''
    '''Line Equation: y = t*x + b '''
    b = end_point_1
    x = (end_point_2 - end_point_1)
    t_max = np.linalg.norm(x)# length of line
    x = x / t_max # normalised direction
    normal = np.random.rand(len(end_point_1)) - .5
    normal -= normal.dot(x) * x
    normal /= np.linalg.norm(normal)

    # Find where in the line the point projects
    t = (points - b[None, :]).dot(x)
    # Fix t to be between t and t_max. Has buffer to ensure that 
    # nodes are move even if they lie outside of range, but close to end points of line.
    out_range = (t < -2) + (t > t_max + 2)
    t[t < 0] = 0.
    t[t > t_max] = t_max
    projected = t[:, None] * x + b
    dist = np.linalg.norm(points - projected, axis=1)
    # Return distance and closest point in line, and the side if the line
    # The side is arbitrary
    orthogonal = points - projected
    side = np.sign(orthogonal.dot(normal))
    side[dist < tol] = 0.
    side[out_range] = 0
    return dist, projected, side


def _triangle_baricenters(triangles, nodes):
    return np.average(nodes[triangles], axis=1)


def _edge_list(triangles):
    '''Gets the list of edges and adjacencies. This is a 2D version of the get_faces()
    method in mesh_io '''
    edges = triangles[:, [[0, 1], [0, 2], [1, 2]]]
    edges = edges.reshape(-1, 2)
    hash_array = _hash_rows(edges)
    unique, idx, inv, count = np.unique(hash_array, return_index=True,
                                        return_inverse=True, return_counts=True)
    edges = edges[idx]
    edge_adjacency_list = -np.ones((len(unique), 2), dtype=int)
    edge_adjacency_list[:, 0] = idx // 3

    if np.any(count > 2):
        warnings.warn(
            'Found an edge with more than 2 adjacent triangles!'
        )

    # Remove the edges already seen from consideration
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
    edges_repeated = np.searchsorted(unique, hash_array_reduced)
    # Finally, fill out the second column in the adjacency list
    edge_adjacency_list[edges_repeated, 1] = hash_array_repeated_pos // 3

    return edges, inv.reshape(-1, 3), edge_adjacency_list


def _baricentric_coordinates(point, triangle):
    '''here the triangle is a 3x2 array '''
    A = triangle[1:] - triangle[0]
    b = point - triangle[0]
    x = np.linalg.solve(A.T, b)
    bari = np.zeros(3)
    bari[:2] = x
    bari[2] = 1 - np.sum(bari)
    return bari

def _orientation(point, triangle):
    '''here the triangle is a 3x2 array '''
    #edges = triangle[[[0, 1], [2, 0], [1, 2]], :]
    edges = triangle[[[1, 0], [0, 2], [2, 1]], :]
    orientation = np.zeros(3, dtype=int)
    for i, e in enumerate(edges):
        orientation[i] = np.sign(np.linalg.det(e - point[None, :]))
    return orientation

def _calc_kdtree(triangles, nodes):
    tr_baricenters = _triangle_baricenters(triangles, nodes)
    return scipy.spatial.cKDTree(tr_baricenters)


def _triangle_with_points(points, triangles, nodes, eps=1e-5,
                          edge_list=None, kdtree=None):
    if edge_list is None:
        edges, tr_edges, adjacency_list = _edge_list(triangles)
    else:
        edges, tr_edges, adjacency_list = edge_list
    # Find triangles with points using a walking algorithm
    if kdtree is None:
        kdtree = _calc_kdtree(triangles, nodes)
    # Starting old_position for walking algorithm: the closest baricenter
    _, closest_tr = kdtree.query(points)
    #TODO: Walk in only a subset of the triangles (20 closest?)
    tr_with_points = np.zeros(len(points), dtype=int)
    # Walk
    for i, p, t in zip(range(len(points)), points, closest_tr):
        stop = False
        previous_t = 1e10
        nr_cycles = 0
        outside = False
        edge_order = np.arange(3)
        while not stop:
            stop = True
            np.random.shuffle(edge_order)
            orient = _orientation(p, nodes[triangles[t]])
            for e in edge_order:
                if orient[e] < 0:
                    adjacent = adjacency_list[tr_edges[t, e], 0]
                    if adjacent == t:
                        adjacent = adjacency_list[tr_edges[t, e], 1]
                    # if the edge has no adjacent tetrahedra
                    if adjacent == -1:
                        outside = True
                    # If this is not the triangle we just came from
                    elif adjacent != previous_t:
                        # Move to the adjacent triangle
                        # The walks are getting out here
                        previous_t = t
                        t = adjacent
                        stop = False
                        break
            # this is only here for ensure that it will not loop forever
            nr_cycles += 1
            if nr_cycles >= 1000:
                outside = True
                break

        if not outside:
            tr_with_points[i] = t
        if outside:
            tr_with_points[i] = -1

    return tr_with_points

def _calc_gamma(affected_nodes):
    ''' gamma value of tetrahedra 
    (see Parthasarathy et al., Finite Elements in Analysis and Design, 1994)

    Note: negative ...
    
    Parameters
    ------------
    affected_nodes: 
        node coordinates for each node in each triangle
    '''

    M = affected_nodes[:, 1:] - affected_nodes[:, 0, None]
    vol = np.linalg.det(M) / 6.
    if np.any(vol == 0.0):
        return np.inf

    # unique combination of nodes -> six edges
    unique_pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    edge_len = np.zeros((affected_nodes.shape[0], 6, affected_nodes.shape[2]))
    for e, (i, j) in enumerate(unique_pairs):
        edge_len[:, e, :] = affected_nodes[:, i, :] - affected_nodes[:, j, :]
    edge_rms = np.sqrt(np.sum(edge_len**2, axis=(1, 2))/ 6.0)
    gamma = edge_rms**3/vol
    gamma /= 8.479670
    # if volume is negative, th has been inverted ->
    # return high gamma to penalize
    gamma[np.signbit(vol)] = 100
    
    return gamma

def _move_point(new_position, to_be_moved, tr_nodes, th_nodes, roi_tr_nodes, roi_th_nodes, triangles, tetrahedra, poly, inv_affine, holes, el_layer, loop_count,
                edge_list=None, kdtree=None):
    '''
    Moves one point to new old_position on the line of the electrode. The gamma metric [source] is calculated for all adjacent tetrahedra,
    including the electrode tetrahedra built on the scalp. If the gamma metric becomes less than the minimal value,
    the movement of the current node is rejected.

    The triangle roi (tr roi) and tetrahedra roi (th roi) are indexes differently. The node to_be_moved in the tr roi
    is mapped to the original mesh indices and then to the th roi to determine the teterahedra that need to be checked for quality.
    
    new_position: 
        position in electrode plane, where node should be moved
    to_be_moved:
        index of node to be moved in tr roi
    tr_nodes:
        coordinates of triangles in tr roi
    th_nodes:
        coordinates of tetrahdera in th roi
    roi_tr_nodes:
        original node indices of triangles in roi
    roi_th_nodes:   
        original node indices of tetrahedra in roi
    triangles:
        node inidces belonging to one triangle (roi)
    tetrahedra:
        node inidces belonging to one tetrahedra (roi)
    poly:
        description of lines marking the electrode outline
    affine:
        rigid affine, mapping from head model space to electrode plane
    holes: 
        description of holes, needed to built electrode tetrahedra
    el_layer: float
        thickness of first electrode layer calculated by average edge length of triangles in roi
    
    '''
    # Certify that the new_position is inside the patch = adj_tr_nodes
    tr = _triangle_with_points(np.atleast_2d(new_position), triangles, tr_nodes[:, :2],
                               edge_list=edge_list, kdtree=kdtree)
    tr_with_node = np.where(np.any(triangles == to_be_moved, axis=1))[0]
    adj_tr_nodes = triangles[tr_with_node]
    if not np.in1d(tr, tr_with_node):
        return None, None
    
    # node movement in 2D electrode plane
    old_position = tr_nodes[to_be_moved] 
    d = new_position - old_position[:2]
    
    # The triangle node is indexed in the tr roi. To extract all adjacent tetrahedra,
    # the original index of to_be_moved is extracted from roi_tr_nodes[to_be_moved].
    th_with_node = np.any(roi_th_nodes[tetrahedra] == roi_tr_nodes[to_be_moved], axis=1)
    adj_th_nodes = roi_th_nodes[tetrahedra[th_with_node]]
    roi_adj_th_nodes = tetrahedra[th_with_node]
    adj_th_node_coords = th_nodes[roi_adj_th_nodes]  

    # which node in which tetrahedron is moved
    th_moved, node_moved = np.where(adj_th_nodes == roi_tr_nodes[to_be_moved])
    
    moved_tr_nodes = np.copy(tr_nodes)
    moved_th_nodes = adj_th_node_coords.copy()

    gamma_threshold = 7.0
    for t in np.linspace(loop_count/6, 0, num=10):

        moved_tr_nodes[to_be_moved, :2] = old_position[:2] + t*d
        # node movement in tetrahedra is done in original 3D space to assess th quality ->
        # apply inverse affine on node movement to map to original space
        moved_th_nodes[th_moved, node_moved, :2] = _apply_affine(inv_affine, (old_position + np.append(t*d, 0)).reshape(-1,3))[:, :2]
        # build electrode tetrahedra according to current movement
        el_nodes, el_tetrahedra_, _, _, _, _, _ = _build_electrode(poly, el_layer, moved_tr_nodes[:, :2], adj_tr_nodes, holes=holes, test_gamma=True) 
        el_nodes = el_nodes[el_tetrahedra_]
        # reorder elec th nodes, so volume is positive
        MM = el_nodes[:, 1:] - el_nodes[:, 0, None]
        vol = np.linalg.det(MM) / 6.
        neg_vol = np.signbit(vol)
        if np.any(neg_vol):
            tmp = el_nodes[neg_vol, 1]
            el_nodes[neg_vol, 1] = el_nodes[neg_vol, 0]
            el_nodes[neg_vol, 0] = tmp
            del tmp

        max_gamma = np.max(_calc_gamma(np.vstack((moved_th_nodes, el_nodes))))

        if max_gamma <= gamma_threshold:
            break

    return moved_tr_nodes, max_gamma

def _make_line(line, tr_nodes, th_nodes, roi_tr_nodes, roi_th_nodes, triangles, tetrahedra, poly, inv_affine, holes, el_layer, loop_count, ends=True):
    moved_tr_nodes = np.copy(tr_nodes)
    moved_th_nodes = np.copy(th_nodes)
    moved = []
    edge_list = _edge_list(triangles)
    if ends:
        # Find the triangles containing the end points
        tr_w_points = _triangle_with_points(line, triangles, tr_nodes[:, :2], edge_list=edge_list)
        for p, t in zip(line, tr_w_points):
            # If not outside
            if t != -1:
                # Try to move the nodes, Choose the movement that maximizes the minimum angle
                gammas = []
                moved_tr_nodes_list = []
                for n in triangles[t]:
                    mn, mg = _move_point(p, n, moved_tr_nodes, moved_th_nodes, roi_tr_nodes, roi_th_nodes, triangles, tetrahedra, poly, inv_affine, holes, el_layer, loop_count,
                                        edge_list=edge_list)
                    if mn is not None:
                        gammas.append(mg)
                        moved_tr_nodes_list.append(mn)
                if moved_tr_nodes_list != []:
                    gammas = np.array(gammas)
                    moved_tr_nodes = moved_tr_nodes_list[gammas.argmin()]
                    moved_th_nodes[np.where(np.isin(roi_th_nodes, roi_tr_nodes))[0]] = _apply_affine(inv_affine, moved_tr_nodes_list[gammas.argmin()])
                    moved.append(triangles[t, gammas.argmin()])
    # Finds the edges that are cut by the line segment
    edges, tr_edges, adjacency_list = _edge_list(triangles)

    _, closest, side = _point_line_distance(line[0], line[1], moved_tr_nodes[:, :2])
    edge_s = side[edges]
    edges_crossing = np.where(np.prod(edge_s, axis=1) < -.1)[0]
    kdtree = _calc_kdtree(triangles, moved_tr_nodes[:, :2])
    # For each edge crossing the line segment, try to move each node. choose the movement
    # that minimizes the maximum gamma value
    for e in edges_crossing:
        gammas = []
        moved_tr_nodes_list = []
        for i, n in enumerate(edges[e]):
            mn, mg = _move_point(closest[n], n, moved_tr_nodes, moved_th_nodes, roi_tr_nodes, roi_th_nodes, triangles, tetrahedra, poly, inv_affine, holes, el_layer, loop_count,
                                edge_list=edge_list, kdtree=kdtree)
            if mn is not None:
                gammas.append(mg)
                moved_tr_nodes_list.append(mn)
        if moved_tr_nodes_list != []:
            gammas = np.array(gammas)
            moved_tr_nodes = moved_tr_nodes_list[gammas.argmin()]
            moved_th_nodes[np.where(np.isin(roi_th_nodes, roi_tr_nodes))[0]] = _apply_affine(inv_affine, moved_tr_nodes_list[gammas.argmin()])
            moved.append(edges[e, gammas.argmin()])
    
    return moved_tr_nodes, moved_th_nodes, moved

def _draw_polygon_2D(poly, tr_nodes, th_nodes, roi_tr_nodes, roi_th_nodes, triangles, tetrahedra, inv_affine, holes, el_layer, loop_count, ends=True):
    moved_tr_nodes = np.copy(tr_nodes)
    moved_th_nodes = np.copy(th_nodes)
    moved = []
    for v in range(len(poly)):
        if v == len(poly) - 1:
            next_v = 0
        else:
            next_v = v+1
        vertices = poly[[v, next_v]]
        moved_tr_nodes, moved_th_nodes, m = _make_line(vertices, moved_tr_nodes, moved_th_nodes, roi_tr_nodes, roi_th_nodes, triangles, tetrahedra, poly, inv_affine, holes, el_layer, loop_count, ends=ends)
        moved.append(m)

    moved = [n for m in moved for n in m]

    return moved_tr_nodes, moved_th_nodes, moved


def _inside_complex_polygon(poly, nodes, triangles, holes=[], tol=1e-2):
    ''' Determines the triangles inside a complex polygon '''
    avg_l = np.average(
        np.linalg.norm(nodes[triangles[:, 0]] - nodes[triangles[:, 1]], axis=1))
    bar = np.mean(nodes[triangles], axis=1)
    tr_inside = np.where(_point_inside_polygon(poly, bar, tol=tol*avg_l))[0]
    if len(holes) > 0:
        tr_hole = []
        for h in holes:
            tr_hole.append(
                np.where(_point_inside_polygon(h, bar, tol=tol*avg_l))[0])
        tr_hole = np.unique(np.hstack(tr_hole))
        tr_inside = np.setdiff1d(tr_inside, tr_hole)
    return tr_inside


def _build_electrode(poly, height, nodes, triangles, holes=[], plug=None,
                     middle_layer=None, test_gamma=False):
    '''Builds the 3d structure of the electrode '''
    h = np.atleast_1d(height)
    assert np.all(np.array(height) > 0)
    assert len(h) <= 3
    if middle_layer is not None and len(h) < 3:
        raise ValueError('Can only define a middle layer with sandwich'
                         ' electrodes')
    
    if not test_gamma: 
        # Average length of an egde, determined by sampling
        avg_l = np.average(
            np.linalg.norm(nodes[triangles[:, 0]] - nodes[triangles[:, 1]], axis=1))
        
        tr_inside = _inside_complex_polygon(poly, nodes, triangles, holes=holes, tol=1e-2)
        if len(tr_inside) == 0:
            return None
        # Find the nodes in the area
        nodes_inside_num, tr = np.unique(triangles[tr_inside], return_inverse=True)

    else:
        # Average length of an egdes in roi equals the height in order to build one layer
        avg_l = height

        nodes_inside_num, tr = np.unique(triangles, return_inverse=True)
    
    tr = tr.reshape(-1, 3)
    nodes_inside = nodes[nodes_inside_num]

    # Determine the number of layers
    n_th_layers = np.rint(h/avg_l).astype(int)
    n_th_layers[n_th_layers < 1] = 1
    total_layers = np.sum(n_th_layers)
    layer_height = h/n_th_layers

    # Create the nodes in the main electrode volume (and also copy the surface ones)
    nodes_3d = np.zeros(((total_layers + 1)*nodes_inside.shape[0], 3), dtype=float)
    nodes_3d[:, :2] = np.tile(nodes_inside, (total_layers + 1, 1))

    n_in_layer = nodes_inside.shape[0]
    curr_height = 0
    curr_layer = 1
    for l_h, n_l in zip(layer_height, n_th_layers):
        for _ in range(1, n_l + 1):
            nodes_3d[curr_layer*n_in_layer:(curr_layer+1)*n_in_layer, 2] =\
                    curr_height + l_h
            curr_height += l_h
            curr_layer += 1

    # Vectors with the corresponding surface node of each electrode node
    corresponding = np.tile(nodes_inside_num, (total_layers + 1, ))
    # Vector with information on whether or not the point is in the bottom surface
    in_surface = np.zeros_like(corresponding, dtype=bool)
    in_surface[:n_in_layer] = 1

    min_vertex = tr.argmin(axis=1)
    tr_reordered = np.copy(tr)

    tr_reordered[min_vertex == 1, 0] = tr[min_vertex == 1, 1]
    tr_reordered[min_vertex == 1, 1] = tr[min_vertex == 1, 2]
    tr_reordered[min_vertex == 1, 2] = tr[min_vertex == 1, 0]

    tr_reordered[min_vertex == 2, 0] = tr[min_vertex == 2, 2]
    tr_reordered[min_vertex == 2, 1] = tr[min_vertex == 2, 0]
    tr_reordered[min_vertex == 2, 2] = tr[min_vertex == 2, 1]

    tetrahedra = []
    tetrahedra_idx = []
    
    count = 0
    c1 = tr_reordered[:, 1] < tr_reordered[:, 2]
    c2 = ~c1    

    for n, n_l in enumerate(n_th_layers):
        for nn in range(n_l):
            # we should vary the type of prism we get to that the faces are
            # properly shared
            # http://ldc.usb.ve/~vtheok/cursos/ci6322/escogidos/how%20to%20divide%20in%20tetrahedra.pdf
            l = count
            tetrahedra.append(
                np.vstack([l * n_in_layer + tr_reordered[:, 0],
                           (l+1) * n_in_layer + tr_reordered[:, 1],
                           (l+1) * n_in_layer + tr_reordered[:, 2],
                           (l+1) * n_in_layer + tr_reordered[:, 0]]).T)
            c1 = tr_reordered[:, 1] < tr_reordered[:, 2]
            tetrahedra.append(
                np.vstack([l * n_in_layer + tr_reordered[c1, 0],
                           l * n_in_layer + tr_reordered[c1, 1],
                           l * n_in_layer + tr_reordered[c1, 2],
                           (l+1) * n_in_layer + tr_reordered[c1, 2]]).T)
            tetrahedra.append(
                np.vstack([l * n_in_layer + tr_reordered[c1, 0],
                           l * n_in_layer + tr_reordered[c1, 1],
                           (l+1) * n_in_layer + tr_reordered[c1, 2],
                           (l+1) * n_in_layer + tr_reordered[c1, 1]]).T)

            c2 = ~c1
            tetrahedra.append(
                np.vstack([l * n_in_layer + tr_reordered[c2, 0],
                           l * n_in_layer + tr_reordered[c2, 1],
                           l * n_in_layer + tr_reordered[c2, 2],
                           (l+1) * n_in_layer + tr_reordered[c2, 1]]).T)
            tetrahedra.append(
                np.vstack([l * n_in_layer + tr_reordered[c2, 0],
                           (l+1) * n_in_layer + tr_reordered[c2, 1],
                           l * n_in_layer + tr_reordered[c2, 2],
                           (l+1) * n_in_layer + tr_reordered[c2, 2]]).T)
            count += 1
            tetrahedra_idx.append(n * np.ones(3 * len(tr_reordered), dtype=int))

    tetrahedra = np.vstack(tetrahedra)
    tetrahedra_idx = np.hstack(tetrahedra_idx).astype(int)

    el_triangles = []
    triangles_idx = []
    '''
    for c in components:
        faces = tetrahedra[tetrahedra_idx==c][:, [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]]
        faces = faces.reshape(-1, 3)
        hash_array = np.array([hash(f.tobytes()) for f in np.sort(faces, axis=1)])
        unique, idx, inv, count = np.unique(hash_array, return_index=True,
                                            return_inverse=True, return_counts=True)
        faces = faces[idx]
        el_triangles.append(faces[count == 1])
        triangles_idx.append(c * np.ones(len(el_triangles[-1])))
    '''
    if len(h) == 1 or len(h) == 2:
        layer_count = 0
        for n, n_l in enumerate(n_th_layers):
            layer_count += n_l
            el_triangles.append(tr + layer_count * n_in_layer)
            triangles_idx.append(n * np.ones(len(el_triangles[-1])))

    else:
        # Create the middle layer
        if middle_layer is None:
            middle_layer = poly
        tr_mid = _inside_complex_polygon(middle_layer, nodes, triangles, holes=holes, tol=1e-2)
        mid_tr = np.intersect1d(tr_inside, tr_mid)
        mid_nodes = np.unique(triangles[mid_tr])
        mid_th = np.all(np.isin(corresponding[tetrahedra], mid_nodes), axis=1)
        tetrahedra_idx[mid_th * (tetrahedra_idx == 1)] = -1
        tetrahedra_idx[tetrahedra_idx != -1] = 0
        tetrahedra_idx[tetrahedra_idx == -1] = 1
        # Make the electrode surfaces
        layer_count = np.sum(n_th_layers[:2])
        candidates = tr + layer_count * n_in_layer
        mid_tr = np.all(np.isin(corresponding[candidates], mid_nodes), axis=1)
        el_triangles.append(candidates[mid_tr])
        triangles_idx.append(1 * np.ones(len(el_triangles[-1])))

        candidates = tr + total_layers * n_in_layer
        el_triangles.append(tr + total_layers * n_in_layer)
        triangles_idx.append(0 * np.ones(len(el_triangles[-1])))

    el_triangles = np.vstack(el_triangles)
    triangles_idx = np.hstack(triangles_idx).astype(int)
 

    # If it has a plug
    if plug is not None:
        tr_plug = _inside_complex_polygon(plug, nodes, triangles, holes=holes, tol=1e-2)
        plug_tr = np.intersect1d(tr_inside, tr_plug)
        plug_nodes = np.unique(triangles[plug_tr])
        plug_tr = np.all(np.isin(corresponding[el_triangles], plug_nodes), axis=1)
        if len(h) == 1:
            i = 0
        else:
            i = 1
        #triangles_idx[plug_tr * (triangles_idx == i)] = -1
        plug_tr = plug_tr * (triangles_idx == i)
        el_triangles = np.vstack((el_triangles, el_triangles[plug_tr]))
        triangles_idx = np.hstack((triangles_idx,
                                   -np.ones(np.sum(plug_tr), dtype=int)))


    # Remove the triangles in the electrode-skin interface
    #triangles = triangles[np.sum(in_surface[triangles], axis=1) != 3]
    '''
    n_components = len(h)
    triangles = np.tile(tr, (n_components, 1))
    triangles_idx = np.ones(len(triangles), dtype=int)
    for n in n_th_layers:
        triangles[len(tr):(n+1)*len(tr)] += n * n_in_layer
        triangles[len(tr):(n+1)*len(tr)] += (n - 1) * n_in_layer
    '''

    return nodes_3d, tetrahedra, el_triangles, corresponding, in_surface, tetrahedra_idx,triangles_idx


def _build_electrode_on_mesh(center, ydir, poly, h, mesh, el_vol_tag, el_surf_tag,
                             on_top_of=[ElementTags.SCALP, ElementTags.SCALP_TH_SURFACE], holes=[], plug=None, plug_tag=None,
                             middle_layer=None):
    ''' Given the spatial localization of the electrode and the polygon description, add
    it to the mesh '''
    # Get Roi
    h = np.atleast_1d(h).astype(float)
    el_vol_tag = np.atleast_1d(el_vol_tag).astype(int)
    el_surf_tag = np.atleast_1d(el_surf_tag).astype(int)
    assert len(h) <= 3
    if len(h) == 1:
        assert len(el_vol_tag) == 1
        assert len(el_surf_tag) == 1
    else:
        assert len(el_vol_tag) == 2
        assert len(el_surf_tag) == 2

    if plug is not None:
        assert plug_tag is not None

    R = np.linalg.norm(poly, axis=1).max() * 1.2

    roi_triangles, roi_tr_nodes, triangles, _, _, _ = _get_roi(
             center, R, mesh, mesh_surface=on_top_of)
        
    # Sometimes the ROI can include nodes inside the head. Figure out the 
    # connected components and take the group, which includes the center node
    roi_tr_nodes, triangles, roi_triangles = \
    _remove_unconnected_triangles(mesh, roi_triangles, center, roi_tr_nodes, triangles)
    
    tr_nodes = mesh.nodes[roi_tr_nodes]
    affine, _ = _get_transform(
        center, ydir, mesh,
        mesh_surface=on_top_of,
        nodes_roi=roi_tr_nodes
    )
    tr_nodes = _apply_affine(affine, tr_nodes)
    tr_nodes_z = tr_nodes[:, 2]
    tr_nodes = tr_nodes[:, :2]

    # Build the electrode
    out = _build_electrode(
        poly, h, tr_nodes, triangles, holes=holes, plug=plug, middle_layer=middle_layer
    )
    if out is None:
        return mesh
    else:
        moved_tr_nodes, tetrahedra, triangles, corresponding, in_surf, th_tag, tr_tag = out
    
    ### Add electrodes to mesh, transform the nodes back
    inv_affine = np.linalg.inv(affine)

    tr_nodes = np.vstack([tr_nodes[:, :2].T, tr_nodes_z]).T
    tr_nodes = _apply_affine(inv_affine, tr_nodes)

    moved_tr_nodes[:, 2] += tr_nodes_z[corresponding]
    moved_tr_nodes = moved_tr_nodes[~in_surf]
    moved_tr_nodes = _apply_affine(inv_affine, moved_tr_nodes)

    # Add nodes to mesh
    mesh.nodes.node_coord[roi_tr_nodes-1, :] = tr_nodes
    mesh.nodes.node_coord = np.vstack(
        [mesh.nodes.node_coord, moved_tr_nodes])

    # Create dictonary in order to make the elements
    node_dict = np.zeros(len(corresponding), dtype=int)
    node_dict[in_surf] = roi_tr_nodes[corresponding[in_surf]]
    node_dict[~in_surf] = 1 + mesh.nodes.nr - len(moved_tr_nodes) + np.arange(len(moved_tr_nodes))

    # Tetrahedra tags
    tetra_tags = np.ones_like(th_tag, dtype=int)
    if len(h) == 1:
        tetra_tags[:] = el_vol_tag[0]
    else:
        tetra_tags[th_tag == 0] = el_vol_tag[0]
        tetra_tags[th_tag == 1] = el_vol_tag[1]

    # Make the tetrahedra
    mesh.elm.tag1 = np.hstack(
        [mesh.elm.tag1, tetra_tags])
    mesh.elm.tag2 = np.hstack(
        [mesh.elm.tag2, tetra_tags])
    mesh.elm.elm_type = np.hstack(
        [mesh.elm.elm_type, int(4) * np.ones(len(tetrahedra), dtype=int)])

    mesh.elm.node_number_list = np.vstack(
        [mesh.elm.node_number_list, node_dict[tetrahedra]])

    # Make the triangles
    s = np.min(np.where(mesh.elm.elm_type==4)[0])
    tr = node_dict[triangles]
    tr = np.concatenate((tr, -1 * np.ones((len(tr), 1), dtype=int)), axis=1)

    surf_tags = np.ones_like(tr_tag, dtype=int)
    if len(h) == 1:
        surf_tags[:] = el_surf_tag[0]
    else:
        surf_tags[tr_tag == 0] = el_surf_tag[0]
        surf_tags[tr_tag == 1] = el_surf_tag[1]

    if plug is not None:
        if not np.any(tr_tag == -1):
            raise ValueError('Could not find any triangles in the plug region')
        surf_tags[tr_tag == -1] = plug_tag

    mesh.elm.tag1 = np.insert(
        mesh.elm.tag1, s,  surf_tags)

    mesh.elm.tag2 = np.insert(
        mesh.elm.tag2, s,  surf_tags)

    mesh.elm.elm_type = np.insert(
        mesh.elm.elm_type, s,  int(2) * np.ones(len(tr), dtype=int))

    mesh.elm.node_number_list = np.insert(
       mesh.elm.node_number_list, s,  tr, axis=0)

    return mesh


def _create_polygon_from_elec(elec, mesh, skin_tag=[ElementTags.SCALP, ElementTags.SCALP_TH_SURFACE]):
    if elec.definition == 'plane':
        if elec.shape in ['rect', 'rectangle', 'ellipse']:
            if elec.dimensions is None or len(elec.dimensions) == 0:
                raise ValueError('Undefined electrode dimension')
            if np.any(np.array(elec.dimensions) <= 1e-5):
                raise ValueError('Found zero or negative value for electrode dimension')
            try:
                x = float(elec.dimensions[0]) / 2.
                y = float(elec.dimensions[1]) / 2.
            except IndexError:
                x = float(elec.dimensions[0]) / 2.
                y = float(elec.dimensions[0]) / 2.
            except TypeError:
                x = float(elec.dimensions) / 2.
                y = float(elec.dimensions) / 2.

        if elec.shape in ['rect', 'rectangle']:
            poly = np.array([[x, y], [-x, y], [-x, -y], [x, -y]], dtype=float)
        elif elec.shape == 'ellipse':
            angles = np.linspace(0, 2*np.pi, endpoint=False, num=20)
            poly = np.vstack([x*np.cos(angles), y*np.sin(angles)]).T
        elif elec.shape == 'custom':
            assert elec.vertices is not None, 'Undefined electrode vertices'
            poly = np.array(elec.vertices, dtype=float)
            assert poly.shape[0] > 2, 'Electrode must have more than 2 vertices'
            assert poly.shape[1] == 2, 'Vertices in electrode space must have 2 dimensions'
        else:
            raise ValueError('Invalid electrode shape: {0}, valid shapes are "rect",'
                             ' "ellipse" and "custom"'.format(elec.shape))
        center = np.array(elec.centre, dtype=float).squeeze()
        if len(center) == 2:  # If the center is 2-dimensional, it is relative
            poly += center
            if elec.pos_ydir is not None and len(elec.pos_ydir) == 2:
                raise NotImplementedError('Relative y axis not implemented yet')
        elif len(center) == 3:
            center = project_points_on_surface(mesh, center, surface_tags = skin_tag).flatten()
        else:
            raise ValueError('Wrong dimension of electrode centre: it should be 1x3 (or 1x2 for plugs)')

        try:
            pos_y = elec.pos_ydir.tolist()
        except:
            pos_y = elec.pos_ydir

        if pos_y:
            y_axis = np.array(elec.pos_ydir, dtype=float)
            if len(y_axis) == 3:
                y_axis = project_points_on_surface(mesh, y_axis, surface_tags = skin_tag).flatten()
        else:
            y_axis = None

    elif elec.definition == 'conf':
        v = np.array(elec.vertices)
        assert v.shape[0] > 2, 'Electrode must have more than 2 vertices'
        assert v.shape[1] == 3, 'Vertices in conform space must have 3 dimensions'
        # In this case, we disregard the electrode.center argument
        center = np.average(v, axis=0)
        R = np.linalg.norm(center - v, axis=1).max() * 1.2
        _, roi_nodes, _, _, _, _ = _get_roi(
            center, R, mesh, mesh_surface=skin_tag
        )
        transform, c = _get_transform(center, None, mesh, mesh_surface=skin_tag,
                                      nodes_roi=roi_nodes)
        center = c
        y_axis = None
        poly = _apply_affine(transform, v)
        poly = poly[:, :2]

    else:
        raise ValueError('Invalid electrode definition: {0}'.format(elec.definition))

    return poly, center, y_axis

def put_electrode_on_mesh(elec, mesh, elec_tag, skin_tag=[ElementTags.SCALP, ElementTags.SCALP_TH_SURFACE], extra_add=(ElementTags.SALINE_START - ElementTags.ELECTRODE_RUBBER_START),
                          surf_add=ElementTags.TH_SURFACE_START, plug_add=ElementTags.ELECTRODE_PLUG_SURFACE_START):
    ''' Places an electrode in the given mesh

    Parameters
    ---------------
    elec: simnibs.simulation.sim_struct.ELEC
        Structure with electrode information
    mesh: simnibs.msh.mesh_io.mesh
        Mesh where the electrode is to be placed
    elec_tag: int
        Tag for the electrode volume.
    skin_tag: int or list of ints (optional)
        Tags where the electrode is to be placed. Default: [5, 1005]
    extra_add: int
        Number to be added to electrode_tag to designate gel/sponge. Default:400
    surf_add:
        Number to be added to designate electrode surfaces. Default:1000
    plug_add:
        Number to be added to designate plug surfaces. Default:2000

    Returns:
    ---------
    mesh_w_electrode: simnibs.msh.mesh_io
        Mesh structure with the electrode added
    el_surf: int
        Electrode surface tag
    '''

    # Generate the polygon
    elec = copy.deepcopy(elec)
    elec.substitute_positions_from_cap()
    elec_poly, elec_center, ydir = _create_polygon_from_elec(elec, mesh, skin_tag=skin_tag)
    thick = np.atleast_1d(elec.thickness)
    holes_poly = []
    for h in elec.holes:
        holes_poly.append(
            _create_polygon_from_elec(h, mesh, skin_tag=skin_tag)[0])

    if elec.plug in [[], None]:
        plug_poly = None
    else:
        try:
            plug = elec.plug[0]
            if len(elec.plug) > 1:
                raise NotImplementedError
        except (IndexError, TypeError):
            plug = elec.plug
        plug_poly, plug_center, plug_ydir = \
            _create_polygon_from_elec(plug, mesh, skin_tag=skin_tag)

    if elec.dimensions_sponge is not None and len(elec.dimensions_sponge) > 0:
        if elec.definition != 'plane':
            raise NotImplementedError('Sponges are only implemented for plane'
                                      ' electrodes')
        if elec.shape not in ['rect', 'ellipse']:
            raise NotImplementedError('Sponges are only implemented for rectangular'
                                      ' and elliptical electrodes')
        if len(elec.thickness) != 3:
            raise ValueError('Sponge electrodes must have 3 thickness values!')

        has_sponge = True
        sponge_el = copy.deepcopy(elec)
        sponge_el.dimensions = elec.dimensions_sponge
        sponge_poly = _create_polygon_from_elec(sponge_el, mesh, skin_tag=skin_tag)[0]

    else:
        has_sponge = False

    if has_sponge:
        R = np.linalg.norm(sponge_poly, axis=1).max() * 1.2
    else:
        R = np.linalg.norm(elec_poly, axis=1).max() * 1.2


    roi_triangles, roi_tr_nodes, triangles, roi_tetrahedra, roi_th_nodes, tetrahedra =  _get_roi(
        elec_center, R, mesh, mesh_surface=skin_tag)
    
    # Sometimes the ROI can include nodes inside the head. Figure out the 
    # connected components and take the group, which includes the center node
    roi_tr_nodes, triangles, roi_triangles = \
        _remove_unconnected_triangles(mesh, roi_triangles, elec_center, roi_tr_nodes, triangles)

    tr_nodes = mesh.nodes[roi_tr_nodes]
    th_nodes = mesh.nodes[roi_th_nodes]
    mesh.fix_th_node_ordering()

    affine, elec_center = _get_transform(
        elec_center, ydir, mesh,
        mesh_surface=skin_tag,
        nodes_roi=roi_tr_nodes)
    tr_nodes = _apply_affine(affine, tr_nodes)
    tr_nodes_z = tr_nodes[:, 2]
    
    inv_affine = np.linalg.inv(affine)

    # thickness of one electrode layer like defined in _build_electrode() when calculating gamma metric
    # for moving nodes. One layer of electrode is built to check quality of electrode tetrahedra.
    el_layer = np.average(
            np.linalg.norm(tr_nodes[triangles[:, 0]] - tr_nodes[triangles[:, 1]], axis=1))

    # Draw polygon and holes
    moved = []
    ends = elec.shape != 'ellipse'
    for loop_count in range(1,7):
        tr_nodes, th_nodes, m = _draw_polygon_2D(
            elec_poly, tr_nodes, th_nodes, roi_tr_nodes, roi_th_nodes, triangles, tetrahedra, inv_affine, holes_poly, el_layer, loop_count, ends=ends)
    moved.append(m)

    for h_p, h in zip(holes_poly, elec.holes):
        ends = h.shape != 'ellipse'
        for loop_count in range(1,7):
            tr_nodes, th_nodes, m = _draw_polygon_2D(
                h_p, tr_nodes, th_nodes, roi_tr_nodes, roi_th_nodes, triangles, tetrahedra, inv_affine, holes_poly, el_layer, loop_count, ends=ends)
        moved.append(m)

    if plug_poly is not None:
        ends = plug.shape != 'ellipse'
        if len(plug_center) == 3:
            plug_center = _apply_affine(affine, plug_center[None, :])[0, :2]
            plug_poly += plug_center
        for loop_count in range(1,7):
            tr_nodes, th_nodes, m = _draw_polygon_2D(
                plug_poly, tr_nodes, th_nodes, roi_tr_nodes, roi_th_nodes, triangles, tetrahedra, inv_affine, holes_poly, el_layer, loop_count, ends=ends)
        moved.append(m)

    if has_sponge:
        ends = sponge_el.shape != 'ellipse'
        for loop_count in range(1,7):
            tr_nodes, th_nodes, m = _draw_polygon_2D(
                sponge_poly, tr_nodes, th_nodes, roi_tr_nodes, roi_th_nodes, triangles, tetrahedra, inv_affine, holes_poly, el_layer, loop_count, ends=ends)
        moved.append(m)

    # Optimize old_positions
    moved = [n for m in moved for n in m]


    # Change the mesh
    tr_nodes = np.vstack([tr_nodes[: ,:2].T, tr_nodes_z]).T
    tr_nodes = _apply_affine(inv_affine, tr_nodes)
    mesh.nodes.node_coord[roi_tr_nodes-1, :] = tr_nodes

    # Build electrodes
    if plug_poly is None:
        plug_poly = elec_poly

    if len(thick) == 1:
        mesh = _build_electrode_on_mesh(
            elec_center, ydir, elec_poly, thick, mesh,
            elec_tag + extra_add, elec_tag + extra_add + surf_add,
            on_top_of=skin_tag, holes=holes_poly, plug=plug_poly, plug_tag=elec_tag + plug_add)

    elif len(thick) == 2:
        mesh = _build_electrode_on_mesh(
            elec_center, ydir, elec_poly, thick, mesh,
            [elec_tag + extra_add, elec_tag],
            [elec_tag + extra_add + surf_add, elec_tag + surf_add],
            on_top_of=skin_tag, holes=holes_poly, plug=plug_poly, plug_tag=elec_tag + plug_add)

    elif len(thick) == 3:
        if not has_sponge:
            sponge_poly = copy.deepcopy(elec_poly)

        mesh = _build_electrode_on_mesh(
            elec_center, ydir, sponge_poly, thick, mesh,
            [elec_tag + extra_add, elec_tag],
            [elec_tag + extra_add + surf_add, elec_tag + surf_add],
            on_top_of=skin_tag, holes=holes_poly, plug=plug_poly,
            plug_tag=elec_tag + plug_add,
            middle_layer=elec_poly)
    else:
        raise ValueError('Electrodes must have 1, 2, or 3 layers')


    return mesh, elec_tag + plug_add

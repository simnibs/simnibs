import os

import copy
#import matplotlib.pyplot as plt
import numpy as np
import scipy.spatial
import pytest


import simnibs.msh.mesh_io as mesh_io
import simnibs.simulation.electrode_placement as electrode_placement
import simnibs.simulation.sim_struct as sim_struct


@pytest.fixture(scope='module')
def sphere3_msh():
    fn = os.path.join(os.path.dirname(os.path.realpath(
        __file__)), '..', 'testing_files', 'sphere3.msh')
    return mesh_io.read_msh(fn)




class TestDrawElectrodes:
    def test_get_transformation(self, sphere3_msh):
        # Notice, the coordinate system is left-handed for coherence with the GUI
        affine, c = electrode_placement._get_transform(
             [0., 0., 100.], None, sphere3_msh)
        assert np.allclose(c, [0., 0., 95.])
        assert np.allclose(affine[:3, :3].dot([0., 0., 1]), [0., 0., 1.], atol=1e-2)
        assert np.allclose(affine[:3, :3].dot([0., 1., 0]), [0., 1., 0.], atol=1e-2)
        assert np.allclose(affine[:3, :3].dot([1., 0., 0]), [-1., 0., 0.], atol=1e-2)
        assert np.allclose(affine.dot([0., 0., 100., 1.]), [0, 0, 5, 1], atol=1e-1,
                           rtol=1e-1)

        affine, c = electrode_placement._get_transform(
             [0., 100., 0.], [0., 100., 2.], sphere3_msh)
        assert np.allclose(c, [0., 95., 0.])
        assert np.allclose(affine[:3, :3].dot([0., 0., 1]), [0., 1., 0.], atol=1e-2)
        assert np.allclose(affine[:3, :3].dot([0., 1., 0]), [0., 0., 1.], atol=1e-2)
        assert np.allclose(affine[:3, :3].dot([1., 0., 0]), [1., 0., 0.], atol=1e-2)
        assert np.allclose(affine.dot([0., 100., 0., 1.]), [0, 0, 5, 1], atol=1e-1,
                           rtol=1e-1)

        affine, c = electrode_placement._get_transform(
             [0., 100., 0.], [0., 100., 2.], sphere3_msh)
        assert np.allclose(c, [0., 95., 0.])
        assert np.allclose(affine[:3, :3].dot([0., 0., 1]), [0., 1., 0.], atol=1e-2)
        assert np.allclose(affine[:3, :3].dot([0., 1., 0]), [0., 0., 1.], atol=1e-2)
        assert np.allclose(affine[:3, :3].dot([1., 0., 0]), [1., 0., 0.], atol=1e-2)
        assert np.allclose(affine.dot([0., 100., 0., 1.]), [0, 0, 5, 1], atol=1e-1,
                           rtol=1e-1)



    def test_get_roi(self, sphere3_msh):
        roi_triangles, roi_nodes, roi_triangles_reordering = electrode_placement._get_roi(
             [0., 0., 95.], 15, sphere3_msh)
        bar = sphere3_msh.elements_baricenters()
        assert np.all(np.linalg.norm(bar[roi_triangles] - [0., 0., 95.], axis=1) < 21)
        assert np.all(np.linalg.norm(sphere3_msh.nodes[roi_nodes] - [0., 0., 95.],
                                     axis=1) < 35)
        assert np.allclose(roi_nodes[roi_triangles_reordering],
                           sphere3_msh.elm[roi_triangles, :3])

    def test_line_line_intersection(self):
        line1 = np.array([[1., 1.], [-1., -1.]])
        line2 = np.array([[-1., 1.], [1., -1.]])
        p, intersect = electrode_placement._line_line_intersect(line1[0], line1[1],
                                                                line2[0], line2[1])
        assert np.allclose(p, 0)
        assert intersect
        line2 = np.array([[-1., 2.], [1., 0.]])
        p, intersect = electrode_placement._line_line_intersect(line1[0], line1[1],
                                                                line2[0], line2[1])
        assert np.allclose(p, .5)
        assert intersect

        line21 = np.array([[-1., 1.], [-1., 2.]]).T
        line22 = np.array([[1., -1.], [1., 0.]]).T
        p, intersect = electrode_placement._line_line_intersect(line1[0], line1[1],
                                                                line21, line22)
        assert np.allclose(p[:, 1], .5)
        assert np.allclose(p[:, 0], 0)
        assert np.all(intersect)
 
        line21 = np.array([[-1., 1.], [-1., 2.], [-1, 10]]).T
        line22 = np.array([[1., -1.], [1., 0.], [1, 8]]).T
        p, intersect = electrode_placement._line_line_intersect(line1[0], line1[1],
                                                                line21, line22)
        assert np.allclose(p[:, 1], .5)
        assert np.allclose(p[:, 0], 0)
        assert np.allclose(p[:, 2], 4.5)
        assert np.all(intersect[:2])
        assert ~intersect[2]

        line21 = np.array([[-1., 1.], [-1., 2.], [-1, 10], [1., 1.]]).T
        line22 = np.array([[1., -1.], [1., 0.], [1, 8], [2., 2.]]).T
        p, intersect = electrode_placement._line_line_intersect(line1[0], line1[1],
                                                                line21, line22)
        assert np.allclose(p[:, 1], .5)
        assert np.allclose(p[:, 0], 0)
        assert np.allclose(p[:, 2], 4.5)
        assert np.all(np.isnan(p[:, 3]))
        assert np.all(intersect[:2])
        assert np.all(~intersect[2:])

    def test_point_inside_polygon(self):
        vertices = [[5., -5.], [5., 5.], [-5., 5.], [-5., -5.]]
        X, Y = np.meshgrid(np.arange(-10, 10, dtype=float), np.arange(-10, 10, dtype=float))
        points = np.vstack([X.reshape(-1), Y.reshape(-1)]).T
        points *= 1.01
        inside = electrode_placement._point_inside_polygon(vertices, points)
        assert np.all(np.abs(points[inside, 0]) <= 5)
        assert np.all(np.abs(points[inside, 1]) <= 5)
        assert np.all((np.abs(points[~inside, 0]) > 5) + (np.abs(points[~inside, 1]) > 5))

    def test_calc_triangle_angles(self):
        tr = np.array([[[0., 0.], [0., 1.], [1., 0.]]])
        a = electrode_placement._calc_triangle_angles(tr)
        assert np.allclose(a, [np.pi/2, np.pi/4, np.pi/4])

        tr = np.array([[[0., 0.], [0., 1.], [1., 0.]],
                       [[0., 0.], [0., 1.], [np.sqrt(3) / 2., 1./2.]]])

        a = electrode_placement._calc_triangle_angles(tr)
        assert np.allclose(a[0], [np.pi/2, np.pi/4, np.pi/4])
        assert np.allclose(a[1], np.pi/3)

    def test_point_line_distance(self):
        np.random.seed(0)
        line = np.array([[3., 3.], [-3., -3.]])
        p = np.array([[1., -1.]])
        dist, closest_p, side = electrode_placement._point_line_distance(line[0], line[1], p)
        assert np.isclose(dist, np.sqrt(2))
        assert np.allclose(closest_p, 0)
        assert np.allclose(side, -1)

        p = np.array([[1., -1.], [-4, 0]])
        dist, closest_p, side = electrode_placement._point_line_distance(line[0], line[1], p)
        assert np.allclose(dist, [np.sqrt(2), 2*np.sqrt(2)])
        assert np.allclose(closest_p[0], 0)
        assert np.allclose(closest_p[1], -2)
        assert np.allclose(side, [1, -1])

        p = np.array([[1., -1.], [4., 0], [.5, .5]])
        dist, closest_p, side = electrode_placement._point_line_distance(line[0], line[1], p)
        assert np.allclose(dist, [np.sqrt(2), 2*np.sqrt(2), 0])
        assert np.allclose(closest_p[0], 0)
        assert np.allclose(closest_p[1], 2)
        assert np.allclose(closest_p[2], .5)
        assert np.allclose(side, [-1, -1, 0])

        p = np.array([[1., -1.], [4., 0], [.5, .5], [-7, -5]])
        dist, closest_p, side = electrode_placement._point_line_distance(line[0], line[1], p)
        assert np.allclose(dist, [np.sqrt(2), 2*np.sqrt(2), 0, np.linalg.norm([4, 2])])
        assert np.allclose(closest_p[0], 0)
        assert np.allclose(closest_p[1], 2)
        assert np.allclose(closest_p[2], .5)
        assert np.allclose(closest_p[3], -3)
        assert np.allclose(side, [-1, -1, 0, 0])


    def test_edge_list(self):
        X, Y = np.meshgrid(np.arange(-10, 11, dtype=float), np.arange(-10, 11, dtype=float))
        nodes = np.vstack([X.reshape(-1), Y.reshape(-1)]).T
        tri = scipy.spatial.Delaunay(nodes)
        edges, tr_edges, adjacency_list = electrode_placement._edge_list(tri.simplices)
        assert np.sum(np.in1d(tri.simplices[adjacency_list[0, 0]],
                              tri.simplices[adjacency_list[0, 1]])) == 2
        assert np.any(np.in1d(tr_edges[0], tr_edges[adjacency_list[tr_edges[0, 0], 1]]))

    def test_baricentric(self):
        tri = np.array([[0., 0.], [0., 1.], [np.sqrt(3) / 2., 1./2.]])
        p = np.average(tri, axis=0)
        bari = electrode_placement._baricentric_coordinates(p, tri)
        assert np.allclose(bari, 1./3.)

    def test_orientation(self):
        tri = np.array([[0., 0.], [np.sqrt(3) / 2., 1./2.], [0., 1.]])
        p = np.average(tri, axis=0)
        orient = electrode_placement._orientation(p, tri)
        assert np.allclose(orient, 1.)

        p = np.array([-1, 0.])
        orient = electrode_placement._orientation(p, tri)
        assert np.allclose(orient, [1, -1, 1])

    def test_fix_node_ordering(self):
        X, Y = np.meshgrid(np.arange(-10, 11, dtype=float), np.arange(-10, 11, dtype=float))
        nodes = np.vstack([X.reshape(-1), Y.reshape(-1)]).T
        tri = scipy.spatial.Delaunay(nodes)
        triangles = tri.simplices
        nodes = tri.points
        M = nodes[triangles][:, 1:, :] - nodes[triangles][:, 0, None]
        assert np.all(np.linalg.det(M) > 0)

        tmp = np.copy(triangles[:, 1])
        triangles[:, 1] = triangles[:, 0]
        triangles[:, 0] = tmp
        M = nodes[triangles][:, 1:, :] - nodes[triangles][:, 0, None]
        assert np.all(np.linalg.det(M) < 0)

        triangles = electrode_placement._fix_triangle_node_ordering(triangles, nodes)
        M = nodes[triangles][:, 1:, :] - nodes[triangles][:, 0, None]
        assert np.all(np.linalg.det(M) > 0)

 

    def test_triangle_with_points(self):
        X, Y = np.meshgrid(np.arange(-10, 11, dtype=float), np.arange(-10, 11, dtype=float))
        nodes = np.vstack([X.reshape(-1), Y.reshape(-1)]).T
        tri = scipy.spatial.Delaunay(nodes)
        query = np.random.rand(10, 2) * 20 - 10
        triangles = electrode_placement._triangle_with_points(query, tri.simplices,
                                                              tri.points)
        assert np.all(triangles == tri.find_simplex(query))

    def test_move_point(self):
        X, Y = np.meshgrid(np.arange(-10, 11, dtype=float), np.arange(-10, 11, dtype=float))
        nodes = np.vstack([X.reshape(-1), Y.reshape(-1)]).T
        tri = scipy.spatial.Delaunay(nodes)
        new_position = [.5, .5]
        candidate = np.argmin(np.linalg.norm(tri.points - new_position, axis=1))
        min_angle = 0.2
        new_points, angle = electrode_placement._move_point(new_position, candidate, tri.points,
                                                            tri.simplices, min_angle)
        angles = electrode_placement._calc_triangle_angles(new_points[tri.simplices])
        assert np.linalg.norm(new_points[candidate] - new_position) < 1e-5
        assert np.all(angles > 0.2)

        # Move too far out
        new_position = [.9, .9]
        new_points, angle = electrode_placement._move_point(new_position, candidate, tri.points,
                                                            tri.simplices, min_angle)
        angles = electrode_placement._calc_triangle_angles(new_points[tri.simplices])
        assert np.linalg.norm(new_points[candidate] - new_position) < 1
        assert np.all(angles > 0.2)
        # Move node outside
        new_position = [5, 5]
        new_points, angle = electrode_placement._move_point(new_position, candidate, tri.points,
                                                            tri.simplices, min_angle)
        assert new_points is None
        assert angle is None

    def test_make_line(self):
        #X, Y = np.meshgrid(np.arange(-10, 11, dtype=float), np.arange(-10, 11, dtype=float))
        X, Y = np.meshgrid(np.arange(-6, 7, dtype=float), np.arange(-2, 3, dtype=float))
        nodes = np.vstack([X.reshape(-1), Y.reshape(-1)]).T
        tri = scipy.spatial.Delaunay(nodes)
        line = np.array([[-5, .5], [5, .5]])
        new_points, _= electrode_placement._make_line(line, tri.points, tri.simplices,
                                                    ends=True)


        edges, tr_edges, adjacency_list = electrode_placement._edge_list(tri.simplices)
        _, closest, side = electrode_placement._point_line_distance(line[0], line[1], new_points)
        edge_s = side[edges]
        edges_crossing = np.prod(edge_s, axis=1) < -.1
        #bellow = (tri.points <= .5) * segment
        assert np.min(np.linalg.norm(new_points - line[0], axis=1)) < 1e-5
        assert np.min(np.linalg.norm(new_points - line[1], axis=1)) < 1e-5
        assert not np.all(edges_crossing)

        line = np.array([[-5, 0], [5, 0]])
        new_points, _= electrode_placement._make_line(line, tri.points, tri.simplices,
                                                    ends=True)
        edges, tr_edges, adjacency_list = electrode_placement._edge_list(tri.simplices)
        _, closest, side = electrode_placement._point_line_distance(line[0], line[1], new_points)
        edge_s = side[edges]
        edges_crossing = np.prod(edge_s, axis=1) < -.1
        #bellow = (tri.points <= .5) * segment
        assert np.min(np.linalg.norm(new_points - line[0], axis=1)) < 1e-5
        assert np.min(np.linalg.norm(new_points - line[1], axis=1)) < 1e-5
        assert not np.all(edges_crossing)
        assert np.allclose(new_points, tri.points)

    def test_draw_polygon_2D(self):
        X, Y = np.meshgrid(np.linspace(-8, 8, 25, dtype=float), np.linspace(-8, 8, 25, dtype=float))
        nodes = np.vstack([X.reshape(-1), Y.reshape(-1)]).T
        tri = scipy.spatial.Delaunay(nodes)
        poly = np.array([[5.5, 5.5], [5.5, -5.5], [-5.5, -5.5], [-5.5, 5.5]])
        #angles = np.linspace(0, 2*np.pi, endpoint=False, num=20)
        #poly = np.vstack([5*np.sin(angles), 5*np.cos(angles)]).T
        new_points, _ = electrode_placement._draw_polygon_2D(
            poly, tri.points, tri.simplices, ends=True)
        bar = np.mean(new_points[tri.simplices], axis=1)
        m = new_points[tri.simplices[:, 1:]] -\
            new_points[tri.simplices[:, 0]][:, None, :]
        area = .5 * np.linalg.det(m)
        inside = electrode_placement._point_inside_polygon(poly, bar, tol=1e-3)
        #plt.triplot(new_points[:, 0], new_points[:, 1], tri.simplices.copy())
        #plt.show()
        assert np.isclose(np.sum(area[inside]), 11**2)

    def test_inside_complex_polygon(self):
        X, Y = np.meshgrid(np.linspace(-8, 8, 17, dtype=float), np.linspace(-8, 8, 17, dtype=float))
        nodes = np.vstack([X.reshape(-1), Y.reshape(-1)]).T
        tri = scipy.spatial.Delaunay(nodes)
        poly = np.array([[5, 5], [5, -5], [-5, -5], [-5, 5]])
        hole1 = np.array([[2, 2], [2, -2], [-2, -2], [-2, 2]])
        hole2 = np.array([[4, 4], [4, 3], [3, 3], [3, 4]])
        inside = electrode_placement._inside_complex_polygon(poly, tri.points,
                                                             tri.simplices,
                                                             holes=[hole1, hole2],
                                                             tol=1e-3)
        m = tri.points[tri.simplices[inside, 1:]] -\
            tri.points[tri.simplices[inside, 0]][:, None, :]
        area = .5 * np.linalg.det(m)
        assert np.isclose(np.sum(area), 100-16-1)


    def test_build_electrode(self):
        X, Y = np.meshgrid(np.linspace(-8, 8, 17, dtype=float), np.linspace(-8, 8, 17, dtype=float))
        nodes = np.vstack([X.reshape(-1), Y.reshape(-1)]).T
        tri = scipy.spatial.Delaunay(nodes)
        poly = np.array([[5, 5], [5, -5], [-5, -5], [-5, 5]])
        #hole = [np.array([[2, 2], [2, -2], [-2, -2], [-2, 2]])]
        h = 5
        new_points, tetrahedra, triangles, _, _, _, _ = electrode_placement._build_electrode(
            poly, h, tri.points, tri.simplices)
        mesh = mesh_io.Msh()
        mesh.elm = mesh_io.Elements(#triangles=triangles+1)
                                 tetrahedra=tetrahedra + 1)
        mesh.nodes = mesh_io.Nodes(new_points)
        #mesh_io.write_msh(mesh, '~/Tests/electrode.msh')
        m = new_points[tetrahedra[:, 1:]] -\
            new_points[tetrahedra[:, 0]][:, None, :]
        vol = np.linalg.det(m) / 6.0
        assert np.isclose(np.sum(vol), 100 * h)

        sides = new_points[triangles[:, 1:]] - \
            new_points[triangles[:, 0]][:, None, :]

        n = np.cross(sides[:, 0], sides[:, 1])
        area = np.linalg.norm(n, axis=1) * 0.5
        assert np.isclose(np.sum(area), 100)
        #assert np.isclose(np.sum(area), 2 * 100 + 4 * 50)

        h = [3, 1]
        new_points, tetrahedra, triangles, _, _, th_tags, tr_tags = electrode_placement._build_electrode(
            poly, h, tri.points, tri.simplices)
        mesh = mesh_io.Msh()
        mesh.elm = mesh_io.Elements(#triangles=triangles+1)
                                 tetrahedra=tetrahedra + 1)
        mesh.elm.tag1 = th_tags + 1
        mesh.elm.tag2 = th_tags + 1
        mesh.nodes = mesh_io.Nodes(new_points)
        #mesh_io.write_msh(mesh, '~/Tests/electrode.msh')
        m = new_points[tetrahedra[:, 1:]] -\
            new_points[tetrahedra[:, 0]][:, None, :]
        vol = np.linalg.det(m) / 6.0
        assert np.isclose(np.sum(vol[th_tags==0]), 100 * 3)
        assert np.isclose(np.sum(vol[th_tags==1]), 100 * 1)
        mesh.elm = mesh_io.Elements(triangles=triangles+1)
        #mesh_io.write_msh(mesh, '~/Tests/electrode.msh')
        sides = new_points[triangles[:, 1:]] - \
            new_points[triangles[:, 0]][:, None, :]

        n = np.cross(sides[:, 0], sides[:, 1])
        area = np.linalg.norm(n, axis=1) * 0.5
        assert np.isclose(np.sum(area[tr_tags==0]), 100)
        assert np.isclose(np.sum(area[tr_tags==1]), 100)
        #assert np.isclose(np.sum(area[tr_tags==1]), 2 * 100 + 4 * 10)
        #assert np.isclose(np.sum(area[tr_tags==0]), 2 * 100 + 4 * 30)

        plug = np.array([[2, 2], [2, -2], [-2, -2], [-2, 2]])
        new_points, tetrahedra, triangles, _, _, th_tags, tr_tags = electrode_placement._build_electrode(
            poly, h, tri.points, tri.simplices, plug=plug)
        sides = new_points[triangles[:, 1:]] - \
            new_points[triangles[:, 0]][:, None, :]

        n = np.cross(sides[:, 0], sides[:, 1])
        area = np.linalg.norm(n, axis=1) * 0.5
        assert np.isclose(np.sum(area[tr_tags==0]), 100)
        assert np.isclose(np.sum(area[tr_tags==1]), 100)
        #assert np.isclose(np.sum(area[tr_tags==0]), 2 * 100 + 4 * 30)
        #assert np.isclose(np.sum(area[tr_tags==1]), 2 * 100 + 4 * 10 - 16)
        assert np.isclose(np.sum(area[tr_tags==-1]), 16)


        middle_layer = np.array([[2, 2], [2, -2], [-2, -2], [-2, 2]])
        h = [3, 1, 3]
        new_points, tetrahedra, triangles, _, _, th_tags, tr_tags = electrode_placement._build_electrode(
            poly, h, tri.points, tri.simplices, middle_layer=middle_layer)
        m = new_points[tetrahedra[:, 1:]] -\
            new_points[tetrahedra[:, 0]][:, None, :]
        vol = np.linalg.det(m) / 6.0
        assert np.isclose(np.sum(vol[th_tags==0]), 100 * 6 + (100-16) * 1)
        assert np.isclose(np.sum(vol[th_tags==1]), 16 * 1)
        sides = new_points[triangles[:, 1:]] - \
            new_points[triangles[:, 0]][:, None, :]

        n = np.cross(sides[:, 0], sides[:, 1])
        area = np.linalg.norm(n, axis=1) * 0.5
        assert np.isclose(np.sum(area[tr_tags==0]), 100)
        assert np.isclose(np.sum(area[tr_tags==1]), 16)
        #assert np.isclose(np.sum(area[tr_tags==0]), 2 * 100 + 4 * 30)
        #assert np.isclose(np.sum(area[tr_tags==1]), 2 * 100 + 4 * 10 - 16)




    def test_build_electrode_on_mesh(self):
        ''' Idea: Function to build a simple electrode (with holes) '''
        X, Y, Z = np.meshgrid(np.linspace(-8, 8, 17, dtype=float),
                              np.linspace(-8, 8, 17, dtype=float),
                              0)
        nodes = np.vstack([X.reshape(-1), Y.reshape(-1), Z.reshape(-1)]).T
        poly = np.array([[5, 5], [5, -5], [-5, -5], [-5, 5]])
        h = 5
        c = [0, 0, 0]
        ydir = [0, 1, 0]
        tri = scipy.spatial.Delaunay(nodes[:, :2])
        mesh = mesh_io.Msh()
        mesh.elm = mesh_io.Elements(triangles=tri.simplices+1)
        mesh.nodes = mesh_io.Nodes(nodes)
        w_elec = electrode_placement._build_electrode_on_mesh(c, ydir, poly, h, mesh, el_vol_tag=2,
                                                              el_surf_tag=1002, on_top_of=1)

        #mesh_io.write_msh(w_elec, '~/Tests/electrode.msh')
        vols = w_elec.elements_volumes_and_areas().value
        assert np.isclose(vols[w_elec.elm.tag1 == 2].sum(), 100 * h)
        assert np.isclose(vols[w_elec.elm.tag1 == 1002].sum(), 100)

        h = [5, 3]
        mesh = mesh_io.Msh()
        mesh.elm = mesh_io.Elements(triangles=tri.simplices+1)
        mesh.nodes = mesh_io.Nodes(nodes)
        w_elec = electrode_placement._build_electrode_on_mesh(c, ydir, poly, h,
                                                              copy.deepcopy(mesh),
                                                              el_vol_tag=[2, 3],
                                                              el_surf_tag=[1002, 1003], on_top_of=1)

        vols = w_elec.elements_volumes_and_areas().value
        assert np.isclose(vols[w_elec.elm.tag1 == 2].sum(), 100 * 5)
        assert np.isclose(vols[w_elec.elm.tag1 == 3].sum(), 100 * 3)
        assert np.isclose(vols[w_elec.elm.tag1 == 1002].sum(), 100)
        assert np.isclose(vols[w_elec.elm.tag1 == 1003].sum(), 100)

        plug = np.array([[2, 2], [2, -2], [-2, -2], [-2, 2]])
        w_elec = electrode_placement._build_electrode_on_mesh(c, ydir, poly, h,
                                                              copy.deepcopy(mesh),
                                                              el_vol_tag=[2, 3],
                                                              el_surf_tag=[1002, 1003],
                                                              on_top_of=1,
                                                              plug=plug,
                                                              plug_tag=2003)

        vols = w_elec.elements_volumes_and_areas().value
        assert np.isclose(vols[w_elec.elm.tag1 == 2].sum(), 100 * 5)
        assert np.isclose(vols[w_elec.elm.tag1 == 3].sum(), 100 * 3)
        assert np.isclose(vols[w_elec.elm.tag1 == 1002].sum(), 100)
        assert np.isclose(vols[w_elec.elm.tag1 == 1003].sum(), 100)
        assert np.isclose(vols[w_elec.elm.tag1 == 2003].sum(), 16)



        h = [5, 3, 1]
        middle_layer = np.array([[2, 2], [2, -2], [-2, -2], [-2, 2]])
        mesh = mesh_io.Msh()
        mesh.elm = mesh_io.Elements(triangles=tri.simplices+1)
        mesh.nodes = mesh_io.Nodes(nodes)
        w_elec = electrode_placement._build_electrode_on_mesh(c, ydir, poly, h,
                                                              copy.deepcopy(mesh),
                                                              el_vol_tag=[2, 3],
                                                              el_surf_tag=[1002, 1003],
                                                              on_top_of=1,
                                                              middle_layer=middle_layer)

        vols = w_elec.elements_volumes_and_areas().value
        assert np.isclose(vols[w_elec.elm.tag1 == 2].sum(), 100 * 6 + (100 - 16) * 3)
        assert np.isclose(vols[w_elec.elm.tag1 == 3].sum(), 16 * 3)
        assert np.isclose(vols[w_elec.elm.tag1 == 1002].sum(), 100)
        assert np.isclose(vols[w_elec.elm.tag1 == 1003].sum(), 16)


    '''
    def test_create_polygon_from_elec(self):
        X, Y, Z = np.meshgrid(np.linspace(-8, 8, 17, dtype=float),
                              np.linspace(-8, 8, 17, dtype=float),
                              0)
        nodes = np.vstack([X.reshape(-1), Y.reshape(-1), Z.reshape(-1)]).T
        tri = scipy.spatial.Delaunay(nodes[:, :2])
        mesh = mesh_io.Msh()
        mesh.elm = mesh_io.Elements(triangles=tri.simplices+1)
        mesh.nodes = mesh_io.Nodes(nodes)

        elec = sim_struct.ELECTRODE()
        elec.definition = 'plane'
        elec.shape = 'rect'
        elec.centre = [0, 0, 0]
        elec.dimX = 10
        elec.dimY = 8

        poly = electrode_placement._create_polygon_from_elec(elec, mesh)[0]
        assert np.allclose(poly, np.array([[5, 4], [-5, 4], [-5, -4], [5, -4]]))

        elec.centre = [2, 1]
        poly = electrode_placement._create_polygon_from_elec(elec, mesh)[0]
        assert np.allclose(poly, np.array([[7, 5], [-3, 5], [-3, -3], [7, -3]]))

        elec.definition = 'conf'
        elec.vertices = [[7, 5, 0], [-3, 5, 0], [-3, -3, 0], [7, -3, 0]]
        poly, centre, _ = electrode_placement._create_polygon_from_elec(elec, mesh, skin_tag=1)
        assert np.allclose(poly, np.array([[5, 4], [-5, 4], [-5, -4], [5, -4]]))
        assert np.allclose(centre, [2, 1, 0])

    def test_mesh_electrode(self):
        X, Y, Z = np.meshgrid(np.linspace(-8, 8, 17, dtype=float),
                              np.linspace(-8, 8, 17, dtype=float),
                              0)
        nodes = np.vstack([X.reshape(-1), Y.reshape(-1), Z.reshape(-1)]).T
        tri = scipy.spatial.Delaunay(nodes[:, :2])
        mesh = mesh_io.Msh()
        mesh.elm = mesh_io.Elements(triangles=tri.simplices+1)
        mesh.nodes = mesh_io.Nodes(nodes)

        elec = sim_struct.ELECTRODE()
        elec.definition = 'plane'
        elec.shape = 'rect'
        elec.centre = [1, 0, 0]
        elec.pos_ydir = [1, 1, 0]
        elec.dimX = 10
        elec.dimY = 5
        elec.thickness = 5
        w_elec = electrode_placement.put_electrode_on_mesh(elec, mesh, 2, skin_tag=1,
                                                           plug_add=2000)

        #mesh_io.write_msh(w_elec, '~/Tests/electrode1.msh')
        vols = w_elec.elements_volumes_and_areas().value
        assert np.isclose(vols[w_elec.elm.tag1 == 402].sum(), 50 * 5)
        assert np.isclose(vols[w_elec.elm.tag1 == 2002].sum(), 50)

        elec.shape = 'ellipse'
        elec.dimX = 8
        elec.dimY = 8
        elec.thickness = 5
        w_elec = electrode_placement.put_electrode_on_mesh(elec, mesh, 2, skin_tag=1)

        #mesh_io.write_msh(w_elec, '~/Tests/electrode2.msh')
        vols = w_elec.elements_volumes_and_areas().value
        assert np.isclose(vols[w_elec.elm.tag1 == 402].sum(), np.pi * 4**2 * 5, rtol=0.1)
        assert np.isclose(vols[w_elec.elm.tag1 == 2002].sum(), np.pi * 4**2, rtol=0.1)

        elec.shape = 'rect'
        elec.dimX = 10
        elec.dimY = 10
        elec.thickness = 5
        hole = sim_struct.ELECTRODE()
        hole.definition = 'plane'
        hole.centre = [1, 0]
        hole.shape = 'ellipse'
        hole.dimX = 6
        hole.dimY = 6
        elec.holes = [hole]
        w_elec = electrode_placement.put_electrode_on_mesh(elec, mesh, 2, skin_tag=1)

        #mesh_io.write_msh(w_elec, '~/Tests/electrode3.msh')
        vols = w_elec.elements_volumes_and_areas().value
        assert np.isclose(vols[w_elec.elm.tag1 == 402].sum(), (100 - np.pi * 3**2) * 5, rtol=0.1)

        elec.holes = []
        plug = sim_struct.ELECTRODE()
        plug.definition = 'plane'
        plug.centre = [.5, 0]
        plug.shape = 'rect'
        plug.dimX = 4
        plug.dimY = 4
        elec.plug = [plug]
        w_elec = electrode_placement.put_electrode_on_mesh(elec, mesh, 2, skin_tag=1)

        vols = w_elec.elements_volumes_and_areas().value
        #mesh_io.write_msh(w_elec, '~/Tests/electrode4.msh')
        assert np.isclose(vols[w_elec.elm.tag1 == 2002].sum(), 16, rtol=0.1)

        # Tests with 2 layers
        elec = sim_struct.ELECTRODE()
        elec.definition = 'plane'
        elec.shape = 'rect'
        elec.centre = [0, 0, 0]
        elec.pos_ydir = [0, 1, 0]
        elec.dimX = 10
        elec.dimY = 5
        elec.thickness = [5, 2]
        w_elec = electrode_placement.put_electrode_on_mesh(elec, mesh, 2, skin_tag=1)

        #mesh_io.write_msh(w_elec, '~/Tests/electrode5.msh')
        vols = w_elec.elements_volumes_and_areas().value
        assert np.isclose(vols[w_elec.elm.tag1 == 2].sum(), 50 * 2)
        assert np.isclose(vols[w_elec.elm.tag1 == 2002].sum(), 50)
        assert np.isclose(vols[w_elec.elm.tag1 == 402].sum(), 50 * 5)
        assert np.isclose(vols[w_elec.elm.tag1 == 1402].sum(), 50)

        #test with sponge
        elec.sponge_x = 12
        elec.sponge_y = 12
        elec.thickness = [5, 3, 5]
        elec.dimX = 10
        elec.dimY = 10
        w_elec = electrode_placement.put_electrode_on_mesh(elec, mesh, 2, skin_tag=1)

        #mesh_io.write_msh(w_elec, '~/Tests/electrode6.msh')
        vols = w_elec.elements_volumes_and_areas().value
        assert np.isclose(vols[w_elec.elm.tag1 == 2].sum(), 100 * 3)
        assert np.isclose(vols[w_elec.elm.tag1 == 2002].sum(), 100)
        assert np.isclose(vols[w_elec.elm.tag1 == 402].sum(), 144 * 10 + 44 * 3)
        assert np.isclose(vols[w_elec.elm.tag1 == 1402].sum(), 144)
    '''

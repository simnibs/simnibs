#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Surface_mesh.h>

// #include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
// #include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
// #include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point3;
typedef CGAL::Surface_mesh<Point3> Surface_mesh;
typedef Surface_mesh::Vertex_index vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;

typedef std::vector<std::vector<float>> matrix2d_float;
typedef std::vector<std::vector<int>> matrix2d_int;

namespace PMP = CGAL::Polygon_mesh_processing;


Surface_mesh construct_Surface_mesh(matrix2d_float vertices, matrix2d_int faces)
{
    CGAL::Real_timer timer;
    timer.start();

    Surface_mesh mesh;
    int n_vertices = vertices.size();
    int n_faces = faces.size();

    std::vector<vertex_descriptor> vertex_indices(n_vertices);
    for (int i = 0; i < n_vertices; ++i)
    {
        vertex_indices[i] = mesh.add_vertex(
            Point3(vertices[i][0], vertices[i][1], vertices[i][2]));
    }
    for (int i = 0; i < n_faces; ++i)
    {
        mesh.add_face(
            vertex_indices[faces[i][0]],
            vertex_indices[faces[i][1]],
            vertex_indices[faces[i][2]]);
    }
    assert(CGAL::is_triangle_mesh(mesh));

    // std::cout << "Elapsed time (mesh construction): " << timer.time() << std::endl;

    return mesh;
}

matrix2d_float vertices_as_vector(Surface_mesh mesh)
{

    int n_vertices = mesh.number_of_vertices();
    matrix2d_float vertices(n_vertices, std::vector<float>(3));

    for (vertex_descriptor vi : mesh.vertices())
    {
        Point3 p = mesh.point(vi);
        vertices[vi][0] = (float)p.x();
        vertices[vi][1] = (float)p.y();
        vertices[vi][2] = (float)p.z();
    }
    return vertices;
};

matrix2d_int faces_as_vector(Surface_mesh mesh)
{
    int n_faces = mesh.number_of_faces();
    matrix2d_int faces(n_faces, std::vector<int>(3));

    // for each face index, iterate over its halfedges and return all `target` vertices
    int i = 0;
    for (Surface_mesh::Face_index fi : mesh.faces())
    {
        int j = 0;
        Surface_mesh::Halfedge_index h = mesh.halfedge(fi);
        for (Surface_mesh::Halfedge_index hi : mesh.halfedges_around_face(h))
        {
            vertex_descriptor vi = mesh.target(hi);
            faces[i][j] = (int)vi;
            j++;
        }
        j = 0;
        i++;
    }
    return faces;
}

std::pair<matrix2d_float, matrix2d_int> get_vertices_and_faces(Surface_mesh mesh)
{
    auto pair = std::make_pair(vertices_as_vector(mesh), faces_as_vector(mesh));
    return pair;
}

matrix2d_int pmp_self_intersections(matrix2d_float vertices, matrix2d_int faces)
{

    Surface_mesh mesh = construct_Surface_mesh(vertices, faces);

    CGAL::Real_timer timer;
    timer.start();

    // std::cout << "Using parallel mode? " << std::is_same<CGAL::Parallel_if_available_tag, CGAL::Parallel_tag>::value << std::endl;

    std::vector<std::pair<face_descriptor, face_descriptor>> intersecting_tris;
    PMP::self_intersections<CGAL::Parallel_if_available_tag>(mesh.faces(), mesh, std::back_inserter(intersecting_tris));
    // std::cout << intersecting_tris.size() << " pairs of triangles intersect." << std::endl;

    int n_intersections = intersecting_tris.size();
    matrix2d_int intersecting_faces(n_intersections, std::vector<int>(2));
    for (int i = 0; i < n_intersections; i++)
    {
        intersecting_faces[i][0] = (int)intersecting_tris[i].first;
        intersecting_faces[i][1] = (int)intersecting_tris[i].second;
    }

    // std::cout << "Elapsed time (self intersections): " << timer.time() << std::endl;

    return intersecting_faces;
}


std::pair<std::vector<int>,std::vector<int>> pmp_connected_components(matrix2d_float vertices, matrix2d_int faces, std::vector<int> constrained_faces)
{
    Surface_mesh mesh = construct_Surface_mesh(vertices, faces);

    // should be enough to just construct a facelistgraph
    // Surface_mesh mesh = construct_FaceListGraph(faces);

    CGAL::Real_timer timer;
    timer.start();

    // Extract *all* edges of `constrained_faces` and use these as constraints
    // std::set<Surface_mesh::Edge_index> indices;
    // for (auto fi : constrained_faces)
    // {
    //     Surface_mesh::Halfedge_index h = mesh.halfedge((Surface_mesh::Face_index)fi);
    //     for (Surface_mesh::Halfedge_index hi : mesh.halfedges_around_face(h))
    //     {
    //         indices.insert(mesh.edge(hi));
    //     }
    // }
    // CGAL::Boolean_property_map<std::set<Surface_mesh::Edge_index>> constrained_edges_map(indices);
    // std::cout << "constraining " << indices.size() << " edges" << std::endl;

    // Extract the *outer* edges of `constrained_faces` and use these as constraints
    std::map<Surface_mesh::Edge_index, int> indices_with_counts;
    for (auto fi : constrained_faces)
    {
        Surface_mesh::Halfedge_index h = mesh.halfedge((Surface_mesh::Face_index)fi);
        for (Surface_mesh::Halfedge_index hi : mesh.halfedges_around_face(h))
        {
            auto edge = mesh.edge(hi);
            if (indices_with_counts.count(edge) == 0){
                // new edge
                indices_with_counts[edge] = 1;
            } else {
                // already seen edge
                indices_with_counts[edge]++;
            }
        }
    }
    // Keep only edges which occur once (i.e., "outer" edges)
    std::set<Surface_mesh::Edge_index> indices;
    for ( const auto &pair : indices_with_counts ) {
        if (pair.second == 1){
            indices.insert(pair.first);
        }
    }
    CGAL::Boolean_property_map<std::set<Surface_mesh::Edge_index>> constrained_edges_map(indices);
    // std::cout << "constraining " << indices.size() << " edges" << std::endl;

    // face component map (output)
    Surface_mesh::Property_map<face_descriptor, std::size_t> fccmap = mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;

    std::size_t num = PMP::connected_components(
        mesh,
        fccmap,
        CGAL::parameters::edge_is_constrained_map(constrained_edges_map));

    // std::cerr << "- The graph has " << num << " connected components (face connectivity)" << std::endl;

    typedef std::map<std::size_t /*index of CC*/, unsigned int /*nb*/> Components_size;
    Components_size nb_per_cc;
    for (face_descriptor f : mesh.faces())
    {
        nb_per_cc[fccmap[f]]++;
    }
    // for (const Components_size::value_type &cc : nb_per_cc)
    // {
    //     std::cout << "\t CC #" << cc.first
    //               << " is made of " << cc.second << " faces" << std::endl;
    // }
    // std::cout << "Elapsed time (connected components): " << timer.time() << std::endl;

    std::vector<int> cc(mesh.number_of_faces());
    std::vector<int> cc_size(num);
    for (face_descriptor f : mesh.faces())
    {
        cc[f] = (int)fccmap[f];
        cc_size[fccmap[f]]++;
    }
    auto pair = std::make_pair(cc, cc_size);

    return pair;
}

// std::pair<std::vector<int>,std::vector<int>> pmp_volume_connected_components(
//     matrix2d_int faces,
//     bool do_orientation_tests = false,
//     bool do_self_intersection_tests = false)
// {
//     Surface_mesh mesh = construct_FaceListGraph(faces);

//     CGAL::Real_timer timer;
//     timer.start();

//     // face component map (output)
//     Surface_mesh::Property_map<face_descriptor, std::size_t> fccmap = mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;

//     std::size_t num = PMP::volume_connected_components(
//         mesh,
//         fccmap,
//         CGAL::parameters::do_orientation_tests(do_orientation_tests).do_self_intersection_tests(do_self_intersection_tests));

//     // CGAL::parameters::do_orientation_tests(true).do_self_intersection_tests(true).is_cc_outward_oriented(true)

//     std::cerr << "- The graph has " << num << " connected components (face connectivity)" << std::endl;

//     typedef std::map<std::size_t /*index of CC*/, unsigned int /*nb*/> Components_size;
//     Components_size nb_per_cc;
//     for (face_descriptor f : mesh.faces())
//     {
//         nb_per_cc[fccmap[f]]++;
//     }
//     for (const Components_size::value_type &cc : nb_per_cc)
//     {
//         std::cout << "\t CC #" << cc.first
//                   << " is made of " << cc.second << " faces" << std::endl;
//     }
//     std::cout << "Elapsed time (connected components): " << timer.time() << std::endl;

//     std::vector<int> cc(mesh.number_of_faces());
//     std::vector<int> cc_size(num);
//     for (face_descriptor f : mesh.faces())
//     {
//         cc[f] = (int)fccmap[f];
//         cc_size[fccmap[f]]++;
//     }
//     auto pair = std::make_pair(cc, cc_size);

//     return pair;
// }

matrix2d_float pmp_smooth_shape(
    matrix2d_float vertices,
    matrix2d_int faces,
    std::vector<int> constrained_vertices,
    const double time,
    const unsigned int nb_iterations)
{
    Surface_mesh mesh = construct_Surface_mesh(vertices, faces);

    CGAL::Real_timer timer;
    timer.start();

    std::set<vertex_descriptor> indices;
    for (int i : constrained_vertices)
    {
        indices.insert((vertex_descriptor)i);
    }
    // std::cout << "constraining " << indices.size() << " vertices" << std::endl;
    CGAL::Boolean_property_map<std::set<vertex_descriptor>> vcmap(indices);

    PMP::smooth_shape(
        mesh,
        time,
        CGAL::parameters::number_of_iterations(nb_iterations)
            .vertex_is_constrained_map(vcmap));

    // std::cout << "Elapsed time (smoothing): " << timer.time() << std::endl;

    auto vertices_out = vertices_as_vector(mesh);

    return vertices_out;
}

// matrix2d_float pmp_angle_and_area_smoothing(
//     matrix2d_float vertices,
//     matrix2d_int faces,
//     std::vector<int> constrained_vertices,
//     const unsigned int nb_iterations,
//     bool use_angle_smoothing = true,
//     bool use_area_smoothing = true,
//     bool use_delaunay_flips = true,
//     bool use_safety_constraints = false)
// {
//     Surface_mesh mesh = construct_Surface_mesh(vertices, faces);

//     CGAL::Real_timer timer;
//     timer.start();

//     std::set<vertex_descriptor> indices;
//     for (int i : constrained_vertices)
//     {
//         indices.insert((vertex_descriptor)i);
//     }
//     std::cout << "constraining " << indices.size() << " vertices" << std::endl;
//     CGAL::Boolean_property_map<std::set<vertex_descriptor>> vcmap(indices);

//     PMP::angle_and_area_smoothing(mesh, CGAL::parameters::number_of_iterations(nb_iterations)
//                                             .use_angle_smoothing(use_angle_smoothing)
//                                             .use_area_smoothing(use_area_smoothing)
//                                             .use_Delaunay_flips(use_delaunay_flips)
//                                             .use_safety_constraints(use_safety_constraints)
//                                             .vertex_is_constrained_map(vcmap));

//     std::cout << "Elapsed time (smoothing): " << timer.time() << std::endl;

//     auto vertices_out = vertices_as_vector(mesh);

//     return vertices_out;
// }

// matrix2d_float pmp_fair(
//     matrix2d_float vertices,
//     matrix2d_int faces,
//     std::vector<int> indices)
// {
//     Surface_mesh mesh = construct_Surface_mesh(vertices, faces);

//     CGAL::Real_timer timer;
//     timer.start();

//     std::set<vertex_descriptor> vertex_indices;
//     for (int i : indices)
//     {
//         vertex_indices.insert((vertex_descriptor)i);
//     }
//     std::cout << "fairing " << vertex_indices.size() << " vertices" << std::endl;
//     CGAL::Boolean_property_map<std::set<vertex_descriptor>> vcmap(vertex_indices);

//     auto success = PMP::fair(mesh, vertex_indices);

//     std::cout << "Fairing : " << (success ? "succeeded" : "failed") << std::endl;
//     std::cout << "Elapsed time (fairing): " << timer.time() << std::endl;

//     auto vertices_faired = vertices_as_vector(mesh);

//     return vertices_faired;
// }


std::pair<matrix2d_float, matrix2d_int> pmp_isotropic_remeshing(
    matrix2d_float vertices,
    matrix2d_int faces,
    const double target_edge_length,
    const int n_iterations)
{
    Surface_mesh mesh = construct_Surface_mesh(vertices, faces);

    CGAL::Real_timer timer;
    timer.start();

    PMP::isotropic_remeshing(
        mesh.faces(),
        target_edge_length,
        mesh,
        CGAL::parameters::number_of_iterations(n_iterations)
    );

    // std::cout << "Elapsed time (isotropic remeshing): " << timer.time() << std::endl;

    // explicit garbage collection needed as vertices are only *marked* as removed
    //
    //   https://github.com/CGAL/cgal/discussions/6625
    //   https://doc.cgal.org/latest/Surface_mesh/index.html#sectionSurfaceMesh_memory
    mesh.collect_garbage();

    auto pair = get_vertices_and_faces(mesh);

    return pair;
}


// // Compute union between two meshes and refine.
// std::pair<matrix2d_float, matrix2d_int> pmp_corefine_and_union(
//     matrix2d_float vertices1,
//     matrix2d_int faces1,
//     matrix2d_float vertices2,
//     matrix2d_int faces2)
// {
//     Surface_mesh mesh1 = construct_Surface_mesh(vertices1, faces1);
//     Surface_mesh mesh2 = construct_Surface_mesh(vertices2, faces2);
//     Surface_mesh mesh_union;

//     bool valid_union = PMP::corefine_and_compute_union(mesh1, mesh2, mesh_union);
//     if (valid_union)
//     {
//         std::cout << "Union was successfully computed\n";
//         auto pair = get_vertices_and_faces(mesh_union);
//         return pair;
//     }
// }

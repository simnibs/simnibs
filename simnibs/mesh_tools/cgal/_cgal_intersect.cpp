#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>


#include <cstdlib>
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef K::Vector_3 Vector;
typedef FT (Function)(const Point&);
typedef CGAL::Surface_mesh<Point> Surface_mesh;

// Intersection
typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Segment>::Type> Segment_intersection;
typedef Tree::Primitive_id Primitive_id;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

std::pair<std::vector<int>, std::vector<float>> _segment_triangle_intersection(
        float* vertices, int n_vertices, int* tris, int n_faces,
        float* segment_start, float* segment_end, int n_segments)

{
  // Create mesh
  Surface_mesh m;
  std::vector<Surface_mesh::Vertex_index> vertex_indices(n_vertices);
  for(int i = 0; i < n_vertices; ++i) {
    vertex_indices[i] = m.add_vertex(Point(vertices[3*i], vertices[3*i+1], vertices[3*i+2]));
  }
  for(int i = 0; i < n_faces; ++i) {
    m.add_face(
        vertex_indices[tris[3*i]],
        vertex_indices[tris[3*i+1]],
        vertex_indices[tris[3*i+2]]);
  }
  // Create tree
  Tree tree(faces(m).first, faces(m).second, m);
  // Output vectors
  std::vector<int> indices_pairs;
  std::vector<float> intersect_positions;
  for (int i=0; i < n_segments; i++){
      Segment segment(
        Point(segment_start[3*i], segment_start[3*i + 1], segment_start[3*i + 2]),
        Point(segment_end[3*i], segment_end[3*i + 1], segment_end[3*i + 2])
      );
      // Test intersections
      std::list<Segment_intersection> intersections;
      tree.all_intersections(segment, std::back_inserter(intersections));
      for (std::list<Segment_intersection>::iterator it=intersections.begin(); it != intersections.end(); ++it){
          Segment_intersection inter = *it;
          std::size_t id = boost::get<CGAL::SM_Face_index>(inter->second);
          indices_pairs.push_back(i);
          indices_pairs.push_back((int) id);
          Point* p = boost::get<Point>(&(inter->first));
          Segment *s = boost::get<Segment>(&(inter->first));
          if (p) {
              for (int j=0; j<3; j++) intersect_positions.push_back((float) (*p)[j]);
           }
          else if (s) {
              for (int j=0; j<3; j++) intersect_positions.push_back((float) s->source()[j]);
          }
      }
  }
  return std::make_pair(indices_pairs, intersect_positions);
}

class TreeC{
  public:
  Tree tree;
  Surface_mesh m;
  TreeC() {tree = Tree();};
  ~TreeC();
  TreeC(float* , int , int* , int);
  std::pair<std::vector<int>, std::vector<float>> _intersections(float*, float*, int);
};

TreeC::TreeC(float* vertices, int n_vertices, int* tris, int n_faces) {
  // Create mesh
  this->m = Surface_mesh();
  std::vector<Surface_mesh::Vertex_index> vertex_indices(n_vertices);
  for(int i = 0; i < n_vertices; ++i) {
    vertex_indices[i] = m.add_vertex(Point(vertices[3*i], vertices[3*i+1], vertices[3*i+2]));
  }
  for(int i = 0; i < n_faces; ++i) {
    this->m.add_face(
        vertex_indices[tris[3*i]],
        vertex_indices[tris[3*i+1]],
        vertex_indices[tris[3*i+2]]);
  }
  // Create tree
  this->tree = Tree();
  this->tree.insert(faces(m).first, faces(m).second, m);
  this->tree.build();
};

TreeC::~TreeC () {
  this->tree.clear();
  this->m.clear();
}

std::pair<std::vector<int>, std::vector<float>> TreeC::_intersections(
        float* segment_start, float* segment_end, int n_segments) {
  // Output vectors
  std::vector<int> indices_pairs;
  std::vector<float> intersect_positions;
  for (int i=0; i < n_segments; i++){
      Segment segment(
        Point(segment_start[3*i], segment_start[3*i + 1], segment_start[3*i + 2]),
        Point(segment_end[3*i], segment_end[3*i + 1], segment_end[3*i + 2])
      );
      // Test intersections
      std::list<Segment_intersection> intersections;
      this->tree.all_intersections(segment, std::back_inserter(intersections));
      for (std::list<Segment_intersection>::iterator it=intersections.begin(); it != intersections.end(); ++it){
          Segment_intersection inter = *it;
          std::size_t id = boost::get<CGAL::SM_Face_index>(inter->second);
          indices_pairs.push_back(i);
          indices_pairs.push_back((int) id);
          Point* p = boost::get<Point>(&(inter->first));
          Segment *s = boost::get<Segment>(&(inter->first));
          if (p) {
              for (int j=0; j<3; j++) intersect_positions.push_back((float) (*p)[j]);
           }
          else if (s) {
              for (int j=0; j<3; j++) intersect_positions.push_back((float) s->source()[j]);
          }
      }
  }
  return std::make_pair(indices_pairs, intersect_positions);
}
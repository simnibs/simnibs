#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Polyhedral_complex_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>


#include <cstdlib>
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT (Function)(const Point&);
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain_img;
typedef CGAL::Polyhedral_complex_mesh_domain_3<K> Mesh_domain_surf;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain_img,CGAL::Default,Concurrency_tag>::type Tr_img;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr_img> C3t3_img;
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Mesh_triangulation_3<Mesh_domain_surf,CGAL::Default,Concurrency_tag>::type Tr_surf;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr_surf,Mesh_domain_surf::Corner_index,Mesh_domain_surf::Curve_index> C3t3_surf;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr_img> Mesh_criteria_img;
typedef Mesh_criteria_img::Facet_criteria    Facet_criteria_img;
typedef Mesh_criteria_img::Cell_criteria     Cell_criteria_img;

typedef CGAL::Mesh_criteria_3<Tr_surf> Mesh_criteria_surf;
typedef Mesh_criteria_surf::Facet_criteria    Facet_criteria_surf;
typedef Mesh_criteria_surf::Cell_criteria     Cell_criteria_surf;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;
// Function


int _mesh_image(
  char *fn_image, char *fn_out,
  float facet_angle, float facet_size, float facet_distance,
  float cell_radius_edge_ratio, float cell_size,
  bool optimize
)
{
  /// Load image
  CGAL::Image_3 image;
  if(!image.read(fn_image)){
    std::cerr << "Error: Cannot read file " <<  fn_image << std::endl;
    return EXIT_FAILURE;
  }
  // Mesh domain
  Mesh_domain_img domain = Mesh_domain_img::create_labeled_image_mesh_domain(image, 1e-10);

  // Mesh criteria
  Facet_criteria_img facet_criteria(facet_angle, facet_size, facet_distance);
  Cell_criteria_img cell_criteria(cell_radius_edge_ratio, cell_size);
  Mesh_criteria_img criteria(facet_criteria, cell_criteria);
  
  // Mesh generation
  C3t3_img c3t3 = CGAL::make_mesh_3<C3t3_img>(domain, criteria, no_perturb(), no_exude());
 
  if (optimize) CGAL::lloyd_optimize_mesh_3(c3t3, domain);
  CGAL::perturb_mesh_3(c3t3, domain);
  CGAL::exude_mesh_3(c3t3);
  

  // Output
  std::ofstream medit_file(fn_out);
  c3t3.output_to_medit(medit_file);
  return EXIT_SUCCESS;

}


int _mesh_surfaces(
  std::vector<char *>filenames, std::vector<std::pair<int, int> > incident_subdomains,
  char *fn_out,
  float facet_angle, float facet_size, float facet_distance,
  float cell_radius_edge_ratio, float cell_size,
  bool optimize
)

{
  std::vector<Polyhedron> patches(filenames.size());
  for(std::size_t i = 0; i < filenames.size(); ++i) {
    std::ifstream input(filenames[i]);

    std::vector<Point> points;
    std::vector< std::vector<std::size_t> > polygons;
    if (!CGAL::read_OFF(input, points, polygons))
    {
      std::cerr << "Error parsing the OFF file " << std::endl;
      return 1;
    }
    // We might need to apply some fixes to the surfaces
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, patches[i]);
  }
  // Create domain
  Mesh_domain_surf domain(
          patches.begin(), patches.end(),
          incident_subdomains.begin(), incident_subdomains.end()
  );

  //domain.detect_features(); //includes detection of borders

  // Mesh criteria
  Facet_criteria_surf facet_criteria(facet_angle, facet_size, facet_distance);
  Cell_criteria_surf cell_criteria(cell_radius_edge_ratio, cell_size);
  Mesh_criteria_surf criteria(facet_criteria, cell_criteria);

  // Mesh generation
  C3t3_surf c3t3 = CGAL::make_mesh_3<C3t3_surf>(domain, criteria, no_perturb(), no_exude());
  if (optimize) CGAL::lloyd_optimize_mesh_3(c3t3, domain);
  CGAL::perturb_mesh_3(c3t3, domain);
  CGAL::exude_mesh_3(c3t3);


  // Output
  std::ofstream medit_file(fn_out);
  c3t3.output_to_medit(medit_file);

  return EXIT_SUCCESS;
}



int _check_self_intersections(float *vertices, int n_vertices, int *faces, int n_faces)
{
  Surface_mesh m;
  std::vector<Surface_mesh::Vertex_index> vertex_indices(n_vertices);
  for(int i = 0; i < n_vertices; ++i) {
    vertex_indices[i] = m.add_vertex(Point(vertices[3*i], vertices[3*i+1], vertices[3*i+2]));
  }
  for(int i = 0; i < n_faces; ++i) {
    m.add_face(
        vertex_indices[faces[3*i]],
        vertex_indices[faces[3*i+1]],
        vertex_indices[faces[3*i+2]]);
  }
  bool intersecting = CGAL::Polygon_mesh_processing::does_self_intersect(m);
  return intersecting;
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_complex_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#ifdef CGAL_CONCURRENT_MESH_3
#include "tbb/task_arena.h"
#include "tbb/task_group.h"
#endif

#include <cstdlib>
// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef K::Vector_3 Vector;
typedef FT (Function)(const Point&);
typedef CGAL::Polyhedral_complex_mesh_domain_3<K> Mesh_domain_surf;
typedef CGAL::Surface_mesh<Point> Surface_mesh;
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Triangulation
typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
typedef CGAL::Mesh_triangulation_3<Mesh_domain_surf,CGAL::Default,Concurrency_tag>::type Tr_surf;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr_surf,Mesh_domain_surf::Corner_index,Mesh_domain_surf::Curve_index> C3t3_surf;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr_surf> Mesh_criteria_surf;
typedef Mesh_criteria_surf::Facet_criteria    Facet_criteria_surf;
typedef Mesh_criteria_surf::Cell_criteria     Cell_criteria_surf;
// To avoid verbose function and named parameters call
using namespace CGAL::parameters;


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
    if (!CGAL::IO::read_OFF(input, points, polygons))
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
  #ifdef CGAL_CONCURRENT_MESH_3
    tbb::task_arena limited(1);        // No more than 2 threads in this arena.
    tbb::task_group tg;
    limited.execute([&]{ // Use at most 2 threads for this job.
        tg.run([&]{ // run in task group
            if (optimize) CGAL::lloyd_optimize_mesh_3(c3t3, domain);
        });
    });
    // Wait for completion of the task group in the limited arena.
    limited.execute([&]{ tg.wait(); });
  #else
    if (optimize) CGAL::lloyd_optimize_mesh_3(c3t3, domain);
  #endif
  CGAL::perturb_mesh_3(c3t3, domain);
  CGAL::exude_mesh_3(c3t3);


  // Output
  std::ofstream medit_file(fn_out);
  c3t3.output_to_medit(medit_file);

  return EXIT_SUCCESS;
}


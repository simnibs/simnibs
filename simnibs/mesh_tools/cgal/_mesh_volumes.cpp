#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

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
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain_img;
#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif
// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain_img,CGAL::Default,Concurrency_tag>::type Tr_img;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr_img> C3t3_img;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr_img> Mesh_criteria_img;
typedef Mesh_criteria_img::Facet_criteria    Facet_criteria_img;
typedef Mesh_criteria_img::Cell_criteria     Cell_criteria_img;
// To avoid verbose function and named parameters call
using namespace CGAL::parameters;


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
  std::cout << "Began meshing \n";
  CGAL::get_default_random() = CGAL::Random(0);
  C3t3_img c3t3 = CGAL::make_mesh_3<C3t3_img>(domain, criteria, no_perturb(), no_exude());

  std::cout << "Lloyd \n";
  // Run Lloyd optimization using single core as it often fails in parallel
  // https://github.com/CGAL/cgal/issues/4566
  // When the bug gets fixed, please remove this whole block and only keep the simple version below
  // if (optimize) CGAL::lloyd_optimize_mesh_3(c3t3, domain);
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
  std::cout << "Perturb \n";
  CGAL::perturb_mesh_3(c3t3, domain);
  std::cout << "Exude \n";
  CGAL::exude_mesh_3(c3t3);
  std::cout << "Exude done\n";

  // Output
  std::ofstream medit_file(fn_out);
  c3t3.output_to_medit(medit_file);
  return EXIT_SUCCESS;

}

struct Sizing_field
{
    typedef ::FT FT;
    typedef Point Point_3;
    typedef Mesh_domain_img::Index Index;
    double tx, ty, tz;
    double vx, vy, vz;
    std::size_t sx, sy, sz;
    float *sizing_field_image;

    FT operator()(const Point_3 &p, const int, const Index&) const
    {
        // I add 0.5 offset because CGAL assumes centered voxels
        const std::size_t i = (p.x() - tx)/vx + 0.5;
        const std::size_t j = (p.y() - ty)/vy + 0.5;
        const std::size_t k = (p.z() - tz)/vz + 0.5;
        if (i < 0 || j < 0 || k < 0 || i >= sx || j >= sy || k >= sz) {
            std::cerr << "trying to access sizing field out-of-bounds" << std::endl;
            return 1;
        }
	      FT val = sizing_field_image[i + sx*j + sx*sy*k];
        return val;
    };
};

int _mesh_image_sizing_field(
  char *fn_image, char *fn_out,
  float facet_angle, float *facet_size, float *facet_distance,
  float cell_radius_edge_ratio, float *cell_size,
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
  Sizing_field sizing_field_cell = {
     image.tx(), image.ty(), image.tz(),
     image.vx(), image.vy(), image.vz(),
     image.xdim(), image.ydim(), image.zdim(),
     cell_size
  };
  Sizing_field sizing_field_facet_distance = {
     image.tx(), image.ty(), image.tz(),
     image.vx(), image.vy(), image.vz(),
     image.xdim(), image.ydim(), image.zdim(),
     facet_distance
  };
  Sizing_field sizing_field_facet_size = {
     image.tx(), image.ty(), image.tz(),
     image.vx(), image.vy(), image.vz(),
     image.xdim(), image.ydim(), image.zdim(),
     facet_size
  };
  Facet_criteria_img facet_criteria(
          facet_angle,
          sizing_field_facet_size,
          sizing_field_facet_distance
  );
  Cell_criteria_img cell_criteria(
          cell_radius_edge_ratio,
          sizing_field_cell
  );
  Mesh_criteria_img criteria(facet_criteria, cell_criteria);

  // Mesh generation
  std::cout << "Began meshing \n";
  C3t3_img c3t3 = CGAL::make_mesh_3<C3t3_img>(domain, criteria, no_perturb(), no_exude());

  std::cout << "Lloyd \n";
  // Run Lloyd optimization using single core as it often fails in parallel
  // https://github.com/CGAL/cgal/issues/4566
  // When the bug gets fixed, please remove this whole block and only keep the simple version
  // if (optimize) CGAL::lloyd_optimize_mesh_3(c3t3, domain);
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
  std::cout << "Perturb \n";
  CGAL::perturb_mesh_3(c3t3, domain);
  std::cout << "Exude \n";
  CGAL::exude_mesh_3(c3t3);


  // Output
  std::ofstream medit_file(fn_out);
  c3t3.output_to_medit(medit_file);
  return EXIT_SUCCESS;
}

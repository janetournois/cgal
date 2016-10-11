#define CGAL_MESH_3_VERBOSE 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/Mesh_3/sizing.h>
#include <CGAL/Mesh_3/lfs.h>

// IO
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/File_binary_mesh_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

typedef Sizing_grid<K> Sizing;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;



int main(int argc, char*argv[])
{
  const char* fname = (argc>1)?argv[1]:"data/elephant.off";
  const double k = (argc > 2) ? atof(argv[2]) : 1.;

  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input(fname);
  input >> polyhedron;
  if(input.fail()){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }
  input.close();
   
  // Create domain
  Mesh_domain domain(polyhedron);

  // Init sizing
  Sizing size(k);
  compute_lfs_sizing<K>(polyhedron, size);

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(facet_angle=25
                       , facet_size=0.15
                       , facet_distance=0.001
                       , cell_radius_edge_ratio=3
                       , cell_size = size);
  
  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

  // Output
  std::ofstream medit_file("out.mesh");
  c3t3.output_to_medit(medit_file);
  medit_file.close();

  std::ofstream outfile;
  outfile.open("out.binary.cgal",
    std::ios_base::binary | std::ios_base::out);
  CGAL::Mesh_3::save_binary_file(outfile, c3t3, true);

  return 0;
}

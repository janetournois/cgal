#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/Adaptive_remeshing_sizing_field.h>

#include <CGAL/IO/File_medit.h>

// Domain
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using FT = K::FT;
using Polyhedron = CGAL::Polyhedron_3<K>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_3<Polyhedron, K>;

// Triangulation for Meshing
using Tr = CGAL::Mesh_triangulation_3<Mesh_domain>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

// Triangulation for Remeshing
using T3 = CGAL::Triangulation_3<Tr::Geom_traits,
                                 Tr::Triangulation_data_structure>;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elk.off");
  std::ifstream input(fname);
  Polyhedron polyhedron;
  input >> polyhedron;
  if (input.fail()) {
    std::cerr << "Error: Cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  if (!CGAL::is_triangle_mesh(polyhedron)) {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  // Create domain
  Mesh_domain domain(polyhedron);

  std::cout << "Meshing...";
  std::cout.flush();

  // Mesh criteria
  Mesh_criteria criteria(facet_angle = 25,
                         facet_distance = 0.2,
                         cell_radius_edge_ratio = 3);

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb().no_exude());

  std::cout << "\rMeshing done." << std::endl;

  T3 tr = CGAL::convert_to_triangulation_3(std::move(c3t3));
  //note we use the move semantic, with std::move(c3t3),
  //  to avoid a copy of the triangulation by the function
  //  `CGAL::convert_to_triangulation_3()`
  //  After the call to this function, c3t3 is an empty and valid C3t3.
  //It is possible to use :  CGAL::convert_to_triangulation_3(c3t3),
  //  Then the triangulation is copied and duplicated, and c3t3 remains as is.

  std::cout << "Remeshing...";
  std::cout.flush();

  CGAL::Tetrahedral_remeshing::Adaptive_remeshing_sizing_field<T3>
    adaptive_field(tr);

  CGAL::tetrahedral_isotropic_remeshing(tr,
    adaptive_field,
    CGAL::parameters::number_of_iterations(5));

  std::cout << "\rRemeshing done." << std::endl;

  return EXIT_SUCCESS;
}

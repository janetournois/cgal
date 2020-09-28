
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

int main(int argc, char* argv[])
{
  using namespace CGAL;
  using Kernel = Exact_predicates_inexact_constructions_kernel;
  using Mesh = Surface_mesh<Kernel::Point_3>;

  Mesh mesh;
  Mesh::Vertex_index v1 = mesh.add_vertex(Mesh::Point(0, 0, 0));
  Mesh::Vertex_index v2 = mesh.add_vertex(Mesh::Point(0, 0, 0));
  Mesh::Vertex_index v3 = mesh.add_vertex(Mesh::Point(1, 0, 0));
  Mesh::Vertex_index v4 = mesh.add_vertex(Mesh::Point(0, 1, 0));
  mesh.add_face(v1, v2, v4);
  mesh.add_face(v2, v3, v4);
  mesh.add_face(v1, v3, v2);

  Mesh::Property_map<Mesh::Vertex_index, double> dmap
    = mesh.add_property_map<Mesh::Vertex_index, double>("v:dist").first;

  using Mode = Heat_method_3::Intrinsic_Delaunay;
  Heat_method_3::Surface_mesh_geodesic_distances_3<Mesh, Mode> heat(mesh);
  heat.add_source(v1);
  heat.estimate_geodesic_distances(dmap);

  return EXIT_SUCCESS; // this will never happen
}

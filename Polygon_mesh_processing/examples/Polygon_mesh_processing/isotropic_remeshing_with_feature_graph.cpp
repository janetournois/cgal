#define CGAL_PMP_DEBUG_FEATURE_GRAPH

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/feature_graph.h>

#include <boost/iterator/function_output_iterator.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/IO/OM.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;

using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;


int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/pig.off");

  std::unordered_map<vertex_descriptor, bool> sm_vfeature_map;
  auto sm_vfeature_pmap = boost::make_assoc_property_map(sm_vfeature_map);

  std::unordered_map<edge_descriptor, bool> sm_efeature_map;
  auto sm_efeature_pmap = boost::make_assoc_property_map(sm_efeature_map);

  // Try building a surface_mesh
  Mesh sm;
  bool ok = CGAL::IO::read_polygon_mesh(filename, sm,
                                        params::vertex_is_constrained_map(sm_vfeature_pmap)
                                       .edge_is_constrained_map(sm_efeature_pmap));


  const double snap_distance = 2.;

  PMP::snap_endpoints(sm_efeature_pmap, sm, snap_distance, params::min_feature_length(3));

  std::cout << "Snapping done." << std::endl;
  std::cout << "sm size = " << sm.number_of_vertices() << std::endl;

  CGAL::IO::write_OM("out_features_snapped.om", sm,
                      params::edge_is_constrained_map(sm_efeature_pmap)
                     .stream_precision(17));

//  CGAL::IO::write_polygon_mesh("out.off", mesh, CGAL::parameters::stream_precision(17));

//  std::cout << "Remeshing done." << std::endl;

  return 0;
}

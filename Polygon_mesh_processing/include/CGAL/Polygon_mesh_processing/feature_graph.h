// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois
//

#ifndef CGAL_POLYGON_MESH_PROCESSING_FEATURE_GRAPH_H
#define CGAL_POLYGON_MESH_PROCESSING_FEATURE_GRAPH_H

#include <CGAL/license/Polygon_mesh_processing/feature_graph.h>

#include <CGAL/Mesh_3/polylines_to_protect.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

#include <CGAL/boost/graph/shortest_path.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/graph/adjacency_list.hpp>

#include <unordered_set>
#include <unordered_map>

namespace CGAL {
namespace Polygon_mesh_processing {

  namespace internal
  {
    template<typename V, typename G>
    struct Features_visitor
    {
      // V is for vertex_descriptor of PolygonMesh
      std::vector<std::vector<V> >& polylines;
      G& graph;
      Features_visitor(typename std::vector<std::vector<V> >& lines, G& p_graph)
        : polylines(lines), graph(p_graph)
      {
      }
      void start_new_polyline()
      {
        std::vector<V> polyline;
        polylines.push_back(polyline);
      }
      void add_node(typename boost::graph_traits<G>::vertex_descriptor vd)
      {
        std::vector<V>& polyline = polylines.back();
        polyline.push_back(graph[vd]);
      }
      void end_polyline()
      {
        // ignore degenerated polylines
        if (polylines.back().size() < 2)
          polylines.resize(polylines.size() - 1);
      }
    };

    template <typename Graph>
    struct Less_for_Graph_vertex_descriptors
    {
      const Graph& graph;
      Less_for_Graph_vertex_descriptors(const Graph& graph) : graph(graph) {}

      template <typename vertex_descriptor>
      bool operator()(vertex_descriptor v1, vertex_descriptor v2) const {
        return graph[v1] < graph[v2];
      }
    };
  }

  /*!
  * GT
  * VPMap
  * min_feature_length (shorter polylines are removed)
  */
  template <typename EdgeIsConstrainedMap,
            typename PolygonMesh,
            typename NamedParameters = parameters::Default_named_parameters>
  void snap_endpoints(EdgeIsConstrainedMap& ecmap,
                      const PolygonMesh& mesh,
                      const double snap_distance,
                      const NamedParameters& np = parameters::default_values())
  {
    using PM = PolygonMesh;
    using vertex_descriptor = typename boost::graph_traits<PM>::vertex_descriptor;
    using edge_descriptor = typename boost::graph_traits<PM>::edge_descriptor;
    using face_descriptor = typename boost::graph_traits<PM>::face_descriptor;

    using GT = typename GetGeomTraits<PM, NamedParameters>::type;
    using FT = typename GT::FT;
    using Point_3 = typename GT::Point_3;
    using VPMap = typename GetVertexPointMap<PM, NamedParameters>::type;

    using Tree_traits_base = CGAL::Search_traits_3<GT>;
    using Tree_traits = CGAL::Search_traits_adapter<vertex_descriptor, VPMap, Tree_traits_base>;
    using K_neighbor_search = CGAL::Orthogonal_k_neighbor_search<Tree_traits>;

    using Tree = typename K_neighbor_search::Tree;
    using Splitter = typename Tree::Splitter;
    using Distance = typename K_neighbor_search::Distance;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
    VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                   get_const_property_map(vertex_point, mesh));
    std::size_t min_feature_length = choose_parameter(get_parameter(np, internal_np::min_feature_length),
                                                      1);

    using Polyline = std::vector<vertex_descriptor>;
    using Polylines = std::vector<Polyline>;

    // collect all the constrained edges
    Polylines constrained_edges;
    for (edge_descriptor e : edges(mesh))
    {
      if (get(ecmap, e))
      {
        vertex_descriptor v0 = source(halfedge(e, mesh), mesh);
        vertex_descriptor v1 = target(halfedge(e, mesh), mesh);
        constrained_edges.push_back({v0, v1});
      }
    }
#ifdef CGAL_PMP_DEBUG_FEATURE_GRAPH
    std::cout << "# constrained edges in input: " << constrained_edges.size() << std::endl;
#endif

    // organize them into a graph
    using Graph = boost::adjacency_list<boost::setS, boost::vecS,
                                        boost::undirectedS,
                                        vertex_descriptor>;
    using Graph_v = typename boost::graph_traits<Graph>::vertex_descriptor;
    using Polyline_v = typename std::iterator_traits<typename Polylines::iterator>::value_type;

    Graph feature_graph;
    std::unordered_map<vertex_descriptor, Graph_v> mv2gv; //mesh vertex to graph vertex
    auto find_or_add_vertex = [&feature_graph, &mv2gv](vertex_descriptor v)
      {
        auto it = mv2gv.find(v);
        if (it != mv2gv.end())
          return it->second;
        Graph_v gv = add_vertex(feature_graph);
        mv2gv[v] = gv;
        feature_graph[gv] = v;
        return gv;
      };

    for (const Polyline& polyline : constrained_edges)
    {
      if (polyline.size() < 2)
        continue;

      auto vit = polyline.begin();
      vertex_descriptor v = *vit;
      Graph_v gv = find_or_add_vertex(v);

      while (std::next(vit) != polyline.end())
      {
        Graph_v gw = find_or_add_vertex(*std::next(vit));
        add_edge(gv, gw, feature_graph);
        gv = gw;
        ++vit;
      }
    }

    Polylines polylines;
    internal::Features_visitor<vertex_descriptor, Graph> visitor(polylines, feature_graph);
    PMP::internal::Less_for_Graph_vertex_descriptors<Graph> less(feature_graph);
    const Graph& const_graph = feature_graph;
    split_graph_into_polylines(const_graph, visitor, CGAL::internal::IsTerminalDefault(), less);

#ifdef CGAL_PMP_DEBUG_FEATURE_GRAPH
    std::cout << "# constrained polylines after split_graph_into_polylines() : " << polylines.size() << std::endl;
#endif

    auto add_to_map = [](const vertex_descriptor& v,
                         std::unordered_map<vertex_descriptor, unsigned int>& map)
      {
        auto it = map.find(v);
        if (it != map.end())
          map[v]++;
        else
          map.insert({ v, 1 });
      };

    // keep only the endpoints of the polylines
    std::unordered_map<vertex_descriptor, unsigned int> endpoints;
    for (const Polyline& polyline : polylines)
    {
      add_to_map(polyline.front(), endpoints);
      add_to_map(polyline.back(), endpoints);
    }

#ifdef CGAL_PMP_DEBUG_FEATURE_GRAPH
    std::cout << "# endpoints : " << endpoints.size() << std::endl;
    std::ofstream out("endpoints.xyz");
    for (const auto& [v, n] : endpoints)
    {
      const Point_3& p = get(vpmap, v);
      out << p.x() << " " << p.y() << " " << p.z() << std::endl;
    }
    out.close();
#endif

    // build a Kd-tree of endpoints
    std::vector<vertex_descriptor> tree_points;
    for (const auto& [v, n] : endpoints)
      tree_points.push_back(v);
    Tree tree(tree_points.begin(), tree_points.end(), Splitter(), Tree_traits(vpmap));
    Distance tr_dist(vpmap);

#ifdef CGAL_PMP_DEBUG_FEATURE_GRAPH
    unsigned int nb_features_added = 0;
#endif

    std::unordered_set<vertex_descriptor> removed_from_corners;
    std::unordered_set<vertex_descriptor> to_be_snapped_later;

    // for each endpoint, find closest endpoint
    for (const auto& [endpoint, constraint_valence] : endpoints)
    {
      if (constraint_valence > 1)
        continue;

      if (removed_from_corners.find(endpoint) != removed_from_corners.end())
        continue;

      // search for close endpoints
      const Point_3& query = get(vpmap, endpoint);
      const unsigned knn = 5;
      const FT epsilon = 0.0;
      K_neighbor_search search(tree, query, knn, epsilon, true, tr_dist);

      for (auto it = search.begin(); it != search.end(); it++)
      {
        const vertex_descriptor closest_endpoint = it->first;
        if (closest_endpoint == endpoint)
          continue;

        else if (removed_from_corners.find(closest_endpoint) != removed_from_corners.end())
        {
          to_be_snapped_later.insert(endpoint);
          continue;
        }

        // if the distance is too big, skip it
        if (tr_dist.inverse_of_transformed_distance(it->second) > snap_distance)
        {
          to_be_snapped_later.insert(endpoint);
          break; //too far from source
        }

        // otherwise, find the shortest path between the endpoint and the close endpoint
        std::vector<halfedge_descriptor> halfedge_sequence;
        CGAL::shortest_path_between_two_vertices(endpoint, closest_endpoint, mesh,
                                                 std::back_inserter(halfedge_sequence));

        // connect the two endpoints
        for (const halfedge_descriptor he : halfedge_sequence)
        {
          put(ecmap, edge(he, mesh), true);

#ifdef CGAL_PMP_DEBUG_FEATURE_GRAPH
          nb_features_added++;
#endif //CGAL_PMP_DEBUG_FEATURE_GRAPH
        }

        //removed_from_corners.insert(endpoint);//keep endpoint for star-shaped connections
        removed_from_corners.insert(closest_endpoint);
        break;
      }
    }
#ifdef CGAL_PMP_DEBUG_FEATURE_GRAPH
    std::cout << "# Edges added to feature graph to connect corners : " << nb_features_added << std::endl;
    nb_features_added = 0;
#endif //CGAL_PMP_DEBUG_FEATURE_GRAPH

    // the set of endpoints is now reduced, because of connections that have been made
    // some endpoints may not have been connected, we need to snap them later

    if (to_be_snapped_later.empty())
      return;

    // Build a kd-tree of all feature vertices,
    // to snap remaining endpoints to the closest feature vertex,
    // and create a T-shaped feature graph
    tree.clear();
    std::unordered_set<vertex_descriptor> feature_vertices;
    for (const edge_descriptor e : edges(mesh))
    {
      if (get(ecmap, e))
      {
        feature_vertices.insert(source(halfedge(e, mesh), mesh));
        feature_vertices.insert(target(halfedge(e, mesh), mesh));
      }
    }
    tree.insert(feature_vertices.begin(), feature_vertices.end());

    for (const vertex_descriptor endpoint : to_be_snapped_later)
    {
      // search for close endpoints
      const Point_3& query = get(vpmap, endpoint);
      const unsigned knn = 5;
      const FT epsilon = 0.0;
      K_neighbor_search search(tree, query, knn, epsilon, true, tr_dist);

      for (auto it = search.begin(); it != search.end(); it++)
      {
        const vertex_descriptor closest_pt = it->first;
        if (closest_pt == endpoint)
          continue;

        // if the distance is too big, skip it
        if (tr_dist.inverse_of_transformed_distance(it->second) > snap_distance)
          break; //too far from source

        // otherwise, find the shortest path between the endpoint and the close endpoint
        std::vector<halfedge_descriptor> halfedge_sequence;
        CGAL::shortest_path_between_two_vertices(endpoint, closest_pt, mesh,
                                                 std::back_inserter(halfedge_sequence));

        // connect the two endpoints
        for (const halfedge_descriptor he : halfedge_sequence)
        {
          put(ecmap, edge(he, mesh), true);

#ifdef CGAL_PMP_DEBUG_FEATURE_GRAPH
          nb_features_added++;
#endif //CGAL_PMP_DEBUG_FEATURE_GRAPH
        }
        break;
      }
    }

#ifdef CGAL_PMP_DEBUG_FEATURE_GRAPH
    std::cout << "# Edges added to feature graph as T-junctions: " << nb_features_added << std::endl;
#endif //CGAL_PMP_DEBUG_FEATURE_GRAPH
  }

} // end namespace PMP
} // end namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_FEATURE_GRAPH_H

#ifndef _LFS_
#define _LFS_

#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Bbox_3.h>
#include <set>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>


template < class Kernel>
class DT3 : public CGAL::Delaunay_triangulation_3<Kernel>
{
public:
  typedef DT3<Kernel> Dt;

  typedef typename Kernel::FT            FT;
  typedef typename Kernel::Point_3       Point;
  typedef typename Kernel::Vector_3      Vector;
  typedef typename Kernel::Segment_3     Segment;
  typedef typename Kernel::Line_3        Line;
  typedef typename Kernel::Triangle_3    Triangle;
  typedef typename Kernel::Tetrahedron_3 Tetrahedron;

  typedef typename Dt::Vertex                   Vertex;
  typedef typename Dt::Vertex_handle            Vertex_handle;
  typedef typename Dt::Vertex_iterator          Vertex_iterator;
  typedef typename Dt::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename Dt::Edge                  Edge;
  typedef typename Dt::Edge_iterator         Edge_iterator;
  typedef typename Dt::Finite_edges_iterator Finite_edges_iterator;

  typedef typename Dt::Facet                   Facet;
  typedef typename Dt::Facet_iterator          Face_iterator;
  typedef typename Dt::Facet_circulator        Face_circulator;
  typedef typename Dt::Finite_facets_iterator  Finite_facets_iterator;

  typedef typename Dt::Cell                  Cell;
  typedef typename Dt::Cell_handle           Cell_handle;
  typedef typename Dt::Cell_iterator         Cell_iterator;
  typedef typename Dt::Cell_circulator       Cell_circulator;
  typedef typename Dt::Finite_cells_iterator Finite_cells_iterator;

  // build KD tree for searching K nearest neighbors
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;


public:
  DT3() {}
  ~DT3() {}

private:
  std::set<Point> m_poles;
  Tree m_tree;

public:

  Vertex_handle get_source_vertex(const Edge& edge) const
  {
    return edge.first->vertex(edge.second);
  }

  Vertex_handle get_target_vertex(const Edge& edge) const
  {
    return edge.first->vertex(edge.third);
  }

  // RENDERING
  void gl_vertex(const Point& p)
  {
    ::glVertex3d(p.x(), p.y(), p.z());
  }

  double lfs(const Point& query)
  {
    Neighbor_search search(m_tree, query, 1);
    typedef typename Neighbor_search::iterator Search_iterator;
    Search_iterator search_iterator = search.begin();
    const Point p = search_iterator->first;
    return std::sqrt(CGAL::squared_distance(query, p));
  }

  template <class InputIterator>
  void init(InputIterator begin, InputIterator end)
  {
    Dt::clear();
    Dt::insert(begin, end);
    extract_poles();
    m_tree.clear();
    m_tree.insert(m_poles.begin(), m_poles.end());
  }

  void extract_poles()
  {
    std::cout << "extract poles...";
    std::ofstream outpoles("poles.points.txt");
    m_poles.clear();

    Finite_vertices_iterator v;
    for (v = Dt::finite_vertices_begin();
      v != Dt::finite_vertices_end();
      v++)
    {
      const Point& p = v->point();

      // first pole (farthest Voronoi vertex)
      Point pole1;
      if (farthest_voronoi_point_from(v, pole1))
      {
        m_poles.insert(pole1);
        outpoles << pole1 << std::endl;

        // second pole (farthest Voronoi vertex in opposite half-space)
        Point pole2;
        if (farthest_voronoi_point_halfspace_from(v, pole1, pole2))
        {
          m_poles.insert(pole2);
          outpoles << pole2 << std::endl;
        }
      }
    }
    std::cout << "done (" << m_poles.size() << " poles)" << std::endl;
  }

  bool farthest_voronoi_point_from(Vertex_handle v,
    Point& pole)
  {
    const Point& p1 = v->point();

    // get all cells incident to v
    std::list<Cell_handle> cells;
    Dt::incident_cells(v, std::back_inserter(cells));

    if (cells.size() == 0)
      return false;

    // pick first finite cell
    FT max_sq_distance = 0.0;
    bool success = false;
    typename std::list<Cell_handle>::iterator it;
    for (it = cells.begin();
      it != cells.end();
      it++)
    {
      Cell_handle cell = *it;
      if (!Dt::is_infinite(cell))
      {
        success = true;
        Point p2 = this->dual(cell);
        FT sq_distance = CGAL::squared_distance(p1, p2);
        if (it == cells.begin() ||
          (sq_distance > max_sq_distance))
        {
          pole = p2;
          max_sq_distance = sq_distance;
        }
      }
    }
    return success;
  }

  bool farthest_voronoi_point_halfspace_from(Vertex_handle v,
    const Point& pole1,
    Point& pole2)
  {
    const Point& p1 = v->point();
    Vector vector1 = pole1 - p1;

    // get all cells incident to v
    std::list<Cell_handle> cells;
    this->incident_cells(v, std::back_inserter(cells));

    if (cells.size() == 0)
      return false;

    // pick first finite cell
    FT max_sq_distance = 0;
    bool success = false;
    typename std::list<Cell_handle>::iterator it;
    for (it = cells.begin();
      it != cells.end();
      it++)
    {
      Cell_handle cell = *it;
      if (!Dt::is_infinite(cell))
      {
        Point p2 = Dt::dual(cell);
        Vector vector2 = p2 - p1;

        // check half-space 
        if (vector2 * vector1 < 0.0)
        {
          success = true;
          FT sq_distance = CGAL::squared_distance(p1, p2);
          if (it == cells.begin() ||
            (sq_distance > max_sq_distance))
          {
            pole2 = p2;
            max_sq_distance = sq_distance;
          }
        }
      }
    }
    return success;
  }
};

namespace CGAL {

template<typename PointIterator>
CGAL::Bbox_3 bbox_3_points(PointIterator pbegin, PointIterator pend)
{
  if (pbegin == pend)
    return CGAL::Bbox_3();

  CGAL::Bbox_3 bb = (*pbegin).bbox();
  for (PointIterator it = pbegin; it != pend; ++it)
    bb = bb + (*it).bbox();

  return bb;
}

template <typename Kernel
        , typename Polyhedron
        , typename SizingField>
void compute_lfs_sizing(const Polyhedron& poly
                      , SizingField& sizing)
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

  DT3<Kernel> lfs;

  // extract poles and such
  lfs.init(poly.points_begin(), poly.points_end());

  CGAL::Bbox_3 bbox = PMP::bbox_3(poly);

  unsigned int nb_samples = static_cast<unsigned int>(1e6);

  sizing.init(bbox.xmin(), bbox.xmax(),
              bbox.ymin(), bbox.ymax(),
              bbox.zmin(), bbox.zmax(),
              nb_samples);

  std::cout << "set lfs constraints...";
  BOOST_FOREACH(vertex_descriptor v, vertices(poly))
  {
    const typename Polyhedron::Point_3& query = v->point();
    const double size = lfs.lfs(query);
    sizing.add_constraint(query, size);
  }
  std::cout << "done" << std::endl;

  std::cout << "update field...";
  sizing.update();
  std::cout << "done" << std::endl;
}


template <typename Kernel
        , typename PointIterator
        , typename SizingField>
void compute_lfs_sizing(PointIterator pbegin
                      , PointIterator pend
                      , SizingField& sizing)
{
  DT3<Kernel> lfs;

  // extract poles and such
  lfs.init(pbegin, pend);

  CGAL::Bbox_3 bbox = CGAL::bbox_3_points(pbegin, pend);

  unsigned int nb_samples = static_cast<unsigned int>(1e6);

  sizing.init(bbox.xmin(), bbox.xmax(),
              bbox.ymin(), bbox.ymax(),
              bbox.zmin(), bbox.zmax(),
              nb_samples);

  std::cout << "set lfs constraints...";
  for (PointIterator it = pbegin; it != pend; ++it)
  {
    const Kernel::Point_3& query = *it;
    const double size = lfs.lfs(query);
    sizing.add_constraint(query, size);
  }
  std::cout << "done" << std::endl;

  std::cout << "update field...";
  sizing.update();
  std::cout << "done" << std::endl;

}


}//end namespace CGAL


#endif // _LFS_

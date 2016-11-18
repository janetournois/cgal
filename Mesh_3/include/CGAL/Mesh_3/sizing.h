#ifndef _SIZING_GRID_
#define _SIZING_GRID_

#include <queue>
#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>
#include <limits>


template <class Kernel>
class Sizing_grid_node
{
public:
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;

  FT m_init_size;
  FT m_size;
  //FT m_distance;
  bool m_done;
  Point m_point;
  std::array<int, 3> m_indices;
  Sizing_grid_node* m_pRef_node;

public:
  Sizing_grid_node()
  {
    m_done = false;
    m_size = (std::numeric_limits<double>::max)();
    m_pRef_node = NULL;
    m_indices[0] = m_indices[1] = m_indices[2] = 0;
  }

  ~Sizing_grid_node() {}

  const Point& point() const { return m_point; }
  Point& point() { return m_point; }

  Sizing_grid_node* ref_node() { return m_pRef_node; }
  void ref_node(Sizing_grid_node* pRef_node) { m_pRef_node = pRef_node; }

  bool& done() { return m_done; }
  const bool& done() const { return m_done; }

  FT& size() { return m_size; }
  const FT& size() const { return m_size; }

  FT& init_size() { return m_init_size; }
  const FT& init_size() const { return m_init_size; }

  const int& i() const { return m_indices[0]; }
  const int& j() const { return m_indices[1]; }
  const int& k() const { return m_indices[2]; }
  int& i() { return m_indices[0]; }
  int& j() { return m_indices[1]; }
  int& k() { return m_indices[2]; }

  void indices(const int i,
    const int j,
    const int k)
  {
    m_indices[0] = i;
    m_indices[1] = j;
    m_indices[2] = k;
  }
};


// functor for priority queue
template<class c>
struct less_candidate_size
{
  bool operator()(const c& c1, const c& c2) const
  {
    return c1.size() > c2.size();
  }
};


template <class Kernel>
class Sizing_grid
{
public:
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef Sizing_grid_node<Kernel> Node;
  typedef typename std::pair<Point, FT> Constraint;

private:
  // grid
  boost::shared_ptr<boost::multi_array<Node, 3> > m_pppNodes;
  FT m_k;
  FT m_ds;
  FT m_dv;
  std::array<FT, 3> m_xrange;
  std::array<FT, 3> m_yrange;
  std::array<FT, 3> m_zrange;
  unsigned int m_nx, m_ny, m_nz;
  FT m_max_size;
  bool m_updated;
  std::list<Constraint> m_constraints;
  FT m_factor;

  class Candidate_size
  {
  private:
    FT m_size;
    Node* m_pNode;
    Node* m_pRef_node;

  public:

    Candidate_size(Node* pNode,
      Node* pRef_node,
      const FT k)
    {
      m_pNode = pNode;
      m_pRef_node = pRef_node;

      // size = K d(node,v) + init_size(v)
      const Point& p1 = pNode->point();
      const Point& p2 = m_pRef_node->point();
      FT distance = (FT)CGAL::sqrt(CGAL_NTS to_double(CGAL::squared_distance(p1, p2)));
      m_size = k * distance + m_pRef_node->init_size();
    }

    ~Candidate_size() {}

  public:
    Candidate_size(const Candidate_size& c)
    {
      m_pNode = c.node();
      m_pRef_node = c.ref_node();
      m_size = c.size();
    }

    Node* node() const { return m_pNode; }
    Node* ref_node() const { return m_pRef_node; }

    FT& size() { return m_size; }
    const FT& size() const { return m_size; }
  };


  typedef typename std::priority_queue<Candidate_size,
    std::vector<Candidate_size>,
    less_candidate_size<Candidate_size> > PQueue;


public:
  Sizing_grid(const FT& k, const FT& factor = 1.)
    : m_pppNodes(new boost::multi_array<Node, 3>(boost::extents[0][0][0]))
    , m_k(k)
    , m_ds(0)
    , m_dv(0)
    , m_xrange({ 0,0,0 })
    , m_yrange({ 0,0,0 })
    , m_zrange({ 0,0,0 })
    , m_nx(0)
    , m_ny(0)
    , m_nz(0)
    , m_max_size(0)
    , m_updated(false)
    , m_constraints()
    , m_factor(factor)
  {}

  ~Sizing_grid()
  {
    std::cout << "delete sizing grid..."<< std::endl;
  }

  FT& k() { return m_k; }
  const FT& k() const { return m_k; }

  template<typename MeshDomainIndex>
  FT operator()(const Point& p, const int, const MeshDomainIndex&) const
  {
    return m_factor*this->value(p);
  }

  void cleanup()
  {
    m_constraints.clear();
    m_pppNodes->resize(boost::extents[0][0][0]);
    m_nx = m_ny = m_nz = 0;
    m_updated = false;
  }

  bool alloc(const unsigned int nx,
             const unsigned int ny,
             const unsigned int nz)
  {
    // cleanup
    cleanup();

    // alloc
    m_pppNodes->resize(boost::extents[nx][ny][nz]);

    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
    return true;
  }

  const FT value(const Point& query) const
  {
    // return constant value if not updated
    if (!m_updated)
      return 1.0;

    const Node *pNode = node(query);
    if (pNode == NULL)
      return (std::numeric_limits<double>::max)(); // outside bounding box
    else
      return pNode->size();
  }

  const FT density_1d(const Point& query)
  {
    if (!m_updated)
      return 1.0;

    const FT s = size_trilinear(query);
    return 1.0 / (s*s*s); // s^(d+2)
  }
  const FT density_2d(const Point& query)
  {
    if (!m_updated)
      return 1.0;

    const FT s = size_trilinear(query);
    return 1.0 / (s*s*s*s); // s^(d+2)
  }

  const FT density_3d(const Point& query)
  {
    if (!m_updated)
      return 1.0;

    const FT s = size_trilinear(query);
    return 1.0 / (s*s*s*s*s); // s^(d+2)
  }

  const FT measure_trilinear(const Point& query,
    int dim) const
  {
    if (!m_updated)
      return 1.0;

    const FT s = size_trilinear(query);
    switch (dim)
    {
    case 3:
      return 1.0 / (s*s*s);
    case 2:
      return 1.0 / (s*s);
    default: // 1
      return 1.0 / s;
    }
  }

  // trilinear interpolation of size
  // http://en.wikipedia.org/wiki/Trilinear_interpolation
  FT size_trilinear(const Point& query) const
  {
    const FT px = query.x();
    const FT py = query.y();
    const FT pz = query.z();
    if (px >= m_xrange[0] && px < m_xrange[1] &&
      py >= m_yrange[0] && py < m_yrange[1] &&
      pz >= m_zrange[0] && pz < m_zrange[1])
    {
      FT x = ((px - m_xrange[0]) / m_xrange[2] * (FT)m_nx);
      FT y = ((py - m_yrange[0]) / m_yrange[2] * (FT)m_ny);
      FT z = ((pz - m_zrange[0]) / m_zrange[2] * (FT)m_nz);
      int xl = (int)floor(x);
      int yl = (int)floor(y);
      int zl = (int)floor(z);
      int xh = (int)ceil(x);
      int yh = (int)ceil(y);
      int zh = (int)ceil(z);
      FT xd = x - xl;
      FT yd = y - yl;
      FT zd = z - zl;
      assert(xd >= 0.0 && xd <= 1.0);
      assert(yd >= 0.0 && yd <= 1.0);
      assert(zd >= 0.0 && zd <= 1.0);
      FT i1 = node(xl, yl, zl)->size() * (1.0 - zd) + node(xl, yl, zh)->size() * zd;
      FT i2 = node(xl, yh, zl)->size() * (1.0 - zd) + node(xl, yh, zh)->size() * zd;
      FT j1 = node(xh, yl, zl)->size() * (1.0 - zd) + node(xh, yl, zh)->size() * zd;
      FT j2 = node(xh, yh, zl)->size() * (1.0 - zd) + node(xh, yh, zh)->size() * zd;
      FT w1 = i1 * (1.0 - yd) + i2 * yd;
      FT w2 = j1 * (1.0 - yd) + j2 * yd;
      return w1 * (1.0 - xd) + w2 * xd;
    }
    assert(false);
    return 1.0;
  }

  Node *node(const Point& query) const
  {
    const FT x = query.x();
    const FT y = query.y();
    const FT z = query.z();
    if (x >= m_xrange[0] && x < m_xrange[1] &&
      y >= m_yrange[0] && y < m_yrange[1] &&
      z >= m_zrange[0] && z < m_zrange[1])
    {
      int i = (int)((x - m_xrange[0]) / m_xrange[2] * (FT)m_nx);
      int j = (int)((y - m_yrange[0]) / m_yrange[2] * (FT)m_ny);
      int k = (int)((z - m_zrange[0]) / m_zrange[2] * (FT)m_nz);
      return node(i, j, k);
    }
    return NULL;
  }

  Node* node(const int i,
             const int j,
             const int k) const
  {
    if (m_pppNodes == NULL)
      return NULL;

    if (i < 0 ||
      j < 0 ||
      k < 0 ||
      i >= (int)m_nx ||
      j >= (int)m_ny ||
      k >= (int)m_nz)
      return NULL;

    return &((*m_pppNodes)[i][j][k]);
  }

  Node* neighbor(Node* n,
    unsigned int index)
  {
    if (n == NULL)
      return NULL;

    switch (index)
    {
    case 0:
      return node(n->i() - 1, n->j(), n->k());
    case 1:
      return node(n->i() + 1, n->j(), n->k());
    case 2:
      return node(n->i(), n->j() - 1, n->k());
    case 3:
      return node(n->i(), n->j() + 1, n->k());
    case 4:
      return node(n->i(), n->j(), n->k() - 1);
    case 5:
      return node(n->i(), n->j(), n->k() + 1);
    default:
      return NULL;
    }
  }

  bool init(const FT xmin,
    const FT xmax,
    const FT ymin,
    const FT ymax,
    const FT zmin,
    const FT zmax,
    const unsigned int nb_samples,
    const FT ratio = 1.3)
  {
    reset();
    if (!init_range_and_alloc(xmin, xmax, ymin, ymax, zmin, zmax, nb_samples, ratio))
      return false;
    init_positions_and_indices();
    return true;
  }

  static FT len(const Vector &v)
  {
    return (FT)CGAL::sqrt(CGAL_NTS to_double(v*v));
  }
  static FT sqlen(const Vector &v)
  {
    return v*v;
  }

  void flood(PQueue& priority_queue)
  {
    m_max_size = 0.0;
    // flood using priority queue
    while (!priority_queue.empty())
    {
      //std::cerr << priority_queue.size() << "\t";
      // pop candidate out of the queue
      Candidate_size candidate = priority_queue.top();
      priority_queue.pop();

      Node* pNode = candidate.node();
      if (pNode->done() == true)
        continue;

      pNode->ref_node(candidate.ref_node());
      pNode->size() = candidate.size();
      pNode->done() = true;
      m_max_size = std::max(m_max_size, pNode->size());

      // explore neighbors
      for (unsigned int index_neighbor = 0;
        index_neighbor < 6;
        index_neighbor++)
      {
        // TODO: change size of seeds
        Node *pNeighbor = neighbor(pNode, index_neighbor);
        if (pNeighbor != NULL &&
          pNeighbor->done() == false)
          priority_queue.push(Candidate_size(pNeighbor, pNode->ref_node(), m_k));
      }
    }
  }

  void init_pqueue(PQueue& priority_queue)
  {
    // insert sizing constraints and init priority queue
    typename std::list<Constraint>::iterator it;
    for (it = m_constraints.begin();
      it != m_constraints.end();
      it++)
    {
      const Point& p = (*it).first;
      const FT init_size = (*it).second;

      // get node at position p
      Node *pNode = node(p);
      if (pNode != NULL && !pNode->done()) // specific
      {
        pNode->point() = p;
        pNode->init_size() = init_size;
        pNode->done() = true;
        pNode->size() = init_size;
        pNode->ref_node(pNode);

        // look for min size on the boundary
        //FT min_size;
        //Boundary_vertex_handle vertex_min = boundary_vertex_min_size(p,mesh,k,s,min_size);
        //pNode->size() = min_size + c;
        //pNode->ref_node(node(vertex_min->point()));

        // insert all valid neighbors to the priority queue
        for (unsigned int index_neighbor = 0;
          index_neighbor < 6;
          index_neighbor++)
        {
          Node *pNeighbor = neighbor(pNode, index_neighbor);
          if (pNeighbor != NULL &&
            pNeighbor->done() == false)
          {
            Candidate_size candidate(pNeighbor, pNode->ref_node(), m_k);
            priority_queue.push(candidate);
          }
        }
      }
    }
  }

  bool init_range_and_alloc(const FT xmin,
    const FT xmax,
    const FT ymin,
    const FT ymax,
    const FT zmin,
    const FT zmax,
    const unsigned int nb_samples,
    const FT ratio)
  {
    m_xrange[2] = ratio*(xmax - xmin);
    m_yrange[2] = ratio*(ymax - ymin);
    m_zrange[2] = ratio*(zmax - zmin);

    if (m_xrange[2] == 0.0 ||
        m_yrange[2] == 0.0 ||
        m_zrange[2] == 0.0)
      return false;

    FT xmid = 0.5*(xmin + xmax);
    FT ymid = 0.5*(ymin + ymax);
    FT zmid = 0.5*(zmin + zmax);

    m_xrange[0] = xmid - 0.5*m_xrange[2];
    m_yrange[0] = ymid - 0.5*m_yrange[2];
    m_zrange[0] = zmid - 0.5*m_zrange[2];

    m_xrange[1] = xmid + 0.5*m_xrange[2];
    m_yrange[1] = ymid + 0.5*m_yrange[2];
    m_zrange[1] = zmid + 0.5*m_zrange[2];

    // deduce nx, ny, nz
    FT volume = m_xrange[2] * m_yrange[2] * m_zrange[2];
    m_dv = volume / (FT)nb_samples;
    m_ds = pow(m_dv, 1.0 / 3.0);

    unsigned nx = (unsigned)(m_xrange[2] / m_ds);
    unsigned ny = (unsigned)(m_yrange[2] / m_ds);
    unsigned nz = (unsigned)(m_zrange[2] / m_ds);

    // alloc (and set m_nx.. variables)
    if (!alloc(nx, ny, nz))
      return false;
    return true;
  }

  void init_positions_and_indices()
  {
    // init positions and tags
    FT x = m_xrange[0] + m_ds / 2;
    for (int i = 0; i < (signed)m_nx; i++)
    {
      FT y = m_yrange[0] + m_ds / 2;
      for (int j = 0; j < (signed)m_ny; j++)
      {
        FT z = m_zrange[0] + m_ds / 2;
        for (int k = 0; k < (signed)m_nz; k++)
        {
          Node* pNode = node(i, j, k);
          pNode->indices(i, j, k);
          pNode->point() = Point(x, y, z);
          z += m_ds;
        }
        y += m_ds;
      }
      x += m_ds;
    }
  }


  double closest(const double x)
  {
    if ((ceil(x) - x) < 0.5)
      return ceil(x);
    else
      return floor(x);
  }

  void tag_done(const bool done)
  {
    for (unsigned int i = 0; i < m_nx; i++)
    for (unsigned int j = 0; j < m_ny; j++)
    for (unsigned int k = 0; k < m_nz; k++)
      node(i, j, k)->done() = done;
  }

  // total reset
  void reset(FT k)
  {
    m_k = k;
    reset();
  }

  // reset where m_k is kept
  void reset()
  {
    cleanup();
  }

  void set_k(const FT k)
  {
    if (k != m_k)
    {
      m_k = k;
      update();
    }
  }

  void add_constraint(const Point& point,
    const FT size)
  {
    Constraint constraint(point, size);
    m_constraints.push_back(constraint);
  }

  // update sizing grid
  void update()
  {
    tag_done(false);
    PQueue priority_queue;
    init_pqueue(priority_queue);
    flood(priority_queue);
    m_updated = true;
  }

};

#endif // _SIZING_GRID_

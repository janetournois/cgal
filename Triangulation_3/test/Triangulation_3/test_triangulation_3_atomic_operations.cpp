// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Random.h>

#include <cassert>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using T3 = CGAL::Triangulation_3<K>;

int main()
{
  std::size_t nbv = 1234;

  std::cout << "CGAL::Random seed = " << CGAL::get_default_random().get_seed() << std::endl;
  CGAL::Random rng;

  T3 tr;
  for(std::size_t i = 0; i < nbv; ++i)
  {
    K::Point_3 p(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.));
    tr.insert(p);
  }

  std::size_t count_edge_flips = 0;
  std::size_t count_facet_flips = 0;

  T3::Edge e = *tr.finite_edges_begin();
  if(tr.flip(e))
    ++count_edge_flips;

  T3::Facet f = *tr.finite_facets_begin();
  if(tr.flip(f))
    ++count_facet_flips;

  std::cout << "Number of edge flips: " << count_edge_flips << std::endl;
  std::cout << "Number of facet flips: " << count_facet_flips << std::endl;

  return EXIT_SUCCESS;
}


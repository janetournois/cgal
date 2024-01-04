// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_UNIFORM_SIZING_FIELD_H
#define CGAL_UNIFORM_SIZING_FIELD_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Tetrahedral_remeshing/Sizing_field.h>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
template <class Kernel>
class Uniform_sizing_field
  : public Sizing_field<Kernel>
{
private:
  typedef Sizing_field<Kernel>        Base;
public:
  typedef typename Base::FT         FT;
  typedef typename Base::Point_3    Point_3;

  Uniform_sizing_field(const FT& size)
    : m_size(size)
  {}

  template<typename Index>
  FT operator()(const Point_3&, const int, const Index&) const
  {
    return m_size;
  }

private:
  const FT m_size;
};

}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif //CGAL_UNIFORM_SIZING_FIELD_H

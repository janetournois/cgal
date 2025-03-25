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

#ifndef CGAL_INTERNAL_PEEL_SLIVERS_H
#define CGAL_INTERNAL_PEEL_SLIVERS_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/property_maps.h>

#include <optional>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
template <typename C3t3>
bool is_peelable(const C3t3& c3t3,
                 const typename C3t3::Cell_handle ch,
                 std::array<bool, 4>& facets_on_surface)
{
  typedef typename C3t3::Triangulation::Geom_traits::FT FT;
  typedef typename C3t3::Facet Facet;

  if(!c3t3.is_in_complex(ch))
    return false;

  bool on_surface = false;
  for(int i = 0; i < 4; ++i) {
    facets_on_surface[i] = !c3t3.is_in_complex(ch->neighbor(i));
    on_surface = on_surface || facets_on_surface[i];
  }
  if(!on_surface)
    return false;

  FT area_on_surface = 0.;
  FT area_inside = 0.;
  for(int i = 0; i < 4; ++i) {
    Facet f(ch, i);
    const FT facet_area = CGAL::approximate_sqrt(c3t3.triangulation().triangle(f).squared_area());
    if(facets_on_surface[i])
      area_on_surface += facet_area;
    else
      area_inside += facet_area;
  }

  return (area_inside < 1.5 * area_on_surface);
}

template<typename C3t3>
bool facet_on_convex_hull(const typename C3t3::Facet& f, const C3t3& c3t3)
{
  for(int i = 0; i < 4; ++i)
  {
    if(c3t3.is_infinite(f.first) != c3t3.is_infinite(f.first->neighbor(f.second)))
     return true;
  }
  return false;
}

template <typename C3t3>
bool is_peelable_from_convex_hull(const C3t3& c3t3,
                                  const typename C3t3::Cell_handle ch,
                                  std::array<bool, 4>& facets_to_be_kept)
{
  typedef typename C3t3::Triangulation::Geom_traits::FT FT;
  typedef typename C3t3::Facet Facet;
  const auto& tr = c3t3.triangulation();

  if(c3t3.is_in_complex(ch))
    return false;

  FT area_kept = 0.;
  FT area_peeled = 0.;
  for(int i = 0; i < 4; ++i)
  {
    facets_to_be_kept[i] = !tr.is_infinite(ch->neighbor(i));

    Facet f(ch, i);
    const FT facet_area = CGAL::approximate_sqrt(c3t3.triangulation().triangle(f).squared_area());
    if(!facets_to_be_kept[i])
      area_peeled += facet_area;
    else
      area_kept += facet_area;
  }

  return (area_peeled < 1.5 * area_kept);
}

}// end namespace internal

template<typename C3T3>
std::size_t peel_slivers(C3T3& c3t3,
                         const typename C3T3::Triangulation::Geom_traits::FT& sliver_angle);

template<typename C3T3, typename CellSelector>
std::size_t peel_slivers(C3T3& c3t3,
                         const typename C3T3::Triangulation::Geom_traits::FT& sliver_angle,
                         const CellSelector& cell_selector)
{
  using FT = typename C3T3::Triangulation::Geom_traits::FT;
  using Cell_handle = typename C3T3::Triangulation::Cell_handle;
  using Surface_patch_index = typename C3T3::Surface_patch_index;

  auto& tr = c3t3.triangulation();

  std::size_t nb_slivers_peel = 0;
  std::vector<std::pair<Cell_handle, std::array<bool, 4> > > peelable_cells;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  FT mindh = FT(180);
#endif
  for (Cell_handle cit : c3t3.cells_in_complex())
  {
    const bool selected = get(cell_selector, cit);
    if(!selected)
      continue;

    std::array<bool, 4> facets_on_surface;

    const FT dh = min_dihedral_angle(tr, cit);
    if (dh < sliver_angle && internal::is_peelable(c3t3, cit, facets_on_surface))
      peelable_cells.push_back(std::make_pair(cit, facets_on_surface));

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    mindh = (std::min)(dh, mindh);
#endif
  }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  std::cout << "Min dihedral angle : " << mindh << std::endl;
  std::cout << "Peelable cells : " << peelable_cells.size() << std::endl;
#endif

  for (auto c_i : peelable_cells)
  {
    Cell_handle c = c_i.first;
    const std::array<bool, 4>& f_on_surface = c_i.second;

    std::optional<Surface_patch_index> patch;
    for (int i = 0; i < 4; ++i)
    {
      if (f_on_surface[i])
      {
        Surface_patch_index spi = c3t3.surface_patch_index(c, i);
        if (patch.has_value() && patch.value() != spi)
        {
          //there are 2 different patches
          patch.reset();
          break;
        }
        else
        {
          patch = spi;
        }
      }
    }
    if (!patch.has_value())
      continue;

    for (int i = 0; i < 4; ++i)
    {
      if (f_on_surface[i])
        c3t3.remove_from_complex(c, i);
      else
        c3t3.add_to_complex(c, i, patch.value());
    }

    c3t3.remove_from_complex(c);
    ++nb_slivers_peel;
  }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  mindh = FT(180);
  for (Cell_handle cit : c3t3.cells_in_complex())
  {
    if(get(cell_selector, cit))
      mindh = (std::min)(min_dihedral_angle(tr, cit), mindh);
  }
  std::cout << "Peeling done (removed " << nb_slivers_peel << " slivers, "
    << "min dihedral angle = " << mindh << ")." << std::endl;
#endif

  return nb_slivers_peel;
}

template<typename C3T3>
std::size_t peel_slivers(C3T3 & c3t3,
  const typename C3T3::Triangulation::Geom_traits::FT & sliver_angle)
{
  using Tr = typename C3T3::Triangulation;
  return peel_slivers(c3t3, sliver_angle,
    CGAL::Tetrahedral_remeshing::internal::All_cells_selected<Tr>());
}


template<typename C3T3, typename FacetsKeptIndices>
void peel_cell_from_convex_hull(C3T3& c3t3,
                                typename C3T3::Cell_handle sliver,
                                const FacetsKeptIndices& f_kept)
{
  using Surface_patch_index = typename C3T3::Surface_patch_index;
  using Cell_handle = typename C3T3::Cell_handle;
  using Vertex_handle = typename C3T3::Vertex_handle;
  using Facet = typename C3T3::Facet;

  auto& tr = c3t3.triangulation();

  // update C3t3
  std::optional<Surface_patch_index> patch;
  for(int i = 0; i < 4; ++i)
  {
    if(f_kept[i]) {
      Surface_patch_index spi = sliver->surface_patch_index(i);
      if(patch.has_value() && patch.value() != spi) {
        // there are 2 different patches
        patch.reset();
        break;
      } else {
        patch = spi;
      }
    }
  }

  if(patch.has_value())
  {
    for(int i = 0; i < 4; ++i)
    {
      if(f_kept[i])
      {
        if(!c3t3.is_in_complex(sliver, i) && patch.value() != 0)
          c3t3.add_to_complex(sliver, i, patch.value());
      }
      else // f not kept
      {
        if(c3t3.is_in_complex(sliver, i))
          c3t3.remove_from_complex(sliver, i);
      }
    }
  }
  if(c3t3.is_in_complex(sliver))
    c3t3.remove_from_complex(sliver);

  // collect infinite cells that are concerned to delete them later
  std::vector<Cell_handle> to_be_deleted = {sliver};
  for(int i = 0; i < 4; ++i)
  {
    if(!f_kept[i])
      to_be_deleted.push_back(sliver->neighbor(i));
  }

  // collect facets seen from outside the (future) cavity
  using FV = CGAL::Triple<Vertex_handle, Vertex_handle, Vertex_handle>;
  struct hash_FV
  {
    std::size_t operator()(const FV& fv) const { return hash_value(fv); }
  };

  std::unordered_map<FV, Facet, hash_FV> outer_map;
  for(int i = 0; i < 4; ++i)
  {
    Cell_handle ni = sliver->neighbor(i);
    if(f_kept[i])
    {
      outer_map.insert(
        std::make_pair(make_vertex_triple(Facet(sliver, i)),
                       Facet(ni, ni->index(sliver))));
    }
    else
    {
      CGAL_assertion(tr.is_infinite(ni));
      CGAL_assertion_code(FV fvi = make_vertex_triple(Facet(ni, ni->index(sliver))));
      CGAL_assertion(outer_map.find(fvi) == outer_map.end());

      for(int j = 0; j < 4; ++j)
      {
        Cell_handle nij = ni->neighbor(j);
        if(nij == sliver)
          continue;

        Facet f_nij(nij, nij->index(ni));
        FV fvij = make_vertex_triple(f_nij);

        auto it_ij = outer_map.find(fvij);
        if(it_ij == outer_map.end())
          outer_map.insert(std::make_pair(fvij, tr.mirror_facet(f_nij)));
        else
          outer_map.erase(it_ij);//already met, remove duplicates
      }
    }
  }

  // create the new infinite cells
  std::unordered_map<FV, Facet> inner_map;
  for(int i = 0; i < 4; ++i)
  {
    if(f_kept[i])
    {
      Cell_handle new_inf_c = (i == 0 || i == 2)
                            ? tr.tds().create_cell(tr.infinite_vertex(),
                                                   sliver->vertex((i + 1) % 4),
                                                   sliver->vertex((i + 2) % 4),
                                                   sliver->vertex((i + 3) % 4))
                            : tr.tds().create_cell(tr.infinite_vertex(),
                                                   sliver->vertex((i + 1) % 4),
                                                   sliver->vertex((i + 3) % 4),
                                                   sliver->vertex((i + 2) % 4));
      for(int j = 0; j < 4; ++j)
      {
        Facet f_inf_j(new_inf_c, j);
        FV fv_inf_j = make_vertex_triple(f_inf_j);

        auto itout = outer_map.find(fv_inf_j);
        if(itout == outer_map.end())
        {
          auto itin = inner_map.find(fv_inf_j);
          if(itin == inner_map.end())
            inner_map.insert(std::make_pair(fv_inf_j, f_inf_j));
          else
          {
            Facet f_inf_j_from_inside = itin->second;
            f_inf_j_from_inside.first->set_neighbor(f_inf_j_from_inside.second, new_inf_c);
            new_inf_c->set_neighbor(j, f_inf_j_from_inside.first);
          }
        }
        else
        {
          Facet f_inf_j_from_outside = itout->second;
          f_inf_j_from_outside.first->set_neighbor(f_inf_j_from_outside.second, new_inf_c);
          new_inf_c->set_neighbor(j, f_inf_j_from_outside.first);
        }
      }
    }
  }

  // delete cells
  tr.tds().delete_cells(to_be_deleted.begin(), to_be_deleted.end());
  CGAL_assertion(tr.tds().is_valid(true));
}


template<typename C3T3>
std::size_t peel_slivers_from_convex_hull(C3T3& c3t3,
  const typename C3T3::Triangulation::Geom_traits::FT& sliver_angle)
{
  using Tr = typename C3T3::Triangulation;
  using FT = typename Tr::Geom_traits::FT;
  using Vertex_handle = typename Tr::Vertex_handle;
  using Cell_handle = typename Tr::Cell_handle;
  using Facet = typename Tr::Facet;
  using Surface_patch_index = typename C3T3::Surface_patch_index;

  auto& tr = c3t3.triangulation();

  std::size_t nb_slivers_peel = 0;
  std::vector<std::pair<Cell_handle, std::array<bool, 4> >> peelable_cells;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  FT mindh = FT(180);
#endif

  std::vector<Cell_handle> ch_cells;
  tr.incident_cells(tr.infinite_vertex(), std::back_inserter(ch_cells));

  for(const Cell_handle chc : ch_cells)
  {
    const Cell_handle c = chc->neighbor(chc->index(tr.infinite_vertex()));
    if(c3t3.is_in_complex(c))
      continue;// will be treated later by remeshing

    const FT dh = min_dihedral_angle(tr, c);
    if(dh > sliver_angle)
      continue;

    std::array<bool, 4> facets_to_be_kept;
    if(internal::is_peelable_from_convex_hull(c3t3, c, facets_to_be_kept))
      peelable_cells.push_back(std::make_pair(c, facets_to_be_kept));

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    mindh = (std::min)(dh, mindh);
#endif
  }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  std::cout << "Min dihedral angle : " << mindh << std::endl;
  std::cout << "Peelable cells : " << peelable_cells.size() << std::endl;
#endif

  for(auto c_i : peelable_cells)
  {
    Cell_handle sliver = c_i.first;
    const std::array<bool, 4>& f_kept = c_i.second;
    peel_cell_from_convex_hull(c3t3, sliver, f_kept);
    ++nb_slivers_peel;
  }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  std::cout << "Peeling done (removed " << nb_slivers_peel << " slivers)." << std::endl;
#endif

  return nb_slivers_peel;
}

} // end namespace Tetrahedral_remeshing
} // end namespace CGAL

#endif // CGAL_INTERNAL_PEEL_SLIVERS_H

// Copyright (C) 2007-2019 Garth N. Wells
//
// This file is part of DOLFIN (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#pragma once

#include <array>
#include <dolfin/la/SparsityPattern.h>
#include <dolfin/fem/MultiPointConstraint.h>

namespace dolfin
{
namespace la
{
class SparsityPattern;
}

namespace mesh
{
class Mesh;
}

namespace fem
{
class DofMap;
 class MultiPointConstraint;

/// This class provides functions to compute the sparsity pattern
/// based on DOF maps

class SparsityPatternBuilder
{
public:
  /// Iterate over cells and insert entries into sparsity pattern
  static void cells(la::SparsityPattern& pattern, const mesh::Mesh& mesh,
                    const std::array<const fem::DofMap*, 2> dofmaps);

  /// Iterate over interior facets and insert entries into sparsity pattern
  static void interior_facets(la::SparsityPattern& pattern,
                              const mesh::Mesh& mesh,
                              const std::array<const fem::DofMap*, 2> dofmaps);

  /// Iterate over exterior facets and insert entries into sparsity pattern
  static void exterior_facets(la::SparsityPattern& pattern,
                              const mesh::Mesh& mesh,
                              const std::array<const fem::DofMap*, 2> dofmaps);

  // Iterate over the master-slave cells and insert entries into sparsity pattern
  /// @param[in] pattern The sparsity pattern which is modified
  /// @param[in] mesh The relevant mesh
  /// @param[in] dofmaps, the dofmaps for the test and trial function (assuming same sapce
  /// @param[in] mpc the multipointconstraint
  static void MultiPointConstraint(la::SparsityPattern& pattern,
								   const mesh::Mesh& mesh,
								   const std::array<const fem::DofMap*, 2> dofmaps,
								   fem::MultiPointConstraint& mpc);
};
} // namespace fem
} // namespace dolfin

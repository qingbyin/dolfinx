// Copyright (C) 2020 Garth N. Wells
//
// This file is part of DOLFINX (https://www.fenicsproject.org)
//
// SPDX-License-Identifier:    LGPL-3.0-or-later

#pragma once

#include <Eigen/Dense>
#include <dolfinx/common/types.h>
#include <dolfinx/fem/DofMap.h>
#include <dolfinx/fem/FiniteElement.h>
#include <dolfinx/la/PETScVector.h>
#include <functional>
#include <memory>
#include <petscsys.h>
#include <petscvec.h>
#include <string>
#include <vector>

namespace dolfinx::function
{
/// Interpolate a finite element Function into this function space,
/// filling the array of expansion coefficients associated with this
/// function space
/// @param[in,out] u The expansion coefficients. It must be correctly
///   sized by the calling function.
/// @param[in] v The function to be interpolated
void interpolate(Function& u, const Function& v)
{
  u.function_space
}

// /// Interpolate an expression into this function space, filling the
// /// array of expansion coefficients associated with this function
// /// space.
// /// @param[in,out] u The expansion coefficients. It must be correctly
// ///   sized by the calling function.
// /// @param[in] f The function to be interpolated
// template <typename T>
// void interpolate(
//     Function& u,
//     const std::function<
//         Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(
//             const Eigen::Ref<const Eigen::Array<double, 3, Eigen::Dynamic,
//                                                 Eigen::RowMajor>>&)>& f);

// /// Interpolate an expression into this function space, filling the
// /// array of expansion coefficients associated with this function
// /// space.
// /// @note This interface is not intended for general use. It supports
// ///   the use of an expression function with a C-signature; it is
// ///   typically used by compiled Numba functions with C interface.
// /// @param[in,out] u The expansion coefficients to be filled.
// ///   It must be correctly sized by the calling function.
// /// @param[in] f The function to be interpolated
// template <typename T>
// void interpolate_c(
//     Function& u,
//     const std::function<void(
//         Eigen::Ref<
//             Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>,
//         const Eigen::Ref<const Eigen::Array<double, Eigen::Dynamic, 3,
//                                             Eigen::RowMajor>>&)>& f);

// void interpolate_from_any(Function& u, const Function& v);

// void interpolate_private(
//     Eigen::Ref<Eigen::Array<PetscScalar, Eigen::Dynamic, 1>> coefficients,
//     const Eigen::Ref<const Eigen::Array<
//         PetscScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>& values);

// inline void interpolate(Function& u, const Function& v)
// {
//   auto mesh = u.function_space()->mesh();
//   auto element = u.function_space()->element();
//   assert(mesh);
//   assert(element);

//   // Check that function ranks match
//   if (int rank_v = v.function_space()->element()->value_rank();
//       element->value_rank() != rank_v)
//   {
//     throw std::runtime_error("Cannot interpolate function into function space. "
//                              "Rank of function ("
//                              + std::to_string(rank_v)
//                              + ") does not match rank of function space ("
//                              + std::to_string(element->value_rank()) + ")");
//   }

//   // Check that function dimension match
//   for (int i = 0; i < element->value_rank(); ++i)
//   {
//     if (int v_dim = v.function_space()->element()->value_dimension(i);
//         element->value_dimension(i) != v_dim)
//     {
//       throw std::runtime_error(
//           "Cannot interpolate function into function space. "
//           "Dimension "
//           + std::to_string(i) + " of function (" + std::to_string(v_dim)
//           + ") does not match dimension " + std::to_string(i)
//           + " of function space(" + std::to_string(element->value_dimension(i))
//           + ")");
//     }
//   }

//   // interpolate_from_any(expansion_coefficients, v);
// }

// template <typename T>
// void interpolate(
//     Function& u,
//     const std::function<
//         Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(
//             const Eigen::Ref<const Eigen::Array<double, 3, Eigen::Dynamic,
//                                                 Eigen::RowMajor>>&)>& f)
// {
//   // Evaluate expression at dof points
//   const Eigen::Array<double, 3, Eigen::Dynamic, Eigen::RowMajor> x
//       = u.function_space()->tabulate_dof_coordinates().transpose();
//   const Eigen::Array<PetscScalar, Eigen::Dynamic, Eigen::Dynamic,
//                      Eigen::RowMajor>
//       values = f(x);

//   auto element = u.function_space()->element();
//   assert(element);
//   std::vector<int> vshape(element->value_rank(), 1);
//   for (std::size_t i = 0; i < vshape.size(); ++i)
//     vshape[i] = element->value_dimension(i);
//   const int value_size = std::accumulate(std::begin(vshape), std::end(vshape),
//                                          1, std::multiplies<>());

//   // Note: pybind11 maps 1D NumPy arrays to column vectors for
//   // Eigen::Array<PetscScalar, Eigen::Dynamic,Eigen::Dynamic, Eigen::RowMajor>
//   // types, therefore we need to handle vectors as a special case.
//   if (values.cols() == 1 and values.rows() != 1)
//   {
//     if (values.rows() != x.cols())
//     {
//       throw std::runtime_error("Number of computed values is not equal to the "
//                                "number of evaluation points. (1)");
//     }
//     interpolate(u, values);
//   }
//   else
//   {
//     if (values.rows() != value_size)
//       throw std::runtime_error("Values shape is incorrect. (2)");
//     if (values.cols() != x.cols())
//     {
//       throw std::runtime_error("Number of computed values is not equal to the "
//                                "number of evaluation points. (2)");
//     }
//     interpolate(u, values.transpose());
//   }
// }

// template <typename T>
// void interpolate_c(
//     Function& u,
//     const std::function<void(
//         Eigen::Ref<
//             Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>,
//         const Eigen::Ref<const Eigen::Array<double, Eigen::Dynamic, 3,
//                                             Eigen::RowMajor>>&)>& f)
// {
//   // Build list of points at which to evaluate the Expression
//   const Eigen::Array<double, Eigen::Dynamic, 3, Eigen::RowMajor> x
//       = u.function_space()->tabulate_dof_coordinates();

//   // Evaluate expression at points
//   assert(_element);
//   std::vector<int> vshape(_element->value_rank(), 1);
//   for (std::size_t i = 0; i < vshape.size(); ++i)
//     vshape[i] = _element->value_dimension(i);
//   const int value_size = std::accumulate(std::begin(vshape), std::end(vshape),
//                                          1, std::multiplies<>());
//   Eigen::Array<PetscScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
//       values(x.rows(), value_size);
//   f(values, x);

//   interpolate(u, values);
// }

// void FunctionSpace::interpolate_from_any(Function& u, const Function& v)
// {
//   assert(v.function_space());
//   if (!v.function_space()->has_element(*u.function_space()->element()))
//   {
//     throw std::runtime_error("Restricting finite elements function in "
//                              "different elements not supported.");
//   }

//   assert(u.function_space()->mesh());
//   assert(v.function_space()->mesh());
//   if (u.function_space()->mesh()->id() != v.function_space()->mesh()->id())
//   {
//     throw std::runtime_error(
//         "Interpolation on different meshes not supported (yet).");
//   }

//   const int tdim = u.function_space()->mesh()->topology().dim();

//   // Get dofmaps
//   assert(v.function_space());
//   std::shared_ptr<const fem::DofMap> dofmap_v = v.function_space()->dofmap();
//   assert(dofmap_v);
//   auto map = _mesh->topology().index_map(tdim);
//   assert(map);

//   // Iterate over mesh and interpolate on each cell
//   assert(_dofmap);
//   la::VecReadWrapper v_vector_wrap(v.vector().vec());
//   Eigen::Map<const Eigen::Matrix<PetscScalar, Eigen::Dynamic, 1>> v_array
//       = v_vector_wrap.x;
//   const int num_cells = map->size_local() + map->num_ghosts();
//   for (int c = 0; c < num_cells; ++c)
//   {
//     auto dofs_v = dofmap_v->cell_dofs(c);
//     auto cell_dofs = _dofmap->cell_dofs(c);
//     assert(dofs_v.size() == cell_dofs.size());
//     for (Eigen::Index i = 0; i < dofs_v.size(); ++i)
//       expansion_coefficients[cell_dofs[i]] = v_array[dofs_v[i]];
//   }
// }
//-----------------------------------------------------------------------------

} // namespace dolfinx::function

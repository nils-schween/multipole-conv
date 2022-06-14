/*
 * Copyright (C) 2022 Schween, Nils W. and Reville, Brian
 *
 * This file is part of multipole-conv.
 *
 * multipole-conv is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * multipole-conv is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with multipole-conv. If not, see <https://www.gnu.org/licenses/>.
 */

#include "matrix_operations.h"

#include <complex>
#include <cstddef>

#include "matrix.h"

using std::size_t;
namespace multipole_conv {
Matrix<double> inverse_diagonal(const Matrix<double>& mat) {
  // TODO: no error handling, namely it is not checked if mat is a square
  // matrix (mat.rows() == mat.columns())
  Matrix<double> res(mat.rows());
  for (size_t i = 0; i < mat.rows(); ++i) res(i, i) = 1. / mat(i, i);
  return res;
}

// Convert a matrix with double entries to a matrix complex<double> entries
Matrix<std::complex<double>> convert_real_to_complex(
    const Matrix<double>& mat) {
  Matrix<std::complex<double>> cplx_mat(mat.rows(), mat.columns());
  for (size_t i = 0; i < mat.rows(); ++i)
    for (size_t j = 0; j < mat.columns(); ++j)
      cplx_mat(i, j) = mat(i, j);
  return cplx_mat;
}
}  // namespace multipole_conv

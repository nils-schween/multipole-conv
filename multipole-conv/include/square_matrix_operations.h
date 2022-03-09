#ifndef SQUARE_MATRIX_OPERATIONS_H
#define SQUARE_MATRIX_OPERATIONS_H

#include "square_matrix.h"
namespace multipole_conv {
// Computes the inverse of a diagonal square matrix
template <typename T>
SquareMatrix<T> inverse_diagonal(const SquareMatrix<T>& mat) {
  SquareMatrix<T> res(mat.dim());
  for (std::size_t i = 0; i < mat.dim(); ++i) res(i, i) = 1. / mat(i, i);
  return res;
}
}  // namespace mutlipole_conv
#endif

#ifndef SQUARE_MATRIX_OPERATIONS_H
#define SQUARE_MATRIX_OPERATIONS_H

#include "matrix.h"
namespace multipole_conv {
// Computes the inverse of a diagonal square matrix
template <typename T>
Matrix<T> inverse_diagonal(const Matrix<T>& mat) {
  // TODO: no error handling, namely it is not checked if mat is a square
  // matrix (mat.rows() == mat.columns())
  Matrix<T> res(mat.rows());
  for (std::size_t i = 0; i < mat.dim(); ++i) res(i, i) = 1. / mat(i, i);
  return res;
}
}  // namespace mutlipole_conv
#endif

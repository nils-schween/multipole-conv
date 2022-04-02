#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include "matrix.h"

namespace multipole_conv {
// Computes the inverse of a diagonal square matrix
Matrix<double> inverse_diagonal(const Matrix<double>& mat);

// Convert a matrix with double entries to a matrix complex<double> entries
Matrix<std::complex<double>> convert_real_to_complex(const Matrix<double>& mat);
}  // namespace multipole_conv
#endif

#ifndef MULTIPOLE_CONV_FUNCS_H
#define MULTIPOLE_CONV_FUNCS_H

#include <complex>
#include <cstddef>

#include "square_matrix.h"

namespace multipole_conv {
// Compute the normalisations of all the real spherical harmonics of degree l
SquareMatrix<double> norms_real_sph(std::size_t degree);

// Compute the normalisations of all the complex spherical harmonics of degree l
SquareMatrix<double> norms_complex_sph(std::size_t degree);

// Compute the transformation matrix from the real spherical harmonics to the
// complex spherical harmonics
SquareMatrix<std::complex<double>> real_sph_to_complex(std::size_t degree);

// Compute the transformation matrix from the complex spherical harmonics to the
// real spherical harmonics
SquareMatrix<std::complex<double>> complex_sph_to_real(std::size_t degree);

// Computes the coefficient in the basis transformation from the real spherical
// harmonics without normalisation to multipole basis functions.
// p_index = 0 belongs to the basis functions M_0qr
// p_index = 1 belongs to the basis functions M_1qr
double coefficient(std::size_t degree, std::size_t order, bool p,
                   std::size_t sum_start);

// Computes the basis transformation between the real spherical harmonics
// without normalisation and multipole basis function
SquareMatrix<double> basis_transformation(std::size_t degree);

// Create a permutation matrix which, when multiplied with the basis
// transformation matrix, results in a matrix with four blocks of triangular
// matrices
SquareMatrix<double> permutation(std::size_t degree);

// Invert the basis transformation matrix
SquareMatrix<double> invert_basis_transformation(
    const SquareMatrix<double>& trans_mat);
}  // namespace multipole_conv
#endif

#ifndef MULTIPOLE_CONV_FUNCS_H
#define MULTIPOLE_CONV_FUNCS_H

#include <cmath>
#include <complex>
#include <cstddef>

#include "math.h"
#include "matrix.h"
#include "options.h"

namespace multipole_conv {
// Compute the Condon Shortley Phase for all the spherical harmonics of degree l
Matrix<double> condon_shortley_phase(std::size_t degree);

// Compute the normalisations of all the real spherical harmonics of degree l
Matrix<double> norms_real_sph(std::size_t degree);

// Compute the normalisations of all the complex spherical harmonics of degree l
Matrix<double> norms_complex_sph(std::size_t degree);

// Use options to determine if you norms_real_sph or norms_complex_sph is called
Matrix<double> norms_sph(size_t degree, MPOptions options);

// Due to the Addition theoreom of spherical harmonics, there are additional
// factors in the expansion of the potential in terms of spherical harmonics:
// The factor (2l + 1)/4pi for (real/complex) spherical harmonics with
// normalisation. And the factor (2/(1 + \delta_0m)*(l - m)!/(l + m)! for real
// spherical harmonics without normalisation.
Matrix<double> addition_theorem_factor(std::size_t degree, MPOptions options);

// Compute the transformation matrix from the real spherical harmonics to the
// complex spherical harmonics
Matrix<std::complex<double>> real_sph_to_complex(std::size_t degree);

// Compute the transformation matrix from the complex spherical harmonics to the
// real spherical harmonics
Matrix<std::complex<double>> complex_sph_to_real(std::size_t degree);

// Computes the coefficient in the basis transformation from the real spherical
// harmonics without normalisation to multipole basis functions.
// p_index = 0 belongs to the basis functions M_0qr
// p_index = 1 belongs to the basis functions M_1qr
double coefficient(std::size_t degree, std::size_t order, size_t s_index,
                   std::size_t sum_start);

// Computes the basis transformation between the real spherical harmonics
// without normalisation and multipole basis function
Matrix<double> basis_transformation(std::size_t degree);

// Create a permutation matrix which, when multiplied with the basis
// transformation matrix, results in a matrix with four blocks of triangular
// matrices
Matrix<double> permutation(std::size_t degree);

// Invert the basis transformation matrix
Matrix<double> invert_basis_transformation(const Matrix<double>& trans_mat);

// Compute the dependent components of Cartesian multipole moment
Matrix<double> dependent_components(Matrix<double>& multipole_basis_mat);
}  // namespace multipole_conv
#endif

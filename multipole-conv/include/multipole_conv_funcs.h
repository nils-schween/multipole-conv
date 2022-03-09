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
SquareMatrix<std::complex<double>> transform_real_sph_to_complex(
    std::size_t degree);

// Compute the transformation matrix from the complex spherical harmonics to the
// real spherical harmonics
SquareMatrix<std::complex<double>> transform_complex_sph_to_real(
    std::size_t degree);

// Computes the common coefficient in front of the multipole basis functions.
// p_index = 0 belongs to the basis functions M_0qr
// p_index = 1 belongs to the basis functions M_1qr
double coefficient(std::size_t limit, std::size_t index, bool p);
}  // namespace multipole_conv
#endif

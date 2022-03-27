#ifndef MULTIPOLE_CONV_FUNCS_H
#define MULTIPOLE_CONV_FUNCS_H

#include <cmath>
#include <complex>
#include <cstddef>

#include "math.h"
#include "matrix.h"

namespace multipole_conv {
// Compute the Condon Shortley Phase for all the spherical harmonics of degree l
Matrix<double> condon_shortley_phase(std::size_t degree);

// Compute the normalisations of all the real spherical harmonics of degree l
Matrix<double> norms_real_sph(std::size_t degree);

// Compute the normalisations of all the complex spherical harmonics of degree l
Matrix<double> norms_complex_sph(std::size_t degree);

// Due to the Addition theoreom of spherical harmonics, there are additional
// factors in the expansion of the potential in terms of spherical harmonics:
// The factor (2l + 1)/4pi for (real/complex) spherical harmonics with
// normalisation. And the factor 2*(l - m)!/(l + m)! for real spherical
// harmonics without normalisation.
Matrix<double> johnston_factor(std::size_t degree);
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

// Compute the dependent components of Cartesian multipole moment. Since the
// multipole basis may be expressed as a linear combination of real or complex
// spherical harmonics, I need to use a template here. And the definiton of the
// template must go in the header file. This breaks with the previous
// programming style (separation of declaration and definiton). Maybe there is a
// better way to solve this.
template <typename T>  // T either double or std::complex<double>
Matrix<T> dependent_components(Matrix<T>& multipole_basis_mat) {
  std::size_t degree = (multipole_basis_mat.rows() - 1) / 2;
  Matrix<double> dep_comps(((degree - 1) * degree) / 2, 2 * degree + 1);
  std::size_t dep_comp_idx = 0;
  for (std::size_t order = degree, i = 0; i < (degree - 1); ++i, --order) {
    for (std::size_t p_idx = order; p_idx > 1; --p_idx) {
      // To each combination order,p_index belongs one dependent component. The
      // dependent components are a linear combination of the multipole basis
      // functions (which are represented) by the rows of the
      // multipole_basis_mat (pattern of dependence is Pascal's triangle).

      // Loop through the correspoding layer of Pascal's triangle
      // Distinguish between even and odd p_index
      int sign = (p_idx % 2 ? -1 : 1);
      size_t steps = std::floor(p_idx / 2);  // steps in Pascal's triangle
      // NOTE: steps + 1 = number of functions on which the dependent component
      // depends on
      int sign_in_pascals_triangle = (steps % 2 ? -1 : 1);
      for (std::size_t row_idx = (p_idx % 2) * 2 * degree + sign * i, j = 0;
           j <= steps; ++j, row_idx += sign * 2) {
        // add up the rows of multipole_basis_mat
        for (size_t k = 0; k < 2 * degree + 1; ++k) {
          dep_comps(dep_comp_idx, k) += sign_in_pascals_triangle *
                                        binomial(steps, j) *
                                        multipole_basis_mat(row_idx, k);
        }
      }
      ++dep_comp_idx;
    }
  }
  return dep_comps;
}

}  // namespace multipole_conv
#endif

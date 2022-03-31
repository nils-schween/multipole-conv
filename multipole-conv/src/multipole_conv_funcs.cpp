#include "multipole_conv_funcs.h"

#include <cmath>
#include <cstddef>
#include <vector>

#include "constants.h"
#include "math.h"
#include "matrix.h"
#include "options.h"

using std::size_t;

multipole_conv::Matrix<double> multipole_conv::condon_shortley_phase(
    size_t degree) {
  Matrix<double> csp(2 * degree + 1);
  for (size_t order = degree, i = 0; i < degree; ++i, --order) {
    int sign = order % 2 ? -1 : 1;
    csp(i, i) = sign;
    csp(2 * degree - i, 2 * degree - i) = sign;
  }
  csp(degree, degree) = 1;
  return csp;
}

multipole_conv::Matrix<double> multipole_conv::norms_real_sph(size_t degree) {
  Matrix<double> res{2 * degree + 1};
  double factor = {std::sqrt((2 * degree + 1) / (2 * pi))};
  for (size_t order = degree, i = 0; i < degree; ++i, --order) {
    double norm_degree_order =
        factor *
        std::sqrt((factorial(degree - order) / factorial(degree + order)));
    res(i, i) = norm_degree_order;
    res(2 * degree - i, 2 * degree - i) = norm_degree_order;
  }
  res(degree, degree) = factor / std::sqrt(2);  // (1 + delta_0m) included
  return res;
}

multipole_conv::Matrix<double> multipole_conv::norms_complex_sph(
    size_t degree) {
  Matrix<double> res{2 * degree + 1};
  double factor = {std::sqrt((2 * degree + 1) / (4 * pi))};
  for (size_t order = degree, i = 0; i < degree; ++i, --order) {
    // positive order (m > 0)
    res(i, i) =
        factor *
        std::sqrt((factorial(degree - order) / factorial(degree + order)));
    // negative order (m < 0)
    res(2 * degree - i, 2 * degree - i) =
        factor *
        std::sqrt((factorial(degree + order) / factorial(degree - order)));
    //  m = 0
    res(degree, degree) = factor;
  }
  return res;
}

multipole_conv::Matrix<double> multipole_conv::addition_theorem_factor(
    std::size_t degree, multipole_conv::MPOptions options) {
  Matrix<double> add_thm_factor(2 * degree + 1);
  double factor = 1.;
  if ((options & MPOptions::normalisation) != MPOptions::none) {
    // If splitted, take the square of the factor
    if ((options & MPOptions::split_addition_theroem) != MPOptions::none)
      factor = std::sqrt((2 * degree + 1) / (4 * pi));
    else
      factor = (2 * degree + 1) / (4 * pi);
    // fill diagonal with factor
    for (size_t i = 0; i < 2 * degree + 1; ++i) add_thm_factor(i, i) = factor;

  } else {
    for (size_t order = degree, i = 0; i < degree; ++i, --order) {
      if ((options & MPOptions::split_addition_theroem) != MPOptions::none)
        factor = std::sqrt(2 * factorial(degree - order) /
                           factorial(degree + order));
      else
        factor = 2 * factorial(degree - order) / factorial(degree + order);

      add_thm_factor(i, i) = factor;
      add_thm_factor(2 * degree - i, 2 * degree - i) = add_thm_factor(i, i);
    }
    add_thm_factor(degree, degree) = 1;
  }
  return add_thm_factor;
}

multipole_conv::Matrix<std::complex<double>>
multipole_conv::real_sph_to_complex(size_t degree) {
  using namespace std::complex_literals;
  Matrix<std::complex<double>> transformation_matrix{2 * degree + 1};
  double factor = 1 / std::sqrt(2);
  for (size_t order = degree, i = 0; i < degree; ++i, --order) {
    // (-1)^order term
    double sign = (order % 2) ? -1 : 1;
    // Upper part of the transformation matrix (s = 0)
    transformation_matrix(i, i) = factor;
    transformation_matrix(i, 2 * degree - i) = factor * 1.i;
    // Lower part of the transformation matrix (s = 1)
    transformation_matrix(2 * degree - i, i) = sign * factor;
    transformation_matrix(2 * degree - i, 2 * degree - i) =
        -sign * factor * 1.i;
  }
  // m = 0 and s = 0 term
  transformation_matrix(degree, degree) = 1.;
  return transformation_matrix;
}

multipole_conv::Matrix<std::complex<double>>
multipole_conv::complex_sph_to_real(size_t degree) {
  // The transformation between the real and the complex spherical harmonics is
  // a unitary transformation, i.e its inverse is the transpose conjugate of the
  // corresponding transformation matrix
  Matrix<std::complex<double>> transformation_matrix =
      real_sph_to_complex(degree).transpose().conjugate();
  return transformation_matrix;
}

double multipole_conv::coefficient(size_t degree, size_t order, size_t s_index,
                                   size_t sum_start) {
  double sign = 1.;  // (-1)^(f(order))
  size_t sum_limit = 0;
  unsigned int p = 0;
  if (order % 2 == 0 && s_index == 1) {  // m even, s = 1
    sum_limit = order / 2 - 1;
    sign = std::pow(-1., order / 2 + 1);
    p = 1;
  } else if (order % 2 == 0 && s_index == 0) {  // m even, s = 0
    sum_limit = order / 2;
    sign = std::pow(-1., sum_limit);
    p = 0;
  } else if (order % 2 == 1 && s_index == 1) {  // m odd, s = 1
    sum_limit = (order - 1) / 2;
    sign = std::pow(-1., (order + 1) / 2);
    p = 0;
  } else {  // m odd, s = 0
    sum_limit = (order - 1) / 2;
    sign = std::pow(-1., (order + 1) / 2);
    p = 1;
  }
  double res = 0.;
  for (size_t i = sum_start; i <= sum_limit; ++i) {
    res += multipole_conv::binomial(order, 2 * i + p) *
           multipole_conv::binomial(i, sum_start);
  }
  double factor = 1 / factorial(degree - order);
  res *= sign * factor;
  return res;
}

multipole_conv::Matrix<double> multipole_conv::basis_transformation(
    size_t degree) {
  // NOTE: The implementation of of this function strongly depends on the chosen
  // ordering of the real spherical harmonics and the multipole basis functions
  Matrix<double> basis_transformation{2 * degree + 1};
  for (size_t order = degree, i = 0; i < degree; ++i, --order) {
    int sign_upper_part = (order % 2 ? -1 : 1);
    int sign_lower_part = (order % 2 ? 1 : -1);
    for (size_t k = 0; k <= std::floor(order / 2); ++k) {
      // Upper part (s = 0)
      basis_transformation(
          i, (degree - sign_upper_part * order) + sign_upper_part * 2 * k) =
          coefficient(degree, order, 0, k);
      // Lower part (s = 1)
      if (order % 2 == 0 && k == order / 2)
        continue;  // the sum for even m ends at k = order/2 - 1
      basis_transformation(2 * degree - i, (degree - sign_lower_part * order) +
                                               sign_lower_part * 2 * k) =
          coefficient(degree, order, 1, k);
    }
    // m = 0 and s = 0
    basis_transformation(degree, degree) = coefficient(degree, 0, 0, 0);
  }
  return basis_transformation;
}

multipole_conv::Matrix<double> multipole_conv::permutation(size_t degree) {
  Matrix<double> permutation_mat(2 * degree + 1);
  // l
  // upper part: r^l Y_ll0, r^l Y_ll-20 ...
  // lower part: r^l Y_ll1, r^l Y_ll-21 ...
  size_t num_permutations_l = std::floor(degree / 2);
  permutation_mat(0, 0) = 1.;                    // first row upper part (s = 0)
  permutation_mat(degree + 1, 2 * degree) = 1.;  // firsr row lower part (s = 1)
  for (size_t i = 1; i <= num_permutations_l; ++i) {
    permutation_mat(i, 2 * i) = 1.;  // upper part
    // in the lower part (s = 1) Y_l01 does not exist
    if (degree % 2 == 0 && i == num_permutations_l) continue;
    permutation_mat(degree + 1 + i, 2 * degree - 2 * i) = 1.;  // lower part
  }
  // l-1
  // upper part: r^l Y_l-1l-10, r^l Y_l-1l-30 ...
  // lower part: r^l Y_l-1l-11, r^l Y_l-1l-31 ...
  size_t num_permutations_l_minus_one = std::floor((degree - 1) / 2);
  // first row upper part
  permutation_mat(num_permutations_l + 1, 1) = 1.;
  // first row lower part
  permutation_mat(degree + 1 + num_permutations_l_minus_one + 1,
                  2 * degree - 1) = 1.;
  for (size_t i = 1; i <= num_permutations_l_minus_one; ++i) {
    permutation_mat(num_permutations_l + 1 + i, 2 * i + 1) = 1.;  // upper part
    // in the lower part (s = 1) Y_l-101 does not exist
    if (degree % 2 == 1 && i == num_permutations_l_minus_one) continue;
    permutation_mat(degree + 1 + num_permutations_l_minus_one + 1 + i,
                    2 * degree - 1 - 2 * i) = 1.;  // lower part
  }
  return permutation_mat;
}

// Auxiliary functions for invert_basis_transformation
void inv_upper_triangular_mat(size_t row_idx, size_t column_idx, size_t size,
                              const multipole_conv::Matrix<double>& mat,
                              multipole_conv::Matrix<double>& inverse) {
  size_t last_row_idx = row_idx + size - 1;
  size_t last_column_idx = column_idx + size - 1;

  // Position of the inverted triangular matrix in the inverse matrix
  size_t inv_row_idx = column_idx;
  size_t inv_column_idx = row_idx;
  size_t inv_last_row_idx = last_column_idx;
  size_t inv_last_column_idx = last_row_idx;

  // Backward substitution (for "size" times unit vectors)
  // last row of the upper triangular matrix corresponds to one non zero entry
  // in the inverse
  inverse(inv_last_row_idx, inv_last_column_idx) =
      1. / mat(last_row_idx, last_column_idx);
  for (size_t i = 0; i < size; ++i) {  // iterating over rhs = e_1, ... , e_size
    // k = 1 because we already dealt with last row of the upper triangular
    // matrix
    for (size_t k = 1; k < size; ++k) {
      double sum = 0.;
      for (size_t j = 0; j < k; ++j) {
        sum += mat(last_row_idx - k, last_column_idx - j) *
               inverse(inv_last_row_idx - j, inv_column_idx + i);
      }
      inverse(inv_last_row_idx - k, inv_column_idx + i) =
          (((last_row_idx - k == row_idx + i) ? 1. : 0) - sum) /
          mat(last_row_idx - k, last_column_idx - k);
    }
  }
}

multipole_conv::Matrix<double> multipole_conv::invert_basis_transformation(
    const Matrix<double>& trans_mat) {
  Matrix<double> inverse(trans_mat.rows());
  size_t degree = (trans_mat.rows() - 1) / 2;
  Matrix<double> preconditioned =
      permutation(degree) * trans_mat * permutation((degree)).transpose();
  // Since the preconditioned transformation matrix has a different structure
  // for even and odd degrees, we have to distinguish these cases.
  if (degree % 2) {  // degree odd
    // s = 0 and m odd (l part)
    // size ((degree - 1)/2 + 1)^2
    // row index = 0, column index = degree + 1
    inv_upper_triangular_mat(0, degree + 1, (degree - 1) / 2 + 1,
                             preconditioned, inverse);

    // s = 0 and m even (l-1 part)
    // size ((degree - 1)/2 + 1)^2
    // row index = (degree - 1)/2 + 1, column index = (degree - 1)/2 + 1
    inv_upper_triangular_mat((degree - 1) / 2 + 1, (degree - 1) / 2 + 1,
                             (degree - 1) / 2 + 1, preconditioned, inverse);

    // s = 1 and m odd (l part)
    // size ((degree-1)/2 + 1)^2
    // row index = degree + 1, column index = 0
    inv_upper_triangular_mat(degree + 1, 0, (degree - 1) / 2 + 1,
                             preconditioned, inverse);

    // s = 1 and m even (l-1 part)
    // size ((degree - 1)/2)^2
    // row_index = 3 * (degree + 1)/2, column index = 3 * (degree + 1)/2
    inv_upper_triangular_mat(3 * (degree + 1) / 2, 3 * (degree + 1) / 2,
                             (degree - 1) / 2, preconditioned, inverse);
  } else {  // degree even
    // s = 0 and m even (l part)
    inv_upper_triangular_mat(0, 0, degree / 2 + 1, preconditioned, inverse);

    // s = 0 and m odd (l-1 part)
    inv_upper_triangular_mat(degree / 2 + 1, (3 * degree) / 2 + 1, degree / 2,
                             preconditioned, inverse);

    // s = 1 and m even (l part)
    // size triangular matrix (degree/2)^2 (no +1, since there is no Yl01)
    inv_upper_triangular_mat(degree + 1, degree + 1, degree / 2, preconditioned,
                             inverse);

    // s = 1 and m odd (l - 1 part)
    inv_upper_triangular_mat((3 * degree) / 2 + 1, degree / 2 + 1, degree / 2,
                             preconditioned, inverse);
  }
  return permutation(degree).transpose() * inverse * permutation(degree);
}

multipole_conv::Matrix<double> multipole_conv::dependent_components(
    Matrix<double>& multipole_basis_mat) {
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

#include "multipole_conv_funcs.h"

#include <cmath>
#include <cstddef>

#include "constants.h"
#include "math.h"
#include "square_matrix.h"

using std::size_t;

multipole_conv::SquareMatrix<double> multipole_conv::norms_real_sph(
    size_t degree) {
  SquareMatrix<double> res{2 * degree + 1};
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

multipole_conv::SquareMatrix<double> multipole_conv::norms_complex_sph(
    size_t degree) {
  SquareMatrix<double> res{2 * degree + 1};
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

multipole_conv::SquareMatrix<std::complex<double>>
multipole_conv::real_sph_to_complex(size_t degree) {
  using namespace std::complex_literals;
  SquareMatrix<std::complex<double>> transformation_matrix{2 * degree + 1};
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

multipole_conv::SquareMatrix<std::complex<double>>
multipole_conv::complex_sph_to_real(size_t degree) {
  // The transformation between the real and the complex spherical harmonics is
  // a unitary transformation, i.e its inverse is the transpose conjugate of the
  // corresponding transformation matrix
  SquareMatrix<std::complex<double>> transformation_matrix =
      real_sph_to_complex(degree).transpose().conjugate();
  return transformation_matrix;
}

double multipole_conv::coefficient(size_t degree, size_t order, bool p_index,
                                   size_t sum_start) {
  double sign = 1.;  // (-1)^(f(order))
  size_t sum_limit = 0;
  if (order % 2 == 0 && p_index) {  // m even, p = 1
    sum_limit = order / 2 - 1;
    sign = std::pow(-1., order / 2 + 1);
  } else if (order % 2 == 0 && !p_index) {  // m even, p = 0
    sum_limit = order / 2;
    sign = std::pow(-1., sum_limit);
  } else {  // m odd, p = 0, p = 1
    sum_limit = (order - 1) / 2;
    sign = std::pow(-1., (order + 1) / 2);
  }

  unsigned int p = static_cast<unsigned int>(p_index);
  double res = 0.;
  for (size_t i = sum_start; i <= sum_limit; ++i) {
    res += multipole_conv::binomial(order, 2 * i + p) *
           multipole_conv::binomial(i, sum_start);
  }
  double factor = 1 / factorial(degree - order);
  res *= sign * factor;
  return res;
}

multipole_conv::SquareMatrix<double> multipole_conv::basis_transformation(
    size_t degree) {
  // NOTE: The implementation of of this function strongly depends on the chosen
  // ordering of the real spherical harmonics and the multipole basis functions
  SquareMatrix<double> basis_transformation{2 * degree + 1};
  for (size_t order = degree, i = 0; i < degree; ++i, --order) {
    unsigned int p_upper_part =
        (order % 2 ? 1 : 0);  // p_index for the s = 0 case
    unsigned int p_lower_part =
        (order % 2 ? 0 : 1);  // p_index for the s = 1 case
    for (size_t k = 0; k <= std::floor(order / 2); ++k) {
      // Upper part (s = 0)
      basis_transformation(
          i, p_upper_part * (degree + 1) + order - (2 * k + p_upper_part)) =
          coefficient(degree, order, p_upper_part, k);
      // Lower part (s = 1)
      if (order % 2 == 0 && k == order / 2)
        continue;  // the sum for even m ends at k = order/2 - 1
      basis_transformation(2 * degree - i, p_lower_part * (degree + 1) + order -
                                               (2 * k + p_lower_part)) =
          coefficient(degree, order, p_lower_part, k);
    }
    // m = 0 and s = 0
    basis_transformation(degree, 0) = coefficient(degree, 0, 0, 0);
  }
  return basis_transformation;
}

multipole_conv::SquareMatrix<double> multipole_conv::permutation(
    size_t degree) {
  SquareMatrix<double> permutation_mat(2 * degree + 1);
  // l
  // upper part: r^l Y_ll0, r^l Y_ll-20 ...
  // lower part: r^l Y_ll1, r^l Y_ll-21 ...
  size_t num_permutations_l = std::floor(degree / 2);
  permutation_mat(0, 0) = 1.;                    // first row upper part (s = 0)
  permutation_mat(degree + 1, 2 * degree) = 1.;  // firsr row lower part (s = 1)
  for (size_t i = 1; i <= num_permutations_l; ++i) {
    permutation_mat(i, 2 * i) = 1.;                            // upper part
    // in the lower part (s = 1) Y_l01 does not exist
    if ( degree % 2 == 0 && i == num_permutations_l)
      continue;
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
    permutation_mat(num_permutations_l + 1 + i, 2 * i + 1) = 1.; // upper part
    // in the lower part (s = 1) Y_l-101 does not exist
    if ( degree % 2 == 1 && i == num_permutations_l_minus_one)
      continue;
    permutation_mat(degree + 1 + num_permutations_l_minus_one + 1 + i, 
                    2 * degree - 1 - 2 * i) = 1.; // lower part
  }
  return permutation_mat;
}

multipole_conv::SquareMatrix<double>
multipole_conv::invert_basis_transformation(const SquareMatrix<double>
					    &trans_mat) {
  SquareMatrix<double> inverse (trans_mat.dim());
  size_t degree = (trans_mat.dim() - 1)/2;
  SquareMatrix<double> preconditioned = permutation(degree) * trans_mat *
  permutation((degree)).transpose();
  // Since the preconditioned transformation matrix has a different structure
  // for even and odd degrees, we have to distinguish these cases. This leads to
  // duplication of code. TODO: improve when you have got time
  if (degree % 2) { 		// degree odd

  } else {			// degree even
    // s = 0 and m even
    // size triangular matrix (degree/2 + 1)/2 x (degree/2 + 1)/2
    // rhs b = unit vectors, e_1 ... e_degree/2
    for (size_t i = 0 ; i < degree/2 + 1; ++i) { // loop over the rhs'
      // b_k = delta_ik
      inverse(0, i) = (i == (degree/2)) ? 1/preconditioned(degree/2, 0) : 0;
      // int was necessary because I could not decrement size_t through zero
      for(int k = degree/2 - 1; k >= 0; --k) {
	double sum = 0.; // sum_j mat_kj x_j
	for(size_t j = 0; j < degree/2 - k; ++j)
	  sum += preconditioned(k, j) * inverse(j,i);
	// std::cout << k << "\n";
	// std::cout << preconditioned(degree/2 - k, degree/2 - k)  << "\n";
	// std::cout << preconditioned(2,2)  << "\n";

	 // b_k - sum_j mat_kj x_j
	inverse(degree/2 - k, i) = ((k == i ? 1 : 0) -
  sum)/preconditioned(k, degree/2 - k);
      }
    }
  }
  return inverse;
}

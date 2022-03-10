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

double multipole_conv::coefficient(std::size_t order, std::size_t k,
                                   bool p_index) {
  std::size_t limit = 0;
  if (order % 2 == 0 && p_index)  // m even, p = 1
    limit = order / 2 - 1;
  else if (order % 2 == 0 && !p_index)  // m even, p = 0
    limit = order / 2;
  else  // m odd, p = 0, p = 1
    limit = (order - 1) / 2;

  double res = 0.;
  unsigned int p = static_cast<unsigned int>(p_index);
  for (std::size_t i = k; i <= limit; ++i) {
    res += multipole_conv::binomial(order, 2 * i + p) *
           multipole_conv::binomial(i, k);
  }
  return res;
}

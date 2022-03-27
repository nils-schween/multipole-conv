#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "math.h"
#include "multipole_conv_funcs.h"
#include "options.h"
#include "matrix.h"
#include "square_matrix_operations.h"

int main(int argc, char *argv[]) {
  using namespace std::complex_literals;
  using namespace multipole_conv;

  std::size_t degree = 3;
  // Johnston
  Matrix<double> basis_transformation_mat = basis_transformation(degree);
  basis_transformation_mat.print();

  Matrix<double> multipole_basis =
      invert_basis_transformation(johnston_factor(degree) *
                                  condon_shortley_phase(degree) *
                                  basis_transformation_mat) /
      factorial(degree);
  multipole_basis.print();

  // Jackson
  Matrix<std::complex<double>> spherical_multipole_moments =
      (real_sph_to_complex(degree) * norms_real_sph(degree) *
       basis_transformation_mat)
          .conjugate();
  spherical_multipole_moments.print();

  return 0;
}

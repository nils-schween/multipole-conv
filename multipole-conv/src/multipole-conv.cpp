#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "math.h"
#include "matrix.h"
#include "multipole_conv_funcs.h"
#include "options.h"
#include "print_multipole_moments.h"
#include "square_matrix_operations.h"

using namespace multipole_conv;

template <typename T>
void pipeline(Matrix<T> &mat, MPOptions options, std::size_t degree) {
  // Include normalisation
  if ((options & MPOptions::normalisation) != MPOptions::none)
    mat = norms_sph(degree, options) * mat;

  // Remove Condon Shortley Phase
  if ((options & MPOptions::remove_condon_shortley_phase) != MPOptions::none)
    mat = condon_shortley_phase(degree) * mat;

  // Include (and split) the addition theorem factor
  if ((options & MPOptions::include_addition_theorem) != MPOptions::none)
    mat = addition_theorem_factor(degree, options) * mat;

  if ((options & MPOptions::cartesian) == MPOptions::none) {
    // End of pipeline for the spherical multipole moments
    print_spherical_multipole_moments(mat, options);
  } else {  // Cartesian multipole moments
    // Express the multipole basis functions as a linear combination of the
    // (real/complex) spherical harmonics
    mat = invert_basis_transformation(mat);
    Matrix<T> dep_comps(degree);
    // Compute the rest of the components of Cartesian multipole moments
    dep_comps = dependent_components(mat);

    // End of pipeline Cartesian multipole moments
    print_cartesian_multipole_moments(mat, dep_comps, options);
  }
}

int main(int argc, char *argv[]) {
  MPOptions options = MPOptions::complex | MPOptions::normalisation;
  std::size_t degree = 3;
  if ((options & MPOptions::complex) != MPOptions::none) {
    Matrix<std::complex<double>> mat =
        convert_real_mat_to_complex_mat(basis_transformation(degree));
    pipeline(mat, options, degree);
  } else {
    Matrix<double> mat = basis_transformation(degree);
    pipeline(mat, options, degree);
  }

  return 0;
}

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "math.h"
#include "matrix.h"
#include "matrix_operations.h"
#include "multipole_conv_funcs.h"
#include "options.h"
#include "print_multipole_moments.h"

using namespace multipole_conv;

template <typename T>
void pipeline_spherical_multipole_moments(Matrix<T> &mat, MPOptions options,
                                          std::size_t degree) {
  if ((options & MPOptions::normalisation) != MPOptions::none)
    mat = norms_sph(degree, options) * mat;

  if ((options & MPOptions::remove_condon_shortley_phase) != MPOptions::none)
    mat = condon_shortley_phase(degree) * mat;

  if ((options & MPOptions::include_addition_theorem) != MPOptions::none)
    mat = addition_theorem_factor(degree, options) * mat;
}

template <typename T>
void pipeline_cartesian_multipole_moments(Matrix<T> &inv_mat, MPOptions options,
                                          std::size_t degree) {
  // CMP = Cartesian multipole moment; SMP = spherical multipole moment
  // CMP = INV_BASIS_TRANSFORMATION * U^dagger * N^(-1) (N*U*SMP)
  if ((options & MPOptions::normalisation) != MPOptions::none)
    inv_mat = inv_mat * inverse_diagonal(norms_sph(degree, options));

  // CONDON_SHORTLEY_PHASE^(-1) = CONDON_SHORTLEY_PHASE
  if ((options & MPOptions::remove_condon_shortley_phase) != MPOptions::none)
    inv_mat = inv_mat * condon_shortley_phase(degree);

  if ((options & MPOptions::include_addition_theorem) != MPOptions::none)
    inv_mat =
        inv_mat * inverse_diagonal(addition_theorem_factor(degree, options));

  if ((options & MPOptions::include_l_factorial) != MPOptions::none)
    inv_mat = inv_mat / factorial(degree);
}

int main(int argc, char *argv[]) {
  MPOptions options = MPOptions::normalisation;
  // Johnston
  // MPOptions options =
  //     MPOptions::cartesian | MPOptions::include_addition_theorem |
  //     MPOptions::remove_condon_shortley_phase |
  //     MPOptions::include_l_factorial | MPOptions::dependent_components;
  std::cout << options << "\n";
  std::size_t degree = 30;

  if ((options & MPOptions::cartesian) != MPOptions::none) {  // CMPs
    Matrix<double> basis_trans = basis_transformation(degree);
    Matrix<double> inv_trans = invert_basis_transformation(basis_trans);
    basis_trans = 0;  // free memory
    if ((options & MPOptions::complex) != MPOptions::none) { // complex
      if ((options & MPOptions::dependent_components) != MPOptions::none) {
        Matrix<std::complex<double>> dep_comps =
            convert_real_to_complex(dependent_components(inv_trans));
        // Convert to complex solid harmonics
        dep_comps =
            dep_comps * real_sph_to_complex(degree).transpose().conjugate();
        pipeline_cartesian_multipole_moments(dep_comps, options, degree);
        print_dependent_components(dep_comps, options);
      }
      Matrix<std::complex<double>> cplx_inv_trans =
          convert_real_to_complex(inv_trans);
      // Transform the real solid harmonics to complex solid harmonics
      // CMP = INV_BASIS_TRANSFORMATION * U^dagger (U*SMP)
      cplx_inv_trans =
          cplx_inv_trans * real_sph_to_complex(degree).transpose().conjugate();
      inv_trans = 0;  // free memory
      pipeline_cartesian_multipole_moments(cplx_inv_trans, options, degree);
      print_cartesian_multipole_moments(cplx_inv_trans, options);
    } else { 			// real
      if ((options & MPOptions::dependent_components) != MPOptions::none) {
        Matrix<double> dep_comps = dependent_components(inv_trans);
        pipeline_cartesian_multipole_moments(dep_comps, options, degree);
        print_dependent_components(dep_comps, options);
      }
      pipeline_cartesian_multipole_moments(inv_trans, options, degree);
      print_cartesian_multipole_moments(inv_trans, options);
    }
  } else {  // SMPs
    if ((options & MPOptions::complex) != MPOptions::none) {
      Matrix<std::complex<double>> basis_trans =
          convert_real_to_complex(basis_transformation(degree));
      // Transform the real solid harmonics to complex solid harmonics
      basis_trans = real_sph_to_complex(degree) * basis_trans;
      pipeline_spherical_multipole_moments(basis_trans, options, degree);
      print_spherical_multipole_moments(basis_trans, options);
    } else {
      Matrix<double> basis_trans = basis_transformation(degree);
      pipeline_spherical_multipole_moments(basis_trans, options, degree);
      print_spherical_multipole_moments(basis_trans, options);
    }
  }

  return 0;
}

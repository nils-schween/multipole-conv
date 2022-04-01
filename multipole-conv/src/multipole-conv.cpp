#include <cmath>
#include <complex>
#include <iostream>
#include <utility>
#include <vector>

#include "command_line_parser.h"
#include "conventions.h"
#include "math.h"
#include "matrix.h"
#include "matrix_operations.h"
#include "multipole_conv_funcs.h"
#include "options.h"
#include "print_multipole_moments.h"

using namespace multipole_conv;

int main(int argc, char *argv[]) {
  // CMP = Cartesian multipole moment; SMP = spherical multipole moment
  std::pair<std::size_t, MPOptions> cmd_line_arguments = cmd_parser(argc, argv);
  MPOptions options = cmd_line_arguments.second;
  std::size_t degree = cmd_line_arguments.first;
  if (options == MPOptions::none)
    std::cout << "No options were set. Default: Real solid harmonics (without "
                 "normalisation) are computed.\n\n";
  else
    std::cout << options << "\n";

  if ((options & MPOptions::cartesian) != MPOptions::none) {  // CMPs
    if ((options & MPOptions::complex) != MPOptions::none) {
      Matrix<std::complex<double>> cmp =
          pipeline_cmp_cmplx_smp(options, degree);
      print_cartesian_multipole_moments(cmp, options);
      if ((options & MPOptions::dependent_components) != MPOptions::none) {
        Matrix<std::complex<double>> dep_comps = dependent_components(cmp);
        print_dependent_components(dep_comps, options);
      }
    } else {
      Matrix<double> cmp = pipeline_cmp_real_smp(options, degree);
      print_cartesian_multipole_moments(cmp, options);
      if ((options & MPOptions::dependent_components) != MPOptions::none) {
        Matrix<double> dep_comps = dependent_components(cmp);
        print_dependent_components(dep_comps, options);
      }
    }
  } else {  // SMPsp
    if ((options & MPOptions::complex) != MPOptions::none) {
      Matrix<std::complex<double>> smp = pipeline_cmplx_spm(options, degree);
      print_spherical_multipole_moments(smp, options);
    } else {
      Matrix<double> smp = pipeline_real_smp(options, degree);
      print_spherical_multipole_moments(smp, options);
    }
  }
  return 0;
}

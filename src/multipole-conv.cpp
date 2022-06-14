/*
 * Copyright (C) 2022 Schween, Nils W. and Reville, Brian
 *
 * This file is part of multipole-conv.
 *
 * multipole-conv is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * multipole-conv is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with multipole-conv. If not, see <https://www.gnu.org/licenses/>.
 */

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
  if (options == MPOptions::none) {
    std::cout << R"(Type "multipole-conv -h" to see the available options. )"
              << "\n";
    return 0;
  } else {
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
}

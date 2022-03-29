#ifndef PRINT_MULTIPOLE_MOMENTS_H
#define PRINT_MULTIPOLE_MOMENTS_H

#include <cstddef>
#include <iomanip>

#include "matrix.h"
#include "options.h"

template <typename T>
void print_spherical_multipole_moments(multipole_conv::Matrix<T> mat,
                                       multipole_conv::MPOptions options) {
  std::size_t degree = (mat.rows() - 1) / 2;
  if ((options & multipole_conv::MPOptions::complex) !=
      multipole_conv::MPOptions::none) {
    for (int order = degree, i = 0; i < mat.rows(); ++i, --order) {
      std::cout << "rho^" << std::setw(3) << order << "_" << degree << " = ";
      for (std::size_t j = 0; j < mat.columns(); ++j) {
      }
      std::cout << "\n";
    }
  } else {
    for (int order = degree, i = 0; i < mat.rows(); ++i, --order) {
      std::size_t s = (i > degree) ? 1 : 0;
      std::size_t m = ((i > degree) ? -1 : 1) * order;
      std::cout << "rho_" << std::setw(2) << degree << "_" << std::setw(2) << m
                << "_" << s << " = ";
      for (std::size_t j = 0; j < mat.columns(); ++j) {
      }
      std::cout << "\n";
    }
  }
}

template <typename T>
void print_cartesian_multipole_moments(multipole_conv::Matrix<T> mat,
                                       multipole_conv::MPOptions options) {
  mat.print();
}

template <typename T>
void print_dependent_components(multipole_conv::Matrix<T> mat,
                                multipole_conv::MPOptions options) {
  mat.print();
}

#endif

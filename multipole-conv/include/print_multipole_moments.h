#ifndef PRINT_MULTIPOLE_MOMENTS_H
#define PRINT_MULTIPOLE_MOMENTS_H

#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>

#include "matrix.h"
#include "options.h"

namespace multipole_conv {
template <typename T>
void print_spherical_multipole_moments(const Matrix<T>& mat,
                                       const MPOptions& options) {
  std::size_t degree = (mat.rows() - 1) / 2;
  if ((options & MPOptions::complex) !=
      MPOptions::none) {
    std::cout << std::showpos;
    for (int order = degree, i = 0; i < mat.rows(); ++i, --order) {
      // std::size_t q_index = ((i > degree) ? -1 : 1) * order;
      std::cout << "rho^" << std::setw(3) << std::left << order << "_" << degree
                << " =";
      for (std::size_t j = 0; j < mat.columns(); ++j) {
        // it is not possible to call mat(i,j).real/imag
        // I have to create a copy and then I can call it. I do not know why.
        std::complex<double> current_element{mat(i, j)};
        std::size_t p_index = (j > degree ? 1 : 0);
        std::size_t q_index = (j > degree ? j - degree : degree - j) - p_index;
        std::size_t r_index = degree - p_index - q_index;
        std::cout << " " << current_element.real() << current_element.imag()
                  << "i"
                  << " M_" << p_index << "," << q_index << "," << r_index;
      }
      std::cout << "\n";
    }
  } else {
    std::cout << std::showpoint;
    for (int order = degree, i = 0; i < mat.rows(); ++i, --order) {
      std::size_t s = (i > degree) ? 1 : 0;
      // no negative order for real solid harmonics (i.e. m >= 0)
      std::size_t m = ((i > degree) ? -1 : 1) * order;
      std::cout << "rho_" << degree << "," << m << "," << s << " = ";
      for (std::size_t j = 0; j < mat.columns(); ++j) {
        std::size_t p_index = (j > degree ? 1 : 0);
        std::size_t q_index = (j > degree ? j - degree : degree - j) - p_index;
        std::size_t r_index = degree - p_index - q_index;
        std::cout << " " << mat(i, j) << " M_" << p_index << "," << q_index
                  << "," << r_index;
      }
      std::cout << "\n";
    }
  }
}

template <typename T>
void print_cartesian_multipole_moments(const Matrix<T>& mat,
                                       const MPOptions& options) {
  std::size_t degree = (mat.rows() - 1) / 2;

  for (std::size_t i = 0; i < mat.rows(); ++i ) {
    std::size_t p_index = (i > degree ? 1 : 0);
    std::size_t q_index = (i > degree ? i - degree : degree - i) - p_index;
    std::size_t r_index = degree - p_index - q_index;
    std::cout << std::showpos;
    std::cout << "M_" << p_index << "," << q_index << "," << r_index << " =";
    if((options &  MPOptions::complex) != MPOptions::none) {
      for(int order = degree, j = 0; j < mat.columns(); ++j, --order) {
	std::complex<double> current_element {mat(i,j)};
	std::cout << " " << current_element.real() << current_element.imag()
		  << "i"
		  << " rho^" << order << "_" << degree;
      }
    }
    else {
      for(int order = degree, j = 0; j < mat.columns(); ++j, --order) {
	std::size_t m = ((j > degree) ? -1 : 1) * order;
	std::size_t s = (j > degree) ? 1 : 0;
      }
    }
    std::cout << "\n";
  }
}




template <typename T>
void print_dependent_components(const Matrix<T>& mat,
                                const MPOptions& options) {
  mat.print();
}
}
#endif

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
  std::cout << "Spherical multipole moments\n\n";
  if ((options & MPOptions::complex) != MPOptions::none) {
    std::cout << std::showpos;
    for (int order = degree, i = 0; i < mat.rows(); ++i, --order) {
      // std::size_t q_index = ((i > degree) ? -1 : 1) * order;
      std::cout << "q^" << std::setw(3) << std::left << order << "_" << degree
                << " =";
      for (std::size_t j = 0; j < mat.columns(); ++j) {
        // it is not possible to call mat(i,j).real/imag
        // I have to create a copy and then I can call it. I do not know why.
        std::complex<double> current_element{mat(i, j)};
        std::size_t p_index = (j > degree ? 1 : 0);
        std::size_t q_index = (j > degree ? j - degree : degree - j) - p_index;
        std::size_t r_index = degree - p_index - q_index;
	// only print elements unequal zero
	if (current_element.real() != 0. || current_element.imag() != 0.) {
        std::cout << " " << current_element.real() << current_element.imag()
                  << "i"
                  << " Q_" << p_index << "," << q_index << "," << r_index;
	}
      }
      std::cout << "\n";
    }
  } else {
    std::cout << std::showpoint;
    for (int order = degree, i = 0; i < mat.rows(); ++i, --order) {
      std::size_t s = (i > degree) ? 1 : 0;
      // no negative order for real solid harmonics (i.e. m >= 0)
      std::size_t m = ((i > degree) ? -1 : 1) * order;
      std::cout << "q_" << degree << "," << m << "," << s << " = ";
      for (std::size_t j = 0; j < mat.columns(); ++j) {
        std::size_t p_index = (j > degree ? 1 : 0);
        std::size_t q_index = (j > degree ? j - degree : degree - j) - p_index;
        std::size_t r_index = degree - p_index - q_index;
	if( mat(i,j) != 0.) {
        std::cout << " " << mat(i, j) << " Q_" << p_index << "," << q_index
                  << "," << r_index;
	}
      }
      std::cout << "\n";
    }
  }
  std::cout << "\n";
}

template <typename T>
void print_cartesian_multipole_moments(const Matrix<T>& mat,
                                       const MPOptions& options) {
  std::size_t degree = (mat.rows() - 1) / 2;
  std::cout << "Cartesian multipole moments (independent components = "
               "multipole basis functions)\n\n";
  for (std::size_t i = 0; i < mat.rows(); ++i) {
    std::size_t p_index = (i > degree ? 1 : 0);
    std::size_t q_index = (i > degree ? i - degree : degree - i) - p_index;
    std::size_t r_index = degree - p_index - q_index;
    std::cout << "Q_" << p_index << "," << q_index << "," << r_index << " =";
    if ((options & MPOptions::complex) != MPOptions::none) {
      std::cout << std::showpos;
      for (int order = degree, j = 0; j < mat.columns(); ++j, --order) {
        std::complex<double> current_element{mat(i, j)};
	// only print non zero elements
	if (current_element.real() != 0. || current_element.imag() != 0.) {
        std::cout << " " << current_element.real() << current_element.imag()
                  << "i"
                  << " q^" << order << "_" << degree;
	}
      }
    } else {
      for (int order = degree, j = 0; j < mat.columns(); ++j, --order) {
        std::size_t m = ((j > degree) ? -1 : 1) * order;
        std::size_t s = (j > degree) ? 1 : 0;
	if( mat(i,j) != 0.) {
        std::cout << " " << mat(i, j) << " q_" << degree << "," << m << ","
                  << s;
	}
      }
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

template <typename T>
void print_dependent_components(const Matrix<T>& mat,
                                const MPOptions& options) {
  std::cout << "Cartesian multipole moments (dependent components)\n\n";
  std::size_t degree = (mat.columns() - 1) / 2;
  // Because of the ordering of the dependent components the for loops are more
  // complicated
  for (std::size_t order = degree, r_index = 0, row_idx = 0; order > 1;
       --order, ++r_index) {
    for (std::size_t p_index = order; p_index > 1; --p_index, ++row_idx) {
      std::size_t q_index = degree - p_index - r_index;
      std::cout << "Q_" << p_index << "," << q_index << "," << r_index << " =";
      if ((options & MPOptions::complex) != MPOptions::none) {
        std::cout << std::showpos;
        for (int order_spm = degree, column_idx = 0; column_idx < mat.columns();
             ++column_idx, --order_spm) {
          std::complex<double> current_element{mat(row_idx, column_idx)};
	  // only print non-zero elements
	  if (current_element.real() != 0. || current_element.imag() != 0.) {
          std::cout << " " << current_element.real() << current_element.imag()
                    << "i"
                    << " q^" << order_spm << "_" << degree;
	  }
        }
      } else {
        for (int order_spm = degree, column_idx = 0; column_idx < mat.columns();
             ++column_idx, --order_spm) {
          std::size_t m = ((column_idx > degree) ? -1 : 1) * order_spm;
          std::size_t s = (column_idx > degree) ? 1 : 0;
	  if ( mat(row_idx, column_idx) != 0.) {
          std::cout << " " << mat(row_idx, column_idx) << " q_" << degree
                    << "," << m << "," << s;
	  }
        }
      }
      std::cout << "\n";
    }
  }
  std::cout << "\n";
}

}  // namespace multipole_conv
#endif

#ifndef PRINT_MULTIPOLE_MOMENTS_H
#define PRINT_MULTIPOLE_MOMENTS_H

#include "matrix.h"
#include "options.h"

template<typename T>
void print_spherical_multipole_moments(multipole_conv::Matrix<T> mat,
				       multipole_conv::MPOptions options) {
  mat.print();
}

template<typename T>
void print_cartesian_multipole_moments(multipole_conv::Matrix<T> mat,
				       multipole_conv::Matrix<T> dep_comps,
				       multipole_conv::MPOptions options) {
  mat.print();
}

#endif

#include <cmath>
#include <complex>
#include <iostream>

#include "options.h"
#include "square_matrix.h"
#include "square_matrix_operations.h"
#include "multipole_conv_funcs.h"

int main(int argc, char *argv[]) {
  using namespace std::complex_literals;
  using namespace multipole_conv;

  SquareMatrix<double> basis_transformation_mat = basis_transformation(4);
  basis_transformation_mat.print();

  return 0;
}

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

  std::size_t degree = 6;
  SquareMatrix<double> basis_transformation_mat = basis_transformation(degree);
  // basis_transformation_mat.print();

  // SquareMatrix<double> permutation_mat = permutation(5);
  // permutation_mat.print();
  SquareMatrix<double> temp = permutation(degree)*basis_transformation_mat*permutation(degree).transpose();
  temp.print();

  SquareMatrix<double> inverse = invert_basis_transformation(basis_transformation_mat);
  (inverse * temp).print();
  return 0;
}

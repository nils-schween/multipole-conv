#include <cmath>
#include <complex>
#include <iostream>

#include "math.h"
#include "options.h"
#include "square_matrix.h"
#include "square_matrix_operations.h"
#include "multipole_conv_funcs.h"

int main(int argc, char *argv[]) {
  using namespace std::complex_literals;
  using namespace multipole_conv;

  MPOptions test = MPOptions::complex | MPOptions::normalisation;

  std::complex<double> a = 3. + 1.i;
  // std::cout << is_complex<int>() ;
  SquareMatrix<double> normalisations {norms_real_sph(3)};
  // transform_real_sph_to_complex(3);
  normalisations.print();
  // inverse_diagonal(normalisations).print();
  SquareMatrix<double> A{3};
  A(1, 1) = 1;
  A(2, 2) = 2;
  A(0, 2) = 3;
  A(1, 0) = 1;
  
  std::cout << "A \n";
  A.print();
  std::cout << "\n";
  auto B{A * 2.};
  std::cout << "B \n";
  B.print();
  B.transpose().print();

  std::cout << "\n";
  std::cout << "A * B \n";
  ((3.) * A).print();
  auto C = A * B;
  C.print();

  SquareMatrix<std::complex<double>> D(3);
  D(1, 1) = 0.3 + 2i;
  (D / (2. + 1i)).print();
  std::cout << factorial(4) << "\n";
  std::cout << factorial(3) << "\n";
  std::cout << factorial(2) << "\n";
  std::cout << factorial(1) << "\n";
  std::cout << factorial(0) << "\n";
    
  (D * A).print();
  std::cout << "\n";
  D.conjugate().print();
  std::cout << test;
  std::cout << binomial(5, 3) << "\n";

  SquareMatrix<std::complex<double>> transformation =
  transform_complex_sph_to_real(2);
  transformation.print();
  return 0;
}

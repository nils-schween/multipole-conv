#include <iostream>
#include <complex>
#include <type_traits>
#include "math.h"
#include "options.h"
#include "square_matrix.h"

int main(int argc, char *argv[])
{
  using namespace std::complex_literals;
  MPOptions test = MPOptions::complex | MPOptions::normalisation;

  std::complex<float> a = 3;
  // std::cout << is_complex<int>() ;
    
  SquareMatrix<double> A(3);
  A(1,1) = 1;
  A(2,2) = 2;
  A(0,2) = 3;
  A(1,0) = 1;
  std::cout << "A \n";
  A.print();
  std::cout << "\n";
  auto B {(2i + 0.3)*A};
  std::cout << "B \n";
  B.print();
  (B/2.).transpose().print();
  
  std::cout << "\n";
  std::cout << "A * B \n";

  auto C = A*B;
  C.print();

  SquareMatrix<std::complex<double>> D(3);
  D(1,1) = 0.3 + 2i;
  std::cout << "\n";
  D.conjugate().print();
  std::cout << test;
  std::cout << binomial(3, 3) << "\n";
  return 0;
}

#include <iostream>
#include <complex>
#include <type_traits>
#include "math.h"
#include "options.h"
#include "square_matrix.h"

// Type traits
// template<typename T>
// struct is_complex : std::false_type {};

// template<typename T>
// struct is_complex<std::complex<T>> : public std::true_type {};

int main(int argc, char *argv[])
{
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
  SquareMatrix<double> B {A/2};
  std::cout << "B \n";
  B.print();
  B.transpose().print();
  
  std::cout << "\n";
  std::cout << "A * B \n";

  SquareMatrix<double> C(3);
  
  C = A*B;
  C.print();
  
  std::cout << test;
  std::cout << binomial(3, 3) << "\n";
  return 0;
}

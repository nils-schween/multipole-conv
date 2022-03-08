#include <iostream>
#include <complex>
#include <cmath>

#include "math.h"
#include "options.h"
#include "square_matrix.h"

constexpr double pi = 3.14159265358979;

template <typename T>
class Normalisation : public SquareMatrix<T> {
   public:
  Normalisation(std::size_t l, MPOptions opt) : SquareMatrix<T> {2 * l + 1}, options{opt}
    {
      double factor {0.};
      if ((MPOptions::complex & opt) != MPOptions::none) {
	factor = {std::sqrt((2 * l + 1)/(4*pi))};
	
      }
      else {
	factor = {std::sqrt((2 * l + 1)/(2*pi))};
	for(std::size_t m = 1; m <= l; ++m) {
	  double norm_l_m = factor * std::sqrt((std::tgamma(l - m + 1)/std::tgamma(l + m + 1)));
	  *this->operator()(m, m) = norm_l_m;
	  *this->operator()(2*l - m, 2*l - m) = norm_l_m;
	}
	*this->operator()(0,0) = factor;
      }
    }
  void real_or_complex() {
    if ((MPOptions::complex & options) != MPOptions::none)
      std::cout << "Normalisation of complex spherical harmonics" << "\n";
    else
      std::cout << "Normalisation of real spherical harmonics" << "\n";
  }
  
  Normalisation inverse();
   private:
  MPOptions options;
};

template <typename T>
Normalisation<T> Normalisation<T>::inverse() {
  Normalisation<T> res(*this->dim());
  for(std::size_t i = 0; i < *this->dim(); ++i)
    res(i, i) = 1./(*this->operator()(i,i));
  return res;
}

int main(int argc, char *argv[])
{
  using namespace std::complex_literals;
  MPOptions test = MPOptions::complex | MPOptions::normalisation;

  std::complex<double> a = 3. + 1.i;
  // std::cout << is_complex<int>() ;
    
  SquareMatrix<double> A {3};
  A(1,1) = 1;
  A(2,2) = 2;
  A(0,2) = 3;
  A(1,0) = 1;
  std::cout << "A \n";
  A.print();
  std::cout << "\n";
  auto B {A*2.};
  std::cout << "B \n";
  B.print();
  B.transpose().print();
  
  std::cout << "\n";
  std::cout << "A * B \n";
  ((3.)*A).print();
  auto C = A*B;
  C.print();

  SquareMatrix<std::complex<double>> D(3);
  D(1,1) = 0.3 + 2i;
  (D/(2. + 1i)).print();

  (D*A).print();
  std::cout << "\n";
  D.conjugate().print();
  std::cout << test;
  std::cout << binomial(5, 3) << "\n";
  return 0;
}

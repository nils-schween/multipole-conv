#include "math.h"

// Computes the binomial coefficient
double multipole_conv::binomial(unsigned int n, unsigned int k) {
  double r = 1.;
  for (unsigned int d = 1.; d <= k; ++d) {
    r *= n--;
    r /= d;
  }
  return r;
}

// Computes n factorial
double multipole_conv::factorial(unsigned int n) {
  double res = (n > 1) ? n : 1;
  if (n > 2)
    for (unsigned int i = n - 1; i > 1; --i) res *= i;
  return res;
}

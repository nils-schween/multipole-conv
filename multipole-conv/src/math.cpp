#include "math.h"

unsigned int binomial(unsigned int n, unsigned int k) {
  unsigned r = 1.;
  for (unsigned int d = 1; d <= k; ++d) {
    r *= n--;
    r /= d;
  }
  return r;
}

#include <iostream>
#include "math.h"
#include "options.h"

int main(int argc, char *argv[])
{
  MPOptions test = MPOptions::complex | MPOptions::normalisation;
  std::cout << test;
  std::cout << binomial(3, 3) << "\n";
  return 0;
}

/*
 * Copyright (C) 2022 Schween, Nils W. and Reville, Brian
 *
 * This file is part of multipole-conv.
 *
 * multipole-conv is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * multipole-conv is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with multipole-conv. If not, see <https://www.gnu.org/licenses/>.
 */

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

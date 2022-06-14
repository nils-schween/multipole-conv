# Copyright (C) 2022 Schween, Nils W. and Reville, Brian
#
# This file is part of multipole-conv.
#
# multipole-conv is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# multipole-conv is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with multipole-conv. If not, see <https://www.gnu.org/licenses/>.


import math
import numpy as np

def norm_sh(degree, order):
    """Computes the normalisation  Laplace's spherical harmonic
    of degree and order.
    """
    norm = np.sqrt(
        (2 * degree + 1) / (4 * np.pi) * (math.factorial(degree - order)) /
        (math.factorial(degree + order)))
    return norm


def norm_real_sh(degree, order):
    """Computes the normalisation of real Laplace's spherical harmonic
    corresponding to degree and order.
    """
    norm = np.sqrt(
        (2 * degree + 1) / (2 * np.pi) * (math.factorial(degree - order)) /
        (math.factorial(degree + order)))
    return norm if order > 0 else norm * 1/np.sqrt(2)

# N
def normalisations(degree):
    """Computes the normalisaton of all real Laplace's spherical harmonics with
    degree and returns them in a NumPy array, namely [N^degree_degree, ...,
    N^(-degree)_degree].
    """
    norms = np.zeros((2 * degree + 1, 2*degree+1))
    for i in range(degree):
        order = degree - i
        norm = norm_real_sh(degree, order)
        norms[i, i] = norm
        norms[2 * degree - i, 2*degree - i] = norm
    norms[degree, degree] = norm_real_sh(degree, 0)
    return norms

# S
def complex_sh_to_real_sh(degree):
    """ Creates the transformation matrix to convert from complex to real
    Laplace's spherical harmonics.
    """
    transformation_matrix = np.zeros((2 * degree + 1, 2 * degree + 1), complex)
    for i in range(degree):
        order = degree - i
        factor = np.sqrt(1 / 2)
        # Upper part of the transformation matrix
        transformation_matrix[i, i] = factor
        transformation_matrix[i, 2 * degree - i] = factor * (-1)**order
        # Lower part of the transformation matrix
        transformation_matrix[2 * degree - i, i] = factor * (0 + -1j)
        transformation_matrix[2 * degree - i, 2 * degree - i] = \
            factor * (-1)**order * (0 + 1j)

    transformation_matrix[degree, degree] = 1
    return transformation_matrix


def sum_in_coefficients_with_p_zero(order, k):
    """A helper function to compute the basis transformation matrix
    """
    temp = 0
    for n in range(k, int(order / 2) + 1):
        temp += math.comb(order, 2 * n) * math.comb(n, k)
    return temp


def sum_in_coefficients_with_p_one(order, k):
    """A helper function to compute the basis transformation matrix
    """
    temp = 0
    for n in range(k, int(order / 2) + order % 2):
        temp += math.comb(order, 2 * n + 1) * math.comb(n, k)
    return temp


def factor_in_coefficients(degree, order):
    """ Compute the factor in the coefficients
    """
    temp = np.sqrt((2 * degree + 1)/(4 * np.pi)
                   * 1/(math.factorial(degree - order)
                        * math.factorial(degree + order)))
    if order % 2 == 0:
        a_m_l = temp
        i_and_minus_one = (-1)**(abs(order / 2))
    else:
        i_and_minus_one = (-1)**(abs(int(order / 2))) * (0 + 1j)
        a_m_l = -temp
    return i_and_minus_one * a_m_l

# D_1 (or D_3)
def addition_theorem(degree, normalised=False, real=True):
    """ Creates a diagonal matrix with the addition theorem factor. The addition
    theorem factor is different if the (real) spherical harmonics are not
    normalised. """

    # factorial = (-1)**order * (0 + 1j)**order * (math.factorial(degree)) \
        #     * (math.factorial(degree - order))**2
    #  * math.prod(range(2 * degree - 1, 0, -2)
    addition_theorem_factor = np.zeros((2*degree + 1, 2*degree + 1))
    if (normalised):
        factor = (4*np.pi)/(2*degree + 1)
        np.fill_diagonal(addition_theorem_factor, factor)
        return addition_theorem_factor
    else:
        if (real):
            for i in range(degree):
                order = degree - i
                factor = 2 * math.factorial(degree - order) \
                    / math.factorial(degree + order)
                addition_theorem_factor[i, i] = factor
                addition_theorem_factor[2*degree - i, 2*degree - i] = factor

            addition_theorem_factor[degree, degree] = 1
            return addition_theorem_factor
        else:
            for i in range(2*degree + 1):
                order = degree - i
                factor =  math.factorial(degree - order) / math.factorial(degree + order)
                addition_theorem_factor[i, i] = factor
            return addition_theorem_factor


# D_2
def l_factorial(degree):
    """ Creates a diagonal with a l factorial.
    """
    return math.factorial(degree) * np.eye(2*degree + 1)

# C
def condon_shortley_phase(degree):
    """ Creates a diagonal matrix with the Condon-Shortley phase.
    """
    csp = np.zeros((2*degree + 1, 2*degree + 1))
    for i in range(degree):
        order = degree - i
        sign = 1 if order % 2 == 0 else -1
        csp[i, i] = sign
        csp[2*degree - i, 2*degree - i] = sign

    csp[degree, degree] = 1
    return csp

# P
def permutation_matrix(degree):
    """ Create a permuation matrix, which is used to transform the basis
    transformation matrix in a way that it consists of four upper triangular
    matrices.
    """
    permutation_mat = np.zeros((2*degree + 1, 2*degree + 1))
    # Notation: Y_lm
    # l
    # upper part: r^l Y_ll, r^l Y_ll-2 ...
    # lower part: r^l Y_l(-l), r^l Y_l(-(l-2)) ...
    num_permutations_l = math.floor(degree / 2)

    permutation_mat[0, 0] = 1.    # first row upper part (m >= 0)
    permutation_mat[degree + 1, 2 * degree] = 1.  # first row lower part ( m < 0)

    for i in range(num_permutations_l + 1):
            permutation_mat[i, 2 * i] = 1.  # upper part
            # in the lower part ( m < 0) Y_l0 does not exist
            if degree % 2 == 0 and i == num_permutations_l:
                continue
            permutation_mat[degree + 1 + i, 2 * degree - 2 * i] = 1. # lower part

    # l-1
    # upper part: r^l Y_l-1l-1, r^l Y_l-1l-3 ...
    # lower part: r^l Y_l-1(-(l-1)), r^l Y_l-1(-(l-3)) ...
    num_permutations_l_minus_one = math.floor((degree - 1) / 2)

    # first row upper part
    permutation_mat[num_permutations_l + 1, 1] = 1.
    # first row lower part
    permutation_mat[degree + 1 + num_permutations_l_minus_one + 1,
                    2 * degree - 1] = 1.
    for i in range(num_permutations_l_minus_one + 1):
        permutation_mat[num_permutations_l + 1 + i, 2 * i + 1] = 1.  # upper part
        # in the lower part (m < 0) Y_l-10 does not exist
        if degree % 2 == 1 and i == num_permutations_l_minus_one:
            continue
        permutation_mat[degree + 1 + num_permutations_l_minus_one + 1 + i,
                        2 * degree - 1 - 2 * i] = 1.  # lower part

    return permutation_mat


# B
def beta_matrix(degree):
    """ Computes the basis transformation between Laplace's spherical harmonics
    and the basis obtained by an application of Efimov's ladder operator.
    """
    basis_transformation = np.zeros((2 * degree + 1, 2 * degree + 1), complex)
    for i in range(degree + 1):
        order = degree - i
        factor = factor_in_coefficients(degree, order)
        for k in range(int(order / 2) + 1):
            # p = 0
            sum_p_zero = sum_in_coefficients_with_p_zero(order, k)
            # sum_p_zero = 1
            # Upper part (m >= 0)
            basis_transformation[i, i + 2 * k] = factor * sum_p_zero
            # Lower part (m < 0)
            # Do not overwrite what was written for the m = 0 case
            if i != degree:
                basis_transformation[2 * degree - i,
                                      i + 2*k] = (-1)**order \
                                         * factor.conjugate() * sum_p_zero

            # p = 1
            # If the order is even, then the k = int(order/2) then p = 1 basis
            # function is too much, and we must not store an entry in the
            # matrix
            if order % 2 == 0 and k == int(order / 2):
                continue
            sum_p_one = sum_in_coefficients_with_p_one(order, k)
            # sum_p_one = 1
            # Upper part (m >= 0)
            basis_transformation[i,
                                 degree + order - 2*k] = \
                                     (0 + -1j) * factor * sum_p_one
            # Lower part (m < 0)
            basis_transformation[2 * degree - i,
                                 degree + order - 2*k] = (-1)**order \
                                     * (0 + 1j) * factor.conjugate() * sum_p_one

    return basis_transformation


def helper_function(degree, order, j):
    """ Helper function
    """
    result = 1/(2**degree) * (-1)**(order + j) * math.comb(degree, j) \
        * math.comb(2*degree - 2*j, degree)  \
        * math.factorial(degree - 2*j)/math.factorial(degree - 2*j - order)
    return result

def sum_in_alphas(degree, p, q, k):
    """ A helper function to compute the sum in the alpha coefficients
    """
    order = q - 2*k + p
    sum_limit = math.floor((degree - order)/2)
    result = 0
    for j in range(k, sum_limit + 1):
        result += math.comb(j, k) * helper_function(degree, order, j)

    return result


def alpha(degree, p, q, k):
    """ Compute the alpha
    """
    i_and_minus_one = (0. + -1.j)**(q - 2*k)
    order = q - 2*k + p
    norm = norm_sh(degree, order)
    sum_in_alpha = sum_in_alphas(degree, p, q, k)
    result = i_and_minus_one * norm * math.factorial(q) \
        * math.factorial(degree - q - p) * sum_in_alpha \
        * (4*np.pi)/(2*degree + 1)
    if p == 1:
        result *= order

    return result


# A (NOTE: not A^{*})
def alpha_matrix(degree):
    """ Assemble the alpha matrix
    """
    matrix = np.zeros((2*degree + 1, 2*degree + 1), complex)
    for q in reversed(range(degree + 1)):
        for k in range(math.floor(q/2) + 1):
            p = 0
            order = q - 2*k + p
            # m >= 0
            matrix[degree - q, degree - order] = alpha(degree, p, q, k)
            # m < 0 (NOTE: I am writing the m = 0 case twice)
            matrix[degree - q, degree + order] = alpha(degree, p, q, k)

            p = 1
            order = q - 2*k + p
            if q == degree:
                continue
            # m >= 0
            matrix[degree + q + 1, degree - order] = alpha(degree, p, q, k)
            # m < 0
            matrix[degree + q + 1, degree + order] = -alpha(degree, p, q, k)
    return matrix

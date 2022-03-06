#ifndef SQUARE_MATRIX_H
#define SQUARE_MATRIX_H

#include <cstddef>
#include <complex>
#include <iomanip>
#include <ios>
#include <iostream>
#include <type_traits>
#include <vector>

// Type traits
template<typename T>
struct is_complex : public std::false_type {};

template<typename T>
struct is_complex<std::complex<T>> : public std::true_type {};

// Square matrix
template <typename T>
class SquareMatrix {
  // class invariant: creates a dense square matrix
 public:
  using value_type = T;
  
  SquareMatrix(std::size_t n)
      : matrix_dim{n}, matrix_elements(matrix_dim * matrix_dim) {}

  T& operator()(std::size_t i, std::size_t j) {
    return matrix_elements[i * matrix_dim + j];
  }

  const T& operator()(std::size_t i, std::size_t j) const {
    return matrix_elements[i * matrix_dim + j];
  }

  void print();
  std::size_t dim() const { return matrix_dim; };
  std::size_t size() const { return matrix_elements.size(); }

  // Scalar arithmetic operations
  // Should these be member functions?
  SquareMatrix operator*(const T& scalar);
  SquareMatrix operator/(const T& scalar);

  // Transpose
  SquareMatrix& transpose();
  // Conjugate
  std::enable_if<std::is_class<T>(), T&> operator->();
 private:
  // data members
  std::size_t matrix_dim;  // quadratic matrix n = m
  std::vector<T> matrix_elements;
};

// Definition of print()
template <typename T>
void SquareMatrix<T>::print() {
  for (std::size_t i = 0; i < matrix_dim; ++i) {
    std::cout << std::left;
    for (std::size_t j = 0; j < matrix_dim; ++j) {
      std::cout << std::setw(6) << matrix_elements[i * matrix_dim + j] << " ";
    }
    std::cout << "\n";
  }
}
// Definiton transpose
template <typename T>
SquareMatrix<T>& SquareMatrix<T>::transpose() {
  std::vector<T> temp_mat_elements(matrix_dim * matrix_dim);
  for (std::size_t i = 0; i < matrix_dim; ++i) {
    for (std::size_t j = 0; j < matrix_dim; ++j) {
      temp_mat_elements[j * matrix_dim + i] =
	  matrix_elements[i * matrix_dim + j];
    }
  }
  matrix_elements = temp_mat_elements;
  return *this;
}

// Mutliplication from right
template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator*(const T& scalar) {
  SquareMatrix<T> temp_mat{*this};
  for (auto& element : temp_mat.matrix_elements) element *= scalar;
  return temp_mat;
}

// Multiplication from left
template <typename T>
SquareMatrix<T> operator*(T scalar, const SquareMatrix<T>& mat) {
  return mat * scalar;
}

// Division from right
template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator/(const T& scalar) {
  SquareMatrix<T> temp_mat{*this};
  for (auto& element : temp_mat.matrix_elements) element /= scalar;
  return temp_mat;
}

// Division from left
template <typename T>
SquareMatrix<T> operator/(T scalar, const SquareMatrix<T>& mat) {
  return mat / scalar;
}

// Matrix multiplication
template <typename T>
SquareMatrix<T> operator*(const SquareMatrix<T>& matA,
                          const SquareMatrix<T>& matB) {
  SquareMatrix<T> res(matA.dim());
  for (std::size_t i = 0; i < matA.dim(); ++i) {
    for (std::size_t j = 0; j < matA.dim(); ++j) {
      for (std::size_t k = 0; k < matA.dim(); ++k) {
        res(i, j) += matA(i, k) * matB(k, j);
      }
    }
  }
  return res;
}

#endif

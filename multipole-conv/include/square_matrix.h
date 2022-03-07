#ifndef SQUARE_MATRIX_H
#define SQUARE_MATRIX_H

#include <complex>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <type_traits>
#include <vector>

// Type traits
template <typename T>
struct is_complex : public std::false_type {};

template <typename T>
struct is_complex<std::complex<T>> : public std::true_type {};

template <typename DataTypeA, typename DataTypeB>
struct MatrixReturn {
  // None of the types is complex or both. Both types can be the return type
  typedef DataTypeA type;
};

// Specialisation Type A is complex
template <typename DataTypeA, typename DataTypeB>
struct MatrixReturn<std::complex<DataTypeA>, DataTypeB> {
  typedef std::complex<DataTypeA> type;
};

// Specialisation Type B is complex
template <typename DataTypeA, typename DataTypeB>
struct MatrixReturn<DataTypeA, std::complex<DataTypeB>> {
  typedef std::complex<DataTypeB> type;
};

// Square matrix
template <typename T>
class SquareMatrix {
  // class invariant: creates a dense square matrix. Only double and float are
  // allowed. Allows to mix operations between complex<T> and T. But does not
  // allow to mix double and float.
 public:
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

  // Transpose
  SquareMatrix transpose() const;
  // Conjugate
  template <typename Q = T>
  // Only create this function if T is complex
  std::enable_if_t<is_complex<Q>::value, SquareMatrix<Q>> conjugate() const;

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
      std::cout << std::setw(8) << matrix_elements[i * matrix_dim + j] << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

// Definiton transpose
template <typename T>
SquareMatrix<T> SquareMatrix<T>::transpose() const {
  SquareMatrix<T> temp_matrix(matrix_dim);
  for (std::size_t i = 0; i < matrix_dim; ++i) {
    for (std::size_t j = 0; j < matrix_dim; ++j) {
      temp_matrix.matrix_elements[j * matrix_dim + i] =
          matrix_elements[i * matrix_dim + j];
    }
  }
  return temp_matrix;
}

// Definition conjugate
template <typename T>
template <typename Q>
std::enable_if_t<is_complex<Q>::value, SquareMatrix<Q>>
SquareMatrix<T>::conjugate() const {
  SquareMatrix<Q> temp_matrix(matrix_dim);
  for (std::size_t i = 0; i < matrix_elements.size(); ++i)
    temp_matrix.matrix_elements[i] = std::conj(matrix_elements[i]);
  return temp_matrix;
}

// Scalar mutliplication from right
template <typename T, typename S>
SquareMatrix<typename MatrixReturn<T, S>::type> operator*(
    const SquareMatrix<T>& mat, const S& scalar) {
  SquareMatrix<typename MatrixReturn<T, S>::type> res(mat.dim());
  for (std::size_t i = 0; i < mat.dim(); ++i) {
    for (std::size_t j = 0; j < mat.dim(); ++j) {
      res(i, j) = mat(i, j) * scalar;
    }
  }
  return res;
}

// Scalar multiplication from left
template <typename S, typename T>
SquareMatrix<typename MatrixReturn<S, T>::type> operator*(
    const S& scalar, const SquareMatrix<T>& mat) {
  return mat * scalar;
}

// Scalar division from right
template <typename T, typename S>
SquareMatrix<typename MatrixReturn<T, S>::type> operator/(
    const SquareMatrix<T>& mat, const S& scalar) {
  SquareMatrix<typename MatrixReturn<T, S>::type> res(mat.dim());
  for (std::size_t i = 0; i < mat.dim(); ++i) {
    for (std::size_t j = 0; j < mat.dim(); ++j) {
      res(i, j) = mat(i, j) / scalar;
    }
  }
  return res;
}
// Scalar division from left is from a mathematical perspective strange
template <typename T, typename S>
SquareMatrix<typename MatrixReturn<T, S>::type> operator*(
    const SquareMatrix<T>& matA, const SquareMatrix<S>& matB) {
  SquareMatrix<typename MatrixReturn<T, S>::type> res(matA.dim());
  for (std::size_t i = 0; i < matA.dim(); ++i) {
    for (std::size_t j = 0; j < matA.dim(); ++j) {
      for (std::size_t k = 0; k < matB.dim(); ++k) {
        res(i, j) += matA(i, k) * matB(k, j);
      }
    }
  }
  return res;
}

#endif

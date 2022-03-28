#ifndef MATRIX_H
#define MATRIX_H

#include <complex>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <iostream>
#include <type_traits>
#include <vector>

namespace multipole_conv {
// Type traits
template <typename T>
struct is_complex : public std::false_type {};

template <typename T>
struct is_complex<std::complex<T>> : public std::true_type {};

template <typename DataTypeA, typename DataTypeB>
struct MatrixReturn {
  // None of the types is complex
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

// Specialisation both are complex
template <typename DataTypeA, typename DataTypeB>
struct MatrixReturn<std::complex<DataTypeA>, std::complex<DataTypeB>> {
  typedef std::complex<DataTypeA> type;
};

template <typename T>
class Matrix {
  // class invariant: creates a dense matrix. Only double and float are
  // allowed. Allows to mix operations between complex<T> and T. But does not
  // allow to mix double and float.
 public:
  using value_type = T;

  Matrix(std::size_t n)
      : n_rows{n}, n_columns(n), matrix_elements(n_rows * n_columns) {}

  Matrix(std::size_t n, std::size_t m)
      : n_rows{n}, n_columns{m}, matrix_elements(n_rows * n_columns) {}

  // TODO: No error handling for copy constructor, copy assignment operator and
  // move constructor and move assignment operator. Dimensions of the matrix are
  // not checked.
  T& operator()(std::size_t i, std::size_t j) {
    return matrix_elements[i * n_columns + j];
  }

  const T& operator()(std::size_t i, std::size_t j) const {
    return matrix_elements[i * n_columns + j];
  }

  // Allow to convert a matrix with double entries to complex matrix in an
  // assignment operation
  template <typename Q = T>
  std::enable_if_t<is_complex<Q>::value, Matrix<Q>> operator=(
      const Matrix<double>& rhs) {
    
  }

  void print();
  std::size_t rows() const { return n_rows; };
  std::size_t columns() const { return n_columns; };
  std::size_t size() const { return matrix_elements.size(); }

  // Transpose
  Matrix transpose() const;
  // Conjugate
  template <typename Q = T>
  // Only create this function if T is complex
  std::enable_if_t<is_complex<Q>::value, Matrix<Q>> conjugate() const;

 private:
  // data members
  std::size_t n_rows;
  std::size_t n_columns;
  std::vector<T> matrix_elements;
};

// Definition of print()
template <typename T>
void Matrix<T>::print() {
  for (std::size_t i = 0; i < n_rows; ++i) {
    std::cout << std::left;
    for (std::size_t j = 0; j < n_columns; ++j) {
      std::cout << std::setw(9) << std::setprecision(3)
                << matrix_elements[i * n_columns + j] << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

// Definiton transpose
template <typename T>
Matrix<T> Matrix<T>::transpose() const {
  Matrix<T> temp_matrix(n_columns, n_rows);
  for (std::size_t i = 0; i < n_rows; ++i) {
    for (std::size_t j = 0; j < n_columns; ++j) {
      temp_matrix.matrix_elements[j * n_rows + i] =
          matrix_elements[i * n_columns + j];
    }
  }
  return temp_matrix;
}

// Definition conjugate
template <typename T>
template <typename Q>
std::enable_if_t<is_complex<Q>::value, Matrix<Q>> Matrix<T>::conjugate() const {
  Matrix<Q> temp_matrix(n_rows, n_columns);
  for (std::size_t i = 0; i < matrix_elements.size(); ++i)
    temp_matrix.matrix_elements[i] = std::conj(matrix_elements[i]);
  return temp_matrix;
}

// Scalar mutliplication from right
template <typename T, typename S>
Matrix<typename MatrixReturn<T, S>::type> operator*(const Matrix<T>& mat,
                                                    const S& scalar) {
  Matrix<typename MatrixReturn<T, S>::type> res(mat.rows(), mat.columns());
  for (std::size_t i = 0; i < mat.rows(); ++i) {
    for (std::size_t j = 0; j < mat.columns(); ++j) {
      res(i, j) = mat(i, j) * scalar;
    }
  }
  return res;
}

// Scalar multiplication from left
template <typename S, typename T>
Matrix<typename MatrixReturn<S, T>::type> operator*(const S& scalar,
                                                    const Matrix<T>& mat) {
  return mat * scalar;
}

// Scalar division from right
template <typename T, typename S>
Matrix<typename MatrixReturn<T, S>::type> operator/(const Matrix<T>& mat,
                                                    const S& scalar) {
  Matrix<typename MatrixReturn<T, S>::type> res(mat.rows(), mat.columns());
  for (std::size_t i = 0; i < mat.rows(); ++i) {
    for (std::size_t j = 0; j < mat.columns(); ++j) {
      res(i, j) = mat(i, j) / scalar;
    }
  }
  return res;
}
// Scalar division from left is from a mathematical a scalar divided by a
// matrix, i.e. an undefined operation. Hence, I do not define it.

// MatrixMatrix mult
template <typename T, typename S>
Matrix<typename MatrixReturn<T, S>::type> operator*(const Matrix<T>& matA,
                                                    const Matrix<S>& matB) {
  // TODO: No error handling, namely it is not checked if the dimension of matA
  // and matB are suited for MatrixMatrix multiplication (matA.columns() ==
  // matB.rows())
  Matrix<typename MatrixReturn<T, S>::type> res(matA.rows(), matB.columns());
  for (std::size_t i = 0; i < matA.rows(); ++i) {
    for (std::size_t j = 0; j < matB.columns(); ++j) {
      for (std::size_t k = 0; k < matB.rows(); ++k) {
        res(i, j) += matA(i, k) * matB(k, j);
      }
    }
  }
  return res;
}

// MatrixVector mult
template <typename T, typename S>
std::vector<typename MatrixReturn<T, S>::type> operator*(
    const Matrix<T>& mat, const std::vector<S>& vec) {
  // TODO: No error handling, namely it is not checked if dimension of matA and
  // the dimension of vec are suited for MatrixVector multiplication
  // (matA.rows() == vec.size())
  std::vector<typename MatrixReturn<T, S>::type> res(vec.size());
  for (size_t i = 0; i < mat.rows(); ++i) {
    for (size_t j = 0; j < vec.size(); ++j) {
      res[i] += mat(i, j) * vec[j];
    }
  }
  return res;
}
}  // namespace multipole_conv
#endif

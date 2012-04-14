/* Copyright (C) 2003-2005  Martin Förg

 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact:
 *   martin.o.foerg@googlemail.com
 *
 */

#ifndef linear_algebrad_h
#define linear_algebrad_h

#include "square_matrix.h"
#include "vector.h"
#include "row_vector.h"
#include "diagonal_matrix.h"
#include "symmetric_matrix.h"
#include "band_matrix.h"
#include <complex>

//-------------------------------------
// Matrix operations
//-------------------------------------

namespace fmatvec {

  /*! \brief Matrix-scalar multiplication.
   *
   * This function multiplies a diagonal matrix by a scalar. 
   * \param A A diagonal matrix. 
   * \param alpha A scalar. 
   * \return A reference to the diagonal matrix. 
   * */
  Matrix<Diagonal, double >& operator*=(const Matrix<Diagonal, double >& A, double alpha);

  /*! \brief Matrix-scalar division.
   *
   * This function divides a diagonal matrix by a scalar. 
   * \param A A diagonal matrix. 
   * \param alpha A scalar. 
   * \return A reference to the diagonal matrix. 
   * */
  Matrix<Diagonal, double >& operator/=(const Matrix<Diagonal, double >& A, double alpha);

  /*! \brief Matrix-scalar multiplication.
   *
   * This function multiplies a square matrix by a scalar. 
   * \param A A square matrix. 
   * \param alpha A scalar. 
   * \return A reference to the square matrix. 
   * */
  SquareMatrix<double>& operator*=(const SquareMatrix<double>& A, double alpha);

  /*! \brief Matrix-scalar division.
   *
   * This function divides a square matrix by a scalar. 
   * \param A A square matrix. 
   * \param alpha A scalar. 
   * \return A reference to the square matrix. 
   * */
  SquareMatrix<double>& operator/=(const SquareMatrix<double>& A, double alpha);

  /*! \brief Matrix-scalar multiplication.
   *
   * This function computes the product of a square matrix and a scalar. 
   * \param A A square matrix. 
   * \param alpha A scalar. 
   * \return A new square matrix containig the result.
   * */
  SquareMatrix<double> operator*(const SquareMatrix<double> &A, double alpha);

  /*! \brief Scalar-matrix multiplication.
   *
   * \see operator*(const SquareMatrix<double> &A, double alpha).
   * */
  SquareMatrix<double> operator*(double alpha, const SquareMatrix<double> &A);

  /*! \brief Matrix-scalar division.
   *
   * This function computes the division of a square matrix and a scalar. 
   * \param A A square matrix. 
   * \param alpha A scalar. 
   * \return A new square matrix containig the result.
   * */
  SquareMatrix<double> operator/(const SquareMatrix<double> &A, double alpha);

 /*! \brief Matrix-scalar multiplication.
   *
   * This function multiplies a general matrix by a scalar. 
   * \param A A general matrix. 
   * \param alpha A scalar. 
   * \return A reference to the general matrix. 
   * */
  Matrix<General, double >& operator*=(const Matrix<General, double > &A, const double &alpha);

  /*! \brief Matrix-scalar division.
   *
   * This function divides a general matrix by a scalar. 
   * \param A A general matrix. 
   * \param alpha A scalar. 
   * \return A reference to the general matrix. 
   * */
  Matrix<General, double >& operator/=(const Matrix<General, double > &A, const double &alpha);

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two general matrices. 
   * \param A A general matrix. 
   * \param B A general matrix. 
   * \return A new general matrix containig the result.
   * */
  Matrix<General, double > operator+(const Matrix<General, double > &A, const Matrix<General, double > &B);

  /*! \brief Matrix-matrix subtraction.
   *
   * This function computes the difference of two general matrices
   * according to \f[\boldsymbol{A}-\boldsymbol{B} \f]
   * \param A A general matrix. 
   * \param B A general matrix. 
   * \return A new general matrix containig the result.
   * */
  Matrix<General, double> operator-(const Matrix<General, double> &A, const Matrix<General, double> &B);

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two square matrices. 
   * \param A A square matrix. 
   * \param B A square matrix. 
   * \return A new square matrix containig the result.
   * */
  SquareMatrix<double> operator+(const SquareMatrix<double> &A, const SquareMatrix<double> &B);

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two symmetric matrices. 
   * \param A A symmetric matrix. 
   * \param B A symmetric matrix. 
   * \return A new symmetric matrix containig the result.
   * */
  Matrix<Symmetric, double> operator+(const Matrix<Symmetric, double > &A, const Matrix<Symmetric, double > &B);

  /*! \brief Matrix-matrix subtraction.
   *
   * This function computes the difference of two square matrices. 
   * \param A A square matrix. 
   * \param B A square matrix. 
   * \return A new square matrix containig the result.
   * */
  SquareMatrix<double> operator-(const SquareMatrix<double> &A, const SquareMatrix<double> &B);

   /*! \brief Matrix-matrix subtraction.
   *
   * This function computes the sum of two symmetric matrices. 
   * \param A A symmetric matrix. 
   * \param B A symmetric matrix. 
   * \return A new symmetric matrix containig the result.
   * */
  Matrix<Symmetric, double> operator-(const Matrix<Symmetric, double > &A, const Matrix<Symmetric, double > &B);

  /*! \brief Matrix-scalar multiplication.
   *
   * This function computes the product of a general matrix and a scalar. 
   * \param A A general matrix. 
   * \param alpha A scalar. 
   * \return A new general matrix containig the result.
   * */
  Matrix<General, double> operator*(const Matrix<General, double> &A, double alpha);

  /*! \brief Scalar-matrix multiplication.
   *
   * \see operator*(const Matrix<General, double> &A, double alpha).
   * */
  Matrix<General, double> operator*(double alpha, const Matrix<General, double> &A);

  /*! \brief Matrix-scalar division.
   *
   * This function computes the division of a general matrix and a scalar. 
   * \param A A general matrix. 
   * \param alpha A scalar. 
   * \return A new general matrix containig the result.
   * */
  Matrix<General, double> operator/(const Matrix<General, double> &A, double alpha);

  /*! \brief Matrix-scalar multiplication.
   *
   * This function computes the product of a symmetric matrix and a scalar. 
   * \param A A symmetric matrix. 
   * \param alpha A scalar. 
   * \return A new symmetric matrix containig the result.
   * */
  Matrix<Symmetric, double> operator*(const Matrix<Symmetric, double> &A, double alpha);

   /*! \brief Scalar-matrix multiplication.
   *
   * \see operator*(const Matrix<Symmetric, double> &A, double alpha).
   * */
  Matrix<Symmetric, double> operator*(double alpha, const Matrix<Symmetric, double> &A);

  /*! \brief Matrix-scalar division.
   *
   * This function computes the division of a symmetric matrix and a scalar. 
   * \param A A symmetric matrix. 
   * \param alpha A scalar. 
   * \return A new symmetric matrix containig the result.
   * */
  Matrix<Symmetric, double> operator/(const Matrix<Symmetric, double> &A, double alpha);

  /*! \brief Matrix-scalar multiplication.
   *
   * This function computes the product of a diagonal matrix and a scalar. 
   * \param A A diagonal matrix. 
   * \param alpha A scalar. 
   * \return A new diagonal matrix containig the result.
   * */
  Matrix<Diagonal, double> operator*(const Matrix<Diagonal, double> &A, double alpha);

 /*! \brief Scalar-matrix multiplication.
   *
   * \see operator*(const Matrix<Diagonal, double> &D, double alpha).
   * */
  Matrix<Diagonal, double> operator*(double alpha, const Matrix<Diagonal, double> &A);

   /*! \brief Matrix-scalar division.
   *
   * This function computes the division of a diagonal matrix and a scalar. 
   * \param A A diagonal matrix. 
   * \param alpha A scalar. 
   * \return A new diagonal matrix containig the result.
   * */
  Matrix<Diagonal, double> operator/(const Matrix<Diagonal, double> &A, double alpha);

   /*! \brief Special multiplication.
   *
   * This function computes the product of a general matrix 
   * according to \f[ \boldsymbol{A}^T\,\boldsymbol{A} \f]
   * \param A A general matrix. 
   * \return A new symmetric matrix containig the result.
   * */
  Matrix<Symmetric, double> JTJ(const Matrix<General, double> &A);

   /*! \brief Special multiplication.
   *
   * This function computes the product of a symmetric and a general matrix
   * according to \f[ \boldsymbol{A}^T\,\boldsymbol{B}\,\boldsymbol{A} \f]
   * \param A A general matrix. 
   * \param B A symmetric matrix. 
   * \return A new symmetric matrix containig the result.
   * */
  Matrix<Symmetric, double> JTMJ(const Matrix<Symmetric, double> &B, const Matrix<General, double> &A);

   /*! \brief Special multiplication.
   *
   * This function computes the product of a diagonal and a symmetric matrix
   * according to \f[ \boldsymbol{A}^T\,\boldsymbol{B}\,\boldsymbol{A} \f]
   * \param A A general matrix. 
   * \param B A diagonal matrix. 
   * \return A new symmetric matrix containig the result.
   * */
  Matrix<Symmetric, double> JTMJ(const Matrix<Diagonal, double> &B, const Matrix<General, double> &A);

 /*! \brief Matrix-matrix addition.
   *
   * This function adds a diagonal matrix to a diagonal matrix. 
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}+\boldsymbol{B} \f]
   * \param A A diagonal matrix. 
   * \param B A diagonal matrix. 
   * \return A reference to the first diagonal matrix. 
   * */
  Matrix<Diagonal, double >& operator+=(const Matrix<Diagonal, double >& A, const Matrix<Diagonal, double > &B);

 /*! \brief Matrix-matrix subtraction.
   *
   * This function subtracts a diagonal matrix from a diagonal matrix. 
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}-\boldsymbol{B} \f]
   * \param A A diagonal matrix. 
   * \param B A diagonal matrix. 
   * \return A reference to the first matrix.
   * */
  Matrix<Diagonal, double >& operator-=(const Matrix<Diagonal, double >& A, const Matrix<Diagonal, double > &B);

 /*! \brief Matrix-matrix addition.
   *
   * This function adds a symmetric matrix to a symmetric matrix. 
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}+\boldsymbol{B} \f]
   * \param A A symmetric matrix. 
   * \param B A symmetric matrix. 
   * \return A reference to the first matrix. 
   * */
  Matrix<Symmetric, double>& operator+=(const Matrix<Symmetric, double >& A, const Matrix<Symmetric, double> &B);

 /*! \brief Matrix-matrix subtraction.
   *
   * This function subtracts a symmetric matrix from a symmetric matrix. 
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}-\boldsymbol{B} \f]
   * \param A A symmetric matrix. 
   * \param B A symmetric matrix. 
   * \return A reference to the first matrix.
   * */
  Matrix<Symmetric, double>& operator-=(const Matrix<Symmetric, double >& A, const Matrix<Symmetric, double> &B);

 /*! \brief Matrix-matrix addition.
   *
   * This function adds a square matrix to a square matrix. 
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}+\boldsymbol{B} \f]
   * \param A A square matrix. 
   * \param B A square matrix. 
   * \return A reference to the first matrix. 
   * */
  SquareMatrix<double>& operator+=(const SquareMatrix<double>& A, const SquareMatrix<double> &B);

 /*! \brief Matrix-matrix subtraction.
   *
   * This function subtracts a square matrix from a square matrix. 
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}-\boldsymbol{B} \f]
   * \param A A square matrix. 
   * \param B A square matrix. 
   * \return A reference to the first matrix.
   * */
  SquareMatrix<double>& operator-=(const SquareMatrix<double>& A, const SquareMatrix<double> &B);

 /*! \brief Matrix-matrix addition.
   *
   * This function adds a general matrix to a general matrix. 
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}+\boldsymbol{B} \f]
   * \param A A general matrix. 
   * \param B A general matrix. 
   * \return A reference to the first matrix. 
   * */
  Matrix<General, double >& operator+=(const Matrix<General, double > &A, const Matrix<General, double > &B);

 /*! \brief Matrix-matrix subtraction.
   *
   * This function subtracts a general matrix from a general matrix
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}-\boldsymbol{B} \f]
   * \param A A general matrix. 
   * \param B A general matrix. 
   * \return A reference to the first matrix.
   * */
  Matrix<General, double >& operator-=(const Matrix<General, double > &A, const Matrix<General, double > &B);

  /*! \brief Matrix-matrix multiplication.
   *
   * This function computes the product of two general matrices
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}\,\boldsymbol{B} \f]
   * \param A A general matrix. 
   * \param B A general matrix. 
   * \return A new general matrix containig the result.
   * */
  Matrix<General, double> operator*(const Matrix<General, double> &A, const Matrix<General, double> &B);

  /*! \brief Matrix-matrix multiplication.
   *
   * This function computes the product of two square matrices
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}\,\boldsymbol{B} \f]
   * \param A A square matrix. 
   * \param B A square matrix. 
   * \return A new square matrix containig the result.
   * */
  SquareMatrix<double> operator*(const SquareMatrix<double> &A, const SquareMatrix<double> &B);

  /*! \brief Matrix-matrix multiplication.
   *
   * This function computes the product of a general and a symmetric matrix.
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}\,\boldsymbol{B} \f]
   * \param A A general matrix. 
   * \param B A symmetric matrix. 
   * \return A new general matrix containig the result.
   * */
  Matrix<General, double> operator*(const Matrix<General, double> &A, const Matrix<Symmetric, double> &B);

  /*! \brief Scalar-matrix multiplication.
   *
   * \see operator*(const Matrix<General, double> &A, const Matrix<Symmetric, double> &B)
   * */
  Matrix<General, double> operator*(const Matrix<Symmetric, double> &A, const Matrix<General, double> &B);

  /*! \brief Matrix-matrix multiplication.
   *
   * This function computes the product of a general and a diagonal matrix.
   * according to \f[\boldsymbol{A} \rightarrow \boldsymbol{A}\,\boldsymbol{B} \f]
   * \param A A general matrix. 
   * \param B A diagonal matrix. 
   * \return A new general matrix containig the result.
   * */
  Matrix<General, double> operator*(const Matrix<General, double> &A,const Matrix<Diagonal, double> &B);

  /*! \brief Scalar-matrix multiplication.
   *
   * \see operator*(const Matrix<General, double> &A,const Matrix<Diagonal, double> &B)
   * */
  Matrix<General, double> operator*(const Matrix<Diagonal, double> &B, const Matrix<General, double> &A);

/*! \brief Vector-scalar multiplication.
   *
   * This function multiplies a vector by a scalar. 
   * \param x A vector.
   * \param alpha A scalar. 
   * \return A reference to the vector.
   * */
  Vector<double>& operator*=(const Vector<double> &x, const double& alpha);

 /*! \brief Vector-scalar multiplication.
   *
   * This function computes the product of a vector and a scalar. 
   * \param x A vector.
   * \param alpha A scalar. 
   * \return A new vector containig the result.
   * */
  Vector<double> operator*(const Vector<double> &x, double alpha);

  /*! \brief Scalar-vector multiplication.
   *
   * \see operator*(const Vector<double> &x, double alpha)
   * */
  Vector<double> operator*(double alpha, const Vector<double> &x);

  /*! \brief Vector-scalar division.
   *
   * This function divides a vector by a scalar. 
   * \param x A vector.
   * \param alpha A scalar. 
   * \return A reference to the vector.
   * */
  Vector<double>& operator/=(const Vector<double> &x, const double &alpha);

 /*! \brief Vector-scalar division.
   *
   * This function computes the division of a vector and a scalar. 
   * \param x A vector.
   * \param alpha A scalar. 
   * \return A new vector containig the result.
   * */
  Vector<double> operator/(const Vector<double> &x, double alpha);

  /*! \brief Vector-vector addition.
   *
   * This function adds two vectors 
   * according to \f[\boldsymbol{x} \rightarrow \boldsymbol{x}+\boldsymbol{y} \f]
   * \param x A vector.
   * \param y A vector.
   * \return A reference to the first vector.
   * */
  Vector<double>& operator+=(const Vector<double> &x, const Vector<double> &y);

  /*! \brief Vector-vector subtraction.
   *
   * This function subtracts two vectors 
   * according to \f[\boldsymbol{x} \rightarrow \boldsymbol{x}-\boldsymbol{y} \f]
   * \param x A vector.
   * \param y A vector.
   * \return A reference to the first vector.
   * */
  Vector<double>& operator-=(const Vector<double> &x, const Vector<double> &y);

  /*! \brief Vector-vector addition.
   *
   * This function computes the sum of two vectors
   * according to \f[\boldsymbol{x}+\boldsymbol{y} \f]
   * \param x A vector.
   * \param y A vector.
   * \return A new vector containig the result.
   * */
  //Vector<double> operator+(const Vector<double> &x, const Vector<double> &y);

  /*! \brief Vector-vector subtraction.
   *
   * This function computes the difference of two vectors
   * according to \f[\boldsymbol{x}-\boldsymbol{y} \f]
   * \param x A vector.
   * \param y A vector.
   * \return A new vector containig the result.
   * */
  Vector<double> operator-(const Vector<double> &x, const Vector<double> &y);

/*! \brief vector-scalar multiplication.
   *
   * This function multiplies a row vector by a scalar. 
   * \param x A row vector.
   * \param alpha A scalar. 
   * \return A reference to the row vector.
   * */
  RowVector<double>& operator*=(const RowVector<double> &x, const double& alpha);

  /*! \brief Vector-scalar division.
   *
   * This function divides a row vector by a scalar. 
   * \param x A row vector.
   * \param alpha A scalar. 
   * \return A reference to the vector.
   * */
  RowVector<double>& operator/=(const RowVector<double> &x, const double &alpha);

  /*! \brief Vector-scalar multiplication.
   *
   * This function computes the product of a vector and a scalar. 
   * \param x A row vector.
   * \param alpha A scalar. 
   * \return A new row vector containig the result.
   * */
  RowVector<double> operator*(double alpha, const RowVector<double> &x);

  /*! \brief Scalar-vector multiplication.
   *
   * \see operator*(const RowVector<double> &x, double alpha)
   * */
  RowVector<double> operator*(const RowVector<double> &x, double alpha);

  /*! \brief Vector-scalar division.
   *
   * This function computes the division of a row vector and a scalar. 
   * \param x A row vector.
   * \param alpha A scalar. 
   * \return A new row vector containig the result.
   * */
  RowVector<double> operator/(const RowVector<double> &x, double alpha);

  /*! \brief Vector-vector addition.
   *
   * This function adds two row vectors 
   * according to \f[\boldsymbol{x} \rightarrow \boldsymbol{x}+\boldsymbol{y} \f]
   * \param x A row vector.
   * \param y A row vector.
   * \return A reference to the first row vector.
   * */
  RowVector<double>& operator+=(const RowVector<double> &x, const RowVector<double> &y);

  /*! \brief Vector-vector subtraction.
   *
   * This function subtracts two row vectors 
   * according to \f[\boldsymbol{x} \rightarrow \boldsymbol{x}-\boldsymbol{y} \f]
   * \param x A row vector.
   * \param y A row vector.
   * \return A reference to the first row vector.
   * */
  RowVector<double>& operator-=(const RowVector<double> &x, const RowVector<double> &y);

  /*! \brief Vector-vector addition.
   *
   * This function computes the sum of two row vectors
   * according to \f[\boldsymbol{x}+\boldsymbol{y} \f]
   * \param x A row vector.
   * \param y A row vector.
   * \return A new vector containig the result.
   * */
  RowVector<double> operator+(const RowVector<double> &x, const RowVector<double> &y);

  /*! \brief Vector-vector subtraction.
   *
   * This function computes the difference of two row vectors
   * according to \f[\boldsymbol{x}-\boldsymbol{y} \f]
   * \param x A row vector.
   * \param y A row vector.
   * \return A new vector containig the result.
   * */
  RowVector<double> operator-(const RowVector<double> &x, const RowVector<double> &y);

  /*! \brief Matrix-vector multiplication.
   *
   * This function computes the product of a general matrix and a vector
   * according to \f[\boldsymbol{A}\,\boldsymbol{x} \f]
   * \param A A general matrix.
   * \param x A vector. 
   * \return A new vector containig the result.
   * */
  //Vector<double> operator*(const Matrix<General, double> &A, const Vector<double> &x);

  /*! \brief Matrix-vector multiplication.
   *
   * This function computes the product of a symmetric matrix and a vector
   * according to \f[\boldsymbol{A}\,\boldsymbol{x} \f]
   * \param A A symmetric matrix.
   * \param x A vector. 
   * \return A new vector containig the result.
   * */
  Vector<double> operator*(const Matrix<Symmetric, double> &A, const Vector<double> &x);

  /*! \brief Matrix-vector multiplication.
   *
   * This function computes the product of a diagonal matrix and a vector
   * according to \f[\boldsymbol{A}\,\boldsymbol{x} \f]
   * \param A A diagonal matrix.
   * \param x A vector. 
   * \return A new vector containig the result.
   * */
  Vector<double> operator*(const Matrix<Diagonal, double> &A, const Vector<double> &x);

  /*! \brief Vector-matrix multiplication.
   *
   * This function computes the product of a row vector and a general matrix
   * according to \f[\boldsymbol{x}\,\boldsymbol{A} \f]
   * \param x A row vector. 
   * \param A A general matrix.
   * \return A new row vector containig the result.
   * */
  RowVector<double> operator*(const RowVector<double> &x, const Matrix<General, double> &A); 

  /*! \brief Vector-matrix multiplication.
   *
   * This function computes the product of a row vector and a symmetric matrix
   * according to \f[\boldsymbol{x}\,\boldsymbol{A} \f]
   * \param x A row vector. 
   * \param A A symmetric matrix.
   * \return A new row vector containig the result.
   * */
  RowVector<double> operator*(const RowVector<double> &x, const Matrix<Symmetric, double> &A);

  /*! \brief Vector-matrix multiplication.
   *
   * This function computes the product of a row vector and a diagonal matrix
   * according to \f[\boldsymbol{x}\,\boldsymbol{A} \f]
   * \param x A row vector. 
   * \param A A diagonal matrix.
   * \return A new row vector containig the result.
   * */
  RowVector<double> operator*(const RowVector<double> &x, const Matrix<Diagonal, double> &A);

  /*! \brief Vector-vector multiplication.
   *
   * This function computes the product of a row vector and a vector 
   * according to \f[\boldsymbol{x}\,\boldsymbol{y} \f]
   * \param x A row vector. 
   * \param y A vector. 
   * \return A scalar containig the result.
   * */
  double operator*(const RowVector<double> &x, const Vector<double> &y); 

  /*! \brief Eigenvalues
   *
   * This function computes the complex eigenvalues of a square matrix.
   * \param A A square matrix. 
   * \return A vector containig the eigenvalues.
   * */
  Vector<std::complex<double> > eigval(const SquareMatrix<double> &A); 

  /*! \brief Eigenvectors and Eigenvalues
   *
   * This function computes all the eigenvectors and the eigenvalues of a real generalized symmetric-definite eigenproblem, of the form A*x=(lambda)*B*x.
   * Here A and B are assumed to be symmetric and B is also positive definite.
   * \param A A symmetric matrix. 
   * \param B A symmetric, positive definite matrix.
   * \param eigenvector A square matrix in the dimension of A, containing the normalized Eigenvectors at the end of the function. 
   * \param eigenvalues A vector in the size of A, containing the Eigenvalues at the end of the function
   * \return void
   * */
  int eigvec(const Matrix<Symmetric, double> &A, const Matrix<Symmetric, double> &B, SquareMatrix<double> &eigenvectors, Vector<double> &eigenvalues); 

  /*! \brief Eigenvalues
   *
   * This function computes the eigenvalues of a symmetric matrix.
   * \param A A symmetric matrix. 
   * \return A vector containig the eigenvalues.
   * */
  Vector<double> eigval(const Matrix<Symmetric, double> &A); 

  /*! \brief Eigenvalues
   *
   * This function computes selected eigenvalues of a symmetric matrix.
   * The eigenvalues can be selected by specifying either a range of values or a 
   * range of indices for the desired eigenvalues.
   * \param A A symmetric matrix. 
   * \param il The index of the smallest eigenvalue to be returned
   * \param iu The index of the largest eigenvalue to be returned
   * \param abstol The absolute error tolerance for the eigenvalues
   * \return A vector containig the eigenvalues.
   * */
  Vector<double> eigvalSel(const Matrix<Symmetric, double> &A, int il, int iu, double abstol=0);

  /*! \brief Systems of linear equations
   *
   * This function solves systems of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]
   * by a LU decompostion.
   * \param A A square matrix. 
   * \param B A general matrix containing the right hand sides.
   * \return A general matrix containig the solution.
   * */
  Matrix<General, double> slvLU(const SquareMatrix<double> &A, const Matrix<General, double> &B);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a LU decompostion.
   * \param A A square matrix. 
   * \param b A vector containing the right hand side.
   * \return A vector containig the solution.
   * */
  Vector<double> slvLU(const SquareMatrix<double> &A, const Vector<double> &b);

  /*! \brief Systems of linear equations
   *
   * This function solves systems of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]. 
   * The matrix is already decomposed by a LU decompostion.
   * \param A A square matrix decomposed by a LU decompostion. 
   * \param B A general matrix containing the right hand sides.
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A general matrix containig the solution.
   * */
  Matrix<General, double> slvLUFac(const SquareMatrix<double> &A, const Matrix<General, double> &B, const Vector<int> &ipiv);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]. 
   * The matrix is already decomposed by a LU decompostion.
   * \param A A square matrix decomposed by a LU decompostion. 
   * \param b A vector containing the right hand side.
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A vector containig the solution.
   * */
  Vector<double> slvLUFac(const SquareMatrix<double> &A, const Vector<double> &b, const Vector<int> &ipiv);

  /*! \brief Systems of linear equations
   *
   * This function solves systems of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]
   * by a LL decompostion.
   * \param A A symmetric matrix. 
   * \param B A general matrix containing the right hand sides.
   * \return A general matrix containig the solution.
   * */
  Matrix<General, double> slvLL(const Matrix<Symmetric, double> &A, const Matrix<General, double> &B);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a LL decompostion.
   * \param A A symmetric matrix. 
   * \param b A vector containing the right hand side.
   * \return A vector containig the solution.
   * */
  Vector<double> slvLL(const Matrix<Symmetric, double> &A, const Vector<double> &b);


  /*! \brief Systems of linear equations
   *
   * This function solves systems of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]
   * by a QR decompostion.
   * \param A A square matrix. 
   * \param B A general matrix containing the right hand sides.
   * \return A general matrix containig the solution.
   * */
  Matrix<General, double> slvQR(const SquareMatrix<double> &A, const Matrix<General, double> &B);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a QR decompostion.
   * \param A A symmetric matrix. 
   * \param b A vector containing the right hand side.
   * \return A vector containig the solution.
   * */
  Vector<double> slvQR(const SquareMatrix<double> &A, const Vector<double> &b);

  /*! \brief Inverse
   *
   * This function computes the inverse of a general matrix 
   * according to \f[\boldsymbol{A}^{-1} \f]
   * \param A A square matrix. 
   * \return A square matrix containig the result.
   * */
  SquareMatrix<double> inv(const SquareMatrix<double> &A);

  /*! \brief Inverse
   *
   * This function computes the inverse of a symmetric matrix 
   * according to \f[\boldsymbol{A}^{-1} \f]
   * \param A A symmetric matrix. 
   * \return A symmetric matrix containig the result.
   * */
  Matrix<Symmetric, double> inv(const Matrix<Symmetric, double> &A);

  /*! \brief Inverse
   *
   * This function computes the inverse of a diagonal matrix 
   * according to \f[\boldsymbol{A}^{-1} \f]
   * \param A A diagonal matrix. 
   * \return A diagonal matrix containig the result.
   * */
  Matrix<Diagonal, double> inv(const Matrix<Diagonal, double> &A);

  /*! \brief LU decomposition
   *
   * This function computes the LU decomposition of a general matrix 
   * according to \f[\boldsymbol{A}=\boldsymbol{P}\,\boldsymbol{L}\,\boldsymbol{U} \f]
   * \param A A general matrix. 
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A general matrix containig the result.
   * */
  Matrix<General, double> facLU(const Matrix<General, double> &A, Vector<int> &ipiv);

  /*! \brief LU decomposition
   *
   * This function computes the LU decomposition of a square matrix 
   * according to \f[\boldsymbol{A}=\boldsymbol{P}\,\boldsymbol{L}\,\boldsymbol{U} \f]
   * \param A A square matrix. 
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A square matrix containig the result.
   * */
  SquareMatrix<double> facLU(const SquareMatrix<double> &A, Vector<int> &ipiv);

  /*! \brief LL decomposition
   *
   * This function computes the Cholesky decomposition of a symmetric matrix 
   * according to \f[\boldsymbol{A}=\boldsymbol{L}\,\boldsymbol{L} \f]
   * \param A A symmetric matrix. 
   * \return A symmetric matrix containig the result.
   * */
  Matrix<Symmetric, double> facLL(const Matrix<Symmetric, double> &A);

  /*! \brief 1-norm
   *
   * This function computes the sum of the absolute values of a vector.
   * \param x A vector.
   * \return A scalar containig the result.
   * */
  double nrm1(const Vector<double> &x);

  /*! \brief 2-norm
   *
   * This function computes the the Euclidean norm of a vector.
   * \param x A vector.
   * \return A scalar containig the result.
   * */
  double nrm2(const Vector<double> &x);

  /*! \brief Infinity-norm
   *
   * This function computes the the largest absolute value of a vector.
   * \param x A vector.
   * \return A scalar containig the result.
   * */
  double nrmInf(const Vector<double> &x);

  /*! \brief 1-norm
   *
   * This function computes the largest column sum of the absolute 
   * values of a general matrix.
   * \param A A general matrix.
   * \return A scalar containig the result.
   * */
  double nrm1(const Matrix<General,double> &A);

  /*! \brief 2-norm
   *
   * This function computes the largest singular value of a general matrix.
   * \param A A general matrix.
   * \return A scalar containig the result.
   * */
  double nrm2(const Matrix<General,double> &A);

  /*! \brief Infinity-norm
   *
   * This function computes the largest row sum of the absolute 
   * values of a general matrix.
   * \param A A general matrix.
   * \return A scalar containig the result.
   * */
  double nrmInf(const Matrix<General,double> &A);

  /*! \brief Frobenius-norm
   *
   * This function computes the Frobenius norm of a general matrix.
   * \param A A general matrix.
   * \return A scalar containig the result.
   * */
  double nrmFro(const Matrix<General,double> &A);

  /*! \brief Spectral radius
   *
   * This function computes the spectral radius of a square matrix.
   * \param A A square matrix.
   * \return A scalar containig the result.
   * */
  double rho(const SquareMatrix<double> &A);

  /*! \brief Spectral radius
   *
   * This function computes the spectral radius of a symmetric matrix.
   * \param A A symmetric matrix.
   * \return A scalar containig the result.
   * */
  double rho(const Matrix<Symmetric, double> &A);
  
  Vector<double> slvLLFac(const Matrix<Symmetric, double> &A, const Vector<double> &x);

  Matrix<General, double> slvLLFac(const Matrix<Symmetric, double> &A, const Matrix<General, double> &X);

  Matrix<General, double> slvLS(const Matrix<General, double> &A, const Matrix<General, double> &B, double rcond=-1);

  Vector<double> slvLS(const Matrix<General,double> &A, const Vector<double> &b, double rcond=-1);

  //Matrix<General, double> slvLU(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_DIAG unit, const SquareMatrix<double> &A, const Matrix<General, double> &X, const Vector<int> &ipiv );

  /*! \brief Row interchanges
   *
   * This function performs a series of row interchanges on a general matrix.
   * \param A A general matrix. 
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A Matrix containig the result.
   * */
  //Matrix<General, double> swap(const Matrix<General, double> &A, const Vector<int> &ipiv );
}

#endif

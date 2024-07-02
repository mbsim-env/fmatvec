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

  /*!
   * This funtion computes the SVD
   * */
  FMATVEC_EXPORT int svd(Matrix<General, Ref, Ref, double> &A,    Matrix<General, Ref, Ref, double> &S,  SquareMatrix<Ref, double> &U, SquareMatrix<Ref, double> &VT, int Rueckgabe);

  /*! \brief Eigenvalues
   *
   * This function computes the complex eigenvalues of a square matrix.
   * \param A A square matrix. 
   * \return A vector containig the eigenvalues.
   * */
  FMATVEC_EXPORT Vector<Ref, std::complex<double>> eigval(const SquareMatrix<Ref, double> &A);

  /*! \brief Eigenvectors and Eigenvalues
   *
   * This function computes the complex eigenvalues and eigenvectors of a square matrix.
   * \param A A square matrix. 
   * \param eigenvector A square matrix in the dimension of A, containing the normalized eigenvectors. 
   * \param eigenvalues A vector in the size of A, containing the eigenvalues.
   * \return If 0, successful exit. If -i, the i-th argument had an illegal value. If i, the QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed.
   * */
  FMATVEC_EXPORT int eigvec(const SquareMatrix<Ref, double> &A, SquareMatrix<Ref, std::complex<double>> &V, Vector<Ref, std::complex<double>> &w);

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
  FMATVEC_EXPORT int eigvec(const Matrix<Symmetric, Ref, Ref, double> &A, const Matrix<Symmetric, Ref, Ref, double> &B, SquareMatrix<Ref, double> &eigenvectors, Vector<Ref, double> &eigenvalues);

  /*! \brief Eigenvalues
   *
   * This function computes the eigenvalues of a symmetric matrix.
   * \param A A symmetric matrix. 
   * \return A vector containig the eigenvalues.
   * */
  FMATVEC_EXPORT Vector<Ref, double> eigval(const Matrix<Symmetric, Ref, Ref, double> &A);

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
  FMATVEC_EXPORT Vector<Ref, double> eigvalSel(const Matrix<Symmetric, Ref, Ref, double> &A, int il, int iu, double abstol = 0);

  /*!
   * \brief solve linear system with LU decomposed matrix
   * \param *LU       Pointer to the first element of the matrix whose elements form the lower and upper triangular matrix factors of A
   * \param *B        Pointer to the column vector, (n x 1) matrix, B.
   * \param  pivot[]  The i-th element is the pivot row interchanged with row i
   * \param *x        Solution to the equation Ax = B
   * \param  n        The number of rows or columns of the matrix LU
   * \return success of routine (1 = Success, -1 Failure (matrix LU is singular)
   *
   * This routine uses Doolittle's method to solve the linear equation
   * Ax = B.  This routine is called after the matrix A has been decomposed
   * into a product of a unit lower triangular matrix L and an upper
   * triangular matrix U with pivoting.  The argument LU is a pointer to the
   * matrix the subdiagonal part of which is L and the superdiagonal
   * together with the diagonal part is U. (The diagonal part of L is 1 and
   * is not stored.)   The matrix A = LU.
   * The solution proceeds by solving the linear equation Ly = B for y and
   * subsequently solving the linear equation Ux = y for x.
   *
   * REMARK: Copied and fitted on 28.01.2013 of http://www.mymathlib.com/c_source/matrices/linearsystems/doolittle_pivot.c
   */
  FMATVEC_EXPORT int Doolittle_LU_with_Pivoting_Solve(double *A, double B[], const int pivot[], double x[], int n);

  /*! \brief Systems of linear equations
   *
   * This function solves systems of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]
   * by a LU decompostion.
   * \param A A square matrix. 
   * \param B A general matrix containing the right hand sides.
   * \return A general matrix containig the solution.
   * */
  FMATVEC_EXPORT Matrix<General, Ref, Ref, double> slvLU(const SquareMatrix<Ref, double> &A, const Matrix<General, Ref, Ref, double> &X);
  FMATVEC_EXPORT Matrix<General, Var, Var, double> slvLU(const SquareMatrix<Var, double> &A, const Matrix<General, Var, Var, double> &X, int & info);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a LU decompostion.
   * \param A A square matrix. 
   * \param b A vector containing the right hand side.
   * \return A vector containig the solution.
   * */
  FMATVEC_EXPORT Vector<Ref, double> slvLU(const SquareMatrix<Ref, double> &A, const Vector<Ref, double> &x);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a LU decompostion.
   * \param A    A square matrix.
   * \param b    A vector containing the right hand side.
   * \param info Information about the success of the routine (0 = success)
   * \return A vector containig the solution.
   * */
  FMATVEC_EXPORT Vector<Ref, double> slvLU(const SquareMatrix<Ref, double> &A, const Vector<Ref, double> &x, int & info);
  FMATVEC_EXPORT Vector<Var, double> slvLU(const SquareMatrix<Var, double> &A, const Vector<Var, double> &x, int & info);

  /*! \brief Systems of linear equations
   *
   * This function solves systems of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]. 
   * The matrix is already decomposed by a LU decompostion using facLU.
   * \param A A square matrix decomposed by a LU decompostion using facLU. 
   * \param B A general matrix containing the right hand sides.
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A general matrix containig the solution.
   * */
  FMATVEC_EXPORT Matrix<General, Ref, Ref, double> slvLUFac(const SquareMatrix<Ref, double> &A, const Matrix<General, Ref, Ref, double> &X, const Vector<Ref, int> &ipiv);
  FMATVEC_EXPORT Matrix<General, Var, Var, double> slvLUFac(const SquareMatrix<Var, double> &A, const Matrix<General, Var, Var, double> &X, const Vector<Var, int> &ipiv);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]. 
   * The matrix is already decomposed by a LU decompostion using facLU.
   * \param A A square matrix decomposed by a LU decompostion using facLU. 
   * \param b A vector containing the right hand side.
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A vector containig the solution.
   * */
  FMATVEC_EXPORT Vector<Ref, double> slvLUFac(const SquareMatrix<Ref, double> &A, const Vector<Ref, double> &x, const Vector<Ref, int> &ipiv);
  FMATVEC_EXPORT Vector<Var, double> slvLUFac(const SquareMatrix<Var, double> &A, const Vector<Var, double> &x, const Vector<Var, int> &ipiv);

  /*! \brief Systems of linear equations
   *
   * This function solves systems of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]
   * by a LL decompostion.
   * \param A A symmetric matrix. 
   * \param B A general matrix containing the right hand sides.
   * \return A general matrix containig the solution.
   * */
  FMATVEC_EXPORT Matrix<General, Ref, Ref, double> slvLL(const Matrix<Symmetric, Ref, Ref, double> &A, const Matrix<General, Ref, Ref, double> &X);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a LL decompostion.
   * \param A A symmetric matrix. 
   * \param b A vector containing the right hand side.
   * \return A vector containig the solution.
   * */
  FMATVEC_EXPORT Vector<Ref, double> slvLL(const Matrix<Symmetric, Ref, Ref, double> &A, const Vector<Ref, double> &x);

  /*! \brief Systems of linear equations
   *
   * This function solves systems of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]
   * by a QR decompostion.
   * \param A A square matrix. 
   * \param B A general matrix containing the right hand sides.
   * \return A general matrix containig the solution.
   * */
  FMATVEC_EXPORT Matrix<General, Ref, Ref, double> slvQR(const SquareMatrix<Ref, double> &A, const Matrix<General, Ref, Ref, double> &X);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a QR decompostion.
   * \param A A symmetric matrix. 
   * \param b A vector containing the right hand side.
   * \return A vector containig the solution.
   * */
  FMATVEC_EXPORT Vector<Ref, double> slvQR(const SquareMatrix<Ref, double> &A, const Vector<Ref, double> &x);

  /*! \brief Inverse
   *
   * This function computes the inverse of a general matrix 
   * according to \f[\boldsymbol{A}^{-1} \f]
   * \param A A square matrix. 
   * \return A square matrix containig the result.
   * */
  FMATVEC_EXPORT SquareMatrix<Ref, double> inv(const SquareMatrix<Ref, double> &A);

  /*! \brief Inverse
   *
   * This function computes the inverse of a symmetric matrix 
   * according to \f[\boldsymbol{A}^{-1} \f]
   * \param A A symmetric matrix. 
   * \return A symmetric matrix containig the result.
   * */
  FMATVEC_EXPORT Matrix<Symmetric, Ref, Ref, double> inv(const Matrix<Symmetric, Ref, Ref, double> &A);

  /*! \brief Inverse
   *
   * This function computes the inverse of a diagonal matrix 
   * according to \f[\boldsymbol{A}^{-1} \f]
   * \param A A diagonal matrix. 
   * \return A diagonal matrix containig the result.
   * */
  FMATVEC_EXPORT Matrix<Diagonal, Ref, Ref, double> inv(const Matrix<Diagonal, Ref, Ref, double> &A);

  /*! \brief LU decomposition
   *
   * This function computes the LU decomposition of a general matrix 
   * according to \f[\boldsymbol{A}=\boldsymbol{P}\,\boldsymbol{L}\,\boldsymbol{U} \f]
   * \param A A general matrix. 
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A general matrix containig the result.
   * */
  FMATVEC_EXPORT Matrix<General, Ref, Ref, double> facLU(const Matrix<General, Ref, Ref, double> &A, Vector<Ref, int> &ipiv);

  /*! \brief LU decomposition
   *
   * This function computes the LU decomposition of a square matrix 
   * according to \f[\boldsymbol{A}=\boldsymbol{P}\,\boldsymbol{L}\,\boldsymbol{U} \f]
   * \param A A square matrix. 
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A square matrix containig the result.
   * */
  FMATVEC_EXPORT SquareMatrix<Ref, double> facLU(const SquareMatrix<Ref, double> &A, Vector<Ref, int> &ipiv);

  FMATVEC_EXPORT SquareMatrix<Var, double> facLU(const SquareMatrix<Var, double> &A, Vector<Var, int> &ipiv);

  /*! \brief LU decomposition
   *
   * This function computes the LU decomposition of a square matrix
   * according to \f[\boldsymbol{A}=\boldsymbol{P}\,\boldsymbol{L}\,\boldsymbol{U} \f]
   * \param A Pointer to first element of a square matrix.
   * \param ipiv A vector of integers containing the pivot indices.
   * \param info Information about success of routine
   * \return A square matrix containig the result.
   *
   * Copied and fitted on 28.01.2013 of http://www.mymathlib.com/c_source/matrices/linearsystems/doolittle_pivot.c
   * */
  FMATVEC_EXPORT int facLU(double *A, int pivot[], int n);
  
  /*! \brief LL decomposition
   *
   * This function computes the Cholesky decomposition of a symmetric matrix 
   * according to \f[\boldsymbol{A}=\boldsymbol{L}\,\boldsymbol{L} \f]
   * \param A A symmetric matrix. 
   * \return A symmetric matrix containig the result.
   * */
  FMATVEC_EXPORT Matrix<Symmetric, Ref, Ref, double> facLL(const Matrix<Symmetric, Ref, Ref, double> &A);

  /*! \brief 1-norm
   *
   * This function computes the sum of the absolute values of a vector.
   * \param x A vector.
   * \return A scalar containig the result.
   * */
  FMATVEC_EXPORT double nrm1(const Vector<Ref, double> &x);

  /*! \brief 2-norm
   *
   * This function computes the the Euclidean norm of a vector.
   * \param x A vector.
   * \return A scalar containig the result.
   * */
  FMATVEC_EXPORT double nrm2(const Vector<Ref, double> &x);

  /*! \brief Infinity-norm
   *
   * This function computes the the largest absolute value of a vector.
   * \param x A vector.
   * \return A scalar containig the result.
   * */
  FMATVEC_EXPORT double nrmInf(const Vector<Ref, double> &x);
  
  /*! \brief 1-norm
   *
   * This function computes the largest column sum of the absolute 
   * values of a general matrix.
   * \param A A general matrix.
   * \return A scalar containig the result.
   * */
  FMATVEC_EXPORT double nrm1(const Matrix<General, Ref, Ref, double> &A);

  /*! \brief 2-norm
   *
   * This function computes the largest singular value of a general matrix.
   * \param A A general matrix.
   * \return A scalar containig the result.
   * */
  FMATVEC_EXPORT double nrm2(const Matrix<General, Ref, Ref, double> &A);

  /*! \brief Infinity-norm
   *
   * This function computes the largest row sum of the absolute 
   * values of a general matrix.
   * \param A A general matrix.
   * \return A scalar containing the result.
   * */
  FMATVEC_EXPORT double nrmInf(const Matrix<General, Ref, Ref, double> &A);

  /*! \brief Infinity-norm
   *
   * This function computes the largest row sum of the absolute
   * values of a symmetric matrix.
   * \param A A general matrix.
   * \return A scalar containing the result.
   * */
  FMATVEC_EXPORT double nrmInf(const Matrix<Symmetric, Ref, Ref, double> &A);

  /*! \brief Frobenius-norm
   *
   * This function computes the Frobenius norm of a general matrix.
   * \param A A general matrix.
   * \return A scalar containig the result.
   * */
  FMATVEC_EXPORT double nrmFro(const Matrix<General, Ref, Ref, double> &A);

  /*! \brief Spectral radius
   *
   * This function computes the spectral radius of a square matrix.
   * \param A A square matrix.
   * \return A scalar containig the result.
   * */
  FMATVEC_EXPORT double rho(const SquareMatrix<Ref, double> &A);

  /*! \brief Spectral radius
   *
   * This function computes the spectral radius of a symmetric matrix.
   * \param A A symmetric matrix.
   * \return A scalar containig the result.
   * */
  FMATVEC_EXPORT double rho(const Matrix<Symmetric, Ref, Ref, double> &A);

  FMATVEC_EXPORT Vector<Ref, double> slvLLFac(const Matrix<Symmetric, Ref, Ref, double> &A, const Vector<Ref, double> &x);

  FMATVEC_EXPORT Matrix<General, Ref, Ref, double> slvLLFac(const Matrix<Symmetric, Ref, Ref, double> &A, const Matrix<General, Ref, Ref, double> &X);

  FMATVEC_EXPORT Matrix<General, Ref, Ref, double> slvLS(const Matrix<General, Ref, Ref, double> &A, const Matrix<General, Ref, Ref, double> &B, double rcond = -1);
  
  FMATVEC_EXPORT Vector<Ref, double> slvLS(const Matrix<General, Ref, Ref, double> &A, const Vector<Ref, double> &b, double rcond = -1);

//Matrix<General,Ref,Ref,double> slvLU(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_DIAG unit, const SquareMatrix<Ref,double> &A, const Matrix<General,Ref,Ref,double> &X, const Vector<Ref,int> &ipiv );

  /*! \brief Row interchanges
   *
   * This function performs a series of row interchanges on a general matrix.
   * \param A A general matrix.
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A Matrix containig the result.
   * */
//Matrix<General,Ref,Ref,double> swap(const Matrix<General,Ref,Ref,double> &A, const Vector<Ref,int> &ipiv );

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a LU decompostion.
   * \param A    A square matrix.
   * \param b    A vector containing the right hand side.
   * \param info Information about the success of the routine (0 = success)
   * \return A vector containig the solution.
   * */
  template <int size>
  Vector<Fixed<size>, double> slvLU(const SquareMatrix<Fixed<size>, double> &A, const Vector<Fixed<size>, double> &b, int & info) {

    if (size == 0)
      return b;

    int ipiv[size];

    SquareMatrix<Fixed<size>, double> LU = A;

    info = facLU(LU(), ipiv, size);

    Vector<Fixed<size>, double> x = b;
    Vector<Fixed<size>, double> y = b;

    info = Doolittle_LU_with_Pivoting_Solve(LU(), y(), ipiv, x(), size);

    return x;

  }

}

#endif

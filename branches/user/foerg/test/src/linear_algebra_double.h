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


  /*! \brief Eigenvalues
   *
   * This function computes the complex eigenvalues of a square matrix.
   * \param A A square matrix. 
   * \return A vector containig the eigenvalues.
   * */
  Vector<Ref,std::complex<double> > eigval(const SquareMatrix<Ref,double> &A); 

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
  int eigvec(const Matrix<Symmetric,Ref,Ref,double> &A, const Matrix<Symmetric,Ref,Ref,double> &B, SquareMatrix<Ref,double> &eigenvectors, Vector<Ref,double> &eigenvalues); 

  /*! \brief Eigenvalues
   *
   * This function computes the eigenvalues of a symmetric matrix.
   * \param A A symmetric matrix. 
   * \return A vector containig the eigenvalues.
   * */
  Vector<Ref,double> eigval(const Matrix<Symmetric,Ref,Ref,double> &A); 

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
  Vector<Ref,double> eigvalSel(const Matrix<Symmetric,Ref,Ref,double> &A, int il, int iu, double abstol=0);

  /*! \brief Systems of linear equations
   *
   * This function solves systems of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]
   * by a LU decompostion.
   * \param A A square matrix. 
   * \param B A general matrix containing the right hand sides.
   * \return A general matrix containig the solution.
   * */
  Matrix<General,Ref,Ref,double> slvLU(const SquareMatrix<Ref,double> &A, const Matrix<General,Ref,Ref,double> &B);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a LU decompostion.
   * \param A A square matrix. 
   * \param b A vector containing the right hand side.
   * \return A vector containig the solution.
   * */
  Vector<Ref,double> slvLU(const SquareMatrix<Ref,double> &A, const Vector<Ref,double> &b);

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
  Matrix<General,Ref,Ref,double> slvLUFac(const SquareMatrix<Ref,double> &A, const Matrix<General,Ref,Ref,double> &B, const Vector<Ref,int> &ipiv);

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
  Vector<Ref,double> slvLUFac(const SquareMatrix<Ref,double> &A, const Vector<Ref,double> &b, const Vector<Ref,int> &ipiv);

  /*! \brief Systems of linear equations
   *
   * This function solves systems of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]
   * by a LL decompostion.
   * \param A A symmetric matrix. 
   * \param B A general matrix containing the right hand sides.
   * \return A general matrix containig the solution.
   * */
  Matrix<General,Ref,Ref,double> slvLL(const Matrix<Symmetric,Ref,Ref,double> &A, const Matrix<General,Ref,Ref,double> &B);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a LL decompostion.
   * \param A A symmetric matrix. 
   * \param b A vector containing the right hand side.
   * \return A vector containig the solution.
   * */
  Vector<Ref,double> slvLL(const Matrix<Symmetric,Ref,Ref,double> &A, const Vector<Ref,double> &b);


  /*! \brief Systems of linear equations
   *
   * This function solves systems of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{X}=\boldsymbol{B} \f]
   * by a QR decompostion.
   * \param A A square matrix. 
   * \param B A general matrix containing the right hand sides.
   * \return A general matrix containig the solution.
   * */
  Matrix<General,Ref,Ref,double> slvQR(const SquareMatrix<Ref,double> &A, const Matrix<General,Ref,Ref,double> &B);

  /*! \brief System of linear equations
   *
   * This function solves a system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a QR decompostion.
   * \param A A symmetric matrix. 
   * \param b A vector containing the right hand side.
   * \return A vector containig the solution.
   * */
  Vector<Ref,double> slvQR(const SquareMatrix<Ref,double> &A, const Vector<Ref,double> &b);

  /*! \brief Inverse
   *
   * This function computes the inverse of a general matrix 
   * according to \f[\boldsymbol{A}^{-1} \f]
   * \param A A square matrix. 
   * \return A square matrix containig the result.
   * */
  SquareMatrix<Ref,double> inv(const SquareMatrix<Ref,double> &A);

  /*! \brief Inverse
   *
   * This function computes the inverse of a symmetric matrix 
   * according to \f[\boldsymbol{A}^{-1} \f]
   * \param A A symmetric matrix. 
   * \return A symmetric matrix containig the result.
   * */
  Matrix<Symmetric,Ref,Ref,double> inv(const Matrix<Symmetric,Ref,Ref,double> &A);

  /*! \brief Inverse
   *
   * This function computes the inverse of a diagonal matrix 
   * according to \f[\boldsymbol{A}^{-1} \f]
   * \param A A diagonal matrix. 
   * \return A diagonal matrix containig the result.
   * */
  Matrix<Diagonal,Ref,Ref,double> inv(const Matrix<Diagonal,Ref,Ref,double> &A);

  /*! \brief LU decomposition
   *
   * This function computes the LU decomposition of a general matrix 
   * according to \f[\boldsymbol{A}=\boldsymbol{P}\,\boldsymbol{L}\,\boldsymbol{U} \f]
   * \param A A general matrix. 
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A general matrix containig the result.
   * */
  Matrix<General,Ref,Ref,double> facLU(const Matrix<General,Ref,Ref,double> &A, Vector<Ref,int> &ipiv);

  /*! \brief LU decomposition
   *
   * This function computes the LU decomposition of a square matrix 
   * according to \f[\boldsymbol{A}=\boldsymbol{P}\,\boldsymbol{L}\,\boldsymbol{U} \f]
   * \param A A square matrix. 
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A square matrix containig the result.
   * */
  SquareMatrix<Ref,double> facLU(const SquareMatrix<Ref,double> &A, Vector<Ref,int> &ipiv);

  /*! \brief LL decomposition
   *
   * This function computes the Cholesky decomposition of a symmetric matrix 
   * according to \f[\boldsymbol{A}=\boldsymbol{L}\,\boldsymbol{L} \f]
   * \param A A symmetric matrix. 
   * \return A symmetric matrix containig the result.
   * */
  Matrix<Symmetric,Ref,Ref,double> facLL(const Matrix<Symmetric,Ref,Ref,double> &A);

  /*! \brief 1-norm
   *
   * This function computes the sum of the absolute values of a vector.
   * \param x A vector.
   * \return A scalar containig the result.
   * */
  double nrm1(const Vector<Ref,double> &x);

  /*! \brief 2-norm
   *
   * This function computes the the Euclidean norm of a vector.
   * \param x A vector.
   * \return A scalar containig the result.
   * */
  double nrm2(const Vector<Ref,double> &x);

  /*! \brief Infinity-norm
   *
   * This function computes the the largest absolute value of a vector.
   * \param x A vector.
   * \return A scalar containig the result.
   * */
  double nrmInf(const Vector<Ref,double> &x);

  /*! \brief 1-norm
   *
   * This function computes the largest column sum of the absolute 
   * values of a general matrix.
   * \param A A general matrix.
   * \return A scalar containig the result.
   * */
  double nrm1(const Matrix<General,Ref,Ref,double> &A);

  /*! \brief 2-norm
   *
   * This function computes the largest singular value of a general matrix.
   * \param A A general matrix.
   * \return A scalar containig the result.
   * */
  double nrm2(const Matrix<General,Ref,Ref,double> &A);

  /*! \brief Infinity-norm
   *
   * This function computes the largest row sum of the absolute 
   * values of a general matrix.
   * \param A A general matrix.
   * \return A scalar containig the result.
   * */
  double nrmInf(const Matrix<General,Ref,Ref,double> &A);

  /*! \brief Frobenius-norm
   *
   * This function computes the Frobenius norm of a general matrix.
   * \param A A general matrix.
   * \return A scalar containig the result.
   * */
  double nrmFro(const Matrix<General,Ref,Ref,double> &A);

  /*! \brief Spectral radius
   *
   * This function computes the spectral radius of a square matrix.
   * \param A A square matrix.
   * \return A scalar containig the result.
   * */
  double rho(const SquareMatrix<Ref,double> &A);

  /*! \brief Spectral radius
   *
   * This function computes the spectral radius of a symmetric matrix.
   * \param A A symmetric matrix.
   * \return A scalar containig the result.
   * */
  double rho(const Matrix<Symmetric,Ref,Ref,double> &A);
  
  Vector<Ref,double> slvLLFac(const Matrix<Symmetric,Ref,Ref,double> &A, const Vector<Ref,double> &x);

  Matrix<General,Ref,Ref,double> slvLLFac(const Matrix<Symmetric,Ref,Ref,double> &A, const Matrix<General,Ref,Ref,double> &X);

  Matrix<General,Ref,Ref,double> slvLS(const Matrix<General,Ref,Ref,double> &A, const Matrix<General,Ref,Ref,double> &B, double rcond=-1);

  Vector<Ref,double> slvLS(const Matrix<General,Ref,Ref,double> &A, const Vector<Ref,double> &b, double rcond=-1);

  //Matrix<General,Ref,Ref,double> slvLU(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_DIAG unit, const SquareMatrix<Ref,double> &A, const Matrix<General,Ref,Ref,double> &X, const Vector<Ref,int> &ipiv );

  /*! \brief Row interchanges
   *
   * This function performs a series of row interchanges on a general matrix.
   * \param A A general matrix. 
   * \param ipiv A vector of integers containing the pivot indices.
   * \return A Matrix containig the result.
   * */
  //Matrix<General,Ref,Ref,double> swap(const Matrix<General,Ref,Ref,double> &A, const Vector<Ref,int> &ipiv );
}

#endif

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
 *   martin.foerg@web.de
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
  Matrix<Diagonal, double >& operator*=(const Matrix<Diagonal, double >& A, double alpha);

  Matrix<Diagonal, double >& operator/=(const Matrix<Diagonal, double >& A, double alpha);

  SquareMatrix<double>& operator*=(const SquareMatrix<double>& A, double alpha);

  SquareMatrix<double>& operator/=(const SquareMatrix<double>& A, double alpha);

  SquareMatrix<double> operator*(const SquareMatrix<double> &A, double alpha);

  SquareMatrix<double> operator*(double alpha, const SquareMatrix<double> &A);

  SquareMatrix<double> operator/(const SquareMatrix<double> &A, double alpha);

  /*! \brief Matrix multiplying.
   *
   * Multiplies the calling matrix by the scalar a. This is equivalent to A<<A*a but may be slightly faster.
   * \param a The scalar the calling matrix will be multiplied by. 
   * \return The calling matrix.
   * */
  Matrix<General, double >& operator*=(const Matrix<General, double > &A, const double &a);

  /*! \brief Matrix division.
   *
   * Divides the calling matrix by the scalar a. This is equivalent to A<<A/a but may be slightly faster.
   * \param a The scalar the calling matrix will be divided by.
   * \return The calling matrix.
   * */
  Matrix<General, double >& operator/=(const Matrix<General, double > &A, const double &a);

  Matrix<General, double > operator+(const Matrix<General, double > &A, const Matrix<General, double > &B);

  SquareMatrix<double> operator+(const SquareMatrix<double> &A, const SquareMatrix<double> &B);

  Matrix<Symmetric, double > operator+(const Matrix<Symmetric, double > &A, const Matrix<Symmetric, double > &B);

  SquareMatrix<double> operator-(const SquareMatrix<double> &A, const SquareMatrix<double> &B);

  Matrix<General, double> operator*(const Matrix<General, double> &A, double alpha);

  Matrix<General, double> operator/(const Matrix<General, double> &A, double alpha);

  Matrix<General, double> operator*(double alpha, const Matrix<General, double> &A);

  Matrix<Symmetric, double> operator*(const Matrix<Symmetric, double> &A, double alpha);

  Matrix<Symmetric, double> operator/(const Matrix<Symmetric, double> &A, double alpha);

  Matrix<Symmetric, double> operator*(double alpha, const Matrix<Symmetric, double> &A);

  Matrix<Diagonal, double> operator*(const Matrix<Diagonal, double> &D, double alpha);

  Matrix<Diagonal, double> operator*(double alpha, const Matrix<Diagonal, double> &D);

  Matrix<Diagonal, double> operator/(const Matrix<Diagonal, double> &D, double alpha);

  //-------------------------------------
  // Matrix/Matrix operations
  //-------------------------------------

  Matrix<Symmetric, double> JTJ(const Matrix<General, double> &J);

  Matrix<Symmetric, double> JTMJ(const Matrix<Symmetric, double> &M, const Matrix<General, double> &J);

  Matrix<Symmetric, double> JTMJ(const Matrix<Diagonal, double> &M, const Matrix<General, double> &J);

  Matrix<Diagonal, double >& operator+=(const Matrix<Diagonal, double >& A, const Matrix<Diagonal, double > &A);

  Matrix<Symmetric, double>& operator+=(const Matrix<Symmetric, double >& A_, const Matrix<Symmetric, double> &B);

  Matrix<Diagonal, double >& operator-=(const Matrix<Diagonal, double >& A, const Matrix<Diagonal, double > &A);

  SquareMatrix<double>& operator+=(const SquareMatrix<double>& A, const SquareMatrix<double> &A);

  SquareMatrix<double>& operator-=(const SquareMatrix<double>& A, const SquareMatrix<double> &A);

  /*! \brief Matrix addition.
   *
   * Adds the matrix A on the calling matrix. This is equivalent to A<<A+B but may be slightly faster.
   * \param A The matrix, that will be added on the calling matrix. 
   * \return The calling matrix.
   * */
  Matrix<General, double >& operator+=(const Matrix<General, double > &A, const Matrix<General, double > &A);

  /*! \brief Matrix subtraction.
   *
   * Subtracts the matrix A from the calling matrix. This is equivalent to A<<A-B but may be slightly faster.
   * \param A The matrix, that will be subtracted from the calling matrix. 
   * \return The calling matrix.
   * */
  Matrix<General, double >& operator-=(const Matrix<General, double > &A, const Matrix<General, double > &A);

  Matrix<General, double> operator-(const Matrix<General, double> &A, const Matrix<General, double> &B);

  Matrix<General, double> operator*(const Matrix<General, double> &A, const Matrix<General, double> &B);

  SquareMatrix<double> operator*(const SquareMatrix<double> &A, const SquareMatrix<double> &B);

  Matrix<General, double> operator*(const Matrix<General, double> &A, const Matrix<Symmetric, double> &B);

  Matrix<General, double> operator*(const Matrix<Symmetric, double> &A, const Matrix<General, double> &B);

  Matrix<General, double> operator*(const Matrix<General, double> &A,const Matrix<Diagonal, double> &D);

  Matrix<General, double> operator*(const Matrix<Diagonal, double> &D, const Matrix<General, double> &B);

  //-------------------------------------
  // Vector operations
  //-------------------------------------

  /*! \brief Vector multiplying.
   *
   * Multiplies the calling vectortor by the scalar a. This is equivalent to x<<x*a but may be slightly faster.
   * \param a The scalar the calling vectortor will be multiplied by. 
   * \return The calling vectortor.
   * */
  Vector<double>& operator*=(const Vector<double> &x, const double& a);

  /*! \brief Vector division.
   *
   * Divides the calling vectortor by the scalar a. This is equivalent to A<<A/a but may be slightly faster.
   * \param a The scalar the calling vectortor will be divided by.
   * \return The calling vectortor.
   * */
  Vector<double>& operator/=(const Vector<double> &x, const double &a);

  /*! \brief Vector addition.
   *
   * Adds the vectortor x on the calling vectortor. This is equivalent to a<<a+x but may be slightly faster.
   * \param x The vectortor, that will be added on the calling vectortor. 
   * \return The calling vectortor.
   * */

  Vector<double> operator/(const Vector<double> &x, double alpha);

  Vector<double> operator*(double alpha, const Vector<double> &x);

  Vector<double> operator*(const Vector<double> &x, double alpha);

  //-------------------------------------
  // Vector/Vector operations
  //-------------------------------------

  Vector<double>& operator+=(const Vector<double> &x, const Vector<double> &x);

  /*! \brief Vector subtraction.
   *
   * Subtracts the vectortor x from the calling vectortor. This is equivalent to a<<a-x but may be slightly faster.
   * \param x The vectortor, that will be subtracted from the calling vectortor. 
   * \return The calling vectortor.
   * */
  Vector<double>& operator-=(const Vector<double> &x, const Vector<double> &x);

  Vector<double> operator+(const Vector<double> &x, const Vector<double> &y);

  Vector<double> operator-(const Vector<double> &x, const Vector<double> &y);


  //-------------------------------------
  // RowVector operations
  //-------------------------------------

  /*! \brief RowVector multiplying.
   *
   * Multiplies the calling vectortor by the scalar a. This is equivalent to x<<x*a but may be slightly faster.
   * \param a The scalar the calling vectortor will be multiplied by. 
   * \return The calling vectortor.
   * */
  RowVector<double>& operator*=(const RowVector<double> &x, const double& a);

  /*! \brief RowVector division.
   *
   * Divides the calling vectortor by the scalar a. This is equivalent to A<<A/a but may be slightly faster.
   * \param a The scalar the calling vectortor will be divided by.
   * \return The calling vectortor.
   * */
  RowVector<double>& operator/=(const RowVector<double> &x, const double &a);

  RowVector<double> operator*(double alpha, const RowVector<double> &x);

  RowVector<double> operator*(const RowVector<double> &x, double alpha);

  RowVector<double> operator/(const RowVector<double> &x, double alpha);

  //-------------------------------------
  // RowVector/RowVector operations
  //-------------------------------------

  /*! \brief RowVector addition.
   *
   * Adds the vectortor x on the calling vectortor. This is equivalent to a<<a+x but may be slightly faster.
   * \param x The vectortor, that will be added on the calling vectortor. 
   * \return The calling vectortor.
   * */
  RowVector<double>& operator+=(const RowVector<double> &x, const RowVector<double> &x);

  /*! \brief RowVector subtraction.
   *
   * Subtracts the vectortor x from the calling vectortor. This is equivalent to a<<a-x but may be slightly faster.
   * \param x The vectortor, that will be subtracted from the calling vectortor. 
   * \return The calling vectortor.
   * */
  RowVector<double>& operator-=(const RowVector<double> &x, const RowVector<double> &x);

  RowVector<double> operator+(const RowVector<double> &x, const RowVector<double> &y);

  RowVector<double> operator-(const RowVector<double> &x, const RowVector<double> &y);

  //-------------------------------------
  // Matrix/Vector operations
  //-------------------------------------

  Vector<double> operator*(const Matrix<General, double> &A, const Vector<double> &x);

  Vector<double> operator*(const Matrix<Symmetric, double> &A, const Vector<double> &x);

  Vector<double> operator*(const Matrix<Diagonal, double> &D, const Vector<double> &x);

  //-------------------------------------
  // RowVector/Matrix operations
  //-------------------------------------

  RowVector<double> operator*(const RowVector<double> &x, const Matrix<General, double> &A); 

  RowVector<double> operator*(const RowVector<double> &x, const Matrix<Diagonal, double> &D);

  RowVector<double> operator*(const RowVector<double> &x, const Matrix<Symmetric, double> &A);

  //-------------------------------------
  // RowVector/Vector operations
  //-------------------------------------

  double operator*(const RowVector<double> &x, const Vector<double> &y); 

  //-------------------------------------
  // operations
  //-------------------------------------
  Vector<complex<double> > eigval(const SquareMatrix<double> &A); 

  Vector<double> eigval(const Matrix<Symmetric, double> &A); 

  Matrix<General, double> slvLU(const SquareMatrix<double> &A, const Matrix<General, double> &X);

  Vector<double> slvLU(const SquareMatrix<double> &A, const Vector<double> &x);

  Matrix<General, double> slvLU(const SquareMatrix<double> &A, const Matrix<General, double> &X, const Vector<int> &ipiv);

  Vector<double> slvLU(const SquareMatrix<double> &A, const Vector<double> &x, const Vector<int> &ipiv);

  //Matrix<General, double> slvLU(CBLAS_SIDE side, CBLAS_UPLO uplo, const SquareMatrix<double> &A, const Matrix<General, double> &X, const Vector<int> &ipiv );
  Matrix<General, double> slvLU(CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_DIAG unit, const SquareMatrix<double> &A, const Matrix<General, double> &X, const Vector<int> &ipiv );

  Vector<double> slvLL(const Matrix<Symmetric, double> &A, const Vector<double> &x);

  Matrix<General, double> slvLL(const Matrix<Symmetric, double> &A, const Matrix<General, double> &X);
  
  Matrix<General, double> swap(const Matrix<General, double> &X, const Vector<int> &ipiv );

  Matrix<General, double> slvQR(const SquareMatrix<double> &A, const Matrix<General, double> &X);

  Vector<double> slvQR(const SquareMatrix<double> &A, const Vector<double> &x);

  SquareMatrix<double> inv(const SquareMatrix<double> &A);

  Matrix<Symmetric, double> inv(const Matrix<Symmetric, double> &A);

  Matrix<Diagonal, double> inv(const Matrix<Diagonal, double> &A);

  Matrix<General, double> lu(const Matrix<General, double> &A, Vector<int> &ipiv);
  SquareMatrix<double> lu(const SquareMatrix<double> &A, Vector<int> &ipiv);

  double nrm1(const Vector<double> &x);
  double nrm2(const Vector<double> &x);
  double nrmInf(const Vector<double> &x);

  double nrmInf(const Matrix<General,double> &A);
  double nrm1(const Matrix<General,double> &A);
  double nrm2(const Matrix<General,double> &A);
  double nrmFro(const Matrix<General,double> &A);

  double rho(const SquareMatrix<double> &A);
  double rho(const Matrix<Symmetric, double> &A);

  Vector<double> eseigval(const Matrix<Symmetric, double> &A, int il, int iu, double abstol=0);
  Vector<double> gelss(const Matrix<General,double> &A, const Vector<double> &b, double rcond=-1);

}

#endif

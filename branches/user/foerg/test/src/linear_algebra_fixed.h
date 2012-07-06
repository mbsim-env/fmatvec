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

#ifndef linear_algebra_fixed_h
#define linear_algebra_fixed_h

#include "fixed_square_matrix.h"
#include "fixed_vector.h"
#include "fixed_row_vector.h"
#include "fixed_symmetric_matrix.h"

namespace fmatvec {

  template <int M, int N, class AT>
    inline Matrix<General<Fixed<M,N> >, AT> operator+(const Matrix<General<Fixed<M,N> >, AT> &A1, const Matrix<General<Fixed<M,N> >, AT> &A2) {
      Matrix<General<Fixed<M,N> >, AT> A3(NONINIT);
      for(int i=0; i<M; i++) 
        for(int j=0; j<N; j++) 
          A3.e(i,j) = A1.e(i,j) + A2.e(i,j);
      return A3;
    }

  template <int M, class AT>
    inline Vector<General<Fixed<M,1> >, AT> operator+(const Vector<General<Fixed<M,1> >, AT> &x1, const Vector<General<Fixed<M,1> >, AT> &x2) {
      Vector<General<Fixed<M,1> >, AT> y(NONINIT);
      for(int i=0; i<M; i++) 
	y.e(i) = x1.e(i) + x2.e(i);
      return y;
    }

  template <int M, class AT>
    inline Vector<General<Fixed<M,1> >, AT> operator-(const Vector<General<Fixed<M,1> >, AT> &x1, const Vector<General<Fixed<M,1> >, AT> &x2) {
      Vector<General<Fixed<M,1> >, AT> y(NONINIT);
      for(int i=0; i<M; i++) 
	y.e(i) = x1.e(i) - x2.e(i);
      return y;
    }

  template <int M, class AT>
    inline RowVector<General<Fixed<1,M> >, AT> operator+(const RowVector<General<Fixed<1,M> >, AT> &x1, const RowVector<General<Fixed<1,M> >, AT> &x2) {
      RowVector<General<Fixed<1,M> >, AT> y(NONINIT);
      for(int i=0; i<M; i++) 
	y.e(i) = x1.e(i) + x2.e(i);
      return y;
    }




  template <int M, class AT>
    inline Vector<General<Fixed<M,1> >, AT> operator*(const AT &a, const Vector<General<Fixed<M,1> >, AT> &x) {
      Vector<General<Fixed<M,1> >, AT> y(NONINIT);
      for(int i=0; i<M; i++)
	y.e(i) = a*x.e(i);
      return y;
    }

  template <int M, class AT>
    inline Vector<General<Fixed<M,1> >, AT> operator*(const Vector<General<Fixed<M,1> >, AT> &x, const AT &a) {
      Vector<General<Fixed<M,1> >, AT> y(NONINIT);
      for(int i=0; i<M; i++)
	y.e(i) = a*x.e(i);
      return y;
    }

  template <int M, class AT>
    inline Vector<General<Fixed<M,1> >, AT> operator/(const Vector<General<Fixed<M,1> >, AT> &x, const AT &a) {
      Vector<General<Fixed<M,1> >, AT> y(NONINIT);
      for(int i=0; i<M; i++)
	y.e(i) = x.e(i)/a;
      return y;
    }

  template <int N, class AT>
    inline RowVector<General<Fixed<1,N> >, AT> operator*(const RowVector<General<Fixed<1,N> >, AT> &x, const AT &a) {
      RowVector<General<Fixed<1,N> >, AT> y(NONINIT);
      for(int i=0; i<N; i++)
	y.e(i) = a*x.e(i);
      return y;
    }

  template <int N, class AT>
    inline RowVector<General<Fixed<1,N> >, AT> operator/(const RowVector<General<Fixed<1,N> >, AT> &x, const AT &a) {
      RowVector<General<Fixed<1,N> >, AT> y(NONINIT);
      for(int i=0; i<N; i++)
	y.e(i) = x.e(i)/a;
      return y;
    }

  template <int M, class AT>
    inline Matrix<Symmetric<Fixed<M,M> >, AT> operator*(const AT &a, const Matrix<Symmetric<Fixed<M,M> >,AT> &A) {
      Matrix<Symmetric<Fixed<M,M> >, AT> B(NONINIT);
      for(int i=0; i<M; i++)
	for(int j=i; j<M; j++)
	  B.ej(i,j) = a*A.ej(i,j);
      return B;
    }

  template <int M, class AT>
    inline Matrix<Symmetric<Fixed<M,M> >, AT> operator*(const Matrix<Symmetric<Fixed<M,M> >,AT> &A, const AT &a) {
      Matrix<Symmetric<Fixed<M,M> >, AT> B(NONINIT);
      for(int i=0; i<M; i++)
	for(int j=i; j<M; j++)
	  B.ej(i,j) = a*A.ej(i,j);
      return B;
    }

  template <int M, int N, class AT>
    inline Matrix<General<Fixed<M,N> >, AT> trans(const Matrix<General<Fixed<M,N> >, AT> &A) {
      Matrix<General<Fixed<M,N> >, AT> B(NONINIT);
      for(int i=0; i<N; i++)
	for(int j=0; j<M; j++)
	  B.e(i,j) = A.e(j,i);
      return B;
    }

  template <int M, class AT>
    inline SquareMatrix<General<Fixed<M,M> >, AT> trans(const SquareMatrix<General<Fixed<M,M> >, AT> &A) {
      SquareMatrix<General<Fixed<M,M> >, AT> B(NONINIT);
      for(int i=0; i<M; i++)
	for(int j=0; j<M; j++)
	  B.e(i,j) = A.e(j,i);
      return B;
    }

  template <int M, class AT>
    inline Vector<General<Fixed<M,1> >, AT>& operator+=(const Vector<General<Fixed<M,1> >, AT> &x_, const Vector<General<Fixed<M,1> >, AT> &y) {

      Vector<General<Fixed<M,1> >, AT> &x = const_cast<Vector<General<Fixed<M,1> >, AT> &>(x_);
      for(int i=0; i<M; i++) 
	x.e(i) += y.e(i);
      return x;
    }

  template <int M, class AT>
    inline Vector<General<Fixed<M,1> >, AT>& operator-=(const Vector<General<Fixed<M,1> >, AT> &x_, const Vector<General<Fixed<M,1> >, AT> &y) {

      Vector<General<Fixed<M,1> >, AT> &x = const_cast<Vector<General<Fixed<M,1> >, AT> &>(x_);
      for(int i=0; i<M; i++) 
	x.e(i) -= y.e(i);
      return x;
    }

  template <int M, class AT>
    inline Vector<General<Fixed<M,1> >, AT> operator*=(const Vector<General<Fixed<M,1> >, AT> &x_, const AT &a) {
      Vector<General<Fixed<M,1> >, AT> &x = const_cast<Vector<General<Fixed<M,1> >, AT> &>(x_);
      for(int i=0; i<M; i++)
	x.e(i) *= a;
      return x;
    }

  template <int M, class AT>
    inline Vector<General<Fixed<M,1> >, AT> operator/=(const Vector<General<Fixed<M,1> >, AT> &x_, const AT &a) {
      Vector<General<Fixed<M,1> >, AT> &x = const_cast<Vector<General<Fixed<M,1> >, AT> &>(x_);
      for(int i=0; i<M; i++)
	x.e(i) /= a;
      return x;
    }

  template <class AT>
    inline Vector<General<Fixed<3,1> >, AT> crossProduct(const Vector<General<Fixed<3,1> >, AT> &x, const Vector<General<Fixed<3,1> >, AT> &y) {

      Vector<General<Fixed<3,1> >, AT> z(NONINIT);

      z.e(0) = x.e(1)*y.e(2) - x.e(2)*y.e(1);
      z.e(1) = x.e(2)*y.e(0) - x.e(0)*y.e(2);
      z.e(2) = x.e(0)*y.e(1) - x.e(1)*y.e(0);

      return z;
    }

  template <class AT>
    inline SquareMatrix<General<Fixed<3,3> >,AT> tilde(const Vector<General<Fixed<3,1> >, AT> &x) {

      SquareMatrix<General<Fixed<3,3> >,AT> B(NONINIT);

      B.e(0,0) =  0;
      B.e(1,1) =  0;
      B.e(2,2) =  0;
      B.e(0,1) = -x.e(2);
      B.e(0,2) =  x.e(1);
      B.e(1,0) =  x.e(2);
      B.e(1,2) = -x.e(0);
      B.e(2,0) = -x.e(1);
      B.e(2,1) =  x.e(0);

      return B;
    }

  template <int M, class AT>
    inline AT operator*(const RowVector<General<Fixed<1,M> >, AT> &x, const Vector<General<Fixed<M,1> >, AT> &y) {
      AT c = 0;
      for(int i=0; i<M; i++)
	c += x.e(i)*y.e(i);
      return c;
    }

  template <int M, class AT>
    inline AT nrm2(const Vector<General<Fixed<M,1> >, AT> &x) {
      AT c = 0;
      for(int i=0; i<M; i++)
	c += pow(x.e(i),2);
      return sqrt(c);
    }

  template <int N, class AT>
    inline AT nrm2(const RowVector<General<Fixed<1,N> >, AT> &x) {
      AT c = 0;
      for(int i=0; i<N; i++)
	c += pow(x.e(i),2);
      return sqrt(c);
    }

  template <int M, int N, class AT>
  inline Matrix<Symmetric<Fixed<N,N> >, AT> JTJ(const Matrix<General<Fixed<M,N> >, AT> &A) { 
    Matrix<Symmetric<Fixed<N,N> >, AT> S(NONINIT);
    for(int i=0; i<N; i++) {
	for(int k=i; k<N; k++) {
	  S.ej(i,k) = 0;
	  for(int j=0; j<M; j++) 
	    S.ej(i,k) += A.e(j,i)*A.e(j,k);
	}
      }
    return S;
  }

  template <int M, int N, class AT>
  inline Matrix<Symmetric<Fixed<N,N> >, AT> JTMJ(const Matrix<Symmetric<Fixed<M,M> >, AT> &B, const Matrix<General<Fixed<M,N> >, AT> &A) {

    Matrix<Symmetric<Fixed<N,N> >, AT> S(NONINIT);
    Matrix<General<Fixed<M,N> >, AT> C = B*A;

    for(int i=0; i<N; i++) {
	for(int k=i; k<N; k++) {
	  S.ej(i,k) = 0;
	  for(int j=0; j<M; j++) 
	    S.ej(i,k) += A.e(j,i)*C.e(j,k);
	}
      }
    return S;
  }

  template <int M, int N, class AT>
  inline Matrix<Symmetric<Fixed<M,M> >, AT> JMJT(const Matrix<General<Fixed<M,N> >, AT> &A, const Matrix<Symmetric<Fixed<N,N> >, AT> &B) {

    Matrix<Symmetric<Fixed<M,M> >, AT> S(NONINIT);
    Matrix<General<Fixed<M,N> >, AT> C = A*B;

    for(int i=0; i<M; i++) {
	for(int k=i; k<M; k++) {
	  S.ej(i,k) = 0;
	  for(int j=0; j<N; j++) 
	    S.ej(i,k) += C.e(k,j)*A.e(i,j);
	}
      }
    return S;
  }

  template <int M, class AT>
    inline Matrix<Symmetric<Fixed<M,M> >, AT> operator+(const Matrix<Symmetric<Fixed<M,M> >, AT> &A1, const Matrix<Symmetric<Fixed<M,M> >, AT> &A2) {
      Matrix<Symmetric<Fixed<M,M> >, AT> A3(NONINIT);
      for(int i=0; i<M; i++) 
	for(int j=i; j<M; j++) 
	  A3.ej(i,j) = A1.ej(i,j) + A2.ej(i,j);
      return A3;
    }

  template <int M, class AT>
    inline Matrix<Symmetric<Fixed<M,M> >, AT> operator-(const Matrix<Symmetric<Fixed<M,M> >, AT> &A1, const Matrix<Symmetric<Fixed<M,M> >, AT> &A2) {
      Matrix<Symmetric<Fixed<M,M> >, AT> A3(NONINIT);
      for(int i=0; i<M; i++) 
	for(int j=i; j<M; j++) 
	  A3.ej(i,j) = A1.ej(i,j) - A2.ej(i,j);
      return A3;
    }

  template <int M, int N, class AT>
    inline Matrix<General<Fixed<M,N> >, AT> operator-(const Matrix<General<Fixed<M,N> >, AT> &A) {

      Matrix<General<Fixed<M,N> >, AT> B(NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=0; j<A.cols(); j++)
          B.e(i,j)=-A.e(i,j);

      return B;
    }

  template <int M, class AT>
    inline Vector<General<Fixed<M,1> >, AT> operator-(const Vector<General<Fixed<M,1> >, AT> &x) {

      Vector<General<Fixed<M,1> >, AT> y(NONINIT);

      for(int i=0; i<x.size(); i++)
	y.e(i)=-x.e(i);

      return y;
    }

  template <int M, class AT>
    inline RowVector<General<Fixed<1,M> >, AT> operator-(const RowVector<General<Fixed<1,M> >, AT> &x) {

      RowVector<General<Fixed<1,M> >, AT> y(NONINIT);

      for(int i=0; i<x.size(); i++)
	y.e(i)=-x.e(i);

      return y;
    }

}

#endif

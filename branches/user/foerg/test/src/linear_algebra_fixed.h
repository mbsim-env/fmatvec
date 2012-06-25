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
    inline Matrix<GeneralFixed<M,N>, AT> operator+(const Matrix<GeneralFixed<M,N>, AT> &A1, const Matrix<GeneralFixed<M,N>, AT> &A2) {
      Matrix<GeneralFixed<M,N>, AT> A3(NONINIT);
      for(int i=0; i<M; i++) 
        for(int j=0; j<N; j++) 
          A3.e(i,j) = A1.e(i,j) + A2.e(i,j);
      return A3;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator+(const Vector<GeneralFixed<M,1>, AT> &x1, const Vector<GeneralFixed<M,1>, AT> &x2) {
      Vector<GeneralFixed<M,1>, AT> y(NONINIT);
      for(int i=0; i<M; i++) 
	y.e(i) = x1.e(i) + x2.e(i);
      return y;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator-(const Vector<GeneralFixed<M,1>, AT> &x1, const Vector<GeneralFixed<M,1>, AT> &x2) {
      Vector<GeneralFixed<M,1>, AT> y(NONINIT);
      for(int i=0; i<M; i++) 
	y.e(i) = x1.e(i) - x2.e(i);
      return y;
    }

  template <int M, class AT>
    inline RowVector<GeneralFixed<1,M>, AT> operator+(const RowVector<GeneralFixed<1,M>, AT> &x1, const RowVector<GeneralFixed<1,M>, AT> &x2) {
      RowVector<GeneralFixed<1,M>, AT> y(NONINIT);
      for(int i=0; i<M; i++) 
	y.e(i) = x1.e(i) + x2.e(i);
      return y;
    }

  template <int M, int N, int K, class AT>
    inline Matrix<GeneralFixed<M,N>, AT> operator*(const Matrix<GeneralFixed<M,N>, AT> &A1, const Matrix<GeneralFixed<N,K>, AT> &A2) {
      Matrix<GeneralFixed<M,K>, AT> A3(NONINIT);
      for(int i=0; i<M; i++) {
	for(int k=0; k<K; k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<N; j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }

  template <int M, class AT>
    inline SquareMatrix<GeneralFixed<M,M>, AT> operator*(const SquareMatrix<GeneralFixed<M,M>, AT> &A1, const SquareMatrix<GeneralFixed<M,M>, AT> &A2) {
      SquareMatrix<GeneralFixed<M,M>, AT> A3(NONINIT);
      for(int i=0; i<M; i++) {
        for(int k=0; k<M; k++) {
          A3.e(i,k) = 0;
          for(int j=0; j<M; j++) 
            A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
        }
      }
      return A3;
    }

  template <int M, int N, class AT>
    inline Matrix<GeneralFixed<M,N>, AT> operator*(const Matrix<SymmetricFixed<M>, AT> &A1, const Matrix<GeneralFixed<M,N>, AT> &A2) {

    Matrix<GeneralFixed<M,N>, AT> A3(NONINIT);

    for(int i=0; i<M; i++) {
	for(int k=0; k<N; k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<i; j++) 
	    A3.e(i,k) += A1.ei(i,j)*A2.e(j,k);
	  for(int j=i; j<M; j++) 
	    A3.e(i,k) += A1.ej(i,j)*A2.e(j,k);
	}
      }

    return A3;
    }

  template <int M, int N, class AT>
    inline Matrix<GeneralFixed<M,N>, AT> operator*( const Matrix<GeneralFixed<M,N>, AT> &A1, const Matrix<SymmetricFixed<N>, AT> &A2) {

    Matrix<GeneralFixed<M,N>, AT> A3(NONINIT);

    for(int i=0; i<M; i++) {
	for(int k=0; k<N; k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<k; j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.ej(j,k);
	  for(int j=k; j<N; j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.ei(j,k);
	}
      }

    return A3;
    }

  template <int M, int N, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator*(const Matrix<GeneralFixed<M,N>, AT> &A, const Vector<GeneralFixed<N,1>, AT> &x) {
      Vector<GeneralFixed<M,1>, AT> y(NONINIT);
      for(int i=0; i<M; i++) {
	y.e(i) = 0;
	for(int j=0; j<N; j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator*(const Matrix<SymmetricFixed<M>, AT> &A, const Vector<GeneralFixed<M,1>, AT> &x) {
      Vector<GeneralFixed<M,1>, AT> y(NONINIT);
      for(int i=0; i<M; i++) {
	y.e(i) = 0;
	for(int j=0; j<i; j++) 
	  y.e(i) += A.ei(i,j)*x.e(j);
	for(int j=i; j<M; j++) 
	  y.e(i) += A.ej(i,j)*x.e(j);
      }
      return y;
    }

  template <int M, int N, class AT>
    inline RowVector<GeneralFixed<1,N>, AT> operator*(const RowVector<GeneralFixed<1,M>, AT> &x, const Matrix<GeneralFixed<M,N>, AT> &A) {
      RowVector<GeneralFixed<1,N>, AT> y(NONINIT);
      for(int i=0; i<N; i++) {
	y.e(i) = 0;
	for(int j=0; j<M; j++) 
	  y.e(i) += x.e(j)*A.e(j,i);
      }
      return y;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator*(const AT &a, const Vector<GeneralFixed<M,1>, AT> &x) {
      Vector<GeneralFixed<M,1>, AT> y(NONINIT);
      for(int i=0; i<M; i++)
	y.e(i) = a*x.e(i);
      return y;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator*(const Vector<GeneralFixed<M,1>, AT> &x, const AT &a) {
      Vector<GeneralFixed<M,1>, AT> y(NONINIT);
      for(int i=0; i<M; i++)
	y.e(i) = a*x.e(i);
      return y;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator/(const Vector<GeneralFixed<M,1>, AT> &x, const AT &a) {
      Vector<GeneralFixed<M,1>, AT> y(NONINIT);
      for(int i=0; i<M; i++)
	y.e(i) = x.e(i)/a;
      return y;
    }

  template <int M, class AT>
    inline Matrix<SymmetricFixed<M>, AT> operator*(const AT &a, const Matrix<SymmetricFixed<M>,AT> &A) {
      Matrix<SymmetricFixed<M>, AT> B(NONINIT);
      for(int i=0; i<M; i++)
	for(int j=i; j<M; j++)
	  B.ej(i,j) = a*A.ej(i,j);
      return B;
    }

  template <int M, class AT>
    inline Matrix<SymmetricFixed<M>, AT> operator*(const Matrix<SymmetricFixed<M>,AT> &A, const AT &a) {
      Matrix<SymmetricFixed<M>, AT> B(NONINIT);
      for(int i=0; i<M; i++)
	for(int j=i; j<M; j++)
	  B.ej(i,j) = a*A.ej(i,j);
      return B;
    }

  template <int M, int N, class AT>
    inline Matrix<GeneralFixed<M,N>, AT> trans(const Matrix<GeneralFixed<M,N>, AT> &A) {
      Matrix<GeneralFixed<M,N>, AT> B(NONINIT);
      for(int i=0; i<N; i++)
	for(int j=0; j<M; j++)
	  B.e(i,j) = A.e(j,i);
      return B;
    }

  template <int M, class AT>
    inline SquareMatrix<GeneralFixed<M,M>, AT> trans(const SquareMatrix<GeneralFixed<M,M>, AT> &A) {
      SquareMatrix<GeneralFixed<M,M>, AT> B(NONINIT);
      for(int i=0; i<M; i++)
	for(int j=0; j<M; j++)
	  B.e(i,j) = A.e(j,i);
      return B;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT>& operator+=(const Vector<GeneralFixed<M,1>, AT> &x_, const Vector<GeneralFixed<M,1>, AT> &y) {

      Vector<GeneralFixed<M,1>, AT> &x = const_cast<Vector<GeneralFixed<M,1>, AT> &>(x_);
      for(int i=0; i<M; i++) 
	x.e(i) += y.e(i);
      return x;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT>& operator-=(const Vector<GeneralFixed<M,1>, AT> &x_, const Vector<GeneralFixed<M,1>, AT> &y) {

      Vector<GeneralFixed<M,1>, AT> &x = const_cast<Vector<GeneralFixed<M,1>, AT> &>(x_);
      for(int i=0; i<M; i++) 
	x.e(i) -= y.e(i);
      return x;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator*=(const Vector<GeneralFixed<M,1>, AT> &x_, const AT &a) {
      Vector<GeneralFixed<M,1>, AT> &x = const_cast<Vector<GeneralFixed<M,1>, AT> &>(x_);
      for(int i=0; i<M; i++)
	x.e(i) *= a;
      return x;
    }

  template <class AT>
    inline Vector<GeneralFixed<3,1>, AT> crossProduct(const Vector<GeneralFixed<3,1>, AT> &x, const Vector<GeneralFixed<3,1>, AT> &y) {

      Vector<GeneralFixed<3,1>, AT> z(NONINIT);

      z.e(0) = x.e(1)*y.e(2) - x.e(2)*y.e(1);
      z.e(1) = x.e(2)*y.e(0) - x.e(0)*y.e(2);
      z.e(2) = x.e(0)*y.e(1) - x.e(1)*y.e(0);

      return z;
    }

  template <class AT>
    inline SquareMatrix<GeneralFixed<3,3>,AT> tilde(const Vector<GeneralFixed<3,1>, AT> &x) {

      SquareMatrix<GeneralFixed<3,3>,AT> B(NONINIT);

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
    inline AT operator*(const RowVector<GeneralFixed<1,M>, AT> &x, const Vector<GeneralFixed<M,1>, AT> &y) {
      AT c = 0;
      for(int i=0; i<M; i++)
	c += x.e(i)*y.e(i);
      return c;
    }

  template <int M, class AT>
    inline AT nrm2(const Vector<GeneralFixed<M,1>, AT> &x) {
      AT c = 0;
      for(int i=0; i<M; i++)
	c += pow(x.e(i),2);
      return sqrt(c);
    }

  template <int M, int N, class AT>
  inline Matrix<SymmetricFixed<N>, AT> JTJ(const Matrix<GeneralFixed<M,N>, AT> &A) { 
    Matrix<SymmetricFixed<N>, AT> S(NONINIT);
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
  inline Matrix<SymmetricFixed<N>, AT> JTMJ(const Matrix<SymmetricFixed<M>, AT> &B, const Matrix<GeneralFixed<M,N>, AT> &A) {

    Matrix<SymmetricFixed<N>, AT> S(NONINIT);
    Matrix<GeneralFixed<M,N>, AT> C = B*A;

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
  inline Matrix<SymmetricFixed<M>, AT> JMJT(const Matrix<GeneralFixed<M,N>, AT> &A, const Matrix<SymmetricFixed<N>, AT> &B) {

    Matrix<SymmetricFixed<M>, AT> S(NONINIT);
    Matrix<GeneralFixed<M,N>, AT> C = A*B;

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
    inline Matrix<SymmetricFixed<M>, AT> operator+(const Matrix<SymmetricFixed<M>, AT> &A1, const Matrix<SymmetricFixed<M>, AT> &A2) {
      Matrix<SymmetricFixed<M>, AT> A3(NONINIT);
      for(int i=0; i<M; i++) 
	for(int j=i; j<M; j++) 
	  A3.ej(i,j) = A1.ej(i,j) + A2.ej(i,j);
      return A3;
    }

  template <int M, class AT>
    inline Matrix<SymmetricFixed<M>, AT> operator-(const Matrix<SymmetricFixed<M>, AT> &A1, const Matrix<SymmetricFixed<M>, AT> &A2) {
      Matrix<SymmetricFixed<M>, AT> A3(NONINIT);
      for(int i=0; i<M; i++) 
	for(int j=i; j<M; j++) 
	  A3.ej(i,j) = A1.ej(i,j) - A2.ej(i,j);
      return A3;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator-(const Vector<GeneralFixed<M,1>, AT> &x) {

      Vector<GeneralFixed<M,1>, AT> y(NONINIT);

      for(int i=0; i<x.size(); i++)
	y.e(i)=-x.e(i);

      return y;
    }
}

#endif

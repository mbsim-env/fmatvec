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
    inline Matrix<FixedGeneral<M,N>, AT> operator+(const Matrix<FixedGeneral<M,N>, AT> &A1, const Matrix<FixedGeneral<M,N>, AT> &A2) {
      Matrix<FixedGeneral<M,N>, AT> A3(NONINIT);
      for(int i=0; i<M*N; i++) 
	A3.e(i) = A1.e(i) + A2.e(i);
      return A3;
    }

  template <int M, class AT>
    inline FixedVector<M, AT> operator+(const FixedVector<M, AT> &x1, const FixedVector<M, AT> &x2) {
      FixedVector<M, AT> y(NONINIT);
      for(int i=0; i<M; i++) 
	y.e(i) = x1.e(i) + x2.e(i);
      return y;
    }

  template <int M, class AT>
    inline FixedVector<M, AT> operator-(const FixedVector<M, AT> &x1, const FixedVector<M, AT> &x2) {
      FixedVector<M, AT> y(NONINIT);
      for(int i=0; i<M; i++) 
	y.e(i) = x1.e(i) - x2.e(i);
      return y;
    }

  template <int M, int N, int K, class AT>
    inline Matrix<FixedGeneral<M,N>, AT> operator*(const Matrix<FixedGeneral<M,N>, AT> &A1, const Matrix<FixedGeneral<N,K>, AT> &A2) {
      Matrix<FixedGeneral<M,K>, AT> A3(NONINIT);
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
    inline FixedSquareMatrix<M,AT> operator*(const FixedSquareMatrix<M,AT> &A1, const FixedSquareMatrix<M,AT> &A2) {
      FixedSquareMatrix<M,AT> A3(NONINIT);
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
    inline Matrix<FixedGeneral<M,N>, AT> operator*(const Matrix<FixedSymmetric<M>, AT> &A1, const Matrix<FixedGeneral<M,N>, AT> &A2) {

    Matrix<FixedGeneral<M,N>, AT> A3(NONINIT);

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
    inline FixedVector<M, AT> operator*(const Matrix<FixedGeneral<M,N>, AT> &A, const FixedVector<N, AT> &x) {
      FixedVector<M, AT> y(NONINIT);
      for(int i=0; i<M; i++) {
	y.e(i) = 0;
	for(int j=0; j<N; j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }

  template <int M, class AT>
    inline FixedVector<M, AT> operator*(const Matrix<FixedSymmetric<M>, AT> &A, const FixedVector<M, AT> &x) {
      FixedVector<M, AT> y(NONINIT);
      for(int i=0; i<M; i++) {
	y.e(i) = 0;
	for(int j=0; j<i; j++) 
	  y.e(i) += A.ei(i,j)*x.e(j);
	for(int j=i; j<M; j++) 
	  y.e(i) += A.ej(i,j)*x.e(j);
      }
      return y;
    }

  template <int M, class AT>
    inline FixedVector<M, AT> operator*(AT a, const FixedVector<M, AT> &x) {
      FixedVector<M, AT> y(NONINIT);
      for(int i=0; i<M; i++)
	y.e(i) = a*x.e(i);
      return y;
    }

  template <int M, class AT>
    inline FixedVector<M, AT> operator*(const FixedVector<M, AT> &x, AT a) {
      FixedVector<M, AT> y(NONINIT);
      for(int i=0; i<M; i++)
	y.e(i) = a*x.e(i);
      return y;
    }

  template <int M, class AT>
    inline FixedVector<M, AT> operator/(const FixedVector<M, AT> &x, AT a) {
      FixedVector<M, AT> y(NONINIT);
      for(int i=0; i<M; i++)
	y.e(i) = x.e(i)/a;
      return y;
    }

  template <int M, class AT>
    inline Matrix<FixedSymmetric<M>, AT> operator*(AT a, const Matrix<FixedSymmetric<M>,AT> &A) {
      Matrix<FixedSymmetric<M>, AT> B(NONINIT);
      for(int i=0; i<M; i++)
	for(int j=i; j<M; j++)
	  B.ej(i,j) = a*A.ej(i,j);
      return B;
    }

  template <int M, class AT>
    inline Matrix<FixedSymmetric<M>, AT> operator*(const Matrix<FixedSymmetric<M>,AT> &A, AT a) {
      Matrix<FixedSymmetric<M>, AT> B(NONINIT);
      for(int i=0; i<M; i++)
	for(int j=i; j<M; j++)
	  B.ej(i,j) = a*A.ej(i,j);
      return B;
    }

  template <int M, int N, class AT>
    inline Matrix<FixedGeneral<M,N>, AT> trans(const Matrix<FixedGeneral<M,N>, AT> &A) {
      Matrix<FixedGeneral<M,N>, AT> B(NONINIT);
      for(int i=0; i<N; i++)
	for(int j=0; j<M; j++)
	  B.e(i,j) = A.e(j,i);
      return B;
    }

  template <int M, class AT>
    inline FixedSquareMatrix<M, AT> trans(const FixedSquareMatrix<M, AT> &A) {
      FixedSquareMatrix<M, AT> B(NONINIT);
      for(int i=0; i<M; i++)
	for(int j=0; j<M; j++)
	  B.e(i,j) = A.e(j,i);
      return B;
    }

  template <int M, class AT>
    inline FixedVector<M, AT>& operator+=(const FixedVector<M, AT> &x_, const FixedVector<M, AT> &y) {

      FixedVector<M, AT> &x = const_cast<FixedVector<M, AT> &>(x_);
      for(int i=0; i<M; i++) 
	x.e(i) += y.e(i);
      return x;
    }

  template <int M, class AT>
    inline FixedVector<M, AT> operator*=(const FixedVector<M, AT> &x_, AT a) {
      FixedVector<M, AT> &x = const_cast<FixedVector<M, AT> &>(x_);
      for(int i=0; i<M; i++)
	x.e(i) *= a;
      return x;
    }

  template <class AT>
    inline FixedVector<3,AT> crossProduct(const FixedVector<3,AT> &x, const FixedVector<3,AT> &y) {

      FixedVector<3,AT> z(NONINIT);

      z.e(0) = x.e(1)*y.e(2) - x.e(2)*y.e(1);
      z.e(1) = x.e(2)*y.e(0) - x.e(0)*y.e(2);
      z.e(2) = x.e(0)*y.e(1) - x.e(1)*y.e(0);

      return z;
    }

  template <class AT>
    inline FixedSquareMatrix<3,AT> tilde(const FixedVector<3,AT > &x) {

      FixedSquareMatrix<3,AT> B(NONINIT);

      B.e(0) =  0;
      B.e(4) =  0;
      B.e(8) =  0;
      B.e(3) = -x.e(2);
      B.e(6) =  x.e(1);
      B.e(1) =  x.e(2);
      B.e(7) = -x.e(0);
      B.e(2) = -x.e(1);
      B.e(5) =  x.e(0);

      return B;
    }

  template <int M, class AT>
    inline AT operator*(const FixedRowVector<M,AT> &x, const FixedVector<M,AT> &y) {
      AT c = 0;
      for(int i=0; i<M; i++)
	c += x.e(i)*y.e(i);
      return c;
    }

  template <int M, class AT>
    inline AT nrm2(const FixedVector<M,AT> &x) {
      AT c = 0;
      for(int i=0; i<M; i++)
	c += pow(x.e(i),2);
      return sqrt(c);
    }

  template <int M, int N, class AT>
  inline Matrix<FixedSymmetric<N>, AT> JTJ(const Matrix<FixedGeneral<M,N>, AT> &A) { 
    Matrix<FixedSymmetric<N>, AT> S(NONINIT);
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
  inline Matrix<FixedSymmetric<N>, AT> JTMJ(const Matrix<FixedSymmetric<M>, AT> &B, const Matrix<FixedGeneral<M,N>, AT> &A) {

    Matrix<FixedSymmetric<N>, AT> S(NONINIT);
    Matrix<FixedGeneral<M,N>, AT> C = B*A;

    for(int i=0; i<N; i++) {
	for(int k=i; k<N; k++) {
	  S.ej(i,k) = 0;
	  for(int j=0; j<M; j++) 
	    S.ej(i,k) += A.e(j,i)*C.e(j,k);
	}
      }
    return S;
  }

  template <int M, class AT>
    inline Matrix<FixedSymmetric<M>, AT> operator+(const Matrix<FixedSymmetric<M>, AT> &A1, const Matrix<FixedSymmetric<M>, AT> &A2) {
      Matrix<FixedSymmetric<M>, AT> A3(NONINIT);
      for(int i=0; i<M; i++) 
	for(int j=i; j<M; j++) 
	  A3.ej(i,j) = A1.ej(i,j) + A2.ej(i,j);
      return A3;
    }

  template <int M, class AT>
    inline Matrix<FixedSymmetric<M>, AT> operator-(const Matrix<FixedSymmetric<M>, AT> &A1, const Matrix<FixedSymmetric<M>, AT> &A2) {
      Matrix<FixedSymmetric<M>, AT> A3(NONINIT);
      for(int i=0; i<M; i++) 
	for(int j=i; j<M; j++) 
	  A3.ej(i,j) = A1.ej(i,j) - A2.ej(i,j);
      return A3;
    }

  template <int M, class AT>
    FixedVector<M,AT> operator-(const FixedVector<M,AT> &x) {

      FixedVector<M,AT> y(NONINIT);

      for(int i=0; i<x.size(); i++)
	y.e(i)=-x.e(i);

      return y;
    }
}

#endif

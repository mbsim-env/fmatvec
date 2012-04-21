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

#ifndef linear_algebra_fixed_var_h
#define linear_algebra_fixed_var_h

#include "fixed_var_general_matrix.h"
//#include "var_square_matrix.h"
#include "fixed_vector.h"
#include "var_vector.h"
//#include "var_row_vector.h"
#include "var_symmetric_matrix.h"

namespace fmatvec {

  template <int M, class AT>
    inline Matrix<GeneralFixedVar<M>, AT> operator+(const Matrix<GeneralFixedVar<M>, AT> &A1, const Matrix<GeneralFixedVar<M>, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.cols() == A2.cols());
#endif
      Matrix<GeneralFixedVar<M>, AT> A3(A1.cols(),NONINIT);
      for(int i=0; i<A3.rows()*A3.cols(); i++) 
	A3.e(i) = A1.e(i) + A2.e(i);
      return A3;
    }

  template <int M, class AT>
    inline Matrix<GeneralFixedVar<M>, AT> operator-(const Matrix<GeneralFixedVar<M>, AT> &A1, const Matrix<GeneralFixedVar<M>, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.cols() == A2.cols());
#endif
      Matrix<GeneralFixedVar<M>, AT> A3(A1.cols(),NONINIT);
      for(int i=0; i<A3.rows()*A3.cols(); i++) 
	A3.e(i) = A1.e(i) - A2.e(i);
      return A3;
    }

  template <int M, class AT>
    inline Matrix<GeneralFixedVar<M>, AT>& operator+=(const Matrix<GeneralFixedVar<M>, AT> &A1_, const Matrix<GeneralFixedVar<M>, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1_.cols() == A2.cols());
#endif
      Matrix<GeneralFixedVar<M>, AT> &A1 = const_cast<Matrix<GeneralFixedVar<M>, AT> &>(A1_);
      for(int i=0; i<A1.rows()*A1.cols(); i++) 
	A1.e(i) += A2.e(i);
      return A1;
    }

  template <int M, class AT>
    inline Matrix<GeneralFixedVar<M>, AT>& operator-=(const Matrix<GeneralFixedVar<M>, AT> &A1_, const Matrix<GeneralFixedVar<M>, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1_.cols() == A2.cols());
#endif
      Matrix<GeneralFixedVar<M>, AT> &A1 = const_cast<Matrix<GeneralFixedVar<M>, AT> &>(A1_);
      for(int i=0; i<A1.rows()*A1.cols(); i++) 
	A1.e(i) -= A2.e(i);
      return A1;
    }

  template <int M, int N, class AT>
    inline Matrix<GeneralFixedVar<M>, AT> operator*(const Matrix<GeneralFixedVar<M>, AT> &A1, const Matrix<GeneralFixedVar<N>, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.cols() == N);
#endif
      Matrix<GeneralFixedVar<M>, AT> A3(A2.cols(),NONINIT);
      for(int i=0; i<M; i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<N; j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }

  template <int M, class AT>
    inline Matrix<GeneralFixedVar<M>, AT> operator*(const Matrix<GeneralFixedVar<M>, AT> &A1, const Matrix<GeneralVar, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.cols() == A2.rows());
#endif
      Matrix<GeneralFixedVar<M>, AT> A3(A2.cols(),NONINIT);
      for(int i=0; i<M; i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<A2.rows(); j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }

  template <int M, int N, class AT>
    inline Matrix<GeneralFixedVar<M>, AT> operator*(const Matrix<GeneralFixed<M,N>, AT> &A1, const Matrix<GeneralFixedVar<N>, AT> &A2) {

      Matrix<GeneralFixedVar<M>, AT> A3(A2.cols(),NONINIT);
      for(int i=0; i<M; i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<N; j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }


  template <int M, class AT>
    inline Matrix<GeneralFixedVar<M>, AT> operator*(const Matrix<SymmetricFixed<M>, AT> &A1, const Matrix<GeneralFixedVar<M>, AT> &A2) {

    Matrix<GeneralFixedVar<M>, AT> A3(A2.cols(),NONINIT);

    for(int i=0; i<M; i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<i; j++) 
	    A3.e(i,k) += A1.ei(i,j)*A2.e(j,k);
	  for(int j=i; j<A1.cols(); j++) 
	    A3.e(i,k) += A1.ej(i,j)*A2.e(j,k);
	}
      }

    return A3;
    }

  template <int M, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator*(const Matrix<GeneralFixedVar<M>, AT> &A, const Vector<GeneralVar, AT> &x) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.cols() == x.size());
#endif
      Vector<GeneralFixed<M,1>, AT> y(NONINIT);
      for(int i=0; i<M; i++) {
	y.e(i) = 0;
	for(int j=0; j<A.cols(); j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }

  template <int M, int N, class AT>
    inline Vector<GeneralFixed<M,1>, AT> operator*(const Matrix<GeneralFixedVar<M>, AT> &A, const Vector<GeneralFixed<N,1>, AT> &x) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.cols() == N);
#endif
      Vector<GeneralFixed<M,1>, AT> y(NONINIT);
      for(int i=0; i<M; i++) {
	y.e(i) = 0;
	for(int j=0; j<N; j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }

  template <int M, class AT>
    inline Matrix<SymmetricVar, AT> JTJ(const Matrix<GeneralFixedVar<M>, AT> &A) { 
      Matrix<SymmetricVar, AT> S(A.cols(),NONINIT);
      for(int i=0; i<A.cols(); i++) {
	for(int k=i; k<A.cols(); k++) {
	  S.ej(i,k) = 0;
	  for(int j=0; j<M; j++) 
	    S.ej(i,k) += A.e(j,i)*A.e(j,k);
	}
      }
      return S;
    }

  template <int M, class AT>
  inline Matrix<SymmetricVar, AT> JTMJ(const Matrix<SymmetricFixed<M>, AT> &B, const Matrix<GeneralFixedVar<M>, AT> &A) {

    Matrix<SymmetricVar, AT> S(A.cols(),NONINIT);
    Matrix<GeneralFixedVar<M>, AT> C = B*A;

    for(int i=0; i<A.cols(); i++) {
	for(int k=i; k<A.cols(); k++) {
	  S.ej(i,k) = 0;
	  for(int j=0; j<M; j++) 
	    S.ej(i,k) += A.e(j,i)*C.e(j,k);
	}
      }
    return S;
  }

  template <int M, class AT>
    inline Matrix<GeneralFixedVar<M>, AT> operator-(const Matrix<GeneralFixedVar<M>, AT> &A) {

      Matrix<GeneralFixedVar<M>, AT> B(A.cols(),NONINIT);

      for(int i=0; i<M*A.cols(); i++)
	B.e(i)=-A.e(i);

      return B;
    }
}

#endif

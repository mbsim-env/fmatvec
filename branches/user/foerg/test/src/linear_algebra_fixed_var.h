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
//#include "var_symmetric_matrix.h"

namespace fmatvec {

  template <int M, class AT>
    inline Matrix<FixedVarGeneral<M>, AT> operator+(const Matrix<FixedVarGeneral<M>, AT> &A1, const Matrix<FixedVarGeneral<M>, AT> &A2) {
#ifdef FMATVEC_SIZE_CHECK
      assert(A1.cols() == A2.cols());
#endif
      Matrix<FixedVarGeneral<M>, AT> A3(A1.cols(),NONINIT);
      for(int i=0; i<A3.rows()*A3.cols(); i++) 
	A3.e(i) = A1.e(i) + A2.e(i);
      return A3;
    }

  template <int M, class AT>
    inline Matrix<FixedVarGeneral<M>, AT> operator-(const Matrix<FixedVarGeneral<M>, AT> &A1, const Matrix<FixedVarGeneral<M>, AT> &A2) {
#ifdef FMATVEC_SIZE_CHECK
      assert(A1.cols() == A2.cols());
#endif
      Matrix<FixedVarGeneral<M>, AT> A3(A1.cols(),NONINIT);
      for(int i=0; i<A3.rows()*A3.cols(); i++) 
	A3.e(i) = A1.e(i) - A2.e(i);
      return A3;
    }

  template <int M, int N, class AT>
    inline Matrix<FixedVarGeneral<M>, AT> operator*(const Matrix<FixedVarGeneral<M>, AT> &A1, const Matrix<FixedVarGeneral<N>, AT> &A2) {
#ifdef FMATVEC_SIZE_CHECK
      assert(A1.cols() == N);
#endif
      Matrix<FixedVarGeneral<M>, AT> A3(A2.cols(),NONINIT);
      for(int i=0; i<M; i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<N; j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }

  template <int M, int N, class AT>
    inline Matrix<FixedVarGeneral<M>, AT> operator*(const Matrix<FixedGeneral<M,N>, AT> &A1, const Matrix<FixedVarGeneral<N>, AT> &A2) {

      Matrix<FixedVarGeneral<M>, AT> A3(A2.cols(),NONINIT);
      for(int i=0; i<M; i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<N; j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }


//  template <class AT>
//    inline FixedSquareMatrix<M,AT> operator*(const FixedSquareMatrix<M,AT> &A1, const FixedSquareMatrix<M,AT> &A2) {
//      FixedSquareMatrix<M,AT> A3(NONINIT);
//      for(int i=0; i<A3.size(); i++) {
//        for(int k=0; k<A3.size(); k++) {
//          A3.e(i,k) = 0;
//          for(int j=0; j<A1.size(); j++) 
//            A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
//        }
//      }
//      return A3;
//    }

//  template <class AT>
//    inline Matrix<VarGeneral, AT> operator*(const Matrix<FixedSymmetric<M>, AT> &A1, const Matrix<VarGeneral, AT> &A2) {
//
//    Matrix<VarGeneral, AT> A3(NONINIT);
//
//    for(int i=0; i<A.rows(); i++) {
//	for(int k=0; k<A.cols(); k++) {
//	  A3.e(i,k) = 0;
//	  for(int j=0; j<i; j++) 
//	    A3.e(i,k) += A1.ei(i,j)*A2.e(j,k);
//	  for(int j=i; j<A1.cols(); j++) 
//	    A3.e(i,k) += A1.ej(i,j)*A2.e(j,k);
//	}
//      }
//
//    return A3;
//    }

  template <int M, class AT>
    inline FixedVector<M,AT> operator*(const Matrix<FixedVarGeneral<M>, AT> &A, const VarVector<AT> &x) {
#ifdef FMATVEC_SIZE_CHECK
      assert(A.cols() == x.size());
#endif
      FixedVector<M,AT> y(NONINIT);
      for(int i=0; i<M; i++) {
	y.e(i) = 0;
	for(int j=0; j<A.cols(); j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }

  template <int M, int N, class AT>
    inline FixedVector<M,AT> operator*(const Matrix<FixedVarGeneral<M>, AT> &A, const FixedVector<N,AT> &x) {
#ifdef FMATVEC_SIZE_CHECK
      assert(A.cols() == N);
#endif
      FixedVector<M,AT> y(NONINIT);
      for(int i=0; i<M; i++) {
	y.e(i) = 0;
	for(int j=0; j<N; j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }

//  template <class AT>
//    inline VarVector<AT> operator*(const Matrix<FixedSymmetric<M>, AT> &A, const VarVector<AT> &x) {
//      VarVector<AT> y(NONINIT);
//      for(int i=0; i<M; i++) {
//	y.e(i) = 0;
//	for(int j=0; j<i; j++) 
//	  y.e(i) += A.ei(i,j)*x.e(j);
//	for(int j=i; j<M; j++) 
//	  y.e(i) += A.ej(i,j)*x.e(j);
//      }
//      return y;
//    }
//
//  template <class AT>
//    inline VarVector<AT> operator*(double a, const VarVector<AT> &x) {
//      VarVector<AT> y(NONINIT);
//      for(int i=0; i<M; i++)
//	y.e(i) = a*x.e(i);
//      return y;
//    }
//
//  template <class AT>
//    inline VarVector<AT> operator*(const VarVector<AT> &x, double a) {
//      VarVector<AT> y(NONINIT);
//      for(int i=0; i<M; i++)
//	y.e(i) = a*x.e(i);
//      return y;
//    }
//
//  template <class AT>
//    inline VarVector<AT> operator/(const VarVector<AT> &x, double a) {
//      VarVector<AT> y(NONINIT);
//      for(int i=0; i<M; i++)
//	y.e(i) = x.e(i)/a;
//      return y;
//    }
//
//  template <class AT>
//    inline Matrix<VarGeneral, AT> trans(const Matrix<VarGeneral, AT> &A) {
//      Matrix<VarGeneral, AT> B(NONINIT);
//      for(int i=0; i<N; i++)
//	for(int j=0; j<M; j++)
//	  B.e(i,j) = A.e(j,i);
//      return B;
//    }
//
//  template <class AT>
//    inline FixedSquareMatrix<M, AT> trans(const FixedSquareMatrix<M, AT> &A) {
//      FixedSquareMatrix<M, AT> B(NONINIT);
//      for(int i=0; i<M; i++)
//	for(int j=0; j<M; j++)
//	  B.e(i,j) = A.e(j,i);
//      return B;
//    }
//
//  template <class AT>
//    inline VarVector<AT>& operator+=(const VarVector<AT> &x_, const VarVector<AT> &y) {
//
//      VarVector<AT> &x = const_cast<VarVector<AT> &>(x_);
//      for(int i=0; i<M; i++) 
//	x.e(i) += y.e(i);
//      return x;
//    }
//
//  template <class AT>
//    inline VarVector<AT> operator*=(const VarVector<AT> &x_, double a) {
//      VarVector<AT> &x = const_cast<VarVector<AT> &>(x_);
//      for(int i=0; i<M; i++)
//	x.e(i) *= a;
//      return x;
//    }
//
//  template <class AT>
//    inline VarVector<AT> crossProduct(const VarVector<AT> &x, const VarVector<AT> &y) {
//
//      VarVector<AT> z(NONINIT);
//
//      z.e(0) = x.e(1)*y.e(2) - x.e(2)*y.e(1);
//      z.e(1) = x.e(2)*y.e(0) - x.e(0)*y.e(2);
//      z.e(2) = x.e(0)*y.e(1) - x.e(1)*y.e(0);
//
//      return z;
//    }
//
//  template <class AT>
//    inline FixedSquareMatrix<3,AT> tilde(const VarVector<AT > &x) {
//
//      FixedSquareMatrix<3,AT> B(NONINIT);
//
//      B.e(0) =  0;
//      B.e(4) =  0;
//      B.e(8) =  0;
//      B.e(3) = -x.e(2);
//      B.e(6) =  x.e(1);
//      B.e(1) =  x.e(2);
//      B.e(7) = -x.e(0);
//      B.e(2) = -x.e(1);
//      B.e(5) =  x.e(0);
//
//      return B;
//    }
//
//  template <class AT>
//    inline AT operator*(const FixedRowVector<M,AT> &x, const VarVector<AT> &y) {
//      AT c = 0;
//      for(int i=0; i<M; i++)
//	c += x.e(i)*y.e(i);
//      return c;
//    }
//
//  template <class AT>
//    inline AT nrm2(const VarVector<AT> &x) {
//      AT c = 0;
//      for(int i=0; i<M; i++)
//	c += pow(x.e(i),2);
//      return sqrt(c);
//    }
//
  template <int M, class AT>
    inline Matrix<VarSymmetric, AT> JTJ(const Matrix<FixedVarGeneral<M>, AT> &A) { 
      Matrix<VarSymmetric, AT> S(A.cols(),NONINIT);
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
  inline Matrix<VarSymmetric, AT> JTMJ(const Matrix<FixedSymmetric<M>, AT> &B, const Matrix<FixedVarGeneral<M>, AT> &A) {

    Matrix<VarSymmetric, AT> S(A.cols(),NONINIT);
    Matrix<FixedVarGeneral<M>, AT> C = B*A;

    for(int i=0; i<A.cols(); i++) {
	for(int k=i; k<A.cols(); k++) {
	  S.ej(i,k) = 0;
	  for(int j=0; j<M; j++) 
	    S.ej(i,k) += A.e(j,i)*C.e(j,k);
	}
      }
    return S;
  }
//
//  template <class AT>
//    inline Matrix<FixedSymmetric<M>, AT> operator+(const Matrix<FixedSymmetric<M>, AT> &A1, const Matrix<FixedSymmetric<M>, AT> &A2) {
//      Matrix<FixedSymmetric<M>, AT> A3(NONINIT);
//      for(int i=0; i<M; i++) 
//	for(int j=i; j<M; j++) 
//	  A3.ej(i,j) = A1.ej(i,j) + A2.ej(i,j);
//      return A3;
//    }
//
//  template <class AT>
//    inline Matrix<FixedSymmetric<M>, AT> operator-(const Matrix<FixedSymmetric<M>, AT> &A1, const Matrix<FixedSymmetric<M>, AT> &A2) {
//      Matrix<FixedSymmetric<M>, AT> A3(NONINIT);
//      for(int i=0; i<M; i++) 
//	for(int j=i; j<M; j++) 
//	  A3.ej(i,j) = A1.ej(i,j) - A2.ej(i,j);
//      return A3;
//    }
//
//  template <class AT>
//    VarVector<AT> operator-(const VarVector<AT> &x) {
//
//      VarVector<AT> y(NONINIT);
//
//      for(int i=0; i<x.size(); i++)
//	y.e(i)=-x.e(i);
//
//      return y;
//    }
}

#endif

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

#ifndef linear_algebra_var_fixed_h
#define linear_algebra_var_fixed_h

#include "var_fixed_general_matrix.h"
#include "fixed_vector.h"
#include "var_vector.h"

namespace fmatvec {

  template <int N, class AT, class Type>
    inline Vector<GeneralVar, AT> operator*(const Matrix<GeneralVarFixed<N>, AT> &A, const Vector<Type, AT> &x) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.cols() == x.size());
#endif
      Vector<GeneralVar, AT> y(A.rows(),NONINIT);
      for(int i=0; i<A.rows(); i++) {
	y.e(i) = 0;
	for(int j=0; j<N; j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }

//  template <int N, class AT>
//    inline Vector<GeneralVar, AT> operator*(const Matrix<GeneralVarFixed<N>, AT> &A, const Vector<GeneralVar, AT> &x) {
//#ifndef FMATVEC_NO_SIZE_CHECK
//      assert(N == x.size());
//#endif
//      Vector<GeneralVar, AT> y(A.rows(),NONINIT);
//      for(int i=0; i<A.rows(); i++) {
//	y.e(i) = 0;
//	for(int j=0; j<N; j++) 
//	  y.e(i) += A.e(i,j)*x.e(j);
//      }
//      return y;
//    }
//
//  template <int N, class AT>
//    inline Vector<GeneralVar, AT> operator*(const Matrix<GeneralVarFixed<N>, AT> &A, const Vector<GeneralFixed<N,1>, AT> &x) {
//      Vector<GeneralVar, AT> y(A.rows(),NONINIT);
//      for(int i=0; i<A.rows(); i++) {
//	y.e(i) = 0;
//	for(int j=0; j<N; j++) 
//	  y.e(i) += A.e(i,j)*x.e(j);
//      }
//      return y;
//    }

  template <int N, class AT>
    inline Matrix<GeneralVar, AT> operator*(const Matrix<GeneralVarFixed<N>, AT> &A1, const Matrix<GeneralFixedVar<N>, AT> &A2) {
      Matrix<GeneralVar, AT> A3(A1.rows(), A2.cols(), NONINIT);
      for(int i=0; i<A3.rows(); i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<N; j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }

  template <int N, int K, class AT>
    inline Matrix<GeneralVarFixed<K>, AT> operator*(const Matrix<GeneralVarFixed<N>, AT> &A1, const Matrix<GeneralVarFixed<K>, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.cols() == A2.rows());
#endif
      Matrix<GeneralVarFixed<K>, AT> A3(A1.rows(), A2.cols(), NONINIT);
      for(int i=0; i<A3.rows(); i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<N; j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }

  template <int N, int K, class AT>
    inline Matrix<GeneralVarFixed<K>, AT> operator*(const Matrix<GeneralVarFixed<N>, AT> &A1, const Matrix<GeneralFixed<N,K>, AT> &A2) {
      Matrix<GeneralVarFixed<K>, AT> A3(A1.rows(), A2.cols(), NONINIT);
      for(int i=0; i<A3.rows(); i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<N; j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }

  template <int N, class AT, class Type>
    inline Matrix<GeneralVar, AT> operator*(const Matrix<GeneralVarFixed<N>, AT> &A1, const Matrix<Type, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.cols() == A2.rows());
#endif
      Matrix<GeneralVar, AT> A3(A1.rows(), A2.cols(), NONINIT);
      for(int i=0; i<A3.rows(); i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<N; j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }
}

#endif

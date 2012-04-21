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

#ifndef linear_algebra_var_h
#define linear_algebra_var_h

#include "var_general_matrix.h"
//#include "var_square_matrix.h"
#include "var_vector.h"
//#include "var_row_vector.h"
//#include "var_symmetric_matrix.h"

namespace fmatvec {

  template <class AT>
    inline Matrix<VarGeneral, AT> operator+(const Matrix<VarGeneral, AT> &A1, const Matrix<VarGeneral, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      Matrix<VarGeneral, AT> A3(A1.rows(),A1.cols(),NONINIT);
      for(int i=0; i<A3.rows()*A3.cols(); i++) 
	A3.e(i) = A1.e(i) + A2.e(i);
      return A3;
    }

  template <class AT>
    inline Matrix<VarSymmetric, AT> operator+(const Matrix<VarSymmetric, AT> &A1, const Matrix<VarSymmetric, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.size() == A2.size());
#endif
      Matrix<VarSymmetric, AT> A3(A1.size(),NONINIT);
      for(int i=0; i<A3.size(); i++) 
	for(int j=i; j<A3.size(); j++) 
	  A3.ej(i,j) = A1.ej(i,j) + A2.ej(i,j);
      return A3;
    }

  template <class AT>
    inline Vector<VarGeneral, AT> operator+(const Vector<VarGeneral, AT> &x1, const Vector<VarGeneral, AT> &x2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x1.size() == x2.size());
#endif
      Vector<VarGeneral, AT> y(x1.size(),NONINIT);
      for(int i=0; i<y.size(); i++) 
	y.e(i) = x1.e(i) + x2.e(i);
      return y;
    }

  template <class AT>
    inline Vector<VarGeneral, AT> operator-(const Vector<VarGeneral, AT> &x1, const Vector<VarGeneral, AT> &x2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x1.size() == x2.size());
#endif
      Vector<VarGeneral, AT> y(x1.size(),NONINIT);
      for(int i=0; i<y.size(); i++) 
	y.e(i) = x1.e(i) - x2.e(i);
      return y;
    }

  template <class AT>
    inline Matrix<VarGeneral, AT> operator*(const Matrix<VarGeneral, AT> &A1, const Matrix<VarGeneral, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.cols() == A2.rows());
#endif
      Matrix<VarGeneral, AT> A3(A1.rows(),A2.cols(),NONINIT);
      for(int i=0; i<A3.rows(); i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<A1.cols(); j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }

  template <class AT>
    inline Vector<VarGeneral, AT> operator*(const Matrix<VarGeneral, AT> &A, const Vector<VarGeneral, AT> &x) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.cols() == x.size());
#endif
      Vector<VarGeneral, AT> y(A.rows(),NONINIT);
      for(int i=0; i<y.size(); i++) {
	y.e(i) = 0;
	for(int j=0; j<A.cols(); j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }

  template <class AT>
    inline Vector<VarGeneral, AT> operator*(const AT &a, const Vector<VarGeneral, AT> &x) {
      Vector<VarGeneral, AT> y(x.size(),NONINIT);
      for(int i=0; i<x.size(); i++)
	y.e(i) = a*x.e(i);
      return y;
    }

  template <class AT>
    inline Vector<VarGeneral, AT> operator*(const Vector<VarGeneral, AT> &x, const AT &a) {
      Vector<VarGeneral, AT> y(x.size(),NONINIT);
      for(int i=0; i<x.size(); i++)
	y.e(i) = a*x.e(i);
      return y;
    }

  template <class AT>
    inline Matrix<VarSymmetric,AT> operator*(const AT &a, const Matrix<VarSymmetric,AT> &A) {
      Matrix<VarSymmetric,AT> B(A.size(),NONINIT);
      for(int i=0; i<A.size(); i++)
	for(int j=i; j<A.size(); j++)
	  B.ej(i,j) = a*A.ej(i,j);
      return B;
    }

  template <class AT>
    inline Matrix<VarSymmetric,AT> operator*(const Matrix<VarSymmetric,AT> &A, const AT &a) {
      Matrix<VarSymmetric,AT> B(A.size(),NONINIT);
      for(int i=0; i<A.size(); i++)
	for(int j=i; j<A.size(); j++)
	  B.ej(i,j) = a*A.ej(i,j);
      return B;
    }

}

#endif

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
    inline Matrix<GeneralVar, AT> operator+(const Matrix<GeneralVar, AT> &A1, const Matrix<GeneralVar, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      Matrix<GeneralVar, AT> A3(A1.rows(),A1.cols(),NONINIT);
      for(int i=0; i<A3.rows()*A3.cols(); i++) 
	A3.e(i) = A1.e(i) + A2.e(i);
      return A3;
    }

  template <class AT>
    inline Matrix<SymmetricVar, AT> operator+(const Matrix<SymmetricVar, AT> &A1, const Matrix<SymmetricVar, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.size() == A2.size());
#endif
      Matrix<SymmetricVar, AT> A3(A1.size(),NONINIT);
      for(int i=0; i<A3.size(); i++) 
	for(int j=i; j<A3.size(); j++) 
	  A3.ej(i,j) = A1.ej(i,j) + A2.ej(i,j);
      return A3;
    }

  template <class AT>
    inline Vector<GeneralVar, AT> operator+(const Vector<GeneralVar, AT> &x1, const Vector<GeneralVar, AT> &x2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x1.size() == x2.size());
#endif
      Vector<GeneralVar, AT> y(x1.size(),NONINIT);
      for(int i=0; i<y.size(); i++) 
	y.e(i) = x1.e(i) + x2.e(i);
      return y;
    }

  template <class AT>
    inline Vector<GeneralVar, AT> operator-(const Vector<GeneralVar, AT> &x1, const Vector<GeneralVar, AT> &x2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x1.size() == x2.size());
#endif
      Vector<GeneralVar, AT> y(x1.size(),NONINIT);
      for(int i=0; i<y.size(); i++) 
	y.e(i) = x1.e(i) - x2.e(i);
      return y;
    }

  template <class AT>
    inline Matrix<GeneralVar, AT> operator*(const Matrix<GeneralVar, AT> &A1, const Matrix<GeneralVar, AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.cols() == A2.rows());
#endif
      Matrix<GeneralVar, AT> A3(A1.rows(),A2.cols(),NONINIT);
      for(int i=0; i<A3.rows(); i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<A1.cols(); j++) 
	    A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
	}
      }
      return A3;
    }

  template <class AT, class Type>
    inline Vector<GeneralVar, AT> operator*(const Matrix<GeneralVar, AT> &A, const Vector<Type, AT> &x) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.cols() == x.size());
#endif
      Vector<GeneralVar, AT> y(A.rows(),NONINIT);
      for(int i=0; i<y.size(); i++) {
	y.e(i) = 0;
	for(int j=0; j<A.cols(); j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }

//  template <class AT>
//    inline Vector<GeneralVar, AT> operator*(const Matrix<GeneralVar, AT> &A, const Vector<GeneralVar, AT> &x) {
//#ifndef FMATVEC_NO_SIZE_CHECK
//      assert(A.cols() == x.size());
//#endif
//      Vector<GeneralVar, AT> y(A.rows(),NONINIT);
//      for(int i=0; i<y.size(); i++) {
//	y.e(i) = 0;
//	for(int j=0; j<A.cols(); j++) 
//	  y.e(i) += A.e(i,j)*x.e(j);
//      }
//      return y;
//    }

  template <class AT>
    inline Vector<GeneralVar, AT> operator*(const AT &a, const Vector<GeneralVar, AT> &x) {
      Vector<GeneralVar, AT> y(x.size(),NONINIT);
      for(int i=0; i<x.size(); i++)
	y.e(i) = a*x.e(i);
      return y;
    }

  template <class AT>
    inline Vector<GeneralVar, AT> operator*(const Vector<GeneralVar, AT> &x, const AT &a) {
      Vector<GeneralVar, AT> y(x.size(),NONINIT);
      for(int i=0; i<x.size(); i++)
	y.e(i) = a*x.e(i);
      return y;
    }

  template <class AT>
    inline Matrix<SymmetricVar,AT> operator*(const AT &a, const Matrix<SymmetricVar,AT> &A) {
      Matrix<SymmetricVar,AT> B(A.size(),NONINIT);
      for(int i=0; i<A.size(); i++)
	for(int j=i; j<A.size(); j++)
	  B.ej(i,j) = a*A.ej(i,j);
      return B;
    }

  template <class AT>
    inline Matrix<SymmetricVar,AT> operator*(const Matrix<SymmetricVar,AT> &A, const AT &a) {
      Matrix<SymmetricVar,AT> B(A.size(),NONINIT);
      for(int i=0; i<A.size(); i++)
	for(int j=i; j<A.size(); j++)
	  B.ej(i,j) = a*A.ej(i,j);
      return B;
    }

  template <class AT>
    inline Matrix<GeneralVar, AT> operator*(const Matrix<GeneralVar, AT> &A, const AT &a) {
      Matrix<GeneralVar, AT> B(A.rows(),A.cols(),NONINIT);

      for(int i=0; i<A.rows(); i++)
	for(int j=0; j<A.cols(); j++)
	  B.e(i,j) = a*A.e(i,j);
      return B;
    }

  /*! \brief Negation.
   *
   * This function computes the negation of a matrix
   * \return The negation. 
   * */
  template <class AT>
    Matrix<GeneralVar, AT> operator-(const Matrix<GeneralVar, AT> &A) {

      Matrix<GeneralVar, AT> C(A.rows(),A.cols(),NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=0; j<A.cols(); j++)
          C.e(i,j)=-A.e(i,j);

      return C;
    }

  /*! \brief Negation.
   *
   * This function computes the negation of a matrix
   * \return The negation. 
   * */
  template <class AT>
    Vector<GeneralVar, AT> operator-(const Vector<GeneralVar, AT> &a) {

      Vector<GeneralVar, AT> c(a.size(),NONINIT);

      for(int i=0; i<a.size(); i++)
          c.e(i)=-a.e(i);

      return c;
    }
}

#endif

/* Copyright (C) 2003-2005  Martin Förg, Rober Huber

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

#ifndef linear_algebra_h
#define linear_algebra_h

#include "square_matrix.h"
#include "fixed_vector.h"
namespace fmatvec {

/////////////////////////////////// vecvecadd //////////////////////////////

  /*! \brief Vector-vector addition.
   *
   * This function computes the sum of two vectors
   * \return The sum.
   * */
  template <class AT>
    inline Vector<General<Ref,Fixed<1> >, AT> operator+(const Vector<General<Ref,Fixed<1> >, AT> &x, const Vector<General<Ref,Fixed<1> >, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      Vector<General<Ref,Fixed<1> >, AT> z(x.size(),NONINIT);

      if(y.transposed()) {
        if(x.transposed())
          for(int i=0; i<x.size(); i++)
            z.er(i) = x.et(i) + y.et(i);
        else
          for(int i=0; i<x.size(); i++)
            z.er(i) = x.er(i) + y.et(i);
      }
      else {
        if(x.transposed())
          for(int i=0; i<x.size(); i++)
            z.er(i) = x.et(i) + y.er(i);
        else 
          for(int i=0; i<x.size(); i++)
            z.er(i) = x.er(i) + y.er(i);
      } 

      return z;
    }

  /*! \brief Vector-vector addition.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, int M, class Type>
    Vector<General<Fixed<M>,Fixed<1> >, AT> operator+(const Vector<General<Fixed<M>,Fixed<1> >, AT> &x, const Vector<Type, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      Vector<General<Fixed<M>,Fixed<1> >, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) + y.e(i);

      return z;
    }

  /*! \brief Vector-vector addition.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, class Type, int M> 
    Vector<General<Fixed<M>,Fixed<1> >, AT> operator+(const Vector<Type, AT> &x, const Vector<General<Fixed<M>,Fixed<1> >, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      Vector<General<Fixed<M>,Fixed<1> >, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) + y.e(i);

      return z;
    }


  /*! \brief Vector-vector addition.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, class Type1, class Type2>
    Vector<General<Var,Fixed<1> >, AT> operator+(const Vector<Type1, AT> &x, const Vector<Type2, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      Vector<General<Var,Fixed<1> >, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) + y.e(i);

      return z;
    }

  /*! \brief Vector-vector subtraction.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, int M, class Type>
    Vector<General<Fixed<M>,Fixed<1> >, AT> operator-(const Vector<General<Fixed<M>,Fixed<1> >, AT> &x, const Vector<Type, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      Vector<General<Fixed<M>,Fixed<1> >, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) - y.e(i);

      return z;
    }

  /*! \brief Vector-vector subtraction.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, class Type, int M> 
    Vector<General<Fixed<M>,Fixed<1> >, AT> operator-(const Vector<Type, AT> &x, const Vector<General<Fixed<M>,Fixed<1> >, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      Vector<General<Fixed<M>,Fixed<1> >, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) - y.e(i);

      return z;
    }

  /*! \brief Vector-vector subtraction.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, class Type1, class Type2>
    Vector<General<Var,Fixed<1> >, AT> operator-(const Vector<Type1, AT> &x, const Vector<Type2, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      Vector<General<Var,Fixed<1> >, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) - y.e(i);

      return z;
    }

  /*! \brief Vector-vector addition.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, int M, class Type>
    RowVector<General<Fixed<1>,Fixed<M> >, AT> operator+(const RowVector<General<Fixed<1>,Fixed<M> >, AT> &x, const Vector<Type, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      RowVector<General<Fixed<1>,Fixed<M> >, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) + y.e(i);

      return z;
    }

  /*! \brief Vector-vector addition.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, class Type, int M> 
    RowVector<General<Fixed<1>,Fixed<M> >, AT> operator+(const Vector<Type, AT> &x, const RowVector<General<Fixed<1>,Fixed<M> >, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      RowVector<General<Fixed<1>,Fixed<M> >, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) + y.e(i);

      return z;
    }


  /*! \brief Vector-vector addition.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, class Type1, class Type2>
    RowVector<General<Fixed<1>,Var>, AT> operator+(const Vector<Type1, AT> &x, const Vector<Type2, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      RowVector<General<Fixed<1>,Var>, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) + y.e(i);

      return z;
    }

  /*! \brief Vector-vector subtraction.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, int M, class Type>
    RowVector<General<Fixed<1>,Fixed<M> >, AT> operator-(const RowVector<General<Fixed<1>,Fixed<M> >, AT> &x, const Vector<Type, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      RowVector<General<Fixed<1>,Fixed<M> >, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) - y.e(i);

      return z;
    }

  /*! \brief Vector-vector subtraction.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, class Type, int M> 
    RowVector<General<Fixed<1>,Fixed<M> >, AT> operator-(const Vector<Type, AT> &x, const RowVector<General<Fixed<1>,Fixed<M> >, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      RowVector<General<Fixed<1>,Fixed<M> >, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) - y.e(i);

      return z;
    }

  /*! \brief Vector-vector subtraction.
   *
   * This function computes the difference of two vectors
   * \return The difference.
   * */
  template <class AT, class Type1, class Type2>
    RowVector<General<Fixed<1>,Var>, AT> operator-(const Vector<Type1, AT> &x, const Vector<Type2, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      RowVector<General<Fixed<1>,Var>, AT> z(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        z.e(i) = x.e(i) - y.e(i);

      return z;
    }

  template <class Type1, class Type2, class AT>
    inline Vector<Type1, AT>& operator+=(const Vector<Type1, AT> &x_, const Vector<Type2, AT> &y) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x_.size() == y.size());
#endif
      Vector<Type1, AT> &x = const_cast<Vector<Type1, AT> &>(x_);
      for(int i=0; i<x.size(); i++) 
        x.e(i) += y.e(i);
      return x;
    }

  template <class Type1, class Type2, class AT>
    inline Vector<Type1, AT>& operator-=(const Vector<Type1, AT> &x_, const Vector<Type2, AT> &y) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x_.size() == y.size());
#endif
      Vector<Type1, AT> &x = const_cast<Vector<Type1, AT> &>(x_);
      for(int i=0; i<x.size(); i++) 
        x.e(i) -= y.e(i);
      return x;
    }

  template <class Type1, class Type2, class AT>
    inline RowVector<Type1, AT>& operator+=(const RowVector<Type1, AT> &x_, const RowVector<Type2, AT> &y) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x_.size() == y.size());
#endif
      RowVector<Type1, AT> &x = const_cast<RowVector<Type1, AT> &>(x_);
      for(int i=0; i<x.size(); i++) 
        x.e(i) += y.e(i);
      return x;
    }

  template <class Type1, class Type2, class AT>
    inline RowVector<Type1, AT>& operator-=(const RowVector<Type1, AT> &x_, const RowVector<Type2, AT> &y) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x_.size() == y.size());
#endif
      RowVector<Type1, AT> &x = const_cast<RowVector<Type1, AT> &>(x_);
      for(int i=0; i<x.size(); i++) 
        x.e(i) -= y.e(i);
      return x;
    }

/////////////////////////////////// end vecvecadd //////////////////////////////
  template <class Type1, class Type2, class Type3, class AT>
    inline void add(const Matrix<Type1, AT> &A1, const Matrix<Type2, AT> &A2, Matrix<Type3, AT> &A3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      for(int i=0; i<A1.rows(); i++)
        for(int j=0; j<A2.cols(); j++)
          A3.e(i,j) = A1.e(i,j) + A2.e(i,j);
    }

  template <class Type1, class Type2, class Type3, class AT>
    inline void sub(const Matrix<Type1, AT> &A1, const Matrix<Type2, AT> &A2, Matrix<Type3, AT> &A3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      for(int i=0; i<A1.rows(); i++)
        for(int j=0; j<A1.cols(); j++)
          A3.e(i,j) = A1.e(i,j) - A2.e(i,j);
    }

 /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, int M, int N, class Type> 
    Matrix<General<Fixed<M>,Fixed<N> >, AT> operator+(const Matrix<General<Fixed<M>,Fixed<N> >, AT> &A, const Matrix<Type, AT > &B) {
      Matrix<General<Fixed<M>,Fixed<N> >, AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, int M, int N, class Type> 
    Matrix<General<Fixed<M>,Fixed<N> >, AT> operator-(const Matrix<General<Fixed<M>,Fixed<N> >, AT> &A, const Matrix<Type, AT > &B) {
      Matrix<General<Fixed<M>,Fixed<N> >, AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, class Type, int M, int N>
    Matrix<General<Fixed<M>,Fixed<N> >, AT> operator+(const Matrix<Type, AT > &A, const Matrix<General<Fixed<M>,Fixed<N> >, AT> &B) {
      Matrix<General<Fixed<M>,Fixed<N> >, AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, class Type, int M, int N>
    Matrix<General<Fixed<M>,Fixed<N> >, AT> operator-(const Matrix<Type, AT > &A, const Matrix<General<Fixed<M>,Fixed<N> >, AT> &B) {
      Matrix<General<Fixed<M>,Fixed<N> >, AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, class Row, int N, class Type> 
    Matrix<General<Row,Fixed<N> >, AT> operator+(const Matrix<General<Row,Fixed<N> >, AT> &A, const Matrix<Type, AT > &B) {
      Matrix<General<Row,Fixed<N> >, AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, class Row, int N, class Type> 
    Matrix<General<Row,Fixed<N> >, AT> operator-(const Matrix<General<Row,Fixed<N> >, AT> &A, const Matrix<Type, AT > &B) {
      Matrix<General<Row,Fixed<N> >, AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, class Type, class Row, int N>
    Matrix<General<Row,Fixed<N> >, AT> operator+(const Matrix<Type, AT > &A, const Matrix<General<Row,Fixed<N> >, AT> &B) {
      Matrix<General<Row,Fixed<N> >, AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, class Type, class Row, int N>
    Matrix<General<Row,Fixed<N> >, AT> operator-(const Matrix<Type, AT > &A, const Matrix<General<Row,Fixed<N> >, AT> &B) {
      Matrix<General<Row,Fixed<N> >, AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, int M, class Col, class Type> 
    Matrix<General<Fixed<M>,Col>, AT> operator+(const Matrix<General<Fixed<M>,Col>, AT> &A, const Matrix<Type, AT > &B) {
      Matrix<General<Fixed<M>,Col>, AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, int M, class Col, class Type> 
    Matrix<General<Fixed<M>,Col>, AT> operator-(const Matrix<General<Fixed<M>,Col>, AT> &A, const Matrix<Type, AT > &B) {
      Matrix<General<Fixed<M>,Col>, AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, class Type, int M, class Col>
    Matrix<General<Fixed<M>,Col>, AT> operator+(const Matrix<Type, AT > &A, const Matrix<General<Fixed<M>,Col>, AT> &B) {
      Matrix<General<Fixed<M>,Col>, AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, class Type, int M, class Col>
    Matrix<General<Fixed<M>,Col>, AT> operator-(const Matrix<Type, AT > &A, const Matrix<General<Fixed<M>,Col>, AT> &B) {
      Matrix<General<Fixed<M>,Col>, AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <template <class Row, class Col> class Type1, template <class Row, class Col> class Type2, class AT> 
    Matrix<General<Ref,Ref>, AT> operator+(const Matrix<Type1<Ref,Ref>, AT> &A, const Matrix<Type2<Ref,Ref>, AT > &B) {
      Matrix<General<Ref,Ref>, AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <template <class Row, class Col> class Type1, template <class Row, class Col> class Type2, class AT> 
    Matrix<General<Ref,Ref>, AT> operator-(const Matrix<Type1<Ref,Ref>, AT> &A, const Matrix<Type2<Ref,Ref>, AT > &B) {
      Matrix<General<Ref,Ref>, AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, class Type1, class Type2> 
    Matrix<General<Var,Var>, AT> operator+(const Matrix<Type1, AT> &A, const Matrix<Type2, AT > &B) {
      Matrix<General<Var,Var>, AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }

  /*! \brief Matrix-matrix addition.
   *
   * This function computes the sum of two matrices
   * \return The sum.
   * */
  template <class AT, class Type1, class Type2> 
    Matrix<General<Var,Var>, AT> operator-(const Matrix<Type1, AT> &A, const Matrix<Type2, AT > &B) {
      Matrix<General<Var,Var>, AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }

  template <class Type1, class Type2, class AT>
    inline Matrix<Type1, AT>& operator+=(const Matrix<Type1, AT> &A_, const Matrix<Type2, AT> &B) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A_.rows() == B.rows());
      assert(A_.cols() == B.cols());
#endif
      Matrix<Type1, AT> &A = const_cast<Matrix<Type1, AT> &>(A_);
      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++) 
          A.e(i,j) += B.e(i,j);
      return A;
    }

  template <class Type1, class Type2, class AT>
    inline Matrix<Type1, AT>& operator-=(const Matrix<Type1, AT> &A_, const Matrix<Type2, AT> &B) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A_.rows() == B.rows());
      assert(A_.cols() == B.cols());
#endif
      Matrix<Type1, AT> &A = const_cast<Matrix<Type1, AT> &>(A_);
      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++) 
          A.e(i,j) -= B.e(i,j);
      return A;
    }

  template <class Row1, class Row2, class AT>
    inline Matrix<Symmetric<Row1,Row1>, AT>& operator+=(const Matrix<Symmetric<Row1,Row1>, AT> &A_, const Matrix<Symmetric<Row2,Row2>, AT> &B) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A_.size() == B.size());
#endif
      Matrix<Symmetric<Row1,Row1>, AT> &A = const_cast<Matrix<Symmetric<Row1,Row1>, AT> &>(A_);
      for(int i=0; i<A.size(); i++) 
        for(int j=i; j<A.size(); j++) 
          A.ej(i,j) += B.ej(i,j);
      return A;
    }

  template <class Row1, class Row2, class AT>
    inline Matrix<Symmetric<Row1,Row1>, AT>& operator-=(const Matrix<Symmetric<Row1,Row1>, AT> &A_, const Matrix<Symmetric<Row2,Row2>, AT> &B) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A_.size() == B.size());
#endif
      Matrix<Symmetric<Row1,Row1>, AT> &A = const_cast<Matrix<Symmetric<Row1,Row1>, AT> &>(A_);
      for(int i=0; i<A.size(); i++) 
        for(int j=i; j<A.size(); j++) 
          A.ej(i,j) -= B.ej(i,j);
      return A;
    }

/////////////////////////////////// end vecvecadd //////////////////////////////

/////////////////////////////////// matmatmult //////////////////////////////

  template <class Type1, class Type2, class Type3, class AT>
    inline void mult(const Matrix<Type1, AT> &A1, const Matrix<Type2, AT> &A2, Matrix<Type3, AT> &A3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.cols() == A2.rows());
#endif
      for(int i=0; i<A3.rows(); i++) {
        for(int k=0; k<A3.cols(); k++) {
          A3.e(i,k) = 0;
          for(int j=0; j<A1.cols(); j++) 
            A3.e(i,k) += A1.e(i,j)*A2.e(j,k);
        }
      }
    }

  template <class Row, class Type2, class Type3, class AT>
    inline void mult(const Matrix<Symmetric<Row,Row>, AT> &A1, const Matrix<Type2, AT> &A2, Matrix<Type3, AT> &A3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.cols() == A2.rows());
#endif
    for(int i=0; i<A3.rows(); i++) {
	for(int k=0; k<A3.cols(); k++) {
	  A3.e(i,k) = 0;
	  for(int j=0; j<i; j++) 
	    A3.e(i,k) += A1.ei(i,j)*A2.e(j,k);
	  for(int j=i; j<A1.cols(); j++) 
	    A3.e(i,k) += A1.ej(i,j)*A2.e(j,k);
	}
      }
    }

  template <class Type1, class Row, class Type3, class AT>
    inline void mult(const Matrix<Type1, AT> &A1, const Matrix<Symmetric<Row,Row>, AT> &A2, Matrix<Type3, AT> &A3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.cols() == A2.rows());
#endif
      for(int i=0; i<A3.rows(); i++) {
        for(int k=0; k<A3.cols(); k++) {
          A3.e(i,k) = 0;
          for(int j=0; j<k; j++) 
            A3.e(i,k) += A1.e(i,j)*A2.ej(j,k);
          for(int j=k; j<A1.cols(); j++) 
            A3.e(i,k) += A1.e(i,j)*A2.ei(j,k);
        }
      }
    }

  template <class Type1, class Type2, class AT>
    inline Matrix<General<Var,Var>, AT> operator*(const Matrix<Type1, AT> &A1, const Matrix<Type2, AT> &A2) {
      Matrix<General<Var,Var>, AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }

  template <class Row1, class Col1, class Row2, class Col2, class AT>
    inline Matrix<General<Row1,Col2>, AT> operator*(const Matrix<General<Row1,Col1>, AT> &A1, const Matrix<General<Row2,Col2>, AT> &A2) {
      Matrix<General<Row1,Col2>, AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }

  template <class Row1, class Col1, class AT>
    inline Matrix<General<Row1,Var>, AT> operator*(const Matrix<General<Row1,Col1>, AT> &A1, const Matrix<General<Ref,Ref>, AT> &A2) {
      Matrix<General<Row1,Var>, AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }

  template <class Row2, class Col2, class AT>
    inline Matrix<General<Var,Col2>, AT> operator*(const Matrix<General<Ref,Ref>, AT> &A1, const Matrix<General<Row2,Col2>, AT> &A2) {
      Matrix<General<Var,Col2>, AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }

  template <class Row1, class Col1, class Row2, class Col2, class AT>
    inline SquareMatrix<General<Row1,Col2>, AT> operator*(const SquareMatrix<General<Row1,Col1>, AT> &A1, const SquareMatrix<General<Row2,Col2>, AT> &A2) {
      SquareMatrix<General<Row1,Col2>, AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }

  template <class Row1, class Row2, class Col2, class AT>
    inline Matrix<General<Row1,Col2>, AT> operator*(const Matrix<Symmetric<Row1,Row1>, AT> &A1, const Matrix<General<Row2,Col2 >, AT> &A2) {
      Matrix<General<Row1,Col2>, AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }

  template <class Row1, class Col1, class Row2, class AT>
    inline Matrix<General<Row1,Row2>, AT> operator*(const Matrix<General<Row1,Col1 >, AT> &A1, const Matrix<Symmetric<Row2,Row2>, AT> &A2) {
      Matrix<General<Row1,Row2>, AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }

  /////////////////////////////////// end matmatmult //////////////////////////////
  
  /////////////////////////////////// matvecmult //////////////////////////////
  template <class Type1, class Type2, class Type3, class AT>
    inline void mult(const Matrix<Type1, AT> &A, const Vector<Type2, AT> &x, Vector<Type3, AT> &y) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.cols() == x.size());
#endif
      for(int i=0; i<y.size(); i++) {
        y.e(i) = 0;
        for(int j=0; j<x.size(); j++) 
          y.e(i) += A.e(i,j)*x.e(j);
      }
    }

  template <class AT, class Type1, class Type2>
    inline Vector<General<Var,Fixed<1> >, AT> operator*(const Matrix<Type1, AT> &A, const Vector<Type2, AT> &x) {
      Vector<General<Var,Fixed<1> >, AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }

  template <class Row1, class Col1, class Row2, class AT>
    inline Vector<General<Row1,Fixed<1> >, AT> operator*(const Matrix<General<Row1,Col1>, AT> &A, const Vector<General<Row2,Fixed<1> >, AT> &x) {
      Vector<General<Row1,Fixed<1> >, AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }

  template <class Row1, class Row2, class AT>
    inline Vector<General<Row1,Fixed<1> >, AT> operator*(const Matrix<Symmetric<Row1,Row1>, AT> &A, const Vector<General<Row2,Fixed<1> >, AT> &x) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.cols() == x.size());
#endif
      Vector<General<Row1,Fixed<1> >, AT> y(A.rows(),NONINIT);
      for(int i=0; i<y.size(); i++) {
	y.e(i) = 0;
	for(int j=0; j<i; j++) 
	  y.e(i) += A.ei(i,j)*x.e(j);
	for(int j=i; j<A.cols(); j++) 
	  y.e(i) += A.ej(i,j)*x.e(j);
      }
      return y;
    }

  /*! \brief Matrix-vector multiplication.
   *
   * This function computes the product of a matrix 
   * and a vector.
   * \return The product.
   * */
  template <class AT>
    inline Vector<General<Ref,Fixed<1> >, AT> operator*(const Matrix<General<Ref,Ref>, AT> &A, const Vector<General<Ref,Fixed<1> >, AT> &x) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.cols() == x.size());
#endif

      Vector<General<Ref,Fixed<1> >, AT> z(A.rows(),NONINIT);

      if(x.transposed()) {
        if(A.transposed()) {
          for(int i=0; i<z.size(); i++) {
            z.er(i)=0;
            for(int j=0; j<x.size(); j++)
              z.er(i) += A.et(i,j) * x.et(j);
          }
        } 
        else {
          for(int i=0; i<z.size(); i++) {
            z.er(i)=0;
            for(int j=0; j<x.size(); j++)
              z.er(i) += A.er(i,j) * x.et(j);
          }
        }
      }
      else {
        if(A.transposed()) {
          for(int i=0; i<z.size(); i++) {
            z.er(i)=0;
            for(int j=0; j<x.size(); j++)
              z.er(i) += A.et(i,j) * x.er(j);
          }
        } 
        else {
          for(int i=0; i<z.size(); i++) {
            z.er(i)=0;
            for(int j=0; j<x.size(); j++)
              z.er(i) += A.er(i,j) * x.er(j);
          }
        }
      } 

      return z;
    }

  /*! \brief Matrix-vector multiplication.
   *
   * This function computes the product of a sparse matrix
   * and a vector.
   * \return The product.
   * */
  template <class AT, class Type> 
    Vector<General<Ref,Fixed<1> >, AT> operator*(const Matrix<Sparse<Ref,Ref>, AT> &A, const Vector<Type, AT> &x) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.cols() == x.size());
#endif

      Vector<General<Ref,Fixed<1> >, AT> z(A.rows(),NONINIT);

      const AT *a = A();
      const int *ia = A.Ip();
      const int *ja = A.Jp();
      for(int i=0; i<A.rows(); i++) {
        z.e(i) = 0;
        for(int j=ia[i]; j<ia[i+1]; j++)
          z.e(i) += a[j]*x.e(ja[j]);
      }

      return z;
    }

  /////////////////////////////////// end matvecmult //////////////////////////////
  //
  /////////////////////////////////// rowvecmatmult //////////////////////////////
  template <class Type1, class Type2, class Type3, class AT>
    inline void mult(const Vector<Type1, AT> &x, const Matrix<Type2, AT> &A, RowVector<Type3, AT> &y) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x.size() == A.rows());
#endif
      for(int i=0; i<y.size(); i++) {
	y.e(i) = 0;
	for(int j=0; j<x.size(); j++) 
	  y.e(i) += x.e(j)*A.e(j,i);
      }
    }

  template <class AT, class Type1, class Type2>
    inline RowVector<General<Fixed<1>,Var>, AT> operator*(const RowVector<Type1, AT> &x, const Matrix<Type2, AT> &A) {
      RowVector<General<Fixed<1>,Var>, AT> y(A.cols(),NONINIT);
      mult(x,A,y);
      return y;
    }

  template <class Col1, class Row2, class Col2, class AT>
    inline RowVector<General<Fixed<1>,Col2>, AT> operator*(const RowVector<General<Fixed<1>,Col1>, AT> &x, const Matrix<General<Row2,Col2>, AT> &A) {
      RowVector<General<Fixed<1>,Col2>, AT> y(A.cols(),NONINIT);
      mult(x,A,y);
      return y;
    }

//  /*! \brief Vector-matrix multiplication.
//   *
//   * This function computes the product of a vector 
//   * and a matrix.
//   * \return The product.
//   * */
//  template <class AT, class Type> 
//    RowVector<General<Fixed<1>,Ref>, AT> operator*(const RowVector<General<Fixed<1>,Ref>, AT> &x, const Matrix<Type, AT> &A) {
//
//#ifndef FMATVEC_NO_SIZE_CHECK
//      assert(x.size() == A.rows());
//#endif
//
//      RowVector<General<Fixed<1>,Ref>, AT> z(A.cols(),NONINIT);
//
//      for(int i=0; i<z.size(); i++) {
//        z(i)=0;
//        for(int j=0; j<x.size(); j++)
//          z.e(i) += x.e(j) * A.e(j,i);
//      }
//
//      return z;
//    }

  /////////////////////////////////// end rowvecmatmult //////////////////////////////
  //
  /////////////////////////////////// end vecscalmult //////////////////////////////

  /*! \brief Vector-scalar multiplication.
   *
   * This function computes the product of a vector 
   * and a scalar.
   * \return The product.
   * */
  template <class Row, class AT>
    Vector<General<Row,Fixed<1> >, AT> operator*(const Vector<General<Row,Fixed<1> >, AT> &x, const AT& alpha) {

      Vector<General<Row,Fixed<1> >, AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
        y.e(i) = x.e(i)*alpha;

      return y;
    }

  /*! \brief Scalar-vector multiplication.
   *
   * \see operator*(const Vector<General<Ref,Fixed<1> >, AT>&x,const AT&).
   * */
  template <class Row, class AT>
    Vector<General<Row,Fixed<1> >, AT> operator*(const AT& alpha, const Vector<General<Row,Fixed<1> >, AT> &x) {

      Vector<General<Row,Fixed<1> >, AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
        y.e(i) = x.e(i)*alpha;

      return y;
    }

  /*! \brief Rowvector-scalar multiplication.
   *
   * This function computes the product of a rowvector 
   * and a scalar.
   * \return The product.
   * */
  template <class Col, class AT>
    RowVector<General<Fixed<1>,Col>, AT> operator*(const RowVector<General<Fixed<1>,Col>, AT> &x, const AT& alpha) {

      RowVector<General<Fixed<1>,Col>, AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
        y.e(i) = x.e(i)*alpha;

      return y;
    }

  /*! \brief Scalar-rowvector multiplication.
   *
   * \see operator*(const RowVector<General<Fixed<1>,Col>, AT>&, const AT&).
   * */
  template <class Col, class AT>
    RowVector<General<Fixed<1>,Col>, AT> operator*(const AT &alpha, const RowVector<General<Fixed<1>,Col>, AT> &x) {

      RowVector<General<Fixed<1>,Col>, AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
        y.e(i) = x.e(i)*alpha;

      return y;
    }

  /////////////////////////////////// end vecscalmult //////////////////////////////
  
  /////////////////////////////////// end vecvecmult //////////////////////////////

   /*! \brief Scalar product.
   *
   * This function computes the product of a vector and a rowvector. The result is
   * a scalar.
   * \return The product.
   * */
  template <class Col, class Row, class AT> 
    AT operator*(const RowVector<General<Fixed<1>,Col>, AT> &x, const Vector<General<Row,Fixed<1> >, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x.size() == y.size());
#endif

      AT res=0;

      for(int i=0; i<x.size(); i++)
        res+=x.e(i)*y.e(i);

      return res;
    }

  /*! \brief Scalar product.
   *
   * This function computes the scalar product of two vectors.
   * \return The scalar product.
   * */
  template <class Row, class AT>
    AT scalarProduct(const Vector<General<Row,Fixed<1> >, AT> &x, const Vector<General<Row,Fixed<1> >, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x.size()==y.size());
#endif

      AT res = 0;

      for(int i=0; i<x.size(); i++)
        res += x.e(i)*y.e(i);

      return res;
    }

  /*! \brief Cross product.
   *
   * This function computes the cross product of two vectors.
   * \return The cross product.
   * */
  template <class Row, class AT>
    Vector<General<Row,Fixed<1> >, AT> crossProduct(const Vector<General<Row,Fixed<1> >, AT> &x, const Vector<General<Row,Fixed<1> >, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x.size()==3);
      assert(y.size()==3);
#endif

      Vector<General<Row,Fixed<1> >, AT> z(3,NONINIT);

      z.e(0) = x.e(1)*y.e(2) - x.e(2)*y.e(1);
      z.e(1) = x.e(2)*y.e(0) - x.e(0)*y.e(2);
      z.e(2) = x.e(0)*y.e(1) - x.e(1)*y.e(0);

      return z;
    }

  /*! \brief Triple product.
   *
   * This function computes the triple product of three vectors.
   * \return The triple product.
   * */
  template <class Row, class AT>
    double tripleProduct(const Vector<General<Row,Fixed<1> >, AT> &a, const Vector<General<Row,Fixed<1> >, AT> &x, const Vector<General<Row,Fixed<1> >, AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(a.size()==3);
      assert(x.size()==3);
      assert(y.size()==3);
#endif

      return a.e(0)*(x.e(1)*y.e(2)-x.e(2)*y.e(1)) + a.e(1)*(x.e(2)*y.e(0)-x.e(0)*y.e(2)) + a.e(2)*(x.e(0)*y.e(1)-x.e(1)*y.e(0));
    }

  /////////////////////////////////// end vecvecmult //////////////////////////////
  //
  /////////////////////////////////// end matscalmult //////////////////////////////

  /*! \brief Matrix-scalar multiplication.
   *
   * This function computes the product of a matrix 
   * and a scalar.
   * \param A The matrix. 
   * \param alpha The scalar. 
   * \return The product.
   * */
  template <class Row, class Col, class AT>
    Matrix<General<Row,Col>, AT> operator*(const Matrix<General<Row,Col>, AT > &A, const AT &alpha) {

      Matrix<General<Row,Col>, AT> B(A.rows(),A.cols(),NONINIT);

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++) 
          B.e(i,j) = A.e(i,j)*alpha;

      return B;
    }

  /*! \brief Scalar-matrix multiplication.
   *
   * \see operator*(const Matrix<Type, AT >&, const AT&).
   * */
  template <class Row, class Col, class AT>
    Matrix<General<Row,Col>, AT> operator*(const AT &alpha, const Matrix<General<Row,Col>, AT > &A) {

      Matrix<General<Row,Col>, AT> B(A.rows(),A.cols(),NONINIT);

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++) 
          B.e(i,j) = A.e(i,j)*alpha;

      return B;
    }

  /////////////////////////////////// end matscalmult //////////////////////////////
 
  /////////////////////////////////// negation //////////////////////////////

  /*! \brief Negation.
   *
   * This function computes the negation of a vector.
   * \return The negation.
   * */
  template <class AT>
    Vector<General<Ref,Fixed<1> >, AT> operator-(const Vector<General<Ref,Fixed<1> >, AT> &x) {

      Vector<General<Ref,Fixed<1> >, AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        y.e(i)=-x.e(i);

      return y;
    }

  /*! \brief Negation.
   *
   * This function computes the negation of a rowvector.
   * \return The negation.
   * */
  template <class AT>
    RowVector<General<Fixed<1>,Ref>, AT> operator-(const RowVector<General<Fixed<1>,Ref>, AT> &a) {

      RowVector<General<Fixed<1>,Ref>, AT> c(a.size(),NONINIT);

      for(int i=0; i<a.size(); i++)
        c.e(i)=-a.e(i);

      return c;
    }

 /*! \brief Negation.
   *
   * This function computes the negation of a vector.
   * \return The negation.
   * */
  template <class AT>
    SquareMatrix<General<Ref,Ref>, AT> operator-(const SquareMatrix<General<Ref,Ref>, AT> &A) {

      SquareMatrix<General<Ref,Ref>, AT> B(A.size(),NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=0; j<A.cols(); j++)
          B.e(i,j)=-A.e(i,j);

      return B;
    }

  /////////////////////////////////// end negation //////////////////////////////
  /////////////////////////////////// transpose //////////////////////////////

  /*! \brief Transpose of a vector.
   *
   * This function computes the transpose of a vector. The result is a rowvector.
   * \return The transpose.
   * */
  template <class AT>
    RowVector<General<Fixed<1>,Ref>, AT> trans(const Vector<General<Ref,Fixed<1> >, AT> &x) {

      return RowVector<General<Fixed<1>,Ref>, AT>(x.m,x.lda,x.tp?false:true,x.memory,x.ele).copy();
    }


  /*! \brief Transpose of a rowvector.
   *
   * This function computes the transpose of a rowvector.
   * The result is a vector.
   * \param x A vector.
   * \return The transpose.
   * */
  template <class AT>
    Vector<General<Ref,Fixed<1> >, AT> trans(const RowVector<General<Fixed<1>,Ref>, AT> &x) {

      return Vector<General<Ref,Fixed<1> >, AT>(x.n,x.lda,x.tp?false:true,x.memory,x.ele).copy();
    }

  /*! \brief Transpose of a matrix.
   *
   * This function computes the transpose of a matrix.
   * \f[ \boldsymbol{A} \rightarrow \boldsymbol{A}^T  \f]
   * The result is a matrix.
   * \param A A general matrix.
   * \return The transpose.
   * */
  template <class AT>
    Matrix<General<Ref,Ref>, AT> trans(const Matrix<General<Ref,Ref>, AT> &A) {

      return Matrix<General<Ref,Ref>, AT>(A.n,A.m,A.lda,A.tp?false:true,A.memory,A.ele).copy();
    }

  /*! \brief Transpose of a matrix.
   *
   * This function computes the transpose of a square matrix.
   * \f[ \boldsymbol{A} \rightarrow \boldsymbol{A}^T  \f]
   * The result is a square matrix.
   * \param A A square matrix.
   * \return The transpose.
   * */
  template <class AT>
    SquareMatrix<General<Ref,Ref>, AT> trans(const SquareMatrix<General<Ref,Ref>, AT> &A) {

      return SquareMatrix<General<Ref,Ref>, AT>(A.n, A.lda, A.tp?false:true, A.memory, A.ele).copy();
    }

  /////////////////////////////////// end transpose //////////////////////////////


 
  ///////  /*! \brief Negation.
  ///////   *
  ///////   * This function computes the negation of a matrix
  ///////   * \return The negation. 
  ///////   * */
  ///////  template <class AT, class Type>
  ///////    Matrix<Type, AT> operator-(const Matrix<Type, AT> &A) {
  ///////
  ///////      Matrix<Type, AT> C(A.rows(),A.cols(),NONINIT);
  ///////
  ///////      for(int i=0; i<A.rows(); i++)
  ///////	for(int j=0; j<A.cols(); j++)
  ///////	  C(i,j)=-A(i,j);
  ///////
  ///////      return C;
  ///////    }


  /*! \brief Tilde operator
   *
   *  Example:
   * \code 
   * tilde(x);
   * \endcode 
   *  \f[
   *  x=\begin{pmatrix}x_1 \\ x_2 \\ x_3\end{pmatrix} \quad \Rightarrow \quad
   *  \tilde x = \begin{pmatrix} 0 & -x_3 & x_2 \\ x_3 &  0 & -x_1 \\ -x_2 & x_1 & 0\end{pmatrix}
   *  \f]
   * */
  template <class AT>
    SquareMatrix<General<Ref,Ref>, AT> tilde(const Vector<General<Ref,Fixed<1> >, AT> &x) {

      SquareMatrix<General<Ref,Ref>, AT> B(x.size(),NONINIT);

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


  ////////
  ////////   /*! \brief Matrix-matrix comparison.
  ////////   *
  ////////   * This function compares two matrices
  ////////   * \return True, if the matrices are identical, false otherwise.
  ////////   * */
  ////////  template <class AT, class Type1, class Type2> 
  ////////    bool operator==(const Matrix<Type1, AT> &A, const Matrix<Type2, AT > &B) {
  ////////
  ////////      if(&A == &B)
  ////////	return true;
  ////////
  ////////      if(A.rows() != B.rows() || A.cols() != B.cols())
  ////////	return false;
  ////////
  ////////      for(int i=0; i<A.rows(); i++) 
  ////////	for(int j=0; j<A.cols(); j++)
  ////////	  if(A(i,j)!=B(i,j))
  ////////	    return false;
  ////////
  ////////      return true;
  ////////    }
  ////////
  ////////   /*! \brief Matrix-matrix comparison.
  ////////   *
  ////////   * This function compares two matrices.
  ////////   * \return True, if the matrices are different, false otherwise.
  ////////   * */
  ////////  template <class AT, class Type1, class Type2> 
  ////////    bool operator!=(const Matrix<Type1, AT> &A, const Matrix<Type2, AT > &B) {
  ////////
  ////////      if(A.rows() != B.rows() || A.cols() != B.cols())
  ////////	return true;
  ////////
  ////////      for(int i=0; i<A.rows(); i++) 
  ////////	for(int j=0; j<A.cols(); j++)
  ////////	  if(A(i,j)!=B(i,j))
  ////////	    return true;
  ////////
  ////////      return false;
  ////////    }


  /*! \brief Maximum value.
   *
   * This function computes the maximum value of a vector.
   * \return The maximum value.
   * */
  template <class AT>
    AT max(const Vector<General<Ref,Fixed<1> >, AT> &x) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x.size() > 0);
#endif
      AT maximum = x.e(0);
      for(int i=1; i<x.size(); i++) {
        if(x.e(i) > maximum)
          maximum = x.e(i);
      }
      return maximum;
    }

  /*! \brief Index of maximum value.
   *
   * This function computes the index of the maximum value of a vector
   * \return The index of the maximum value.
   * */
  template <class AT>
    int maxIndex(const Vector<General<Ref,Fixed<1> >, AT> &x) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x.size() > 0);
#endif
      AT maximum = x.e(0);
      int index = 0;
      for(int i=1; i<x.size(); i++) {
        if(x.e(i) > maximum) {
          maximum = x.e(i);
          index = i;
        }
      }
      return index;
    }

  /*! \brief Minimum value.
   *
   * This function computes the minimum value of a vector.
   * \return The minimum value.
   * */
  template <class AT>
    AT min(const Vector<General<Ref,Fixed<1> >, AT> &x) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x.size() > 0);
#endif
      AT minimum = x.e(0);
      for(int i=1; i<x.size(); i++) {
        if(x.e(i) < minimum)
          minimum = x.e(i);
      }
      return minimum;
    }

  /*! \brief Index of minimum value.
   *
   * This function computes the index of the minimum value of a vector
   * \return The index of the minimum value.
   * */
  template <class AT>
    int minIndex(const Vector<General<Ref,Fixed<1> >, AT> &x) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x.size() > 0);
#endif
      AT minimum = x.e(0);
      int index = 0;
      for(int i=1; i<x.size(); i++) {
        if(x.e(i) < minimum) {
          minimum = x.e(i);
          index = i;
        }
      }
      return index;
    }

  // HR 28.09.2006
  /*! \brief Bubble Sort Algorithm (stable sorting Algorithm )
   *
   * Values of rowvectors of Matrix A are sorted in ascending order 
   * \param A_ Matrix to be sorted
   * \param PivotCol Column of A used as sorting index 
   */
  template <class AT>
    Matrix<General<Ref,Ref>, AT> bubbleSort(const Matrix<General<Ref,Ref>, AT> &A_, int PivotCol) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A_.rows()>0);
      assert(A_.cols()>0);
      assert(A_.cols()> PivotCol);
#endif
      Matrix<General<Ref,Ref>, AT> A = A_.copy();
      int i,j,N;
      RowVector<General<Fixed<1>,Ref>, AT> tmp;
      N = A.rows();
      for (i=1; i<= N-1; i++) {
        for (j=0; j< N-1; j++) {
          if (A(j,PivotCol) >  A(j+1,PivotCol)) {
            tmp = A.row(j);  
            A.row(j) = A.row(j+1);
            A.row(j+1) = tmp;
          }
        }
      }
      return A;
    }

  // HR 28.09.2006
  /*! internal function of QuickSortMedian */
  template <class AT>
    void quicksortmedian_intern(Matrix<General<Ref,Ref>, AT> &A, int PivotCol, RowVector<General<Fixed<1>,Ref>, AT> &tmp, int l, int r){

      if(r>l){
        int i=l-1, j=r;
        //RowVec tmp;     
        if(r-l > 3){ //Median of three
          int m=l+(r-l)/2;
          if(A(l,PivotCol)>A(m,PivotCol))
          { tmp= A.row(l);
            A.row(l) = A.row(m);
            A.row(m) =tmp; }
          if(A(l,PivotCol)> A(r,PivotCol))
          { tmp= A.row(l);
            A.row(l) = A.row(r);
            A.row(r) =tmp; }
          else if(A(r,PivotCol) > A(m,PivotCol) )
          { tmp= A.row(r);
            A.row(r) = A.row(m);
            A.row(m) =tmp; }
        }

        for(;;){
          while(A(++i,PivotCol) < A(r,PivotCol));
          while(A(--j,PivotCol) > A(r,PivotCol) && j>i);
          if(i>=j) break;
          tmp= A.row(i);
          A.row(i) = A.row(j);
          A.row(j) =tmp;
        }
        tmp= A.row(i);
        A.row(i) = A.row(r);
        A.row(r) =tmp;
        quicksortmedian_intern(A, PivotCol,tmp, l, i-1);
        quicksortmedian_intern(A, PivotCol,tmp, i+1, r);
      }
    }

  // HR 28.09.2006
  /*! \brief Modified QuickSort Algorithm (unstable sorting Algorithm )
   *
   * Values of rowvectors of Matrix A are sorted in ascending order 
   * unstabel but very quick sorting Algorithm
   * Pivot Elements as 'Median of three'
   * \param A_ Matrix to be sorted
   * \param PivotCol Column of A used as sorting index 
   */
  template <class AT>
    Matrix<General<Ref,Ref>, AT> quickSortMedian(const Matrix<General<Ref,Ref>, AT> &A_, int PivotCol) {
      Matrix<General<Ref,Ref>, AT> A = A_.copy();
      int N = A.rows();
      RowVector<General<Fixed<1>,Ref>, AT> tmp;
      quicksortmedian_intern(A, PivotCol,tmp, 0, N-1);
      return A;
    }

  /*! \brief Count nonzero elements.
   *
   * This function counts the nonzero elements of a matrix. ALL diagonal
   * elements are treated as NONZERO!!! (See the storage format for sparse matrix)
   * \return The number of nonzero or diagonal elemets.
   * */
  template <class AT> int countElements(const SquareMatrix<General<Ref,Ref>, AT> &A) { 
    int k=0;
    for(int i=0; i<A.size(); i++) {
      for(int j=0; j<A.size(); j++) {
        if(fabs(A(i,j))>1e-16 || i==j) {
          k++;
        }
      }
    }
    return k;
  }


}

#endif

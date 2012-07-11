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
#include "vector.h"
#include "row_vector.h"
namespace fmatvec {

/////////////////////////////////// vecvecadd //////////////////////////////

  /*! \brief Vector-vector addition.
   *
   * This function computes the sum of two vectors
   * \return The sum.
   * */
  template <class AT>
    inline Vector<General,Ref,Fixed<1>,AT> operator+(const Vector<General,Ref,Fixed<1>,AT> &x, const Vector<General,Ref,Fixed<1>,AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(y.size() == x.size());
#endif

      Vector<General,Ref,Fixed<1>,AT> z(x.size(),NONINIT);

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

  // Vector-Vector
  template <class Type1, class Row1, class Type2, class Row2, class Type3, class Row3, class AT> 
    inline void add(const Vector<Type1,Row1,Fixed<1>,AT> &a1, const Vector<Type2,Row2,Fixed<1>,AT> &a2, Vector<Type3,Row3,Fixed<1>,AT> &a3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(a1.size() == a2.size());
#endif
      for(int i=0; i<a1.size(); i++)
          a3.e(i) = a1.e(i) + a2.e(i);
    }

  template <class Type1, class Row1, class Type2, class Row2, class AT> 
    inline void add(Vector<Type1,Row1,Fixed<1>,AT> &a1, const Vector<Type2,Row2,Fixed<1>,AT> &a2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(a1.size() == a2.size());
#endif
      for(int i=0; i<a1.size(); i++)
          a1.e(i) += a2.e(i);
    }

  template <class Type1, class Row1, class Type2, class Row2, class Type3, class Row3, class AT> 
    inline void sub(const Vector<Type1,Row1,Fixed<1>,AT> &a1, const Vector<Type2,Row2,Fixed<1>,AT> &a2, Vector<Type3,Row3,Fixed<1>,AT> &a3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(a1.size() == a2.size());
#endif
      for(int i=0; i<a1.size(); i++)
          a3.e(i) = a1.e(i) - a2.e(i);
    }

  template <class Type1, class Row1, class Type2, class Row2, class AT> 
    inline void sub(Vector<Type1,Row1,Fixed<1>,AT> &a1, const Vector<Type2,Row2,Fixed<1>,AT> &a2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(a1.size() == a2.size());
#endif
      for(int i=0; i<a1.size(); i++)
          a1.e(i) -= a2.e(i);
    }

  // Addition
  template <class AT, class Type1, class Type2, class Row1, class Row2> 
    inline Vector<General,Var,Fixed<1>,AT> operator+(const Vector<Type1,Row1,Fixed<1>,AT> &a, const Vector<Type2,Row2,Fixed<1>,AT> &b) {
      Vector<General,Var,Fixed<1>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, class Type1, class Type2, class Row> 
    inline Vector<General,Row,Fixed<1>,AT> operator+(const Vector<Type1,Row,Fixed<1>,AT> &a, const Vector<Type2,Row,Fixed<1>,AT> &b) {
      Vector<General,Row,Fixed<1>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, int M, class Type1, class Type2, class Row2> 
    inline Vector<General,Fixed<M>,Fixed<1>,AT> operator+(const Vector<Type1,Fixed<M>,Fixed<1>,AT> &a, const Vector<Type2,Row2,Fixed<1>,AT> &b) {
      Vector<General,Fixed<M>,Fixed<1>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, int M, class Type1, class Type2> 
    inline Vector<General,Fixed<M>,Fixed<1>,AT> operator+(const Vector<Type1,Fixed<M>,Fixed<1>,AT> &a, const Vector<Type2,Fixed<M>,Fixed<1>,AT> &b) {
      Vector<General,Fixed<M>,Fixed<1>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, class Type, class Row1, class Row2> 
    inline Vector<General,Var,Fixed<1>,AT> operator+(const Vector<Type,Row1,Fixed<1>,AT> &a, const Vector<Type,Row2,Fixed<1>,AT> &b) {
      Vector<General,Var,Fixed<1>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, class Type, class Row> 
    inline Vector<General,Row,Fixed<1>,AT> operator+(const Vector<Type,Row,Fixed<1>,AT> &a, const Vector<Type,Row,Fixed<1>,AT> &b) {
      Vector<General,Row,Fixed<1>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, int M, class Type, class Row2> 
    inline Vector<General,Fixed<M>,Fixed<1>,AT> operator+(const Vector<Type,Fixed<M>,Fixed<1>,AT> &a, const Vector<Type,Row2,Fixed<1>,AT> &b) {
      Vector<General,Fixed<M>,Fixed<1>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, int M, class Type> 
    inline Vector<General,Fixed<M>,Fixed<1>,AT> operator+(const Vector<Type,Fixed<M>,Fixed<1>,AT> &a, const Vector<Type,Fixed<M>,Fixed<1>,AT> &b) {
      Vector<General,Fixed<M>,Fixed<1>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, class Type1, class Type2, class Row1, class Row2> 
    inline Vector<General,Row1,Fixed<1>,AT> operator+=(const Vector<Type1,Row1,Fixed<1>,AT> &a_, const Vector<Type2,Row2,Fixed<1>,AT> &b) {
      Vector<Type1,Row1,Fixed<1>,AT> &a = const_cast<Vector<Type1,Row1,Fixed<1>,AT> &>(a_);
      add(a,b);
      return a;
    }

///  // Subtraction
  template <class AT, class Type1, class Type2, class Row1, class Row2> 
    inline Vector<General,Var,Fixed<1>,AT> operator-(const Vector<Type1,Row1,Fixed<1>,AT> &a, const Vector<Type2,Row2,Fixed<1>,AT> &b) {
      Vector<General,Var,Fixed<1>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, class Type1, class Type2, class Row> 
    inline Vector<General,Row,Fixed<1>,AT> operator-(const Vector<Type1,Row,Fixed<1>,AT> &a, const Vector<Type2,Row,Fixed<1>,AT> &b) {
      Vector<General,Row,Fixed<1>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, int M, class Type1, class Type2, class Row2> 
    inline Vector<General,Fixed<M>,Fixed<1>,AT> operator-(const Vector<Type1,Fixed<M>,Fixed<1>,AT> &a, const Vector<Type2,Row2,Fixed<1>,AT> &b) {
      Vector<General,Fixed<M>,Fixed<1>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, int M, class Type1, class Type2> 
    inline Vector<General,Fixed<M>,Fixed<1>,AT> operator-(const Vector<Type1,Fixed<M>,Fixed<1>,AT> &a, const Vector<Type2,Fixed<M>,Fixed<1>,AT> &b) {
      Vector<General,Fixed<M>,Fixed<1>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, class Type, class Row1, class Row2> 
    inline Vector<General,Var,Fixed<1>,AT> operator-(const Vector<Type,Row1,Fixed<1>,AT> &a, const Vector<Type,Row2,Fixed<1>,AT> &b) {
      Vector<General,Var,Fixed<1>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, class Type, class Row> 
    inline Vector<General,Row,Fixed<1>,AT> operator-(const Vector<Type,Row,Fixed<1>,AT> &a, const Vector<Type,Row,Fixed<1>,AT> &b) {
      Vector<General,Row,Fixed<1>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, int M, class Type, class Row2> 
    inline Vector<General,Fixed<M>,Fixed<1>,AT> operator-(const Vector<Type,Fixed<M>,Fixed<1>,AT> &a, const Vector<Type,Row2,Fixed<1>,AT> &b) {
      Vector<General,Fixed<M>,Fixed<1>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, int M, class Type> 
    inline Vector<General,Fixed<M>,Fixed<1>,AT> operator-(const Vector<Type,Fixed<M>,Fixed<1>,AT> &a, const Vector<Type,Fixed<M>,Fixed<1>,AT> &b) {
      Vector<General,Fixed<M>,Fixed<1>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, class Type1, class Type2, class Row1, class Row2> 
    inline Vector<General,Row1,Fixed<1>,AT> operator-=(const Vector<Type1,Row1,Fixed<1>,AT> &a_, const Vector<Type2,Row2,Fixed<1>,AT> &b) {
      Vector<Type1,Row1,Fixed<1>,AT> &a = const_cast<Vector<Type1,Row1,Fixed<1>,AT> &>(a_);
      sub(a,b);
      return a;
    }

  // RowVector-Vector
  template <class Type1, class Col1, class Type2, class Col2, class Type3, class Col3, class AT> 
    inline void add(const RowVector<Type1,Fixed<1>,Col1,AT> &a1, const RowVector<Type2,Fixed<1>,Col2,AT> &a2, RowVector<Type3,Fixed<1>,Col3,AT> &a3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(a1.size() == a2.size());
#endif
      for(int i=0; i<a1.size(); i++)
          a3.e(i) = a1.e(i) + a2.e(i);
    }

  template <class Type1, class Col1, class Type2, class Col2, class AT> 
    inline void add(RowVector<Type1,Fixed<1>,Col1,AT> &a1, const RowVector<Type2,Fixed<1>,Col2,AT> &a2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(a1.size() == a2.size());
#endif
      for(int i=0; i<a1.size(); i++)
          a1.e(i) += a2.e(i);
    }

  template <class Type1, class Col1, class Type2, class Col2, class Type3, class Col3, class AT> 
    inline void sub(const RowVector<Type1,Fixed<1>,Col1,AT> &a1, const RowVector<Type2,Fixed<1>,Col2,AT> &a2, RowVector<Type3,Fixed<1>,Col3,AT> &a3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(a1.size() == a2.size());
#endif
      for(int i=0; i<a1.size(); i++)
          a3.e(i) = a1.e(i) - a2.e(i);
    }

  template <class Type1, class Col1, class Type2, class Col2, class AT> 
    inline void sub(RowVector<Type1,Fixed<1>,Col1,AT> &a1, const RowVector<Type2,Fixed<1>,Col2,AT> &a2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(a1.size() == a2.size());
#endif
      for(int i=0; i<a1.size(); i++)
          a1.e(i) -= a2.e(i);
    }

  // Addition
  template <class AT, class Type1, class Type2, class Col1, class Col2> 
    inline RowVector<General,Var,Fixed<1>,AT> operator+(const RowVector<Type1,Fixed<1>,Col1,AT> &a, const RowVector<Type2,Fixed<1>,Col2,AT> &b) {
      RowVector<General,Fixed<1>,Var,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, class Type1, class Type2, class Col> 
    inline RowVector<General,Fixed<1>,Col,AT> operator+(const RowVector<Type1,Fixed<1>,Col,AT> &a, const RowVector<Type2,Fixed<1>,Col,AT> &b) {
      Vector<General,Fixed<1>,Col,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, int N, class Type1, class Type2, class Col2> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator+(const RowVector<Type1,Fixed<1>,Fixed<N>,AT> &a, const RowVector<Type2,Fixed<1>,Col2,AT> &b) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, int N, class Type1, class Type2> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator+(const RowVector<Type1,Fixed<1>,Fixed<N>,AT> &a, const RowVector<Type2,Fixed<1>,Fixed<N>,AT> &b) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, class Type, class Col1, class Col2> 
    inline RowVector<General,Var,Fixed<1>,AT> operator+(const RowVector<Type,Fixed<1>,Col1,AT> &a, const RowVector<Type,Fixed<1>,Col2,AT> &b) {
      RowVector<General,Fixed<1>,Var,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, class Type, class Col> 
    inline RowVector<General,Fixed<1>,Col,AT> operator+(const RowVector<Type,Fixed<1>,Col,AT> &a, const RowVector<Type,Fixed<1>,Col,AT> &b) {
      Vector<General,Fixed<1>,Col,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, int N, class Type, class Col2> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator+(const RowVector<Type,Fixed<1>,Fixed<N>,AT> &a, const RowVector<Type,Fixed<1>,Col2,AT> &b) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, int N, class Type> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator+(const RowVector<Type,Fixed<1>,Fixed<N>,AT> &a, const RowVector<Type,Fixed<1>,Fixed<N>,AT> &b) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> c(a.size(),NONINIT);
      add(a,b,c);
      return c;
    }
  template <class AT, class Type1, class Type2, class Col1, class Col2> 
    inline RowVector<General,Fixed<1>,Col1,AT> operator+=(const RowVector<Type1,Fixed<1>,Col1,AT> &a_, const RowVector<Type2,Fixed<1>,Col2,AT> &b) {
      RowVector<Type1,Fixed<1>,Col1,AT> &a = const_cast<RowVector<Type1,Fixed<1>,Col1,AT> &>(a_);
      add(a,b);
      return a;
    }

///  // Subtraction
  template <class AT, class Type1, class Type2, class Col1, class Col2> 
    inline RowVector<General,Var,Fixed<1>,AT> operator-(const RowVector<Type1,Fixed<1>,Col1,AT> &a, const RowVector<Type2,Fixed<1>,Col2,AT> &b) {
      RowVector<General,Fixed<1>,Var,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, class Type1, class Type2, class Col> 
    inline RowVector<General,Fixed<1>,Col,AT> operator-(const RowVector<Type1,Fixed<1>,Col,AT> &a, const RowVector<Type2,Fixed<1>,Col,AT> &b) {
      Vector<General,Fixed<1>,Col,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, int N, class Type1, class Type2, class Col2> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator-(const RowVector<Type1,Fixed<1>,Fixed<N>,AT> &a, const RowVector<Type2,Fixed<1>,Col2,AT> &b) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, int N, class Type1, class Type2> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator-(const RowVector<Type1,Fixed<1>,Fixed<N>,AT> &a, const RowVector<Type2,Fixed<1>,Fixed<N>,AT> &b) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, class Type, class Col1, class Col2> 
    inline RowVector<General,Var,Fixed<1>,AT> operator-(const RowVector<Type,Fixed<1>,Col1,AT> &a, const RowVector<Type,Fixed<1>,Col2,AT> &b) {
      RowVector<General,Fixed<1>,Var,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, class Type, class Col> 
    inline RowVector<General,Fixed<1>,Col,AT> operator-(const RowVector<Type,Fixed<1>,Col,AT> &a, const RowVector<Type,Fixed<1>,Col,AT> &b) {
      Vector<General,Fixed<1>,Col,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, int N, class Type, class Col2> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator-(const RowVector<Type,Fixed<1>,Fixed<N>,AT> &a, const RowVector<Type,Fixed<1>,Col2,AT> &b) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, int N, class Type> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator-(const RowVector<Type,Fixed<1>,Fixed<N>,AT> &a, const RowVector<Type,Fixed<1>,Fixed<N>,AT> &b) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> c(a.size(),NONINIT);
      sub(a,b,c);
      return c;
    }
  template <class AT, class Type1, class Type2, class Col1, class Col2> 
    inline RowVector<General,Fixed<1>,Col1,AT> operator-=(const RowVector<Type1,Fixed<1>,Col1,AT> &a_, const RowVector<Type2,Fixed<1>,Col2,AT> &b) {
      RowVector<Type1,Fixed<1>,Col1,AT> &a = const_cast<RowVector<Type1,Fixed<1>,Col1,AT> &>(a_);
      sub(a,b);
      return a;
    }

  // Matrix-Matrix
  template <class Type1, class Row1, class Col1, class Type2, class Row2, class Col2, class Type3, class Row3, class Col3, class AT> 
    inline void add(const Matrix<Type1,Row1,Col1,AT> &A1, const Matrix<Type2,Row2,Col2,AT> &A2, Matrix<Type3,Row3,Col3,AT> &A3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      for(int i=0; i<A1.rows(); i++)
        for(int j=0; j<A2.cols(); j++)
          A3.e(i,j) = A1.e(i,j) + A2.e(i,j);
    }

  template <class Row1, class Row2, class Row3, class AT> 
    inline void add(const Matrix<Symmetric,Row1,Row1,AT> &A1, const Matrix<Symmetric,Row2,Row2,AT> &A2, Matrix<Symmetric,Row3,Row3,AT> &A3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      for(int i=0; i<A1.size(); i++)
        for(int j=0; j<A2.size(); j++)
          A3.ej(i,j) = A1.ej(i,j) + A2.ej(i,j);
    }

  template <class Type1, class Row1, class Col1, class Type2, class Row2, class Col2, class AT> 
    inline void add(Matrix<Type1,Row1,Col1,AT> &A1, const Matrix<Type2,Row2,Col2,AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      for(int i=0; i<A1.rows(); i++)
        for(int j=0; j<A2.cols(); j++)
          A1.e(i,j) += A2.e(i,j);
    }

  template <class Row1, class Row2, class AT> 
    inline void add(Matrix<Symmetric,Row1,Row1,AT> &A1, const Matrix<Symmetric,Row2,Row2,AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      for(int i=0; i<A1.size(); i++)
        for(int j=0; j<A2.size(); j++)
          A1.ej(i,j) += A2.ej(i,j);
    }

  template <class Type1, class Row1, class Col1, class Type2, class Row2, class Col2, class Type3, class Row3, class Col3, class AT> 
    inline void sub(const Matrix<Type1,Row1,Col1,AT> &A1, const Matrix<Type2,Row2,Col2,AT> &A2, Matrix<Type3,Row3,Col3,AT> &A3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      for(int i=0; i<A1.rows(); i++)
        for(int j=0; j<A2.cols(); j++)
          A3.e(i,j) = A1.e(i,j) - A2.e(i,j);
    }

  template <class Row1, class Row2, class Row3, class AT> 
    inline void sub(const Matrix<Symmetric,Row1,Row1,AT> &A1, const Matrix<Symmetric,Row2,Row2,AT> &A2, Matrix<Symmetric,Row3,Row3,AT> &A3) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      for(int i=0; i<A1.size(); i++)
        for(int j=0; j<A2.size(); j++)
          A3.ej(i,j) = A1.ej(i,j) - A2.ej(i,j);
    }

  template <class Type1, class Row1, class Col1, class Type2, class Row2, class Col2, class AT> 
    inline void sub(Matrix<Type1,Row1,Col1,AT> &A1, const Matrix<Type2,Row2,Col2,AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      for(int i=0; i<A1.rows(); i++)
        for(int j=0; j<A2.cols(); j++)
          A1.e(i,j) -= A2.e(i,j);
    }

  template <class Row1, class Row2, class AT> 
    inline void sub(Matrix<Symmetric,Row1,Row1,AT> &A1, const Matrix<Symmetric,Row2,Row2,AT> &A2) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A1.rows() == A2.rows());
      assert(A1.cols() == A2.cols());
#endif
      for(int i=0; i<A1.size(); i++)
        for(int j=0; j<A2.size(); j++)
          A1.ej(i,j) -= A2.ej(i,j);
    }


  // Addition
  // Type1 Type2
  template <class AT, class Type1, class Type2, class Row1, class Col1, class Row2, class Col2> 
    inline Matrix<General,Var,Var,AT> operator+(const Matrix<Type1,Row1,Col1,AT> &A, const Matrix<Type2,Row2,Col2,AT> &B) {
      Matrix<General,Var,Var,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, class Type1, class Type2, class Row, class Col>
    inline Matrix<General,Row,Col,AT> operator+(const Matrix<Type1,Row,Col,AT> &A, const Matrix<Type2,Row,Col,AT> &B) {
      Matrix<General,Row,Col,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, class Type1, class Type2, class Row, class Col> 
    inline Matrix<General,Fixed<M>,Var,AT> operator+(const Matrix<Type1,Fixed<M>,Var,AT> &A, const Matrix<Type2,Row,Col,AT> &B) {
      Matrix<General,Fixed<M>,Var,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Var,AT> operator+(const Matrix<Type1,Fixed<M>,Var,AT> &A, const Matrix<Type2,Fixed<M>,Var,AT> &B) {
      Matrix<General,Fixed<M>,Var,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int N, class Type1, class Type2, class Row, class Col> 
    inline Matrix<General,Var,Fixed<N>,AT> operator+(const Matrix<Type1,Var,Fixed<N>,AT> &A, const Matrix<Type2,Row,Col,AT> &B) {
      Matrix<General,Var,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int N, class Type1, class Type2> 
    inline Matrix<General,Var,Fixed<N>,AT> operator+(const Matrix<Type1,Var,Fixed<N>,AT> &A, const Matrix<Type2,Var,Fixed<N>,AT> &B) {
      Matrix<General,Var,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2, class Row, class Col> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type2,Row,Col,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type2,Fixed<M>,Var,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type2,Var,Fixed<N>,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type2,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2, class Row, class Col> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type2,Row,Col,AT> &A, const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type2,Fixed<M>,Var,AT> &A, const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type2,Var,Fixed<N>,AT> &A, const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  // Type
  template <class AT, class Type, class Row1, class Col1, class Row2, class Col2> 
    inline Matrix<Type,Var,Var,AT> operator+(const Matrix<Type,Row1,Col1,AT> &A, const Matrix<Type,Row2,Col2,AT> &B) {
      Matrix<Type,Var,Var,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, class Type, class Row, class Col> 
    inline Matrix<Type,Row,Col,AT> operator+(const Matrix<Type,Row,Col,AT> &A, const Matrix<Type,Row,Col,AT> &B) {
      Matrix<Type,Row,Col,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, class Type, class Row, class Col> 
    inline Matrix<Type,Fixed<M>,Var,AT> operator+(const Matrix<Type,Fixed<M>,Var,AT> &A, const Matrix<Type,Row,Col,AT> &B) {
      Matrix<Type,Fixed<M>,Var,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, class Type> 
    inline Matrix<Type,Fixed<M>,Var,AT> operator+(const Matrix<Type,Fixed<M>,Var,AT> &A, const Matrix<Type,Fixed<M>,Var,AT> &B) {
      Matrix<Type,Fixed<M>,Var,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int N, class Type, class Row, class Col> 
    inline Matrix<Type,Var,Fixed<N>,AT> operator+(const Matrix<Type,Var,Fixed<N>,AT> &A, const Matrix<Type,Row,Col,AT> &B) {
      Matrix<Type,Var,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int N, class Type> 
    inline Matrix<Type,Var,Fixed<N>,AT> operator+(const Matrix<Type,Var,Fixed<N>,AT> &A, const Matrix<Type,Var,Fixed<N>,AT> &B) {
      Matrix<Type,Var,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type, class Row, class Col> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type,Row,Col,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type,Fixed<M>,Var,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type,Var,Fixed<N>,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type, class Row, class Col> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type,Row,Col,AT> &A, const Matrix<Type,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type,Fixed<M>,Var,AT> &A, const Matrix<Type,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator+(const Matrix<Type,Var,Fixed<N>,AT> &A, const Matrix<Type,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      add(A,B,C);
      return C;
    }
  template <class AT, class Type1, class Row1, class Col1, class Type2, class Row2, class Col2> 
    inline Matrix<Type1,Row1,Col1,AT>& operator+=(const Matrix<Type1,Row1,Col1,AT> &A_, const Matrix<Type2,Row2,Col2,AT> &B) {
      Matrix<Type1,Row1,Col1,AT> &A = const_cast<Matrix<Type1,Row1,Col1,AT> &>(A_);
      add(A,B);
      return A;
    }

  // Subtraction
  // Type1 Type2
  template <class AT, class Type1, class Type2, class Row1, class Col1, class Row2, class Col2> 
    inline Matrix<General,Var,Var,AT> operator-(const Matrix<Type1,Row1,Col1,AT> &A, const Matrix<Type2,Row2,Col2,AT> &B) {
      Matrix<General,Var,Var,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, class Type1, class Type2, class Row, class Col>
    inline Matrix<General,Row,Col,AT> operator-(const Matrix<Type1,Row,Col,AT> &A, const Matrix<Type2,Row,Col,AT> &B) {
      Matrix<General,Row,Col,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, class Type1, class Type2, class Row, class Col> 
    inline Matrix<General,Fixed<M>,Var,AT> operator-(const Matrix<Type1,Fixed<M>,Var,AT> &A, const Matrix<Type2,Row,Col,AT> &B) {
      Matrix<General,Fixed<M>,Var,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Var,AT> operator-(const Matrix<Type1,Fixed<M>,Var,AT> &A, const Matrix<Type2,Fixed<M>,Var,AT> &B) {
      Matrix<General,Fixed<M>,Var,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int N, class Type1, class Type2, class Row, class Col> 
    inline Matrix<General,Var,Fixed<N>,AT> operator-(const Matrix<Type1,Var,Fixed<N>,AT> &A, const Matrix<Type2,Row,Col,AT> &B) {
      Matrix<General,Var,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int N, class Type1, class Type2> 
    inline Matrix<General,Var,Fixed<N>,AT> operator-(const Matrix<Type1,Var,Fixed<N>,AT> &A, const Matrix<Type2,Var,Fixed<N>,AT> &B) {
      Matrix<General,Var,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2, class Row, class Col> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type2,Row,Col,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type2,Fixed<M>,Var,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type2,Var,Fixed<N>,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type2,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2, class Row, class Col> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type2,Row,Col,AT> &A, const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type2,Fixed<M>,Var,AT> &A, const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type2,Var,Fixed<N>,AT> &A, const Matrix<Type1,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  // Type
  template <class AT, class Type, class Row1, class Col1, class Row2, class Col2> 
    inline Matrix<Type,Var,Var,AT> operator-(const Matrix<Type,Row1,Col1,AT> &A, const Matrix<Type,Row2,Col2,AT> &B) {
      Matrix<Type,Var,Var,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, class Type, class Row, class Col> 
    inline Matrix<Type,Row,Col,AT> operator-(const Matrix<Type,Row,Col,AT> &A, const Matrix<Type,Row,Col,AT> &B) {
      Matrix<Type,Row,Col,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, class Type, class Row, class Col> 
    inline Matrix<Type,Fixed<M>,Var,AT> operator-(const Matrix<Type,Fixed<M>,Var,AT> &A, const Matrix<Type,Row,Col,AT> &B) {
      Matrix<Type,Fixed<M>,Var,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, class Type> 
    inline Matrix<Type,Fixed<M>,Var,AT> operator-(const Matrix<Type,Fixed<M>,Var,AT> &A, const Matrix<Type,Fixed<M>,Var,AT> &B) {
      Matrix<Type,Fixed<M>,Var,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int N, class Type, class Row, class Col> 
    inline Matrix<Type,Var,Fixed<N>,AT> operator-(const Matrix<Type,Var,Fixed<N>,AT> &A, const Matrix<Type,Row,Col,AT> &B) {
      Matrix<Type,Var,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int N, class Type> 
    inline Matrix<Type,Var,Fixed<N>,AT> operator-(const Matrix<Type,Var,Fixed<N>,AT> &A, const Matrix<Type,Var,Fixed<N>,AT> &B) {
      Matrix<Type,Var,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type, class Row, class Col> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type,Row,Col,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type,Fixed<M>,Var,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type,Var,Fixed<N>,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type,Fixed<M>,Fixed<N>,AT> &A, const Matrix<Type,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type, class Row, class Col> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type,Row,Col,AT> &A, const Matrix<Type,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type,Fixed<M>,Var,AT> &A, const Matrix<Type,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, int M, int N, class Type> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator-(const Matrix<Type,Var,Fixed<N>,AT> &A, const Matrix<Type,Fixed<M>,Fixed<N>,AT> &B) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> C(A.rows(),A.cols(),NONINIT);
      sub(A,B,C);
      return C;
    }
  template <class AT, class Type1, class Row1, class Col1, class Type2, class Row2, class Col2> 
    inline Matrix<Type1,Row1,Col1,AT>& operator-=(const Matrix<Type1,Row1,Col1,AT> &A_, const Matrix<Type2,Row2,Col2,AT> &B) {
      Matrix<Type1,Row1,Col1,AT> &A = const_cast<Matrix<Type1,Row1,Col1,AT> &>(A_);
      sub(A,B);
      return A;
    }

//  SquareMatrix-SquareMatrix
  template <class AT, class Row1, class Row2> 
    inline SquareMatrix<General,Var,Var,AT> operator+(const SquareMatrix<General,Row1,Row1,AT> &A1, const SquareMatrix<General,Row2,Row2,AT> &A2) {
      SquareMatrix<General,Var,Var,AT> A3(A1.rows(),A1.cols(),NONINIT);
      add(A1,A2,A3);
      return A3;
    }
  template <class AT, class Row> 
    inline SquareMatrix<General,Row,Row,AT> operator+(const SquareMatrix<General,Row,Row,AT> &A1, const SquareMatrix<General,Row,Row,AT> &A2) {
      SquareMatrix<General,Row,Row,AT> A3(A1.rows(),A1.cols(),NONINIT);
      add(A1,A2,A3);
      return A3;
    }

//  SquareMatrix-SquareMatrix
  template <class AT, class Row1, class Row2> 
    inline SquareMatrix<General,Var,Var,AT> operator-(const SquareMatrix<General,Row1,Row1,AT> &A1, const SquareMatrix<General,Row2,Row2,AT> &A2) {
      SquareMatrix<General,Var,Var,AT> A3(A1.rows(),A1.cols(),NONINIT);
      sub(A1,A2,A3);
      return A3;
    }
  template <class AT, class Row> 
    inline SquareMatrix<General,Row,Row,AT> operator-(const SquareMatrix<General,Row,Row,AT> &A1, const SquareMatrix<General,Row,Row,AT> &A2) {
      SquareMatrix<General,Row,Row,AT> A3(A1.rows(),A1.cols(),NONINIT);
      sub(A1,A2,A3);
      return A3;
    }

  //////      
  //////      /////////////////////////////////// end vecvecadd //////////////////////////////
  //////      
  //////      /////////////////////////////////// matmatmult //////////////////////////////
  //////      
  template <class Type1, class Row1, class Col1, class Type2, class Row2, class Type3, class Row3, class AT> 
    inline void mult(const Matrix<Type1,Row1,Col1,AT> &A, const Vector<Type2,Row2,Fixed<1>,AT> &x, Vector<Type3,Row3,Fixed<1>,AT> &y) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A.cols() == x.size());
#endif
      for(int i=0; i<y.size(); i++) {
        y.e(i) = 0;
        for(int j=0; j<x.size(); j++) 
          y.e(i) += A.e(i,j)*x.e(j);
      }
    }

   // Matrix-Vector
  template <class AT, class Type1, class Type2, class Row1, class Col1, class Row2> 
    inline Vector<General,Var,Fixed<1>,AT> operator*(const Matrix<Type1,Row1,Col1,AT> &A, const Vector<Type2,Row2,Fixed<1>,AT> &x) {
      Vector<General,Var,Fixed<1>,AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }
  template <class AT, class Type1, class Type2, class Row, class Col> 
    inline Vector<General,Row,Fixed<1>,AT> operator*(const Matrix<Type1,Row,Col,AT> &A, const Vector<Type2,Row,Fixed<1>,AT> &x) {
      Vector<General,Row,Fixed<1>,AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }
  template <class AT, int M, class Type1, class Type2, class Col1, class Row2> 
    inline Vector<General,Fixed<M>,Fixed<1>,AT> operator*(const Matrix<Type1,Fixed<M>,Col1,AT> &A, const Vector<Type2,Row2,Fixed<1>,AT> &x) {
      Vector<General,Fixed<M>,Fixed<1>,AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }
  template <class AT, int M, class Type1, class Type2> 
    inline Vector<General,Fixed<M>,Fixed<1>,AT> operator*(const Matrix<Type1,Fixed<M>,Fixed<M>,AT> &A, const Vector<Type2,Fixed<M>,Fixed<1>,AT> &x) {
      Vector<General,Fixed<M>,Fixed<1>,AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }
  template <class AT, class Type, class Row1, class Col1, class Row2> 
    inline Vector<Type,Var,Fixed<1>,AT> operator*(const Matrix<Type,Row1,Col1,AT> &A, const Vector<Type,Row2,Fixed<1>,AT> &x) {
      Vector<Type,Var,Fixed<1>,AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }
  template <class AT, class Type, class Row, class Col> 
    inline Vector<Type,Row,Fixed<1>,AT> operator*(const Matrix<Type,Row,Col,AT> &A, const Vector<Type,Row,Fixed<1>,AT> &x) {
      Vector<Type,Row,Fixed<1>,AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }
  template <class AT, int M, class Type, class Col1, class Row2> 
    inline Vector<Type,Fixed<M>,Fixed<1>,AT> operator*(const Matrix<Type,Fixed<M>,Col1,AT> &A, const Vector<Type,Row2,Fixed<1>,AT> &x) {
      Vector<Type,Fixed<M>,Fixed<1>,AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }
  template <class AT, int M, int N, class Type, class Row2> 
    inline Vector<Type,Fixed<M>,Fixed<1>,AT> operator*(const Matrix<Type,Fixed<M>,Fixed<N>,AT> &A, const Vector<Type,Row2,Fixed<1>,AT> &x) {
      Vector<Type,Fixed<M>,Fixed<1>,AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }
  template <class AT, int M, int N, class Type> 
    inline Vector<Type,Fixed<M>,Fixed<1>,AT> operator*(const Matrix<Type,Fixed<M>,Fixed<N>,AT> &A, const Vector<Type,Fixed<N>,Fixed<1>,AT> &x) {
      Vector<Type,Fixed<M>,Fixed<1>,AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }
  template <class AT, int M, class Type> 
    inline Vector<Type,Fixed<M>,Fixed<1>,AT> operator*(const Matrix<Type,Fixed<M>,Fixed<M>,AT> &A, const Vector<Type,Fixed<M>,Fixed<1>,AT> &x) {
      Vector<Type,Fixed<M>,Fixed<1>,AT> y(A.rows(),NONINIT);
      mult(A,x,y);
      return y;
    }

  // RowVector-Matrix
  template <class Type1, class Col1, class Type2, class Row2, class Col2, class Type3, class Col3, class AT> 
    inline void mult(const RowVector<Type1,Fixed<1>,Col1,AT> &x, const Matrix<Type2,Row2,Col2,AT> &A, RowVector<Type3,Fixed<1>,Col3,AT> &y) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x.size() == A.rows());
#endif
      for(int i=0; i<y.size(); i++) {
        y.e(i) = 0;
        for(int j=0; j<x.size(); j++) 
          y.e(i) += x.e(j)*A.e(j,i);
      }
    }

  template <class AT, class Type1, class Type2, class Col2, class Row1, class Col1> 
    inline RowVector<General,Fixed<1>,Var,AT> operator*(const RowVector<Type1,Fixed<1>,Col2,AT> &x, const Matrix<Type2,Row1,Col1,AT> &A) {
      RowVector<General,Fixed<1>,Var,AT> y(A.rows(),NONINIT);
      mult(x,A,y);
      return y;
    }
  template <class AT, class Type1, class Type2, class Row, class Col> 
    inline RowVector<General,Fixed<1>,Col,AT> operator*(const RowVector<Type1,Fixed<1>,Col,AT> &x, const Matrix<Type2,Row,Col,AT> &A) {
      RowVector<General,Fixed<1>,Col,AT> y(A.rows(),NONINIT);
      mult(x,A,y);
      return y;
    }
  template <class AT, int N, class Type1, class Type2, class Col1, class Row2> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator*(const RowVector<Type1,Fixed<1>,Col1,AT> &x, const Matrix<Type2,Row2,Fixed<N>,AT> &A) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> y(A.rows(),NONINIT);
      mult(x,A,y);
      return y;
    }
  template <class AT, int N, class Type1, class Type2> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator*(const RowVector<Type1,Fixed<1>,Fixed<N>,AT> &x, const Matrix<Type2,Fixed<N>,Fixed<N>,AT> &A) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> y(A.rows(),NONINIT);
      mult(x,A,y);
      return y;
    }
  template <class AT, class Type, class Col2, class Row1, class Col1> 
    inline RowVector<General,Fixed<1>,Var,AT> operator*(const RowVector<Type,Fixed<1>,Col2,AT> &x, const Matrix<Type,Row1,Col1,AT> &A) {
      RowVector<General,Fixed<1>,Var,AT> y(A.rows(),NONINIT);
      mult(x,A,y);
      return y;
    }
  template <class AT, class Type, class Row, class Col> 
    inline RowVector<General,Fixed<1>,Col,AT> operator*(const RowVector<Type,Fixed<1>,Col,AT> &x, const Matrix<Type,Row,Col,AT> &A) {
      RowVector<General,Fixed<1>,Col,AT> y(A.rows(),NONINIT);
      mult(x,A,y);
      return y;
    }
  template <class AT, int N, class Type, class Col1, class Row2> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator*(const RowVector<Type,Fixed<1>,Col1,AT> &x, const Matrix<Type,Row2,Fixed<N>,AT> &A) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> y(A.rows(),NONINIT);
      mult(x,A,y);
      return y;
    }
  template <class AT, int N, class Type> 
    inline RowVector<General,Fixed<1>,Fixed<N>,AT> operator*(const RowVector<Type,Fixed<1>,Fixed<N>,AT> &x, const Matrix<Type,Fixed<N>,Fixed<N>,AT> &A) {
      RowVector<General,Fixed<1>,Fixed<N>,AT> y(A.rows(),NONINIT);
      mult(x,A,y);
      return y;
    }

  // Matrix-Matrix
  template <class Type1, class Row1, class Col1, class Type2, class Row2, class Col2, class Type3, class Row3, class Col3, class AT> 
    inline void mult(const Matrix<Type1,Row1,Col1,AT> &A1, const Matrix<Type2,Row2,Col2,AT> &A2, Matrix<Type3,Row3,Col3,AT> &A3) {
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
  template <class Row1, class Type2, class Row2, class Col2, class Type3, class Row3, class Col3, class AT> 
    inline void mult(const Matrix<Symmetric,Row1,Row1,AT> &A1, const Matrix<Type2,Row2,Col2,AT> &A2, Matrix<Type3,Row3,Col3,AT> &A3) {
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

  template <class AT, class Type1, class Type2, class Row1, class Col1, class Row2, class Col2> 
    inline Matrix<General,Var,Var,AT> operator*(const Matrix<Type1,Row1,Col1,AT> &A1, const Matrix<Type2,Row2,Col2,AT> &A2) {
      Matrix<General,Var,Var,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  template <class AT, class Type1, class Type2, class Row, class Col> 
    inline Matrix<General,Row,Col,AT> operator*(const Matrix<Type1,Row,Col,AT> &A1, const Matrix<Type2,Row,Col,AT> &A2) {
      Matrix<General,Row,Col,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  template <class AT, int M, class Type1, class Type2, class Col1, class Row2, class Col2> 
    inline Matrix<General,Fixed<M>,Var,AT> operator*(const Matrix<Type1,Fixed<M>,Col1,AT> &A1, const Matrix<Type2,Row2,Col2,AT> &A2) {
      Matrix<General,Fixed<M>,Var,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  template <class AT, int N, class Type1, class Row1, class Col1, class Type2, class Row2> 
    inline Matrix<General,Var,Fixed<N>,AT> operator*(const Matrix<Type1,Row1,Col1,AT> &A1, const Matrix<Type2,Row2,Fixed<N>,AT> &A2) {
      Matrix<General,Var,Fixed<N>,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  template <class AT, int M, int N, class Type1, class Type2, class Col1, class Row2> 
    inline Matrix<General,Fixed<M>,Fixed<N>,AT> operator*(const Matrix<Type1,Fixed<M>,Col1,AT> &A1, const Matrix<Type2,Row2,Fixed<N>,AT> &A2) {
      Matrix<General,Fixed<M>,Fixed<N>,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  template <class AT, int M, class Type1, class Type2> 
    inline Matrix<General,Fixed<M>,Fixed<M>,AT> operator*(const Matrix<Type1,Fixed<M>,Fixed<M>,AT> &A1, const Matrix<Type2,Fixed<M>,Fixed<M>,AT> &A2) {
      Matrix<General,Fixed<M>,Fixed<M>,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  template <class AT, class Type, class Row1, class Col1, class Row2, class Col2> 
    inline Matrix<Type,Var,Var,AT> operator*(const Matrix<Type,Row1,Col1,AT> &A1, const Matrix<Type,Row2,Col2,AT> &A2) {
      Matrix<Type,Var,Var,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  template <class AT, class Type, class Row, class Col> 
    inline Matrix<Type,Row,Col,AT> operator*(const Matrix<Type,Row,Col,AT> &A1, const Matrix<Type,Row,Col,AT> &A2) {
      Matrix<Type,Row,Col,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  template <class AT, int M, class Type, class Col1, class Row2, class Col2> 
    inline Matrix<Type,Fixed<M>,Var,AT> operator*(const Matrix<Type,Fixed<M>,Col1,AT> &A1, const Matrix<Type,Row2,Col2,AT> &A2) {
      Matrix<Type,Fixed<M>,Var,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  template <class AT, int N, class Type, class Row1, class Col1, class Row2> 
    inline Matrix<Type,Var,Fixed<N>,AT> operator*(const Matrix<Type,Row1,Col1,AT> &A1, const Matrix<Type,Row2,Fixed<N>,AT> &A2) {
      Matrix<Type,Var,Fixed<N>,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  template <class AT, int M, int N, class Type, class Col1, class Row2> 
    inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator*(const Matrix<Type,Fixed<M>,Col1,AT> &A1, const Matrix<Type,Row2,Fixed<N>,AT> &A2) {
      Matrix<Type,Fixed<M>,Fixed<N>,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
 // template <class AT, int M, int K, int N, class Type> 
 //   inline Matrix<Type,Fixed<M>,Fixed<N>,AT> operator*(const Matrix<Type,Fixed<M>,Fixed<K>,AT> &A1, const Matrix<Type,Fixed<K>,Fixed<N>,AT> &A2) {
 //     Matrix<Type,Fixed<M>,Fixed<N>,AT> A3(A1.rows(),A2.cols(),NONINIT);
 //     mult(A1,A2,A3);
 //     return A3;
 //   }
  template <class AT, int M, class Type> 
    inline Matrix<Type,Fixed<M>,Fixed<M>,AT> operator*(const Matrix<Type,Fixed<M>,Fixed<M>,AT> &A1, const Matrix<Type,Fixed<M>,Fixed<M>,AT> &A2) {
      Matrix<Type,Fixed<M>,Fixed<M>,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  
//  SquareMatrix-SquareMatrix
  template <class AT, class Row1, class Row2> 
    inline SquareMatrix<General,Var,Var,AT> operator*(const SquareMatrix<General,Row1,Row1,AT> &A1, const SquareMatrix<General,Row2,Row2,AT> &A2) {
      SquareMatrix<General,Var,Var,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  template <class AT, class Row> 
    inline SquareMatrix<General,Row,Row,AT> operator*(const SquareMatrix<General,Row,Row,AT> &A1, const SquareMatrix<General,Row,Row,AT> &A2) {
      SquareMatrix<General,Row,Row,AT> A3(A1.rows(),A2.cols(),NONINIT);
      mult(A1,A2,A3);
      return A3;
    }
  /////////////////////////////////// end vecscalmult //////////////////////////////

  /*! \brief Vector-scalar multiplication.
   *
   * This function computes the product of a vector 
   * and a scalar.
   * \return The product.
   * */
  template <class Row, class AT>
    Vector<General,Row,Fixed<1>,AT> operator*(const Vector<General,Row,Fixed<1>,AT> &x, const AT& alpha) {

      Vector<General,Row,Fixed<1>,AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
        y.e(i) = x.e(i)*alpha;

      return y;
    }

  /*! \brief Scalar-vector multiplication.
   *
   * \see operator*(const Vector<General,Ref,Fixed<1>,AT>&x,const AT&).
   * */
  template <class Row, class AT>
    Vector<General,Row,Fixed<1>,AT> operator*(const AT& alpha, const Vector<General,Row,Fixed<1>,AT> &x) {

      Vector<General,Row,Fixed<1>,AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
        y.e(i) = x.e(i)*alpha;

      return y;
    }

  template <class Row, class AT>
    inline Vector<General,Row,Fixed<1>,AT> operator*=(const Vector<General,Row,Fixed<1>,AT> &x_, const AT &a) {
      Vector<General,Row,Fixed<1>,AT> &x = const_cast<Vector<General,Row,Fixed<1>,AT> &>(x_);
      for(int i=0; i<x.size(); i++)
	x.e(i) *= a;
      return x;
    }

  template <class Row, class AT>
    inline Vector<General,Row,Fixed<1>,AT> operator/(const Vector<General,Row,Fixed<1>,AT> &x, const AT &a) {
      Vector<General,Row,Fixed<1>,AT> y(x.size(),NONINIT);
      for(int i=0; i<x.size(); i++)
	y.e(i) = x.e(i)/a;
      return y;
    }

  template <class Row, class AT>
    inline Vector<General,Row,Fixed<1>,AT> operator/=(const Vector<General,Row,Fixed<1>,AT> &x_, const AT &a) {
      Vector<General,Row,Fixed<1>,AT> &x = const_cast<Vector<General,Row,Fixed<1>,AT> &>(x_);
      for(int i=0; i<x.size(); i++)
	x.e(i) /= a;
      return x;
    }
  /*! \brief Rowvector-scalar multiplication.
   *
   * This function computes the product of a rowvector 
   * and a scalar.
   * \return The product.
   * */
  template <class Col, class AT>
    RowVector<General,Fixed<1>,Col,AT> operator*(const RowVector<General,Fixed<1>,Col,AT> &x, const AT& alpha) {

      RowVector<General,Fixed<1>,Col,AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
        y.e(i) = x.e(i)*alpha;

      return y;
    }

  /*! \brief Scalar-rowvector multiplication.
   *
   * \see operator*(const RowVector<General<Fixed<1>,Col>, AT>&, const AT&).
   * */
  template <class Col, class AT>
    RowVector<General,Fixed<1>,Col,AT> operator*(const AT &alpha, const RowVector<General,Fixed<1>,Col,AT> &x) {

      RowVector<General,Fixed<1>,Col,AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++) 
        y.e(i) = x.e(i)*alpha;

      return y;
    }

  template <class Col, class AT>
    inline RowVector<General,Fixed<1>,Col,AT> operator*=(const RowVector<General,Fixed<1>,Col,AT> &x_, const AT &a) {
      RowVector<General,Fixed<1>,Col,AT> &x = const_cast<RowVector<General,Fixed<1>,Col,AT> &>(x_);
      for(int i=0; i<x.size(); i++)
	x.e(i) *= a;
      return x;
    }

  template <class Col, class AT>
    inline RowVector<General,Fixed<1>,Col,AT> operator/(const RowVector<General,Fixed<1>,Col,AT> &x, const AT &a) {
      RowVector<General,Fixed<1>,Col,AT> y(x.size(),NONINIT);
      for(int i=0; i<x.size(); i++)
	y.e(i) = x.e(i)/a;
      return y;
    }

  template <class Col, class AT>
    inline RowVector<General,Fixed<1>,Col,AT> operator/=(const RowVector<General,Fixed<1>,Col,AT> &x_, const AT &a) {
      RowVector<General,Fixed<1>,Col,AT> &x = const_cast<RowVector<General,Fixed<1>,Col,AT> &>(x_);
      for(int i=0; i<x.size(); i++)
	x.e(i) /= a;
      return x;
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
    AT operator*(const RowVector<General,Fixed<1>,Col,AT> &x, const Vector<General,Row,Fixed<1>,AT> &y) {

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
    AT scalarProduct(const Vector<General,Row,Fixed<1>,AT> &x, const Vector<General,Row,Fixed<1>,AT> &y) {

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
    Vector<General,Row,Fixed<1>,AT> crossProduct(const Vector<General,Row,Fixed<1>,AT> &x, const Vector<General,Row,Fixed<1>,AT> &y) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(x.size()==3);
      assert(y.size()==3);
#endif

      Vector<General,Row,Fixed<1>,AT> z(3,NONINIT);

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
    double tripleProduct(const Vector<General,Row,Fixed<1>,AT> &a, const Vector<General,Row,Fixed<1>,AT> &x, const Vector<General,Row,Fixed<1>,AT> &y) {

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
  template <class Type, class Row, class Col, class AT>
    Matrix<Type,Row,Col,AT> operator*(const Matrix<Type,Row,Col,AT > &A, const AT &alpha) {

      Matrix<Type,Row,Col,AT> B(A.rows(),A.cols(),NONINIT);

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++) 
          B.e(i,j) = A.e(i,j)*alpha;

      return B;
    }

  template <class Type, class Row, class Col, class AT>
    Matrix<Type,Row,Col,AT> operator*(const AT &alpha, const Matrix<Type,Row,Col,AT > &A) {

      Matrix<Type,Row,Col,AT> B(A.rows(),A.cols(),NONINIT);

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++) 
          B.e(i,j) = A.e(i,j)*alpha;

      return B;
    }
  
  template <class Row, class AT>
   inline Matrix<Symmetric,Row,Row,AT> operator*(const AT &a, const Matrix<Symmetric,Row,Row,AT> &A) {
     Matrix<Symmetric,Row,Row,AT> B(A.size(),NONINIT);
     for(int i=0; i<A.size(); i++)
       for(int j=i; j<A.size(); j++)
         B.ej(i,j) = a*A.ej(i,j);
     return B;
   }

  template <class Row, class AT>
    inline Matrix<Symmetric,Row,Row,AT> operator*(const Matrix<Symmetric,Row,Row,AT> &A, const AT &a) {
      Matrix<Symmetric,Row,Row,AT> B(A.size(),NONINIT);
      for(int i=0; i<A.size(); i++)
        for(int j=i; j<A.size(); j++)
          B.ej(i,j) = a*A.ej(i,j);
      return B;
    }

  /////////////////////////////////// end matscalmult //////////////////////////////

  /////////////////////////////////// negation //////////////////////////////

  /*! \brief Negation.
   *
   * This function computes the negation of a vector.
   * \return The negation.
   * */
  template <class Type, class Row, class AT>
    Vector<Type,Row,Fixed<1>,AT> operator-(const Vector<Type,Row,Fixed<1>,AT> &x) {

      Vector<Type,Row,Fixed<1>,AT> y(x.size(),NONINIT);

      for(int i=0; i<x.size(); i++)
        y.e(i)=-x.e(i);

      return y;
    }

  /*! \brief Negation.
   *
   * This function computes the negation of a rowvector.
   * \return The negation.
   * */
  template <class Type, class Col, class AT>
    RowVector<Type,Fixed<1>,Col,AT> operator-(const RowVector<Type,Fixed<1>,Col,AT> &a) {

      RowVector<Type,Fixed<1>,Col,AT> c(a.size(),NONINIT);

      for(int i=0; i<a.size(); i++)
        c.e(i)=-a.e(i);

      return c;
    }

  /*! \brief Negation.
   *
   * This function computes the negation of a vector.
   * \return The negation.
   * */
  template <class Type, class Row, class AT>
    SquareMatrix<Type,Row,Row,AT> operator-(const SquareMatrix<Type,Row,Row,AT> &A) {

      SquareMatrix<Type,Row,Row,AT> B(A.size(),NONINIT);

      for(int i=0; i<A.rows(); i++)
        for(int j=0; j<A.cols(); j++)
          B.e(i,j)=-A.e(i,j);

      return B;
    }

  /*! \brief Negation.
   *
   * This function computes the negation of a vector.
   * \return The negation.
   * */
  template <class Type, class Row, class Col, class AT>
    Matrix<Type,Row,Col,AT> operator-(const Matrix<Type,Row,Col,AT> &A) {

      Matrix<Type,Row,Col,AT> B(A.rows(),A.cols(),NONINIT);

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
    RowVector<General,Fixed<1>,Ref,AT> trans(const Vector<General,Ref,Fixed<1>,AT> &x) {

      return RowVector<General,Fixed<1>,Ref,AT>(x.m,x.lda,x.tp?false:true,x.memory,x.ele).copy();
    }

  /*! \brief Transpose of a rowvector.
   *
   * This function computes the transpose of a rowvector.
   * The result is a vector.
   * \param x A vector.
   * \return The transpose.
   * */
  template <class AT>
    Vector<General,Ref,Fixed<1>,AT> trans(const RowVector<General,Fixed<1>,Ref,AT> &x) {

      return Vector<General,Ref,Fixed<1>,AT>(x.n,x.lda,x.tp?false:true,x.memory,x.ele).copy();
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
    Matrix<General,Ref,Ref,AT> trans(const Matrix<General,Ref,Ref,AT> &A) {

      return Matrix<General,Ref,Ref,AT>(A.n,A.m,A.lda,A.tp?false:true,A.memory,A.ele).copy();
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
    SquareMatrix<General,Ref,Ref,AT> trans(const SquareMatrix<General,Ref,Ref,AT> &A) {

      return SquareMatrix<General,Ref,Ref,AT>(A.n, A.lda, A.tp?false:true, A.memory, A.ele).copy();
    }

  template <class Type, class Row, class Col, class AT>
    inline Matrix<Type,Col,Row,AT> trans(const Matrix<Type,Row,Col,AT> &A) {
      Matrix<Type,Col,Row,AT> B(NONINIT);
      for(int i=0; i<B.rows(); i++)
        for(int j=0; j<B.cols(); j++)
          B.e(i,j) = A.e(j,i);
      return B;
    }

  template <class Row, class AT>
    inline SquareMatrix<General,Row,Row,AT> trans(const SquareMatrix<General,Row,Row,AT> &A) {
      SquareMatrix<General,Row,Row,AT> B(NONINIT);
      for(int i=0; i<B.size(); i++)
        for(int j=0; j<B.size(); j++)
          B.e(i,j) = A.e(j,i);
      return B;
    }
  /////////////////////////////////// end transpose //////////////////////////////

  template <class Row, class AT>
    inline AT nrm2(const Vector<General,Row,Fixed<1>,AT> &x) {
      AT c = 0;
      for(int i=0; i<x.size(); i++)
        c += pow(x.e(i),2);
      return sqrt(c);
    }

  template <class Col, class AT>
    inline AT nrm2(const RowVector<General,Fixed<1>,Col,AT> &x) {
      AT c = 0;
      for(int i=0; i<x.size(); i++)
        c += pow(x.e(i),2);
      return sqrt(c);
    }

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
  template <class Type, class Row, class AT>
    SquareMatrix<Type,Row,Row,AT> tilde(const Vector<Type,Row,Fixed<1>,AT> &x) {

      SquareMatrix<Type,Row,Row,AT> B(x.size(),NONINIT);

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

  /*! \brief Matrix-matrix comparison.
   *
   * This function compares two matrices
   * \return True, if the matrices are identical, false otherwise.
   * */
  template <class AT, class Type1, class Row1, class Col1, class Type2, class Row2, class Col2> 
    bool operator==(const Matrix<Type1,Row1,Col1,AT> &A, const Matrix<Type2,Row2,Col2,AT > &B) {

      if(&A == &B)
        return true;

      if(A.rows() != B.rows() || A.cols() != B.cols())
        return false;

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++)
          if(A(i,j)!=B(i,j))
            return false;

      return true;
    }

  /*! \brief Matrix-matrix comparison.
   *
   * This function compares two matrices.
   * \return True, if the matrices are different, false otherwise.
   * */
  template <class AT, class Type1, class Row1, class Col1, class Type2, class Row2, class Col2> 
    bool operator!=(const Matrix<Type1,Row1,Col1,AT> &A, const Matrix<Type2,Row2,Col2,AT > &B) {

      if(A.rows() != B.rows() || A.cols() != B.cols())
        return true;

      for(int i=0; i<A.rows(); i++) 
        for(int j=0; j<A.cols(); j++)
          if(A(i,j)!=B(i,j))
            return true;

      return false;
    }


  /*! \brief Maximum value.
   *
   * This function computes the maximum value of a vector.
   * \return The maximum value.
   * */
  template <class Type, class Row, class AT>
    AT max(const Vector<Type,Row,Fixed<1>,AT> &x) {

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
  template <class Type, class Row, class AT>
    int maxIndex(const Vector<Type,Row,Fixed<1>,AT> &x) {

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
  template <class Type, class Row, class AT>
    AT min(const Vector<Type,Row,Fixed<1>,AT> &x) {

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
  template <class Type, class Row, class AT>
    int minIndex(const Vector<Type,Row,Fixed<1>,AT> &x) {

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
    Matrix<General,Ref,Ref,AT> bubbleSort(const Matrix<General,Ref,Ref,AT> &A_, int PivotCol) {

#ifndef FMATVEC_NO_SIZE_CHECK
      assert(A_.rows()>0);
      assert(A_.cols()>0);
      assert(A_.cols()> PivotCol);
#endif
      Matrix<General,Ref,Ref,AT> A = A_.copy();
      int i,j,N;
      RowVector<General,Fixed<1>,Ref,AT> tmp;
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
    void quicksortmedian_intern(Matrix<General,Ref,Ref,AT> &A, int PivotCol, RowVector<General,Fixed<1>,Ref,AT> &tmp, int l, int r){

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
    Matrix<General,Ref,Ref,AT> quickSortMedian(const Matrix<General,Ref,Ref,AT> &A_, int PivotCol) {
      Matrix<General,Ref,Ref,AT> A = A_.copy();
      int N = A.rows();
      RowVector<General,Fixed<1>,Ref,AT> tmp;
      quicksortmedian_intern(A, PivotCol,tmp, 0, N-1);
      return A;
    }

  /*! \brief Count nonzero elements.
   *
   * This function counts the nonzero elements of a matrix. ALL diagonal
   * elements are treated as NONZERO!!! (See the storage format for sparse matrix)
   * \return The number of nonzero or diagonal elemets.
   * */
  template <class AT> int countElements(const SquareMatrix<General,Ref,Ref,AT> &A) { 
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

  template <class Row, class Col, class AT>
    inline Matrix<Symmetric,Col,Col,AT> JTJ(const Matrix<General,Row,Col,AT> &A) { 
      Matrix<Symmetric,Col,Col,AT> S(A.cols(),NONINIT);
      for(int i=0; i<A.cols(); i++) {
        for(int k=i; k<A.cols(); k++) {
          S.ej(i,k) = 0;
          for(int j=0; j<A.rows(); j++) 
            S.ej(i,k) += A.e(j,i)*A.e(j,k);
        }
      }
      return S;
    }

  template <class Row, class Col, class AT>
    inline Matrix<Symmetric,Col,Col,AT> JTMJ(const Matrix<Symmetric,Row,Row,AT> &B, const Matrix<General,Row,Col,AT> &A) {

      Matrix<Symmetric,Col,Col,AT> S(A.cols(),NONINIT);
      Matrix<General,Row,Col,AT> C = B*A;

      for(int i=0; i<A.cols(); i++) {
        for(int k=i; k<A.cols(); k++) {
          S.ej(i,k) = 0;
          for(int j=0; j<A.rows(); j++) 
            S.ej(i,k) += A.e(j,i)*C.e(j,k);
        }
      }
      return S;
    }

  template <class Row, class Col, class AT>
    inline Matrix<Symmetric,Row,Row,AT> JMJT(const Matrix<General,Row,Col,AT> &A, const Matrix<Symmetric,Col,Col,AT> &B) {

      Matrix<Symmetric,Row,Row,AT> S(A.rows(),NONINIT);
      Matrix<General,Row,Col,AT> C = A*B;

      for(int i=0; i<S.size(); i++) {
        for(int k=i; k<S.size(); k++) {
          S.ej(i,k) = 0;
          for(int j=0; j<A.cols(); j++) 
            S.ej(i,k) += C.e(k,j)*A.e(i,j);
        }
      }
      return S;
    }


}

#endif

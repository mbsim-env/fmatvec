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

#include <cmath>

#include "square_matrix.h"
#include "vector.h"
#include "row_vector.h"
#include "fixed_square_matrix.h"
#include "fixed_vector.h"
#include "fixed_row_vector.h"
#include "fixed_symmetric_matrix.h"
#include "var_square_matrix.h"
#include "var_symmetric_matrix.h"
#include "var_vector.h"
#include "var_row_vector.h"
#include "var_fixed_general_matrix.h"
#include "fixed_var_general_matrix.h"
#include <cmath>

namespace fmatvec {

/////////////////////////////////// vecvecadd //////////////////////////////

  // Vector-Vector
  template <class Row1, class Row2, class Row3, class AT1, class AT2, class AT3>
  inline void add(const Vector<Row1, AT1> &a1, const Vector<Row2, AT2> &a2, Vector<Row3, AT3> &a3) {
    assert(a1.size() == a2.size());
    for (int i = 0; i < a1.size(); i++)
      a3.e(i) = a1.e(i) + a2.e(i);
  }

  template <class Row1, class Row2, class AT>
  inline void add(Vector<Row1, AT> &a1, const Vector<Row2, AT> &a2) {
    assert(a1.size() == a2.size());
    for (int i = 0; i < a1.size(); i++)
      a1.e(i) += a2.e(i);
  }

  template <class Row1, class Row2, class Row3, class AT1, class AT2>
  inline void sub(const Vector<Row1, AT1> &a1, const Vector<Row2, AT2> &a2, Vector<Row3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &a3) {
    assert(a1.size() == a2.size());
    for (int i = 0; i < a1.size(); i++)
      a3.e(i) = a1.e(i) - a2.e(i);
  }

  template <class Row1, class Row2, class AT1, class AT2>
  inline void sub(Vector<Row1, AT1> &a1, const Vector<Row2, AT2> &a2) {
    assert(a1.size() == a2.size());
    for (int i = 0; i < a1.size(); i++)
      a1.e(i) -= a2.e(i);
  }

  // Addition
  template <class AT1, class AT2, class Row1, class Row2>
  inline Vector<Var, typename OperatorResult<AT1, AT2>::Type> operator+(const Vector<Row1, AT1> &a, const Vector<Row2, AT2> &b) {
    Vector<Var, typename OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    add(a, b, c);
    return c;
  }
  template <class AT1, class AT2, class Row>
  inline Vector<Row, typename OperatorResult<AT1, AT2>::Type> operator+(const Vector<Row, AT1> &a, const Vector<Row, AT2> &b) {
    Vector<Row, typename OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    add(a, b, c);
    return c;
  }
  template <class AT1, class AT2, int M, class Row>
  inline Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> operator+(const Vector<Fixed<M>, AT1> &a, const Vector<Row, AT2> &b) {
    Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    add(a, b, c);
    return c;
  }
  template <class AT1, class AT2, int M>
  inline Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> operator+(const Vector<Fixed<M>, AT1> &a, const Vector<Fixed<M>, AT2> &b) {
    Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    add(a, b, c);
    return c;
  }
  template <class AT, class Row1, class Row2>
  inline Vector<Row1, AT> operator+=(const Vector<Row1, AT> &a_, const Vector<Row2, AT> &b) {
    auto &a = const_cast<Vector<Row1, AT> &>(a_);
    add(a, b);
    return a;
  }

///  // Subtraction
  template <class AT1, class AT2, class Row1, class Row2>
  inline Vector<Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Vector<Row1, AT1> &a, const Vector<Row2, AT2> &b) {
    Vector<Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    sub(a, b, c);
    return c;
  }
  template <class AT1, class AT2, class Row>
  inline Vector<Row, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Vector<Row, AT1> &a, const Vector<Row, AT2> &b) {
    Vector<Row, typename fmatvec::OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    sub(a, b, c);
    return c;
  }
  template <class AT1, class AT2, int M, class Row2>
  inline Vector<Fixed<M>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Vector<Fixed<M>, AT1> &a, const Vector<Row2, AT2> &b) {
    Vector<Fixed<M>, typename fmatvec::OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    sub(a, b, c);
    return c;
  }
  template <class AT1, class AT2, int M>
  inline Vector<Fixed<M>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Vector<Fixed<M>, AT1> &a, const Vector<Fixed<M>, AT2> &b) {
    Vector<Fixed<M>, typename fmatvec::OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    sub(a, b, c);
    return c;
  }
  template <class AT, class Row1, class Row2>
  inline Vector<Row1, AT> operator-=(const Vector<Row1, AT> &a_, const Vector<Row2, AT> &b) {
    auto &a = const_cast<Vector<Row1, AT> &>(a_);
    sub(a, b);
    return a;
  }

  // RowVector-Vector
  template <class Col1, class Col2, class Col3, class AT>
  inline void add(const RowVector<Col1, AT> &a1, const RowVector<Col2, AT> &a2, RowVector<Col3, AT> &a3) {
    assert(a1.size() == a2.size());
    for (int i = 0; i < a1.size(); i++)
      a3.e(i) = a1.e(i) + a2.e(i);
  }

  template <class Col1, class Col2, class AT>
  inline void add(RowVector<Col1, AT> &a1, const RowVector<Col2, AT> &a2) {
    assert(a1.size() == a2.size());
    for (int i = 0; i < a1.size(); i++)
      a1.e(i) += a2.e(i);
  }

  template <class Col1, class Col2, class Col3, class AT1, class AT2>
  inline void sub(const RowVector<Col1, AT1> &a1, const RowVector<Col2, AT2> &a2, RowVector<Col3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &a3) {
    assert(a1.size() == a2.size());
    for (int i = 0; i < a1.size(); i++)
      a3.e(i) = a1.e(i) - a2.e(i);
  }

  template <class Col1, class Col2, class AT1, class AT2>
  inline void sub(RowVector<Col1, AT1> &a1, const RowVector<Col2, AT2> &a2) {
    assert(a1.size() == a2.size());
    for (int i = 0; i < a1.size(); i++)
      a1.e(i) -= a2.e(i);
  }

  // Addition
  template <class AT1, class AT2, class Col1, class Col2>
  inline RowVector<Var, typename OperatorResult<AT1, AT2>::Type> operator+(const RowVector<Col1, AT1> &a, const RowVector<Col2, AT2> &b) {
    RowVector<Var, typename OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    add(a, b, c);
    return c;
  }
  template <class AT1, class AT2, class Col>
  inline RowVector<Col, typename OperatorResult<AT1, AT2>::Type> operator+(const RowVector<Col, AT1> &a, const RowVector<Col, AT2> &b) {
    RowVector<Col, typename OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    add(a, b, c);
    return c;
  }
  template <class AT1, class AT2, int N, class Col2>
  inline RowVector<Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const RowVector<Fixed<N>, AT1> &a, const RowVector<Col2, AT2> &b) {
    RowVector<Fixed<N>, typename OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    add(a, b, c);
    return c;
  }
  template <class AT1, class AT2, int N>
  inline RowVector<Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const RowVector<Fixed<N>, AT1> &a, const RowVector<Fixed<N>, AT2> &b) {
    RowVector<Fixed<N>, typename OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    add(a, b, c);
    return c;
  }
  template <class AT, class Col1, class Col2>
  inline RowVector<Col1, AT> operator+=(const RowVector<Col1, AT> &a_, const RowVector<Col2, AT> &b) {
    auto &a = const_cast<RowVector<Col1, AT> &>(a_);
    add(a, b);
    return a;
  }

///  // Subtraction
  template <class AT1, class AT2, class Col1, class Col2>
  inline RowVector<Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const RowVector<Col1, AT1> &a, const RowVector<Col2, AT2> &b) {
    RowVector<Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    sub(a, b, c);
    return c;
  }
  template <class AT1, class AT2, class Col>
  inline RowVector<Col, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const RowVector<Col, AT1> &a, const RowVector<Col, AT2> &b) {
    RowVector<Col, typename fmatvec::OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    sub(a, b, c);
    return c;
  }
  template <class AT1, class AT2, int N, class Col2>
  inline RowVector<Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const RowVector<Fixed<N>, AT1> &a, const RowVector<Col2, AT2> &b) {
    RowVector<Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    sub(a, b, c);
    return c;
  }
  template <class AT1, class AT2, int N>
  inline RowVector<Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const RowVector<Fixed<N>, AT1> &a, const RowVector<Fixed<N>, AT2> &b) {
    RowVector<Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> c(a.size(), NONINIT);
    sub(a, b, c);
    return c;
  }
  template <class AT, class Col1, class Col2>
  inline RowVector<Col1, AT> operator-=(const RowVector<Col1, AT> &a_, const RowVector<Col2, AT> &b) {
    auto &a = const_cast<RowVector<Col1, AT> &>(a_);
    sub(a, b);
    return a;
  }

  // Matrix-Matrix
  template <class Type1, class Row1, class Col1, class Type2, class Row2, class Col2, class Type3, class Row3, class Col3, class AT1, class AT2, class AT3>
  inline void add(const Matrix<Type1, Row1, Col1, AT1> &A1, const Matrix<Type2, Row2, Col2, AT2> &A2, Matrix<Type3, Row3, Col3, AT3> &A3) {
    assert(A1.rows() == A2.rows());
    assert(A1.cols() == A2.cols());
    for (int i = 0; i < A1.rows(); i++)
      for (int j = 0; j < A2.cols(); j++)
        A3.e(i, j) = A1.e(i, j) + A2.e(i, j);
  }

  template <class Row1, class Row2, class Row3, class AT>
  inline void add(const Matrix<Symmetric, Row1, Row1, AT> &A1, const Matrix<Symmetric, Row2, Row2, AT> &A2, Matrix<Symmetric, Row3, Row3, AT> &A3) {
    assert(A1.size() == A2.size());
    for (int i = 0; i < A1.size(); i++)
      for (int j = i; j < A2.size(); j++)
        A3.ej(i, j) = A1.ej(i, j) + A2.ej(i, j);
  }

  template <class Row1, class Row2, class Row3, class AT>
  inline void add(const Matrix<Diagonal, Row1, Row1, AT> &A1, const Matrix<Diagonal, Row2, Row2, AT> &A2, Matrix<Diagonal, Row3, Row3, AT> &A3) {
    assert(A1.size() == A2.size());
    for (int i = 0; i < A3.size(); i++)
      A3.e(i) = A1.e(i) + A2.e(i);
  }

  template <class Type1, class Row1, class Col1, class Type2, class Row2, class Col2, class AT>
  inline void add(Matrix<Type1, Row1, Col1, AT> &A1, const Matrix<Type2, Row2, Col2, AT> &A2) {
    assert(A1.rows() == A2.rows());
    assert(A1.cols() == A2.cols());
    for (int i = 0; i < A1.rows(); i++)
      for (int j = 0; j < A2.cols(); j++)
        A1.e(i, j) += A2.e(i, j);
  }

  template <class Row1, class Row2, class AT>
  inline void add(Matrix<Symmetric, Row1, Row1, AT> &A1, const Matrix<Symmetric, Row2, Row2, AT> &A2) {
    assert(A1.size() == A2.size());
    for (int i = 0; i < A1.size(); i++)
      for (int j = i; j < A2.size(); j++)
        A1.ej(i, j) += A2.ej(i, j);
  }

  template <class Type1, class Row1, class Col1, class Type2, class Row2, class Col2, class Type3, class Row3, class Col3, class AT1, class AT2>
  inline void sub(const Matrix<Type1, Row1, Col1, AT1> &A1, const Matrix<Type2, Row2, Col2, AT2> &A2, Matrix<Type3, Row3, Col3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &A3) {
    assert(A1.rows() == A2.rows());
    assert(A1.cols() == A2.cols());
    for (int i = 0; i < A1.rows(); i++)
      for (int j = 0; j < A2.cols(); j++)
        A3.e(i, j) = A1.e(i, j) - A2.e(i, j);
  }

  template <class Row1, class Row2, class Row3, class AT1, class AT2>
  inline void sub(const Matrix<Symmetric, Row1, Row1, AT1> &A1, const Matrix<Symmetric, Row2, Row2, AT2> &A2, Matrix<Symmetric, Row3, Row3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &A3) {
    assert(A1.size() == A2.size());
    for (int i = 0; i < A1.size(); i++)
      for (int j = i; j < A2.size(); j++)
        A3.ej(i, j) = A1.ej(i, j) - A2.ej(i, j);
  }

  template <class Row1, class Row2, class Row3, class AT1, class AT2>
  inline void sub(const Matrix<Diagonal, Row1, Row1, AT1> &A1, const Matrix<Diagonal, Row2, Row2, AT2> &A2, Matrix<Diagonal, Row3, Row3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &A3) {
    assert(A1.size() == A2.size());
    for (int i = 0; i < A3.size(); i++)
      A3.e(i) = A1.e(i) - A2.e(i);
  }

  template <class Type1, class Row1, class Col1, class Type2, class Row2, class Col2, class AT1, class AT2>
  inline void sub(Matrix<Type1, Row1, Col1, AT1> &A1, const Matrix<Type2, Row2, Col2, AT2> &A2) {
    assert(A1.rows() == A2.rows());
    assert(A1.cols() == A2.cols());
    for (int i = 0; i < A1.rows(); i++)
      for (int j = 0; j < A2.cols(); j++)
        A1.e(i, j) -= A2.e(i, j);
  }

  template <class Row1, class Row2, class AT1, class AT2>
  inline void sub(Matrix<Symmetric, Row1, Row1, AT1> &A1, const Matrix<Symmetric, Row2, Row2, AT2> &A2) {
    assert(A1.size() == A2.size());
    for (int i = 0; i < A1.size(); i++)
      for (int j = i; j < A2.size(); j++)
        A1.ej(i, j) -= A2.ej(i, j);
  }

  // Addition
  // Type1 Type2
  template <class AT1, class AT2, class Type1, class Type2, class Row1, class Col1, class Row2, class Col2>
  inline Matrix<General, Var, Var, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type1, Row1, Col1, AT1> &A, const Matrix<Type2, Row2, Col2, AT2> &B) {
    Matrix<General, Var, Var, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, class Type1, class Type2, class Row, class Col>
  inline Matrix<General, Row, Col, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type1, Row, Col, AT1> &A, const Matrix<Type2, Row, Col, AT2> &B) {
    Matrix<General, Row, Col, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, class Type1, class Type2, class Row, class Col>
  inline Matrix<General, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type1, Fixed<M>, Var, AT1> &A, const Matrix<Type2, Row, Col, AT2> &B) {
    Matrix<General, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type1, Fixed<M>, Var, AT1> &A, const Matrix<Type2, Fixed<M>, Var, AT2> &B) {
    Matrix<General, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int N, class Type1, class Type2, class Row, class Col>
  inline Matrix<General, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type1, Var, Fixed<N>, AT1> &A, const Matrix<Type2, Row, Col, AT2> &B) {
    Matrix<General, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int N, class Type1, class Type2>
  inline Matrix<General, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type1, Var, Fixed<N>, AT1> &A, const Matrix<Type2, Var, Fixed<N>, AT2> &B) {
    Matrix<General, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2, class Row, class Col>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type1, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type2, Row, Col, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type1, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type2, Fixed<M>, Var, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type1, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type2, Var, Fixed<N>, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type1, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type2, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2, class Row, class Col>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type2, Row, Col, AT1> &A, const Matrix<Type1, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type2, Fixed<M>, Var, AT1> &A, const Matrix<Type1, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type2, Var, Fixed<N>, AT1> &A, const Matrix<Type1, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  // Type
  template <class AT1, class AT2, class Type, class Row1, class Col1, class Row2, class Col2>
  inline Matrix<Type, Var, Var, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Row1, Col1, AT1> &A, const Matrix<Type, Row2, Col2, AT2> &B) {
    Matrix<Type, Var, Var, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, class Type, class Row, class Col>
  inline Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Row, Col, AT1> &A, const Matrix<Type, Row, Col, AT2> &B) {
    Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, class Type, class Row, class Col>
  inline Matrix<Type, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Fixed<M>, Var, AT1> &A, const Matrix<Type, Row, Col, AT2> &B) {
    Matrix<Type, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, class Type>
  inline Matrix<Type, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Fixed<M>, Var, AT1> &A, const Matrix<Type, Fixed<M>, Var, AT2> &B) {
    Matrix<Type, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int N, class Type, class Row, class Col>
  inline Matrix<Type, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Var, Fixed<N>, AT1> &A, const Matrix<Type, Row, Col, AT2> &B) {
    Matrix<Type, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int N, class Type>
  inline Matrix<Type, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Var, Fixed<N>, AT1> &A, const Matrix<Type, Var, Fixed<N>, AT2> &B) {
    Matrix<Type, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type, class Row, class Col>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type, Row, Col, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type, Fixed<M>, Var, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type, Var, Fixed<N>, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type, class Row, class Col>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Row, Col, AT1> &A, const Matrix<Type, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Fixed<M>, Var, AT1> &A, const Matrix<Type, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator+(const Matrix<Type, Var, Fixed<N>, AT1> &A, const Matrix<Type, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    add(A, B, C);
    return C;
  }
  template <class AT, class Type1, class Row1, class Col1, class Type2, class Row2, class Col2>
  inline Matrix<Type1, Row1, Col1, AT>& operator+=(const Matrix<Type1, Row1, Col1, AT> &A_, const Matrix<Type2, Row2, Col2, AT> &B) {
    auto &A = const_cast<Matrix<Type1, Row1, Col1, AT> &>(A_);
    add(A, B);
    return A;
  }

  // Subtraction
  // Type1 Type2
  template <class AT1, class AT2, class Type1, class Type2, class Row1, class Col1, class Row2, class Col2>
  inline Matrix<General, Var, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type1, Row1, Col1, AT1> &A, const Matrix<Type2, Row2, Col2, AT2> &B) {
    Matrix<General, Var, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, class Type1, class Type2, class Row, class Col>
  inline Matrix<General, Row, Col, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type1, Row, Col, AT1> &A, const Matrix<Type2, Row, Col, AT2> &B) {
    Matrix<General, Row, Col, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, class Type1, class Type2, class Row, class Col>
  inline Matrix<General, Fixed<M>, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type1, Fixed<M>, Var, AT1> &A, const Matrix<Type2, Row, Col, AT2> &B) {
    Matrix<General, Fixed<M>, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type1, Fixed<M>, Var, AT1> &A, const Matrix<Type2, Fixed<M>, Var, AT2> &B) {
    Matrix<General, Fixed<M>, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int N, class Type1, class Type2, class Row, class Col>
  inline Matrix<General, Var, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type1, Var, Fixed<N>, AT1> &A, const Matrix<Type2, Row, Col, AT2> &B) {
    Matrix<General, Var, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int N, class Type1, class Type2>
  inline Matrix<General, Var, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type1, Var, Fixed<N>, AT1> &A, const Matrix<Type2, Var, Fixed<N>, AT2> &B) {
    Matrix<General, Var, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2, class Row, class Col>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type1, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type2, Row, Col, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type1, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type2, Fixed<M>, Var, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type1, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type2, Var, Fixed<N>, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type1, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type2, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2, class Row, class Col>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type2, Row, Col, AT1> &A, const Matrix<Type1, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type2, Fixed<M>, Var, AT1> &A, const Matrix<Type1, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type2, Var, Fixed<N>, AT1> &A, const Matrix<Type1, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<General, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  // Type
  template <class AT1, class AT2, class Type, class Row1, class Col1, class Row2, class Col2>
  inline Matrix<Type, Var, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Row1, Col1, AT1> &A, const Matrix<Type, Row2, Col2, AT2> &B) {
    Matrix<Type, Var, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, class Type, class Row, class Col>
  inline Matrix<Type, Row, Col, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Row, Col, AT1> &A, const Matrix<Type, Row, Col, AT2> &B) {
    Matrix<Type, Row, Col, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, class Type, class Row, class Col>
  inline Matrix<Type, Fixed<M>, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Fixed<M>, Var, AT1> &A, const Matrix<Type, Row, Col, AT2> &B) {
    Matrix<Type, Fixed<M>, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, class Type>
  inline Matrix<Type, Fixed<M>, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Fixed<M>, Var, AT1> &A, const Matrix<Type, Fixed<M>, Var, AT2> &B) {
    Matrix<Type, Fixed<M>, Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int N, class Type, class Row, class Col>
  inline Matrix<Type, Var, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Var, Fixed<N>, AT1> &A, const Matrix<Type, Row, Col, AT2> &B) {
    Matrix<Type, Var, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int N, class Type>
  inline Matrix<Type, Var, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Var, Fixed<N>, AT1> &A, const Matrix<Type, Var, Fixed<N>, AT2> &B) {
    Matrix<Type, Var, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type, class Row, class Col>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type, Row, Col, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type, Fixed<M>, Var, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type, Var, Fixed<N>, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Fixed<M>, Fixed<N>, AT1> &A, const Matrix<Type, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type, class Row, class Col>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Row, Col, AT1> &A, const Matrix<Type, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Fixed<M>, Var, AT1> &A, const Matrix<Type, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT1, class AT2, int M, int N, class Type>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const Matrix<Type, Var, Fixed<N>, AT1> &A, const Matrix<Type, Fixed<M>, Fixed<N>, AT2> &B) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename fmatvec::OperatorResult<AT1, AT2>::Type> C(A.rows(), A.cols(), NONINIT);
    sub(A, B, C);
    return C;
  }
  template <class AT, class Type1, class Row1, class Col1, class Type2, class Row2, class Col2>
  inline Matrix<Type1, Row1, Col1, AT>& operator-=(const Matrix<Type1, Row1, Col1, AT> &A_, const Matrix<Type2, Row2, Col2, AT> &B) {
    auto &A = const_cast<Matrix<Type1, Row1, Col1, AT> &>(A_);
    sub(A, B);
    return A;
  }

//  SquareMatrix-SquareMatrix
  template <class AT1, class AT2, class Row1, class Row2>
  inline SquareMatrix<Var, typename OperatorResult<AT1, AT2>::Type> operator+(const SquareMatrix<Row1, AT1> &A1, const SquareMatrix<Row2, AT2> &A2) {
    SquareMatrix<Var, typename OperatorResult<AT1, AT2>::Type> A3(A1.size(), NONINIT);
    add(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, class Row>
  inline SquareMatrix<Row, typename OperatorResult<AT1, AT2>::Type> operator+(const SquareMatrix<Row, AT1> &A1, const SquareMatrix<Row, AT2> &A2) {
    SquareMatrix<Row, typename OperatorResult<AT1, AT2>::Type> A3(A1.size(), NONINIT);
    add(A1, A2, A3);
    return A3;
  }

//  SquareMatrix-SquareMatrix
  template <class AT1, class AT2, class Row1, class Row2>
  inline SquareMatrix<Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const SquareMatrix<Row1, AT1> &A1, const SquareMatrix<Row2, AT2> &A2) {
    SquareMatrix<Var, typename fmatvec::OperatorResult<AT1, AT2>::Type> A3(A1.size(), NONINIT);
    sub(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, class Row>
  inline SquareMatrix<Row, typename fmatvec::OperatorResult<AT1, AT2>::Type> operator-(const SquareMatrix<Row, AT1> &A1, const SquareMatrix<Row, AT2> &A2) {
    SquareMatrix<Row, typename fmatvec::OperatorResult<AT1, AT2>::Type> A3(A1.size(), NONINIT);
    sub(A1, A2, A3);
    return A3;
  }

  //////      
  //////      /////////////////////////////////// end vecvecadd //////////////////////////////
  //////      
  //////      /////////////////////////////////// matmatmult //////////////////////////////
  //////      
  template <class Type1, class Row1, class Col1, class Row2, class Row3, class AT1, class AT2>
  inline void mult(const Matrix<Type1, Row1, Col1, AT1> &A, const Vector<Row2, AT2> &x, Vector<Row3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &y) {
    assert(A.cols() == x.size());
    for (int i = 0; i < y.size(); i++) {
      y.e(i) = 0;
      for (int j = 0; j < x.size(); j++)
        y.e(i) += A.e(i, j) * x.e(j);
    }
  }

  template <class Row1, class Row2, class Row3, class AT1, class AT2>
  inline void mult(const Matrix<Symmetric, Row1, Row1, AT1> &A, const Vector<Row2, AT2> &x, Vector<Row3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &y) {
    assert(A.cols() == x.size());
    for (int i = 0; i < y.size(); i++) {
      y.e(i) = 0;
      for (int j = 0; j < i; j++)
        y.e(i) += A.ei(i, j) * x.e(j);
      for (int j = i; j < x.size(); j++)
        y.e(i) += A.ej(i, j) * x.e(j);
    }
  }

  // Matrix-Vector
  template <class AT1, class AT2, class Type1, class Row1, class Col1, class Row2>
  inline Vector<Var, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type1, Row1, Col1, AT1> &A, const Vector<Row2, AT2> &x) {
    Vector<Var, typename OperatorResult<AT1, AT2>::Type> y(A.rows(), NONINIT);
    mult(A, x, y);
    return y;
  }
  template <class AT1, class AT2, class Type1, class Row, class Col>
  inline Vector<Row, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type1, Row, Col, AT1> &A, const Vector<Row, AT2> &x) {
    Vector<Row, typename OperatorResult<AT1, AT2>::Type> y(A.rows(), NONINIT);
    mult(A, x, y);
    return y;
  }
  template <class AT1, class AT2, int M, class Type1, class Col1, class Row2>
  inline Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type1, Fixed<M>, Col1, AT1> &A, const Vector<Row2, AT2> &x) {
    Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> y(A.rows(), NONINIT);
    mult(A, x, y);
    return y;
  }
  template <class AT1, class AT2, int M, class Type1>
  inline Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type1, Fixed<M>, Fixed<M>, AT1> &A, const Vector<Fixed<M>, AT2> &x) {
    Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> y(A.rows(), NONINIT);
    mult(A, x, y);
    return y;
  }
  template <class AT1, class AT2, int M, int N, class Type, class Row2>
  inline Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type, Fixed<M>, Fixed<N>, AT1> &A, const Vector<Row2, AT2> &x) {
    Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> y(A.rows(), NONINIT);
    mult(A, x, y);
    return y;
  }
  template <class AT1, class AT2, int M, int N, class Type>
  inline Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type, Fixed<M>, Fixed<N>, AT1> &A, const Vector<Fixed<N>, AT2> &x) {
    Vector<Fixed<M>, typename OperatorResult<AT1, AT2>::Type> y(A.rows(), NONINIT);
    mult(A, x, y);
    return y;
  }

  // RowVector-Matrix
  template <class Col1, class Type2, class Row2, class Col2, class Type3, class Col3, class AT1, class AT2>
  inline void mult(const RowVector<Col1, AT1> &x, const Matrix<Type2, Row2, Col2, AT2> &A, RowVector<Col3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &y) {
    assert(x.size() == A.rows());
    for (int i = 0; i < y.size(); i++) {
      y.e(i) = 0;
      for (int j = 0; j < x.size(); j++)
        y.e(i) += x.e(j) * A.e(j, i);
    }
  }

  template <class Col1, class Row2, class Col3, class AT1, class AT2>
  inline void mult(const RowVector<Col1, AT1> &x, const Matrix<Symmetric, Row2, Row2, AT2> &A, RowVector<Col3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &y) {
    assert(x.size() == A.rows());
    for (int i = 0; i < y.size(); i++) {
      y.e(i) = 0;
      for (int j = 0; j < i; j++)
        y.e(i) += x.e(j) * A.ej(j, i);
      for (int j = i; j < x.size(); j++)
        y.e(i) += x.e(j) * A.ei(j, i);
    }
  }

  template <class AT1, class AT2, class Type2, class Col2, class Row1, class Col1>
  inline RowVector<Var, typename OperatorResult<AT1, AT2>::Type> operator*(const RowVector<Col2, AT1> &x, const Matrix<Type2, Row1, Col1, AT2> &A) {
    RowVector<Var, typename OperatorResult<AT1, AT2>::Type> y(A.cols(), NONINIT);
    mult(x, A, y);
    return y;
  }
  template <class AT1, class AT2, class Type2, class Row, class Col>
  inline RowVector<Col, typename OperatorResult<AT1, AT2>::Type> operator*(const RowVector<Col, AT1> &x, const Matrix<Type2, Row, Col, AT2> &A) {
    RowVector<Col, typename OperatorResult<AT1, AT2>::Type> y(A.cols(), NONINIT);
    mult(x, A, y);
    return y;
  }
  template <class AT1, class AT2, int N, class Type2, class Col1, class Row2>
  inline RowVector<Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator*(const RowVector<Col1, AT1> &x, const Matrix<Type2, Row2, Fixed<N>, AT2> &A) {
    RowVector<Fixed<N>, typename OperatorResult<AT1, AT2>::Type> y(A.cols(), NONINIT);
    mult(x, A, y);
    return y;
  }
  template <class AT1, class AT2, int N, class Type2>
  inline RowVector<Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator*(const RowVector<Fixed<N>, AT1> &x, const Matrix<Type2, Fixed<N>, Fixed<N>, AT2> &A) {
    RowVector<Fixed<N>, typename OperatorResult<AT1, AT2>::Type> y(A.cols(), NONINIT);
    mult(x, A, y);
    return y;
  }

  // Matrix-Matrix
  template <class Type1, class Row1, class Col1, class Type2, class Row2, class Col2, class Type3, class Row3, class Col3, class AT1, class AT2>
  inline void mult(const Matrix<Type1, Row1, Col1, AT1> &A1, const Matrix<Type2, Row2, Col2, AT2> &A2, Matrix<Type3, Row3, Col3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &A3) {
    assert(A1.cols() == A2.rows());
    for (int i = 0; i < A3.rows(); i++) {
      for (int k = 0; k < A3.cols(); k++) {
        A3.e(i, k) = 0;
        for (int j = 0; j < A1.cols(); j++)
          A3.e(i, k) += A1.e(i, j) * A2.e(j, k);
      }
    }
  }
  template <class Row1, class Type2, class Row2, class Col2, class Type3, class Row3, class Col3, class AT1, class AT2>
  inline void mult(const Matrix<Symmetric, Row1, Row1, AT1> &A1, const Matrix<Type2, Row2, Col2, AT2> &A2, Matrix<Type3, Row3, Col3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &A3) {
    assert(A1.size() == A2.rows());
    for (int i = 0; i < A3.rows(); i++) {
      for (int k = 0; k < A3.cols(); k++) {
        A3.e(i, k) = 0;
        for (int j = 0; j < i; j++)
          A3.e(i, k) += A1.ei(i, j) * A2.e(j, k);
        for (int j = i; j < A1.cols(); j++)
          A3.e(i, k) += A1.ej(i, j) * A2.e(j, k);
      }
    }
  }
  template <class Row1, class Row2, class Row3, class AT1, class AT2>
  inline void mult(const Matrix<Diagonal, Row1, Row1, AT1> &A1, const Matrix<Diagonal, Row2, Row2, AT2> &A2, Matrix<Diagonal, Row3, Row3, typename fmatvec::OperatorResult<AT1, AT2>::Type> &A3) {
    assert(A1.size() == A2.size());
    for (int i = 0; i < A3.size(); i++)
      A3.e(i) = A1.e(i) * A2.e(i);
  }

  template <class AT1, class AT2, class Type1, class Type2, class Row1, class Col1, class Row2, class Col2>
  inline Matrix<General, Var, Var, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type1, Row1, Col1, AT1> &A1, const Matrix<Type2, Row2, Col2, AT2> &A2) {
    Matrix<General, Var, Var, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, class Type1, class Type2, class Row, class Col>
  inline Matrix<General, Row, Col, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type1, Row, Col, AT1> &A1, const Matrix<Type2, Row, Col, AT2> &A2) {
    Matrix<General, Row, Col, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, int M, class Type1, class Type2, class Col1, class Row2, class Col2>
  inline Matrix<General, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type1, Fixed<M>, Col1, AT1> &A1, const Matrix<Type2, Row2, Col2, AT2> &A2) {
    Matrix<General, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, int N, class Type1, class Row1, class Col1, class Type2, class Row2>
  inline Matrix<General, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type1, Row1, Col1, AT1> &A1, const Matrix<Type2, Row2, Fixed<N>, AT2> &A2) {
    Matrix<General, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, int M, int N, class Type1, class Type2, class Col1, class Row2>
  inline Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type1, Fixed<M>, Col1, AT1> &A1, const Matrix<Type2, Row2, Fixed<N>, AT2> &A2) {
    Matrix<General, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, int M, class Type1, class Type2>
  inline Matrix<General, Fixed<M>, Fixed<M>, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type1, Fixed<M>, Fixed<M>, AT1> &A1, const Matrix<Type2, Fixed<M>, Fixed<M>, AT2> &A2) {
    Matrix<General, Fixed<M>, Fixed<M>, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, class Type, class Row1, class Col1, class Row2, class Col2>
  inline Matrix<Type, Var, Var, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type, Row1, Col1, AT1> &A1, const Matrix<Type, Row2, Col2, AT2> &A2) {
    Matrix<Type, Var, Var, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, class Type, class Row, class Col>
  inline Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type, Row, Col, AT1> &A1, const Matrix<Type, Row, Col, AT2> &A2) {
    Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, int M, class Type, class Col1, class Row2, class Col2>
  inline Matrix<Type, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type, Fixed<M>, Col1, AT1> &A1, const Matrix<Type, Row2, Col2, AT2> &A2) {
    Matrix<Type, Fixed<M>, Var, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, int N, class Type, class Row1, class Col1, class Row2>
  inline Matrix<Type, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type, Row1, Col1, AT1> &A1, const Matrix<Type, Row2, Fixed<N>, AT2> &A2) {
    Matrix<Type, Var, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, int M, int N, class Type, class Col1, class Row2>
  inline Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type, Fixed<M>, Col1, AT1> &A1, const Matrix<Type, Row2, Fixed<N>, AT2> &A2) {
    Matrix<Type, Fixed<M>, Fixed<N>, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, int M, class Type>
  inline Matrix<Type, Fixed<M>, Fixed<M>, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type, Fixed<M>, Fixed<M>, AT1> &A1, const Matrix<Type, Fixed<M>, Fixed<M>, AT2> &A2) {
    Matrix<Type, Fixed<M>, Fixed<M>, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  
//  SquareMatrix-SquareMatrix
  template <class AT1, class AT2, class Row1, class Row2>
  inline SquareMatrix<Var, typename OperatorResult<AT1, AT2>::Type> operator*(const SquareMatrix<Row1, AT1> &A1, const SquareMatrix<Row2, AT2> &A2) {
    SquareMatrix<Var, typename OperatorResult<AT1, AT2>::Type> A3(A1.size(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, class Row>
  inline SquareMatrix<Row, typename OperatorResult<AT1, AT2>::Type> operator*(const SquareMatrix<Row, AT1> &A1, const SquareMatrix<Row, AT2> &A2) {
    SquareMatrix<Row, typename OperatorResult<AT1, AT2>::Type> A3(A1.size(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }

  // SymmetricMatrix-SymmetricMatrix
  template <class AT1, class AT2, class Row1, class Row2>
  inline Matrix<General, Var, Var, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Symmetric, Row1, Row1, AT1> &A1, const Matrix<Symmetric, Row2, Row2, AT2> &A2) {
    Matrix<General, Var, Var, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }
  template <class AT1, class AT2, class Row>
  inline Matrix<General, Row, Row, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Symmetric, Row, Row, AT1> &A1, const Matrix<Symmetric, Row, Row, AT2> &A2) {
    Matrix<General, Row, Row, typename OperatorResult<AT1, AT2>::Type> A3(A1.rows(), A2.cols(), NONINIT);
    mult(A1, A2, A3);
    return A3;
  }

  /////////////////////////////////// end vecscalmult //////////////////////////////

  /*! \brief Vector-scalar multiplication.
   *
   * This function computes the product of a vector 
   * and a scalar.
   * \return The product.
   * */
  template <class Row, class AT1, class AT2>
  Vector<Row, typename OperatorResult<AT1, AT2>::Type> operator*(const Vector<Row, AT1> &x, const AT2& alpha) {

    Vector<Row, typename OperatorResult<AT1, AT2>::Type> y(x.size(), NONINIT);

    for (int i = 0; i < x.size(); i++)
      y.e(i) = x.e(i) * alpha;

    return y;
  }

  /*! \brief Scalar-vector multiplication.
   *
   * \see operator*(const Vector<Ref,AT1>&x,const AT&).
   * */
  template <class Row, class AT1, class AT2>
  Vector<Row, typename OperatorResult<AT1, AT2>::Type> operator*(const AT1& alpha, const Vector<Row, AT2> &x) {

    Vector<Row, typename OperatorResult<AT1, AT2>::Type> y(x.size(), NONINIT);

    for (int i = 0; i < x.size(); i++)
      y.e(i) = x.e(i) * alpha;

    return y;
  }

  template <class Row, class AT1, class AT2>
  inline Vector<Row, AT1> operator*=(const Vector<Row, AT1> &x_, const AT2 &alpha) {
    auto &x = const_cast<Vector<Row, AT1> &>(x_);
    for (int i = 0; i < x.size(); i++)
      x.e(i) *= alpha;
    return x;
  }

  template <class Row, class AT1, class AT2>
  inline Vector<Row, typename OperatorResult<AT1, AT2>::Type> operator/(const Vector<Row, AT1> &x, const AT2 &alpha) {
    Vector<Row, typename OperatorResult<AT1, AT2>::Type> y(x.size(), NONINIT);
    for (int i = 0; i < x.size(); i++)
      y.e(i) = x.e(i) / alpha;
    return y;
  }

  template <class Row, class AT1, class AT2>
  inline Vector<Row, AT1> operator/=(const Vector<Row, AT1> &x_, const AT2 &alpha) {
    auto &x = const_cast<Vector<Row, AT1> &>(x_);
    for (int i = 0; i < x.size(); i++)
      x.e(i) /= alpha;
    return x;
  }
  /*! \brief Rowvector-scalar multiplication.
   *
   * This function computes the product of a rowvector 
   * and a scalar.
   * \return The product.
   * */
  template <class Col, class AT1, class AT2>
  RowVector<Col, typename OperatorResult<AT1, AT2>::Type> operator*(const RowVector<Col, AT1> &x, const AT2& alpha) {

    RowVector<Col, typename OperatorResult<AT1, AT2>::Type> y(x.size(), NONINIT);

    for (int i = 0; i < x.size(); i++)
      y.e(i) = x.e(i) * alpha;

    return y;
  }

  /*! \brief Scalar-rowvector multiplication.
   *
   * \see operator*(const RowVector<Col>, AT2>&, const AT&).
   * */
  template <class Col, class AT1, class AT2>
  RowVector<Col, typename OperatorResult<AT1, AT2>::Type> operator*(const AT1 &alpha, const RowVector<Col, AT2> &x) {

    RowVector<Col, typename OperatorResult<AT1, AT2>::Type> y(x.size(), NONINIT);

    for (int i = 0; i < x.size(); i++)
      y.e(i) = x.e(i) * alpha;

    return y;
  }

  template <class Col, class AT1, class AT2>
  inline RowVector<Col, AT1> operator*=(const RowVector<Col, AT1> &x_, const AT2 &alpha) {
    auto &x = const_cast<RowVector<Col, AT1> &>(x_);
    for (int i = 0; i < x.size(); i++)
      x.e(i) *= alpha;
    return x;
  }

  template <class Col, class AT1, class AT2>
  inline RowVector<Col, typename OperatorResult<AT1, AT2>::Type> operator/(const RowVector<Col, AT1> &x, const AT2 &a) {
    RowVector<Col, typename OperatorResult<AT1, AT2>::Type> y(x.size(), NONINIT);
    for (int i = 0; i < x.size(); i++)
      y.e(i) = x.e(i) / a;
    return y;
  }

  template <class Col, class AT1, class AT2>
  inline RowVector<Col, AT1> operator/=(const RowVector<Col, AT1> &x_, const AT2 &a) {
    auto &x = const_cast<RowVector<Col, AT1> &>(x_);
    for (int i = 0; i < x.size(); i++)
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
  template <class Col, class Row, class AT1, class AT2>
  typename OperatorResult<AT1, AT2>::Type operator*(const RowVector<Col, AT1> &x, const Vector<Row, AT2> &y) {

    assert(x.size() == y.size());

    typename OperatorResult<AT1, AT2>::Type res = 0;

    for (int i = 0; i < x.size(); i++)
      res += x.e(i) * y.e(i);

    return res;
  }

  /*! \brief Scalar product.
   *
   * This function computes the scalar product of two vectors.
   * \return The scalar product.
   * */
  template <class Row, class AT>
  AT scalarProduct(const Vector<Row, AT> &x, const Vector<Row, AT> &y) {

    assert(x.size() == y.size());

    AT res = 0;

    for (int i = 0; i < x.size(); i++)
      res += x.e(i) * y.e(i);

    return res;
  }

  /*! \brief Cross product.
   *
   * This function computes the cross product of two vectors.
   * \return The cross product.
   * */
  template <class Row1, class Row2, class AT1, class AT2>
  Vector<Fixed<3>, typename OperatorResult<AT1, AT2>::Type> crossProduct(const Vector<Row1, AT1> &x, const Vector<Row2, AT2> &y) {

    assert(x.size() == 3);
    assert(y.size() == 3);

    Vector<Fixed<3>, typename OperatorResult<AT1, AT2>::Type> z(3, NONINIT);

    z.e(0) = x.e(1) * y.e(2) - x.e(2) * y.e(1);
    z.e(1) = x.e(2) * y.e(0) - x.e(0) * y.e(2);
    z.e(2) = x.e(0) * y.e(1) - x.e(1) * y.e(0);

    return z;
  }

  /*! \brief Triple product.
   *
   * This function computes the triple product of three vectors.
   * \return The triple product.
   * */
  template <class Row, class AT>
  double tripleProduct(const Vector<Row, AT> &a, const Vector<Row, AT> &x, const Vector<Row, AT> &y) {

    assert(a.size() == 3);
    assert(x.size() == 3);
    assert(y.size() == 3);

    return a.e(0) * (x.e(1) * y.e(2) - x.e(2) * y.e(1)) + a.e(1) * (x.e(2) * y.e(0) - x.e(0) * y.e(2)) + a.e(2) * (x.e(0) * y.e(1) - x.e(1) * y.e(0));
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
  template <class Type, class Row, class Col, class AT1, class AT2>
  Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Type, Row, Col, AT1> &A, const AT2 &alpha) {

    Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> B(A.rows(), A.cols(), NONINIT);

    for (int i = 0; i < A.rows(); i++)
      for (int j = 0; j < A.cols(); j++)
        B.e(i, j) = A.e(i, j) * alpha;

    return B;
  }

  template <class Type, class Row, class Col, class AT1, class AT2>
  Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> operator*(const AT1 &alpha, const Matrix<Type, Row, Col, AT2> &A) {

    Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> B(A.rows(), A.cols(), NONINIT);

    for (int i = 0; i < A.rows(); i++)
      for (int j = 0; j < A.cols(); j++)
        B.e(i, j) = A.e(i, j) * alpha;

    return B;
  }

  template <class Row, class AT1, class AT2>
  inline Matrix<Symmetric, Row, Row, typename OperatorResult<AT1, AT2>::Type> operator*(const AT1 &alpha, const Matrix<Symmetric, Row, Row, AT2> &A) {
    Matrix<Symmetric, Row, Row, typename OperatorResult<AT1, AT2>::Type> B(A.size(), A.size(), NONINIT);
    for (int i = 0; i < A.size(); i++)
      for (int j = i; j < A.size(); j++)
        B.ej(i, j) = A.ej(i, j) * alpha;
    return B;
  }

  template <class Row, class AT1, class AT2>
  inline Matrix<Symmetric, Row, Row, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Symmetric, Row, Row, AT1> &A, const AT2 &alpha) {
    Matrix<Symmetric, Row, Row, typename OperatorResult<AT1, AT2>::Type> B(A.size(), A.size(), NONINIT);
    for (int i = 0; i < A.size(); i++)
      for (int j = i; j < A.size(); j++)
        B.ej(i, j) = A.ej(i, j) * alpha;
    return B;
  }

  template <class AT1, class AT2>
  inline Matrix<Diagonal, Ref, Ref, typename OperatorResult<AT1, AT2>::Type> operator*(const AT1 &alpha, const Matrix<Diagonal, Ref, Ref, AT2> &A) {
    Matrix<Diagonal, Ref, Ref, typename OperatorResult<AT1, AT2>::Type> B(A.rows(), A.cols(), NONINIT);
    for (int i = 0; i < A.rows(); i++)
      B.e(i) = A.e(i) * alpha;
    return B;
  }

  template <class AT1, class AT2>
  inline Matrix<Diagonal, Ref, Ref, typename OperatorResult<AT1, AT2>::Type> operator*(const Matrix<Diagonal, Ref, Ref, AT1> &A, const AT2 &alpha) {
    Matrix<Diagonal, Ref, Ref, typename OperatorResult<AT1, AT2>::Type> B(A.rows(), A.cols(), NONINIT);
    for (int i = 0; i < A.size(); i++)
      B.e(i) = A.e(i) * alpha;
    return B;
  }

  template <class Type, class Row, class Col, class AT1, class AT2>
  Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> operator/(const Matrix<Type, Row, Col, AT1> &A, const AT2 &alpha) {

    Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> B(A.rows(), A.cols(), NONINIT);

    for (int i = 0; i < A.rows(); i++)
      for (int j = 0; j < A.cols(); j++)
        B.e(i, j) = A.e(i, j) / alpha;

    return B;
  }

  template <class Row, class AT1, class AT2>
  inline Matrix<Symmetric, Row, Row, typename OperatorResult<AT1, AT2>::Type> operator/(const Matrix<Symmetric, Row, Row, AT1> &A, const AT2 &alpha) {
    Matrix<Symmetric, Row, Row, typename OperatorResult<AT1, AT2>::Type> B(A.size(), A.size(), NONINIT);
    for (int i = 0; i < A.size(); i++)
      for (int j = i; j < A.size(); j++)
        B.ej(i, j) = A.ej(i, j) / alpha;
    return B;
  }

  template <class Row, class AT1, class AT2>
  inline Matrix<Diagonal, Row, Row, typename OperatorResult<AT1, AT2>::Type> operator/(const Matrix<Diagonal, Row, Row, AT1> &A, const AT2 &alpha) {
    Matrix<Diagonal, Row, Row, typename OperatorResult<AT1, AT2>::Type> B(A.size(), A.size(), NONINIT);
    for (int i = 0; i < A.size(); i++)
      B.e(i) = A.e(i) / alpha;
    return B;
  }

  template <class Row, class AT1, class AT2>
  SquareMatrix<Row, typename OperatorResult<AT1, AT2>::Type> operator*(const SquareMatrix<Row, AT1> &A, const AT2 &alpha) {

    SquareMatrix<Row, typename OperatorResult<AT1, AT2>::Type> B(A.size(), NONINIT);

    for (int i = 0; i < A.size(); i++)
      for (int j = 0; j < A.size(); j++)
        B.e(i, j) = A.e(i, j) * alpha;

    return B;
  }

  template <class Row, class AT1, class AT2>
  SquareMatrix<Row, typename OperatorResult<AT1, AT2>::Type> operator*(const AT1 &alpha, const SquareMatrix<Row, AT2> &A) {

    SquareMatrix<Row, typename OperatorResult<AT1, AT2>::Type> B(A.size(), NONINIT);

    for (int i = 0; i < A.size(); i++)
      for (int j = 0; j < A.size(); j++)
        B.e(i, j) = A.e(i, j) * alpha;

    return B;
  }

  template <class Row, class AT1, class AT2>
  SquareMatrix<Row, typename OperatorResult<AT1, AT2>::Type> operator/(const SquareMatrix<Row, AT1> &A, const AT2 &alpha) {

    SquareMatrix<Row, typename OperatorResult<AT1, AT2>::Type> B(A.size(), NONINIT);

    for (int i = 0; i < A.size(); i++)
      for (int j = 0; j < A.size(); j++)
        B.e(i, j) = A.e(i, j) / alpha;

    return B;
  }

  template <class Type, class Row, class Col, class AT1, class AT2>
  Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> operator*=(const Matrix<Type, Row, Col, AT1> &A_, const AT2 &alpha) {
    auto &A = const_cast<Matrix<Type, Row, Col, typename OperatorResult<AT1, AT2>::Type> &>(A_);
    for (int i = 0; i < A.rows(); i++)
      for (int j = 0; j < A.cols(); j++)
        A.e(i, j) *= alpha;

    return A;
  }

  template <class Row, class AT1, class AT2>
  inline Matrix<Symmetric, Row, Row, AT1> operator*=(const Matrix<Symmetric, Row, Row, AT1> &A_, const AT2 &alpha) {
    auto &A = const_cast<Matrix<Symmetric, Row, Row, AT1> &>(A_);
    for (int i = 0; i < A.size(); i++)
      for (int j = i; j < A.size(); j++)
        A.ej(i, j) *= alpha;
    return A;
  }

  template <class Row, class AT1, class AT2>
  inline Matrix<Diagonal, Row, Row, AT1> operator*=(const Matrix<Diagonal, Row, Row, AT1> &A_, const AT2 &alpha) {
    auto &A = const_cast<Matrix<Diagonal, Row, Row, AT1> &>(A_);
    for (int i = 0; i < A.size(); i++)
      A.e(i) *= alpha;
    return A;
  }

  template <class Type, class Row, class Col, class AT1, class AT2>
  Matrix<Type, Row, Col, AT1> operator/=(const Matrix<Type, Row, Col, AT1> &A_, const AT2 &alpha) {
    Matrix<Type, Row, Col, AT1> A = const_cast<Matrix<Type, Row, Col, AT1> &>(A_);
    for (int i = 0; i < A.rows(); i++)
      for (int j = 0; j < A.cols(); j++)
        A.e(i, j) /= alpha;

    return A;
  }

  template <class Row, class AT1, class AT2>
  inline Matrix<Symmetric, Row, Row, AT1> operator/=(const Matrix<Symmetric, Row, Row, AT1> &A_, const AT2 &alpha) {
    Matrix<Symmetric, Row, Row, AT1> A = const_cast<Matrix<Symmetric, Row, Row, AT1> &>(A_);
    for (int i = 0; i < A.size(); i++)
      for (int j = i; j < A.size(); j++)
        A.ej(i, j) /= alpha;
    return A;
  }

  template <class Row, class AT1, class AT2>
  inline Matrix<Diagonal, Row, Row, AT1> operator/=(const Matrix<Diagonal, Row, Row, AT1> &A_, const AT2 &alpha) {
    Matrix<Diagonal, Row, Row, AT1> A = const_cast<Matrix<Diagonal, Row, Row, AT1> &>(A_);
    for (int i = 0; i < A.size(); i++)
      A.e(i) /= alpha;
    return A;
  }

  /////////////////////////////////// end matscalmult //////////////////////////////
  
  /////////////////////////////////// negation //////////////////////////////

  /*! \brief Negation.
   *
   * This function computes the negation of a vector.
   * \return The negation.
   * */
  template <class Row, class AT>
  Vector<Row, typename OperatorResult<AT, AT>::Type> operator-(const Vector<Row, AT> &x) {

    Vector<Row, typename OperatorResult<AT, AT>::Type> y(x.size(), NONINIT);

    for (int i = 0; i < x.size(); i++)
      y.e(i) = -x.e(i);

    return y;
  }

  /*! \brief Negation.
   *
   * This function computes the negation of a rowvector.
   * \return The negation.
   * */
  template <class Col, class AT>
  RowVector<Col, typename OperatorResult<AT, AT>::Type> operator-(const RowVector<Col, AT> &a) {

    RowVector<Col, typename OperatorResult<AT, AT>::Type> c(a.size(), NONINIT);

    for (int i = 0; i < a.size(); i++)
      c.e(i) = -a.e(i);

    return c;
  }

  /*! \brief Negation.
   *
   * This function computes the negation of a vector.
   * \return The negation.
   * */
  template <class Row, class AT>
  SquareMatrix<Row, typename OperatorResult<AT, AT>::Type> operator-(const SquareMatrix<Row, AT> &A) {

    SquareMatrix<Row, typename OperatorResult<AT, AT>::Type> B(A.size(), NONINIT);

    for (int i = 0; i < A.size(); i++)
      for (int j = 0; j < A.size(); j++)
        B.e(i, j) = -A.e(i, j);

    return B;
  }

  /*! \brief Negation.
   *
   * This function computes the negation of a vector.
   * \return The negation.
   * */
  template <class Type, class Row, class Col, class AT>
  Matrix<Type, Row, Col, typename OperatorResult<AT, AT>::Type> operator-(const Matrix<Type, Row, Col, AT> &A) {

    Matrix<Type, Row, Col, typename OperatorResult<AT, AT>::Type> B(A.rows(), A.cols(), NONINIT);

    for (int i = 0; i < A.rows(); i++)
      for (int j = 0; j < A.cols(); j++)
        B.e(i, j) = -A.e(i, j);

    return B;
  }

  /////////////////////////////////// end negation //////////////////////////////
  /////////////////////////////////// transpose //////////////////////////////

  /*! \brief Transpose of a vector.
   *
   * This function computes the transpose of a vector. The result is a rowvector.
   * \return The transpose.
   * */
  template <class Row, class AT>
  RowVector<Row,AT> trans(const Vector<Row,AT> &x) { return x.T(); }

  /*! \brief Transpose of a rowvector.
   *
   * This function computes the transpose of a rowvector.
   * The result is a vector.
   * \param x A vector.
   * \return The transpose.
   * */
  template <class Col, class AT>
  Vector<Col,AT> trans(const RowVector<Col,AT> &x) { return x.T(); }

  /*! \brief Transpose of a matrix.
   *
   * This function computes the transpose of a matrix.
   * \f[ \boldsymbol{A} \rightarrow \boldsymbol{A}^T  \f]
   * The result is a matrix.
   * \param A A general matrix.
   * \return The transpose.
   * */
  template <class Type, class Row, class Col, class AT>
  Matrix<Type,Col,Row,AT> trans(const Matrix<Type,Row,Col,AT> &A) { return A.T(); }

  /*! \brief Transpose of a matrix.
   *
   * This function computes the transpose of a square matrix.
   * \f[ \boldsymbol{A} \rightarrow \boldsymbol{A}^T  \f]
   * The result is a square matrix.
   * \param A A square matrix.
   * \return The transpose.
   * */
  template <class Row, class AT>
  SquareMatrix<Row,AT> trans(const SquareMatrix<Row,AT> &A) { return A.T(); }

//  template <class Type, class Row, class Col, class AT>
//  inline Matrix<Type, Col, Row, AT> trans(const Matrix<Type, Row, Col, AT> &A) {
//    Matrix<Type, Col, Row, AT> B(A.cols(), A.rows(), NONINIT);
//    for (int i = 0; i < B.rows(); i++)
//      for (int j = 0; j < B.cols(); j++)
//        B.e(i, j) = A.e(j, i);
//    return B;
//  }
//
//  template <class Row, class AT>
//  inline SquareMatrix<Row, AT> trans(const SquareMatrix<Row, AT> &A) {
//    SquareMatrix<Row, AT> B(A.size(), NONINIT);
//    for (int i = 0; i < B.size(); i++)
//      for (int j = 0; j < B.size(); j++)
//        B.e(i, j) = A.e(j, i);
//    return B;
//  }
//
//  template <class Row, class AT>
//  inline RowVector<Row, AT> trans(const Vector<Row, AT> &x) {
//    RowVector<Row, AT> y(x.size(), NONINIT);
//    for (int i = 0; i < y.size(); i++)
//      y.e(i) = x.e(i);
//    return y;
//  }
//
//  template <class Row, class AT>
//  inline Vector<Row, AT> trans(const RowVector<Row, AT> &x) {
//    Vector<Row, AT> y(x.size(), NONINIT);
//    for (int i = 0; i < y.size(); i++)
//      y.e(i) = x.e(i);
//    return y;
//  }

  /////////////////////////////////// end transpose //////////////////////////////

  /////////////////////////////////// norms //////////////////////////////////////
  template <class Row, class AT>
  inline AT nrmInf(const Vector<Row, AT> &x) {
    AT c = 0;
    for (int i = 0; i < x.size(); i++)
      if (c < fabs(x.e(i)))
        c = fabs(x.e(i));
    return c;
  }

  template <class Row, class AT>
  inline AT nrmInf(const RowVector<Row, AT> &x) {
    AT c = 0;
    for (int i = 0; i < x.size(); i++)
      if (c < fabs(x.e(i)))
        c = fabs(x.e(i));
    return c;
  }

  template <class Row, class AT>
  inline typename OperatorResult<AT, AT>::Type nrm2(const Vector<Row, AT> &x) {
    typename OperatorResult<AT, AT>::Type c=typename OperatorResult<AT, AT>::Type();
    for (int i = 0; i < x.size(); i++)
      c += pow(x.e(i), 2);
    return sqrt(c);
  }

  template <class Col, class AT>
  inline typename OperatorResult<AT, AT>::Type nrm2(const RowVector<Col, AT> &x) {
    typename OperatorResult<AT, AT>::Type c=typename OperatorResult<AT, AT>::Type();
    for (int i = 0; i < x.size(); i++)
      c += pow(x.e(i), 2);
    return sqrt(c);
  }

  /////////////////////////////////// end norms //////////////////////////////////

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
  template <class Row, class AT>
  SquareMatrix<Row, AT> tilde(const Vector<Row, AT> &x) {

    SquareMatrix<Row, AT> B(x.size(), NONINIT);

    B.e(0, 0) = 0;
    B.e(1, 1) = 0;
    B.e(2, 2) = 0;
    B.e(0, 1) = -x.e(2);
    B.e(0, 2) = x.e(1);
    B.e(1, 0) = x.e(2);
    B.e(1, 2) = -x.e(0);
    B.e(2, 0) = -x.e(1);
    B.e(2, 1) = x.e(0);

    return B;
  }

  /*! \brief Matrix-matrix comparison.
   *
   * This function compares two matrices
   * \return True, if the matrices are identical, false otherwise.
   * */
  template <class AT, class Type1, class Row1, class Col1, class Type2, class Row2, class Col2>
  bool operator==(const Matrix<Type1, Row1, Col1, AT> &A, const Matrix<Type2, Row2, Col2, AT> &B) {

    if (&A == &B)
      return true;

    if (A.rows() != B.rows() || A.cols() != B.cols())
      return false;

    for (int i = 0; i < A.rows(); i++)
      for (int j = 0; j < A.cols(); j++)
        if (A(i, j) != B(i, j))
          return false;

    return true;
  }

  /*! \brief Matrix-matrix comparison.
   *
   * This function compares two matrices.
   * \return True, if the matrices are different, false otherwise.
   * */
  template <class AT, class Type1, class Row1, class Col1, class Type2, class Row2, class Col2>
  bool operator!=(const Matrix<Type1, Row1, Col1, AT> &A, const Matrix<Type2, Row2, Col2, AT> &B) {

    if (A.rows() != B.rows() || A.cols() != B.cols())
      return true;

    for (int i = 0; i < A.rows(); i++)
      for (int j = 0; j < A.cols(); j++)
        if (A(i, j) != B(i, j))
          return true;

    return false;
  }

  /*! \brief Maximum value.
   *
   * This function computes the maximum value of a vector.
   * \return The maximum value.
   * */
  template <class Row, class AT>
  AT max(const Vector<Row, AT> &x) {

    assert(x.size() > 0);
    AT maximum = x.e(0);
    for (int i = 1; i < x.size(); i++) {
      if (x.e(i) > maximum)
        maximum = x.e(i);
    }
    return maximum;
  }

  /*! \brief Index of maximum value.
   *
   * This function computes the index of the maximum value of a vector
   * \return The index of the maximum value.
   * */
  template <class Row, class AT>
  int maxIndex(const Vector<Row, AT> &x) {

    assert(x.size() > 0);
    AT maximum = x.e(0);
    int index = 0;
    for (int i = 1; i < x.size(); i++) {
      if (x.e(i) > maximum) {
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
  template <class Row, class AT>
  AT min(const Vector<Row, AT> &x) {

    assert(x.size() > 0);
    AT minimum = x.e(0);
    for (int i = 1; i < x.size(); i++) {
      if (x.e(i) < minimum)
        minimum = x.e(i);
    }
    return minimum;
  }

  /*! \brief Index of minimum value.
   *
   * This function computes the index of the minimum value of a vector
   * \return The index of the minimum value.
   * */
  template <class Row, class AT>
  int minIndex(const Vector<Row, AT> &x) {

    assert(x.size() > 0);
    AT minimum = x.e(0);
    int index = 0;
    for (int i = 1; i < x.size(); i++) {
      if (x.e(i) < minimum) {
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
  Matrix<General, Ref, Ref, AT> bubbleSort(const Matrix<General, Ref, Ref, AT> &A_, int PivotCol) {

    assert(A_.rows() > 0);
    assert(A_.cols() > 0);
    assert(A_.cols() > PivotCol);
    Matrix<General, Ref, Ref, AT> A = A_;
    int i, j, N;
    RowVector<Ref, AT> tmp(A.cols(), NONINIT);
    N = A.rows();
    for (i = 1; i <= N - 1; i++) {
      for (j = 0; j < N - 1; j++) {
        if (A(j, PivotCol) > A(j + 1, PivotCol)) {
          tmp = A.row(j);
          A.row(j) = A.row(j + 1);
          A.row(j + 1) = tmp;
        }
      }
    }
    return A;
  }

  // HR 28.09.2006
  /*! internal function of QuickSortMedian */
  template <class AT>
  void quicksortmedian_intern(Matrix<General, Ref, Ref, AT> &A, int PivotCol, RowVector<Ref, AT> &tmp, int l, int r) {

    if (r > l) {
      int i = l - 1, j = r;
      //RowVec tmp;
      if (r - l > 3) { //Median of three
        int m = l + (r - l) / 2;
        if (A(l, PivotCol) > A(m, PivotCol)) {
          tmp = A.row(l);
          A.row(l) = A.row(m);
          A.row(m) = tmp;
        }
        if (A(l, PivotCol) > A(r, PivotCol)) {
          tmp = A.row(l);
          A.row(l) = A.row(r);
          A.row(r) = tmp;
        }
        else if (A(r, PivotCol) > A(m, PivotCol)) {
          tmp = A.row(r);
          A.row(r) = A.row(m);
          A.row(m) = tmp;
        }
      }

      for (;;) {
          while(A(++i,PivotCol) < A(r,PivotCol));
          while(A(--j,PivotCol) > A(r,PivotCol) && j>i);
          if(i>=j) break;
        tmp = A.row(i);
        A.row(i) = A.row(j);
        A.row(j) = tmp;
      }
      tmp = A.row(i);
      A.row(i) = A.row(r);
      A.row(r) = tmp;
      quicksortmedian_intern(A, PivotCol, tmp, l, i - 1);
      quicksortmedian_intern(A, PivotCol, tmp, i + 1, r);
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
  Matrix<General, Ref, Ref, AT> quickSortMedian(const Matrix<General, Ref, Ref, AT> &A_, int PivotCol) {
    Matrix<General, Ref, Ref, AT> A = A_;
    int N = A.rows();
    RowVector<Ref, AT> tmp(A.cols(), NONINIT);
    quicksortmedian_intern(A, PivotCol, tmp, 0, N - 1);
    return A;
  }

  template <class Row, class Col, class AT>
  inline Matrix<Symmetric, Col, Col, AT> JTJ(const Matrix<General, Row, Col, AT> &A) {
    Matrix<Symmetric, Col, Col, AT> S(A.cols(), A.cols(), NONINIT);
    for (int i = 0; i < A.cols(); i++) {
      for (int k = i; k < A.cols(); k++) {
        S.ej(i, k) = 0;
        for (int j = 0; j < A.rows(); j++)
          S.ej(i, k) += A.e(j, i) * A.e(j, k);
      }
    }
    return S;
  }

  template <class Row, class Col, class AT>
  inline Matrix<Symmetric, Col, Col, AT> JTMJ(const Matrix<Symmetric, Row, Row, AT> &B, const Matrix<General, Row, Col, AT> &A) {

    Matrix<Symmetric, Col, Col, AT> S(A.cols(), A.cols(), NONINIT);
    Matrix<General, Row, Col, AT> C = B * A;

    for (int i = 0; i < A.cols(); i++) {
      for (int k = i; k < A.cols(); k++) {
        S.ej(i, k) = 0;
        for (int j = 0; j < A.rows(); j++)
          S.ej(i, k) += A.e(j, i) * C.e(j, k);
      }
    }
    return S;
  }

  template <class Row, class Col, class AT>
  inline Matrix<Symmetric, Row, Row, AT> JMJT(const Matrix<General, Row, Col, AT> &A, const Matrix<Symmetric, Col, Col, AT> &B) {

    Matrix<Symmetric, Row, Row, AT> S(A.rows(), A.rows(), NONINIT);
    Matrix<General, Row, Col, AT> C = A * B;

    for (int i = 0; i < S.size(); i++) {
      for (int k = i; k < S.size(); k++) {
        S.ej(i, k) = 0;
        for (int j = 0; j < A.cols(); j++)
          S.ej(i, k) += C.e(k, j) * A.e(i, j);
      }
    }
    return S;
  }

}

#endif

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

#ifndef fmatrixvector_h
#define fmatrixvector_h

#include "linear_algebra.h"
#include "linear_algebra_double.h"

namespace fmatvec {

  class General;
  class Ref;
  class Diagonal;
  class Symmetric;
  class GeneralBand;
  class Sparse;
  class Rotation;
  class Var;
  template<int N> class Fixed;

  template<class Type, class Row, class Col, class AT>
  class Matrix;

  template<class Type, class AT>
  class SquareMatrix;

  template<class Row, class AT>
  class Vector;

  template<class Col, class AT>
  class RowVector;

  /*! 
   *  \brief General matrix
   *  
   *  Mat is an abbreviation for general matrices of type double.
   * */
  typedef Matrix<General, Ref, Ref, double> Mat;

  /*! 
   *  \brief Vector
   *  
   *  Vec is an abbreviation for vectors of type double.
   * */
  typedef Vector<Ref, double> Vec;

  /*! 
   *  \brief Row vector
   *  
   *  RowVec is an abbreviation for row vectors of type double.
   * */
  typedef RowVector<Ref, double> RowVec;

  /*! 
   *  \brief Square matrix
   *  
   *  SqrMat is an abbreviation for square matrices of type double.
   * */
  typedef SquareMatrix<Ref, double> SqrMat;

  /*! 
   *  \brief Diagonal matrix
   *  
   *  DiagMat is an abbreviation for diagonal matrices of type double.
   * */
  typedef Matrix<Diagonal, Ref, Ref, double> DiagMat;

  /*! 
   *  \brief Symmetric matrix
   *  
   *  SymMat is an abbreviation for symmetric matrices of type double.
   * */
  typedef Matrix<Symmetric, Ref, Ref, double> SymMat;

  /*! 
   *  \brief Band matrix
   *  
   *  BandMat is an abbreviation for band matrices of type double.
   * */
  typedef Matrix<GeneralBand, Ref, Ref, double> BandMat;

  typedef Matrix<Sparse, Ref, Ref, double> SparseMat;

  typedef Vector<Ref, int> VecInt;
  typedef Vector<Var, int> VecVI;

  /* Typedefs for Fixed types */

  //Column Vectors
  typedef Vector<Fixed<1>, double> Vec1;
  typedef Vector<Fixed<2>, double> Vec2;
  typedef Vector<Fixed<2>, int> Vec2I;
  typedef Vector<Fixed<3>, double> Vec3;
  typedef Vector<Fixed<4>, double> Vec4;

  typedef Vector<Fixed<11>, double> Vec11;
  typedef Vector<Fixed<16>, double> Vec16;

  //Row Vectors
  typedef RowVector<Fixed<1>, double> RowVec1;
  typedef RowVector<Fixed<2>, double> RowVec2;
  typedef RowVector<Fixed<3>, double> RowVec3;
  typedef RowVector<Fixed<4>, double> RowVec4;

  typedef RowVector<Fixed<11>, double> RowVec11;
  typedef RowVector<Fixed<16>, double> RowVec16;

  //General Matrices
  typedef Matrix<General, Fixed<3>, Fixed<3>, double> Mat3x3;
  typedef Matrix<General, Fixed<3>, Fixed<2>, double> Mat3x2;
  typedef Matrix<General, Fixed<3>, Fixed<11>, double> Mat3x11;
  typedef Matrix<General, Fixed<3>, Fixed<16>, double> Mat3x16;

  typedef Matrix<General, Fixed<4>, Fixed<16>, double> Mat4x16;

  typedef Matrix<General, Fixed<16>, Fixed<3>, double> Mat16x3;
  typedef Matrix<General, Fixed<16>, Fixed<4>, double> Mat16x4;

  //Square Matrices
  typedef SquareMatrix<Fixed<1>, double> SqrMat1;
  typedef SquareMatrix<Fixed<2>, double> SqrMat2;
  typedef SquareMatrix<Fixed<3>, double> SqrMat3;
  typedef SquareMatrix<Fixed<11>, double> SqrMat11;
  typedef SquareMatrix<Fixed<16>, double> SqrMat16;

  //Symmetric Matrices
  typedef Matrix<Symmetric, Fixed<3>, Fixed<3>, double> SymMat3;
  typedef Matrix<Symmetric, Fixed<4>, Fixed<4>, double> SymMat4;

  typedef Matrix<Symmetric, Var, Var, double> SymMatV;

  typedef Matrix<General, Var, Var, double> MatV;

  typedef Matrix<General, Var, Var, int> MatVI;

  typedef SquareMatrix<Var, double> SqrMatV;

  typedef Vector<Var, double> VecV;

  typedef RowVector<Var, double> RowVecV;
  typedef RowVector<Var, int> RowVecVI;

  //Rotation Matrices
  typedef Matrix<Rotation, Fixed<3>, Fixed<3>, double> RotMat3;

  /*Typedefs for mixed size */

  typedef Matrix<General, Fixed<2>, Var, double> Mat2xV;
  typedef Matrix<General, Var, Fixed<2>, int> MatVx2I;
  typedef Matrix<General, Fixed<3>, Var, double> Mat3xV;
  typedef Matrix<General, Fixed<3>, Var, int> Mat3xVI;

  typedef Matrix<General, Var, Fixed<2>, double> MatVx2;
  typedef Matrix<General, Var, Fixed<3>, double> MatVx3;
  typedef Matrix<General, Var, Fixed<3>, int> MatVx3I;



}

#endif


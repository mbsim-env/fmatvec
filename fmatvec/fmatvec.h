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
  class SymmetricSparse;
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
  using Mat = Matrix<General, Ref, Ref, double>;

  /*! 
   *  \brief Vector
   *  
   *  Vec is an abbreviation for vectors of type double.
   * */
  using Vec = Vector<Ref, double>;

  /*! 
   *  \brief Row vector
   *  
   *  RowVec is an abbreviation for row vectors of type double.
   * */
  using RowVec = RowVector<Ref, double>;

  /*! 
   *  \brief Square matrix
   *  
   *  SqrMat is an abbreviation for square matrices of type double.
   * */
  using SqrMat = SquareMatrix<Ref, double>;

  /*! 
   *  \brief Diagonal matrix
   *  
   *  DiagMat is an abbreviation for diagonal matrices of type double.
   * */
  using DiagMat = Matrix<Diagonal, Ref, Ref, double>;

  /*! 
   *  \brief Symmetric matrix
   *  
   *  SymMat is an abbreviation for symmetric matrices of type double.
   * */
  using SymMat = Matrix<Symmetric, Ref, Ref, double>;

  /*! 
   *  \brief Band matrix
   *  
   *  BandMat is an abbreviation for band matrices of type double.
   * */
  using BandMat = Matrix<GeneralBand, Ref, Ref, double>;

  using SparseMat = Matrix<Sparse, Ref, Ref, double>;

  using SymSparseMat = Matrix<SymmetricSparse, Ref, Ref, double>;

  using VecInt = Vector<Ref, int>;
  using VecVI = Vector<Var, int>;

  /* Typedefs for Fixed types */

  //Column Vectors
  using Vec1 = Vector<Fixed<1>, double>;
  using Vec2 = Vector<Fixed<2>, double>;
  using Vec2I = Vector<Fixed<2>, int>;
  using Vec3 = Vector<Fixed<3>, double>;
  using Vec4 = Vector<Fixed<4>, double>;

  using Vec11 = Vector<Fixed<11>, double>;
  using Vec16 = Vector<Fixed<16>, double>;

  //Row Vectors
  using RowVec1 = RowVector<Fixed<1>, double>;
  using RowVec2 = RowVector<Fixed<2>, double>;
  using RowVec3 = RowVector<Fixed<3>, double>;
  using RowVec4 = RowVector<Fixed<4>, double>;

  using RowVec11 = RowVector<Fixed<11>, double>;
  using RowVec16 = RowVector<Fixed<16>, double>;

  //General Matrices
  using Mat3x3 = Matrix<General, Fixed<3>, Fixed<3>, double>;
  using Mat3x2 = Matrix<General, Fixed<3>, Fixed<2>, double>;
  using Mat3x11 = Matrix<General, Fixed<3>, Fixed<11>, double>;
  using Mat3x16 = Matrix<General, Fixed<3>, Fixed<16>, double>;

  using Mat4x16 = Matrix<General, Fixed<4>, Fixed<16>, double>;

  using Mat16x3 = Matrix<General, Fixed<16>, Fixed<3>, double>;
  using Mat16x4 = Matrix<General, Fixed<16>, Fixed<4>, double>;

  //Square Matrices
  using SqrMat1 = SquareMatrix<Fixed<1>, double>;
  using SqrMat2 = SquareMatrix<Fixed<2>, double>;
  using SqrMat3 = SquareMatrix<Fixed<3>, double>;
  using SqrMat11 = SquareMatrix<Fixed<11>, double>;
  using SqrMat16 = SquareMatrix<Fixed<16>, double>;

  //Symmetric Matrices
  using SymMat3 = Matrix<Symmetric, Fixed<3>, Fixed<3>, double>;
  using SymMat4 = Matrix<Symmetric, Fixed<4>, Fixed<4>, double>;

  using SymMatV = Matrix<Symmetric, Var, Var, double>;

  using MatV = Matrix<General, Var, Var, double>;

  using MatVI = Matrix<General, Var, Var, int>;

  using SqrMatV = SquareMatrix<Var, double>;

  using VecV = Vector<Var, double>;

  using RowVecV = RowVector<Var, double>;
  using RowVecVI = RowVector<Var, int>;

  //Rotation Matrices
  using RotMat3 = Matrix<Rotation, Fixed<3>, Fixed<3>, double>;

  /*Typedefs for mixed size */

  using Mat2xV = Matrix<General, Fixed<2>, Var, double>;
  using MatVx2I = Matrix<General, Var, Fixed<2>, int>;
  using Mat3xV = Matrix<General, Fixed<3>, Var, double>;
  using Mat3xVI = Matrix<General, Fixed<3>, Var, int>;

  using MatVx2 = Matrix<General, Var, Fixed<2>, double>;
  using MatVx3 = Matrix<General, Var, Fixed<3>, double>;
  using MatVx3I = Matrix<General, Var, Fixed<3>, int>;

  using RangeV = Range<Var, Var>;

}

#endif

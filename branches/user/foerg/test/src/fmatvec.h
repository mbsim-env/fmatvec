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
#include "linear_algebra_fixed.h"
#include "linear_algebra_var.h"
#include "linear_algebra_fixed_var.h"
#include "linear_algebra_var_fixed.h"
#include "linear_algebra_double.h"

namespace fmatvec {

  /*! 
   *  \brief General matrix
   *  
   *  Mat is an abbreviation for general matrices of type double.
   * */
  typedef  Matrix<General<Ref,Ref>, double> Mat;

  /*! 
   *  \brief Vector
   *  
   *  Vec is an abbreviation for vectors of type double.
   * */
  typedef  Vector<General<Ref,Fixed<1> >, double> Vec;

  /*! 
   *  \brief Row vector
   *  
   *  RowVec is an abbreviation for row vectors of type double.
   * */
  typedef  RowVector<General<Fixed<1>,Ref>, double> RowVec;

  /*! 
   *  \brief Square matrix
   *  
   *  SqrMat is an abbreviation for square matrices of type double.
   * */
  typedef  SquareMatrix<General<Ref,Ref>, double> SqrMat;

  /*! 
   *  \brief Diagonal matrix
   *  
   *  DiagMat is an abbreviation for diagonal matrices of type double.
   * */
  typedef  Matrix<Diagonal<Ref,Ref>, double> DiagMat;

  /*! 
   *  \brief Symmetric matrix
   *  
   *  SymMat is an abbreviation for symmetric matrices of type double.
   * */
  typedef  Matrix<Symmetric<Ref,Ref>, double> SymMat;

  /*! 
   *  \brief Band matrix
   *  
   *  BandMat is an abbreviation for band matrices of type double.
   * */
  typedef  Matrix<GeneralBand<Ref,Ref>, double> BandMat;

  typedef  Vector<General<Ref,Fixed<1> >, int> VecInt;

  typedef  Matrix<General<Fixed<3>,Fixed<3> >,double> Mat33;

  typedef  Matrix<General<Fixed<3>,Fixed<2> >,double> Mat32;

  typedef  SquareMatrix<General<Fixed<3>,Fixed<3> >,double> SqrMat3;

  typedef  Vector<General<Fixed<3>,Fixed<1> >,double> Vec3;

  typedef  Vector<General<Fixed<2>,Fixed<1> >,double> Vec2;

  typedef  RowVector<General<Fixed<1>,Fixed<3> >,double> RowVec3;

  typedef  Matrix<Symmetric<Fixed<3>,Fixed<3> >,double> SymMat3;

  typedef  Matrix<Symmetric<Var,Var>,double> SymMatV;

  typedef  Matrix<General<Var,Var>,double> MatV;

  typedef  Vector<General<Var,Fixed<1> >,double> VecV;

  typedef  Matrix<General<Fixed<3>, Var>,double> Mat3V;
}

#endif


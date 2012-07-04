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
  typedef  Matrix<General, double> Mat;

  /*! 
   *  \brief Vector
   *  
   *  Vec is an abbreviation for vectors of type double.
   * */
  typedef  Vector<General, double> Vec;

  /*! 
   *  \brief Row vector
   *  
   *  RowVec is an abbreviation for row vectors of type double.
   * */
  typedef  RowVector<General, double> RowVec;

  /*! 
   *  \brief Square matrix
   *  
   *  SqrMat is an abbreviation for square matrices of type double.
   * */
  typedef  SquareMatrix<General, double> SqrMat;

  /*! 
   *  \brief Diagonal matrix
   *  
   *  DiagMat is an abbreviation for diagonal matrices of type double.
   * */
  typedef  Matrix<Diagonal, double> DiagMat;

  /*! 
   *  \brief Symmetric matrix
   *  
   *  SymMat is an abbreviation for symmetric matrices of type double.
   * */
  typedef  Matrix<Symmetric, double> SymMat;

  /*! 
   *  \brief Band matrix
   *  
   *  BandMat is an abbreviation for band matrices of type double.
   * */
  typedef  Matrix<GeneralBand, double> BandMat;

  typedef  Vector<General, int> VecInt;

  typedef  Matrix<GeneralFixed<3,3>,double> FMat;
  
  typedef  Matrix<GeneralFixed<3,3>,double> Mat33;

  typedef  Matrix<GeneralFixed<3,2>,double> Mat32;

  typedef  SquareMatrix<GeneralFixed<3,3>,double> SqrMat3;

  typedef  SquareMatrix<GeneralFixed<3,3>,double> FSqrMat;

  typedef  Vector<GeneralFixed<3,1>,double> FVec;

  typedef  Vector<GeneralFixed<3,1>,double> Vec3;

  typedef  Vector<GeneralFixed<2,1>,double> Vec2;

  typedef  RowVector<GeneralFixed<1,3>,double> FRowVec;

  typedef  Matrix<SymmetricFixed<3>,double> FSymMat;

  typedef  Matrix<SymmetricVar,double> VSymMat;

  typedef  Matrix<GeneralVar,double> VMat;

  typedef  Vector<GeneralVar,double> VVec;

  typedef  Matrix<GeneralFixedVar<3>,double> FVMat;
}

#endif


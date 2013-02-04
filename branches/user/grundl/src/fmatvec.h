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

  /* Typedefs for Fixed types */

  //Column Vectors
  typedef Vector<Fixed<1>, double> Vec1;
  typedef Vector<Fixed<2>, double> Vec2;
  typedef Vector<Fixed<3>, double> Vec3;
  typedef Vector<Fixed<4>, double> Vec4;
  typedef Vector<Fixed<6>, double> Vec6;
  typedef Vector<Fixed<8>, double> Vec8;
  typedef Vector<Fixed<11>, double> Vec11;
  typedef Vector<Fixed<12>, double> Vec12;
  typedef Vector<Fixed<16>, double> Vec16;

  //Row Vectors
  typedef RowVector<Fixed<1>, double> RowVec1;
  typedef RowVector<Fixed<2>, double> RowVec2;
  typedef RowVector<Fixed<3>, double> RowVec3;
  typedef RowVector<Fixed<4>, double> RowVec4;

  typedef RowVector<Fixed<11>, double> RowVec11;
  typedef RowVector<Fixed<16>, double> RowVec16;

  //General Matrices
  typedef Matrix<General, Fixed<2>, Fixed<3>, double> Mat2x3;

  typedef Matrix<General, Fixed<3>, Fixed<2>, double> Mat3x2;
  typedef Matrix<General, Fixed<3>, Fixed<3>, double> Mat3x3;
  typedef Matrix<General, Fixed<3>, Fixed<6>, double> Mat3x6;
  typedef Matrix<General, Fixed<3>, Fixed<11>, double> Mat3x11;
  typedef Matrix<General, Fixed<3>, Fixed<16>, double> Mat3x16;

  typedef Matrix<General, Fixed<4>, Fixed<16>, double> Mat4x16;

  typedef Matrix<General, Fixed<8>, Fixed<3>, double> Mat8x3;
  typedef Matrix<General, Fixed<8>, Fixed<9>, double> Mat8x9;

  typedef Matrix<General, Fixed<11>, Fixed<3>, double> Mat11x3;
  typedef Matrix<General, Fixed<11>, Fixed<16>, double> Mat11x16;
  typedef Matrix<General, Fixed<11>, Fixed<27>, double> Mat11x27;

  typedef Matrix<General, Fixed<16>, Fixed<3>, double> Mat16x3;
  typedef Matrix<General, Fixed<16>, Fixed<4>, double> Mat16x4;
  typedef Matrix<General, Fixed<16>, Fixed<6>, double> Mat16x6;
  typedef Matrix<General, Fixed<16>, Fixed<8>, double> Mat16x8;

  //Square Matrices
  typedef SquareMatrix<Fixed<1>, double> SqrMat1;
  typedef SquareMatrix<Fixed<2>, double> SqrMat2;
  typedef SquareMatrix<Fixed<3>, double> SqrMat3;
  typedef SquareMatrix<Fixed<4>, double> SqrMat4;
  typedef SquareMatrix<Fixed<8>, double> SqrMat8;
  typedef SquareMatrix<Fixed<11>, double> SqrMat11;
  typedef SquareMatrix<Fixed<16>, double> SqrMat16;

  //Symmetric Matrices
  typedef Matrix<Symmetric, Fixed<3>, Fixed<3>, double> SymMat3;
  typedef Matrix<Symmetric, Fixed<4>, Fixed<4>, double> SymMat4;
  typedef Matrix<Symmetric, Fixed<8>, Fixed<8>, double> SymMat8;
  typedef Matrix<Symmetric, Fixed<16>, Fixed<16>, double> SymMat16;

  typedef Matrix<Symmetric, Var, Var, double> SymMatV;

  typedef Matrix<General, Var, Var, double> MatV;

  typedef SquareMatrix<Var, double> SqrMatV;

  typedef Vector<Var, double> VecV;

  typedef RowVector<Var, double> RowVecV;

  /*Typedefs for mixed size */

  typedef Matrix<General, Fixed<3>, Var, double> Mat3V;
  typedef Matrix<General, Var, Fixed<3>, double> MatV3;

  typedef Matrix<General, Fixed<11>, Var, double> Mat11V;
  typedef Matrix<General, Var, Fixed<11>, double> MatV11;

  typedef Matrix<General, Var, Fixed<16>, double> MatV16;
  typedef Matrix<General, Fixed<16>, Var, double> Mat16V;

  typedef Matrix<General, Ref, Fixed<3>, double> MatRef3;

  /*Typedefs for Ranges*/
  typedef Range<Fixed<0>, Fixed<1> > Range0x1;
  typedef Range<Fixed<0>, Fixed<2> > Range0x2;
  typedef Range<Fixed<0>, Fixed<4> > Range0x4;
  typedef Range<Fixed<0>, Fixed<9> > Range0x9;
  typedef Range<Fixed<3>, Fixed<5> > Range3x5;
  typedef Range<Fixed<3>, Fixed<10> > Range3x10;
  typedef Range<Fixed<5>, Fixed<7> > Range5x7;
  typedef Range<Fixed<6>, Fixed<8> > Range6x8;
  typedef Range<Fixed<7>, Fixed<14> > Range7x14;
  typedef Range<Fixed<9>, Fixed<11> > Range9x11;
  typedef Range<Fixed<10>, Fixed<12> > Range10x12;
  typedef Range<Fixed<10>, Fixed<15> > Range10x15;
  typedef Range<Fixed<13>, Fixed<15> > Range13x15;


  /*HELPER Function --> TODO: Move or better: avoid!!*/
  template<class Col>
  std::string getVectorTemplateType(Vector<Col, double> *instance) {
      if(dynamic_cast<Vector<Ref, double> *>(instance))
        return "Ref";
      else if(dynamic_cast<Vector<Fixed<8>, double> *>(instance)) {
        return "Fixed<8>";
      }
      else if(dynamic_cast<Vector<Fixed<16>, double> *>(instance)) {
        return "Fixed<16>";
      }
      else {
        throw "ERROR: getVectorTemplateType --> Not known Columns-Type";
      }

  }
}

#endif


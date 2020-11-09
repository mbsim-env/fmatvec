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

#include "stream_impl.h"
#include "general_matrix.h"
#include "var_general_matrix.h"
#include "fixed_var_general_matrix.h"
#include "var_fixed_general_matrix.h"
#include "fixed_general_matrix.h"
#include "symmetric_matrix.h"
#include "var_symmetric_matrix.h"
#include "fixed_symmetric_matrix.h"
#include "diagonal_matrix.h"
#include "ast.h"

namespace fmatvec {

  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Ref     ,Ref     ,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Ref     ,Ref     ,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Var     ,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Var     ,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Var     ,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Var     ,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<2>,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<2>,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<3>,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<3>,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<4>,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<4>,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<5>,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<5>,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<6>,Fixed<1>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<6>,Fixed<1>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<2>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<2>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<3>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<3>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<4>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<4>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<5>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<5>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<General  ,Fixed<1>,Fixed<6>,IndependentVariable > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<General  ,Fixed<1>,Fixed<6>,IndependentVariable > &A);
  template std::istream& operator>>(std::istream &s,       Matrix<Diagonal ,Ref     ,Ref     ,double              > &A);
  template std::ostream& operator<<(std::ostream &s, const Matrix<Diagonal ,Ref     ,Ref     ,double              > &A);

}

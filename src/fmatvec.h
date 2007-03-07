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
 *   mfoerg@users.berlios.de
 *
 */

#ifndef fmatrixvector_h
#define fmatrixvector_h

#include "linear_algebra.h"
#include "linear_algebra_double.h"

namespace fmatvec {
  typedef  Matrix<General, double> Mat;
  typedef  Vector<double> Vec;
  typedef  RowVector<double> RowVec;
  typedef  SquareMatrix<double> SqrMat;
  typedef  Matrix<Diagonal, double> DiagMat;
  typedef  Matrix<Symmetric, double> SymMat;
  typedef  Matrix<GeneralBand, double> BandMat;
}

#endif


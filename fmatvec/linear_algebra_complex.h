/* Copyright (C) 2003-2014  Martin Förg

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

#ifndef linear_algebraz_h
#define linear_algebraz_h

#include "square_matrix.h"
#include "vector.h"
#include "row_vector.h"
#include "diagonal_matrix.h"
#include "symmetric_matrix.h"
#include "band_matrix.h"
#include <complex>

//-------------------------------------
// Matrix operations
//-------------------------------------

namespace fmatvec {

  /*! \brief System of linear equations
   *
   * This function solves a complex system of linear equations 
   * according to \f[\boldsymbol{A}\,\boldsymbol{x}=\boldsymbol{b} \f]
   * by a LU decompostion.
   * \param A A complex square matrix. 
   * \param b A complex vector containing the right hand side.
   * \return A complex vector containig the solution.
   * */
  FMATVEC_EXPORT Vector<Ref, std::complex<double>> slvLU(const SquareMatrix<Ref, std::complex<double>> &A, const Vector<Ref, std::complex<double>> &b);

  // Because template argument deduction does not consider implicit conversions, the complex standard operators
  // cannot be used for mixed integer/complex arithmetic. Hence, we define these here for mixed
  // integer/complex arithmetic.
  template<class T> std::complex<T> operator+(const std::complex<T> &x, int y) { return x+static_cast<T>(y); }
  template<class T> std::complex<T> operator-(const std::complex<T> &x, int y) { return x+static_cast<T>(y); }
  template<class T> std::complex<T> operator*(const std::complex<T> &x, int y) { return x+static_cast<T>(y); }
  template<class T> std::complex<T> operator/(const std::complex<T> &x, int y) { return x+static_cast<T>(y); }
  template<class T> std::complex<T> operator+(int x, const std::complex<T> &y) { return static_cast<T>(x)+y; }
  template<class T> std::complex<T> operator-(int x, const std::complex<T> &y) { return static_cast<T>(x)+y; }
  template<class T> std::complex<T> operator*(int x, const std::complex<T> &y) { return static_cast<T>(x)+y; }
  template<class T> std::complex<T> operator/(int x, const std::complex<T> &y) { return static_cast<T>(x)+y; }

}

#endif

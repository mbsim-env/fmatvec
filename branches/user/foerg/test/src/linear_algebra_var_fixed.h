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

#ifndef linear_algebra_var_fixed_h
#define linear_algebra_var_fixed_h

#include "var_fixed_general_matrix.h"
#include "fixed_vector.h"
#include "var_vector.h"

namespace fmatvec {

  template <int N, class AT>
    inline VarVector<AT> operator*(const Matrix<VarFixedGeneral<N>, AT> &A, const VarVector<AT> &x) {
#ifndef FMATVEC_NO_SIZE_CHECK
      assert(N == x.size());
#endif
      VarVector<AT> y(A.rows(),NONINIT);
      for(int i=0; i<A.rows(); i++) {
	y.e(i) = 0;
	for(int j=0; j<N; j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }

  template <int N, class AT>
    inline VarVector<AT> operator*(const Matrix<VarFixedGeneral<N>, AT> &A, const FixedVector<N,AT> &x) {
      VarVector<AT> y(A.rows(),NONINIT);
      for(int i=0; i<A.rows(); i++) {
	y.e(i) = 0;
	for(int j=0; j<N; j++) 
	  y.e(i) += A.e(i,j)*x.e(j);
      }
      return y;
    }


}

#endif

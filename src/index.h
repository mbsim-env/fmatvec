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

#ifndef index_h
#define index_h

#include "matrix.h"

namespace fmatvec {

  /*! 
   * \brief This is an index class for creating submatrices.
   *
   * The index class contains indices defining the first and the last element.
   * */
  class Index {
    private:
      int i1, i2;
    public:

      /*! \brief Standard constructor
       *
       * Constructs an index object with both elements set to zero.
       */
      inline Index() { i1 = i2 = 0;}

      /*! \brief Regular constructor
       *
       * Constructs an index object with equal first and last element.
       * \param i1_ First and last element. 
       */
      inline Index(int i1_) : i1(i1_), i2(i1_) { 
#ifndef FMATVEC_NO_BOUNDS_CHECK
	assert(i1>=0);
#endif
      }

      /*! \brief Regular constructor
       *
       * Constructs an index object with first and last element.
       * \param i1_ First index. 
       * \param i2_ Last index. 
       */
      inline Index(int i1_, int i2_) : i1(i1_), i2(i2_) {
#ifndef FMATVEC_NO_BOUNDS_CHECK
	assert(i1>=0);
	assert(i2>=i1-1);
#endif
      }

      /*! \brief First element 
       *
       * \return Returns the first element. 
       */
      int start() const {return i1;}

      /*! \brief Last element 
       *
       * \return Returns the last element. 
       */
      int end() const {return i2;}
  };

  /*! \brief Equality operator for indices.
   *
   * Checks for equality of two indices.
   * \return true if two index objects are equal, false otherwise. 
   */
  bool operator==(const Index &I, const Index &J);

}

#endif

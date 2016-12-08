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

#ifndef range_h
#define range_h

//#include "matrix.h"
#include "types.h"

namespace fmatvec {

  /*! 
   * \brief This is an index class for creating submatrices.
   *
   * The index class contains indices defining the first and the last element.
   * */
  template <class Row, class Col> class Range {
    private:
  };

  /*! 
   * \brief This is an index class for creating submatrices.
   *
   * The index class contains indices defining the first and the last element.
   * */
  template <> class Range<Var, Var> {
    private:
      
   /// @cond NO_SHOW

    int i1, i2;

   /// @endcond

    public:

      /*! \brief Standard constructor
       *
       * Constructs an index object with both elements set to zero.
       */
      inline Range() { i1 = i2 = 0;}

      /*! \brief Regular constructor
       *
       * Constructs an index object with equal first and last element.
       * \param i1_ First and last element. 
       */
      inline Range(int i1_) : i1(i1_), i2(i1_) { 
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
      inline Range(int i1_, int i2_) : i1(i1_), i2(i2_) {
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

      /*! \brief Size
       *
       * \return Returns the Size. 
       */
      int size() const {return i2-i1+1;}
  };


  /*! 
   * \brief This is an index class for creating submatrices.
   *
   * The index class contains indices defining the first and the last element.
   * */
  template <int I1, int I2> class Range<Fixed<I1>, Fixed<I2> > {
    private:
      
    public:

      /*! \brief First element 
       *
       * \return Returns the first element. 
       */
      int start() const {return I1;}

      /*! \brief Last element 
       *
       * \return Returns the last element. 
       */
      int end() const {return I2;}
  };

  /*! \brief Equality operator for indices.
   *
   * Checks for equality of two indices.
   * \return true if two index objects are equal, false otherwise. 
   */
  inline bool operator==(const Range<Var,Var> &I, const Range<Var,Var> &J) {
    if(I.start() == J.start() && I.end() == J.end())
      return true;
    else
      return false;
  }

}

#endif

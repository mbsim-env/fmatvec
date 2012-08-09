/* Copyright (C) 2003-2005  Jan Clauberg

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
 *   jan.p.clauberg@googlemail.com
 *
 */

#ifndef RANDOMGENERATOR
#define RANDOMGENERATOR

#include "fmatvec.h"

using namespace std;

namespace fmatvec {

  class Randomgenerator {
  
    private:
      int Range;
      int DecimalPlaces;
      time_t t;

    public:
      /*! \brief Standard constructor
       *
       * Constructs a Randomgenerator. 
       * */
      Randomgenerator() {}

      /*! \brief Standard destructor
       *
       * Destructs a Randomgenerator. 
       * */
      ~Randomgenerator() {}
  
      /*! \brief Setter: range of random numbers
       *
       * Sets the range from ]-Range/2;Range/2[
       * \param Range_ range for random numbers.
       * */      
      void setRange(int Range_) {Range=Range_;}

      /*! \brief Setter: number of decimal places
       *
       * Sets the number of decimal places
       * \param DecimalPlaces_ number of decimal places.
       * */ 
      void setDecimalPlaces(int DecimalPlaces_) {DecimalPlaces=DecimalPlaces_;}
  
      /*! \brief Initializing the RandomGenerator
       *
       * */      
      void init() {time(&t); srand((unsigned int)t);}

      /*! \brief Operator: Generates a random number
       *
       * Generates a random number
       * \param norm_ norm of random number.
       * \return A random number
       * */       
      double operator()(bool norm_=false);
  
  };
}

#endif /* RANDOMGENERATOR */

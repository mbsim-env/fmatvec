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

#include <randomgenerator.h>
#include <config.h>
#include <iostream>
#include "math.h"

namespace fmatvec {


  double Randomgenerator::operator()(bool norm_) {
    
    if(norm_) {
      return ((rand() % Range) + static_cast<double>(rand() % static_cast<int>(pow(10.,DecimalPlaces)))/(pow(10.,DecimalPlaces)))/2.;
    }
    else {
      return (rand() % Range) + static_cast<double>(rand() % static_cast<int>(pow(10.,DecimalPlaces)))/(pow(10.,DecimalPlaces))-static_cast<double>(Range)/2.;
    }
  }
}

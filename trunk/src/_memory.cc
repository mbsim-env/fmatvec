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

#include "config.h"
#include "_memory.h"
#include <complex>

using namespace std;

namespace fmatvec {

  template <> Memory<char>::allocator_type Memory<char>::ms = Memory<char>::allocator_type();
  template <> Memory<int>::allocator_type Memory<int>::ms = Memory<int>::allocator_type();
  template <> Memory<double>::allocator_type Memory<double>::ms = Memory<double>::allocator_type();
  template <> Memory<complex<double> >::allocator_type Memory<complex<double> >::ms = Memory<complex<double> >::allocator_type();

}



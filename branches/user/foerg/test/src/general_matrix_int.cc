
/* Copyright (C) 2003-2005  Martin Förg, Robert Huber

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
#include "general_matrix.h"
#include "symmetric_matrix.h"
#include "blas_extensions_double.h"
#include <iostream>
#include <sstream>

using namespace std;

#define FMATVEC_NO_INITIALIZATION
#define FMATVEC_NO_BOUNDS_CHECK

namespace fmatvec {

  template <>
  Matrix<General<Ref>, int>::Matrix(const char *strs) {
    istringstream iss(strs);
    int x;
    char c;
    m = 0,n=0;
    int buf=0;
    iss >> c;
    do {
      iss >> c;
      if(c==';') {
	if(buf)
	  assert(buf == n);

	//if(buf && buf!=n)
	//	cout << "Mist" << endl;
	buf=n;
	n=0;
	m++;
      }
      else if(c==',')
	n++;
      c='0';
    } while(iss);

    n++; m++;
    lda=m;
    tp = false;
    memory.resize(m*n);
    ele = (int*)memory.get();
    iss.clear();
    iss.seekg(0);
    iss >> c;
    for(int i=0; i<m; i++)
      for(int j=0; j<n; j++) {
	iss >> x;
	operator()(i,j)=x;
	iss >> c;
      }
  }


}

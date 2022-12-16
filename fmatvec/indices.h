/* Copyright (C) 2022  Martin Förg

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

#ifndef indices_h
#define indices_h

#include <cassert>
#include <utility>
#include <vector>

namespace fmatvec {

  class Indices {
    public:
      Indices() = default;
      Indices(int size) { ind.resize(size); }
      Indices(std::vector<int> ind_) : ind(std::move(ind_)) {
	assert(min()>=0);
      }
      Indices(std::initializer_list<int> il) : ind(il) {
	assert(min()>=0);
      }
      int size() const { return ind.size(); }
      int max() const {
	assert(ind.size());
	int m = ind[0];
	for(size_t i=1; i<ind.size(); i++) {
	  if(ind[i]>m)
	    m = ind[i];
	}
	return m;
      }
      int min() const {
	assert(ind.size());
	int m = ind[0];
	for(size_t i=1; i<ind.size(); i++) {
	  if(ind[i]<m)
	    m = ind[i];
	}
	return m;
      }
      void add(int ind_) { 
	assert(ind_>=0);
	ind.push_back(ind_);
      }
      void set(int i, int ind_) {
	assert(ind_>=0);
       	ind[i] = ind_;
      }
      const int& operator[](int i) const { return ind[i]; }
    private:
      std::vector<int> ind;
  };

}

#endif

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

#ifndef memory_h
#define memory_h

#define MSTACKSIZE  20
#define MAXSZ 10000

#include <iostream>

using namespace std;

namespace fmatvec {

 /// @cond NO_SHOW

  template <class AT> class MemoryStack {
      private:
	AT ***ele;
      public:
	MemoryStack() {
	  ele=new AT**[MAXSZ];
	  for(int i=0; i<MAXSZ; i++)
	    ele[i]=0;
	}
	~MemoryStack() {
	  for(int i=0; i<MAXSZ; i++) {
	    AT **&sp=ele[i];
	    if(sp) {
	      while(*sp) {
		sp--;
	      }
	      for(int j=0; j<MSTACKSIZE; j++)
		if(sp[j]!= (AT *) 1)
		  delete [] ele[i][j];
	      delete [] sp;
	    }
	  }
	  delete [] ele; 
	}

	AT** renewstack(AT ** sp) {
	  AT **sp_neu,**sp_merk=sp;
	  int alte_laenge,neue_laenge,i;
	  while(*(--sp)); 
	  alte_laenge=sp_merk-sp+1;
	  neue_laenge=alte_laenge+MSTACKSIZE;
	  sp_neu = new AT*[neue_laenge];
	  for(i=0;i<alte_laenge-1;i++) sp_neu[i]=sp[i];
	  delete[] sp;
	  for(i=alte_laenge;i<neue_laenge-1;i++) sp_neu[i] = (AT *) 1;
	  sp_neu[neue_laenge-1]=NULL;
	  return sp_neu+alte_laenge-1;
	}

	AT** newstack() {
	  AT** sp;
	  AT** stackptr;
	  sp = stackptr = new AT*[MSTACKSIZE];
	  for (int i=0; i<MSTACKSIZE; i++) *(sp++) = (AT *) 1;
	  stackptr[0]=NULL; stackptr[MSTACKSIZE-1]=NULL;
	  return stackptr;
	}

	AT * getMemory(size_t sz) { 
	  AT **&sp=ele[sz];
	  if(sz<MAXSZ) {
	    if(!sp)
	      sp = newstack();
	    if (*sp) 
	      return *sp--;
	    else 
	      return new AT[sz]; 
	  } else
	    return new AT[sz];  
	}
	void freeMemory(size_t sz, AT *p) {
	  AT **&sp=ele[sz];
	  if(sz<MAXSZ && sz) {
	    if(*(++sp)) 
	      *sp = p; 
	    else {
	      sp=renewstack(sp);
	      *sp=p;
	    }
	  } else delete [] p;
	}
    };

  template <class AT> class Memory {
      private:
	size_t sz;
	AT *ele0;
	static MemoryStack<AT> ms;
	size_t *ref;
	void lock() {(*ref)++;};
	void unlock() {
	  if(!--(*ref)) {
	    ms.freeMemory(sz,ele0);
	    delete ref;
	  }
	};

      public:
	Memory() : sz(0), ele0(0), ref(new size_t(1)){ 
	};

	Memory(size_t n) : sz(n), ele0(ms.getMemory(sz)), ref(new size_t(1)){
	};

	Memory(const Memory &memory) : sz(memory.sz), ele0(memory.ele0), ref(memory.ref)  {
	  lock();
	};

	~Memory() {
	  unlock();
	}; 

	Memory& operator=(const Memory &memory) {
	  if(this == &memory)
	    return *this;
	  unlock(); sz=memory.sz; ref=memory.ref; ele0=memory.ele0; lock();
	  return *this;
	}

	void resize(size_t n) {
	  unlock();
	  ref=new size_t(1);
	  sz = n;
	  ele0 = ms.getMemory(sz); 
	};

	AT* get() const {return ele0;};
    };
}

 /// @endcond
 
#endif

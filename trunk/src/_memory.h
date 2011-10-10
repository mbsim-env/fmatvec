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

#ifndef _memory_h
#define _memory_h

// includes
#include <stdlib.h>
#include <vector>

// include header given by ALLOCATORHEADER
#ifdef ALLOCATORHEADER
#  define CPPBEGININCLUDE <
#  define CPPENDINCLUDE >
#  include CPPBEGININCLUDE./ALLOCATORHEADER CPPENDINCLUDE
#  undef CPPBEGININCLUDE
#  undef CPPENDINCLUDE
#endif
// use fmatvec::PoolAllocator if no allocator is defined by ALLOCATORCLASS
#ifndef ALLOCATORCLASS
#  define ALLOCATORCLASS fmatvec::PoolAllocator
#endif

// create some locking and atomic operations dependent on the thread safety method
//
// FMATVEC_LOCKVAR(l) : create a locking/mutex named l
// FMATVEC_LOCKINIT(l) : initialize the locking/mutex named l
// FMATVEC_LOCK(l) : lock the locking/mutex l
// FMATVEC_UNLOCK(l) : unlock the locking/mutex l
//
// FMATVEC_SYNCVAR(l) : create a locking/mutex named l for a syncronized operation
// FMATVEC_SYNCINIT(l) : initialize a locking/mutex named l for a syncronized operation
// FMATVEC_SYNCINC(var, l) : var++ (syncronized by l)
// FMATVEC_SYNCPREDEC(ret, var, l) : ret=--(var) (syncronized by l)
#if defined FMATVEC_THREADSAFE_GCCBUILTIN
  #define FMATVEC_LOCKVAR(l) short int l;
  #define FMATVEC_LOCKINIT(l) l=0;
  #define FMATVEC_LOCK(l) while(__sync_bool_compare_and_swap(&l, 0, 0xffff)==false);
  #define FMATVEC_UNLOCK(l) __sync_bool_compare_and_swap(&l, 0xffff, 0);

  #define FMATVEC_SYNCVAR(l)
  #define FMATVEC_SYNCINIT(l)
  #define FMATVEC_SYNCINC(var, l) __sync_add_and_fetch(&(var), 1);
  #define FMATVEC_SYNCPREDEC(ret, var, l) ret=__sync_sub_and_fetch(&(var), 1);
#elif defined FMATVEC_THREADSAFE_PTHREAD
  #include <pthread.h> // pthread required a header file

  #define FMATVEC_LOCKVAR(l) pthread_mutex_t l;
  #define FMATVEC_LOCKINIT(l) pthread_mutex_init(&l, NULL);
  #define FMATVEC_LOCK(l) pthread_mutex_lock(&l);
  #define FMATVEC_UNLOCK(l) pthread_mutex_unlock(&l);

  #define FMATVEC_SYNCVAR(l)  pthread_mutex_t l;
  #define FMATVEC_SYNCINIT(l) pthread_mutex_init(&l, NULL);
  #define FMATVEC_SYNCINC(var, l) pthread_mutex_lock(&l); (var)++; pthread_mutex_unlock(&l);
  #define FMATVEC_SYNCPREDEC(ret, var, l) pthread_mutex_lock(&l); ret=--(var); pthread_mutex_unlock(&l);
#elif defined FMATVEC_THREADSAFE_OPENMP
  #define FMATVEC_PRAGMA(x) _Pragma(#x) // openmp required the _Pragma operator (defined by C99)

  #define FMATVEC_LOCKVAR(l)
  #define FMATVEC_LOCKINIT(l)
  #define FMATVEC_LOCK(l) FMATVEC_PRAGMA(omp critical (l)) {
  #define FMATVEC_UNLOCK(l) }

  #define FMATVEC_SYNCVAR(l)
  #define FMATVEC_SYNCINIT(l)
  #define FMATVEC_SYNCINC(var, l) FMATVEC_PRAGMA(omp critical (l)) { (var)++; }
  #define FMATVEC_SYNCPREDEC(ret, var, l) FMATVEC_PRAGMA(omp critical (l)) { ret=--(var); }
#else
  #define FMATVEC_PRAGMA(l)
  #define FMATVEC_LOCKVAR(l)
  #define FMATVEC_LOCKINIT(l)
  #define FMATVEC_LOCK(l)
  #define FMATVEC_UNLOCK(l)

  #define FMATVEC_SYNCVAR(l)
  #define FMATVEC_SYNCINIT(l)
  #define FMATVEC_SYNCINC(var, l) (var)++;
  #define FMATVEC_SYNCPREDEC(ret, var, l) ret=--(var);
#endif

/**
 * \brief memory managment and reference counting 
 * \author Martin Foerg
 * \date 2010-07-07 pragma omp critical can be deactivated by --disable-pragma_omp_critical (Robert Huber)
 */

namespace fmatvec {

 /// @cond NO_SHOW

  // A STL-like pool memory allocator: fast but NOT thread safe till now!
  template <class AT>
  class PoolAllocator {
    private:
      // The minimal allocation size in power of 2: MINSIZE=2^MINEXP
      static const size_t MINEXP=4;
      // A pool of memory of allocation size 2^(MINEXP+0), 2^(MINEXP+1), 2^(MINEXP+2), ...
      // stored at index 0, 1, 2, ... in memoryPool.
      // The pointers stored here are already allocated but currently not used.
      std::vector<std::vector<AT*> > memoryPool;
      // Just a dummy variable used for zero size allocations.
      AT sizeZero;
      FMATVEC_LOCKVAR(fmatvec_memoryPoolLock)
    public:
      // Constructor: alloc memoryPool with at least size 1, since the size is doubled each time the size ran out.
      PoolAllocator() : memoryPool(1) {
        FMATVEC_LOCKINIT(fmatvec_memoryPoolLock)
      }
      // Destructor: delete all cached memory in the pool.
      ~PoolAllocator() {
        // loop over all pool sizes and free memory
        for(typename std::vector<std::vector<AT*> >::iterator i=memoryPool.begin(); i!=memoryPool.end(); ++i)
          for(typename std::vector<AT*>::iterator j=i->begin(); j!=i->end(); ++j)
            delete[](*j);
      }
      // Allocate new memory of at least size size.
      AT *allocate(size_t size) { 
        // do nothing for zero size, but return a valid memory pointer
        if(size==0) return &sizeZero;
        // calculate size index: size 0 to 2^MINEXP = index 0; size 2^MINEXP+1 to 2*2^MINEXP = index 1; ...
        size_t index=0;
        --size>>=MINEXP-1;
        while(size>>=1) ++index;
        AT* ret=NULL;
        FMATVEC_LOCK(fmatvec_memoryPoolLock)
          // increase (double) memoryPool if nessasary
          if(index>=memoryPool.size()) memoryPool.resize(index<<1);
          // try to use cached memory
          std::vector<AT*> &memoryPoolIndex=memoryPool[index];
          if(memoryPoolIndex.size()>0) {
            ret=memoryPoolIndex.back(); // get last cached memory
            memoryPoolIndex.pop_back(); // remove last cached memory
          }
        FMATVEC_UNLOCK(fmatvec_memoryPoolLock)
        if(ret) return ret; // return chaced memory
        // return new memory
        return new AT[1<<(index+MINEXP)]; // allocate 2^(index+MINEXP) elements (max size for this index)
      }
      // Free memory and release it to the pool.
      void deallocate(AT *p, size_t size) {
        // do nothing for zero size, a dummy but valid pointer was returned by allocate
        if(size==0) return;
        // calculate size index: size 0 to 2^MINEXP = index 0; size 2^MINEXP+1 to 2*2^MINEXP = index 1; ...
        size_t index=0;
        --size>>=MINEXP-1;
        while(size>>=1) ++index;
        FMATVEC_LOCK(fmatvec_memoryPoolLock)
          // release memory
          memoryPool[index].push_back(p); // add to cached memory
        FMATVEC_UNLOCK(fmatvec_memoryPoolLock)
      }
  };

  template <class AT> class Memory {
      private:
        typedef ALLOCATORCLASS<AT> allocator_type;
	size_t sz;
	AT *ele0;
        static allocator_type ms;
	size_t *ref;
        FMATVEC_SYNCVAR(fmatvec_refLock)
	void lock() {
          FMATVEC_SYNCINC(*ref, fmatvec_refLock)
        };
	void unlock() {
          size_t refLocal;
          FMATVEC_SYNCPREDEC(refLocal, *ref, fmatvec_refLock)
	  if(!refLocal) {
            ms.deallocate(ele0, sz);
	    delete ref;
	  }
	};

      public:
	Memory() : sz(0), ele0(0), ref(new size_t(1)) {
          FMATVEC_SYNCINIT(fmatvec_refLock);
	};

	Memory(size_t n) : sz(n), ele0(ms.allocate(sz)), ref(new size_t(1)){
          FMATVEC_SYNCINIT(fmatvec_refLock);
	};

	Memory(const Memory &memory) : sz(memory.sz), ele0(memory.ele0), ref(memory.ref)  {
          FMATVEC_SYNCINIT(fmatvec_refLock);
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
	  ele0 = ms.allocate(sz);
	};

	AT* get() const {return ele0;};
    };
}

 /// @endcond
 
#endif

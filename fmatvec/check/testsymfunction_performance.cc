#define FMATVEC_LRUCACHE_SMALLSIZE_WORSTPERFORMANCE
#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#  undef __STRICT_ANSI__ // to define _controlfp which is not part of ANSI and hence not defined in mingw
#  include <cfloat>
#  define __STRICT_ANSI__
#endif
#include <cfenv>
#include <cassert>
#include "fmatvec/symbolic.h"
#include "fmatvec/lrucache.h"
#include <fmatvec/symbolic_function.h>
#include <boost/lexical_cast.hpp>
#include <chrono>

using namespace std;
using namespace fmatvec;

int main() {
#ifdef _WIN32
  SetConsoleCP(CP_UTF8);
  SetConsoleOutputCP(CP_UTF8);
  _controlfp(~(_EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW), _MCW_EM);
#else
  assert(feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)!=-1);
#endif

#ifdef NDEBUG
  constexpr int N=100000000;
#else
  constexpr int N=10;
#endif

  auto func=[](double x) {
    // the length of the expression is parameter 1 for the performance test
    return sin(1*x)+sin(2*x)+sin(3*x)+sin(4*x)+sin(5*x);
  };
  // the LRUCache size is parameters 2 for the performance test
  size_t cacheSize=1;
  // the small LRUCache size is parameters 3 for the performance test
  // (its default 45 is a good value; for smaller size the std::array approach is faster for larger size the std::unordered_map is faster)
  size_t cacheSmallSize=45;

  cout<<setprecision(12);
  Atom::setCurrentMessageStream(Atom::Debug);

  // performance test of LRUCache
  double allNoCache=0;
  {
    auto start=chrono::high_resolution_clock::now();
    for(size_t i=0; i<N; ++i) {
      double a=static_cast<double>(i)/55;
      auto v=func(a);
      allNoCache+=v;
    }
    auto end=chrono::high_resolution_clock::now();
    cout<<"NOCACHE  t="<<static_cast<chrono::duration<double>>(end-start).count()<<" (reference)"<<endl;
  }
  double allNoMatchCache=0;
  {
    LRUCache<double, double> cache(cacheSize, cacheSmallSize);
    cache.setName("NOMATCH");
    auto start=chrono::high_resolution_clock::now();
    for(size_t i=0; i<N; ++i) {
      double a=static_cast<double>(i)/55;
      auto [v,c]=cache(a);
      if(c)
        v=func(a);
      allNoMatchCache+=v;
    }
    auto end=chrono::high_resolution_clock::now();
    cout<<"NOMATCH  t="<<static_cast<chrono::duration<double>>(end-start).count()<<" (cache overhead, compared to NOCACHE)"<<endl;
  }
  double allAllMatchCache=0;
  {
    LRUCache<double, double> cache(cacheSize, cacheSmallSize);
    cache.setName("ALLMATCH");
    auto start=chrono::high_resolution_clock::now();
    for(size_t i=0; i<N; ++i) {
      double a=static_cast<double>(3847)/55;
      auto [v,c]=cache(a);
      if(c)
        v=func(a);
      allAllMatchCache+=v;
    }
    auto end=chrono::high_resolution_clock::now();
    cout<<"ALLMATCH t="<<static_cast<chrono::duration<double>>(end-start).count()<<" (should be faster than NOCACHE, else cache is useless)"<<endl;
  }
//1 no      cache t=5.748068602 157.199032262
//1 new     cache t=5.914504742 157.199032262
//1 Cache hit rate 0%
//1 exiting cache t=0.357286128 131343020.212
//1 Cache hit rate 100%
//
//> no      cache t=6.465329132 157.199032262
//> new     cache t=16.056331865 157.199032262
//> Cache hit rate 0%
//> exiting cache t=4.735711698 131343020.212
//> Cache hit rate 100%

  return allNoMatchCache!=allNoCache;
}

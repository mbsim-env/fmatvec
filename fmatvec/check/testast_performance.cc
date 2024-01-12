#include <chrono>
#include <iostream>
#include "fmatvec/symbolic.h"

using namespace std;
using namespace fmatvec;

int main() {

  cout<<setprecision(12);
#ifdef NDEBUG
  constexpr int N=100000000;
#else
  constexpr int N=10;
#endif

  double rSumSym=0;
  double rSumNative=0;
  {
    IndependentVariable a;
    IndependentVariable b;
    IndependentVariable c;
    SymbolicExpression r=3*a+sin(b)*pow(c,3)+fmatvec::min(a,c)-condition(2*a+5, atan2(4*a, b*c), cos(4*b+c));
    Eval rEval{r};
    auto start=std::chrono::high_resolution_clock::now();
    for(int i=0; i<N; ++i)
    {
      a^=4.7+static_cast<double>(i)/N;
      b^=8.3+static_cast<double>(i)/N;
      c^=2.9+static_cast<double>(i)/N;
      rSumSym+=rEval()/N;
    }
    auto end=std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> delta=end-start;
    cout<<"r.eval.SUM = "<<rSumSym<<endl;
    cout<<"execution time of symbolic performance test with N="<<N<<": "<<delta.count()<<"sec (about factor 2 slower)"<<endl;
  }
  {
    double a;
    double b;
    double c;
    auto rEval=[&a,&b,&c]() -> double {
      return 3*a+sin(b)*pow(c,3)+min(a,c)-(2*a+5>0 ? atan2(4*a, b*c) : cos(4*b+c));
    };
    auto start=std::chrono::high_resolution_clock::now();
    for(int i=0; i<N; ++i)
    {
      a=4.7+static_cast<double>(i)/N;
      b=8.3+static_cast<double>(i)/N;
      c=2.9+static_cast<double>(i)/N;
      rSumNative+=rEval()/N;
    }
    auto end=std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> delta=end-start;
    cout<<"r.eval.SUM = "<<rSumNative<<endl;
    cout<<"execution time of native performance test with N="<<N<<": "<<delta.count()<<"sec (about factor 2 faster)"<<endl;
  }

  return rSumSym != rSumNative;
}

#include <chrono>
#include <iostream>
#include "fmatvec/symbolic.h"

using namespace std;
using namespace fmatvec;

int main() {
  IndependentVariable a;
  IndependentVariable b;
  IndependentVariable c;
  SymbolicExpression r=3*a+sin(b)*pow(c,3)+min(a,c)-condition(2*a+5, cos(4*b+c), atan2(4*a, b*c));
  Eval rEval{r};
  double rSum=0;
#ifdef NDEBUG
  constexpr int N=100000000;
#else
  constexpr int N=10;
#endif
  auto start=std::chrono::high_resolution_clock::now();
  for(int i=0; i<N; ++i)
  {
    a^=4.7+static_cast<double>(i)/N;
    b^=8.3+static_cast<double>(i)/N;
    c^=2.9+static_cast<double>(i)/N;
    rSum+=rEval()/N;
  }
  auto end=std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> delta=end-start;
  cout<<"r.eval.SUM = "<<rSum<<endl;
  cout<<"execution time of performance test with N="<<N<<": "<<delta.count()<<"sec"<<endl;
}


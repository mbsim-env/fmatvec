#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  undef __STRICT_ANSI__ // to define _controlfp which is not part of ANSI and hence not defined in mingw
#  include <cfloat>
#  define __STRICT_ANSI__
#endif
#include <cfenv>
#include <cassert>
#include <iostream>

#include "fmatvec/linear_algebra_double.h"


//#include "fmatvec/fmatvec.h"

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

  array<double,4> in{2.0, 30.0, 20.0, 4.0};
  const SquareMatrix<Ref, double> A(2, in.data());
  Vector<Ref, std::complex<double>> eigenValues = eigval(A);

  // https://www.arndt-bruenner.de/mathe/scripts/engl_eigenwert2.htm
  // (-21.5153,0)
  // (27.5153, 0)
  std::cout << "Example eigval: " << std::endl;
  for (const complex<double> ev : eigenValues)
  {
      std::cout << ev << std::endl;
  }

  return 0;  
}

#include <cfenv>
#include <cassert>
#include <iostream>

#include "fmatvec/linear_algebra_double.h"


//#include "fmatvec/fmatvec.h"

using namespace std;
using namespace fmatvec;


int main() {
#ifndef _WIN32
  assert(feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)!=-1);
#endif

  double in[4]{2.0, 30.0, 20.0, 4.0};
  const SquareMatrix<Ref, double> A(2, in);
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

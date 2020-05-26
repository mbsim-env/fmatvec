#include <complex>

#include "fmatvec/linear_algebra_double.h"
#include "wrapper.h"

using namespace std;
using namespace fmatvec;

void eigval(std::size_t s, double* in, double* out) {
  const SquareMatrix<Ref, double> A(static_cast<int>(s), in);
  Vector<Ref, std::complex<double>> eigenValues = eigval(A);
  
  for (const complex<double> ev : eigenValues)
  {
      *out = ev.real();
      out++;
  }
}

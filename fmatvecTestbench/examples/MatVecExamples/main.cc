#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  define NOMINMAX
#  include <windows.h>
#  undef __STRICT_ANSI__ // to define _controlfp which is not part of ANSI and hence not defined in mingw
#  include <cfloat>
#  define __STRICT_ANSI__
#endif
#include <config.h>
#include <cassert>
#include <cfenv>
#include "fmatvec.h"
#include <iostream>

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

  cout << "----------------------------------------"<< endl;
  cout << "create an empty 3x2 matrix" << endl<<endl;
  Mat A1(3,2); 
  cout << "A1 = " << A1 << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "create a 3x2 matrix from string" << endl<<endl;
  Mat A2("[1.3, -0.2; 2.4, 0; 4, -1.9]"); 
  cout << "A2 = " << A2 << endl<<endl;
  
  cout << "----------------------------------------"<< endl;
  cout << "initialize A2 to 1" << endl<<endl;
  A2.init(1);
  cout << "A2 = " << A2 << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "copy constructor (A3 references A2)" << endl<<endl;
  Mat A3(A2); // same as Mat A3 = A2;
  A2(2,1) = 10;
  cout << "A2 = " << A2 << endl<<endl;
  cout << "A3 = " << A3 << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "assignment operator (deep copy)" << endl<<endl;
  Mat A4;
  A4 = A3;  // same as A4 << A3
  A4(0,0) = -5;
  cout << "A3 = " << A3 << endl<<endl;
  cout << "A4 = " << A4 << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "reference operator (A5 references A4)" << endl<<endl;
  Mat A5;
  A5 >> A4; 
  A5(0,1) = 0;
  cout << "A4 = " << A4 << endl<<endl;
  cout << "A5 = " << A5 << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "create a 2x2 submatrix of A2" << endl<<endl;
  Mat A6 = A2(Range<Var,Var>(1,2), Range<Var,Var>(0,1));
  cout << "A6 = " << A6 << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "edit A2 by A6" << endl<<endl;
  A6.init(-2);
  cout << "A2 = " << A2 << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "resize A1 to a 3x3 matrix" << endl<<endl;
  A1.resize(3,3);
  cout << "A1 = " << A1 << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "resize A1 to a dimensionless matrix" << endl<<endl;
  A1.resize();
  cout << "A1 = " << A1 << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "A1 is dimensionless and can be assigned again" << endl<<endl;
  A1 = A6;
  cout << "A1 = " << A1 << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "create empty vector" << endl<<endl;
  Vec x(3); 
  cout << "x = " << x << endl<<endl;
  cout << "----------------------------------------"<< endl;

  cout << "transpose of a matrix" << endl<<endl;
  SqrMat A("[1, -2, 0; 0, 1, 4; -1, 3, -1"); 
  cout << "A = " << A << endl<<endl;
  cout << "A' = " << A.T() << endl<<endl;
  cout << "trans(A) = " << trans(A) << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "matrix-vector product" << endl<<endl;
  x = "[1;2;1]"; 
  cout << "A = " << A << endl<<endl;
  cout << "x = " << x << endl<<endl;
  cout << "A*x = " << A*x << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "matrix-matrix addition" << endl<<endl;
  x = "[1;2;1]"; 
  cout << "A = " << A << endl<<endl;
  cout << "A+A = " << A+A << endl<<endl;
  
  cout << "----------------------------------------"<< endl;
  cout << "matrix-matrix product" << endl<<endl;
  x = "[1;2;1]"; 
  cout << "A = " << A << endl<<endl;
  cout << "A*A = " << A*A << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "solve linear system of equations" << endl<<endl;
  Vec b = "[0;1;2]"; 
  cout << "A = " << A << endl<<endl;
  cout << "b = " << b << endl<<endl;
  cout << "A*x = b, x = " << slvLU(A,b) << endl<<endl;
  cout << "inv(A)*x = " << inv(A)*b << endl<<endl;

  cout << "----------------------------------------"<< endl;
  cout << "eigenvalues of a matrix" << endl<<endl;
  cout << "A = " << A << endl<<endl;
  cout << "lambda(A) = " << eigval(A) << endl<<endl;

  return 0;
}

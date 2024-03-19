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
#include "fmatvecTestbench/testbench.h"
#include <iostream>

using namespace std;
using namespace fmatvec;

int main (int argc, char* argv[])
{
#ifdef _WIN32
  SetConsoleCP(CP_UTF8);
  SetConsoleOutputCP(CP_UTF8);
  _controlfp(~(_EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW), _MCW_EM);
#else
  assert(feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)!=-1);
#endif

  // Creating Testbench
  Testbench *testbench = new Testbench();
  
  // setting range for random numbers ]-Range/2;Range/2[  
  testbench->setRange(50);
  
  // setting number of decimal places
  testbench->setDecimalPlaces(4);
  
  // setting dimension of operation
  testbench->setDim(4);

  // initializing testbench
  testbench->init_Matgenerator();
  testbench->init_MatVec();

  // set number of runs
  testbench->setNumRuns(1);

  cout << "Testbench created and initialized" << endl;

  return 0;
}


#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
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

/*
includes from mbsim !!!
*/
#include "mbsim/utils/stopwatch.h"


using namespace std;
using namespace fmatvec;
using namespace MBSim;

int main (int argc, char* argv[])
{
#ifdef _WIN32
  SetConsoleCP(CP_UTF8);
  SetConsoleOutputCP(CP_UTF8);
  _controlfp(~(_EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW), _MCW_EM);
#else
  assert(feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)!=-1);
#endif

  // global parameter
  int NumRuns = 1;
  int DimEnd = 500;
  
  // Creating Timer
  StopWatch Timer1,Timer2;

  // Creating Testbench
  Testbench *testbench = new Testbench();
  testbench->setRange(1000);
  testbench->setDecimalPlaces(2);
  testbench->setNumRuns(NumRuns);
  testbench->init_Matgenerator();

  // variables for calculation
  int Dim=0;
  int Zaehler=0;
  double time=0;

  // Matrix for results (size statically 100)
  Mat Erg(200,2,INIT,0.0);

  do {

/////////////////////////////////////////////////////////
// testing individual function, e.g. operate_GenMatGenMat

    // setting dimension
    if (Dim<50) Dim=Dim+1;
    else if (Dim>=50 && Dim<100) Dim=Dim+2;
    else if (Dim>=100 && Dim<200) Dim=Dim+3;
    else if (Dim>=200 && Dim<300) Dim=Dim+20;
    else if (Dim>=300 && Dim<400) Dim=Dim+40;
    else Dim=Dim+60;

    // setting the dimension and initializing the Testbench
    testbench->setDim(Dim);
    //testbench->init_MatVec();

    Timer1.start();
      testbench->init_MatVec();
    cout << "Dimension was: " << Dim << " and init-time was: " << Timer2.stop() << " sec!" << endl;

    Timer2.start();
      testbench->operate_GenMatGenMat();
      //testbench->operate_SymMatSymMat();
      //testbench->operate_SqrMatSqrMat();
      //testbench->operate_DiagMatDiagMat();
      //testbench->operate_VecVec();
      //testbench->operate_inv();
      //testbench->operate_eigval();
      //testbench->operate_facLL();
      //testbench->operate_slvLLFac();
      //testbench->operate_slvLLFac_MatMat();
      //testbench->operate_slvLU();
      //testbench->operate_slvLU_MatMat();
    time = Timer2.stop(true)/NumRuns;

    cout << "Dimension was: " << Dim << " and time was: " << time << " sec!" << endl;
    Erg(Zaehler,0)=Dim;
    Erg(Zaehler,1)=time;   

    Zaehler++;
  } while (Dim < DimEnd);  
/////////////////////////////////////////////////////////


  cout << Erg << endl;

  return 0;
}


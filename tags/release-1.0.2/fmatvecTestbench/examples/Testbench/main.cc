#include "fmatvec.h"
#include "fmatvecTestbench/testbench.h"
#include <iostream>

using namespace std;
using namespace fmatvec;

int main (int argc, char* argv[])
{

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


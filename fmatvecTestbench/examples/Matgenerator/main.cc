#include <config.h>
#include <cassert>
#include <cfenv>
#include "fmatvec.h"
#include "fmatvecTestbench/matgenerator.h"
#include <iostream>

using namespace std;
using namespace fmatvec;

int main (int argc, char* argv[])
{
#ifndef _WIN32
  assert(feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW)!=-1);
#endif

  // Creating Matgenerator
  Matgenerator *matgenerator = new Matgenerator();
  
  // setting range for random numbers ]-Range/2;Range/2[  
  matgenerator->setRange(5000);
  
  // setting number of decimal places
  matgenerator->setDecimalPlaces(4);
  
  // initializing random generator
  matgenerator->init();

  cout << "Mat(3,3) " << matgenerator->getMat(3,3) << endl;
  cout << "SymMat(3,3) " << matgenerator->getSymMat(3,3) << endl;
  cout << "SqrMat(3,3) " << matgenerator->getSqrMat(3) << endl;
  cout << "DiagMat(4) " << matgenerator->getDiagMat(4) << endl;
  cout << "Vec(4) " << matgenerator->getVec(4) << endl;
  cout << "RowVec(4) " << matgenerator->getRowVec(4) << endl;
  cout << "hposdefSymMat(4) " << matgenerator->hposdefSymMat(7) << endl;
  cout << "hposdefSymFacLLMat(4) " << matgenerator->hposdefSymFacLLMat(4) << endl;


  return 0;
}


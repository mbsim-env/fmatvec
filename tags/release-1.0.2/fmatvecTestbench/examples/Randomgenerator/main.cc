#include "fmatvec.h"
#include "fmatvecTestbench/randomgenerator.h"
#include <iostream>

using namespace std;
using namespace fmatvec;

int main (int argc, char* argv[])
{

  // Creating Randomgenerator
  Randomgenerator *generator = new Randomgenerator();

  // setting range for random numbers ]-Range/2;Range/2[
  generator->setRange(500);

  // setting number of decimal places
  generator->setDecimalPlaces(3);

  // initializing random generator
  generator->init();

  for(int i=0; i<100; i++) {
    cout << "Zufallszahl= " << (*generator)() << endl;
  }

  return 0;
}


//#############################################################################
//
// ProtonXsec.cpp
//
//  A desription of what this class does
//
// January 17, 2018
//
//#############################################################################

#include "ProtonXsec.h"
#include "../Selection/EventSelector.h"
#include "../Selection/BeamSelector.h"

#include <iostream>


ProtonXsec::ProtonXsec(){}


int ProtonXsec::testFunction(){

  std::cout << "I calculate the p-Ar cross section! \n";
  
  EventSelector *ES = new EventSelector();
  std::cout << ES->classifyEvent( 4 ) << std::endl;
  std::cout << ES->classifyEvent( 7 ) << std::endl;
  return 0;

}


int main(){

  ProtonXsec protonXsec;
  protonXsec.testFunction();
  
  return 0;
  
}

//#############################################################################
//
// ProtonAnalyzerMC.cpp
//
//  A desription of what this class does
//
// January 17, 2018
//
//#############################################################################

#include "ProtonAnalyzerMC.h"
#include "../Selection/EventSelector.h"

#include <iostream>


ProtonAnalyzerMC::ProtonAnalyzerMC(){}


int ProtonAnalyzerMC::testFunction(){

  std::cout << "I analyze protons in LArIAT! \n";
  
  EventSelector *ES = new EventSelector();
  std::cout << ES->classifyEvent( 4 ) << std::endl;
  std::cout << ES->classifyEvent( 7 ) << std::endl;
  return 0;

}


int main(){

  ProtonAnalyzerMC protonAna;
  protonAna.testFunction();
  
  return 0;
  
}

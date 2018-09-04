
#ifndef BENDSTUDY_H
#define BENDSTUDY_H


#include "../Utilities/LArIATAnalysis.h"
#include "../Utilities/UtilityFunctions.h"
#include "math.h"
#include <stdio.h>
#include <string.h>



class BendStudy : public LArIATAnalysis {

 public:

  BendStudy( );
  BendStudy( char* jobOptionsFile );
 
  //------------------------------------
  // function signatures
  //------------------------------------

  void AnalyzeFromNtuples( );


  //TH1D *BeamSelHistMC;
  //TH1D *BeamSelHistData;
 private:

  TChain *tuple;
  TFile *outputFile;
  TPostScript *ps;
  ofstream IDfile;
  int psPage;

  long numEventsToProcess;

};

#endif

 

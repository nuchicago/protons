
#ifndef PROTONXSEC_H
#define PROTONXSEC_H


#include "../Utilities/LArIATAnalysis.h"
#include "../Utilities/UtilityFunctions.h"
#include "math.h"



class ProtonXsec : public LArIATAnalysis {

 public:

  ProtonXsec( );
  ProtonXsec( char* jobOptionsFile );
 
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
  int psPage;

  long numEventsToProcess;

};

#endif

 


#ifndef PROTONXSEC_H
#define PROTONXSEC_H


#include "../Utilities/LArIATAnalysis.h"


class ProtonAnalyzerMC : public LArIATAnalysis {

 public:

  ProtonAnalyzerMC( );
  ProtonAnalyzerMC( char* jobOptionsFile );
 
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

 

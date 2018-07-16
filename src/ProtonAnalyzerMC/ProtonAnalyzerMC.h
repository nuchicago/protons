
#ifndef PROTONXSEC_H
#define PROTONXSEC_H


#include "../Utilities/LArIATAnalysis.h"
#include <regex> 


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
  TFile *correctFile;
  TPostScript *ps;
  int psPage;

  long numEventsToProcess;

};

#endif

 

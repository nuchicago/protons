
#ifndef PROTONXSEC_H
#define PROTONXSEC_H


#include "../Utilities/LArIATAnalysis.h"


class ProtonXsec : public LArIATAnalysis {

 public:

  ProtonXsec( );
  ProtonXsec( char* jobOptionsFile );
 
  //------------------------------------
  // function signatures
  //------------------------------------

  void AnalyzeFromNtuples( );


 private:

  TChain *tuple;
  TFile *outputFile;
  TPostScript *ps;
  int psPage;

  long numEventsToProcess;

};

#endif

 

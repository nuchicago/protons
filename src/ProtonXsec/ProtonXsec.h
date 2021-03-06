
#ifndef PROTONXSEC_H
#define PROTONXSEC_H


#include "../Utilities/LArIATAnalysis.h"
#include "../Utilities/UtilityFunctions.h"
#include "math.h"
#include <stdio.h>
#include <string.h>



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
  ofstream logFile;
  ofstream IDfile;
  ofstream multiMatchFile;
  TFile *beamPlotFile;
  TFile *haloPlotFile;
  TFile *beamPionFile;
  TFile *beamProtonFile;
  TFile *beamKaonFile;

  //ofstream XSecLog;
  int psPage;

  long numEventsToProcess;

};

#endif

 

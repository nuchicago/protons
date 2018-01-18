
#ifndef PROTONXSEC_H
#define PROTONXSEC_H


#include "../Utilities/LArIATAnalysis.h"


class ProtonXsec : public LArIATAnalysis {

 public:

  ProtonXsec( );
 

  //------------------------------------
  // function signatures
  //------------------------------------

  void AnalyzeFromNtuples( );


 private:

  TChain *tuple;

};

#endif

 

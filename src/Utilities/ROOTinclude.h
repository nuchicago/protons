//#############################################################################
//
// ROOTinclude -                    
//
// D. Schmitz  Feb 8, 2005
//
//#############################################################################


//-----------------------------------------------------------
// Included files needed from ROOT and C++
//-----------------------------------------------------------

#include<iomanip>
#include<iostream>
#include<stdlib.h>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<string>
#include<sstream>
#include<time.h>
#include<zlib.h>
#include<algorithm>
#include<vector>

using std::cout;
using std::cin;
using std::endl;
using std::streamsize;
using std::setprecision;
using std::string;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::sort;
using std::vector;

// ROOT headers
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TMath.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TPostScript.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TText.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TArrow.h"
#include "TLine.h"
#include "TMarker.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMarker.h"
#include "TMultiGraph.h"
#include "TPolyMarker.h"
#include "TPaveText.h"
#include "TPolyMarker3D.h"
#include "TVector.h"
#include "TArray.h"
#include "TLegend.h"
#include "TImage.h"
#endif


 

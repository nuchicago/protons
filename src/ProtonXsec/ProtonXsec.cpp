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



//=============================================================================
// MAIN
//=============================================================================

int main( int argc, char * argv[] ){

  char *jobOptionsFile;

  if (argc == 1 )
  {
    jobOptionsFile    = new char[200];

    cout << endl;
    cout << "Welcome to your favorite LArIAT analysis macro...it's mine too." << endl;
    cout << "First, a question to get to know you and your analysis needs a little better." << endl;
    
    cout << endl;
    cout << "Where is your jobOptions file (full path) : ";
    cin >> jobOptionsFile;
  }
  else if( argc == 2 )
  {
    jobOptionsFile   = argv[1];
  }
  else
  {
    cout << "Program confused by user's input...exiting" << endl;
    exit(0);
  }

  //--------------------------------------------------------------------------
  // Run ROOT in batch mode
  //--------------------------------------------------------------------------

  TROOT theRoot("theRoot","root for HARP Analysis");
  gROOT->SetBatch(true);


  //--------------------------------------------------------------------------
  // ROOT Style Settings
  //--------------------------------------------------------------------------
  
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(0000);
  gStyle->SetOptTitle(0);
  gStyle->SetFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPalette(18,0);
  gStyle->SetPaperSize(20,26);

  //gRandom->SetSeed(0);



  //--------------------------------------------------------------------------
  // Create a ProtonXsec object and call main analysis function
  //--------------------------------------------------------------------------

  ProtonXsec *protonXsec = new ProtonXsec( jobOptionsFile );

  protonXsec->AnalyzeFromNtuples();
  
  return 0;

  
} // end main() 


//=============================================================================
// Constructors
//=============================================================================

ProtonXsec::ProtonXsec( ) : LArIATAnalysis( ) { 

  tuple = NULL;



}


ProtonXsec::ProtonXsec( char* jobOptionsFile ) : LArIATAnalysis( jobOptionsFile ) { 

  tuple = NULL;
  

  //== Open list of ntuple files 
  if( UI->inputFilesSet ){
    openNtupleFiles( UI->inputFiles, tuple );
  }
  else{
    cout << endl << "#### ERROR : Input files not properly specified!!!!" << endl << endl;
    exit(0);
  }

  //== Open output root file and postscript file
  if( UI->rootOutputFileSet && UI->psOutputFileSet ){
    outputFile = new TFile( UI->rootOutputFile, "RECREATE" );
    ps = new TPostScript( UI->psOutputFile, 112 );
    ps->Range(26,18); 
    psPage = 1; 
  }
  else{
    cout << endl << "#### No output files specified!!!!" << endl << endl;
  }

  //== number of events to process
  if( UI->numEventsToProcessSet )
    numEventsToProcess = UI->numEventsToProcess;
  else
    numEventsToProcess = 100000000;
  
  
  verbose = UI->verbose;
  isMC = UI->isMC;
  
}


//=============================================================================
// AnalyzeFromNtuples
//=============================================================================

void ProtonXsec::AnalyzeFromNtuples(){


  std::cout << "I calculate the p-Ar cross section! \n";
  
  EventSelector *ES = new EventSelector();
  std::cout << ES->classifyEvent( 4 ) << std::endl;
  std::cout << ES->classifyEvent( 7 ) << std::endl;


  BeamSelector *BS = new BeamSelector();
  
  bookNtuple( tuple );
  if (tuple == 0) return;

  Long64_t nentries = tuple->GetEntriesFast();
  TCanvas *c = new TCanvas;
  TH1D *BeamSelHistMC = new TH1D("BeamSelHistMC","True False Proton Histogram: MC",100,-3,3);
  BeamSelHistMC->GetYaxis()->SetTitle("Number of Events");
  TH1D *BeamSelHistData = new TH1D("BeamSelHistData","True False Proton Histogram: Data",100,-3,3);
  BeamSelHistData->GetYaxis()->SetTitle("Number of Events");
   
  // ## event loop ##
  for (Long64_t jentry=0; jentry < numEventsToProcess && jentry < nentries; jentry++) {
    std::cout << "loading Entry" << jentry << std::endl;
    Long64_t ientry = tuple->LoadTree(jentry);
    if (ientry < 0){continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);
    
    printEvent();
    bool beam_result = BS->isProton( track_zpos, ntracks_reco, isMC, UI->zBeamCutoff);
    if(beam_result){
      std::cout << "is Proton\n" << std::endl;
      if(isMC){
        BeamSelHistMC->Fill(1);
      
      }
      else{BeamSelHistData->Fill(1);}

    }
    else{
      std::cout << "not Proton\n" << std::endl;
      if(isMC){BeamSelHistMC->Fill(-1);}
      else{BeamSelHistData->Fill(-1);}
    }
  }


  
  
  
  if(isMC){
    c->cd();

    BeamSelHistMC->Draw();


    if(UI->rootOutputFileSet){
       BeamSelHistMC->Write();
    }

    TImage *img = TImage::Create();
    img->FromPad(c);
    img->WriteImage("BeamSelHistMC.png");}
  else if (!(isMC)){BeamSelHistData->Draw();
    c->cd();


    if(UI->rootOutputFileSet){
       BeamSelHistData->Write();
    }

    TImage *img = TImage::Create();
    img->FromPad(c);
    img->WriteImage("BeamSelHistData.png");}




}




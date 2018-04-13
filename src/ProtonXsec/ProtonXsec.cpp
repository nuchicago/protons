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
  //gStyle->SetOptTitle(0);
  //gStyle->SetFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  //gStyle->SetStatColor(0);
  //gStyle->SetTitleColor(0);
  //gStyle->SetPadBorderMode(0);
  //gStyle->SetFrameBorderMode(0);
  //gStyle->SetCanvasBorderMode(0);
  //gStyle->SetPalette(18,0);
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

  // counters for cuts

  double numEventsStart = 0;
  double numWCTrack = 0;
  double xyDeltaCut = 0;
  double angleCut = 0;




  std::cout << "I calculate the p-Ar cross section! \n";
  
  EventSelector *ES = new EventSelector();

  BeamSelector *BS = new BeamSelector();
  
  bookNtuple( tuple , UI->isMC);
  if (tuple == 0) return;

  Long64_t nentries = tuple->GetEntriesFast();
  TCanvas *c = new TCanvas("c","c",1000, 1000);


  TH1D *BeamSelHistMC = new TH1D("BeamSelHistMC","Primary Track Identified for Event: MC",100,-3,3);
  BeamSelHistMC->GetYaxis()->SetTitle("Number of Events");
  
  TH1D *BeamSelHistData = new TH1D("BeamSelHistData","Primary Tracks Identified",100,-3,3);
  BeamSelHistData->GetYaxis()->SetTitle("Number of Events");
  BeamSelHistData->GetYaxis()->SetTitle("Number of Events");

  TH1D *BeamToF = new TH1D("BeamToF","Incoming Particle TOF",20,-2,2);
  BeamToF->GetXaxis()->SetTitle("ns");
  BeamToF->GetYaxis()->SetTitle("Number of Events");

  

  TH2D *TrackPositionXY =  new TH2D("TrackPositionXY","Position of TPC track start",
    20, 0, 47.5, 20, -20, 20);
  TrackPositionXY->GetYaxis()->SetTitle("Y");
  TrackPositionXY->GetXaxis()->SetTitle("X");

   

  TH1D *BeamMomentum = new TH1D("BeamMomentum","Incoming Particle Momentum",100,300,1000);
  BeamMomentum->GetXaxis()->SetTitle("[MeV/c]");
  BeamMomentum->GetYaxis()->SetTitle("Number of Events");

  TH2D * delXYHist =  new TH2D("delXYHist","tpc to wc delta x",100,-40,40,100,-40,40);
  delXYHist->GetXaxis()->SetTitle("X[cm]");
  delXYHist->GetYaxis()->SetTitle("Y[cm]");
  TH1D * delThetaHist =  new TH1D("delThetaHist","tpc to wc delta Theta",100,-1,1);
  delThetaHist->GetXaxis()->SetTitle("[radians]");
  delThetaHist->GetYaxis()->SetTitle("Number of Events");
  TH1D * delPhiHist =  new TH1D("delPhiHist","tpc to wc delta Phi",100,-4,4);
  delPhiHist->GetXaxis()->SetTitle("[radians]");
  delPhiHist->GetYaxis()->SetTitle("Number of Events");


  // ## event loop ##
  for (Long64_t jentry=0; jentry < numEventsToProcess && jentry < nentries; jentry++) {
    
    Long64_t ientry = tuple->LoadTree(jentry);
    if (ientry < 0){continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);

    int reco_primary = -1;
    
    numEventsStart++;
    if(verbose == 2){printEvent();}
    bool found_primary = false; //BS->PrimaryTrack( track_zpos, ntracks_reco, UI->zBeamCutoff,
     //reco_primary, first_reco_z,verbose);


    
    if(!isMC){
      if (num_wctracks < 2){
        numWCTrack++;
        std::vector< std::vector<double> > wctpc_mvect = BS->wcTPCMatchPlots(wctrk_XFace[0],wctrk_YFace[0], wctrk_theta[0], wctrk_phi[0], track_xpos,
         track_ypos, track_zpos, ntracks_reco, UI->zTPCCutoff);

        if(wctpc_mvect.size() > 0 ){
          for(int rtrack = 0; rtrack < wctpc_mvect.size(); rtrack++){
            delXYHist->Fill(wctpc_mvect[rtrack][0],wctpc_mvect[rtrack][1]);
            delThetaHist->Fill(wctpc_mvect[rtrack][2]);
            delPhiHist->Fill(wctpc_mvect[rtrack][3]);
          }
          std::vector <double> MinVector = BS->wcTPCMatch(wctrk_XFace[0],wctrk_YFace[0], wctrk_theta[0], wctrk_phi[0],track_xpos,
          track_ypos, track_zpos, ntracks_reco, UI->zTPCCutoff,reco_primary);

          if(reco_primary > 0 && MinVector[0] < 5){
            xyDeltaCut++;
            if(abs(MinVector[1]) < 1 && abs(MinVector[2]) <1){
              angleCut++;
            }


          }





        }





      }
    }



    if(found_primary){
      if(verbose == 2){std::cout << "Found Primary\n" << std::endl;}
      if(isMC){
        BeamSelHistMC->Fill(1);        
        TrackPositionXY->Fill((*track_xpos)[reco_primary][0],
          (*track_ypos)[reco_primary][0]);
      
      }
      else{//
        BeamSelHistData->Fill(1);
          TrackPositionXY->Fill((*track_xpos)[reco_primary][0],
          (*track_ypos)[reco_primary][0]);

          for (int wctrack = 0 ; wctrack < num_wctracks; wctrack++){
          BeamMomentum->Fill(wctrk_momentum[wctrack]);
          BeamToF->Fill(tofObject[wctrack]);}


      }

    }
    else{
      if(verbose == 2){std::cout << "No Valid Primary \n" << std::endl;}
      if(isMC){BeamSelHistMC->Fill(-1);}
      else{BeamSelHistData->Fill(-1);}

    }
  } //end of event Loop

  if(verbose > 0){
    std::cout << "Total events processed: "<< numEventsStart << std::endl;
    std::cout << "Events with only one wc track: "<< numWCTrack << std::endl;
    std::cout << "Events after XY plane distance cut: "<< xyDeltaCut << std::endl;
    std::cout << "Events after angle cut: "<< angleCut << std::endl;
  }


  
  
  
  if(isMC){
    c->cd();

    BeamSelHistMC->Draw();


    if(UI->rootOutputFileSet){
       BeamSelHistMC->Write();
    }
    /*if(UI->psOutputFileSet){
       BeamSelHistMC->Print(UI->psOutputFile);
       TrackPositionXY->Print(UI->psOutputFile);
    }*/
    TImage *img = TImage::Create();
    img->FromPad(c);
    img->WriteImage("BeamSelHistMC.png");
    TrackPositionXY->Draw("COLZ");
    TImage *BeamXYimg = TImage::Create();
    BeamXYimg->FromPad(c);
    BeamXYimg->WriteImage("BeamXY.png");}

  



  else if (!(isMC)){
    c->cd();
    
    


    if(UI->rootOutputFileSet){



    BeamSelHistData->Draw();
    BeamSelHistData->Write();
    TImage *img = TImage::Create();
    img->FromPad(c);
    img->WriteImage("BeamSelHistData.png");
    
    TrackPositionXY->Draw("COLZ");
    TrackPositionXY->Write();
    TImage *BeamXYimg = TImage::Create();
    BeamXYimg->FromPad(c);
    BeamXYimg->WriteImage("BeamXYData.png");
    
    BeamToF->Draw();
    BeamToF->Write();
    TImage *Beamtofimg = TImage::Create();
    Beamtofimg->FromPad(c);
    Beamtofimg->WriteImage("BeamToFData.png");

    BeamMomentum->Draw();
    BeamMomentum->Write();
    TImage *BeamMomImg = TImage::Create();
    BeamMomImg->FromPad(c);
    BeamMomImg->WriteImage("BeamMomentum.png");


    delXYHist->Draw("COLZ");
    delXYHist->Write();
    //XYcut->Write();
    TImage *delXYimg = TImage::Create();
    delXYimg->FromPad(c);
    delXYimg->WriteImage("delXYHist.png");




    delThetaHist->Draw();
    delThetaHist->Write();
    TImage *delThetaimg = TImage::Create();
    delThetaimg->FromPad(c);
    delThetaimg->WriteImage("delThetaHist.png");


    delPhiHist->Draw();
    delPhiHist->Write();
    TImage *delPhiimg = TImage::Create();
    delPhiimg->FromPad(c);
    delPhiimg->WriteImage("delPhiHist.png");
    }


}




}




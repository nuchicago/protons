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
  
  //gStyle->SetOptStat(0000);
  //gStyle->SetOptFit(0000);
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
  //gStyle->SetPaperSize(20,26);

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
  if( UI->rootOutputFileSet ){
    outputFile = new TFile( UI->rootOutputFile, "RECREATE" );
    //ps = new TPostScript( UI->psOutputFile, 112 );
    //ps->Range(26,18); 
    //psPage = 1; 
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

  // counters for beam cuts

  double numEventsStart = 0;
  double numZcutoff = 0;
  double numWCTrack = 0;
  double xyDeltaCut = 0;
  double angleCut = 0;




  // ### some variables that are needed for the xs calc ###
  double z2 = UI->zSlabSize; //<-- slab size. need to move this to jobOptions



  std::cout << "I calculate the p-Ar cross section! \n";
  
  EventSelector *ES = new EventSelector();

  BeamSelector *BS = new BeamSelector();
  
  bookNtuple( tuple , UI->isMC);
  if (tuple == 0) return;

  Long64_t nentries = tuple->GetEntriesFast();
  TCanvas *c = new TCanvas("c","c",1000, 1000);





  

  TH2D *tpcInTracksXY =  new TH2D("tpcInTracksXY","Position of TPC track start XY",
    200, -100, 100, 200, -100, 100);



  TH2D *wctrkPositionXY =  new TH2D("wctrkPositionXY","Position of Wire Chamber Track",
    200,-100 ,100 , 200, -100, 100);
  TH2D *wctrkSelectedXY =  new TH2D("wctrkSelectedXY","Position of Wire Chamber Track",
  200, -100, 100, 200, -100, 100);


   
  TH1D *tpcInTracksZ =  new TH1D("tpcInTracksZ","Position of TPC track start Z",10, 0, 4);
  TH1D *PrimaryStartZ =  new TH1D("PrimaryStartZ","Selected track start in Z",10, 0, 4);
  TH1D *BadTrackStartZ =  new TH1D("BadTrackStartZ","non-selected track start in Z",10, 0, 4);
  TH1D *InTrackLength =  new TH1D("InTrackLength","Entering Track Length",250, 0, 100);
  TH1D *PrimaryLength =  new TH1D("PrimaryLength","Selected Entering Track Length",250, 0, 100);
  TH1D *BadTrackLength =  new TH1D("BadTrackLength","non-selected Entering Track Length",250, 0, 100);

  TH1D *BeamMomentum = new TH1D("BeamMomentum","Incoming Particle Momentum",100,300,1000);
  TH1D *BeamToF = new TH1D("BeamToF","Incoming Particle Momentum",100,0,200);
  TH1D *inTracksNumHist = new TH1D("inTracksNumHist","Number of Entering Tracks TPC",100,0,10);

  TH1D *wctrkNumHist = new TH1D("wctrkNumHist","Number of Entering Tracks TPC",20,0,5);

  TH2D * delXYHist =  new TH2D("delXYHist","tpc to wc delta x",200,-100,100,200,-100,100);
  TH2D *BadTrackHist =  new TH2D("BadTrackHist","non-selected tracks",200,-100,100,200,-100,100);
  TH2D *delBadTrackHist =  new TH2D("delBadTrackHist","non-selected tracks",200,-100,100,200,-100,100);
  TH1D * delThetaHist =  new TH1D("delThetaHist","tpc to wc delta Theta",100,-1,1);
  TH1D * delPhiHist =  new TH1D("delPhiHist","tpc to wc delta Phi",100,-4,4);

  // ## xs histos ##
  TH1D *hreco_initialKE = new TH1D("hreco_initialKE", "initial ke", 20, 0, 1000);
  TH1D *hreco_intke = new TH1D("hreco_intke", "int ke", 20, 0, 1000);
  TH1D *hreco_incke = new TH1D("hreco_incke", "inc ke", 20, 0, 1000);

  // ## looping once to find Beamline center ##

  if(!isMC){
  for (Long64_t jentry=0; jentry < numEventsToProcess && jentry < nentries; jentry++) {


    //BS->MassCut(wctrk_momentum,tofObject,)
    
    Long64_t ientry = tuple->LoadTree(jentry);
    if (ientry < 0){continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);
    int numEnteringTracks = 0;
    int best_candidate = -1;
    wctrkNumHist->Fill(num_wctracks);

      if (num_wctracks == 1){
        wctrkPositionXY->Fill(wctrk_XFace[0],wctrk_YFace[0]);
        BeamMomentum->Fill(wctrk_momentum[0]);
        BeamToF->Fill(tofObject[0]);

        std::vector <double> wctpc_mvect = BS->wcTPCMatch(wctrk_XFace[0],wctrk_YFace[0], wctrk_theta[0],
          wctrk_phi[0],track_xpos,track_ypos, track_zpos, ntracks_reco, UI->zTPCCutoff,
          best_candidate, numEnteringTracks);

        inTracksNumHist->Fill(numEnteringTracks);

        //phicalctest->Fill(wctrk_phi[0] - atan2(wctrk_XFace[0]-20,wctrk_YFace[0]));

        if(best_candidate != -1){

            delXYHist->Fill(wctpc_mvect[1],wctpc_mvect[2]);
            delThetaHist->Fill(wctpc_mvect[3]);
            delPhiHist->Fill(wctpc_mvect[4]);
            tpcInTracksXY->Fill((*track_xpos)[best_candidate][0],
            (*track_ypos)[best_candidate][0]);
            wctrkSelectedXY->Fill(wctrk_XFace[0],wctrk_YFace[0]);         
          }
        }
      }
    }
  // End of beam centering loop

    // setting beam center values
  
    double xMeanTPCentry = delXYHist->GetMean(1);
    double yMeanTPCentry = delXYHist->GetMean(2);

  // ## event loop ##
  for (Long64_t jentry=0; jentry < numEventsToProcess && jentry < nentries; jentry++){
    
    Long64_t ientry = tuple->LoadTree(jentry);
    if (ientry < 0){continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);
    int numEnteringTracks = 0;
    int reco_primary = -1;
    
    numEventsStart++;
    if(verbose == 2){printEvent();}
    bool found_primary = false; //BS->PrimaryTrack( track_zpos, ntracks_reco, UI->zBeamCutoff,
     //reco_primary, first_reco_z,verbose);

    if(!isMC){
      if (num_wctracks !=1){continue;}
      else{
          numWCTrack++;
          std::vector <double> MinVector = BS->wcTPCMatch(wctrk_XFace[0],wctrk_YFace[0], wctrk_theta[0], wctrk_phi[0],track_xpos,
          track_ypos, track_zpos, ntracks_reco, UI->zTPCCutoff,reco_primary, numEnteringTracks);
          
          if(reco_primary == -1){
            if(verbose == 2){std::cout << "No Valid Primary \n" << std::endl;}
          }
          else{ 
            numZcutoff++;
              if((sqrt(pow((MinVector[1] - xMeanTPCentry),2) + pow((MinVector[2] - yMeanTPCentry),2)))> UI-> rCircleCut) {
                  continue;}
              else{
                  xyDeltaCut++;
                  
                for (int inTrack = 0; inTrack < ntracks_reco; inTrack++ ){
                  if((*track_zpos)[inTrack][0] < UI->zTPCCutoff){
                    tpcInTracksZ->Fill((*track_zpos)[inTrack][0]);
                    InTrackLength->Fill((*track_length)[inTrack]);
                    
                    if (inTrack != reco_primary){
                      BadTrackHist->Fill((*track_xpos)[inTrack][0],(*track_ypos)[inTrack][0]);
                      delBadTrackHist->Fill(wctrk_XFace[0] - (*track_xpos)[inTrack][0],
                          wctrk_YFace[0] - (*track_ypos)[inTrack][0]);
                      BadTrackLength->Fill((*track_length)[inTrack]);
                      BadTrackStartZ->Fill((*track_zpos)[inTrack][0]);

                    }
                    if(inTrack == reco_primary){
                      PrimaryLength->Fill((*track_length)[inTrack]);
                      PrimaryStartZ->Fill((*track_zpos)[inTrack][0]);
                    }
                  }
                }
              }
            }
          }
      }

    else{
      if(isMC){//BeamSelHistMC->Fill(-1);
      }
      else{//BeamSelHistData->Fill(-1);
      }
    }

    // ### porting over the work from ProtonAnalyzerMC module -- ryan ###


    // ## grabbing reco primary ##
    //int reco_primary = -1;
    if(isMC){
      double first_reco_z = 99.;
      reco_primary = BS->isTPCPrimary(track_zpos, ntracks_reco, isMC, UI->zBeamCutoff, 
        reco_primary, first_reco_z, verbose);}
    if(reco_primary == -1){continue;}//<- skipping events that didn't pass isTPCPrimary 
    else{ found_primary = true;}



    // ## grabbing interaction point ##
    double temp[4];
    double* candidate_info = ES->findInt(temp, reco_primary, ntracks_reco, 
                                            ntrack_hits, track_xpos, track_ypos, track_zpos,
                                            track_end_x, track_end_y, track_end_z,
                                            col_track_hits, col_track_dedx, col_track_pitch_hit,
                                            col_track_x, col_track_y, col_track_z);


    // ## grabbing what will be histogram entries ##
    double initial_ke = 99999; //<-- setting to a constant to do dev. needs to be WC info
    std::vector<double> calo_slab_xpos;
    std::vector<double> calo_slab_ypos;
    std::vector<double> calo_slab_zpos;
    std::vector<double> calo_slab_KE;

    if(isMC) {
      initial_ke = BS->getMCInitialKE(initial_ke, geant_list_size, process_primary, 
                                      NTrTrajPts, MidPosX, MidPosY,  MidPosZ, MidPx, MidPy, MidPz); 
    }
    if(!isMC) {
      initial_ke = BS->getDataInitialKE(initial_ke, wctrk_momentum[0]);
    }
    hreco_initialKE->Fill(initial_ke);
    
    int slabPass = ES->getSlabInfo(calo_slab_xpos, calo_slab_ypos, calo_slab_zpos, calo_slab_KE,
                                    reco_primary, z2, initial_ke,
                                    col_track_hits, col_track_dedx, col_track_pitch_hit,
                                    col_track_x, col_track_y, col_track_z);


    
    // ## actually filling histograms ##
    // ## getting which slab will be used for interactions ##
    int calo_int_slab = 999;
    if(candidate_info[0]){
      double int_candidate_x = candidate_info[1];
      double int_candidate_y = candidate_info[2];
      double int_candidate_z = candidate_info[3];
      // loop over slabs to find slab closest to interaction candidate
      double min_dist_int = 99;
      for(unsigned int calo_slab = 0; calo_slab < calo_slab_KE.size(); calo_slab++){
        double calo_slab_x = calo_slab_xpos[calo_slab];  
        double calo_slab_y = calo_slab_ypos[calo_slab];  
        double calo_slab_z = calo_slab_zpos[calo_slab];  
        double dist_int_slab = sqrt( pow(int_candidate_x - calo_slab_x, 2) 
                                   + pow(int_candidate_y - calo_slab_y, 2) 
                                   + pow(int_candidate_z - calo_slab_z, 2));
        if(calo_slab_z < int_candidate_z){
          if(dist_int_slab < min_dist_int){
            min_dist_int = dist_int_slab;
            calo_int_slab = calo_slab;
          }
        }//<--End if this slab is upstream of int
      }//<--End calo slab loop
    }//<---End if interaction candidate


    // ## incident slabs ## 
    for(unsigned int calo_slab = 1; calo_slab < calo_slab_KE.size(); calo_slab++){
      if(calo_slab > calo_int_slab){continue;}//<--stop after interaction slab 
      hreco_incke->Fill(calo_slab_KE[calo_slab]);
      if(calo_slab == calo_int_slab){
        hreco_intke->Fill(calo_slab_KE[calo_slab]);
      }//<-- End if this is the interacting slab
    }//<--End calo slab loop

    // ### end of port work -- ryan ###



   
  
  }//end of event Loop

  if(verbose > 0){
    std::cout << "Total events processed: "<< numEventsStart << std::endl;
    std::cout << "Events with only one wc track: "<< numWCTrack << std::endl;
    std::cout << "Events with best match TPC track Z < Cutoff: "<< numZcutoff << std::endl;
    std::cout << "Events after XY plane distance cut: "<< xyDeltaCut << std::endl;
    //std::cout << "Events after; angle cut: "<< angleCut << std::endl;
  }


  
  
  
  if(isMC){

    if(UI->rootOutputFileSet){
       
      // ## xs histos ##
      hreco_initialKE->Write();
      hreco_incke->Write();
      hreco_intke->Write();
  }


  }

  



  else if (!(isMC)){
    
    if(UI->rootOutputFileSet){
    
    //BeamSelHistData->Write();

    
    PrimaryLength->Write();
    PrimaryStartZ->Write();
    tpcInTracksXY->Write();
    tpcInTracksZ->Write();
    InTrackLength->Write();
    wctrkPositionXY->Write();
    wctrkSelectedXY->Write();
    wctrkNumHist->Write();
    inTracksNumHist->Write();
    BeamMomentum->Write();
    BeamToF->Write();
    delXYHist->Write();
    delThetaHist->Write();
    delPhiHist->Write();
    BadTrackHist->Write();
    BadTrackLength->Write();
    BadTrackStartZ->Write();
    delBadTrackHist->Write();

    // ## xs histos ##
    hreco_initialKE->Write();
    hreco_incke->Write();
    hreco_intke->Write();
    }
  }
}





//#############################################################################
//
// BendStudy.cpp
//
//  Measures the track bending issue in LArIAT data
//
// January 17, 2018
//
//#############################################################################

#include "BendStudy.h"

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

  TROOT theRoot("theRoot","root for XS Analysis");
  gROOT->SetBatch(true);


  //--------------------------------------------------------------------------
  // ROOT Style Settings
  //--------------------------------------------------------------------------
  
  //gStyle->SetOptStat(0000);
  //gStyle->SetOptFit(0000);
  //gStyle->SetOptTitle(0);
  //gStyle->SetFillColor(0);
  //gStyle->SetPadColor(0);
  //gStyle->SetCanvasColor(0);
  //gStyle->SetStatColor(0);
  //gStyle->SetTitleColor(0);
  //gStyle->SetPadBorderMode(0);
  //gStyle->SetFrameBorderMode(0);
  //gStyle->SetCanvasBorderMode(0);
  //gStyle->SetPalette(18,0);
  //gStyle->SetPaperSize(20,26);
  gStyle->SetOptStat("emrou");

  //gRandom->SetSeed(0);



  //--------------------------------------------------------------------------
  // Create a BendStudy object and call main analysis function
  //--------------------------------------------------------------------------

  BendStudy *analyzeBending = new BendStudy( jobOptionsFile );

  analyzeBending->AnalyzeFromNtuples();
  
  return 0;

  
} // end main() 


//=============================================================================
// Constructors
//=============================================================================

BendStudy::BendStudy( ) : LArIATAnalysis( ) { 

  tuple = NULL;



}


BendStudy::BendStudy( char* jobOptionsFile ) : LArIATAnalysis( jobOptionsFile ) { 

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

  //output list of selected primary particle IDs

  if (UI->SelEventListSet){
    
    //IDfile.open(UI->SelEventList, ios::trunc);
  }

  //== number of events to process
  if( UI->numEventsToProcessSet )
    numEventsToProcess = UI->numEventsToProcess;
  else
    numEventsToProcess = 100000000;


  
  
  verbose = UI->verbose;
  isMC = UI->isMC;
  applyMassCut =  UI->applyMassCut;
  
}


//=============================================================================
// AnalyzeFromNtuples
//=============================================================================

void BendStudy::AnalyzeFromNtuples(){

  // counters for beam cuts



  double numEventsStart = 0;
  double numMassCut = 0;
  double numtofvalid = 0;
  double numZcutoff = 0;
  double numHasCandidate = 0; 
  double numTwoTracks = 0;
  double numzProjShort = 0;
  double numBendZcut = 0;
  double numRightSided = 0;
  double numHasLong = 0 ;
  double numXYZMatching = 0;
  double numZmatching = 0;

  double numWCTrack = 0;
  double xyDeltaCut = 0;
  double numThetaCut = 0;
  double numPhiCut = 0;
  int multiple_matches = 0;
  double numInteractions = 0; 
  double NumEventsSelList = 0;
  bool intHistFilled = false;


  const double wc_zpos_val = 0.1;
  const double *wc_zpos;
  wc_zpos = &wc_zpos_val;

  double bendZcut  =  UI->zBeamCutoff; //using this leftover input option to select stubby tracks

  // For individual event plotting mode (plotIndividual)

  //std::vector<double> total_dedx_v;
  //std::vector<double> total_res_v;
  //TGraph gtotal_res_dedx_item;
  //TGraph * gtotal_res_dedx;



  // ### some variables that are needed for the xs calc ###
  double z2 = UI->zSlabSize; //<-- slab size. need to move this to jobOptions

  // EventSelector Options in a vector
  std::vector<double> ESoptions = {static_cast<double> (verbose), UI->dedxNoBraggMax,
    UI->branchMaxDist, UI->clusterMaxDist};

  
  // globals move later...
  double mass = 938.57;
  //double z = 0.03;
  //double z2 = 0.5;

  double rho = 1.3954;                   // ## g/cm3
  double molar_mass = 39.95;                    // ## g/mol
  double N_A = 6.022 * pow(10, 23);      // ## num/mol
  
  double sparse_recip_num_density = molar_mass / (rho * z2 * N_A); // ## cm2/num
  double barn = pow(10, -24);



  std::cout << "I calculate the p-Ar cross section! \n";
  
  EventSelector *ES = new EventSelector();

  BeamSelector *BS = new BeamSelector();
  
  bookNtuple( tuple , UI->isMC);
  if (tuple == 0) return;

  Long64_t nentries = tuple->GetEntriesFast();
  TCanvas *c = new TCanvas("c","c",1500, 750);
  //c->Divide(1,2);

  //TCanvas *csplit = new TCanvas("csplit","csplit",1500, 750);
  //csplit->Divide(2,2);

// ## Histograms for entering tracks
  //TH2D *tpcInTracksXY =  new TH2D("tpcInTracksXY","Position of TPC track start XY",
    //200, -100, 100, 200, -100, 100);
  TH2D *wctrkPositionXY =  new TH2D("wctrkPositionXY","Position of Wire Chamber Track",
    200,-100 ,100 , 200, -100, 100);
  //TH2D *wctrkSelectedXY =  new TH2D("wctrkSelectedXY","Position of Wire Chamber Track",
  //200, -100, 100, 200, -100, 100);
  TH1D *tpcInTracksZ =  new TH1D("tpcInTracksZ","Position of TPC track start Z",25, 0, 10);
  TH1D *ShortStartZ =  new TH1D("ShortStartZ","Earlies track start in Z",25, 0, 10);
  TH1D *ShortEndZ =  new TH1D("ShortEndZ","Earlies track end in Z",25, 0, 10);




  //TH1D *BadTrackStartZ =  new TH1D("BadTrackStartZ","non-selected track start in Z",10, 0, 4);
  //TH1D *InTrackLength =  new TH1D("InTrackLength","Entering Track Length",100, 0, 100);
  TH1D *ShortLength =  new TH1D("ShortLength","Earliest Entering Track Length",200, 0, 100);
  TH1D *LongLength =  new TH1D("LongLength","Second Entering Track Length",200, 0, 100);
  //TH1D *BadTrackLength =  new TH1D("BadTrackLength","non-selected Entering Track Length",250, 0, 100);

  TH1D *inTracksNumHist = new TH1D("inTracksNumHist","Number of Entering Tracks TPC",100,0,10);
  TH1D *CylNumHist = new TH1D("CylNumHist","Number of Tracks in selection cylinder",100,0,10);


  TH1D *BeamMomentum = new TH1D("BeamMomentum","Incoming Particle Momentum",200,0,2000);
  TH1D *wctrkNumHist = new TH1D("wctrkNumHist","Number of Entering Tracks TPC",20,0,5);
  TH1D *BeamToF = new TH1D("BeamToF","Incoming Particle Momentum",100,0,100);

  TH2D * delXYHist =  new TH2D("delXYHist","tpc to wc delta x",200,-100,100,200,-100,100);

  TH2D * BendingZHist =  new TH2D("BendingZHist"," Bending vs Z_{tpc}",50,0,10,50,0,2);
  TH2D * BendingXHist =  new TH2D("BendingXHist","Bending vs X_{tpc}",50,0,48,50,0,2);

  //TH1D * delXYHistPx =  new TH1D("delXYHistPx","tpc to wc delta x",200,-100,100,200,-100,100);
  //TH1D * delXYHistPy =  new TH1D("delXYHistPy","tpc to wc delta x",200,-100,100,200,-100,100);

  //TH2D *BadTrackHist =  new TH2D("BadTrackHist","non-selected tracks",200,-100,100,200,-100,100);
  //TH2D *delBadTrackHist =  new TH2D("delBadTrackHist","non-selected tracks",200,-100,100,200,-100,100);



  
  //TH1D * zProjTrack_out = new TH1D("zProjTrack_out","",100,0,40);
  //TH1D * zProjTrack_in = new TH1D("zProjTrack_in","",100,0,40);
 // TH1D * InTrackLength_out = new TH1D("InTrackLength_out","",100,0,40);
 // TH1D * InTrackLength_in = new TH1D("InTrackLength_in","",100,0,40);
   //z projections of tracks, in/out of circle cut
  TH1D * zProjTrack_short = new TH1D("zProjTrack_short","Earliest Entering Track Length - Z projection",100,0,40);
  //TH1D * InTrackTPCnum_in = new TH1D("InTrackTPCnum_in","",100, 0, 10);
  //TH1D * InTrackTPCnum_out = new TH1D("InTrackTPCnum_out","",100, 0, 10);

  TH1D *BendingProxHist = new TH1D("BendingProxHist","short to long tagged track distance", 50, 0, 10);



  TH1D * numTracksSelHist =  new TH1D("numTracksSelHist","number of Entering Tracks - Selected Events", 10, 0, 5);

  TH2D *tofMomentHist = new TH2D("tofMomentHist","Momentum vs TOF",100,0,2000, 100 , 0,100);
  TH1D *BeamMassHist = new TH1D("BeamMassHist","Beamline particle Mass", 100, 0,3000);
  TH1D *BeamMassCutHist = new TH1D("BeamMassCutHist","Beamline particle Mass - after Cut", 100, 0,3000);
  
  //plots for EventSelection Branches

  //TH1D * BranchDistHist = new TH1D("BranchDistHist", "Inelastic Event Branch Distance",50, 0,25);
  //TH1D * ClusterDistHist = new TH1D("ClusterDistHist", "Additional Branch Distance (Type 4)",50, 0,25);
  // ## plot markers for cuts

  //TLine  *MassMinLine =  new TLine(0,UI->MassCutMin, 100, UI->MassCutMin);
  //TLine  *MassMaxLine =  new TLine(0,UI->MassCutMax, 100, UI->MassCutMax);

  // ## xs histos ##
  //TH1D *hreco_initialKE = new TH1D("hreco_initialKE", "initial ke", 20, 0, 1000);
  //TH1D *hreco_intke = new TH1D("hreco_intke", "int ke", 20, 0, 1000);
  //TH1D *hreco_incke = new TH1D("hreco_incke", "inc ke", 20, 0, 1000);
  //TH1D *hreco_xs = new TH1D("hreco_xs", "P-Ar Inelastic XS", 20, 0, 1000);

  // ## to read in from corrections file
  //TH1D *hreco_intke_background = new TH1D("hreco_intke_background", "int ke background", 20, 0, 1000);
  //TH1D *hreco_incke_background = new TH1D("hreco_incke_background", "inc ke background", 20, 0, 1000);
  //TH2D *hreco_unfolding_matrix_normalized = new TH2D("hreco_unfolding_matrix_normalized", "energy unfolding matrix", 20, 0, 1000, 20, 0, 1000);
  //TH1D *hreco_intke_eff = new TH1D("hreco_intke_eff", "interaction selection efficiency", 20, 0, 1000);
  //TH1D *hreco_incke_eff = new TH1D("hreco_incke_eff", "incident selection efficiency", 20, 0, 1000);

  // ## looping once to find Beamline center ##

std::cout << "num entries in file  = " << nentries << std::endl;
std::cout << "num events to process  = " << numEventsToProcess << std::endl;





  if(!isMC){
  for (Long64_t jentry=0;  jentry < numEventsToProcess && jentry < nentries; jentry++) {
    
    //if(verbose){std::cout << "entry = " << jentry << std::endl;}
    Long64_t ientry = tuple->LoadTree(jentry);
    if (ientry < 0){continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);
     //if(verbose){std::cout << jentry << " loaded succesfully "<< std::endl;}
    
    int numEnteringTracks = 0;
    int best_candidate = -1;
    
    wctrkNumHist->Fill(num_wctracks);

    double xMeanTPCentry = 0;
    double yMeanTPCentry = 0;

    bool skipCircleCut = true;


      if (num_wctracks != 1){continue;}
      else{
          double ParticleMass = -9999999.;

          bool isProton = BS->MassCut(wctrk_momentum[0], tofObject[0], ParticleMass, UI->MassCutMin, UI->MassCutMax);

          if(applyMassCut){
            if(!isProton){continue;}
          }

        wctrkPositionXY->Fill(wctrk_XFace[0],wctrk_YFace[0]);
        //if(verbose){std::cout << "wctrk X,Y= " << wctrk_XFace[0] << " , " << wctrk_YFace[0]  << std::endl;}
        BeamMomentum->Fill(wctrk_momentum[0]);
        BeamToF->Fill(tofObject[0]);

        std::vector<double> matchCandidate = BS->BeamCentering(wctrk_XFace[0], wctrk_YFace[0],
         track_xpos, track_ypos, track_zpos, ntracks_reco, ntrack_hits, UI->zTPCCutoff, best_candidate);

        if(best_candidate != -1){

           //if(verbose){std::cout << "doing static cast int"<< std::endl;}

            //int best_candidate = static_cast <int> (matchCandidates[0][0]);
            int startIndex = static_cast <int> (matchCandidate[1]);


            //if(verbose){std::cout << "Filling delXYHist"<< std::endl;}

            delXYHist->Fill(wctrk_XFace[0] - (*track_xpos)[best_candidate][startIndex],
            wctrk_YFace[0] - (*track_ypos)[best_candidate][startIndex]);
            //tpcInTracksXY->Fill((*track_xpos)[best_candidate][startIndex],
            //(*track_ypos)[best_candidate][startIndex]);
            //wctrkSelectedXY->Fill(wctrk_XFace[0],wctrk_YFace[0]);
            //if(verbose){std::cout << "finished filling histos"<< std::endl;}
                     
          }
        } 
      }
    }
  // End of beam centering loop

    // setting beam center values
    double xMeanTPCentry = delXYHist->GetMean(1);
    double yMeanTPCentry = delXYHist->GetMean(2);

  // ## event loop ##
    if(verbose){std::cout << "Starting main loop" << std::endl;
    std::cout << "x Mean cut " << xMeanTPCentry<< std::endl;
    std::cout << "y Mean cut " << yMeanTPCentry<< std::endl;
  }
  for (Long64_t jentry=0; jentry < numEventsToProcess && jentry < nentries; jentry++){
    
    Long64_t ientry = tuple->LoadTree(jentry); 
    if (ientry < 0){
      continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);
    numEventsStart++;

    int numTracksCircle = 0;
    int numEnteringTracks = 0;
    int reco_primary = -1;
    double ParticleMass;
    
    if(verbose){printEvent();}
    bool found_primary = false;
    bool passed_geo_cuts = true;
    bool skipCircleCut = false;

    if(!isMC){
      if (num_wctracks !=1){if(verbose){std::cout << "no WC track \n" << std::endl;}
        continue;}
      else{
          numWCTrack++;
          if(tofObject[0] < 0){
            if(verbose){std::cout << "no valid ToF \n" << std::endl;}
            continue;}
          numtofvalid++;
          ParticleMass = -9999999. ;
          tofMomentHist->Fill(wctrk_momentum[0], tofObject[0]);

          bool isProton = BS->MassCut(wctrk_momentum[0], tofObject[0], ParticleMass, UI->MassCutMin, UI->MassCutMax);

          BeamMassHist->Fill(ParticleMass);
          if(!isProton){
            if(applyMassCut){continue;}
          }

          if (isProton){
            BeamMassCutHist->Fill(ParticleMass);
            numMassCut++;
          }

          //############## my crappy match finding function (I'll pack it somewhere later)


        //std::vector<std::vector<double>> zIndexes; //vector w/ vectors of indices for point z min, second min, z max.
        std::vector<std::vector<double>> matchCandidates; 
        std::vector<double> centeringMatch;

        double numZtracks = 0;
        double numCylinder = 0; 
        for (int itrack = 0; itrack < ntracks_reco ;  itrack++ ){


          std::vector<int> zIndices = UtilityFunctions::zOrderedTrack(track_zpos,itrack,ntrack_hits);

          int index1 = zIndices[0];
          int index2 = zIndices[1];
          int indexLast = zIndices[2];
          double zmin1 = (*track_zpos)[itrack][zIndices[0]];
          //std::cout << "loading zmin2 : "<< zIndices[1] << std::endl;
          double zmin2 = (*track_zpos)[itrack][zIndices[1]];
          //std::cout << "loading zlast : "<< zIndices[2] << std::endl;
          double zmax = (*track_zpos)[itrack][zIndices[2]];
          

          //zIndexes.push_back(trackInd);



        
          if(zmin1 < UI->zTPCCutoff){
            numZtracks++;

            double delX = wctrk_XFace[0] - (*track_xpos)[itrack][index1];
            double delY = wctrk_YFace[0] - (*track_ypos)[itrack][index1];
            //double delTheta = wc_theta - UtilityFunctions::getTrackTheta(itrack, track_xpos,track_ypos,track_zpos);
            //double delPhi;   ########reminder to reimplement angle.
            double alpha = -999;
            double rValue = sqrt(pow(delX,2) + pow(delY,2));



            
            double adjustedR = (sqrt(pow((delX - xMeanTPCentry),2) + pow((delY - yMeanTPCentry),2)));

              if (adjustedR < UI->rCircleCut){
                numCylinder++;
                std::vector<double> candidate = {1.* itrack, 1.* index1, 1.* index2, 1.* indexLast, rValue, alpha};
                matchCandidates.push_back(candidate);
              
            }
          }
        }

        tpcInTracksZ->Fill(numZtracks);
        CylNumHist->Fill(numCylinder);
        if (numZtracks > 0){numZcutoff++;}
        if (numCylinder > 0){numHasCandidate++;}


        
        if (matchCandidates.size() < 2){continue;} //maybe use this branch to find kinks?
        else{
            numTwoTracks++;
            double minZpos = 9999;

            std::vector<double> ShortCandidate;
            int Short_StartPoint;
            int Short_EndPoint;
            int Short_ID;


            std::vector<double> LongCandidate;
            bool found_bending = false;
            for (int icandidate = 0;  icandidate < matchCandidates.size() ; icandidate++){
              int startPoint = static_cast <int> (matchCandidates[icandidate][1]);
              int genID = static_cast <int> (matchCandidates[icandidate][0]);

              if ((*track_zpos)[genID][startPoint] < minZpos){
                ShortCandidate = matchCandidates[icandidate];
                Short_ID = genID;
                Short_EndPoint = static_cast <int> (matchCandidates[icandidate][3]);
                Short_StartPoint= startPoint;
              }
            }


          
          double zProjLength = (*track_zpos)[Short_ID][Short_EndPoint] - (*track_zpos)[Short_ID][Short_StartPoint];
          zProjTrack_short->Fill(zProjLength);
          ShortLength->Fill((*track_length)[Short_ID]);
          ShortStartZ->Fill((*track_zpos)[Short_ID][Short_StartPoint]);
          ShortEndZ->Fill((*track_zpos)[Short_ID][Short_EndPoint]);


          if(zProjLength > 4) {continue;}
          numzProjShort++;

          if((*track_zpos)[Short_ID][Short_StartPoint] > bendZcut ){ continue;}
          numBendZcut++;

          if((*track_xpos)[Short_ID][Short_StartPoint] < (*track_xpos)[Short_ID][Short_EndPoint]){continue;}

          else{

          numRightSided++;
          
          found_bending = true;
          double XYZMin = 999;
          int Long_ID;
          int Long_StartPoint;

          bool found_long = false;

              for (int icandidate = 0; icandidate < matchCandidates.size(); icandidate++){

                int startPoint = static_cast <int> (matchCandidates[icandidate][1]);
                int genID = static_cast <int> (matchCandidates[icandidate][0]);
                if (genID ==  Short_ID){continue;}
                else{
                    //if  ((*track_zpos)[genID][startPoint] > (*track_zpos)[Short_ID][Short_EndPoint]){
                      found_long = true;

                      double delX = (*track_xpos)[Short_ID][Short_EndPoint] - (*track_xpos)[genID][startPoint];
                      double delY = (*track_ypos)[Short_ID][Short_EndPoint] - (*track_ypos)[genID][startPoint];
                      double delZ = (*track_zpos)[Short_ID][Short_EndPoint] - (*track_zpos)[genID][startPoint];

                      double delMatch =  sqrt( pow(delX,2) + pow(delY,2) + pow(delZ,2));


                      

                      if (delMatch < XYZMin){
                        XYZMin = delMatch;
                        LongCandidate = matchCandidates[icandidate];
                        Long_ID = genID;
                        Long_StartPoint =  startPoint;
                    }
                  //}
                }
              }

          if(!found_long){continue;}
          numHasLong++;
          BendingProxHist->Fill(XYZMin);
          if (XYZMin > UI->branchMaxDist){continue;}
          numXYZMatching++;
          //if( abs((*track_zpos)[Long_ID][Long_StartPoint] - (*track_zpos)[Short_ID][Short_EndPoint]) > 4 ){continue;}
          //numZmatching++;

          //double CandidateDist = UtilityFunctions::pointDistance((*track_xpos)[Short_ID][Short_EndPoint],
                                                      //(*track_ypos)[Short_ID][Short_EndPoint],
                                                      //(*track_zpos)[Short_ID][Short_EndPoint],
                                                      //(*track_xpos)[genID][startPoint],
                                                      //(*track_ypos)[genID][startPoint],
                                                      //(*track_zpos)[genID][startPoint]);
          
          //if (DistanceMin < UI || (*track_zpos)[Long_ID][Long_StartPoint] < (*track_zpos)[Short_ID][Short_EndPoint]){
           // continue;
          //}

      

          ///####### at last, we begin plotting projections!!!

          // ...
          // How do I even do this?


          //TF1 fit, taking only X, Z into account.


          // finding interaction to avoid
          //reco_primary =  Long_ID;
          //double temp[6];
          //double* candidate_info = ES->findInt(temp, reco_primary, ntracks_reco, 
                                            //ntrack_hits, track_xpos, track_ypos, track_zpos,
                                            //track_end_x, track_end_y, track_end_z,
                                            //col_track_hits, col_track_dedx, col_track_pitch_hit,
                                            //col_track_x, col_track_y, col_track_z, ESoptions);

          //double interaction_x = candidate_info[1];
          //double interaction_z =  candidate_info[3];

          double zmaxFit = (*track_zpos)[Long_ID][Long_StartPoint] + 4;

          TGraph *BendGraph =  new TGraph(); // X as a function of Z

          for (int ipoint = 0; ipoint < (*ntrack_hits)[Long_ID] ; ipoint++){
            BendGraph->SetPoint(ipoint,(*track_zpos)[Long_ID][ipoint],(*track_xpos)[Long_ID][ipoint]);
          }
          
          TF1 *FitFz = new TF1("FitFz","[0]+[1]*x",0, zmaxFit);
          FitFz->SetParameters(22,1);
          c->cd(1);

          TGraph *ShortGraph =  new TGraph();

          for (int ipoint = 0; ipoint < (*ntrack_hits)[Short_ID] ; ipoint++){
            ShortGraph->SetPoint(ipoint,(*track_zpos)[Short_ID][ipoint],(*track_xpos)[Short_ID][ipoint]);
          }
          ShortGraph->SetMarkerStyle(17);
          ShortGraph->SetMarkerColor(3);
          BendGraph->SetMarkerStyle(17);
          BendGraph->SetMarkerColor(4);
          

          BendGraph->GetXaxis()->SetLimits(-1,UI->zTPCCutoff);
          BendGraph->GetHistogram()->SetMinimum((wctrk_XFace[0] - xMeanTPCentry - UI->rCircleCut));
          BendGraph->GetHistogram()->SetMaximum((wctrk_XFace[0] - xMeanTPCentry + UI->rCircleCut));
          BendGraph->Draw("AP");
          BendGraph->GetXaxis()->SetTitle("Z[cm]");
          BendGraph->GetYaxis()->SetTitle("X[cm]");


          BendGraph->Fit("FitFz","R");

          ShortGraph->Draw("P SAME");

          TGraph *OthersGraph = new TGraph();
          
          if (matchCandidates.size() > 2){
            int graphPtBuffer = 0;
            for (int itrack = 0; itrack < ntracks_reco; itrack++ ){
              if (itrack != Short_ID && itrack != Long_ID){
                //std::cout << "which track" << itrack << std::endl;
                
                for(int ipoint = 0; ipoint < (*ntrack_hits)[itrack]; ipoint++){
                  //std::cout << ipoint << "point filled" << std::endl;
                  
                  OthersGraph->SetPoint(graphPtBuffer,(*track_zpos)[itrack][ipoint],
                    (*track_ypos)[itrack][ipoint]);
                  graphPtBuffer++;
                }
              }
            }
          }
          OthersGraph->SetMarkerStyle(17);
          OthersGraph->SetMarkerColor(12);
          OthersGraph->Draw("P SAME");

          c->Update();

          
          char graph_title[100];
          sprintf(graph_title,"plotting/images/Bending/Graph%d.png",event);
          c->Update();
          c->Print(graph_title,"png");




          for(int ipoint = 0; ipoint < (*ntrack_hits)[Short_ID]; ipoint++){
            double xval  = (*track_xpos)[Short_ID][ipoint];
            double zval = (*track_zpos)[Short_ID][ipoint];
            double projX = FitFz->Eval(zval);
            double deltaX = xval - projX;
            BendingZHist->Fill(zval,deltaX);
            BendingXHist->Fill(projX,deltaX);
          }



          delete BendGraph;
          delete ShortGraph;


        }

        }
        }
      }
    }


    // ### porting over the work from ProtonAnalyzerMC module -- ryan ###
    // ## grabbing reco primary ##
    //int reco_primary = -1;
 
 

    std::cout << "\n------- Cut Results -------\n"<< std::endl;

    std::cout << "Total events processed: "<< numEventsStart << std::endl;
    std::cout << "Events with only one wc track: "<< numWCTrack << std::endl;
    std::cout << "Events with valid ToF value: "<< numtofvalid << std::endl;
    std::cout << "Events passing mass cut: "<< numMassCut << std::endl;
    std::cout << "Events with at least 1 TPC track Z < "  << UI->zTPCCutoff << ": " << numZcutoff << std::endl;
    std::cout << "Events with at least 1 Candidate in Cylinder cut: " << numHasCandidate << std::endl;
    std::cout << "Events with 2 or more Candidates in Cylinder cut: " << numTwoTracks << std::endl;
    std::cout << "Events with Short track Z proj < 4 :" << numzProjShort << std::endl;
    std::cout << "Events with Short track start < "<< UI->zBeamCutoff << ": " << numBendZcut << std::endl;
    std::cout << "Events with track tilt to the Right in XZ: "<< numRightSided  << std::endl;
    std::cout << "Events with Long track found "<< numHasLong << std::endl;
    std::cout << "Events with Long track matched in XY: "<< numXYZMatching << std::endl;
    //std::cout << "Events with long matched in Z "<< numZmatching << std::endl;




    
    

    




  //if(UI->plotIndividualSet){
    //gtotal_res_dedx_item = TGraph(total_res_v.size(),&total_res_v[0],&total_dedx_v[0]);
    //gtotal_res_dedx_item.SetName("gtotal_res_dedx");
    //gtotal_res_dedx = &gtotal_res_dedx_item;
  //}
  






  if (!(isMC)){
    
    if(UI->rootOutputFileSet){
      if(verbose){std::cout << "Writing to outputFile" << std::endl;}
      outputFile->cd();
        
    
      ShortLength->Write();
      ShortStartZ->Write();
      ShortEndZ->Write();
      CylNumHist->Write();


      //tpcInTracksXY->Write();
      //tpcInTracksZ->Write();
      //InTrackLength->Write();
      //wctrkPositionXY->Write();
      //wctrkSelectedXY->Write();
      //wctrkNumHist->Write();
      //inTracksNumHist->Write();
      //BeamMomentum->Write();
      //BeamToF->Write();
      delXYHist->Write();
      //delThetaHist->Write();
      //delPhiHist->Write();
      //BadTrackHist->Write();
      //BadTrackLength->Write();
      //BadTrackStartZ->Write();
      //delBadTrackHist->Write();
      //BeamMassHist->Write();
      //BeamMassCutHist->Write();
      //tofMomentHist->Write();
      //numTracksSelHist->Write();
      //tpcPhiHist->Write();
      //wcPhiHist->Write();
      //tpcThetaHist->Write();
      //wcThetaHist->Write();
      //primary_dedx->Write();

      //## entrering kink diagnostics

      //InTrackTPCnum_out->Write();
      //InTrackTPCnum_in->Write();
      //zProjTrack_in->Write();
      //zProjTrack_out->Write();
      zProjTrack_short->Write();
      //InTrackLength_in->Write();
      //InTrackLength_out->Write();
      if(verbose){std::cout << "Writing bending histos" << std::endl;}
      BendingXHist->Write();
      BendingZHist->Write();
      BendingProxHist->Write();



      // ## EventSelector histos ##
      //BranchDistHist->Write();
      //ClusterDistHist->Write();
      
      //if(UI->plotIndividualSet){gtotal_res_dedx->Write();}
      
      // ## xs histos ##
      //hreco_initialKE->Write();
     // hreco_incke->Write();
      //hreco_intke->Write();
    }
  }
  //if(UI->SelEventListSet){IDfile.close();}

}




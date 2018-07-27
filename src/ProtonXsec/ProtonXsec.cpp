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

  //output list of selected primary particle IDs

  if (UI->SelEventListSet){
    
    IDfile.open(UI->SelEventList, ios::trunc);
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
  double numMassCut = 0;
  double numtofvalid = 0;
  double numZcutoff = 0;
  double numWCTrack = 0;
  double xyDeltaCut = 0;
  double numThetaCut = 0;
  double numPhiCut = 0;


  const double wc_zpos_val = 0.1;
  const double *wc_zpos;
  wc_zpos = &wc_zpos_val;



  // For individual event plotting mode (plotIndividual)

  std::vector<double> total_dedx_v;
  std::vector<double> total_res_v;
  TGraph gtotal_res_dedx_item;
  TGraph * gtotal_res_dedx;



  // ### some variables that are needed for the xs calc ###
  double z2 = UI->zSlabSize; //<-- slab size. need to move this to jobOptions

  
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
  c->Divide(1,2);

  TCanvas *csplit = new TCanvas("csplit","csplit",1500, 750);
  csplit->Divide(2,2);

// ## Histograms for entering tracks
  TH2D *tpcInTracksXY =  new TH2D("tpcInTracksXY","Position of TPC track start XY",
    200, -100, 100, 200, -100, 100);
  TH2D *wctrkPositionXY =  new TH2D("wctrkPositionXY","Position of Wire Chamber Track",
    200,-100 ,100 , 200, -100, 100);
  TH2D *wctrkSelectedXY =  new TH2D("wctrkSelectedXY","Position of Wire Chamber Track",
  200, -100, 100, 200, -100, 100);
  TH1D *tpcInTracksZ =  new TH1D("tpcInTracksZ","Position of TPC track start Z",10, 0, 4);
  TH1D *PrimaryStartZ =  new TH1D("PrimaryStartZ","Selected track start in Z",10, 0, 4);
  TH1D *BadTrackStartZ =  new TH1D("BadTrackStartZ","non-selected track start in Z",10, 0, 4);
  TH1D *InTrackLength =  new TH1D("InTrackLength","Entering Track Length",100, 0, 100);
  TH1D *PrimaryLength =  new TH1D("PrimaryLength","Selected Entering Track Length",100, 0, 100);
  TH1D *BadTrackLength =  new TH1D("BadTrackLength","non-selected Entering Track Length",250, 0, 100);
  TH1D *inTracksNumHist = new TH1D("inTracksNumHist","Number of Entering Tracks TPC",100,0,10);
  TH1D *BeamMomentum = new TH1D("BeamMomentum","Incoming Particle Momentum",200,0,2000);
  TH1D *wctrkNumHist = new TH1D("wctrkNumHist","Number of Entering Tracks TPC",20,0,5);
  TH1D *BeamToF = new TH1D("BeamToF","Incoming Particle Momentum",100,0,100);

  TH2D * delXYHist =  new TH2D("delXYHist","tpc to wc delta x",200,-100,100,200,-100,100);
  TH2D *BadTrackHist =  new TH2D("BadTrackHist","non-selected tracks",200,-100,100,200,-100,100);
  TH2D *delBadTrackHist =  new TH2D("delBadTrackHist","non-selected tracks",200,-100,100,200,-100,100);
  TH2D *delPileupHist =  new TH2D("delPileupHist","Pileup Track difference",200,-100,100,200,-100,100);
  TH2D *PileupHist =  new TH2D("PileupHist","Pileup Track Position",200,-100,100,200,-100,100);
  TH1D * delThetaHist =  new TH1D("delThetaHist","tpc to wc delta Theta",100,-1,1);
  TH1D * tpcThetaHist =  new TH1D("tpcThetaHist","Tpc #theta",100,-1,1);
  TH1D * wcThetaHist =  new TH1D("wcThetaHist","Wire chamber #theta",100,-1,1);
  TH1D * delPhiHist =  new TH1D("delPhiHist"," Tpc to WC #Delta #phi",100,-4,4);
  TH1D * tpcPhiHist =  new TH1D("tpcPhiHist","Tpc #phi",100,0,7);
  TH1D * wcPhiHist =  new TH1D("wcPhiHist","Wire chamber #phi",100,0,7);

  TH1D * numTracksSelHist =  new TH1D("numTracksSelHist","number of Entering Tracks - Selected Events", 10, 0, 5);

  TH2D *tofMomentHist = new TH2D("tofMomentHist","Momentum vs TOF",100,0,2000, 100 , 0,100);
  TH1D *BeamMassHist = new TH1D("BeamMassHist","Beamline particle Mass", 100, 0,3000);
  TH1D *BeamMassCutHist = new TH1D("BeamMassCutHist","Beamline particle Mass - after Cut", 100, 0,3000);
  TH1D *primary_dedx = new TH1D("primary_dedx","primary track dE/dx", 400, 0,40);
  


  // ## plot markers for cuts

  TLine  *MassMinLine =  new TLine(0,UI->MassCutMin, 100, UI->MassCutMin);
  TLine  *MassMaxLine =  new TLine(0,UI->MassCutMax, 100, UI->MassCutMax);

  // ## xs histos ##
  TH1D *hreco_initialKE = new TH1D("hreco_initialKE", "initial ke", 20, 0, 1000);
  TH1D *hreco_intke = new TH1D("hreco_intke", "int ke", 20, 0, 1000);
  TH1D *hreco_incke = new TH1D("hreco_incke", "inc ke", 20, 0, 1000);
  TH1D *hreco_xs = new TH1D("hreco_xs", "P-Ar Inelastic XS", 20, 0, 1000);

  // ## to read in from corrections file
  TH1D *hreco_intke_background = new TH1D("hreco_intke_background", "int ke background", 20, 0, 1000);
  TH1D *hreco_incke_background = new TH1D("hreco_incke_background", "inc ke background", 20, 0, 1000);
  TH2D *hreco_unfolding_matrix_normalized = new TH2D("hreco_unfolding_matrix_normalized", "energy unfolding matrix", 20, 0, 1000, 20, 0, 1000);
  TH1D *hreco_intke_eff = new TH1D("hreco_intke_eff", "interaction selection efficiency", 20, 0, 1000);
  TH1D *hreco_incke_eff = new TH1D("hreco_incke_eff", "incident selection efficiency", 20, 0, 1000);

  if(UI->plotIndividualSet){

  }
  // ## looping once to find Beamline center ##

  if(!isMC){
  for (Long64_t jentry=0; jentry < numEventsToProcess && jentry < nentries; jentry++) {

    
    Long64_t ientry = tuple->LoadTree(jentry);
    if (ientry < 0){continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);
    int numEnteringTracks = 0;
    int best_candidate = -1;
    wctrkNumHist->Fill(num_wctracks);

      if (num_wctracks == 1){
          double ParticleMass = -9999999. ;
          
          bool isProton = BS->MassCut(wctrk_momentum[0], tofObject[0], ParticleMass, UI->MassCutMin, UI->MassCutMax);
          

          if(!isProton){continue;}

        wctrkPositionXY->Fill(wctrk_XFace[0],wctrk_YFace[0]);
        BeamMomentum->Fill(wctrk_momentum[0]);
        BeamToF->Fill(tofObject[0]);

        std::vector <double> wctpc_mvect = BS->dataTPCPrimary(wctrk_XFace[0],wctrk_YFace[0], wctrk_theta[0],
          wctrk_phi[0],track_xpos,track_ypos, track_zpos, ntracks_reco, UI->zTPCCutoff,
          best_candidate, numEnteringTracks);



        if(best_candidate != -1){

            delXYHist->Fill(wctpc_mvect[1],wctpc_mvect[2]);
            tpcInTracksXY->Fill((*track_xpos)[best_candidate][0],
            (*track_ypos)[best_candidate][0]);
                     
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
    bool found_primary = false; 
    if(!isMC){
      if (num_wctracks !=1){continue;}
      else{
          numWCTrack++;
          if(tofObject[0] < 0){continue;}
          numtofvalid++;
          double ParticleMass = -9999999. ;
          tofMomentHist->Fill(wctrk_momentum[0], tofObject[0]);
          bool isProton = BS->MassCut(wctrk_momentum[0], tofObject[0], ParticleMass, UI->MassCutMin, UI->MassCutMax);
          BeamMassHist->Fill(ParticleMass);
          if(!isProton){
              for(int inTrack = 0; inTrack < ntracks_reco; inTrack++ ){
              delPileupHist->Fill(wctrk_XFace[0] - (*track_xpos)[inTrack][0],
                          wctrk_YFace[0] - (*track_ypos)[inTrack][0]);
              PileupHist->Fill((*track_xpos)[inTrack][0],
                            (*track_ypos)[inTrack][0]);
          }
            continue;}



          BeamMassCutHist->Fill(ParticleMass);
          numMassCut++;
          std::vector <double> MinVector = BS->dataTPCPrimary(wctrk_XFace[0],wctrk_YFace[0], wctrk_theta[0], wctrk_phi[0],track_xpos,
          track_ypos, track_zpos, ntracks_reco, UI->zTPCCutoff,reco_primary, numEnteringTracks);

          inTracksNumHist->Fill(numEnteringTracks);
          
          if(reco_primary == -1){
            if(verbose == 2){std::cout << "No Valid Primary \n" << std::endl;}
          }
          else{ 
            numZcutoff++;
              if((sqrt(pow((MinVector[1] - xMeanTPCentry),2) + pow((MinVector[2] - yMeanTPCentry),2)))> UI-> rCircleCut) {
                  if(verbose == 2){std::cout << "No Valid Primary \n" << std::endl;}
                  continue;}
              else{
                  xyDeltaCut++;
                  int numTracksSel = 0;
                for (int inTrack = 0; inTrack < ntracks_reco; inTrack++ ){
                  if((*track_zpos)[inTrack][0] < UI->zTPCCutoff){
                    numTracksSel++;
                    tpcInTracksZ->Fill((*track_zpos)[inTrack][0]);
                    InTrackLength->Fill((*track_length)[inTrack]);
                    
                    if (inTrack != reco_primary){
                      delPileupHist->Fill(wctrk_XFace[0] - (*track_xpos)[inTrack][0],
                          wctrk_YFace[0] - (*track_ypos)[inTrack][0]);
                      PileupHist->Fill((*track_xpos)[inTrack][0],
                          (*track_ypos)[inTrack][0]);
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
                numTracksSelHist->Fill(numTracksSel);
                tpcPhiHist->Fill(UtilityFunctions::getTrackPhi(reco_primary,track_xpos,track_ypos));
                if (wctrk_phi[0] < 0){
                wcPhiHist->Fill((wctrk_phi[0]+ 8 * atan(1)));
                }
                else{wcPhiHist->Fill(wctrk_phi[0]);}
                tpcThetaHist->Fill(UtilityFunctions::getTrackTheta(reco_primary,track_xpos,track_ypos,track_zpos));
                wcThetaHist->Fill(wctrk_theta[0]);
                delThetaHist->Fill(MinVector[3]);
                delPhiHist->Fill(MinVector[4]);
                wctrkSelectedXY->Fill(wctrk_XFace[0],wctrk_YFace[0]);
                if(MinVector[3] > UI->ThetaCut){continue;}
                else{numThetaCut++;
                  if(MinVector[4] > UI->PhiCut){continue;}
                  else{numPhiCut++;
                    if(UI->SelEventListSet){IDfile << jentry << "," << reco_primary << std::endl;}
                  }
                }
              }
            }
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



    if(found_primary){                    
      if(UI->plotIndividualSet){
        double res_buffer = 0;
        int primary_pts = (*col_track_hits)[reco_primary];
        double_t dedx_graphpts[primary_pts], res_graphpts[primary_pts];
        if (primary_pts > 0){
          for(int ipos = primary_pts - 1; ipos >= 0  ; ipos--){

            //std::cout << ipos << " , "<< (*col_track_dedx)[reco_primary][ipos]<< std::endl;
            primary_dedx->Fill((*col_track_dedx)[reco_primary][ipos]);
            double res_range;
            if (ipos == primary_pts -1){res_range = 0;}
            else{
             res_range =  UtilityFunctions::pointDistance(
              (*col_track_x)[reco_primary][ipos],
              (*col_track_y)[reco_primary][ipos],
              (*col_track_z)[reco_primary][ipos],
              (*col_track_x)[reco_primary][ipos + 1],
              (*col_track_y)[reco_primary][ipos + 1],
              (*col_track_z)[reco_primary][ipos + 1]);}

            res_buffer += res_range;
            res_graphpts[ipos] = res_buffer;
            dedx_graphpts[ipos] = (*col_track_dedx)[reco_primary][ipos];
            total_res_v.push_back(res_buffer);
            total_dedx_v.push_back((*col_track_dedx)[reco_primary][ipos]);

          }
          int primary_size = (*ntrack_hits)[reco_primary];

          std::vector<double> ptracks_vx;
          ptracks_vx.reserve(primary_size);
          std::vector<double> ptracks_vy;
          ptracks_vy.reserve(primary_size);
          std::vector<double> ptracks_vz;
          ptracks_vz.reserve(primary_size);

            ptracks_vx.insert(ptracks_vx.end(), (*track_xpos)[reco_primary].begin(), (*track_xpos)[reco_primary].end());
            ptracks_vy.insert(ptracks_vy.end(), (*track_ypos)[reco_primary].begin(), (*track_ypos)[reco_primary].end());
            ptracks_vz.insert(ptracks_vz.end(), (*track_zpos)[reco_primary].begin(), (*track_zpos)[reco_primary].end());



          TGraph *EventXZprimary = new TGraph(primary_size,
            &ptracks_vz[0], &ptracks_vx[0]);
          EventXZprimary->SetName("EventXZprimary");

          TGraph *EventYZprimary = new TGraph(primary_size,
            &ptracks_vz[0], &ptracks_vy[0]);
          EventYZprimary->SetName("EventYZprimary");


          TGraph2D *Event3dPrimary =  new TGraph2D(primary_size,&ptracks_vz[0],&ptracks_vx[0],&ptracks_vy[0]);
          Event3dPrimary->SetName("Event3dPrimary");



          int secondaries_size = 0;
          for (int itrack = 0; itrack < ntracks_reco && itrack !=reco_primary; itrack++ ){
            
            int ctrack_size = (*ntrack_hits)[itrack];
            secondaries_size += ctrack_size;
          }

          std::vector<double> otracks_vx;
          otracks_vx.reserve(secondaries_size);
          std::vector<double> otracks_vy;
          otracks_vy.reserve(secondaries_size);
          std::vector<double> otracks_vz;
          otracks_vz.reserve(secondaries_size);


          TGraph2D *Event3dOther = new TGraph2D();
          Event3dOther->SetName("Event3dOther");
          if (secondaries_size > 0){
            int graphPtBuffer = 0;
            for (int itrack = 0; itrack < ntracks_reco && itrack !=reco_primary; itrack++ ){
              otracks_vx.insert(otracks_vx.end(), (*track_xpos)[itrack].begin(), (*track_xpos)[itrack].end());
              otracks_vy.insert(otracks_vy.end(), (*track_ypos)[itrack].begin(), (*track_ypos)[itrack].end());
              otracks_vz.insert(otracks_vz.end(), (*track_zpos)[itrack].begin(), (*track_zpos)[itrack].end());
              
              for(int ipoint = 0; ipoint < (*ntrack_hits)[itrack]; ipoint++){
                //std::cout << ipoint <<"point filled" << std::endl;
                Event3dOther->SetPoint(graphPtBuffer,(*track_zpos)[itrack][ipoint],
                  (*track_xpos)[itrack][ipoint],(*track_ypos)[itrack][ipoint]);
                graphPtBuffer++;
              }
            }
          }

          TGraph *EventXZother = new TGraph( secondaries_size, &otracks_vz[0], &otracks_vx[0]);
          EventXZother->SetName("EventXZother");
          TGraph *EventYZother = new TGraph( secondaries_size, &otracks_vz[0], &otracks_vy[0]);
          EventYZother->SetName("EventYZother");
          TGraph *EventXZwc = new TGraph(1);
          EventXZwc->SetName("EventXZwc");
          TGraph *EventYZwc = new TGraph(1);
          EventYZwc->SetName("EventYZwc");
          TGraph2D *Event3dWC = new TGraph2D(1);  
          Event3dWC->SetName("Event3dWC");

          if(!isMC){
            EventXZwc->SetPoint(0,0,wctrk_XFace[0]);
            EventYZwc->SetPoint(0,0,wctrk_YFace[0]);
            Event3dWC->SetPoint(0,0,wctrk_XFace[0],wctrk_YFace[0]);
          }

          EventXZprimary->SetMarkerStyle(7);
          EventXZprimary->SetMarkerColor(2);
          EventYZprimary->SetMarkerStyle(7);
          EventYZprimary->SetMarkerColor(2);

          EventXZother->SetMarkerStyle(7);
          EventXZother->SetMarkerColor(4);
          EventYZother->SetMarkerStyle(7);
          EventYZother->SetMarkerColor(4);
          

          EventXZwc->SetMarkerStyle(9);
          EventXZwc->SetMarkerColor(3);
          EventYZwc->SetMarkerStyle(9);
          EventYZwc->SetMarkerColor(3);


          csplit->cd(1);
          EventXZprimary->SetTitle("");
          EventXZprimary->GetXaxis()->SetTitle("Z [cm]");
          EventXZprimary->GetYaxis()->SetTitle("X [cm]");
          if (primary_size > 0){EventXZprimary->Draw("AP");}
          if (secondaries_size > 0){EventXZother->Draw("psame");}
          if (!isMC){EventXZwc->Draw("psame");}

          csplit->cd(3);
          EventYZprimary->SetTitle("");
          EventYZprimary->GetXaxis()->SetTitle("Z [cm]");
          EventYZprimary->GetYaxis()->SetTitle("Y [cm]");
          if (primary_size > 0){EventYZprimary->Draw("AP");}
          if (secondaries_size > 0){EventYZother->Draw("psame");}
          if (!isMC){EventYZwc->Draw("psame");}

          csplit->cd(2);
          TGraph *gres_dedx = new TGraph(primary_pts,res_graphpts,dedx_graphpts);
          gres_dedx->GetXaxis()->SetTitle("Residual Range [cm]");
          gres_dedx->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
          gres_dedx->SetTitle("");
          gres_dedx->SetMarkerStyle(20);
          gres_dedx->SetMarkerSize(1.5);
          gres_dedx->SetMarkerColor(4);

          gres_dedx->Draw("AP");
          csplit->cd(4);
          
          TGraph2D *Event3dTPC = new TGraph2D();
          Event3dTPC->SetName("Event3dTPC");
          Event3dTPC->SetPoint(0,0,0,-20);
          Event3dTPC->SetPoint(1,90,48,20);
          Event3dTPC->SetTitle("3D Event Display");
          Event3dTPC->GetXaxis()->SetTitle("Z [cm]");
          Event3dTPC->GetYaxis()->SetTitle("X [cm]");
          Event3dTPC->GetZaxis()->SetTitle("Y [cm]");
          Event3dTPC->Draw("p");
          Event3dPrimary->Draw("P SAME");
          
          if(secondaries_size > 0){
            Event3dOther->SetMarkerStyle(7);
            Event3dOther->SetMarkerColor(4);

            Event3dOther->Draw("P SAME");}
          if (!isMC){

            Event3dWC->Draw("P SAME");}

          Event3dPrimary->SetMarkerStyle(7);
          Event3dPrimary->SetMarkerColor(2);
          Event3dWC->SetMarkerStyle(9);
          Event3dWC->SetMarkerColor(3);

          // gPad->Modified(); 
          gPad->Update();

          char eventdisp_title[100];
          sprintf(eventdisp_title,"%s/EventDisplay%d.png",UI->plotIndividual,event);
          
          //char gres_dedx_title[100];
          //sprintf(gres_dedx_title,"%s/dedx_Residuals/gres_dedx%d.png",UI->plotIndividual,event);


          //~gres_dedx;
          csplit->Print(eventdisp_title,"png");
          delete gres_dedx;
          delete EventXZprimary;
          delete EventYZprimary;
          delete EventXZother;
          delete EventYZother;
          delete EventYZwc;
          delete EventXZwc;
          delete Event3dPrimary;
          delete Event3dOther;
          delete Event3dWC;
          delete Event3dTPC;




          
        }
      


      }

    }

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
      }//<-- End if this is the \ slab
    }//<--End calo slab loop

    // ### end of port work -- ryan ###
  
  }//end of event Loop


  // #####################################################################################
  // #####################################################################################
  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
  // ################################### xs calculation ################################## 
  // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- 
  // #####################################################################################
  // #####################################################################################

  // ### Taking raw uncorrected histograms and applying (1) background subtraction
  // ### (2) energy unfolding and (3) efficiency corrections. These corrections are applied
  // ### to both the numerator and denominator. Eventually we should probably use 
  // ### switches for which we apply. For now I'm applying all of them in the order
  // ### described. The default will be to use Bertini model on data, for MC comparisons
  // ### we'll try all of them. -- ryan

  // --- just realized that for now we can only use the background subtraction on MC
  // --- it's a raw number and needs to be scaled for the size of the sample
  // --- for now I'm putting the whole thing in if( isMC ) and we can deal with it later
  if(UI->isMC) {

    // ## should probably incorp this into the joboptions in a clever way -- ryan
    TFile *mc_corrections = new TFile("files/corrections/xsCorrectionBertini.root"); 

    hreco_intke_background = (TH1D*)mc_corrections->Get("hreco_intke_background");
    hreco_incke_background = (TH1D*)mc_corrections->Get("hreco_incke_background");
    hreco_unfolding_matrix_normalized = (TH2D*)mc_corrections->Get("hreco_unfolding_matrix_normalized");
    hreco_intke_eff = (TH1D*)mc_corrections->Get("hreco_int_eff");
    hreco_incke_eff = (TH1D*)mc_corrections->Get("hreco_inc_eff");

    // ## temporary histograms ##
    TH1D *hreco_folded_intke_signal = new TH1D("hreco_folded_intke_signal", "total-back", 20, 0, 1000);
    TH1D *hreco_folded_incke_signal = new TH1D("hreco_folded_incke_signal", "total-back", 20, 0, 1000);
    TH1D *hreco_unfolded_intke_signal = new TH1D("hreco_unfolded_intke_signal", "U(total-back)", 20, 0, 1000);
    TH1D *hreco_unfolded_incke_signal = new TH1D("hreco_unfolded_incke_signal", "U(total-back)", 20, 0, 1000);

    // ## folded signal distributions ##
    for(int iBin = 0; iBin < hreco_intke->GetNbinsX(); iBin++){
      double total = hreco_intke->GetBinContent(iBin);
      double background = hreco_intke_background->GetBinContent(iBin);
      hreco_folded_intke_signal->SetBinContent(iBin, total-background); 
    }
    for(int iBin = 0; iBin < hreco_incke->GetNbinsX(); iBin++){
      double total = hreco_incke->GetBinContent(iBin);
      double background = hreco_incke_background->GetBinContent(iBin);
      hreco_folded_incke_signal->SetBinContent(iBin, total-background); 
    }

    // ## unfolding signal distributions ##
    for(int iBin = 0; iBin < hreco_unfolding_matrix_normalized->GetNbinsX(); iBin++){
      int n_int_entries = hreco_folded_intke_signal->GetBinContent(iBin);
      int n_inc_entries = hreco_folded_incke_signal->GetBinContent(iBin);
      for(int jBin = 0; jBin < hreco_unfolding_matrix_normalized->GetNbinsY(); jBin++){
        double weight = hreco_unfolding_matrix_normalized->GetBinContent(iBin, jBin);
        double int_value = n_int_entries * weight;
        double inc_value = n_inc_entries * weight;
        hreco_unfolded_intke_signal->AddBinContent(jBin, int_value);
        hreco_unfolded_incke_signal->AddBinContent(jBin, inc_value);
      }
    }

    // ## putting the pieces together ##
    for(int iBin = 0; iBin < hreco_unfolded_intke_signal->GetNbinsX(); iBin++){
      if(hreco_unfolded_incke_signal->GetBinContent(iBin) == 0){continue;}
      if(hreco_intke_eff->GetBinContent(iBin) == 0){continue;}
      if(hreco_incke_eff->GetBinContent(iBin) == 0){continue;}
      // num:   (N_int - Background_int)*U_ij * 1/eps
      double num = hreco_unfolded_intke_signal->GetBinContent(iBin) / hreco_intke_eff->GetBinContent(iBin);
      if(num == 0){num = 1;}
      // denom: (N_inc - Background_inc)*U_ij * 1/eps
      double dem = hreco_unfolded_incke_signal->GetBinContent(iBin) / hreco_incke_eff->GetBinContent(iBin);
      if(dem == 0){dem =1;}

      // # ratio #
      double ratio = num / dem;
      double temp_xs = ratio * sparse_recip_num_density;
      double xs = temp_xs / barn;
    
      // # error #
      double num_err = pow(num, .5);
      double term1 = num_err/num;
      //double dem_err = pow(sincke->GetBinContent(iBin), .5);
      double dem_err = pow(dem, .5);
      double term2 = dem_err/dem;
      double totalError = temp_xs*pow(pow(term1,2) + pow(term2,2), 0.5)/barn;//*recip_num_density*barn;

      std::cout<<"iBin: "<<iBin<<std::endl;
      std::cout<<"\tnum: "<<num<<std::endl;
      std::cout<<"\tdem: "<<dem<<std::endl;
      std::cout<<"\tratio: "<<ratio<<std::endl;
      std::cout<<"\txs: "<<xs<<" +- "<<totalError<<std::endl;

      hreco_xs->SetBinContent(iBin, xs); 
      hreco_xs->SetBinError(iBin,totalError);
    }

  }//<--End if is MC for xs calculation
  


  if (!isMC){
    if(verbose > 0){
    std::cout << "\n------- Beam Selection Results -------\n"<< std::endl;
    std::cout << "Total events processed: "<< numEventsStart << std::endl;
    std::cout << "Events with only one wc track: "<< numWCTrack << std::endl;
    std::cout << "Events with valid ToF value: "<< numtofvalid << std::endl;
    std::cout << "Events passing mass cut: "<< numMassCut << std::endl;   
    std::cout << "Events with best match TPC track Z < Cutoff: "<< numZcutoff << std::endl;
    std::cout << "Events after XY plane distance cut: "<< xyDeltaCut << std::endl;
    std::cout << "Events after Phi cut: " << numPhiCut << std::endl;
    std::cout << "Events after Theta cut: " << numThetaCut << std::endl;
    }
  }


  if(UI->plotIndividualSet){
    gtotal_res_dedx_item = TGraph(total_res_v.size(),&total_res_v[0],&total_dedx_v[0]);
    gtotal_res_dedx_item.SetName("gtotal_res_dedx");
    gtotal_res_dedx = &gtotal_res_dedx_item;
  }
  
  
  if(isMC){

    if(UI->rootOutputFileSet){
      outputFile->cd();
       
      // ## xs histos ##
      hreco_initialKE->Write();
      hreco_incke->Write();
      hreco_intke->Write();
      hreco_xs->Write();
  }


  }





  else if (!(isMC)){
    
    if(UI->rootOutputFileSet){
      outputFile->cd();
    
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
      BeamMassHist->Write();
      BeamMassCutHist->Write();
      tofMomentHist->Write();
      delPileupHist->Write();
      PileupHist->Write();
      numTracksSelHist->Write();
      tpcPhiHist->Write();
      wcPhiHist->Write();
      tpcThetaHist->Write();
      wcThetaHist->Write();
      primary_dedx->Write();
      if(UI->plotIndividualSet){gtotal_res_dedx->Write();}


      // ## xs histos ##
      hreco_initialKE->Write();
      hreco_incke->Write();
      hreco_intke->Write();
    }
  }
  if(UI->SelEventListSet){IDfile.close();}

}




//
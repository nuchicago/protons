//#############################################################################
//
// ProtonXsec.cpp
//
//  Calculates Inelastic cross section from either MC or Data
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
  gStyle->SetOptStat("emr");

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
    //char LogTitle[500];
    //sprintf(LogTitle,"%s.txt", UI->rootOutputFile);
    //RunLog.open()
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
    IDfile.close();
  }

  if (UI->MultiMatchEventListSet){
    multiMatchFile.open(UI->MultiMatchEventList, ios::trunc);
    multiMatchFile.close();

  }

  // output for data driven monte carlo

  if (UI->beamCharFileSet){
    beamPlotFile = new TFile( UI->beamCharFile, "RECREATE");
  }
  if(UI->haloCharFileSet){
    haloPlotFile = new TFile( UI->haloCharFile, "RECREATE");
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

void ProtonXsec::AnalyzeFromNtuples(){

  // counters for beam cuts

  double numEventsStart = 0;
  double numMassCut = 0;
  int numtofsingle = 0;
  int  numMatchedMain = 0;
  double numtofvalid = 0;
  double numZcutoff = 0;
  double numWCTrack = 0;
  double num4PointTrack = 0;
  double numPickyTrack = 0;
  double numQualityTrack = 0;
  double xyDeltaCut = 0;
  double numThetaCut = 0;
  double numPhiCut = 0;
  bool intHistFilled = false;




  const double wc_zpos_val = 0.1;
  const double *wc_zpos;
  wc_zpos = &wc_zpos_val;



  int numInteractions = 0;

  // For individual event plotting mode (plotIndividual)

  //std::vector<double> total_dedx_v;
  //std::vector<double> total_res_v;
  //TGraph gtotal_res_dedx_item;
  //TGraph * gtotal_res_dedx;


    if (UI->MultiMatchEventListSet){
    multiMatchFile.open(UI->MultiMatchEventList, ios::out);
    multiMatchFile << "RUN\tSUBRUN\tEVENT\tNUM_MATCHES\n";
    multiMatchFile.close();

  }

  // #### Data Driven MC 


  int beam_run;
  int beam_subrun;
  int beam_event;
  double beam_x;
  double beam_y;
  double beam_z;
  double beam_angle_xz;
  double beam_angle_yz;
  double beam_momentum;

  std::vector< double > halo_pileup_x;
  std::vector< double > halo_pileup_y;
  std::vector< double > halo_pileup_z;
  std::vector< double > halo_pileup_angle_xz;
  std::vector< double > halo_pileup_angle_yz;
  int halo_pileup_run;
  int halo_pileup_subrun;
  int halo_pileup_event;
  //std::vector< double > halo_pileup_momentum;
  int halo_pileup_number_particles;



  TTree * beam_tree;
  if(UI->beamCharFileSet){
    beamPlotFile->cd();
    beam_tree = new TTree("beam", "beam");

    beam_tree->Branch("beam_run", &beam_run, "beam_run/I");
    beam_tree->Branch("beam_subrun", &beam_subrun, "beam_subrun/I");
    beam_tree->Branch("beam_event", &beam_event, "beam_event/I");

    beam_tree->Branch("beam_x", &beam_x, "beam_x/D");
    beam_tree->Branch("beam_y", &beam_y, "beam_y/D");
    beam_tree->Branch("beam_z", &beam_z, "beam_z/D");
    beam_tree->Branch("beam_angle_xz", &beam_angle_xz, "beam_angle_xz/D");
    beam_tree->Branch("beam_angle_yz", &beam_angle_yz, "beam_angle_yz/D");
    beam_tree->Branch("beam_momentum", &beam_momentum, "beam_momentum/D");
  }

  TTree * halo_pileup_tree;
  if(UI->haloCharFileSet){

    haloPlotFile->cd();

    halo_pileup_tree = new TTree("halo_pileup", "halo_pileup");
    halo_pileup_tree->Branch("halo_pileup_run", &halo_pileup_run, "halo_pileup_run/I");
    halo_pileup_tree->Branch("halo_pileup_subrun", &halo_pileup_subrun, "halo_pileup_subrun/I");
    halo_pileup_tree->Branch("halo_pileup_event", &halo_pileup_event, "halo_pileup_event/I");
    halo_pileup_tree->Branch("halo_pileup_x", &halo_pileup_x);
    halo_pileup_tree->Branch("halo_pileup_y", &halo_pileup_y);
    halo_pileup_tree->Branch("halo_pileup_z", &halo_pileup_z);
    halo_pileup_tree->Branch("halo_pileup_angle_xz", &halo_pileup_angle_xz);
    halo_pileup_tree->Branch("halo_pileup_angle_yz", &halo_pileup_angle_yz);
    //halo_pileup_tree->Branch("halo_pileup_momentum", &halo_pileup_momentum);
    halo_pileup_tree->Branch("halo_pileup_number_particles", &halo_pileup_number_particles, "halo_pileup_number_particles/I");


  }



  // ### some variables that are needed for the xs calc ###
  double z2 = UI->zSlabSize; //<-- slab size. need to move this to jobOptions

  // Selector Options in a vector
  std::vector<double> ESoptions = {static_cast<double> (verbose), UI->dedxNoBraggMax,UI->branchMaxDist, UI->clusterMaxDist};
  std::vector<double> BSoptions = {static_cast<double> (verbose), UI->zTPCCutoff, UI->alphaCut, UI->rCircleCut,
    static_cast<double>(UI->pionCuts), UI->numPileupCut, UI->pcPileupDist, UI->numTracksShower, UI->lenTracksShower,static_cast<double>(UI->uniqueMatch) };


  
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
  TH2D *wctrk4XY =  new TH2D("wctrk4XY","Position of Wire Chamber Track on detector 4",
  200, -100, 100, 200, -100, 100);
  TH2D *wctrkTpcXY =  new TH2D("wctrkTpcXY","Position of Wire Chamber Tracks z = 0, after quality cut",
  200, -100, 100, 200, -100, 100);
  TH1D *tpcInTracksZ =  new TH1D("tpcInTracksZ","Position of TPC entering track start Z",25, 0, 10);
  

  TH1D *tpcAllTrackZmin1 = new TH1D("tpcAllTrackZmin1", "Position of TPC track start in Z", 300, -20, 100);
  TH1D *tpcAllTrackZmin2 = new TH1D("tpcAllTrackZmin2", "Position of TPC track start in Z", 300, -20, 100);
  TH1D *tpcAllTrackZmax = new TH1D("tpcAllTrackZmax", "Position of TPC track end in Z", 300, -20, 120);
  TH1D * allAlphaHist = new TH1D("allAlphaHist","Angle between WC and TPC track",360,0,360);
  TH1D * allAlphaHistEnd = new TH1D("allAlphaHistEnd","Angle between WC and TPC track - Using End Point",360,0,360);
  TH1D * alphaHistEnd = new TH1D("alphaHistEnd","Angle between WC and TPC entering track - Using End Point",360,0,360);
  TH1D * dirAngleHist = new TH1D("dirAngleHist","Angle between direction vectors - entering tracks",360,0,360);
  TH1D * allDirAngleHist = new TH1D("allDirAngleHist","Angle between direction vectors - entering tracks",360,0,360);
  TH1D * circleAlphaHist  = new TH1D("circleAlphaHist","Angle between WC and TPC track - Tracks Within Circle Cut",360,0,360);


  TH1D * tpcAllTrackTheta = new TH1D("tpcAllTrackTheta"," #Theta angle with Z axis",360,0,360);
  TH1D * tpcAllTrackNhits = new TH1D("tpcAllTrackNhits"," Number of Hits per Track",600,0,600);
  TH1D * tpcAllTrackLength = new TH1D("tpcAllTrackLength"," Length of All Tracks",600,0,600);
  TH1D * NtracksTotalHist = new TH1D("NtracksTotalHist"," Number of Tracks per event",600,0,600);
  TH1D * tpcAllTrackDelX = new TH1D("tpcAllTrackDelX"," #Delta X For All Tracks",300,-60,60);
  TH1D * tpcAllTrackDelY = new TH1D("tpcAllTrackDelY"," #Delta Y For All Tracks",300,-60,60);
  TH1D * allRdistHist = new TH1D("allRdistHist", "Position of TPC track in #Delta X #DeltaY", 300, -20, 120);




  TH1D *tpcInTrackEndZ = new TH1D("tpcInTrackEndZ","Postion of entering track end in Z", 250, 0, 100);
  TH1D *zProjPrimaryTrack =  new TH1D("zProjPrimaryTrack","Length of primary track - z projection", 250, 0, 100);
  TH1D *zProjBadTrack =  new TH1D("zProjBadTrack","Length of pileup track - z projection", 250, 0, 100);
  TH1D *PrimaryStartZ =  new TH1D("PrimaryStartZ","Selected track start in Z",25, 0, 10);
  TH1D *BadTrackStartZ =  new TH1D("BadTrackStartZ","non-selected track start in Z",25, 0, 10);
  TH1D *InTrackLength =  new TH1D("InTrackLength","Entering Track Length",250, 0, 100);
  TH1D *PrimaryLength =  new TH1D("PrimaryLength","Selected Entering Track Length",250, 0, 100);
  TH1D *BadTrackLength =  new TH1D("BadTrackLength"," Pileup Track Length",250, 0, 100);
  TH1D *BadTrackLength_m =  new TH1D("BadTrackLength_m","Pileup Track Length - match found",250, 0, 100);
  TH1D *inTracksNumHist = new TH1D("inTracksNumHist","Number of Entering Tracks TPC",100,0,10);


  TH1D *electronLifetimeHist = new TH1D("electronLifetimeHist","electron lifetime",300,0, 3000);

  
  TH1D *wctrkNumHist = new TH1D("wctrkNumHist","Number of Entering Tracks TPC",40,0,10);
  TH1D *BeamToF = new TH1D("BeamToF","Incoming Particle ToF",100,0,100);

  TH2D * delXYHist =  new TH2D("delXYHist","tpc to wc delta x",200,-100,100,200,-100,100);
  TH1D * delXYHistPx =  new TH1D("delXYHistPx","tpc to wc delta x",200,-100,100);
  TH1D * delXYHistPy =  new TH1D("delXYHistPy","tpc to wc delta x",200,-100,100);

  TH2D * delXYHistCenter =  new TH2D("delXYHistCenter","tpc to wc delta x",200,-100,100,200,-100,100);
  TH1D * delXYHistPxCenter =  new TH1D("delXYHistPxCenter","tpc to wc delta x",200,-100,100);
  TH1D * delXYHistPyCenter =  new TH1D("delXYHistPyCenter","tpc to wc delta x",200,-100,100);


  TH2D * delXYHistMatch =  new TH2D("delXYHistMatch","tpc to wc delta x",200,-100,100,200,-100,100);
  TH1D * delXYHistPxMatch =  new TH1D("delXYHistPxMatch","tpc to wc delta x",200,-100,100);
  TH1D * delXYHistPyMatch =  new TH1D("delXYHistPyMatch","tpc to wc delta x",200,-100,100);

  TH2D *BadTrackHist =  new TH2D("BadTrackHist","non-selected tracks",200,-100,100,200,-100,100);


  TH2D *delBadTrackHist_mlen =  new TH2D("delBadTrackHist_mlen","pileup - match found, L > 70",200,-100,100,200,-100,100);
  TH2D *delBadTrackHist_len =  new TH2D("delBadTrackHist_len","pileup - L > 70",200,-100,100,200,-100,100);
  TH2D *delBadTrackHist_m =  new TH2D("delBadTrackHist_m","pileup - match found",200,-100,100,200,-100,100);
  TH2D *delBadTrackHist_nocut =  new TH2D("delBadTrackHist_nocut","pileup - no cuts",200,-100,100,200,-100,100);

//  TH1D * delThetaHist =  new TH1D("delThetaHist","tpc to wc delta Theta",100,-1,1);
//  TH1D * tpcThetaHist =  new TH1D("tpcThetaHist","Tpc #theta",100,-1,1);
//  TH1D * wcThetaHist =  new TH1D("wcThetaHist","Wire chamber #theta",100,-1,1);
//  TH1D * delPhiHist =  new TH1D("delPhiHist"," Tpc to WC #Delta #phi",100,-4,4);
//  TH1D * tpcPhiHist =  new TH1D("tpcPhiHist","Tpc #phi",100,0,7);
//  TH1D * wcPhiHist =  new TH1D("wcPhiHist","Wire chamber #phi",100,0,7);

  TH1D * alphaHist = new TH1D("alphaHist","Angle between WC and TPC track",360,0,360);
  TH1D * numPileupTracksHist = new TH1D("numPileupTracksHist","number of Tracks for Pileup Filter",30,0,30);
  TH1D * numShowerCutHist = new TH1D("numShowerCutHist", "number of short tracks - EM shower cut",30,0,30);

  TH1D * numTracksSelHist =  new TH1D("numTracksSelHist","number of Entering Tracks - Selected Events", 10, 0, 5);

  TH2D *tofMomentHist = new TH2D("tofMomentHist","Momentum vs TOF",100,0,2000, 100 , 0,100);
  TH1D *BeamMassHist = new TH1D("BeamMassHist","Beamline particle Mass", 100, -3000,3000);
  TH1D *BeamMassCutHist = new TH1D("BeamMassCutHist","Beamline particle Mass - after Cut", 100, -3000,3000);
  TH1D *BeamToFmassSel =  new  TH1D("BeamToFmassSel","Beamline ToF - after Mass Cut",300, -30, 270);
  TH1D *BeamMomentumMassSel = new TH1D("BeamMomentumMassSel", "Beam Momentum - after Mass Cut", 200, 0, 2000);
  //TH1D *primary_dedx = new TH1D("primary_dedx","primary track dE/dx", 400, 0,40);
  TH1D *BeamMomentum = new TH1D("BeamMomentum","Incoming Particle Momentum",200,0,2000);

  TH1D *beamLengthHist =  new TH1D("beamLengthHist","Calculated beam Length",100,500,900);


  // resolution changes with cuts

  TH1D *BeamMomentumInit = new TH1D("BeamMomentumInit","Incoming Particle Momentum",200,0,2000);
  TH1D *BeamMomentumQual = new TH1D("BeamMomentumQual","Momentum after quality flag",200,0,2000);
  TH1D *BeamMomentumMatch = new TH1D("BeamMomentumMatch","Incoming Particle Momentum",200,0,2000);

  TH1D *BeamMassInit = new TH1D("BeamMassInit","Incoming Particle Mass", 360, -600,3000);
  TH1D *BeamMassQual = new TH1D("BeamMassQual","Mass after quality flag", 360, -600,3000);
  TH1D *BeamMassMatch = new TH1D("BeamMassMatch","Mass after matching", 360, -600,3000);

  //plots for EventSelection Branches

  TH1D * BranchDistHist = new TH1D("BranchDistHist", "Inelastic Event Branch Distance",50, 0,25);
  TH1D * ClusterDistHist = new TH1D("ClusterDistHist", "Additional Branch Distance (Type 4)",50, 0,25);

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

  TH1D *hreco_folded_intke_signal = new TH1D("hreco_folded_intke_signal", "total-back", 20, 0, 1000);
  TH1D *hreco_folded_incke_signal = new TH1D("hreco_folded_incke_signal", "total-back", 20, 0, 1000);
  TH1D *hreco_unfolded_intke_signal = new TH1D("hreco_unfolded_intke_signal", "U(total-back)", 20, 0, 1000);
  TH1D *hreco_unfolded_incke_signal = new TH1D("hreco_unfolded_incke_signal", "U(total-back)", 20, 0, 1000);


  TH1D *beamXtpc0 = new TH1D("beamXtpc0","beam X position, z = 0 cm", 160, -20, 60);
  TH1D *beamYtpc0 = new TH1D("beamYtpc0","beam Y position, z = 0 cm", 160, -40, 40);
  TH1D *beamZtpc0 = new TH1D("beamZtpc0","beam Z at tpc face", 160, -140 , -60);
  // Plots for data driven MC

  //TH1D *beam_x = new TH1D("beam_x","beam X position, z = -100 cm", 160, -20, 60);
  //TH1D *beam_y = new TH1D("beam_y","beam Y position, z = -100 cm", 160, -40, 40);
  //TH1D *beam_z = new TH1D("beam_z","beam Z start", 160, -140 , -60);

  //TH1D *beam_angle_xz =  new TH1D("beam_angle_xz", "beam theta xz", 200, -10, 10);
  //TH1D *beam_angle_yz = new TH1D("beam_angle_yz", "beam theta yz", 200 , -10, 10);
  //TH1D *beam_momentum = new TH1D("beam_momentum","beam momentum after quality flag",200,0,2000);


  //TH1D *halo_pileup_x = new TH1D("halo_pileup_x","halo pileup X position, z = -1 cm", 160, -20, 60);
  //TH1D *halo_pileup_y = new TH1D("halo_pileup_y","halo pileup Y position, z = -1 cm", 160, -40, 40);
  //TH1D *halo_pileup_z = new TH1D("halo_pileup_z","halo pileup Z start", 160, -40 , -40);

  //TH1D *halo_pileup_angle_xz =  new TH1D("halo_pileup_angle_xz", "halo pileup theta xz", 200, -5, 5);
  //TH1D *halo_pileup_angle_yz = new TH1D("halo_pileup_angle_yz", "halo pileup theta yz", 200 , -5, 5);

  //TH1D *halo_pileup_number_particles = new TH1D("halo_pileup_number_particles","Momentum after quality flag",20,0,20);

  // ## looping once to find Beamline center ##

  if(verbose){std::cout << "num entries in file  = " << nentries << std::endl;}
  if(verbose){std::cout << "num events to process  = " << numEventsToProcess << std::endl;}


  if(!isMC){
  for (Long64_t jentry=0;  jentry < UI->numCenteringEvents && jentry < nentries; jentry++) {
    
    
    Long64_t ientry = tuple->LoadTree(jentry);
    if (ientry < 0){continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);
    
    
    int numEnteringTracks = 0;
    int best_candidate = -1;
    double xMeanTPCentry = 0;
    double yMeanTPCentry = 0;
    bool skipCircleCut = true;


      if (num_wctracks != 1){continue;}
      if (wctrk_missed > 0){continue;}
      else{
          if(UI->pickyTracksWC){if (  wctrk_picky != 1){continue;}}
          if(UI->qualityTracksWC){if (wctrk_quality !=1){continue;}}
          

          double ParticleMass = -9999999.;

          bool isProton = BS->MassCut(wctrk_momentum[0], tofObject[0], UI->beamLength, UI->tofOffset, ParticleMass, UI->MassCutMin, UI->MassCutMax);

          if(applyMassCut){
            if(!isProton){continue;}
          }

        std::vector<double> matchCandidate = BS->BeamCentering(wctrk_x_proj_3cm[0], wctrk_y_proj_3cm[0],
         track_xpos, track_ypos, track_zpos, ntracks_reco, ntrack_hits, UI->zTPCCutoff, best_candidate);

        if(best_candidate != -1){
          int startIndex = static_cast <int> (matchCandidate[1]);
          delXYHist->Fill(wctrk_x_proj_3cm[0] - (*track_xpos)[best_candidate][startIndex],
          wctrk_y_proj_3cm[0] - (*track_ypos)[best_candidate][startIndex]);
          delXYHistPx->Fill(wctrk_x_proj_3cm[0] - (*track_xpos)[best_candidate][startIndex]);
          delXYHistPy->Fill(wctrk_y_proj_3cm[0] - (*track_ypos)[best_candidate][startIndex]);
          wctrkPositionXY->Fill(wctrk_x_proj_3cm[0],wctrk_y_proj_3cm[0]);


          if (matchCandidate[4] < 10 && matchCandidate[5] < BSoptions[2]){
            delXYHistCenter->Fill(wctrk_x_proj_3cm[0] - (*track_xpos)[best_candidate][startIndex],
            wctrk_y_proj_3cm[0] - (*track_ypos)[best_candidate][startIndex]);
            delXYHistPxCenter->Fill(wctrk_x_proj_3cm[0] - (*track_xpos)[best_candidate][startIndex]);
            delXYHistPyCenter->Fill(wctrk_y_proj_3cm[0] - (*track_ypos)[best_candidate][startIndex]);
            }
          }
        } 
      }
    }
  // End of beam centering loop



    // setting beam center values
    double xMeanTPCentry = delXYHistCenter->GetMean(1);
    double yMeanTPCentry = delXYHistCenter->GetMean(2);

    BS->SetMeanXY(xMeanTPCentry,yMeanTPCentry);

  // ## event loop ##
    if(verbose){std::cout << "Starting main loop" << std::endl;
    std::cout << "x Entering Mean :" << xMeanTPCentry<< std::endl;
    std::cout << "y Entering Mean :" << yMeanTPCentry<< std::endl;
  }
  for (Long64_t jentry=0; jentry < numEventsToProcess && jentry < nentries; jentry++){
    
    Long64_t ientry = tuple->LoadTree(jentry); 
    if (ientry < 0){
      continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);
    int numEnteringTracks;
    int reco_primary = -1;
    numEventsStart++;
    if(verbose){printEvent();}
    bool found_primary = false;
    bool passed_geo_cuts = true;
    double ParticleMass = -9999999.;
    wctrkNumHist->Fill(num_wctracks);

    double tertiaryLength = 0;

    if(!(UI->beamPlanesSet)){tertiaryLength = UI->beamLength;}
    else{

      std::vector<double> *wc1 = new std::vector<double> {wctrk_wc1_x[0], wctrk_wc1_y[0], wctrk_wc1_z[0]};
      std::vector<double> *wc2 = new std::vector<double> {wctrk_wc2_x[0], wctrk_wc2_y[0], wctrk_wc2_z[0]};
      std::vector<double> *wc3 = new std::vector<double> {wctrk_wc3_x[0], wctrk_wc3_y[0], wctrk_wc3_z[0]};
      std::vector<double> *wc4 = new std::vector<double> {wctrk_wc4_x[0], wctrk_wc4_y[0], wctrk_wc4_z[0]};




      double mg1Zcoef =  cos(UI->mg1Angle * (pi/180))/sin(UI->mg1Angle * (pi/180));
      double mg2Zcoef =  cos(UI->mg2Angle * (pi/180))/sin(UI->mg2Angle * (pi/180));

      //std::cout << "WC1 : x = " << (*wc1)[0] << "; y = " << (*wc1)[1] << "; z = " << (*wc1)[2] << std::endl;
      //std::cout << "WC2 : x = " << (*wc2)[0] << "; y = " << (*wc2)[1] << "; z = " << (*wc2)[2] << std::endl;


      std::vector<double> magnet1 = UtilityFunctions::pointProjector(wc1, wc2, mg1Zcoef, UI->mg1Const);
      //std::cout << "magnet 1 intercept { " << magnet1[0]<< ", " << magnet1[1]<< ", " << magnet1[2] << " }" << std::endl;

      //std::cout << "WC3 : x = " << (*wc3)[0] << "; y = " << (*wc3)[1] << "; z = " << (*wc3)[2] << std::endl;
      //std::cout << "WC4 : x = " << (*wc4)[0] << "; y = " << (*wc4)[1] << "; z = " << (*wc4)[2] << std::endl;


      std::vector<double> magnet2 = UtilityFunctions::pointProjector(wc3, wc4, mg2Zcoef, UI->mg2Const);
      //std::cout << "magnet 2 intercept { " << magnet2[0]<< ", " << magnet2[1]<< ", " << magnet2[2] << " }" << std::endl;

      double d_wc1_mg1 = sqrt(pow((*wc1)[0] - magnet1[0],2) + pow((*wc1)[1] - magnet1[1],2) +pow((*wc1)[2] - magnet1[2],2));
      //std::cout << "dist wc1-mg1: " << d_wc1_mg1 << "cm" << std::endl;  
      double d_mg1_mg2 =  sqrt(pow(magnet2[0] - magnet1[0],2) + pow(magnet2[1] - magnet1[1],2) +pow(magnet2[2] - magnet1[2],2));
      //std::cout << "dist mg1-mg2: " << d_mg1_mg2 << "cm" << std::endl;
      double d_mg2_wc4 =  sqrt(pow((*wc4)[0] - magnet2[0],2) + pow((*wc4)[1] - magnet2[1],2) +pow((*wc4)[2] - magnet2[2],2));
      //std::cout << "dist wc4-mg2: " << d_mg2_wc4 << "cm" << std::endl;

      tertiaryLength = d_wc1_mg1 + d_mg1_mg2 + d_mg2_wc4 + 60;
      beamLengthHist->Fill(tertiaryLength);

    }


    
    bool isProton = BS->MassCut(wctrk_momentum[0], tofObject[0], tertiaryLength, UI->tofOffset, ParticleMass, UI->MassCutMin, UI->MassCutMax);

    BeamMassHist->Fill(ParticleMass);
    electronLifetimeHist->Fill(electron_lifetime);

    if(!isMC){
      if (num_wctracks !=1){continue;}
      numWCTrack++;
      if (wctrk_missed > 0){continue;}
      num4PointTrack++;
      if(UI->pickyTracksWC){if (wctrk_picky != 1){continue;}}
      numPickyTrack++;

      BeamMomentumInit->Fill(wctrk_momentum[0]);
      BeamMassInit->Fill(ParticleMass);


      if(UI->qualityTracksWC){if (wctrk_quality !=1){continue;}}
      numQualityTrack++;
      
      BeamMomentumQual->Fill(wctrk_momentum[0]);
      BeamMassQual->Fill(ParticleMass);

      if(isProton){
        std::vector <double> projVector = BS->backProjections(wctrk_XFace[0],wctrk_YFace[0],wctrk_momentum[0],wctrk_theta[0],wctrk_phi[0]);
      
        //beam_x->Fill(projVector[0]);
        //beam_y->Fill(projVector[1]);
        //beam_z->Fill(projVector[2]);

        //beamXtpc0->Fill(wctrk_XFace[0]);
        //beamYtpc0->Fill(wctrk_YFace[0]);
        //beamZtpc0->Fill(0);

        //wctrk4XY->Fill(projVector[0],projVector[1]);
        //wctrkTpcXY->Fill(wctrk_XFace[0],wctrk_YFace[0]);
        //beam_angle_xz->Fill(projVector[3]);
        //beam_angle_yz->Fill(projVector[4]);}

        /////////

        //beamPlotFile->cd()

        
        beam_x = projVector[0];
        beam_y = projVector[1];
        beam_z = projVector[2];
        beam_angle_xz = projVector[3];
        beam_angle_yz = projVector[4];
        beam_momentum = wctrk_momentum[0];
        beam_run = run;
        beam_subrun = subrun;
        beam_event = event;

        beam_tree->Fill();


      }

      //if (verbose){std::cout << "number of ToF objects : " <<num_tof_objects << std::endl;}
      //if (verbose){std::cout << "ToF object size : " <<  << std::endl;}
      if(num_tof_objects != 1){continue;}
      numtofsingle++;

      if(tofObject[0] < UI->tofMin || tofObject[0] > UI->tofMax){continue;}
      numtofvalid++;
      
      BeamMomentum->Fill(wctrk_momentum[0]);
      BeamToF->Fill(tofObject[0]);
      tofMomentHist->Fill(wctrk_momentum[0], tofObject[0]);

      if(applyMassCut){
        if(!isProton){continue;}
      }
      if (isProton){
        BeamMassCutHist->Fill(ParticleMass);
        numMassCut++;
        BeamMomentumMassSel->Fill(wctrk_momentum[0]);
        BeamToFmassSel->Fill(tofObject[0]);
        double mass_const;

        if (UI->isProton){ mass_const = massProton * 1000;}
        else if(UI->isPion){ mass_const = massPion * 1000;}
        else{mass_const = 999.999 * 1000;}



        double calculatedLength =  sqrt((pow(tofObject[0],2) * pow(c_light,2)) / (pow((mass_const/wctrk_momentum[0]),2) +1 ) );
        beamLengthHist->Fill(calculatedLength);
      }

      //if(verbose){std::cout << "matching to tpc" << std::endl;}

      std::vector <double> matchCandidate = BS->BeamMatching(wctrk_x_proj_3cm[0],wctrk_y_proj_3cm[0], wctrk_theta[0], wctrk_phi[0],
                                                             track_xpos, track_ypos, track_zpos, ntracks_reco, ntrack_hits,
                                                             track_length, reco_primary,BSoptions);

      //std::cout << "matching finished" << std::endl;
      numEnteringTracks = BS->EnteringTrkID.size();
      inTracksNumHist->Fill(numEnteringTracks);
      numPileupTracksHist->Fill(BS->PileupTracksBuffer);
      numShowerCutHist->Fill(BS->ShowerTracksBuffer);

      for (int i = 0; i < ntracks_reco ; i++){
        tpcAllTrackZmin1->Fill(BS->AllTrkZmin1[i]);
        tpcAllTrackZmin2->Fill(BS->AllTrkZmin2[i]);
        tpcAllTrackZmax->Fill(BS->AllTrkZmax[i]);
        tpcAllTrackDelX->Fill(BS->AllTrkDelX[i]);
        tpcAllTrackDelY->Fill(BS->AllTrkDelY[i]);
        tpcAllTrackNhits->Fill(BS->AllTrkNhits[i]);
        tpcAllTrackLength->Fill(BS->AllTrkLength[i]);
        tpcAllTrackTheta->Fill(BS->AllTrkTheta[i]);
        allAlphaHist->Fill(BS->AllTrkAlpha[i]);
        allRdistHist->Fill(BS->AllTrkRdist[i]);
        allAlphaHistEnd->Fill(BS->AllTrkAlphaEnd[i]);
        allDirAngleHist->Fill(BS->AllTrkDirCompare[i]);
      }

      NtracksTotalHist->Fill(ntracks_reco);

      for (int  i = 0; i < BS->CircleTrkAlpha.size() ; i++ ){
        circleAlphaHist->Fill(BS->CircleTrkAlpha[i]);
      }

      if(matchCandidate[8] > 1){
        if(UI->MultiMatchEventListSet){
          multiMatchFile.open(UI->MultiMatchEventList, ios::app);
          if(multiMatchFile){
            multiMatchFile << run << "\t"<< subrun << "\t" << event << "\t" << matchCandidate[8] << "\n";
          }
          multiMatchFile.close();
        }
      }


      //if(verbose){std::cout << "ploting delXY" << std::endl;}


      

      if (matchCandidate[0]){
        numMatchedMain++;
        found_primary = true;
        wctrkSelectedXY->Fill(wctrk_x_proj_3cm[0],wctrk_y_proj_3cm[0]);
        BeamMomentumMatch->Fill(wctrk_momentum[0]);
        BeamMassMatch->Fill(ParticleMass);
        if (verbose){std::cout << "Primary Found" << std::endl;}
      }
      //if(found_primary){halo_pileup_number_particles->Fill(BS->PileupTracksBuffer -1);}
      //else{halo_pileup_number_particles->Fill(BS->PileupTracksBuffer);}

      if (verbose){std::cout << "Num Entering Tracks : " << numEnteringTracks << std::endl;}
      if(!found_primary){
        if(verbose){std::cout << "No Valid Primary \n" << std::endl;}
        //continue;
      }

      int pileup_counter = 0;
      
      for (int i = 0;  i < numEnteringTracks; i++){
          int inTrackID = BS->EnteringTrkID[i];
          int inSecond = BS->EnteringTrkSecond[i];
          int inStart =  BS->EnteringTrkStart[i];
          int inEnd = BS->EnteringTrkEnd[i];
          double trackAlpha = BS->EnteringTrkAlpha[i];
          tpcInTracksZ->Fill((*track_zpos)[inTrackID][inStart]);
          tpcInTrackEndZ->Fill((*track_zpos)[inTrackID][inStart]);
          InTrackLength->Fill((*track_length)[inTrackID]);
          alphaHist->Fill(trackAlpha);
          alphaHistEnd->Fill(BS->EnteringTrkAlphaEnd[i]);
          dirAngleHist->Fill(BS->EnteringTrkDirCompare[i]);

          if(inTrackID != reco_primary){

            delBadTrackHist_nocut->Fill(wctrk_x_proj_3cm[0] - (*track_xpos)[inTrackID][inStart],
                    wctrk_y_proj_3cm[0] - (*track_ypos)[inTrackID][inStart]);
            BadTrackLength->Fill((*track_length)[inTrackID]);
            if((*track_length)[inTrackID] > 70){
              delBadTrackHist_len->Fill(wctrk_x_proj_3cm[0] - (*track_xpos)[inTrackID][inStart],
                    wctrk_y_proj_3cm[0] - (*track_ypos)[inTrackID][inStart]);
            }



          }


          if(found_primary){
            if(inTrackID == reco_primary){
            delXYHistMatch->Fill(wctrk_x_proj_3cm[0] - (*track_xpos)[inTrackID][inStart],
            wctrk_y_proj_3cm[0] - (*track_ypos)[inTrackID][inStart]);
            delXYHistPxMatch->Fill(wctrk_x_proj_3cm[0] - (*track_xpos)[inTrackID][inStart]);
            delXYHistPyMatch->Fill(wctrk_y_proj_3cm[0] - (*track_ypos)[inTrackID][inStart]);

              PrimaryLength->Fill((*track_length)[inTrackID]);
              PrimaryStartZ->Fill((*track_zpos)[inTrackID][inStart]);
              zProjPrimaryTrack->Fill((*track_zpos)[inTrackID][inEnd] -  (*track_zpos)[inTrackID][inStart]);
            }
            else{
              delBadTrackHist_m->Fill(wctrk_x_proj_3cm[0] - (*track_xpos)[inTrackID][inStart],
                    wctrk_y_proj_3cm[0] - (*track_ypos)[inTrackID][inStart]);
              BadTrackLength_m->Fill((*track_length)[inTrackID]);
              

              if((*track_length)[inTrackID] > 70){

                BadTrackHist->Fill((*track_xpos)[inTrackID][inStart],(*track_ypos)[inTrackID][inStart]);
                delBadTrackHist_mlen->Fill(wctrk_x_proj_3cm[0] - (*track_xpos)[inTrackID][inStart],
                    wctrk_y_proj_3cm[0] - (*track_ypos)[inTrackID][inStart]);
                
                BadTrackStartZ->Fill((*track_zpos)[inTrackID][inStart]);
                zProjBadTrack->Fill((*track_zpos)[inTrackID][inEnd] -  (*track_zpos)[inTrackID][inStart]);

                double yproj  = (-1 - (*track_zpos)[inTrackID][inStart]) 
                                *( (*track_ypos)[inTrackID][inSecond] - (*track_ypos)[inTrackID][inStart])
                                /( (*track_zpos)[inTrackID][inSecond] - (*track_zpos)[inTrackID][inStart])
                                + (*track_ypos)[inTrackID][inStart];

                double xproj  = (-1 - (*track_zpos)[inTrackID][inStart]) 
                                *( (*track_xpos)[inTrackID][inSecond] - (*track_xpos)[inTrackID][inStart])
                                /( (*track_zpos)[inTrackID][inSecond] - (*track_zpos)[inTrackID][inStart])
                                + (*track_xpos)[inTrackID][inStart];

                //halo_pileup_x->Fill(xproj);
                //halo_pileup_y->Fill(yproj);
                //halo_pileup_z->Fill(-1);

                double angle_yz = atan2((*track_ypos)[inTrackID][3] - (*track_ypos)[inTrackID][inStart],
                                        (*track_zpos)[inTrackID][3] - (*track_zpos)[inTrackID][inStart]);

                double angle_xz = atan2((*track_xpos)[inTrackID][3] - (*track_xpos)[inTrackID][inStart],
                                        (*track_zpos)[inTrackID][3] - (*track_zpos)[inTrackID][inStart]);

                
                //halo_pileup_angle_xz->Fill(angle_xz);
                //halo_pileup_angle_yz->Fill(angle_yz);

                halo_pileup_x.push_back(xproj);
                halo_pileup_y.push_back(yproj);
                halo_pileup_z.push_back(-1.0);
                halo_pileup_angle_xz.push_back(angle_xz);
                halo_pileup_angle_yz.push_back(angle_yz);
                pileup_counter++;
              }
              //halo_pileup_momentum.push_back(trandom_->Landau(1200, 50));

              
              


            }
          }
      }

      if(found_primary){
      halo_pileup_number_particles =  pileup_counter;
      halo_pileup_run = run;
      halo_pileup_subrun = subrun;
      halo_pileup_event = event;
      halo_pileup_tree->Fill();

      halo_pileup_x.clear();
      halo_pileup_y.clear();
      halo_pileup_z.clear();
      halo_pileup_angle_xz.clear();
      halo_pileup_angle_yz.clear();}
      //halo_pileup_momentum.clear();

      
    }

    // ### porting over the work from ProtonAnalyzerMC module -- ryan ###
    // ## grabbing reco primary ##
    //int reco_primary = -1;
    if(isMC){
      double first_reco_z = 99.;
      reco_primary = BS->isTPCPrimary(track_zpos, ntracks_reco, isMC, UI->zBeamCutoff, 
        reco_primary, first_reco_z, verbose);
    }
    if(reco_primary == -1){continue;}//<- skipping events that didn't pass isTPCPrimary 
    else{ found_primary = true;}

    if(UI->skipEventSelection){continue;}//skipping Event Selection - running in beam diagnostic mode



    // ## grabbing interaction point ##
    double temp[6];
    double* candidate_info = ES->findInt(temp, reco_primary, ntracks_reco, 
                                            ntrack_hits, track_xpos, track_ypos, track_zpos,
                                            track_end_x, track_end_y, track_end_z,
                                            col_track_hits, col_track_dedx, col_track_pitch_hit,
                                            col_track_x, col_track_y, col_track_z, ESoptions);

    // ## grabbing what will be histogram entries ##
    //double initial_ke = 99999; //<-- setting to a constant to do dev. needs to be WC info
    double initial_ke = 0; //<-- setting to a constant to do dev. needs to be WC info
    std::vector<double> calo_slab_xpos;
    std::vector<double> calo_slab_ypos;
    std::vector<double> calo_slab_zpos;
    std::vector<double> calo_slab_KE;

    if(isMC) {
      initial_ke = BS->getMCInitialKE(initial_ke, geant_list_size, process_primary, 
                                      NTrTrajPts, MidPosX, MidPosY,  MidPosZ, MidPx, MidPy, MidPz); 
    }
    if(!isMC) {
      initial_ke = BS->getDataInitialKE(initial_ke, wctrk_momentum[0],ParticleMass);
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
      numInteractions++;
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
                                   + pow(int_candidate_z - calo_slab_z, 2) );
        if(calo_slab_z < int_candidate_z){
          if(dist_int_slab < min_dist_int){
            min_dist_int = dist_int_slab;
            calo_int_slab = calo_slab;
          }
        }//<--End if this slab is upstream of int
      }//<--End calo slab loop
    }//<---End if interaction candidate

    double intKE = -1;
    
    // ## incident slabs ## 
    for(unsigned int calo_slab = 1; calo_slab < calo_slab_KE.size(); calo_slab++){
      if(calo_slab > calo_int_slab){continue;}//<--stop after interaction slab 
      hreco_incke->Fill(calo_slab_KE[calo_slab]);
      if(calo_slab == calo_int_slab){
        hreco_intke->Fill(calo_slab_KE[calo_slab]);
        numInteractions++;
        intKE = calo_slab_KE[calo_int_slab];
      }//<-- End if this is the \ slab
    }//<--End calo slab loop

    // ### end of port work -- ryan ###

  double truthIntKE = -1;

  if(UI->SelEventListSet){


    IDfile.open(UI->SelEventList,ios::app);
    if(IDfile){
      IDfile << jentry << "\t" << reco_primary << "\t" << candidate_info[0] << 
      "\t" << candidate_info[1]  << "\t" << candidate_info[2]  << "\t" << candidate_info[3] << "\t " <<
      candidate_info[4] << "\t"<< candidate_info[5] << "\t" << initial_ke << "\t" << intKE <<"\t" << truthIntKE
      << "\t" << ParticleMass << "\n";}
    IDfile.close();
  }

  
  }//end of event Loop

    std::vector<double> BranchDistVect = ES->BranchDistVect;
    std::vector<double> ClusterDistVect = ES->ClusterDistVect;

    for(int i = 0; i < ES->BranchDistVect.size(); i++){
      BranchDistHist->Fill((ES->BranchDistVect)[i]);
    }
    for(int i = 0; i < ClusterDistVect.size(); i++){
      ClusterDistHist->Fill((ES->ClusterDistVect)[i]);
    }





  
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
    //TFile *mc_corrections = new TFile("files/corrections/xsCorrectionBinary.root"); 
    TFile *mc_corrections = new TFile("files/corrections/xsCorrectionBertini.root"); 

    hreco_intke_background = (TH1D*)mc_corrections->Get("hreco_intke_background");
    hreco_incke_background = (TH1D*)mc_corrections->Get("hreco_incke_background");
    hreco_unfolding_matrix_normalized = (TH2D*)mc_corrections->Get("hreco_unfolding_matrix_normalized");
    hreco_intke_eff = (TH1D*)mc_corrections->Get("hreco_int_eff");
    hreco_incke_eff = (TH1D*)mc_corrections->Get("hreco_inc_eff");


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
    
    std::cout << "\n------- Particle ID Results -------\n"<< std::endl;
    std::cout << "Number of entries in file  = " << nentries << std::endl;
    std::cout << "Total events processed: "<< numEventsStart << std::endl;
    std::cout << "Events with only one wc track: "<< numWCTrack << std::endl;
    std::cout << "Events with 4-point wc track: " << num4PointTrack << std::endl;
    if (UI->pickyTracksWC){std::cout << "Events with picky tracks: " << numPickyTrack << std::endl;}
    if (UI->qualityTracksWC){std::cout << "Events passing quality flag: " << numQualityTrack << std::endl;}
    std::cout << "Events with unique ToF value: " << numtofsingle << std::endl;
    std::cout << "Events with ToF value in range (" << UI->tofMin << " , " << UI->tofMax<< ") [ns]: " << numtofvalid << std::endl;
    std::cout << "Events passing mass cut (" << UI->MassCutMin << " , " << UI->MassCutMax << ") [GeV/c^2]: "<< numMassCut << std::endl;
    //std::cout << "Number of Matched Events (According to ProtonXsec): " << numMatchedMain << std::endl;  

    BS->printSummary(BSoptions);
  }

  if (isMC){
    if(verbose){
      std::cout << "\n------- MC Selection Results -------\n"<< std::endl;
      std::cout << "Number of interacting candidates: " << numInteractions << std::endl;
    }
  }


  //if(UI->plotIndividualSet){
    //gtotal_res_dedx_item = TGraph(total_res_v.size(),&total_res_v[0],&total_dedx_v[0]);
    //gtotal_res_dedx_item.SetName("gtotal_res_dedx");
    //gtotal_res_dedx = &gtotal_res_dedx_item;
  //}
  
  
  if(isMC){

    if(UI->rootOutputFileSet){
      outputFile->cd();

      // ## EventSelector histos ##
      BranchDistHist->Write();
      ClusterDistHist->Write();
       
      // ## xs histos ##
      hreco_initialKE->Write();
      hreco_incke->Write();
      hreco_intke->Write();
      hreco_folded_intke_signal->Write(); 
      hreco_folded_incke_signal->Write();
      hreco_unfolded_intke_signal->Write(); 
      hreco_unfolded_incke_signal->Write();
      hreco_xs->Write();
  }


  }





  else if (!(isMC)){
    
    if(UI->rootOutputFileSet){
      if(verbose){std::cout << "Writing to outputFile" << std::endl;}
      outputFile->cd();

      electronLifetimeHist->Write();
    
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
      delXYHistPx->Write();
      delXYHistPy->Write();

      delXYHistCenter->Write();
      delXYHistPxCenter->Write();
      delXYHistPyCenter->Write();

      delXYHistMatch->Write();
      delXYHistPxMatch->Write();
      delXYHistPyMatch->Write();
//      delThetaHist->Write();
//      delPhiHist->Write();
      BadTrackHist->Write();
      BadTrackLength->Write();
      BadTrackLength_m->Write();

      BadTrackStartZ->Write();

      delBadTrackHist_mlen->Write();
      delBadTrackHist_m->Write();
      delBadTrackHist_len->Write();
      delBadTrackHist_nocut->Write();

      BeamMassHist->Write();
      BeamMassInit->Write();
      BeamMassQual->Write();
      BeamMassMatch->Write();

      BeamMomentumInit->Write();
      BeamMomentumQual->Write();
      BeamMomentumMatch->Write();


      BeamMassCutHist->Write();
      BeamToFmassSel->Write();
      BeamMomentumMassSel->Write();

      tofMomentHist->Write();
      numTracksSelHist->Write();
      alphaHist->Write();

      tpcAllTrackZmin1->Write(); 
      tpcAllTrackZmin2->Write(); 
      tpcAllTrackZmax->Write(); 
      tpcAllTrackTheta->Write(); 
      tpcAllTrackNhits->Write();
      tpcAllTrackLength->Write(); 
      tpcAllTrackDelX->Write();  
      tpcAllTrackDelY->Write();
      NtracksTotalHist->Write();
      allRdistHist->Write();
      allAlphaHist->Write();
      allAlphaHistEnd->Write();
      allDirAngleHist->Write();
      dirAngleHist->Write();
      circleAlphaHist->Write();
      alphaHistEnd->Write();
      tpcInTrackEndZ->Write();
      zProjPrimaryTrack->Write();
      zProjBadTrack->Write();
      numPileupTracksHist->Write();
      numShowerCutHist->Write();
      beamLengthHist->Write();
//      tpcPhiHist->Write();
//      wcPhiHist->Write();
//      tpcThetaHist->Write();
//      wcThetaHist->Write();
//      primary_dedx->Write();

      // ## EventSelector histos ##
      BranchDistHist->Write();
      ClusterDistHist->Write();
      beamXtpc0 ->Write();
      beamYtpc0 ->Write();
      beamZtpc0 ->Write();

      wctrk4XY->Write();
      wctrkTpcXY->Write();
      
      //if(UI->plotIndividualSet){gtotal_res_dedx->Write();}
      if(verbose){std::cout << "Writing xs histos" << std::endl;}
      // ## xs histos ##
      hreco_initialKE->Write();
      hreco_incke->Write();
      hreco_intke->Write();

      //beam_x->Write();
      //beam_y->Write();
      //beam_z->Write();
      
      //beam_angle_xz->Write();
      //beam_angle_yz->Write();
      //beam_momentum->Write();

      //halo_pileup_x->Write();
      //halo_pileup_y->Write();
      //halo_pileup_z->Write();
      
      //halo_pileup_angle_xz->Write();
      //halo_pileup_angle_yz->Write();
      //halo_pileup_number_particles->Write();

    }
    if(UI->beamCharFileSet){
      if(verbose){std::cout << "Writing Beam DDMC file" << std::endl;}

      beamPlotFile->cd();
      beam_tree->Write();
      beamPlotFile->Close();

      //beam_x->Write();
      //beam_y->Write();
      //beam_z->Write();
      
      //beam_angle_xz->Write();
      //beam_angle_yz->Write();
      //beam_momentum->Write();
      }

    if(UI->haloCharFileSet){
      if(verbose){std::cout << "Writing Halo DDMC file" << std::endl;}

      haloPlotFile->cd();
      halo_pileup_tree->Write();
      haloPlotFile->Close();

      //halo_pileup_x->Write();
      //halo_pileup_y->Write();
      //halo_pileup_z->Write();
      
      //halo_pileup_angle_xz->Write();
      //halo_pileup_angle_yz->Write();
      //halo_pileup_number_particles->Write();

      }

  }
  if(verbose){std::cout << "Done" << std::endl;}

}




//

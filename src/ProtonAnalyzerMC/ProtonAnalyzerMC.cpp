//#############################################################################
//
// ProtonAnalyzerMC.cpp
//
//  A desription of what this class does
//
// January 17, 2018
//
//#############################################################################

#include "ProtonAnalyzerMC.h"

#include "../Selection/EventSelector.h"
#include "../Selection/BeamSelector.h"



//=============================================================================
// MAIN
//=============================================================================

int main(int argc, char * argv[]) {

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

  TROOT theRoot("theRoot","root for LArP Analysis");
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
  // Create a ProtonAnalyzerMC object and call main analysis function
  //--------------------------------------------------------------------------

  ProtonAnalyzerMC *protonMC = new ProtonAnalyzerMC( jobOptionsFile );
  protonMC->AnalyzeFromNtuples();
  
  return 0;

  
} // end main() 


//=============================================================================
// Constructors
//=============================================================================

ProtonAnalyzerMC::ProtonAnalyzerMC( ) : LArIATAnalysis( ) { 

  tuple = NULL;

}


ProtonAnalyzerMC::ProtonAnalyzerMC( char* jobOptionsFile ) : LArIATAnalysis( jobOptionsFile ) { 

  tuple = NULL;
  

  //== Open list of ntuple files 
  if( UI->inputFilesSet ){
    openNtupleFiles( UI->inputFiles, tuple );
  }
  else{
    cout << endl << "#### ERROR : Input files not properly specified!!!!" << endl << endl;
    exit(0);
  }

  ////== Open output root file and postscript file
  //if( UI->rootOutputFileSet && UI->psOutputFileSet ){
  //  outputFile = new TFile( UI->rootOutputFile, "RECREATE" );

  //  //ps = new TPostScript( UI->psOutputFile, 112 );
  //  //ps->Range(26,18); 
  //  //psPage = 1; 
  //}
  //else{
  //  cout << endl << "#### No output files specified!!!!" << endl << endl;
  //}

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

void ProtonAnalyzerMC::AnalyzeFromNtuples() {


  // globals move later...
  double mass = 938;
  double z = 0.03;
  double z2 = 0.5;

  // histograms, keeping minimal diagnostics
  TH1D *hdedx = new TH1D("hdedx", "dedx", 500, 0, 50);
  TH1D *hintke = new TH1D("hintke", "int ke", 20, 0, 1000);
  TH1D *hincke = new TH1D("hincke", "inc ke", 20, 0, 1000);
  TH1D *h2incke = new TH1D("h2incke", "inc ke", 20, 0, 1000);
  TH1D *hxs    = new TH1D("hxs",    "xs",     20, 0, 1000);

  TH1D *sdedx = new TH1D("sdedx", "slab dedx", 100, 0,   50);
  TH1D *sincke = new TH1D("sincke", "slab inc ke", 20, 0, 1000);
  TH1D *sintke = new TH1D("sintke", "slab int ke", 20, 0, 1000);
  TH1D *sxs    = new TH1D("sxs",    "slab xs",     20, 0, 1000);


  std::cout<<"welp.\n";
  
  EventSelector *ES = new EventSelector();
  BeamSelector *BS = new BeamSelector();
  
  bookNtuple( tuple );
  if (tuple == 0) return;

  Long64_t nentries = tuple->GetEntriesFast();
   
  // ## event loop ##
  for(Long64_t jentry=0; jentry < numEventsToProcess && jentry < nentries; jentry++) {
    
    Long64_t ientry = tuple->LoadTree(jentry);
    if (ientry < 0){continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);
    printEvent();

    // ### Geant4 Information ###
    bool dense_int = false;
    bool g4_interaction = false;
    double intx = 0;
    double inty = 0;
    double intz = 0;
    double int_ke = -1;

    bool g4_int = false;
    double g4_intx = -99;    
    double g4_inty = -99;
    double g4_intz = -99;

    double initial_ke = 0;
    int num_pts_inTPC = 0;
    int first_pt = 0;
    int first_dense_slab = 0;
    double total_dense_dist = 0;
    int num_entries_dense = 0;

    int inslabs = 0;
    std::vector<double> true_slab_xpos;
    std::vector<double> true_slab_ypos;
    std::vector<double> true_slab_zpos;
    std::vector<double> true_slab_ke;


    // ## looping over g4 particles/filling histograms ##
    for(int g4part = 0; g4part < geant_list_size; g4part++) {
      if((*process_primary)[g4part] != 1) {
        continue;
      }// skipping non primary g4 ID

      for(unsigned int nint = 0; nint < (*InteractionPoint).size(); nint++) {
        if( (*InteractionPointType)[nint] == 13) { // 13 is the ID for protonInelastic in this container
          g4_intx = (*MidPosX)[g4part][(*InteractionPoint)[nint]];
          g4_inty = (*MidPosY)[g4part][(*InteractionPoint)[nint]];
          g4_intz = (*MidPosZ)[g4part][(*InteractionPoint)[nint]];
          std::cout<<"\tg4 int: ("<<g4_intx<<", "<<g4_inty<<", "<<g4_intz<<")\n";
        }//< End if interaction is signal
      }//<- End loop over interaction vector
      
      // ## g4 xs thin ##
      double g4_proj_dist = 0;
      double prev_g4_proj_dist = 0;
      int ndense_slab = 1;
      double last_dense_slabx = -99;
      double last_dense_slaby = -99;
      double last_dense_slabz = -99;
      double last_dense_slabke = -99;

      for(int pt = 1; pt < (*NTrTrajPts)[g4part]; pt++){
        double xpos = (*MidPosX)[g4part][pt];
        double ypos = (*MidPosY)[g4part][pt];
        double zpos = (*MidPosZ)[g4part][pt];
        double prev_xpos = (*MidPosX)[g4part][pt-1];
        double prev_ypos = (*MidPosY)[g4part][pt-1];
        double prev_zpos = (*MidPosZ)[g4part][pt-1];
    
        double p = sqrt(pow(1000*(*MidPx)[g4part][pt], 2)
                      + pow(1000*(*MidPy)[g4part][pt], 2)
                      + pow(1000*(*MidPz)[g4part][pt], 2));
        double prev_p = sqrt(pow(1000*(*MidPx)[g4part][pt-1], 2)
                           + pow(1000*(*MidPy)[g4part][pt-1], 2)
                           + pow(1000*(*MidPz)[g4part][pt-1], 2));
        double ke = sqrt(pow(mass, 2) + pow(p, 2)) - mass; 
        double prev_ke = sqrt(pow(mass, 2) + pow(prev_p, 2)) - mass; 
        if( xpos > 0 && xpos < 47.5 && ypos > -20 && ypos < 20 && zpos > 0 && zpos < 90 ) {
          if(first_pt == 0){initial_ke = ke; first_pt++;}
          num_pts_inTPC++;
          double dist_between_points = sqrt( pow(xpos - prev_xpos, 2) 
                                           + pow(ypos - prev_ypos, 2) 
                                           + pow(zpos - prev_zpos, 2));

          double energy_loss = prev_ke - ke;
          double de_dx = energy_loss / dist_between_points;
          total_dense_dist += dist_between_points;
          if(total_dense_dist < 1.){continue;}
          if(!first_dense_slab){
            first_dense_slab = 1;
            //hfirstx->Fill(xpos);
            //hfirsty->Fill(ypos);
            //hfirstz->Fill(zpos);
          }
          // ### Fill Histos ###
          //hDistanceBetweenPoints->Fill(dist_between_points);
          hdedx->Fill(de_dx);
          hincke->Fill(ke);
          h2incke->Fill(ke, dist_between_points/z);
          //hdistvske->Fill(dist_between_points, ke);
          num_entries_dense++; //TotDenseEntries++;
          for(int nint = 0; nint < (*InteractionPoint).size(); nint++) {
            if(pt == (*InteractionPoint)[nint]) {
              //std::cout << "\t\tinteraction point" << std::endl;
              if((*InteractionPointType)[nint] == 13) {
                //nG4Interactions++;
                g4_interaction = true;
                intx = xpos;inty = ypos;intz = zpos;
                int_ke = prev_ke;
                //std::cout << "\t\t\tinelastic! " << pt << std::endl;
                //std::cout<<"\t\t\t\tx,y,z: "<<prev_xpos<<", "<<prev_ypos<<", "<<prev_zpos<<std::endl;
                //std::cout<<"\t\t\t\tx,y,z: "<<xpos<<", "<<ypos<<", "<<zpos<<std::endl;
                dense_int = true;
                hintke->Fill(prev_ke);
              }//<--End if the interaction is inelastic
            }//<--End if the point we're looking at is the interaction point
          }//<--End loop over all interactions (could be 0!)
        }//<--End if in tpc
      }//<--End spt loop



      // ## g4 xs thick ##
      //int inslabs = 0;
      double total_slab_distance = 0;
      //std::cout << "\tslab loop.." << std::endl;
      for(int slab = 0; slab < (*SlabN).size(); slab++){
        double slab_e = 1000*(*SlabE)[slab];
        double slab_x = (*SlabX)[slab];
        double slab_y = (*SlabY)[slab];
        double slab_z = (*SlabZ)[slab];
        double slab_p = sqrt( pow(1000*(*SlapX)[slab], 2) +
                              pow(1000*(*SlapY)[slab], 2) +
                              pow(1000*(*SlapZ)[slab], 2) );
        double slab_ke = sqrt( pow(slab_p, 2) + pow(mass, 2) ) - mass;
        double zero_protection = 0;
        if(slab_x > 0   && slab_x < 47.5 && 
           slab_y > -20 && slab_y < 20   && 
           slab_z > 0   && slab_z < 90){
          inslabs++;
          if(slab !=0){
            double slab_previous_p = sqrt( pow(1000*(*SlapX)[slab-1], 2) +
                                           pow(1000*(*SlapY)[slab-1], 2) +
                                           pow(1000*(*SlapZ)[slab-1], 2) ); 
            double slab_previous_ke = sqrt(pow(slab_previous_p,2)+pow(mass,2)) - mass;
            zero_protection = slab_previous_ke;
            double slab_de = abs(slab_ke - slab_previous_ke);
            double slab_dx = sqrt( pow(slab_x - (*SlabX)[slab-1], 2) +
                                   pow(slab_y - (*SlabY)[slab-1], 2) +
                                   pow(slab_z - (*SlabZ)[slab-1], 2));
            double slab_dedx = slab_de / slab_dx;
            total_slab_distance += slab_dx;
            //sDistanceBetweenSlabs->Fill(slab_dx);
            sdedx->Fill(slab_dedx);
          }
          //else{sfirstx->Fill(slab_x);sfirsty->Fill(slab_y);sfirstz->Fill(slab_z);}
          true_slab_xpos.push_back(slab_x);
          true_slab_ypos.push_back(slab_y);
          true_slab_zpos.push_back(slab_z);
          //sxpos->Fill(slab_x);
          //sypos->Fill(slab_y);
          //szpos->Fill(slab_z);
          // ## getting ke vars ##
          if(slab_ke){
            sincke->Fill(slab_ke);
            true_slab_ke.push_back(slab_ke);
          }
          if(!slab_ke){
            sincke->Fill(zero_protection);
            true_slab_ke.push_back(zero_protection);
          }
          //TotSlabEntries++;
        }//<---End if in tpc
      }//<---End slab loop
      //snslb->Fill(inslabs);
      //if(fmod(total_dense_dist,(int)total_dense_dist) < .05){}
      //if(geant4_print){
      //  std::cout << "dense total distance: " << total_dense_dist << std::endl;
      //  std::cout << "slab total distance: " << total_slab_distance << std::endl;
      //  std::cout << "num slabs: " << inslabs << std::endl;
      //  std::cout << "\tnum entries dense: " << num_entries_dense << std::endl;
      //  std::cout << "\tratio: " << (.5*inslabs)/(z*num_entries_dense) << std::endl;
      //}
      //ratio_entries->Fill((1.*inslabs)/(z*num_entries_dense));
      //ahh->Fill((1.*inslabs)/(z*num_entries_dense), total_dense_dist);
      
      // ## interaction debugging :/ ##
      if(dense_int){
        //std::cout<<"\nintdebug"<<std::endl;
        //std::cout << "number of slabs: " << inslabs << std::endl;
        //std::cout << "interaction point:  " << intx              << "\t\t"
        //                                    << inty              << "\t\t"
        //                                    << intz              << "\t\t" << std::endl;
        //std::cout << "last slab position: " << (*SlabX)[inslabs-1] << "\t\t" 
        //                                    << (*SlabY)[inslabs-1] << "\t\t"
        //                                    << (*SlabZ)[inslabs-1] << "\t\t" << std::endl;

        //double dist_penult = sqrt( pow(intx - (*SlabX)[inslabs-1], 2) +
        //                           pow(inty - (*SlabY)[inslabs-1], 2) +
        //                           pow(intz - (*SlabZ)[inslabs-1], 2) );
        //if(inslabs && dist_penult < 1){
        //  double slab_int_p = sqrt( pow(1000*(*SlapX)[inslabs-1], 2) +
        //                            pow(1000*(*SlapY)[inslabs-1], 2) +
        //                            pow(1000*(*SlapZ)[inslabs-1], 2) );
        //  double slab_int_ke = sqrt( pow(slab_int_p, 2) + pow(mass, 2)) - mass;
        //  sintke->Fill(slab_int_ke);   
        //}//<---End if inelastic interaction 
        if(true_slab_ke[inslabs-1]){
          sintke->Fill(true_slab_ke[inslabs-1]);
        }
        else{
          sintke->Fill(true_slab_ke[inslabs-2]);
        }
        //slab_vs_dense_intke->Fill(true_slab_ke[inslabs-1], int_ke);
        //if(!true_slab_ke[inslabs-1]){
        //  std::cout<<"dense ke: "<<int_ke<<std::endl;
        //  std::cout<<"slab ke: "<<true_slab_ke[inslabs-1]<<std::endl;
        //  std::cout<<"prev slab ke: "<<true_slab_ke[inslabs-2]<<std::endl;
        //  std::cout<<std::endl;
        //}

      }//<--End if this event had an inelastic interaction
    }//<-- End loop over g4 particles



    // ### Reconstructed Information ###

    // ## grabbing reco primary ##
    int reco_primary = -1;
    double first_reco_z = 99.;
    reco_primary = BS->isTPCPrimary(track_zpos, ntracks_reco, isMC, UI->zBeamCutoff, reco_primary, first_reco_z);
    if(reco_primary == -1) {
      continue;
    }//<- skipping events that didn't pass isTPCPrimary 

    // ## grabbing interaction point ##
    bool signal_candidate = ES->findInt(ntracks_reco, reco_primary);
    if(signal_candidate) {
      std::cout << "found inealstic event.\n";
    }//<- End if we found a possible inelastic event

    // ## filling histograms ##

  }//<-- End eventLoop


  // ## write histos ##
  //if(UI->rootOutputFileSet) {
  //  hdedx->Write();
  //  hintke->Write();
  //  hincke->Write();
  //  h2incke->Write();
  //  hxs->Write();
  //  sdedx->Write();
  //  sintke->Write();
  //  sincke->Write();
  //  sxs->Write();
  //}

}// End AnalyzeFromNtuples


// ### EOF ###

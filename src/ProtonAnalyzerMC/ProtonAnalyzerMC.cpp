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

void ProtonAnalyzerMC::AnalyzeFromNtuples() {


  // globals move later...
  double mass = 938.57;
  double z = 0.03;
  double z2 = 0.5;

  double rho = 1.3954;                   // ## g/cm3
  double molar_mass = 39.95;                    // ## g/mol
  double N_A = 6.022 * pow(10, 23);      // ## num/mol
  
  double recip_num_density = molar_mass / (rho * z * N_A); // ## cm2/num
  double sparse_recip_num_density = molar_mass / (rho * z2 * N_A); // ## cm2/num
  double barn = pow(10, -24);



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

  TH1D *hreco_initialKE = new TH1D("hreco_initialKE", "initial ke", 20, 0, 1000);
  TH1D *hreco_intke = new TH1D("hreco_intke", "int ke", 20, 0, 1000);
  TH1D *hreco_intke_signal = new TH1D("hreco_intke_signal", "int ke signal", 20, 0, 1000);
  TH1D *hreco_folded_intke_signal = new TH1D("hreco_folded_intke_signal", "interacting ke (signal)", 20, 0, 1000);
  TH1D *hreco_unfolded_intke_signal = new TH1D("hreco_unfolded_intke_signal", "interacting ke (signal)", 20, 0, 1000);
  TH1D *hreco_intke_background = new TH1D("hreco_intke_background", "int ke background", 20, 0, 1000);

  TH1D *hreco_incke = new TH1D("hreco_incke", "inc ke", 20, 0, 1000);
  TH1D *hreco_incke_signal = new TH1D("hreco_incke_signal", "inc ke signal", 20, 0, 1000);
  TH1D *hreco_folded_incke_signal = new TH1D("hreco_folded_incke_signal", "incident ke (signal)", 20, 0, 1000);
  TH1D *hreco_unfolded_incke_signal = new TH1D("hreco_unfolded_incke_signal", "incident ke (signal)", 20, 0, 1000);
  TH1D *hreco_incke_background = new TH1D("hreco_incke_background", "inc ke background", 20, 0, 1000);


  TH1D *hreco_intke_eff = new TH1D("hreco_int_eff", "interacting selection efficiency", 20, 0, 1000);
  TH1D *hreco_incke_eff = new TH1D("hreco_inc_eff", "incident selection efficiency", 20, 0, 1000);

  TH2D *hreco_unfolding_matrix = new TH2D("hreco_unfolding_matrix", "unfolding", 20, 0, 1000, 20, 0, 1000);
  TH2D *hreco_unfolding_matrix_normalized=new TH2D("hreco_unfolding_matrix_normalized","energy unfolding matrix",20,0,1000,20,0,1000);
  TH1D *hreco_xs = new TH1D("hreco_xs", "p-ar inelastic xs", 20, 0, 1000);


  std::cout<<"welp.\n";
  
  EventSelector *ES = new EventSelector();
  BeamSelector *BS = new BeamSelector();
  
  bookNtuple( tuple, true );
  if (tuple == 0) return;

  Long64_t nentries = tuple->GetEntriesFast();
  TCanvas *c = new TCanvas("c","c",1000, 1000);
   
  // ## event loop ##
  for(Long64_t jentry=0; jentry < numEventsToProcess && jentry < nentries; jentry++) {
    
    Long64_t ientry = tuple->LoadTree(jentry);
    if (ientry < 0){continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);

    std::cout<<std::endl;
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
      //double g4_proj_dist = 0;
      //double prev_g4_proj_dist = 0;
      //int ndense_slab = 1;
      //double last_dense_slabx = -99;
      //double last_dense_slaby = -99;
      //double last_dense_slabz = -99;
      //double last_dense_slabke = -99;

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
        //double slab_e = 1000*(*SlabE)[slab];
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
    reco_primary = BS->isTPCPrimary(track_zpos, ntracks_reco, isMC, UI->zBeamCutoff, reco_primary, first_reco_z, verbose);
    if(reco_primary == -1) {
      continue;
    }//<- skipping events that didn't pass isTPCPrimary 

    // ## grabbing interaction point ##
    double temp[4];
    double* candidate_info = ES->findInt(temp, reco_primary, ntracks_reco, 
                                            ntrack_hits, track_xpos, track_ypos, track_zpos,
                                            track_end_x, track_end_y, track_end_z,
                                            col_track_hits, col_track_dedx, col_track_pitch_hit,
                                            col_track_x, col_track_y, col_track_z);

    // ## comparing int candidate to g4 info ##
    bool signal = false;
    if(candidate_info[0]){ 
      //nRecoCandidates++;
      double reco_x = candidate_info[1];
      double reco_y = candidate_info[2];
      double reco_z = candidate_info[3];
      //std::cout<<"I found a potential interaction here\n";
      //std::cout<<"\t\t("<<reco_x<<", "<<reco_y<<", "<<reco_z<<")\n";
      // mark as signal or background ...
      if(g4_interaction){
        //std::cout<<"\nthere was an inelastic yeah?"<<std::endl;
        //std::cout<<"\t\t("<<intx<<", "<<inty<<", "<<intz<<")\n";

        double dist_reco_g4 = sqrt( pow(reco_x-intx,2)
                                  + pow(reco_y-inty,2)
                                  + pow(reco_z-intz,2) );
        if(dist_reco_g4 < 2){
          //std::cout<<"signal event!\n";
          signal = true;
          //nRecoSignalEvts++;
          // pass a label to histo filling
        }
        else{
          //std::cout<<"background event!\n";
          signal = false;
          //nRecoBckgrdEvts++;
          // pass a label to histo filling
        }
      }//<--End if flag an interaction
      else{
        //std::cout<<"background event!\n";
        signal = false;
        //nRecoBckgrdEvts++;
        // pass a label to histo filling
      }
    }//<--End if flagged an interaction 


    // ## grabbing what will be histogram entries ##
    std::vector<double> calo_slab_xpos;
    std::vector<double> calo_slab_ypos;
    std::vector<double> calo_slab_zpos;
    std::vector<double> calo_slab_KE;
    std::cout<<"initial_ke: "<<initial_ke<<std::endl;
    hreco_initialKE->Fill(initial_ke);
    int slabPass = ES->getSlabInfo(calo_slab_xpos, calo_slab_ypos, calo_slab_zpos, calo_slab_KE,
                                    reco_primary, z2, initial_ke,
                                    col_track_hits, col_track_dedx, col_track_pitch_hit,
                                    col_track_x, col_track_y, col_track_z);
    if(slabPass) {
      std::cout<<"\ttesting multiple vectors D:\n";
      std::cout<<"\t\tz2: "<<z2<<std::endl;
      std::cout<<"\t\tcalo_slab_xpos.size(): "<<calo_slab_xpos.size()<<std::endl;
    }

    // ## comparing to g4 info ##
    std::vector<int> calo_slab_signal(calo_slab_KE.size(), 0);
    int calo_slab_counter = 0;
    for(int calo_slab = 0; calo_slab < calo_slab_KE.size(); calo_slab++){
      double calo_slab_x = calo_slab_xpos[calo_slab];
      double calo_slab_y = calo_slab_ypos[calo_slab];
      double calo_slab_z = calo_slab_zpos[calo_slab];
      //std::cout<<"\tthis slab pos (x,y,z): ("<<calo_slab_xpos[calo_slab]
      //                                 <<", "<<calo_slab_ypos[calo_slab]
      //                                 <<", "<<calo_slab_zpos[calo_slab]<<")\n";
      double min_dist_slab = 99;
      double min_downstream_slab = 99;
      double min_upstream_slab = 99;
      int closest_downstream_slab = -1;
      int closest_upstream_slab = -1;
      for(int true_slab = 0; true_slab < inslabs; true_slab++){
        double true_slab_x = true_slab_xpos[true_slab];
        double true_slab_y = true_slab_ypos[true_slab];
        double true_slab_z = true_slab_zpos[true_slab];
        double dist_calo_true = sqrt( pow(calo_slab_x - true_slab_x, 2) +
                                      pow(calo_slab_y - true_slab_y, 2) +
                                      pow(calo_slab_z - true_slab_z, 2));
        if(dist_calo_true < min_dist_slab){
          min_dist_slab = dist_calo_true;
        }
        // # true slabs downstream calo #
        if(true_slab_z > calo_slab_z){
          if(dist_calo_true < min_downstream_slab){
            min_downstream_slab = dist_calo_true;
            closest_downstream_slab = true_slab;
          }
        }
        // # true slabs upstream calo #
        if(true_slab_z < calo_slab_z){
          if(dist_calo_true < min_upstream_slab){
            min_upstream_slab = dist_calo_true;
            closest_upstream_slab = true_slab;
          }
        }
      }//<--End loop of true slabs
      //std::cout<<"\t\tdistance to closest true slab: "<<min_dist_slab<<std::endl;
      //std::cout<<"\t\tdistance to closest ds true slab: "<<min_downstream_slab<<std::endl;
      //std::cout<<"\t\tdistance to closest us true slab: "<<min_upstream_slab<<std::endl;
      // ## signal background! ##
      if(closest_downstream_slab != -1 && closest_upstream_slab != -1){
        //std::cout<<"\t\t\tthere's a slab on both side!"<<std::endl;
        //std::cout<<"\t\tdistance to closest true slab: "<<min_dist_slab<<std::endl;
        //std::cout<<"\t\tdistance to closest ds true slab: "<<min_downstream_slab<<std::endl;
        //std::cout<<"\t\tdistance to closest us true slab: "<<min_upstream_slab<<std::endl;
        calo_slab_signal[calo_slab] = 1;
        calo_slab_counter++;
      }//<--End if surrounded!
      else{
        // # check if it's the first slab
          // # check distance to closest downstream true slab
        if(calo_slab == 0){
          //std::cout<<"\t\t\tthis is the first slab :0"<<std::endl;
          calo_slab_signal[calo_slab] = 1;
          calo_slab_counter++;
        }
        if(calo_slab_signal[calo_slab-1] == 1){
          //std::cout<<"\t\t\tthis slab is not surrounded but the previous one was\n";
          if(min_upstream_slab < .5){calo_slab_signal[calo_slab] = 1; calo_slab_counter++;}
        }
      }//<--End if not surrounded
    }//<---End loop on calo level slabs


    // ## actually filling histograms ##
    // ## getting which slab will be used for interactions ##
    int calo_int_slab = 999;
    if(candidate_info[0]){
      double int_candidate_x = candidate_info[1];
      double int_candidate_y = candidate_info[2];
      double int_candidate_z = candidate_info[3];
      // loop over slabs to find slab closest to interaction candidate
      double min_dist_int = 99;
      for(int calo_slab = 0; calo_slab < calo_slab_KE.size(); calo_slab++){
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
      //std::cout<<"\n]]]]]]]]]]]]]]]\n";
      //std::cout<<"found possible interaction\n";
      //std::cout<<"\tint slab: "<<calo_int_slab<<std::endl;
      //std::cout<<"\tdist btwn int candidate and closest slab: "<<min_dist_int<<std::endl;
    }//<---End if interaction candidate

    // ## incident slabs ## 
    int ninc_entries = 0;
    for(int calo_slab = 1; calo_slab < calo_slab_KE.size(); calo_slab++){
      if(calo_slab > calo_int_slab){continue;}//<--stop after interaction slab 
      //std::cout<<"\tincident entry: "<<std::endl;
      //std::cout<<"\t\tke: "<<calo_slab_KE[calo_slab]<<std::endl;
      //std::cout<<"\t\tsignal? "<<calo_slab_signal[calo_slab]<<std::endl;
      hreco_incke->Fill(calo_slab_KE[calo_slab]);
      if(calo_slab_signal[calo_slab]){
        hreco_incke_signal->Fill(calo_slab_KE[calo_slab]); 
        hreco_unfolding_matrix->Fill(calo_slab_KE[calo_slab], true_slab_ke[calo_slab-1]);
        ninc_entries++;
      }
      if(!calo_slab_signal[calo_slab]){
        hreco_incke_background->Fill(calo_slab_KE[calo_slab]);
      }
      if(calo_slab == calo_int_slab){
        //std::cout<<"\tinteraction energy: "<<calo_slab_KE[calo_slab]<<std::endl;
        // ###### need to decide what to do with interaction pts far away from slab
        // ###### likely a feature of differences between reco spts and calo objs :/
        // ## also need to check on non-physical entries (negative ?)
        // ### should probably also just grab any non terminating protons as interactions...
        hreco_intke->Fill(calo_slab_KE[calo_slab]);
        if(signal){
          hreco_intke_signal->Fill(calo_slab_KE[calo_slab]); 
        }
        if(!signal){
          hreco_intke_background->Fill(calo_slab_KE[calo_slab]);
        }
      }//<-- End if this is the interacting slab
    }//<--End calo slab loop

  }//<-- End eventLoop

  // ## truth level xs (can go either here or in Xsec module) ##
  for(int iBin = 0; iBin < hintke->GetNbinsX(); iBin++){
   if(h2incke->GetBinContent(iBin) == 0){continue;}
   double ratio   = 1.*hintke->GetBinContent(iBin) / h2incke->GetBinContent(iBin);  
   double temp_xs = ratio * recip_num_density; 
   double xs      = temp_xs / barn;

   double num_err = pow(hintke->GetBinContent(iBin), .5);
   double num = hintke->GetBinContent(iBin);
   if(num == 0){num = 1;}
   double term1 = num_err/num;
   double dem_err = pow(h2incke->GetBinContent(iBin), .5);
   double dem = h2incke->GetBinContent(iBin);
   if(dem == 0){dem =1;}
   double term2 = dem_err/dem;
   double totalError = temp_xs*pow(pow(term1,2) + pow(term2,2), 0.5)/barn;//*recip_num_density*barn;

   //std::cout <<"xs: " << xs << " +- " << totalError << std::endl;
   hxs->SetBinContent(iBin, xs); 
   hxs->SetBinError(iBin,totalError);
  }

  for(int iBin = 0; iBin < sintke->GetNbinsX(); iBin++){
   if(sincke->GetBinContent(iBin) == 0){continue;}
   double ratio   = 1.*sintke->GetBinContent(iBin) / sincke->GetBinContent(iBin);  
   double temp_xs = ratio * sparse_recip_num_density; 
   double xs      = temp_xs / barn;

   double num_err = pow(sintke->GetBinContent(iBin), .5);
   double num = sintke->GetBinContent(iBin);
   if(num == 0){num = 1;}
   double term1 = num_err/num;
   double dem_err = pow(sincke->GetBinContent(iBin), .5);
   double dem = sincke->GetBinContent(iBin);
   if(dem == 0){dem =1;}
   double term2 = dem_err/dem;
   double totalError = temp_xs*pow(pow(term1,2) + pow(term2,2), 0.5)/barn;//*recip_num_density*barn;

   //std::cout <<"xs: " << xs << " +- " << totalError << std::endl;
   sxs->SetBinContent(iBin, xs); 
   sxs->SetBinError(iBin,totalError);
  }

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

  // ## unfolding matrix normalization ##
  for(int iBin = 0; iBin <hreco_incke->GetNbinsX(); iBin++){
    int column_total = 0;
    for(int jBin = 0; jBin < hreco_unfolding_matrix->GetNbinsX(); jBin++){
      column_total += hreco_unfolding_matrix->GetBinContent(iBin, jBin);
    }
    if(!column_total){continue;}
    for(int jBin = 0; jBin < hreco_unfolding_matrix_normalized->GetNbinsX(); jBin++){
      double normalized_value = hreco_unfolding_matrix->GetBinContent(iBin, jBin) / column_total;
      hreco_unfolding_matrix_normalized->SetBinContent(iBin, jBin, normalized_value);
    }
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

  // ## epsilon curves ##
  for(int iBin = 0; iBin < hreco_intke->GetNbinsX(); iBin++){
    if(sintke->GetBinContent(iBin)){
      double int_eff = hreco_unfolded_intke_signal->GetBinContent(iBin) / sintke->GetBinContent(iBin);
      hreco_intke_eff->SetBinContent(iBin, int_eff);
    }
    if(sincke->GetBinContent(iBin)){
      double inc_eff = hreco_unfolded_incke_signal->GetBinContent(iBin) / sincke->GetBinContent(iBin);
      hreco_incke_eff->SetBinContent(iBin, inc_eff);
    }
  }

    // ## putting the pieces together ##
    for(int iBin = 0; iBin < hreco_unfolded_intke_signal->GetNbinsX(); iBin++){
      if(hreco_unfolded_incke_signal->GetBinContent(iBin) == 0){continue;}
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
      double dem_err = pow(sincke->GetBinContent(iBin), .5);
      double term2 = dem_err/dem;
      double totalError = temp_xs*pow(pow(term1,2) + pow(term2,2), 0.5)/barn;//*recip_num_density*barn;

      //std::cout <<"xs: " << xs << " +- " << totalError << std::endl;
      hreco_xs->SetBinContent(iBin, xs); 
      hreco_xs->SetBinError(iBin,totalError);
    }


  // ## write histos ##
  if(UI->rootOutputFileSet) {
    c->cd();
    //hdedx->Draw();
    hdedx->Write();
    hintke->Write();
    hincke->Write();
    h2incke->Write();
    hxs->Write();
    sdedx->Write();
    sintke->Write();
    sincke->Write();
    sxs->Write();
    hreco_initialKE->Write();
    hreco_intke->Write();
    hreco_intke_signal->Write();
    hreco_folded_intke_signal->Write();
    hreco_unfolded_intke_signal->Write();
    hreco_intke_background->Write();
    hreco_incke->Write();
    hreco_incke_signal->Write();
    hreco_folded_incke_signal->Write();
    hreco_unfolded_incke_signal->Write();
    hreco_incke_background->Write();
    hreco_intke_eff->Write();
    hreco_incke_eff->Write();
    hreco_unfolding_matrix->Write();
    hreco_unfolding_matrix_normalized->Write();
    hreco_xs->Write();
  }

}// End AnalyzeFromNtuples


// ### EOF ###

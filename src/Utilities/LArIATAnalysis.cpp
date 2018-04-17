//#############################################################################
//
// LArIATAnalysis.cpp -- Base Class for all Analysis Software 
//
// D. Schmitz  February 8, 2006
//
//#############################################################################


#include "LArIATAnalysis.h"
//#include "UtilityFunctions.h"


const double LArIATAnalysis::pi = 3.14159;                // the beloved constant  
const double LArIATAnalysis::massProton = 0.938;          // proton mass GeV
const double LArIATAnalysis::massPion = 0.140;            // piplus/minus mass GeV
const double LArIATAnalysis::massElectron = 0.000511;     // electron mass GeV
const double LArIATAnalysis::massKaon = 0.494;            // kplus/kminus mass GeV
const double LArIATAnalysis::c_light = 29.9792458;        // cm/ns - speed of light in vacuum


//#############################################################################
//
// Constructors
//
//#############################################################################


LArIATAnalysis::LArIATAnalysis( ){ 

  UI = NULL;
  verbose = 1;

  setKinematicBinning( );

}


LArIATAnalysis::LArIATAnalysis( char* jobOptionsFile ){


  UI = new UserInputs( jobOptionsFile );  
  verbose = UI->verbose;

  setKinematicBinning( );
  
}


//#############################################################################
//
// Additional Functions
//
//#############################################################################


//=============================================================================
// openNtupleFiles() -- opens input root files and returns TChain object
//                      (just DATA or just MC)
//=============================================================================
void LArIATAnalysis::openNtupleFiles( char* file, TChain* &tuple ){ 

  std::cout << std::endl;
  std::cout << "#### Open Files." << std::endl;

  char buffer[500];
  char words[500];

  ifstream *files = new ifstream( file );
  if( files )
  {
    files->getline( buffer, 500 );

    
    int numFiles = atoi( buffer );
    std::cout << "  -- Opening " << numFiles << " ntuples." << std::endl;
    tuple = new TChain( "tuple", "tuple" );
    
    for( int i = 0; i < numFiles; i++ )
    {
      files->getline( buffer, 500 );
      
      sprintf( words, "%s/anatree/anatree", buffer );
      std::cout << words << std::endl;
      tuple->Add( words, 0 );
      std::cout << "        ntuple " << i+1 << " : " << words << std::endl; 
    }
  }
  else
  {
    std::cout << "No list of ntuple files was supplied !! " << std::endl;
  }

}

//===================================================================
// booNtuple -- set the branch addresses for the current ntuple
//===================================================================
void LArIATAnalysis::bookNtuple( TChain* tuple, bool isMC ){

  if( !tuple ){
    std::cout << "You did not supply an ntuple file!!" << std::endl;
    return;
  }
  
  track_primary = 0;
  track_start_x = 0;
  track_start_y = 0;
  track_start_z = 0;
  track_end_x = 0;
  track_end_y = 0;
  track_end_z = 0;
  track_length = 0;
  ntrack_hits = 0;
  track_xpos = 0;
  track_ypos = 0;
  track_zpos = 0;
  ind_track_hits = 0;
  ind_track_ke = 0;
  ind_track_wire = 0;
  ind_track_dedx = 0;
  ind_track_dqdx = 0;
  ind_track_rr = 0;
  ind_track_pitch_hit = 0;
  col_track_hits = 0;
  col_track_ke = 0;
  col_track_x = 0;
  col_track_y = 0;
  col_track_z = 0;
  col_track_wire = 0;
  col_track_dedx = 0;
  col_track_dqdx = 0;
  col_track_rr = 0;
  col_track_pitch_hit = 0;
  PDG = 0;
  StartEnergy = 0;
  StartPx = 0;
  StartPy = 0;
  StartPz = 0;
  EndEnergy = 0;
  EndPx = 0;
  EndPy = 0;
  EndPz = 0;
  StartPointx = 0;
  StartPointy = 0;
  StartPointz = 0;
  EndPointx = 0;
  EndPointy = 0;
  EndPointz = 0;
  Process = 0;
  NumberDaughters = 0;
  Mother = 0;
  TrackId = 0;
  process_primary = 0;
  G4Process = 0;
  G4FinalProcess = 0;
  NTrTrajPts = 0;
  NDTrTrajPts = 0;
  DPdgCode = 0;
  DStartEnergy = 0;
  DStartP = 0;
  MidPosX = 0;
  MidPosY = 0;
  MidPosZ = 0;
  MidPx = 0;
  MidPy = 0;
  MidPz = 0;
  SlabX = 0;
  SlabY = 0;
  SlabZ = 0;
  SlapX = 0;
  SlapY = 0;
  SlapZ = 0;
  SlabN = 0;
  SlabE = 0;
  DMidPosX = 0;
  DMidPosY = 0;
  DMidPosZ = 0;
  InteractionPoint = 0;
  InteractionPointType = 0;
  num_wctracks = 0;


  
  tuple->SetMakeClass(1);
  
  
  tuple->SetBranchAddress("run", &run, &b_run);
  tuple->SetBranchAddress("subrun", &subrun, &b_subrun);
  tuple->SetBranchAddress("event", &event, &b_event);
  tuple->SetBranchAddress("evttime", &evttime, &b_evttime);
  tuple->SetBranchAddress("efield", efield, &b_efield);
  tuple->SetBranchAddress("t0", &t0, &b_t0);
  tuple->SetBranchAddress("ntracks_reco", &ntracks_reco, &b_ntracks_reco);
  tuple->SetBranchAddress("track_primary", &track_primary, &b_track_primary);
  tuple->SetBranchAddress("track_start_x", &track_start_x, &b_track_start_x);
  tuple->SetBranchAddress("track_start_y", &track_start_y, &b_track_start_y);
  tuple->SetBranchAddress("track_start_z", &track_start_z, &b_track_start_z);
  tuple->SetBranchAddress("track_end_x", &track_end_x, &b_track_end_x);
  tuple->SetBranchAddress("track_end_y", &track_end_y, &b_track_end_y);
  tuple->SetBranchAddress("track_end_z", &track_end_z, &b_track_end_z);
  tuple->SetBranchAddress("track_length", &track_length, &b_track_length);
  tuple->SetBranchAddress("ntrack_hits", &ntrack_hits, &b_ntrack_hits);
  tuple->SetBranchAddress("track_xpos", &track_xpos, &b_track_xpos);
  tuple->SetBranchAddress("track_ypos", &track_ypos, &b_track_ypos);
  tuple->SetBranchAddress("track_zpos", &track_zpos, &b_track_zpos);
  tuple->SetBranchAddress("ind_track_hits", &ind_track_hits, &b_ind_track_hits);
  tuple->SetBranchAddress("ind_track_ke", &ind_track_ke, &b_ind_track_ke);
  tuple->SetBranchAddress("ind_track_wire", &ind_track_wire, &b_ind_track_wire);
  tuple->SetBranchAddress("ind_track_dedx", &ind_track_dedx, &b_ind_track_dedx);
  tuple->SetBranchAddress("ind_track_dqdx", &ind_track_dqdx, &b_ind_track_dqdx);
  tuple->SetBranchAddress("ind_track_rr", &ind_track_rr, &b_ind_track_rr);
  tuple->SetBranchAddress("ind_track_pitch_hit", &ind_track_pitch_hit, &b_ind_track_pitch_hit);
  tuple->SetBranchAddress("col_track_hits", &col_track_hits, &b_col_track_hits);
  tuple->SetBranchAddress("col_track_ke", &col_track_ke, &b_col_track_ke);
  tuple->SetBranchAddress("col_track_x", &col_track_x, &b_col_track_x);
  tuple->SetBranchAddress("col_track_y", &col_track_y, &b_col_track_y);
  tuple->SetBranchAddress("col_track_z", &col_track_z, &b_col_track_z);
  tuple->SetBranchAddress("col_track_wire", &col_track_wire, &b_col_track_wire);
  tuple->SetBranchAddress("col_track_dedx", &col_track_dedx, &b_col_track_dedx);
  tuple->SetBranchAddress("col_track_dqdx", &col_track_dqdx, &b_col_track_dqdx);
  tuple->SetBranchAddress("col_track_rr", &col_track_rr, &b_col_track_rr);
  tuple->SetBranchAddress("col_track_pitch_hit", &col_track_pitch_hit, &b_col_track_pitch_hit);
  tuple->SetBranchAddress("no_primaries", &no_primaries, &b_no_primaries);
  tuple->SetBranchAddress("geant_list_size", &geant_list_size, &b_geant_list_size);
  tuple->SetBranchAddress("primary_p", &primary_p, &b_primary);
  tuple->SetBranchAddress("PDG", &PDG, &b_PDG);
  tuple->SetBranchAddress("StartEnergy", &StartEnergy, &b_StartEnergy);
  tuple->SetBranchAddress("StartPx", &StartPx, &b_StartPx);
  tuple->SetBranchAddress("StartPy", &StartPy, &b_StartPy);
  tuple->SetBranchAddress("StartPz", &StartPz, &b_StartPz);
  tuple->SetBranchAddress("EndEnergy", &EndEnergy, &b_EndEnergy);
  tuple->SetBranchAddress("EndPx", &EndPx, &b_EndPx);
  tuple->SetBranchAddress("EndPy", &EndPy, &b_EndPy);
  tuple->SetBranchAddress("EndPz", &EndPz, &b_EndPz);
  tuple->SetBranchAddress("StartPointx", &StartPointx, &b_StartPointx);
  tuple->SetBranchAddress("StartPointy", &StartPointy, &b_StartPointy);
  tuple->SetBranchAddress("StartPointz", &StartPointz, &b_StartPointz);
  tuple->SetBranchAddress("EndPointx", &EndPointx, &b_EndPointx);
  tuple->SetBranchAddress("EndPointy", &EndPointy, &b_EndPointy);
  tuple->SetBranchAddress("EndPointz", &EndPointz, &b_EndPointz);
  tuple->SetBranchAddress("Process", &Process, &b_Process);
  tuple->SetBranchAddress("NumberDaughters", &NumberDaughters, &b_NumberDaughters);
  tuple->SetBranchAddress("Mother", &Mother, &b_Mother);
  tuple->SetBranchAddress("TrackId", &TrackId, &b_TrackId);
  tuple->SetBranchAddress("process_primary", &process_primary, &b_process_primary);
  tuple->SetBranchAddress("G4Process", &G4Process, &b_G4Process);
  tuple->SetBranchAddress("G4FinalProcess", &G4FinalProcess, &b_G4FinalProcess);
  tuple->SetBranchAddress("NTrTrajPts", &NTrTrajPts, &b_NTrTrajPts);
  tuple->SetBranchAddress("NProtonDaughters", &NProtonDaughters, &b_NProtonDaughters);
  tuple->SetBranchAddress("NNeutronDaughters", &NNeutronDaughters, &b_NNeutronDaughters);
  tuple->SetBranchAddress("NDTrTrajPts", &NDTrTrajPts, &b_NDTrTrajPts);
  tuple->SetBranchAddress("DPdgCode", &DPdgCode, &b_DPdgCode);
  tuple->SetBranchAddress("DStartEnergy", &DStartEnergy, &b_DStartEnergy);
  tuple->SetBranchAddress("DStartP", &DStartP, &b_DStartP);
  tuple->SetBranchAddress("MidPosX", &MidPosX, &b_MidPosX);
  tuple->SetBranchAddress("MidPosY", &MidPosY, &b_MidPosY);
  tuple->SetBranchAddress("MidPosZ", &MidPosZ, &b_MidPosZ);
  tuple->SetBranchAddress("MidPx", &MidPx, &b_MidPx);
  tuple->SetBranchAddress("MidPy", &MidPy, &b_MidPy);
  tuple->SetBranchAddress("MidPz", &MidPz, &b_MidPz);
  tuple->SetBranchAddress("SlabX", &SlabX, &b_SlabX);
  tuple->SetBranchAddress("SlabY", &SlabY, &b_SlabY);
  tuple->SetBranchAddress("SlabZ", &SlabZ, &b_SlabZ);
  tuple->SetBranchAddress("SlapX", &SlapX, &b_SlapX);
  tuple->SetBranchAddress("SlapY", &SlapY, &b_SlapY);
  tuple->SetBranchAddress("SlapZ", &SlapZ, &b_SlapZ);
  tuple->SetBranchAddress("SlabN", &SlabN, &b_SlabN);
  tuple->SetBranchAddress("SlabE", &SlabE, &b_SlabE);
  tuple->SetBranchAddress("LastSlabInt", &LastSlabInt, &b_LastSlabInt);
  tuple->SetBranchAddress("DMidPosX", &DMidPosX, &b_DMidPosX);
  tuple->SetBranchAddress("DMidPosY", &DMidPosY, &b_DMidPosY);
  tuple->SetBranchAddress("DMidPosZ", &DMidPosZ, &b_DMidPosZ);
  tuple->SetBranchAddress("InteractionPoint", &InteractionPoint, &b_InteractionPoint);
  tuple->SetBranchAddress("InteractionPointType", &InteractionPointType, &b_InteractionPointType);
  if(!isMC){
    tuple->SetBranchAddress("wctrk_momentum", &wctrk_momentum, &b_wctrk_momentum);
    tuple->SetBranchAddress("wctrk_XFace", &wctrk_XFace, &b_wctrk_XFace);
    tuple->SetBranchAddress("wctrk_YFace", &wctrk_YFace, &b_wctrk_YFace);
    tuple->SetBranchAddress("wctrk_theta", &wctrk_theta, &b_wctrk_theta);
    tuple->SetBranchAddress("wctrk_phi", &wctrk_phi, &b_wctrk_phi);
    tuple->SetBranchAddress("num_wctracks", &num_wctracks, &b_num_wctracks);
    tuple->SetBranchAddress("tofObject", &tofObject, &b_tofObject); }

  
}


//=============================================================================
// setKinematicBinning( ) -- sets p, th, phi, thx, thy binning for analysis
//=============================================================================
void LArIATAnalysis::setKinematicBinning( ){

  /*
  setpBins( p_bins, p_vector );
  setThBins( th_bins, th_vector );
  setPhiBins( phi_bins, phi_vector );
  setThxBins( thx_bins, thx_vector );
  setThyBins( thy_bins, thy_vector );

  p_min = p_vector[0];
  p_max = p_vector[p_bins];
  th_min = th_vector[0];
  th_max = th_vector[th_bins];
  phi_min = phi_vector[0];
  phi_max = phi_vector[th_bins];
  thx_min = thx_vector[0];
  thx_max = thx_vector[thx_bins];
  thy_min = thy_vector[0];
  thy_max = thy_vector[thy_bins];
  */
}

//=============================================================================
// printEvent( )
//=============================================================================

void LArIATAnalysis::printEvent(){

  


  
  cout << "run = " << run << ", subrun = " << subrun << ", event = " << event << ", Ntracks = " << ntracks_reco << ", wctrk_momentum = " << wctrk_momentum[0] << endl;

  

}


//=============================================================================
// finalize() -- use to close output files, etc
//=============================================================================
void LArIATAnalysis::finalize( ){

  cout << "#### Finalize." << endl;

}


 
//############################################################################
//
// END LArIATAnalysis.cpp
//
//############################################################################

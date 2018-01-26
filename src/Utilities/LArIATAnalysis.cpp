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

  //UI = NULL;
  //verbose = 1;

  setKinematicBinning( );

}


/*
LArIATAnalysis::LArIATAnalysis( char* jobOptionsFile ){


  UI = new UserInputs( jobOptionsFile );  
  verbose = UI->verbose;

  setKinematicBinning( );
  
}
*/

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
    tuple = new TChain("tuple","tuple");
    
    for( int i = 0; i < numFiles; i++ )
    {
      files->getline( buffer, 500 );
      sprintf( words, "%s/analysis", buffer ); 
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
void LArIATAnalysis::bookNtuple( TChain* tuple ){

  if( !tuple ){
    std::cout << "You did not supply an ntuple file!!" << std::endl;
    return;
  }

  tuple->SetBranchAddress("run", &run);
  tuple->SetBranchAddress("subrun", &subrun);
  tuple->SetBranchAddress("event", &event);
  tuple->SetBranchAddress("evttime", &evttime);
  tuple->SetBranchAddress("efield", efield);
  tuple->SetBranchAddress("t0", &t0);
  tuple->SetBranchAddress("ntracks_reco", &ntracks_reco);
  tuple->SetBranchAddress("track_start_x", &track_start_x);
  tuple->SetBranchAddress("track_start_y", &track_start_y);
  tuple->SetBranchAddress("track_start_z", &track_start_z);
  tuple->SetBranchAddress("track_end_x", &track_end_x);
  tuple->SetBranchAddress("track_end_y", &track_end_y);
  tuple->SetBranchAddress("track_end_z", &track_end_z);
  tuple->SetBranchAddress("track_length", &track_length);
  tuple->SetBranchAddress("ntrack_hits", &ntrack_hits);
  tuple->SetBranchAddress("track_xpos", &track_xpos);
  tuple->SetBranchAddress("track_ypos", &track_ypos);
  tuple->SetBranchAddress("track_zpos", &track_zpos);
  tuple->SetBranchAddress("ind_track_hits", &ind_track_hits);
  tuple->SetBranchAddress("ind_track_ke", &ind_track_ke);
  tuple->SetBranchAddress("ind_track_wire", &ind_track_wire);
  tuple->SetBranchAddress("ind_track_dedx", &ind_track_dedx);
  tuple->SetBranchAddress("ind_track_dqdx", &ind_track_dqdx);
  tuple->SetBranchAddress("ind_track_rr", &ind_track_rr);
  tuple->SetBranchAddress("ind_track_pitch_hit", &ind_track_pitch_hit);
  tuple->SetBranchAddress("col_track_hits", &col_track_hits);
  tuple->SetBranchAddress("col_track_ke", &col_track_ke);
  tuple->SetBranchAddress("col_track_wire", &col_track_wire);
  tuple->SetBranchAddress("col_track_dedx", &col_track_dedx);
  tuple->SetBranchAddress("col_track_dqdx", &col_track_dqdx);
  tuple->SetBranchAddress("col_track_rr", &col_track_rr);
  tuple->SetBranchAddress("col_track_pitch_hit", &col_track_pitch_hit);
  tuple->SetBranchAddress("trkpitch", trkpitch);
  tuple->SetBranchAddress("num_hit", &num_hit);
  tuple->SetBranchAddress("hit_channel", &hit_channel);
  tuple->SetBranchAddress("hit_integral", &hit_integral);
  tuple->SetBranchAddress("no_primaries", &no_primaries);
  tuple->SetBranchAddress("geant_list_size", &geant_list_size);
  tuple->SetBranchAddress("primary_p", &primary_p);
  tuple->SetBranchAddress("PDG", &PDG);
  tuple->SetBranchAddress("StartEnergy", &StartEnergy);
  tuple->SetBranchAddress("StartPx", &StartPx);
  tuple->SetBranchAddress("StartPy", &StartPy);
  tuple->SetBranchAddress("StartPz", &StartPz);
  tuple->SetBranchAddress("EndEnergy", &EndEnergy);
  tuple->SetBranchAddress("EndPx", &EndPx);
  tuple->SetBranchAddress("EndPy", &EndPy);
  tuple->SetBranchAddress("EndPz", &EndPz);
  tuple->SetBranchAddress("StartPointx", &StartPointx);
  tuple->SetBranchAddress("StartPointy", &StartPointy);
  tuple->SetBranchAddress("StartPointz", &StartPointz);
  tuple->SetBranchAddress("EndPointx", &EndPointx);
  tuple->SetBranchAddress("EndPointy", &EndPointy);
  tuple->SetBranchAddress("EndPointz", &EndPointz);
  tuple->SetBranchAddress("Process", &Process);
  tuple->SetBranchAddress("NumberDaughters", &NumberDaughters);
  tuple->SetBranchAddress("primary_simChannel_num_voxel", &primary_simChannel_num_voxel);
  tuple->SetBranchAddress("primary_simChannel_voxel_dr", &primary_simChannel_voxel_dr);
  tuple->SetBranchAddress("primary_simChannel_voxel_E", &primary_simChannel_voxel_E);
  tuple->SetBranchAddress("primary_simChannel_voxel_e", &primary_simChannel_voxel_e);
  tuple->SetBranchAddress("primary_simChannel_voxel_x", &primary_simChannel_voxel_x);
  tuple->SetBranchAddress("primary_simChannel_voxel_y", &primary_simChannel_voxel_y);
  tuple->SetBranchAddress("primary_simChannel_voxel_z", &primary_simChannel_voxel_z);
  tuple->SetBranchAddress("primary_num_simChannel", &primary_num_simChannel);
  tuple->SetBranchAddress("primary_simChannel", &primary_simChannel);
  tuple->SetBranchAddress("primary_simChannel_dr", &primary_simChannel_dr);
  tuple->SetBranchAddress("primary_simChannel_E", &primary_simChannel_E);
  tuple->SetBranchAddress("primary_simChannel_e", &primary_simChannel_e);
  tuple->SetBranchAddress("Mother", &Mother);
  tuple->SetBranchAddress("TrackId", &TrackId);
  tuple->SetBranchAddress("process_primary", &process_primary);
  tuple->SetBranchAddress("G4Process", &G4Process);
  tuple->SetBranchAddress("G4FinalProcess", &G4FinalProcess);
  tuple->SetBranchAddress("NTrTrajPts", &NTrTrajPts);
  tuple->SetBranchAddress("NProtonDaughters", &NProtonDaughters);
  tuple->SetBranchAddress("NNeutronDaughters", &NNeutronDaughters);
  tuple->SetBranchAddress("NDTrTrajPts", &NDTrTrajPts);
  tuple->SetBranchAddress("DPdgCode", &DPdgCode);
  tuple->SetBranchAddress("DStartEnergy", &DStartEnergy);
  tuple->SetBranchAddress("MidPosX", &MidPosX);
  tuple->SetBranchAddress("MidPosY", &MidPosY);
  tuple->SetBranchAddress("MidPosZ", &MidPosZ);
  tuple->SetBranchAddress("MidPx", &MidPx);
  tuple->SetBranchAddress("MidPy", &MidPy);
  tuple->SetBranchAddress("MidPz", &MidPz);
  tuple->SetBranchAddress("DPMidPosX", &DPMidPosX);
  tuple->SetBranchAddress("DPMidPosY", &DPMidPosY);
  tuple->SetBranchAddress("DPMidPosZ", &DPMidPosZ);
  tuple->SetBranchAddress("InteractionPoint", &InteractionPoint);
  tuple->SetBranchAddress("InteractionPointType", &InteractionPointType);
  
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
void LArIATAnalysis::printEvent( ){

  cout << endl;
  
  cout << "run = " << run << "subrun = " << subrun << ", event = " << event << ", evt time = " << evttime << ", t0 = " << t0 << endl;

  cout << endl;

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

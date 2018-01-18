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
void openNtupleFiles( char* file, TChain* &tuple ){ 

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

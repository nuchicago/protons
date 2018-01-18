
#ifndef LARIATANALYSIS_H
#define LARIATANALYSIS_H


#include "ROOTinclude.h"
//#include "UserInputs.h"


class LArIATAnalysis {

 public:

  //--------------------------------------------
  // Types of Constructors
  //--------------------------------------------

  LArIATAnalysis( );
  //LArIATAnalysis( char* jobOptionsFile );

  void printEvent( );

 protected:

  void openNtupleFiles( char* file, TChain* &tuple ); 
  void bookNtuple( TChain* tuple ); 
  
  void setKinematicBinning( );

  void finalize( );
 

  //---------------------------------------------
  // UserInput manager
  //---------------------------------------------

  //UserInputs *UI;
  //int verbose;


  //---------------------------------------------
  // utility items
  //---------------------------------------------

  TLatex *label;
  char words[200];
  char morewords[200]; 


  //---------------------------------------------
  // useful constants
  //---------------------------------------------

  static const double pi;                  // the beloved constant  
  static const double massProton;          // proton mass GeV
  static const double massPion;            // piplus/minus mass GeV
  static const double massElectron;        // electron mass GeV
  static const double massKaon;            // kplus/kminus mass GeV
  static const double c_light;             // cm/ns - speed of light in vacuum


  //---------------------------------------------
  // Analysis kinematic binning
  //---------------------------------------------
  
  int p_bins;
  double p_min, p_max;
  double* p_vector;

  int th_bins;
  double th_min, th_max;
  double* th_vector;

  int phi_bins;
  double phi_min, phi_max;
  double* phi_vector;

  int thx_bins;
  double thx_min, thx_max;
  double* thx_vector;

  int thy_bins;
  double thy_min, thy_max;
  double* thy_vector;

  //---------------------------------------------
  // NTuple Variables
  //---------------------------------------------

  Int_t           run;
  Int_t           subrun;
  Int_t           event;
  Double_t        evttime;
  Double_t        efield[3];
  Int_t           t0;
  Int_t           ntracks_reco;

  vector<double>  *track_start_x;
  vector<double>  *track_start_y;
  vector<double>  *track_start_z;
  vector<double>  *track_end_x;
  vector<double>  *track_end_y;
  vector<double>  *track_end_z;
  vector<double>  *track_length;
  vector<int>     *ntrack_hits;

  vector<vector<double> > *track_xpos;
  vector<vector<double> > *track_ypos;
  vector<vector<double> > *track_zpos;

  vector<int>     *ind_track_hits;
  vector<double>  *ind_track_ke;
  vector<vector<double> > *ind_track_wire;
  vector<vector<double> > *ind_track_dedx;
  vector<vector<double> > *ind_track_dqdx;
  vector<vector<double> > *ind_track_rr;
  vector<vector<double> > *ind_track_pitch_hit;
 
  vector<int>     *col_track_hits;
  vector<double>  *col_track_ke;
  vector<vector<double> > *col_track_wire;
  vector<vector<double> > *col_track_dedx;
  vector<vector<double> > *col_track_dqdx;
  vector<vector<double> > *col_track_rr;
  vector<vector<double> > *col_track_pitch_hit;
  
  Double_t        trkpitch[14][2];   //[ntracks_reco]
  Int_t           num_hit;
  vector<double>  *hit_channel;
  vector<double>  *hit_integral;
  Int_t           no_primaries;
  Int_t           geant_list_size;
  Double_t        primary_p;
  vector<int>     *PDG;
  vector<double>  *StartEnergy;
  vector<double>  *StartPx;
  vector<double>  *StartPy;
  vector<double>  *StartPz;
  vector<double>  *EndEnergy;
  vector<double>  *EndPx;
  vector<double>  *EndPy;
  vector<double>  *EndPz;
  vector<double>  *StartPointx;
  vector<double>  *StartPointy;
  vector<double>  *StartPointz;
  vector<double>  *EndPointx;
  vector<double>  *EndPointy;
  vector<double>  *EndPointz;
  vector<int>     *Process;
  vector<int>     *NumberDaughters;
  vector<int>     *primary_simChannel_num_voxel;
  
  vector<vector<double> > *primary_simChannel_voxel_dr;
  vector<vector<double> > *primary_simChannel_voxel_E;
  vector<vector<double> > *primary_simChannel_voxel_e;
  vector<vector<double> > *primary_simChannel_voxel_x;
  vector<vector<double> > *primary_simChannel_voxel_y;
  vector<vector<double> > *primary_simChannel_voxel_z;
  
  Int_t           primary_num_simChannel;
  vector<double>  *primary_simChannel;
  vector<double>  *primary_simChannel_dr;
  vector<double>  *primary_simChannel_E;
  vector<double>  *primary_simChannel_e;
  vector<int>     *Mother;
  vector<int>     *TrackId;
  vector<int>     *process_primary;
  vector<string>  *G4Process;
  vector<string>  *G4FinalProcess;
  vector<int>     *NTrTrajPts;
  Int_t           NProtonDaughters;
  Int_t           NNeutronDaughters;
  vector<int>     *NDTrTrajPts;
  vector<int>     *DPdgCode;
  vector<int>     *DStartEnergy;
  
  vector<vector<double> > *MidPosX;
  vector<vector<double> > *MidPosY;
  vector<vector<double> > *MidPosZ;
  vector<vector<double> > *MidPx;
  vector<vector<double> > *MidPy;
  vector<vector<double> > *MidPz;
  vector<vector<double> > *DPMidPosX;
  vector<vector<double> > *DPMidPosY;
  vector<vector<double> > *DPMidPosZ;
  
  vector<int>     *InteractionPoint;
  vector<int>     *InteractionPointType;
  
};

#endif

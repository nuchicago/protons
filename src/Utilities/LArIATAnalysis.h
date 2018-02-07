
#ifndef LARIATANALYSIS_H
#define LARIATANALYSIS_H

#include <vector>
#include <iostream>
#include "ROOTinclude.h"
#include "UserInputs.h"




class LArIATAnalysis {

 public:

  //--------------------------------------------
  // Types of Constructors
  //--------------------------------------------

  LArIATAnalysis( );
  LArIATAnalysis( char* jobOptionsFile );

  void printEvent(TChain* tuple, Long64_t entry);

 protected:

  void openNtupleFiles( char* file, TChain* &tuple ); 
  void bookNtuple( TChain* tuple ); 
  
  void setKinematicBinning( );

  void finalize( );
 

  //---------------------------------------------
  // UserInput manager
  //---------------------------------------------

  UserInputs *UI;
  int verbose;
  bool isMC;


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
  std::vector<int>     *track_primary;
  std::vector<double>  *track_start_x;
  std::vector<double>  *track_start_y;
  std::vector<double>  *track_start_z;
  std::vector<double>  *track_end_x;
  std::vector<double>  *track_end_y;
  std::vector<double>  *track_end_z;
  std::vector<double>  *track_length;
  std::vector<int>     *ntrack_hits;
  std::vector<std::vector<double> > *track_xpos;
  std::vector<std::vector<double> > *track_ypos;
  std::vector<std::vector<double> > *track_zpos;
  std::vector<int>     *ind_track_hits;
  std::vector<double>  *ind_track_ke;
  std::vector<std::vector<double> > *ind_track_wire;
  std::vector<std::vector<double> > *ind_track_dedx;
  std::vector<std::vector<double> > *ind_track_dqdx;
  std::vector<std::vector<double> > *ind_track_rr;
  std::vector<std::vector<double> > *ind_track_pitch_hit;
  std::vector<int>     *col_track_hits;
  std::vector<double>  *col_track_ke;
  std::vector<std::vector<double> > *col_track_x;
  std::vector<std::vector<double> > *col_track_y;
  std::vector<std::vector<double> > *col_track_z;
  std::vector<std::vector<double> > *col_track_wire;
  std::vector<std::vector<double> > *col_track_dedx;
  std::vector<std::vector<double> > *col_track_dqdx;
  std::vector<std::vector<double> > *col_track_rr;
  std::vector<std::vector<double> > *col_track_pitch_hit;
  Int_t           no_primaries;
  Int_t           geant_list_size;
  Double_t        primary_p;
  std::vector<int>     *PDG;
  std::vector<double>  *StartEnergy;
  std::vector<double>  *StartPx;
  std::vector<double>  *StartPy;
  std::vector<double>  *StartPz;
  std::vector<double>  *EndEnergy;
  std::vector<double>  *EndPx;
  std::vector<double>  *EndPy;
  std::vector<double>  *EndPz;
  std::vector<double>  *StartPointx;
  std::vector<double>  *StartPointy;
  std::vector<double>  *StartPointz;
  std::vector<double>  *EndPointx;
  std::vector<double>  *EndPointy;
  std::vector<double>  *EndPointz;
  std::vector<int>     *Process;
  std::vector<int>     *NumberDaughters;
  std::vector<int>     *Mother;
  std::vector<int>     *TrackId;
  std::vector<int>     *process_primary;
  std::vector<string>  *G4Process;
  std::vector<string>  *G4FinalProcess;
  std::vector<int>     *NTrTrajPts;
  Int_t           NProtonDaughters;
  Int_t           NNeutronDaughters;
  std::vector<double>  *NDTrTrajPts;
  std::vector<int>     *DPdgCode;
  std::vector<int>     *DStartEnergy;
  std::vector<std::vector<double> > *MidPosX;
  std::vector<std::vector<double> > *MidPosY;
  std::vector<std::vector<double> > *MidPosZ;
  std::vector<std::vector<double> > *MidPx;
  std::vector<std::vector<double> > *MidPy;
  std::vector<std::vector<double> > *MidPz;
  std::vector<double>  *SlabX;
  std::vector<double>  *SlabY;
  std::vector<double>  *SlabZ;
  std::vector<double>  *SlapX;
  std::vector<double>  *SlapY;
  std::vector<double>  *SlapZ;
  std::vector<int>     *SlabN;
  std::vector<double>  *SlabE;
  Int_t           LastSlabInt;
  std::vector<std::vector<double> > *DPMidPosX;
  std::vector<std::vector<double> > *DPMidPosY;
  std::vector<std::vector<double> > *DPMidPosZ;
  std::vector<int>     *InteractionPoint;
  std::vector<int>     *InteractionPointType;
  
  // List of branches
  TBranch        *b_run;   //!
  TBranch        *b_subrun;   //!
  TBranch        *b_event;   //!
  TBranch        *b_evttime;   //!
  TBranch        *b_efield;   //!
  TBranch        *b_t0;   //!
  TBranch        *b_ntracks_reco;   //!
  TBranch        *b_track_primary;   //!
  TBranch        *b_track_start_x;   //!
  TBranch        *b_track_start_y;   //!
  TBranch        *b_track_start_z;   //!
  TBranch        *b_track_end_x;   //!
  TBranch        *b_track_end_y;   //!
  TBranch        *b_track_end_z;   //!
  TBranch        *b_track_length;   //!
  TBranch        *b_ntrack_hits;   //!
  TBranch        *b_track_xpos;   //!
  TBranch        *b_track_ypos;   //!
  TBranch        *b_track_zpos;   //!
  TBranch        *b_ind_track_hits;   //!
  TBranch        *b_ind_track_ke;   //!
  TBranch        *b_ind_track_wire;   //!
  TBranch        *b_ind_track_dedx;   //!
  TBranch        *b_ind_track_dqdx;   //!
  TBranch        *b_ind_track_rr;   //!
  TBranch        *b_ind_track_pitch_hit;   //!
  TBranch        *b_col_track_hits;   //!
  TBranch        *b_col_track_ke;   //!
  TBranch        *b_col_track_x;   //!
  TBranch        *b_col_track_y;   //!
  TBranch        *b_col_track_z;   //!
  TBranch        *b_col_track_wire;   //!
  TBranch        *b_col_track_dedx;   //!
  TBranch        *b_col_track_dqdx;   //!
  TBranch        *b_col_track_rr;   //!
  TBranch        *b_col_track_pitch_hit;   //!
  TBranch        *b_no_primaries;   //!
  TBranch        *b_geant_list_size;   //!
  TBranch        *b_primary;   //!
  TBranch        *b_PDG;   //!
  TBranch        *b_StartEnergy;   //!
  TBranch        *b_StartPx;   //!
  TBranch        *b_StartPy;   //!
  TBranch        *b_StartPz;   //!
  TBranch        *b_EndEnergy;   //!
  TBranch        *b_EndPx;   //!
  TBranch        *b_EndPy;   //!
  TBranch        *b_EndPz;   //!
  TBranch        *b_StartPointx;   //!
  TBranch        *b_StartPointy;   //!
  TBranch        *b_StartPointz;   //!
  TBranch        *b_EndPointx;   //!
  TBranch        *b_EndPointy;   //!
  TBranch        *b_EndPointz;   //!
  TBranch        *b_Process;   //!
  TBranch        *b_NumberDaughters;   //!
  TBranch        *b_Mother;   //!
  TBranch        *b_TrackId;   //!
  TBranch        *b_process_primary;   //!
  TBranch        *b_G4Process;   //!
  TBranch        *b_G4FinalProcess;   //!
  TBranch        *b_NTrTrajPts;   //!
  TBranch        *b_NProtonDaughters;   //!
  TBranch        *b_NNeutronDaughters;   //!
  TBranch        *b_NDTrTrajPts;   //!
  TBranch        *b_DPdgCode;   //!
  TBranch        *b_DStartEnergy;   //!
  TBranch        *b_MidPosX;   //!
  TBranch        *b_MidPosY;   //!
  TBranch        *b_MidPosZ;   //!
  TBranch        *b_MidPx;   //!
  TBranch        *b_MidPy;   //!
  TBranch        *b_MidPz;   //!
  TBranch        *b_SlabX;   //!
  TBranch        *b_SlabY;   //!
  TBranch        *b_SlabZ;   //!
  TBranch        *b_SlapX;   //!
  TBranch        *b_SlapY;   //!
  TBranch        *b_SlapZ;   //!
  TBranch        *b_SlabN;   //!
  TBranch        *b_SlabE;   //!
  TBranch        *b_LastSlabInt;   //!
  TBranch        *b_DPMidPosX;   //!
  TBranch        *b_DPMidPosY;   //!
  TBranch        *b_DPMidPosZ;   //!
  TBranch        *b_InteractionPoint;   //!
  TBranch        *b_InteractionPointType;   //!
  
};

#endif

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 17 14:54:30 2017 by ROOT version 6.06/08
// from TTree anatree/analysis tree
// found on file: full_anaTree.root
//////////////////////////////////////////////////////////

#ifndef blah_h
#define blah_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class blah {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_evttime;   //!
   TBranch        *b_efield;   //!
   TBranch        *b_t0;   //!
   TBranch        *b_ntracks_reco;   //!
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
   TBranch        *b_col_track_wire;   //!
   TBranch        *b_col_track_dedx;   //!
   TBranch        *b_col_track_dqdx;   //!
   TBranch        *b_col_track_rr;   //!
   TBranch        *b_col_track_pitch_hit;   //!
   TBranch        *b_trkpitch;   //!
   TBranch        *b_num_hits;   //!
   TBranch        *b_hit_channel;   //!
   TBranch        *b_hit_integral;   //!
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
   TBranch        *b_primary_simChannel_num_voxel;   //!
   TBranch        *b_primary_simChannel_voxel_dr;   //!
   TBranch        *b_primary_simChannel_voxel_E;   //!
   TBranch        *b_primary_simChannel_voxel_e;   //!
   TBranch        *b_primary_simChannel_voxel_x;   //!
   TBranch        *b_primary_simChannel_voxel_y;   //!
   TBranch        *b_primary_simChannel_voxel_z;   //!
   TBranch        *b_primary_num_simChannel;   //!
   TBranch        *b_primary_simChannel;   //!
   TBranch        *b_primary_simChannel_dr;   //!
   TBranch        *b_primary_simChannel_E;   //!
   TBranch        *b_primary_simChannel_e;   //!
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
   TBranch        *b_DPMidPosX;   //!
   TBranch        *b_DPMidPosY;   //!
   TBranch        *b_DPMidPosZ;   //!
   TBranch        *b_InteractionPoint;   //!
   TBranch        *b_InteractionPointType;   //!

   blah(TTree *tree=0);
   virtual ~blah();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   // ### adding functions here ###
   virtual double get_angle(double a[], double b[]); 
   virtual double distance(double a[], double b[]);
   virtual bool   is_in_tpc(int spc_pt);
   virtual int    get_primary(double mc_tpc_x, double mc_tpc_y, double mc_vx, double mc_vy, double mc_vz);
   virtual int    get_primary_track_id();
   virtual int    get_primary_process();
};

#endif

#ifdef blah_cxx
blah::blah(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("full_anaTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("full_anaTree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("full_anaTree.root:/anatree");
      dir->GetObject("anatree",tree);

   }
   Init(tree);
}

blah::~blah()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t blah::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t blah::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void blah::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
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
   col_track_wire = 0;
   col_track_dedx = 0;
   col_track_dqdx = 0;
   col_track_rr = 0;
   col_track_pitch_hit = 0;
   hit_channel = 0;
   hit_integral = 0;
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
   primary_simChannel_num_voxel = 0;
   primary_simChannel_voxel_dr = 0;
   primary_simChannel_voxel_E = 0;
   primary_simChannel_voxel_e = 0;
   primary_simChannel_voxel_x = 0;
   primary_simChannel_voxel_y = 0;
   primary_simChannel_voxel_z = 0;
   primary_simChannel = 0;
   primary_simChannel_dr = 0;
   primary_simChannel_E = 0;
   primary_simChannel_e = 0;
   Mother = 0;
   TrackId = 0;
   process_primary = 0;
   G4Process = 0;
   G4FinalProcess = 0;
   NTrTrajPts = 0;
   NDTrTrajPts = 0;
   DPdgCode = 0;
   DStartEnergy = 0;
   MidPosX = 0;
   MidPosY = 0;
   MidPosZ = 0;
   MidPx = 0;
   MidPy = 0;
   MidPz = 0;
   DPMidPosX = 0;
   DPMidPosY = 0;
   DPMidPosZ = 0;
   InteractionPoint = 0;
   InteractionPointType = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("evttime", &evttime, &b_evttime);
   fChain->SetBranchAddress("efield", efield, &b_efield);
   fChain->SetBranchAddress("t0", &t0, &b_t0);
   fChain->SetBranchAddress("ntracks_reco", &ntracks_reco, &b_ntracks_reco);
   fChain->SetBranchAddress("track_start_x", &track_start_x, &b_track_start_x);
   fChain->SetBranchAddress("track_start_y", &track_start_y, &b_track_start_y);
   fChain->SetBranchAddress("track_start_z", &track_start_z, &b_track_start_z);
   fChain->SetBranchAddress("track_end_x", &track_end_x, &b_track_end_x);
   fChain->SetBranchAddress("track_end_y", &track_end_y, &b_track_end_y);
   fChain->SetBranchAddress("track_end_z", &track_end_z, &b_track_end_z);
   fChain->SetBranchAddress("track_length", &track_length, &b_track_length);
   fChain->SetBranchAddress("ntrack_hits", &ntrack_hits, &b_ntrack_hits);
   fChain->SetBranchAddress("track_xpos", &track_xpos, &b_track_xpos);
   fChain->SetBranchAddress("track_ypos", &track_ypos, &b_track_ypos);
   fChain->SetBranchAddress("track_zpos", &track_zpos, &b_track_zpos);
   fChain->SetBranchAddress("ind_track_hits", &ind_track_hits, &b_ind_track_hits);
   fChain->SetBranchAddress("ind_track_ke", &ind_track_ke, &b_ind_track_ke);
   fChain->SetBranchAddress("ind_track_wire", &ind_track_wire, &b_ind_track_wire);
   fChain->SetBranchAddress("ind_track_dedx", &ind_track_dedx, &b_ind_track_dedx);
   fChain->SetBranchAddress("ind_track_dqdx", &ind_track_dqdx, &b_ind_track_dqdx);
   fChain->SetBranchAddress("ind_track_rr", &ind_track_rr, &b_ind_track_rr);
   fChain->SetBranchAddress("ind_track_pitch_hit", &ind_track_pitch_hit, &b_ind_track_pitch_hit);
   fChain->SetBranchAddress("col_track_hits", &col_track_hits, &b_col_track_hits);
   fChain->SetBranchAddress("col_track_ke", &col_track_ke, &b_col_track_ke);
   fChain->SetBranchAddress("col_track_wire", &col_track_wire, &b_col_track_wire);
   fChain->SetBranchAddress("col_track_dedx", &col_track_dedx, &b_col_track_dedx);
   fChain->SetBranchAddress("col_track_dqdx", &col_track_dqdx, &b_col_track_dqdx);
   fChain->SetBranchAddress("col_track_rr", &col_track_rr, &b_col_track_rr);
   fChain->SetBranchAddress("col_track_pitch_hit", &col_track_pitch_hit, &b_col_track_pitch_hit);
   fChain->SetBranchAddress("trkpitch", trkpitch, &b_trkpitch);
   fChain->SetBranchAddress("num_hit", &num_hit, &b_num_hits);
   fChain->SetBranchAddress("hit_channel", &hit_channel, &b_hit_channel);
   fChain->SetBranchAddress("hit_integral", &hit_integral, &b_hit_integral);
   fChain->SetBranchAddress("no_primaries", &no_primaries, &b_no_primaries);
   fChain->SetBranchAddress("geant_list_size", &geant_list_size, &b_geant_list_size);
   fChain->SetBranchAddress("primary_p", &primary_p, &b_primary);
   fChain->SetBranchAddress("PDG", &PDG, &b_PDG);
   fChain->SetBranchAddress("StartEnergy", &StartEnergy, &b_StartEnergy);
   fChain->SetBranchAddress("StartPx", &StartPx, &b_StartPx);
   fChain->SetBranchAddress("StartPy", &StartPy, &b_StartPy);
   fChain->SetBranchAddress("StartPz", &StartPz, &b_StartPz);
   fChain->SetBranchAddress("EndEnergy", &EndEnergy, &b_EndEnergy);
   fChain->SetBranchAddress("EndPx", &EndPx, &b_EndPx);
   fChain->SetBranchAddress("EndPy", &EndPy, &b_EndPy);
   fChain->SetBranchAddress("EndPz", &EndPz, &b_EndPz);
   fChain->SetBranchAddress("StartPointx", &StartPointx, &b_StartPointx);
   fChain->SetBranchAddress("StartPointy", &StartPointy, &b_StartPointy);
   fChain->SetBranchAddress("StartPointz", &StartPointz, &b_StartPointz);
   fChain->SetBranchAddress("EndPointx", &EndPointx, &b_EndPointx);
   fChain->SetBranchAddress("EndPointy", &EndPointy, &b_EndPointy);
   fChain->SetBranchAddress("EndPointz", &EndPointz, &b_EndPointz);
   fChain->SetBranchAddress("Process", &Process, &b_Process);
   fChain->SetBranchAddress("NumberDaughters", &NumberDaughters, &b_NumberDaughters);
   fChain->SetBranchAddress("primary_simChannel_num_voxel", &primary_simChannel_num_voxel, &b_primary_simChannel_num_voxel);
   fChain->SetBranchAddress("primary_simChannel_voxel_dr", &primary_simChannel_voxel_dr, &b_primary_simChannel_voxel_dr);
   fChain->SetBranchAddress("primary_simChannel_voxel_E", &primary_simChannel_voxel_E, &b_primary_simChannel_voxel_E);
   fChain->SetBranchAddress("primary_simChannel_voxel_e", &primary_simChannel_voxel_e, &b_primary_simChannel_voxel_e);
   fChain->SetBranchAddress("primary_simChannel_voxel_x", &primary_simChannel_voxel_x, &b_primary_simChannel_voxel_x);
   fChain->SetBranchAddress("primary_simChannel_voxel_y", &primary_simChannel_voxel_y, &b_primary_simChannel_voxel_y);
   fChain->SetBranchAddress("primary_simChannel_voxel_z", &primary_simChannel_voxel_z, &b_primary_simChannel_voxel_z);
   fChain->SetBranchAddress("primary_num_simChannel", &primary_num_simChannel, &b_primary_num_simChannel);
   fChain->SetBranchAddress("primary_simChannel", &primary_simChannel, &b_primary_simChannel);
   fChain->SetBranchAddress("primary_simChannel_dr", &primary_simChannel_dr, &b_primary_simChannel_dr);
   fChain->SetBranchAddress("primary_simChannel_E", &primary_simChannel_E, &b_primary_simChannel_E);
   fChain->SetBranchAddress("primary_simChannel_e", &primary_simChannel_e, &b_primary_simChannel_e);
   fChain->SetBranchAddress("Mother", &Mother, &b_Mother);
   fChain->SetBranchAddress("TrackId", &TrackId, &b_TrackId);
   fChain->SetBranchAddress("process_primary", &process_primary, &b_process_primary);
   fChain->SetBranchAddress("G4Process", &G4Process, &b_G4Process);
   fChain->SetBranchAddress("G4FinalProcess", &G4FinalProcess, &b_G4FinalProcess);
   fChain->SetBranchAddress("NTrTrajPts", &NTrTrajPts, &b_NTrTrajPts);
   fChain->SetBranchAddress("NProtonDaughters", &NProtonDaughters, &b_NProtonDaughters);
   fChain->SetBranchAddress("NNeutronDaughters", &NNeutronDaughters, &b_NNeutronDaughters);
   fChain->SetBranchAddress("NDTrTrajPts", &NDTrTrajPts, &b_NDTrTrajPts);
   fChain->SetBranchAddress("DPdgCode", &DPdgCode, &b_DPdgCode);
   fChain->SetBranchAddress("DStartEnergy", &DStartEnergy, &b_DStartEnergy);
   fChain->SetBranchAddress("MidPosX", &MidPosX, &b_MidPosX);
   fChain->SetBranchAddress("MidPosY", &MidPosY, &b_MidPosY);
   fChain->SetBranchAddress("MidPosZ", &MidPosZ, &b_MidPosZ);
   fChain->SetBranchAddress("MidPx", &MidPx, &b_MidPx);
   fChain->SetBranchAddress("MidPy", &MidPy, &b_MidPy);
   fChain->SetBranchAddress("MidPz", &MidPz, &b_MidPz);
   fChain->SetBranchAddress("DPMidPosX", &DPMidPosX, &b_DPMidPosX);
   fChain->SetBranchAddress("DPMidPosY", &DPMidPosY, &b_DPMidPosY);
   fChain->SetBranchAddress("DPMidPosZ", &DPMidPosZ, &b_DPMidPosZ);
   fChain->SetBranchAddress("InteractionPoint", &InteractionPoint, &b_InteractionPoint);
   fChain->SetBranchAddress("InteractionPointType", &InteractionPointType, &b_InteractionPointType);
   Notify();
}

Bool_t blah::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void blah::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t blah::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

// ### functions here ###
// ### some of these are because I'm lazy ###
// ### some of thees are because I like readable code ###
// ### you can test them out if you want or you can trust me ###
double blah::distance(double a[], double b[])
  {
  //double distance = sqrt((b[0] - a[0])**2 + (b[1] - a[1])**2 + (b[2] - a[2])**2);
  double distance = sqrt(pow( (b[0] - a[0]), 2) + pow( (b[1] - a[1]), 2) + pow( (b[2] - a[2]), 2));
  return distance;
  } //<---End distance() Function


bool blah::is_in_tpc(int spc_pt)
  {
  if( ( (*MidPosZ)[0][spc_pt] <   2. || (*MidPosZ)[0][spc_pt] > 88.)
   || ( (*MidPosY)[0][spc_pt] < -19. || (*MidPosY)[0][spc_pt] > 19.)
   || ( (*MidPosX)[0][spc_pt] <   1. || (*MidPosZ)[0][spc_pt] > 46.5) )
   {return false;}
  else
   {return true;}
  }

double blah::get_angle(double a[], double b[])
  {
  //double mag_a  = sqrt( a[0]**2 + a[1]**2 + a[2]**2);
  double mag_a  = sqrt( pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2));
  //double mag_b  = sqrt( b[0]**2 + b[1]**2 + b[2]**2);
  double mag_b  = sqrt( pow(b[0], 2) + pow(b[1], 2) + pow(b[2], 2));
  double dot_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  double theta  = acos( (dot_ab) / (mag_a*mag_b) );
  return theta*(180. /  3.14159265);
  }

int blah::get_primary(double mc_tpc_x, double mc_tpc_y, double mc_vx, double mc_vy, double mc_vz)
  {
  int primary = 999;
  // think about adding earliest space point z as well...
  for(size_t nTPCtrk = 0; nTPCtrk < ntracks_reco; nTPCtrk++)
    {
      //std::cout << "new reco track... " << std::endl;
      double tt = 0;
      double montecarlo_match [] = {mc_vx, mc_vy, mc_vz};
      //for(int i=0;i<3;i++){montecarlo_match[i]=mc_v[i];}
      double rplanar_x = 0; double rplanar_y = 0;
      //double mplanar_x = mc_int[0]; double mplanar_y = mc_int[1]; double mplanar_z = mc_int[2];
      double mplanar_x = mc_tpc_x; double mplanar_y = mc_tpc_y;
      double dx = 999; double dy = 999; double dz = 999; double dt = 0;
      double rp0[] = {(*track_xpos)[nTPCtrk][0], (*track_ypos)[nTPCtrk][0], (*track_zpos)[nTPCtrk][0]};
      double rp1[] = {(*track_xpos)[nTPCtrk][1], (*track_ypos)[nTPCtrk][1], (*track_zpos)[nTPCtrk][1]};
      double reco_match[]  = {rp1[0] - rp0[0],  rp1[1] - rp0[1],  rp1[2] - rp0[2]};
      tt = (0 - rp1[2]) / (reco_match[2]);
      rplanar_x = rp1[0] + tt*reco_match[0];
      rplanar_y = rp1[1] + tt*reco_match[1];
      dx = mplanar_x - rplanar_x;
      dy = mplanar_y - rplanar_y;
      dt = get_angle(montecarlo_match, reco_match);
      //std::cout << "(dx, dy, dtheta): (" << dx << ", " << dy << ", " << dt << ")" << std::endl; 
      if( abs(dx) < 6. && abs(dy) < 6. && dt < 16. )
        {
        //std::cout << "primary match!" << std::endl;
        //std::cout << "(dx, dy, dz, dtheta): (" << dx << ", " << dy << ", " << ", " << dt << ")" << std::endl; 
        primary = nTPCtrk;
        } 
    }//<--- End reco track loop
  //std::cout << primary << std::endl;
  return primary;
  }//<--- End get primary function

int blah::get_primary_track_id()
  {
  int primary_track_id = 0;
  for(size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
    if((*process_primary)[iG4] == 0){continue;}//<--- 0-nonprimary, 1-primary
    primary_track_id = (*TrackId)[iG4];
    }
  return primary_track_id;
  }//<--- End get primary track id function


int blah::get_primary_process()
  {
  int primary = get_primary_track_id();
  int mc_primary_process = 0;
  for(size_t iG4 = 0; iG4 < geant_list_size; iG4++)
    {
    if((*process_primary)[iG4] == 0)//<--- 0-nonprimary, 1-primary
      {
      if((*Mother)[iG4] != primary){continue;} 
      else{mc_primary_process = (*Process)[iG4];}
      //std::cout << G4Process->at(iG4) << std::endl;; 
      }
    }
  return mc_primary_process;
  }//<---End get primary process fucntion

#endif // #ifdef blah_cxx
